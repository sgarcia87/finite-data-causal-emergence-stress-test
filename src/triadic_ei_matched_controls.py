#!/usr/bin/env python3
"""Triadic Control 5: matched-partition and temporal-null audit for Effective Information."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd


def factors(eta: float, a: float = 0.25) -> tuple[np.ndarray, np.ndarray]:
    z = np.exp(-2 * a)
    bit = np.array([[1 + z, 1 - z], [1 - z, 1 + z]]) / 2
    q = np.kron(np.kron(bit, bit), bit)
    j = np.ones((6, 6)) / 6
    v = np.array([1, 1, 1, -1, -1, -1]) / np.sqrt(6)
    vv = np.outer(v, v)
    r = j + eta * (np.eye(6) - j - vv) + eta**2 * vv
    return q, r


def make_models() -> dict[str, np.ndarray]:
    models = {}
    for eta in (0.0, 0.02, 0.2, 0.8):
        q, r = factors(eta)
        models[f"eta_{eta:g}"] = np.kron(q, r)
    p = models["eta_0.02"]
    u = np.repeat(np.array([1, 1, 1, 1, -1, -1, -1, -1]), 6) / np.sqrt(48)
    v = np.kron(np.ones(8) / np.sqrt(8), np.array([1, -1, 0, 0, 0, 0]) / np.sqrt(2))
    theta = 0.01
    o = (
        np.eye(48)
        + (np.cos(theta) - 1) * (np.outer(u, u) + np.outer(v, v))
        + np.sin(theta) * (np.outer(v, u) - np.outer(u, v))
    )
    models["isospectral_broken_8partition"] = o @ p @ o.T
    return models


def membership(labels: np.ndarray, groups: int = 8) -> np.ndarray:
    c = np.zeros((len(labels), groups), dtype=float)
    c[np.arange(len(labels)), labels] = 1.0
    return c


def macro_tpm(p: np.ndarray, labels: np.ndarray) -> np.ndarray:
    c = membership(labels)
    sizes = c.sum(axis=0)
    if np.any(sizes == 0):
        raise ValueError("empty macro group")
    return (c.T @ p @ c) / sizes[:, None]


def closure_tv(p: np.ndarray, labels: np.ndarray) -> float:
    c = membership(labels)
    q = macro_tpm(p, labels)
    return float(np.abs(p @ c - c @ q).sum(axis=1).max() / 2)


def effective_information(p: np.ndarray) -> float:
    output = p.mean(axis=0)
    rows = np.broadcast_to(output, p.shape)
    mask = p > 0
    return float((p[mask] * np.log2(p[mask] / rows[mask])).sum() / p.shape[0])


def balanced_partition(rng: np.random.Generator) -> np.ndarray:
    return rng.permutation(np.repeat(np.arange(8), 6))


def rank_correlation(x: np.ndarray, y: np.ndarray) -> float:
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    return float(np.corrcoef(rx, ry)[0, 1])


def sample_counts(p: np.ndarray, samples: int, rng: np.random.Generator) -> np.ndarray:
    return np.array([rng.multinomial(samples, row / row.sum()) for row in p])


def shuffle_futures_exact(counts: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Destroy X_t--X_t+1 dependence while preserving both margins exactly."""
    remaining = counts.sum(axis=0).astype(np.int64)
    row_totals = counts.sum(axis=1).astype(np.int64)
    out = np.zeros_like(counts)
    for row, total in enumerate(row_totals[:-1]):
        draw = rng.multivariate_hypergeometric(remaining, int(total))
        out[row] = draw
        remaining -= draw
    out[-1] = remaining
    if not np.array_equal(out.sum(axis=0), counts.sum(axis=0)):
        raise AssertionError("future-state margin changed")
    if not np.array_equal(out.sum(axis=1), row_totals):
        raise AssertionError("origin-state margin changed")
    return out


def empirical_tpm(counts: np.ndarray) -> np.ndarray:
    return counts / counts.sum(axis=1, keepdims=True)


def macro_counts(counts: np.ndarray, labels: np.ndarray) -> np.ndarray:
    c = membership(labels)
    return c.T @ counts @ c


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(output: Path, random_partitions: int, repeats: int, seed: int) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    models = make_models()
    truth = np.repeat(np.arange(8), 6)
    population_rows = []
    random_rows = []
    null_rows = []

    for model_index, (name, p) in enumerate(models.items()):
        micro_ei = effective_information(p)
        truth_closure = closure_tv(p, truth)
        truth_macro_ei = effective_information(macro_tpm(p, truth))
        rng = np.random.default_rng(np.random.SeedSequence([seed, model_index, 1]))
        closures = np.empty(random_partitions)
        macro_eis = np.empty(random_partitions)
        for index in range(random_partitions):
            labels = balanced_partition(rng)
            closures[index] = closure_tv(p, labels)
            macro_eis[index] = effective_information(macro_tpm(p, labels))
            random_rows.append({
                "model": name,
                "control": index,
                "closure_tv": closures[index],
                "macro_ei_bits": macro_eis[index],
                "delta_ei_vs_micro_bits": macro_eis[index] - micro_ei,
            })
        best_ei = int(np.argmax(macro_eis))
        best_closure = int(np.argmin(closures))
        rho = rank_correlation(macro_eis, -closures)
        population_rows.append({
            "model": name,
            "micro_ei_bits": micro_ei,
            "designated_closure_tv": truth_closure,
            "designated_macro_ei_bits": truth_macro_ei,
            "designated_delta_ei_bits": truth_macro_ei - micro_ei,
            "designated_closure_better_fraction": float(np.mean(closures > truth_closure)),
            "designated_ei_better_fraction": float(np.mean(macro_eis < truth_macro_ei)),
            "random_median_closure_tv": float(np.median(closures)),
            "random_min_closure_tv": float(closures[best_closure]),
            "random_median_macro_ei_bits": float(np.median(macro_eis)),
            "random_max_macro_ei_bits": float(macro_eis[best_ei]),
            "closure_of_max_ei_random": float(closures[best_ei]),
            "ei_of_min_closure_random": float(macro_eis[best_closure]),
            "spearman_ei_vs_negative_closure": rho,
            "strong_alignment_abs_rho_ge_0_5": bool(abs(rho) >= 0.5),
        })

        for sample_size in (100, 500, 5000):
            for repeat in range(repeats):
                ss = np.random.SeedSequence([seed, model_index, sample_size, repeat, 2])
                data_ss, null_ss = ss.spawn(2)
                counts = sample_counts(p, sample_size, np.random.default_rng(data_ss))
                shuffled = shuffle_futures_exact(counts, np.random.default_rng(null_ss))
                row = {"model": name, "samples_per_state": sample_size, "repeat": repeat}
                for prefix, matrix in (("observed", counts), ("null", shuffled)):
                    micro = effective_information(empirical_tpm(matrix))
                    macro = effective_information(empirical_tpm(macro_counts(matrix, truth)))
                    row[f"{prefix}_micro_ei_bits"] = micro
                    row[f"{prefix}_macro_ei_bits"] = macro
                    row[f"{prefix}_delta_ei_bits"] = macro - micro
                null_rows.append(row)

    write_csv(output / "triadic_ei_population_partition_controls.csv", population_rows)
    write_csv(output / "triadic_ei_random_partition_controls.csv", random_rows)
    write_csv(output / "triadic_ei_finite_null_runs.csv", null_rows)
    null_summary = []
    for name in models:
        for sample_size in (100, 500, 5000):
            subset = [r for r in null_rows if r["model"] == name and r["samples_per_state"] == sample_size]
            null_summary.append({
                "model": name,
                "samples_per_state": sample_size,
                "independent_repeats": len(subset),
                "median_observed_delta_ei_bits": float(np.median([r["observed_delta_ei_bits"] for r in subset])),
                "observed_positive_delta_rate": float(np.mean([r["observed_delta_ei_bits"] > 0 for r in subset])),
                "median_null_delta_ei_bits": float(np.median([r["null_delta_ei_bits"] for r in subset])),
                "null_positive_delta_rate": float(np.mean([r["null_delta_ei_bits"] > 0 for r in subset])),
            })
    write_csv(output / "triadic_ei_finite_null_summary.csv", null_summary)
    result = {
        "status": "Exploratory protocol fixed before first run.",
        "seed": seed,
        "balanced_random_partitions_per_model": random_partitions,
        "finite_null_repeats": repeats,
        "population": population_rows,
        "finite_null_summary": null_summary,
    }
    (output / "triadic_ei_matched_controls_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/triadic"))
    parser.add_argument("--random-partitions", type=int, default=2000)
    parser.add_argument("--repeats", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260914)
    args = parser.parse_args()
    result = run(args.output, args.random_partitions, args.repeats, args.seed)
    print(json.dumps(result["population"], indent=2))

