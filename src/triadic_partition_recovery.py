#!/usr/bin/env python3
"""Triadic control: partition recovery with known macrostate count.

The algorithm is told that eight groups are requested, but not which of the 48
states belong together. It clusters a spectral embedding of a training TPM,
freezes the partition, and evaluates it on an independent test TPM. Hidden
family labels and population matrices are used only by the evaluator.

The exact family and its isospectral control are defined explicitly below.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


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


def population_closure_residual(p: np.ndarray, labels: np.ndarray) -> float:
    c = membership(labels)
    sizes = c.sum(axis=0)
    if np.any(sizes == 0):
        return float("inf")
    macro = (c.T @ p @ c) / sizes[:, None]
    return float(np.abs(p @ c - c @ macro).sum(axis=1).max() / 2)


def sample_counts(p: np.ndarray, samples_per_state: int, rng: np.random.Generator) -> np.ndarray:
    return np.array([rng.multinomial(samples_per_state, row / row.sum()) for row in p])


def empirical_tpm(counts: np.ndarray) -> np.ndarray:
    return counts / counts.sum(axis=1, keepdims=True)


def discover_partition(p_train: np.ndarray, seed: int, groups: int = 8) -> np.ndarray:
    """Simple spectral-clustering baseline conditional on a known group count."""
    symmetric = (p_train + p_train.T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    embedding = vectors[:, np.argsort(values)[-groups:]]
    norms = np.linalg.norm(embedding, axis=1, keepdims=True)
    embedding = embedding / np.maximum(norms, 1e-15)
    return KMeans(n_clusters=groups, n_init=30, random_state=seed).fit_predict(embedding)


def heldout_macro_metrics(
    train_counts: np.ndarray,
    test_counts: np.ndarray,
    labels: np.ndarray,
) -> tuple[float, float]:
    """Return held-out closure TV and macro-transition log loss."""
    c = membership(labels)
    group_counts = c.T @ train_counts @ c
    macro = (group_counts + 0.5) / (group_counts.sum(axis=1, keepdims=True) + 0.5 * c.shape[1])

    p_test = empirical_tpm(test_counts)
    observed_macro_rows = p_test @ c
    predicted_macro_rows = c @ macro
    tv = float(np.abs(observed_macro_rows - predicted_macro_rows).sum(axis=1).mean() / 2)

    test_macro_counts = test_counts @ c
    probabilities = predicted_macro_rows
    log_loss = float(-(test_macro_counts * np.log(np.maximum(probabilities, 1e-15))).sum() / test_macro_counts.sum())
    return tv, log_loss


def random_balanced_partition(rng: np.random.Generator) -> np.ndarray:
    return rng.permutation(np.repeat(np.arange(8), 6))


def run(
    output_dir: Path,
    samples: list[int],
    repeats: int,
    random_controls: int,
    seed: int,
) -> dict:
    models = make_models()
    truth = np.repeat(np.arange(8), 6)
    rows = []

    for model_index, (model_name, p) in enumerate(models.items()):
        for sample_size in samples:
            for repeat in range(repeats):
                seed_sequence = np.random.SeedSequence([seed, model_index, sample_size, repeat])
                train_ss, test_ss, control_ss = seed_sequence.spawn(3)
                train_counts = sample_counts(p, sample_size, np.random.default_rng(train_ss))
                test_counts = sample_counts(p, sample_size, np.random.default_rng(test_ss))
                p_train = empirical_tpm(train_counts)
                discovered = discover_partition(p_train, seed=int(train_ss.generate_state(1)[0]))

                discovered_tv, discovered_loss = heldout_macro_metrics(
                    train_counts, test_counts, discovered
                )
                truth_tv, truth_loss = heldout_macro_metrics(train_counts, test_counts, truth)

                control_rng = np.random.default_rng(control_ss)
                random_metrics = [
                    heldout_macro_metrics(
                        train_counts,
                        test_counts,
                        random_balanced_partition(control_rng),
                    )
                    for _ in range(random_controls)
                ]
                random_tvs = [x[0] for x in random_metrics]
                random_losses = [x[1] for x in random_metrics]

                rows.append({
                    "model": model_name,
                    "samples_per_state": sample_size,
                    "repeat": repeat,
                    "discovered_ari_to_hidden_families": float(adjusted_rand_score(truth, discovered)),
                    "discovered_group_count": int(len(np.unique(discovered))),
                    "discovered_min_group_size": int(np.bincount(discovered).min()),
                    "discovered_max_group_size": int(np.bincount(discovered).max()),
                    "discovered_population_closure_tv": population_closure_residual(p, discovered),
                    "hidden_population_closure_tv": population_closure_residual(p, truth),
                    "discovered_heldout_closure_tv": discovered_tv,
                    "hidden_heldout_closure_tv": truth_tv,
                    "random_median_heldout_closure_tv": float(np.median(random_tvs)),
                    "random_best_heldout_closure_tv": float(np.min(random_tvs)),
                    "discovered_heldout_macro_log_loss": discovered_loss,
                    "hidden_heldout_macro_log_loss": truth_loss,
                    "random_median_heldout_macro_log_loss": float(np.median(random_losses)),
                })

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "triadic_partition_recovery_runs.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for model_name in models:
        for sample_size in samples:
            subset = [
                row for row in rows
                if row["model"] == model_name and row["samples_per_state"] == sample_size
            ]
            summary.append({
                "model": model_name,
                "samples_per_state": sample_size,
                "independent_repeats": len(subset),
                "median_ari": float(np.median([r["discovered_ari_to_hidden_families"] for r in subset])),
                "exact_partition_rate": float(np.mean([
                    r["discovered_ari_to_hidden_families"] == 1.0 for r in subset
                ])),
                "median_discovered_population_closure_tv": float(np.median([
                    r["discovered_population_closure_tv"] for r in subset
                ])),
                "median_hidden_population_closure_tv": float(np.median([
                    r["hidden_population_closure_tv"] for r in subset
                ])),
                "median_discovered_heldout_closure_tv": float(np.median([
                    r["discovered_heldout_closure_tv"] for r in subset
                ])),
                "median_hidden_heldout_closure_tv": float(np.median([
                    r["hidden_heldout_closure_tv"] for r in subset
                ])),
                "median_random_heldout_closure_tv": float(np.median([
                    r["random_median_heldout_closure_tv"] for r in subset
                ])),
                "median_discovered_log_loss": float(np.median([
                    r["discovered_heldout_macro_log_loss"] for r in subset
                ])),
                "median_hidden_log_loss": float(np.median([
                    r["hidden_heldout_macro_log_loss"] for r in subset
                ])),
            })

    with (output_dir / "triadic_partition_recovery_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)

    result = {
        "status": "Exploratory Stage 5; k=8 is supplied, family memberships are hidden.",
        "seed": seed,
        "samples_per_state": samples,
        "independent_repeats": repeats,
        "random_partitions_per_run": random_controls,
        "summary": summary,
    }
    (output_dir / "triadic_partition_recovery_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--samples", type=int, nargs="+", default=[100, 500, 5000])
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--random-controls", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260914)
    args = parser.parse_args()
    run(args.output, args.samples, args.repeats, args.random_controls, args.seed)
