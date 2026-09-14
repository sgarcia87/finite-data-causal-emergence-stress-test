#!/usr/bin/env python3
"""Triadic control: exact hierarchy versus finite-data spectral resolution.

This stage constructs the 48-state triadic Cayley graph, defines the explicit
heat-kernel dynamics P_tau = exp(-tau L), and verifies two nested exact Markov
quotients: 48 -> 8 triadic families -> 2 classes of the invariant t.

The population TPM is full rank for every finite tau.  We then apply the same
split-sample resolved-mode diagnostic used by the finite-data pilot.  The goal
is not to estimate a macrostate count: it is to test whether resolved spectral
modes, exact lumpability, and macro cardinality can be treated as equivalent.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np


def enumerate_triads4() -> np.ndarray:
    states = np.array(list(itertools.product((-1, 1), repeat=12)), dtype=np.int8)
    states = states.reshape(-1, 3, 4)
    balanced = np.all(states.sum(axis=2) == 0, axis=1)
    gram = np.einsum("aik,ajk->aij", states.astype(int), states.astype(int))
    orthogonal = np.all(gram == 4 * np.eye(3, dtype=int), axis=(1, 2))
    return states[balanced & orthogonal]


def adjacency(triads: np.ndarray) -> np.ndarray:
    lookup = {tuple(t.ravel()): i for i, t in enumerate(triads)}
    graph = np.zeros((len(triads), len(triads)), dtype=int)
    for i, triad in enumerate(triads):
        for block in range(3):
            neighbour = triad.copy()
            neighbour[block] *= -1
            graph[i, lookup[tuple(neighbour.ravel())]] = 1
        for left, right in itertools.combinations(range(3), 2):
            neighbour = triad.copy()
            neighbour[[left, right]] = neighbour[[right, left]]
            graph[i, lookup[tuple(neighbour.ravel())]] = 1
    if not np.array_equal(graph, graph.T) or not np.all(graph.sum(axis=1) == 6):
        raise AssertionError("unexpected triadic graph")
    return graph


def triadic_partitions(triads: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    family_keys = [tuple(sorted(map(tuple, triad))) for triad in triads]
    family_map = {key: i for i, key in enumerate(sorted(set(family_keys)))}
    family = np.array([family_map[key] for key in family_keys], dtype=int)
    t_class = np.all(triads == 1, axis=1).sum(axis=1).astype(int)
    return family, t_class


def membership(labels: np.ndarray) -> np.ndarray:
    _, normalized = np.unique(labels, return_inverse=True)
    return np.eye(normalized.max() + 1, dtype=int)[normalized]


def exact_generator_quotient(laplacian: np.ndarray, labels: np.ndarray) -> np.ndarray:
    c = membership(labels)
    lc = laplacian @ c
    quotient = np.array([lc[labels == group][0] for group in range(c.shape[1])])
    if not np.array_equal(lc, c @ quotient):
        raise AssertionError("partition is not an exact generator quotient")
    return quotient


def heat_kernel(laplacian: np.ndarray, tau: float) -> np.ndarray:
    values, vectors = np.linalg.eigh(laplacian)
    p = (vectors * np.exp(-tau * values)) @ vectors.T
    p[np.abs(p) < 1e-15] = 0.0
    p = np.maximum(p, 0.0)
    return p / p.sum(axis=1, keepdims=True)


def quotient_and_closure(p: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, float]:
    c = membership(labels).astype(float)
    sizes = c.sum(axis=0)
    macro = (c.T @ p @ c) / sizes[:, None]
    residual = float(np.abs(p @ c - c @ macro).sum(axis=1).max() / 2)
    return macro, residual


def effective_information(p: np.ndarray) -> float:
    output = p.mean(axis=0)
    output_rows = np.broadcast_to(output, p.shape)
    mask = p > 0
    terms = p[mask] * np.log2(p[mask] / output_rows[mask])
    return float(terms.sum() / p.shape[0])


def sample_counts(p: np.ndarray, samples_per_state: int, rng: np.random.Generator) -> np.ndarray:
    return np.array([rng.multinomial(samples_per_state, row) for row in p])


def resolved_count(counts: np.ndarray, rng: np.random.Generator, splits: int = 60) -> int | None:
    """Split diagnostic transcribed from the published finite-data benchmark.

    It is a conservative resolved-mode count, not an estimator of the number
    of valid macrostates and not a formal confidence bound.
    """
    n = counts.shape[0]
    left = rng.binomial(counts[None, :, :], 0.5, size=(splits, n, n))
    right = counts[None, :, :] - left
    left_n = left.sum(axis=2, keepdims=True)
    right_n = right.sum(axis=2, keepdims=True)
    good = (left_n[:, :, 0] > 0).all(axis=1) & (right_n[:, :, 0] > 0).all(axis=1)
    if good.sum() < 20:
        return None
    left = left[good] / left_n[good]
    right = right[good] / right_n[good]
    difference = (left - right) / 2
    signal = (left.transpose(0, 2, 1) @ right + right.transpose(0, 2, 1) @ left) / 2
    noise_matrix = difference.transpose(0, 2, 1) @ difference
    eigen_signal = np.linalg.eigvalsh(signal)[:, ::-1]
    noise = np.maximum(np.linalg.eigvalsh(noise_matrix)[:, -1], 0)
    margins = (eigen_signal - noise[:, None]) / np.maximum(eigen_signal[:, :1], 1e-15)
    lower = np.quantile(margins, 0.05, axis=0)
    resolved = 0
    for value in lower:
        if value <= 0:
            break
        resolved += 1
    return resolved


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(output_dir: Path, samples: list[int], repeats: int, splits: int, seed: int) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    triads = enumerate_triads4()
    graph = adjacency(triads)
    laplacian = 6 * np.eye(48, dtype=int) - graph
    family, t_class = triadic_partitions(triads)

    family_laplacian = exact_generator_quotient(laplacian, family)
    t_laplacian = exact_generator_quotient(laplacian, t_class)
    if not all(len(set(t_class[family == group])) == 1 for group in range(8)):
        raise AssertionError("the two-state quotient must coarsen the eight-family quotient")

    population_rows: list[dict] = []
    finite_rows: list[dict] = []
    rng = np.random.default_rng(seed)
    for tau in (0.1, 0.4, 1.0):
        p = heat_kernel(laplacian, tau)
        q8, closure8 = quotient_and_closure(p, family)
        q2, closure2 = quotient_and_closure(p, t_class)
        singular = np.linalg.svd(p, compute_uv=False)
        ei_micro = effective_information(p)
        for name, matrix, closure in (
            ("micro48", p, 0.0),
            ("family8", q8, closure8),
            ("t2", q2, closure2),
        ):
            ei = effective_information(matrix)
            population_rows.append({
                "tau": tau,
                "partition": name,
                "states": matrix.shape[0],
                "micro_algebraic_rank": 48,
                "micro_sigma_min": float(singular[-1]),
                "population_closure_max_row_tv": closure,
                "effective_information_bits": ei,
                "delta_ei_vs_micro_bits": ei - ei_micro,
                "effectiveness_ei_over_log2_states": ei / np.log2(matrix.shape[0]),
            })

        for sample_size in samples:
            for repeat in range(repeats):
                counts = sample_counts(p, sample_size, rng)
                value = resolved_count(counts, rng, splits=splits)
                finite_rows.append({
                    "tau": tau,
                    "samples_per_state": sample_size,
                    "repeat": repeat,
                    "resolved_count": "" if value is None else value,
                })

    summary_rows: list[dict] = []
    for tau in (0.1, 0.4, 1.0):
        for sample_size in samples:
            values = [
                int(row["resolved_count"])
                for row in finite_rows
                if row["tau"] == tau
                and row["samples_per_state"] == sample_size
                and row["resolved_count"] != ""
            ]
            summary_rows.append({
                "tau": tau,
                "samples_per_state": sample_size,
                "valid_repeats": len(values),
                "median_resolved": float(np.median(values)) if values else "",
                "q05_resolved": float(np.quantile(values, 0.05)) if values else "",
                "q95_resolved": float(np.quantile(values, 0.95)) if values else "",
                "min_resolved": min(values) if values else "",
                "max_resolved": max(values) if values else "",
            })

    write_csv(output_dir / "triadic_hierarchy_population.csv", population_rows)
    write_csv(output_dir / "triadic_resolved_modes_runs.csv", finite_rows)
    write_csv(output_dir / "triadic_resolved_modes_summary.csv", summary_rows)
    result = {
        "scope": "Exact nested lumpability versus finite-data resolved spectral modes.",
        "microstates": 48,
        "edges": int(graph.sum() // 2),
        "micro_graph": "Q3 Cartesian-product K3,3",
        "exact_nested_hierarchy": [48, 8, 2],
        "family8_quotient_laplacian": family_laplacian.tolist(),
        "t2_quotient_laplacian": t_laplacian.tolist(),
        "population": population_rows,
        "finite_summary": summary_rows,
        "interpretation": [
            "The microscopic TPM is full rank for every finite tau tested.",
            "The 8-state and 2-state partitions are exact for every tau.",
            "Resolved-mode count changes with tau and sample size while the exact hierarchy does not.",
            "Resolved spectral dimension is therefore not a macrostate-count estimator.",
            "Neither exact lumpability nor normalized effectiveness establishes positive causal emergence in absolute EI.",
        ],
    }
    (output_dir / "triadic_hierarchy_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--samples", type=int, nargs="+", default=[100, 500, 5000])
    parser.add_argument("--repeats", type=int, default=100)
    parser.add_argument("--splits", type=int, default=60)
    parser.add_argument("--seed", type=int, default=20260915)
    args = parser.parse_args()
    result = run(args.output, args.samples, args.repeats, args.splits, args.seed)
    print(json.dumps(result["finite_summary"], indent=2))


if __name__ == "__main__":
    main()
