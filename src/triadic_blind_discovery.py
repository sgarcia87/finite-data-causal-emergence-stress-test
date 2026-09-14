#!/usr/bin/env python3
"""Triadic control: exploratory blind discovery of an exact hierarchy.

The baseline receives only a training TPM. It searches k=2,...,12, selects
the largest non-trivial raw eigengap of the symmetrized TPM, clusters the
leading k-dimensional spectral embedding, freezes the result, and evaluates
it on independent data and against the evaluator-only population structure.

Oracle-k=2 and oracle-k=8 runs diagnose whether failure comes from selecting k
or from recovering a partition once k is supplied. This is an exploratory
baseline designed after the exact-hierarchy control, not a preregistered or
universal discovery procedure.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

import triadic_exact_hierarchy as hierarchy


def empirical_tpm(counts: np.ndarray) -> np.ndarray:
    return counts / counts.sum(axis=1, keepdims=True)


def spectral_values_vectors(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    symmetric = (p + p.T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    order = np.argsort(values)[::-1]
    return values[order], vectors[:, order]


def select_k_by_raw_gap(values: np.ndarray, minimum: int = 2, maximum: int = 12) -> int:
    candidates = np.arange(minimum, maximum + 1)
    gaps = values[candidates - 1] - values[candidates]
    return int(candidates[np.argmax(gaps)])


def cluster_embedding(vectors: np.ndarray, k: int, seed: int) -> np.ndarray:
    embedding = vectors[:, :k]
    embedding /= np.maximum(np.linalg.norm(embedding, axis=1, keepdims=True), 1e-15)
    return KMeans(n_clusters=k, n_init=30, random_state=seed).fit_predict(embedding)


def heldout_closure_tv(train_counts: np.ndarray, test_counts: np.ndarray, labels: np.ndarray) -> float:
    c = hierarchy.membership(labels).astype(float)
    group_counts = c.T @ train_counts @ c
    macro = (group_counts + 0.5) / (group_counts.sum(axis=1, keepdims=True) + 0.5 * c.shape[1])
    observed = empirical_tpm(test_counts) @ c
    predicted = c @ macro
    return float(np.abs(observed - predicted).sum(axis=1).mean() / 2)


def one_candidate(
    method: str,
    k: int,
    vectors: np.ndarray,
    seed: int,
    population: np.ndarray,
    train_counts: np.ndarray,
    test_counts: np.ndarray,
    family: np.ndarray,
    t_class: np.ndarray,
    coordinate_cuts: list[np.ndarray],
) -> dict:
    labels = cluster_embedding(vectors.copy(), k, seed)
    closure = hierarchy.quotient_and_closure(population, labels)[1]
    return {
        "method": method,
        "selected_k": k,
        "population_closure_max_row_tv": closure,
        "population_exact_at_1e_minus_10": closure < 1e-10,
        "heldout_closure_mean_row_tv": heldout_closure_tv(train_counts, test_counts, labels),
        "ari_family8": adjusted_rand_score(family, labels) if k == 8 else "",
        "ari_t2": adjusted_rand_score(t_class, labels) if k == 2 else "",
        "max_ari_cube_coordinate2": (
            max(adjusted_rand_score(cut, labels) for cut in coordinate_cuts)
            if k == 2 else ""
        ),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(output_dir: Path, samples: list[int], repeats: int, seed: int) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    triads = hierarchy.enumerate_triads4()
    graph = hierarchy.adjacency(triads)
    laplacian = 6 * np.eye(48, dtype=int) - graph
    family, t_class = hierarchy.triadic_partitions(triads)
    coordinate_cuts = [((family >> bit) & 1) for bit in range(3)]
    rows: list[dict] = []

    for tau_index, tau in enumerate((0.1, 0.4, 1.0)):
        population = hierarchy.heat_kernel(laplacian, tau)
        for sample_size in samples:
            for repeat in range(repeats):
                sequences = np.random.SeedSequence([seed, tau_index, sample_size, repeat]).spawn(2)
                train_counts = hierarchy.sample_counts(
                    population, sample_size, np.random.default_rng(sequences[0])
                )
                test_counts = hierarchy.sample_counts(
                    population, sample_size, np.random.default_rng(sequences[1])
                )
                values, vectors = spectral_values_vectors(empirical_tpm(train_counts))
                blind_k = select_k_by_raw_gap(values)
                for method, k, offset in (
                    ("blind_raw_eigengap", blind_k, 0),
                    ("oracle_k2", 2, 1),
                    ("oracle_k8", 8, 2),
                ):
                    candidate = one_candidate(
                        method,
                        k,
                        vectors,
                        seed + 100000 * tau_index + 1000 * repeat + offset,
                        population,
                        train_counts,
                        test_counts,
                        family,
                        t_class,
                        coordinate_cuts,
                    )
                    candidate.update({
                        "tau": tau,
                        "samples_per_state": sample_size,
                        "repeat": repeat,
                    })
                    rows.append(candidate)

    summary: list[dict] = []
    methods = ("blind_raw_eigengap", "oracle_k2", "oracle_k8")
    for tau in (0.1, 0.4, 1.0):
        for sample_size in samples:
            for method in methods:
                subset = [
                    row for row in rows
                    if row["tau"] == tau
                    and row["samples_per_state"] == sample_size
                    and row["method"] == method
                ]
                ks = np.array([row["selected_k"] for row in subset])
                closures = np.array([row["population_closure_max_row_tv"] for row in subset])
                heldout = np.array([row["heldout_closure_mean_row_tv"] for row in subset])
                ari8 = [float(row["ari_family8"]) for row in subset if row["ari_family8"] != ""]
                ari2 = [float(row["ari_t2"]) for row in subset if row["ari_t2"] != ""]
                coordinate2 = [
                    float(row["max_ari_cube_coordinate2"])
                    for row in subset if row["max_ari_cube_coordinate2"] != ""
                ]
                counts = {str(k): int((ks == k).sum()) for k in sorted(set(ks.tolist()))}
                summary.append({
                    "tau": tau,
                    "samples_per_state": sample_size,
                    "method": method,
                    "runs": len(subset),
                    "selected_k_counts": json.dumps(counts, sort_keys=True),
                    "hierarchy_level_selection_rate": float(np.mean(np.isin(ks, [2, 8]))),
                    "population_exact_partition_rate": float(np.mean(closures < 1e-10)),
                    "median_population_closure_tv": float(np.median(closures)),
                    "median_heldout_closure_tv": float(np.median(heldout)),
                    "median_ari_family8": float(np.median(ari8)) if ari8 else "",
                    "median_ari_t2": float(np.median(ari2)) if ari2 else "",
                    "median_max_ari_cube_coordinate2": (
                        float(np.median(coordinate2)) if coordinate2 else ""
                    ),
                })

    write_csv(output_dir / "triadic_blind_discovery_runs.csv", rows)
    write_csv(output_dir / "triadic_blind_discovery_summary.csv", summary)
    result = {
        "scope": "Exploratory spectral baseline; not preregistered.",
        "hidden_exact_hierarchy": [48, 8, 2],
        "candidate_k_range": [2, 12],
        "selection": "largest raw eigengap of symmetrized training TPM",
        "summary": summary,
        "claim_ceiling": [
            "Failure applies to this baseline, not to all hierarchy-discovery methods.",
            "ARI to designated labels and exact closure are separate outcomes.",
            "Oracle-k diagnostics do not constitute blind discovery.",
        ],
    }
    (output_dir / "triadic_blind_discovery_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--samples", type=int, nargs="+", default=[100, 500, 5000])
    parser.add_argument("--repeats", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260916)
    args = parser.parse_args()
    result = run(args.output, args.samples, args.repeats, args.seed)
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
