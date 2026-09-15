#!/usr/bin/env python3
"""Audit codon partitions against the Kosiol-Holmes-Goldman ECM models."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits


ROOT_SEED = 20260915
HORIZONS = (0.01, 0.05, 0.10, 0.50, 1.00)
TARGET_BY_K = {16: "B16", 20: "AA20", 23: "L23"}
SOURCE_HASHES = {
    "ECMrest.dat": "51f97e4310a41ba3a275446f5f46d3a0dfc40da875e01f62e13b38bb76d7b141",
    "ECMunrest.dat": "4a4701ec245b128c91c38bf36678ad2d6904ebe37c61fbaa6f1cdeaacc246de0",
}

GENETIC_CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"STOP","TAG":"STOP","TGT":"C","TGC":"C","TGA":"STOP","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize_labels(keys: list[object]) -> np.ndarray:
    mapping = {key: i for i, key in enumerate(sorted(set(keys)))}
    return np.array([mapping[key] for key in keys], dtype=int)


def partitions(codons: list[str]) -> dict[str, np.ndarray]:
    if any(GENETIC_CODE[c] == "STOP" for c in codons):
        raise ValueError("ECM input must contain only sense codons")
    return {
        "B16": normalize_labels([c[:2] for c in codons]),
        "AA20": normalize_labels([GENETIC_CODE[c] for c in codons]),
        "L23": normalize_labels([(c[:2], GENETIC_CODE[c]) for c in codons]),
    }


def refines(fine: np.ndarray, coarse: np.ndarray) -> bool:
    return all(len(np.unique(coarse[fine == group])) == 1 for group in np.unique(fine))


def parse_ecm(path: Path) -> tuple[list[str], np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    exchangeability = np.zeros((61, 61), dtype=float)
    for row in range(1, 61):
        values = [float(x) for x in lines[row - 1].split()]
        if len(values) != row:
            raise ValueError(f"Expected {row} lower-triangle entries on line {row}")
        exchangeability[row, :row] = values
        exchangeability[:row, row] = values

    frequencies = None
    codons: list[str] = []
    for line in lines[60:]:
        tokens = line.split()
        if frequencies is None and len(tokens) == 61:
            try:
                candidate = np.array([float(x) for x in tokens], dtype=float)
            except ValueError:
                pass
            else:
                frequencies = candidate
                continue
        codon_tokens = [x for x in tokens if re.fullmatch(r"[ACGT]{3}", x)]
        codons.extend(codon_tokens)
        if len(codons) == 61:
            break

    if frequencies is None or len(codons) != 61:
        raise ValueError("Could not parse 61 frequencies and codons")
    frequencies /= frequencies.sum()
    return codons, exchangeability, frequencies


def build_generator(exchangeability: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
    q = exchangeability * frequencies[None, :]
    np.fill_diagonal(q, 0.0)
    np.fill_diagonal(q, -q.sum(axis=1))
    mean_rate = float(-frequencies @ np.diag(q))
    if mean_rate <= 0:
        raise ValueError("Non-positive equilibrium mean rate")
    return q / mean_rate


def transition_matrix(q: np.ndarray, horizon: float) -> np.ndarray:
    p = expm(horizon * q)
    p[np.abs(p) < 1e-15] = 0.0
    p /= p.sum(axis=1, keepdims=True)
    return p


def closure_tv(p: np.ndarray, labels: np.ndarray) -> float:
    groups = int(labels.max() + 1)
    aggregate = np.column_stack([p[:, labels == g].sum(axis=1) for g in range(groups)])
    macro_rows = np.vstack([aggregate[labels == g].mean(axis=0) for g in range(groups)])
    return float(0.5 * np.abs(aggregate - macro_rows[labels]).sum(axis=1).max())


def jump_profiles(p: np.ndarray) -> np.ndarray:
    profiles = p.copy()
    np.fill_diagonal(profiles, 0.0)
    profiles /= profiles.sum(axis=1, keepdims=True)
    return profiles


def discover_jump(p: np.ndarray, groups: int, seed: int) -> np.ndarray:
    with threadpool_limits(limits=1):
        return KMeans(n_clusters=groups, n_init=20, random_state=seed).fit_predict(jump_profiles(p))


def discover_spectral(p: np.ndarray, frequencies: np.ndarray, groups: int, seed: int) -> np.ndarray:
    root = np.sqrt(frequencies)
    symmetric = root[:, None] * p / root[None, :]
    symmetric = (symmetric + symmetric.T) / 2
    with threadpool_limits(limits=1):
        values, vectors = np.linalg.eigh(symmetric)
        embedding = vectors[:, np.argsort(values)[-groups:]]
        embedding /= np.maximum(np.linalg.norm(embedding, axis=1, keepdims=True), 1e-15)
        return KMeans(n_clusters=groups, n_init=20, random_state=seed).fit_predict(embedding)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def category(median_ari: float, exact_rate: float) -> str:
    if exact_rate >= 0.90:
        return "EXACT"
    if median_ari >= 0.75:
        return "STRONG"
    if median_ari >= 0.40:
        return "PARTIAL"
    return "WEAK_NULL"


def run(source_dir: Path, output: Path, controls: int = 2000, repeats: int = 30,
        seed: int = ROOT_SEED) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    source_rows: list[dict] = []
    closure_rows: list[dict] = []
    recovery_rows: list[dict] = []
    recovery_summary: list[dict] = []
    model_cache: dict[str, tuple[list[str], np.ndarray, np.ndarray]] = {}

    for model_name, filename in (("ECMrest", "ECMrest.dat"), ("ECMunrest", "ECMunrest.dat")):
        path = source_dir / filename
        digest = sha256(path)
        if digest != SOURCE_HASHES[filename]:
            raise ValueError(f"Unexpected SHA256 for {filename}: {digest}")
        codons, exchangeability, frequencies = parse_ecm(path)
        q = build_generator(exchangeability, frequencies)
        truth = partitions(codons)
        detailed_balance = np.max(np.abs(frequencies[:, None] * q - frequencies[None, :] * q.T))
        source_rows.append({
            "model": model_name,
            "filename": filename,
            "sha256": digest,
            "codons": len(codons),
            "frequency_sum": float(frequencies.sum()),
            "generator_max_abs_row_sum": float(np.abs(q.sum(axis=1)).max()),
            "detailed_balance_max_abs": float(detailed_balance),
            "equilibrium_mean_rate": float(-frequencies @ np.diag(q)),
            "positive_offdiagonal_fraction": float(np.mean(exchangeability[np.triu_indices(61, 1)] > 0)),
        })
        model_cache[model_name] = (codons, q, frequencies)

        rng = np.random.default_rng(np.random.SeedSequence([seed, 100 if model_name == "ECMrest" else 200]))
        randomized = {
            name: [rng.permutation(labels) for _ in range(controls)]
            for name, labels in truth.items()
        }

        for horizon_index, horizon in enumerate(HORIZONS):
            p = transition_matrix(q, horizon)
            for target, labels in truth.items():
                observed = closure_tv(p, labels)
                null = np.array([closure_tv(p, x) for x in randomized[target]])
                better_equal = int(np.sum(null <= observed))
                closure_rows.append({
                    "model": model_name,
                    "horizon": horizon,
                    "partition": target,
                    "groups": int(labels.max() + 1),
                    "closure_tv": observed,
                    "exact_closure": bool(observed < 1e-10),
                    "matched_controls": controls,
                    "random_q05": float(np.quantile(null, 0.05)),
                    "random_median": float(np.median(null)),
                    "random_q95": float(np.quantile(null, 0.95)),
                    "random_better_or_equal": better_equal,
                    "one_sided_p": float((1 + better_equal) / (controls + 1)),
                    "exceptionality_percentile": float(100 * np.mean(null >= observed)),
                })

            for groups, target in TARGET_BY_K.items():
                for repeat in range(repeats):
                    sequence = np.random.SeedSequence([
                        seed, 1 if model_name == "ECMrest" else 2,
                        horizon_index, groups, repeat,
                    ])
                    jump_seed, spectral_seed = [int(x.generate_state(1)[0]) for x in sequence.spawn(2)]
                    discovered_by_method = {
                        "jump_profile": discover_jump(p, groups, jump_seed),
                        "spectral": discover_spectral(p, frequencies, groups, spectral_seed),
                    }
                    for method, discovered in discovered_by_method.items():
                        row = {
                            "model": model_name,
                            "horizon": horizon,
                            "requested_k": groups,
                            "matching_target": target,
                            "repeat": repeat,
                            "method": method,
                            "discovered_closure_tv": closure_tv(p, discovered),
                        }
                        for name, labels in truth.items():
                            row[f"ari_{name}"] = float(adjusted_rand_score(labels, discovered))
                        recovery_rows.append(row)

    for model_name in ("ECMrest", "ECMunrest"):
        for horizon in HORIZONS:
            for groups, target in TARGET_BY_K.items():
                for method in ("jump_profile", "spectral"):
                    subset = [r for r in recovery_rows if
                              r["model"] == model_name and r["horizon"] == horizon and
                              r["requested_k"] == groups and r["method"] == method]
                    matching = np.array([r[f"ari_{target}"] for r in subset])
                    exact_rate = float(np.mean(np.isclose(matching, 1.0, atol=1e-12)))
                    median_ari = float(np.median(matching))
                    recovery_summary.append({
                        "model": model_name,
                        "horizon": horizon,
                        "requested_k": groups,
                        "matching_target": target,
                        "method": method,
                        "repeats": repeats,
                        "median_matching_ari": median_ari,
                        "q05_matching_ari": float(np.quantile(matching, 0.05)),
                        "q95_matching_ari": float(np.quantile(matching, 0.95)),
                        "exact_recovery_rate": exact_rate,
                        "category": category(median_ari, exact_rate),
                        "median_discovered_closure_tv": float(np.median([r["discovered_closure_tv"] for r in subset])),
                    })

    closure_witnesses = []
    recovery_witnesses = []
    for horizon in HORIZONS:
        for target in TARGET_BY_K.values():
            rows = [r for r in closure_rows if r["horizon"] == horizon and r["partition"] == target]
            if len(rows) == 2 and all(r["one_sided_p"] <= 0.01 for r in rows):
                closure_witnesses.append({"horizon": horizon, "partition": target})
        for groups, target in TARGET_BY_K.items():
            for method in ("jump_profile", "spectral"):
                rows = [r for r in recovery_summary if r["horizon"] == horizon and
                        r["requested_k"] == groups and r["method"] == method]
                if len(rows) == 2 and all(r["median_matching_ari"] >= 0.40 for r in rows):
                    recovery_witnesses.append({"horizon": horizon, "partition": target, "method": method})

    codons, _, _ = model_cache["ECMrest"]
    truth = partitions(codons)
    partition_rows = [{
        "codon": codon,
        "amino_acid": GENETIC_CODE[codon],
        "B16": int(truth["B16"][i]),
        "AA20": int(truth["AA20"][i]),
        "L23": int(truth["L23"][i]),
    } for i, codon in enumerate(codons)]

    result = {
        "status": "Empirical ECM population audit",
        "root_seed": seed,
        "matched_controls_per_partition": controls,
        "kmeans_repeats": repeats,
        "horizons": list(HORIZONS),
        "partition_cardinalities": {name: int(labels.max() + 1) for name, labels in truth.items()},
        "refinement_checks": {
            "L23_refines_B16": refines(truth["L23"], truth["B16"]),
            "L23_refines_AA20": refines(truth["L23"], truth["AA20"]),
            "B16_refines_AA20": refines(truth["B16"], truth["AA20"]),
            "AA20_refines_B16": refines(truth["AA20"], truth["B16"]),
        },
        "gate": {
            "passed": bool(closure_witnesses or recovery_witnesses),
            "closure_witnesses": closure_witnesses,
            "recovery_witnesses": recovery_witnesses,
        },
        "source_validation": source_rows,
        "closure": closure_rows,
        "recovery_summary": recovery_summary,
    }

    write_csv(output / "source_validation.csv", source_rows)
    write_csv(output / "codon_partitions.csv", partition_rows)
    write_csv(output / "designated_closure.csv", closure_rows)
    write_csv(output / "recovery_runs.csv", recovery_rows)
    write_csv(output / "recovery_summary.csv", recovery_summary)
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("data/empirical_ecm"))
    parser.add_argument("--output", type=Path, default=Path("results/empirical_ecm"))
    parser.add_argument("--controls", type=int, default=2000)
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--seed", type=int, default=ROOT_SEED)
    args = parser.parse_args()
    summary = run(args.source_dir, args.output, args.controls, args.repeats, args.seed)
    print(json.dumps({"refinement_checks": summary["refinement_checks"], "gate": summary["gate"]}, indent=2))
