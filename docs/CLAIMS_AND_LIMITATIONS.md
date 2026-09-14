# Claims and limitations

## Claims supported by this repository

### Claim 1 — representation sensitivity

Exact behavior-preserving state refinement can change rank-based clear-CE status while leaving the exactly lumped dynamics and Effective Information unchanged.

### Claim 2 — exact empirical rank is fragile

A population TPM can have exact low rank while its finite-sample empirical estimate is generically full rank.

### Claim 3 — cutoff location and cutoff existence are different questions

A procedure that is constructed to return some internal cutoff does not establish that the spectrum contains a scientifically meaningful cutoff.

### Claim 4 — resolved modes are a safer one-sided statement

The cross-split diagnostic provides positive evidence for some modes without interpreting unresolved modes as absent.

### Claim 5 — spectrum does not certify a designated partition

Two TPMs can have the same singular values while a fixed partition is exactly
lumpable for one and not for the other.

### Claim 6 — full rank and exact macro-dynamics can coexist

A full-rank 48-state heat-kernel TPM admits the exact nested Markov hierarchy
\(48\rightarrow8\rightarrow2\).

### Claim 7 — resolved dimension is not macrostate cardinality

The number of resolved spectral modes changes with observation time and sample
size while the exact triadic hierarchy remains unchanged.

## Claims explicitly *not* supported

This repository does not establish that:

- SVD causal emergence is wrong;
- causal emergence is purely representational;
- `r_resolved` is a formal confidence lower bound;
- `r_resolved` is optimal;
- all physical refinements that look like state cloning are scientifically equivalent;
- spectral redundancy alone is sufficient for a scientifically useful macroscale;
- exact lumpability implies positive causal emergence;
- a valid macro-partition is unique;
- the exploratory triadic clustering baseline represents all possible
  macrostate-discovery methods.

## Statistical limitations

The repeated count splits quantify stability conditional on the observed finite count matrix. The 5th percentile used for resolved modes is **not** a theorem guaranteeing nominal population coverage.

A formal treatment would require a statistical theory for singular/eigenvalue inference under multinomially estimated Markov transition operators.

## Benchmark limitations

- synthetic TPMs;
- observed dimensions \(N=16\) and \(N=48\) in the two benchmark components;
- row-stratified sampling;
- a small set of full-rank control families;
- no real-world dataset in the final claim;
- one highly symmetric 48-state construction in the triadic extension;
- no theorem of blind macro-partition identifiability.

## Why failed exploratory tests are not headline results

Earlier exploratory work tested automatic epsilon selection, exact rank intervals, and other heuristics. Several failed under stronger controls.

Those failures motivated the final one-sided claim, but they are not presented as successful methods.

See `archive/EXPLORATORY_HISTORY.md`.
