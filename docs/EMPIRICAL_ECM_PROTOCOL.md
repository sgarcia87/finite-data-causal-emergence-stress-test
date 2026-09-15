# Empirical codon structure audit — protocol v0.1

Status: frozen before computing any partition metric or clustering result.

Date: 2026-09-15.

## Question

Do biologically interpretable partitions of the 61 sense codons leave a
detectable signature in empirical codon-substitution dynamics that were not
constructed from those partitions?

This is an empirical-model audit, not a claim about raw molecular trajectories
or the evolutionary origin of the genetic code.

## Authoritative input

Kosiol, Holmes and Goldman (2007), *An empirical codon model for protein
sequence evolution*, Molecular Biology and Evolution 24:1464–1479.

Supplementary archive:

`https://www.ebi.ac.uk/goldman-srv/ECM/SupplMat.tar.gz`

Downloaded archive SHA256:

`44d2dab71d5eca8175e5deb4132ffd3bf1a58380ab92270ba70a28fef426d039`

Two source models are analyzed:

- `ECMrest`: restricted to single-nucleotide changes;
- `ECMunrest`: allows multiple-nucleotide changes.

The files provide 1,830 lower-triangular symmetric exchangeabilities, 61
equilibrium codon frequencies, and the codon order. Following the source
documentation, `q_ij = s_ij * pi_j`; diagonal entries make rows sum to zero;
the generator is rescaled to equilibrium mean substitution rate one; and
`P(t) = exp(tQ)`.

## Designated partitions

- `B16`: 16 boxes sharing the first two nucleotides;
- `AA20`: 20 canonical amino-acid outcomes;
- `L23`: 23 local translation blocks `(first two nucleotides, amino acid)`.

After excluding STOP codons, `L23` refines both `B16` and `AA20`, while `B16`
and `AA20` do not refine one another:

`61 -> L23 -> {B16, AA20}`.

These labels are never supplied to the recovery algorithms.

## Horizons

Analyze population transition matrices at fixed branch lengths
`t = 0.01, 0.05, 0.10, 0.50, 1.00` expected substitutions per codon site.

## Question 1: dynamic closure

For each designated partition and `P(t)`, calculate the maximum total-variation
distance between a microstate's aggregated next-macro distribution and the
uniform-within-block macro TPM. Exact closure is declared only below `1e-10`.

For each partition independently, generate 2,000 random relabelings that
preserve its exact multiset of block sizes. Report its percentile and one-sided
empirical p-value `(1 + count[random <= observed]) / 2001` for unusually low
closure residual.

No exact closure is expected. Exceptional relative performance will not be
called exact closure.

## Question 2: label-blind recovery

Two cardinality-conditioned algorithms receive only `P(t)` and requested
`k in {16, 20, 23}`:

1. `jump_profile`: k-means on the row-normalized conditional destination
   distribution given that the codon changes;
2. `spectral`: k-means on the row-normalized leading-eigenvector embedding of
   the reversible similarity transform of `P(t)`.

Use 30 independently seeded fits per model, horizon, method, and `k`. Evaluate
adjusted Rand index (ARI) against `B16`, `AA20`, and `L23` only after fitting.
This is label-blind but not blind selection of `k`.

Each fit uses k-means++ with `n_init=20`. All randomization derives from the
fixed root seed `20260915`.

## Outcome categories fixed before execution

- `EXACT`: ARI 1.0 in at least 90% of fits at a matching `k`.
- `STRONG`: median ARI at least 0.75 without meeting `EXACT`.
- `PARTIAL`: median ARI at least 0.40 without meeting `STRONG`.
- `WEAK/NULL`: median ARI below 0.40.

These are operational descriptors, not formal significance thresholds.

## Primary gate

The empirical step is considered informative if at least one designated
partition either:

1. beats at least 99% of matched-size random partitions for closure in both ECM
   models at one or more common horizons; or
2. reaches `PARTIAL` or better recovery at its matching `k` in both ECM models
   at one or more common horizons.

Failure is a result and will not trigger threshold or horizon changes.

## Claim ceiling

The ECM matrices summarize substitutions estimated from protein-coding
alignments. They are empirical evolutionary models, not organism-specific
mutation measurements, direct biochemical transition matrices, or tRNA-wobble
models. A positive result shows structural signal in these ECM models only.
