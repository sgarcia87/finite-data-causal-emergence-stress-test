# Triadic controls: spectrum, partitions, and exact macro-dynamics

## Purpose

These controls test distinctions that a singular-value spectrum alone cannot
resolve:

1. whether a designated state partition is dynamically closed;
2. whether a valid partition can be recovered from transition data;
3. whether a resolved spectral dimension is a macrostate count;
4. whether exact closure implies positive causal emergence.

They supplement the finite-sample rank benchmark. They do not propose a new
definition of causal emergence.

## State space

The construction has 48 states arranged as eight unordered triadic families,
each with six orderings. Its undirected Cayley graph is

\[
Q_3\mathbin{\square}K_{3,3},
\]

with 48 vertices, 144 edges, and degree six. The family labels and the classes
of the invariant \(t\) are evaluator-only reference partitions.

## Control 1: isospectral partition breaking

The exact family uses

\[
P_\eta=Q\otimes R_\eta.
\]

For the designated eight-family partition, the models with
\(\eta\in\{0,0.02,0.2,0.8\}\) are exactly lumpable. An orthogonal conjugation
of the \(\eta=0.02\) model preserves its singular values to numerical precision
but breaks closure of that same partition.

| Comparison | Result |
|---|---:|
| maximum singular-value difference | \(2.22\times10^{-16}\) |
| maximum entry of closure residual | 0.00126979 |
| worst-row total-variation residual | 0.00507917 |

Thus an identical singular-value spectrum does not certify the validity of a
particular macro-partition.

## Control 2: full-rank exact hierarchy

For the graph Laplacian \(L\), define the heat dynamics

\[
P_\tau=\exp(-\tau L).
\]

Every singular value is strictly positive for finite \(\tau\), so the
microscopic TPM has algebraic rank 48. Integer generator identities nevertheless
prove the nested exact hierarchy

\[
48\longrightarrow 8\longrightarrow 2.
\]

The eight-state quotient is the cube \(Q_3\). The two-state quotient groups
cube vertices by parity and has generator Laplacian

\[
\begin{pmatrix}3&-3\\-3&3\end{pmatrix}.
\]

Full microscopic rank is therefore compatible with multiple exact autonomous
macrovariables.

## Control 3: resolved modes are not macrostate counts

The exact hierarchy above does not change with diffusion time or sample size.
The finite-data resolved-mode diagnostic does:

| \(\tau\) | 100 samples/state | 500 | 5,000 |
|---:|---:|---:|---:|
| 0.1 | 48 | 48 | 48 |
| 0.4 | 11 | 24 | 40 |
| 1.0 | 1 | 4 | 8 |

Values are medians across 100 independent count matrices per condition. Hence
`r_resolved` describes statistically supported spectral structure at a chosen
observation scale; it is not an estimator of the number of valid macrostates.

## Control 4: recovery, validity, and prediction separate

A spectral clustering baseline is supplied with \(k=8\) but not the family
memberships:

- at \(\eta\le0.2\), it recovers the designated families in all tested runs;
- at \(\eta=0.8\), recovery is 0% with median ARI about -0.127 although the
  designated partition is exactly lumpable;
- for the isospectral broken control, recovery is 100% although the designated
  partition is not exactly lumpable.

An additional exploratory baseline searches \(k=2,\ldots,12\). With 5,000
samples per state it selects \(k=4\) in every run. When forced to use \(k=2\),
it often recovers a different exact two-state macrovariable—a slow cube-coordinate
cut—rather than the designated parity partition.

This distinguishes:

- recovery of designated labels;
- discovery of some valid macrovariable;
- exact dynamical closure;
- held-out predictive performance.

## Effective Information

For the heat dynamics tested here, coarse-graining lowers absolute Effective
Information. Exact lumpability is therefore not, by itself, evidence of a
positive EI gain. Normalized effectiveness can increase for the eight-state
quotient, illustrating why absolute and normalized objectives must be named
separately.

## Claim ceiling

These controls do not show that:

- all spectral causal-emergence methods fail;
- lumpability is the unique definition of a valid macrostate;
- the designated triadic hierarchy is the only exact hierarchy;
- the exploratory eigengap baseline represents an optimal discovery method;
- exact closure implies causal emergence.

## Reproduce

From the repository root:

```bash
python src/triadic_partition_recovery.py --output results/triadic
python src/triadic_exact_hierarchy.py --output results/triadic
python src/triadic_blind_discovery.py --output results/triadic
PYTHONPATH=src python -m unittest discover -s tests -p 'test_triadic_*.py'
python scripts/make_triadic_figure.py
```
