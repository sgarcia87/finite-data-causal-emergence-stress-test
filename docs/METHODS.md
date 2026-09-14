# Methods

## 1. Transition probability matrices

All benchmark systems are row-stochastic Markov transition matrices \(P\).

The final benchmark uses an observed state-space size of:

\[
N=16.
\]

## 2. Exact behavior-preserving refinement

A lower-dimensional TPM \(P\) with \(r\) states is expanded into \(N\) observed states by replacing each original state with an equal number \(m=N/r\) of clones.

If observed clone \(a\) belongs to parent state \(i\), and clone \(b\) belongs to parent state \(j\),

\[
Q_{ab}=\frac{P_{ij}}{m}.
\]

All clones of a parent therefore have identical outgoing behavior at the parent-state level.

The expanded TPM is exactly lumpable back to \(P\).

The experiment records:

- lumping error;
- Effective Information before and after refinement;
- numerical rank;
- clear-CE status under the criterion \(rank(P)<N\).

## 3. Finite transition sampling

For every row \(i\) of a population TPM \(P\), the benchmark draws:

\[
C_i \sim \mathrm{Multinomial}(m,P_i),
\]

where \(m\) is the number of transition samples per row.

The empirical TPM is:

\[
\hat P_{ij}=\frac{C_{ij}}{\sum_j C_{ij}}.
\]

This row-stratified design isolates transition-estimation noise from the separate issue of unvisited states in a single finite trajectory.

## 4. Low-rank ground-truth families

Random full-rank parent TPMs of size 2, 4, and 8 are generated using Dirichlet rows and then equally cloned into the 16-state observed space.

Population algebraic ranks are therefore known.

## 5. Full-rank controls

Two control families are used.

### 5.1 Well-conditioned full-rank

\[
P=(1-\alpha)I+\alpha J/N
\]

with \(\alpha\in\{0.2,0.4\}\).

These matrices are full rank and have no intended low singular-value tail.

### 5.2 Near-low-rank but mathematically full-rank

An exact cloned rank-4 TPM \(Q\) is mixed with a generic full-rank TPM \(R\):

\[
P_\delta=(1-\delta)Q+\delta R
\]

for:

\[
\delta\in\{0.001,0.01,0.05\}.
\]

These systems test the distinction between **algebraic rank** and **currently resolvable modes**.

## 6. Clear causal-emergence stress test

Following the rank criterion in Zhang et al. (2025), the benchmark records whether:

\[
rank(\hat P)<N.
\]

Because finite empirical perturbations generically restore full rank, this test quantifies the gap between an exact population definition and finite-data inference.

## 7. Forced spectral baselines

Three deliberately simple baselines are compared:

1. **Largest log gap** in the singular-value spectrum.
2. **Linear scree elbow** using maximum distance to the endpoint chord.
3. **Generic square-matrix SVHT baseline**, using \(2.858\times\mathrm{median}(\sigma)\).

These methods are stress baselines. They are not claimed to be the official inference procedure of Zhang et al. or TPM-specific optimal methods.

## 8. Cross-split resolved-mode diagnostic

Finite counts are repeatedly split binomially into two independent halves \(A\) and \(B\).

After row normalization:

\[
D=\frac{\hat P_A-\hat P_B}{2}
\]

and

\[
G=
\frac{
\hat P_A^\top\hat P_B+
\hat P_B^\top\hat P_A
}{2}.
\]

Define the split discrepancy:

\[
H=D^\top D
\]

and noise scale:

\[
n=\lambda_{\max}(H).
\]

For ordered eigenvalues \(\lambda_i(G)\), the normalized margin of mode \(i\) is:

\[
m_i=
\frac{\lambda_i(G)-n}{\lambda_1(G)}.
\]

Across repeated count splits, a leading mode is called **resolved** only when the 5th percentile of its margin remains positive.

The reported quantity is the largest leading block satisfying that criterion:

\[
r_{\mathrm{resolved}}.
\]

Interpretation:

\[
r_{\mathrm{true}}\ge r_{\mathrm{resolved}}
\]

is treated as a benchmarked one-sided statement only.

The method does **not** claim:

\[
r_{\mathrm{true}}=r_{\mathrm{resolved}}.
\]

## 9. Reproducibility

Random seeds are fixed in the source code.

Reference outputs are committed under `results/`.
