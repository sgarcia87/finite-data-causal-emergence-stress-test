# Empirical ECM audit

## Question

Do biologically interpretable partitions of the 61 sense codons leave a
detectable signature in empirical codon-substitution dynamics that were not
constructed from those partitions?

This analysis uses the restricted and unrestricted Empirical Codon Models
(`ECMrest` and `ECMunrest`) of Kosiol, Holmes, and Goldman (2007). Their
exchangeabilities and equilibrium frequencies were estimated by maximum
likelihood from reliable coding-DNA columns in 7,332 protein families.

The complete design was fixed before analysis in
[`EMPIRICAL_ECM_PROTOCOL.md`](EMPIRICAL_ECM_PROTOCOL.md). Its pre-execution
SHA256 is retained in
[`EMPIRICAL_ECM_PROTOCOL_PREEXECUTION.sha256`](EMPIRICAL_ECM_PROTOCOL_PREEXECUTION.sha256).

## Source validation

The matrices were downloaded from the official Goldman Group supplementary
archive. The committed source files match these SHA256 values:

```text
51f97e4310a41ba3a275446f5f46d3a0dfc40da875e01f62e13b38bb76d7b141  ECMrest.dat
4a4701ec245b128c91c38bf36678ad2d6904ebe37c61fbaa6f1cdeaacc246de0  ECMunrest.dat
```

Both parsed matrices yield valid reversible continuous-time Markov generators
with equilibrium mean rate one. `ECMrest` has positive exchangeability for
14.37% of codon pairs, consistent with restricting instantaneous changes to a
single nucleotide; `ECMunrest` has 96.50% positive exchangeabilities.

## Partitions

- `B16`: 16 groups sharing the first two nucleotide positions;
- `AA20`: the 20 canonical amino-acid outcomes;
- `L23`: 23 local `(first-two-bases, amino acid)` blocks.

After excluding STOP codons, their refinement geometry is

\[
61\rightarrow L23\rightarrow\{B16,AA20\}.
\]

`L23` refines both branches; neither `B16` nor `AA20` refines the other.

## Matched-partition closure result

For each model, five transition horizons, and all three partitions, closure was
compared with 2,000 random relabelings preserving the exact group-size
multiset. Across all 30 conditions,

\[
\boxed{0/2000\text{ matched controls had equal or lower closure residual}.}
\]

The finite-control one-sided value is therefore `1/2001 = 0.00049975` in every
condition. None of the designated partitions closes exactly.

Representative results:

| Model | Horizon | Partition | Observed residual | Random median | Exact? |
|---|---:|---|---:|---:|---|
| `ECMrest` | 0.10 | `B16` | 0.07459 | 0.16342 | No |
| `ECMrest` | 0.10 | `AA20` | 0.06734 | 0.16057 | No |
| `ECMrest` | 0.10 | `L23` | 0.07790 | 0.15551 | No |
| `ECMunrest` | 0.10 | `B16` | 0.04197 | 0.08301 | No |
| `ECMunrest` | 0.10 | `AA20` | 0.01988 | 0.08318 | No |
| `ECMunrest` | 0.10 | `L23` | 0.03400 | 0.08110 | No |
| `ECMunrest` | 1.00 | `AA20` | 0.09228 | 0.36632 | No |

## Cardinality-conditioned recovery

Two label-blind algorithms receive a population TPM and a supplied value of
`k`: k-means on conditional jump profiles, and k-means on the leading
eigenvector embedding of the reversible transition operator.

Median matching ARI for the spectral method is:

| Model | `B16`, k=16 | `AA20`, k=20 | `L23`, k=23 |
|---|---:|---:|---:|
| `ECMrest` | ~0.61 | ~0.78 | ~0.65 |
| `ECMunrest` | 0.555 | **1.000** | 0.799 |

In `ECMunrest`, spectral clustering recovers `AA20` exactly in all 30 seeded
fits at every tested horizon. Yet the closure residual of that same partition
ranges from 0.00216 to 0.09228. Thus,

\[
\boxed{\text{exact recovery}\not\Rightarrow\text{exact dynamic closure}.}
\]

The jump-profile baseline is substantially weaker, showing that recoverability
also depends on the algorithm's representation of the dynamics.

## Interpretation and novelty ceiling

Kosiol et al. already reported that codon affiliation, encoded amino acid, and
amino-acid physicochemical properties are major factors in codon evolution.
This repository does not present the amino-acid signal as a new biological
discovery.

The narrower methodological observation is that, in externally estimated
empirical dynamics, a canonical biological partition can outperform all
matched random controls and be recovered exactly by a label-blind,
cardinality-conditioned algorithm while still failing exact Markov closure.

## Limitations

- ECM is an average reversible evolutionary model, not event-level molecular
  trajectories or an organism-specific mutation process.
- The training alignments were selected as protein-coding sequences using the
  universal genetic code. Biological selection therefore naturally embeds
  amino-acid structure even though the base ECM did not impose an explicit
  synonymous/nonsynonymous parameter.
- The number of clusters is supplied; `k` is not inferred blindly.
- The five horizons derive from the same generator and are robustness
  conditions, not five independent datasets.
- Random controls preserve group sizes but not every aspect of codon geometry.
- The search is not exhaustive over every partition of 61 states.

## Sources

- Kosiol C, Holmes I, Goldman N. *An empirical codon model for protein sequence
  evolution.* Molecular Biology and Evolution. 2007;24(7):1464–1479.
  https://doi.org/10.1093/molbev/msm064
- Goldman Group source page:
  https://www.ebi.ac.uk/research/goldman/empirical-codon-models/
- Official supplementary archive:
  https://www.ebi.ac.uk/goldman-srv/ECM/SupplMat.tar.gz

