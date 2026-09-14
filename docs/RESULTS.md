# Results

## A. Exact representation refinement

Across 150 random 4-state full-rank chains:

- exact lumping error after equal 4x cloning: **0**
- median absolute EI change: **2.64e-16 bits**
- clear CE before cloning: **0%**
- clear CE after cloning: **100%**

This is the cleanest representation-sensitivity result in the repository.

## B. Finite-data exact-rank recovery

At 5,000 samples per TPM row, aggregated exact-rank recovery over true ranks 2, 4, and 8 was:

| Procedure | Exact-rank recovery |
|---|---:|
| empirical clear CE | 0% |
| largest log gap | 30.7% |
| linear scree elbow | 26.9% |
| generic SVHT | 55.6% |

The methods differ substantially by true rank:

| True rank | log gap | elbow | SVHT | median resolved fraction |
|---:|---:|---:|---:|---:|
| 2 | 56.0% | 22.0% | 92.7% | 100% |
| 4 | 27.3% | 56.7% | 74.0% | 100% |
| 8 | 8.7% | 2.0% | 0% | 87.5% |

## C. Full-rank controls

The benchmark includes both well-conditioned full-rank systems and systems that are mathematically full rank but extremely close to rank 4.

The forced cutoff baselines select an internal lower rank throughout these controls in this benchmark.

By contrast, the conservative diagnostic is interpreted only as the number of currently **resolved** modes:

- well-conditioned full-rank controls resolve all 16 modes at high support;
- near-low-rank full-rank controls typically resolve about 4 strong modes at the tested support.

The latter is not called an exact low-rank result.

## D. Publication endpoints

See:

`results/final_benchmark_publication_endpoints.csv`

for the compact 100 / 500 / 5,000 samples-per-row summary.

## Central interpretation

The robust conclusion is not that a particular cutoff selector is best.

It is that finite data create an inferential asymmetry:

- a strong mode can accumulate positive evidence of being resolved;
- a weak mode being unresolved is not equivalent to evidence that the mode is absent.

This distinction motivates one-sided reporting of resolved modes unless an explicit negligibility criterion for the unresolved spectral tail is supplied.
