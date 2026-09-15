# Effective Information against matched partition controls

## Relation to the existing triadic controls

This control uses the same five 48-state models and the same designated
eight-family partition as the partition-recovery experiment in
[`TRIADIC_CONTROLS.md`](TRIADIC_CONTROLS.md). The earlier experiment asks
whether spectral clustering recovers a partition and whether that partition
predicts held-out transitions. This experiment asks a different question:

> Does Effective Information distinguish the designated partition from
> balanced random partitions of the same cardinality, and does that ranking
> imply exact Markovian closure?

It does not introduce a new definition of causal emergence or a general
three-part validation protocol.

## Design

For each model, the designated partition was compared with 2,000 independently
sampled balanced partitions. Every control partition contains exactly eight
groups of six microstates. We calculated:

- worst-row total-variation closure residual;
- macro-level Effective Information under uniform macro intervention;
- macro EI minus 48-state micro EI;
- rank correlation between macro EI and negative closure residual.

The five models are the exact product systems with
\(\eta\in\{0,0.02,0.2,0.8\}\) and the isospectral rotation of the
\(\eta=0.02\) system that breaks the designated partition.

The protocol was written before the first run of this audit. All reported
results were independently reconstructed from the self-contained code and
matched the committed summary exactly.

## Population results

| Model | Designated closure TV | Designated macro EI | Maximum random EI | Spearman \(\rho(\mathrm{EI},-\mathrm{closure})\) |
|---|---:|---:|---:|---:|
| \(\eta=0\) | \(3.90\times10^{-18}\) | 0.853952 | 0.143924 | -0.720 |
| \(\eta=0.02\) | \(4.51\times10^{-17}\) | 0.853952 | 0.148637 | -0.736 |
| \(\eta=0.2\) | \(1.02\times10^{-16}\) | 0.853952 | 0.180239 | -0.663 |
| \(\eta=0.8\) | \(5.94\times10^{-17}\) | 0.853952 | 0.752335 | -0.399 |
| isospectral broken | 0.005079 | 0.853893 | 0.132857 | -0.742 |

The designated partition had higher macro EI and lower closure residual than
all 2,000 balanced random controls in every model. This establishes that it is
exceptional relative to this control ensemble.

It does **not** establish exact closure. The isospectral broken partition also
outperformed every random control on both rankings while retaining a nonzero
population closure residual.

Within the random-control ensembles, EI and closure did not behave as surrogate
objectives. Higher EI was associated with larger—not smaller—closure residuals.
Therefore, maximizing EI is not a reliable shortcut for minimizing closure,
even though the designated partition is an exceptional outlier on both axes.

## Absolute EI difference

| Model | \(EI_{macro}-EI_{micro}\) (bits) |
|---|---:|
| \(\eta=0\) | approximately 0 |
| \(\eta=0.02\) | -0.001139 |
| \(\eta=0.2\) | -0.104151 |
| \(\eta=0.8\) | -1.507947 |
| isospectral broken | -0.001341 |

Consequently, a partition can be exceptional relative to matched random
partitions—and can even be exactly lumpable—without exhibiting positive causal
emergence in absolute EI.

## Finite-data temporal null

For 100, 500, and 5,000 transitions per state, 100 independent count matrices
were generated per model and condition. The null redistributes observed future
states among origin states without replacement. It preserves every row total
and the pooled future-state counts exactly while destroying their association.

Neither the observed data nor the shuffled null produced a positive
macro-minus-micro EI difference. The null positive rate was 0/100 in all 15
conditions. Its median difference approached zero from below, approximately
-0.371 bits at 100 samples per state, -0.066 at 500, and -0.0065 at 5,000.

This behavior is consistent with unequal finite-sample plug-in bias at 48-state
and eight-state resolutions. In this benchmark, the raw positive-difference
rule is conservative under the null rather than falsely positive.

## Interpretation and limits

Matched random partitions can demonstrate that a candidate mapping is
non-random relative to a specified ensemble. They cannot certify:

- exact Markovian closure;
- uniqueness of the macrovariable;
- positive causal emergence;
- biological or ontological validity;
- validity of the same null transformation for empirical domains.

The correct conclusion is narrower: EI, closure, and performance relative to
matched partition controls are distinct properties and should be reported
separately.

## Reproduce

```bash
python src/triadic_ei_matched_controls.py --output results/triadic
PYTHONPATH=src python -m unittest tests/test_triadic_ei_matched_controls.py
```

