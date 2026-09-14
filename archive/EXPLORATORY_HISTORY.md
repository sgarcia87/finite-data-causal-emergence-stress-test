# Exploratory history

This file exists for transparency. These experiments motivated the final benchmark but should **not** be cited as successful methods.

## Stress 0.1–0.2 — state duplication

Observed that state-space refinement can change EI or rank-based CE depending on the refinement. Equal cloning was retained because it preserves EI exactly and gives the cleanest representation-sensitivity control.

## Stress 0.3 — noise stability

Showed that exact `clear CE` is discontinuous under generic perturbations that restore full rank. `Vague CE` can repair this when a clear spectral gap exists.

## Stress 0.4 — epsilon selectors

Compared several automatic cutoff heuristics. No single simple selector dominated across noise regimes. This was retained as motivation, not as a final estimator comparison claim.

## Stress 0.5 — finite-sample TPM estimation

Showed strong sample-size dependence and motivated moving away from exact empirical rank.

## Stress 0.6 — fixed candidate rank with abstention

A cross-split procedure performed well when the target rank was supplied. This was useful but not a blind solution.

## Stress 0.7 / 0.7b — blind rank discovery

The first scoring used algebraic rank even though a nonzero threshold was estimating a practical thresholded rank. The truth definition was corrected in 0.7b.

This is an example of why the repository emphasizes explicit estimands.

## Stress 0.8 — automatic noise-floor threshold

The automatic threshold controlled strong full-rank false compression well but under-resolved weak true modes. It was rejected as an exact rank estimator.

## Stress 0.9 — attempted rank interval

An attempted \([r_{\min},r_{\max}]\) interval failed because the proposed upper bound excluded weak true modes too often.

The upper-bound claim was discarded.

## Final benchmark

Only the following idea survived:

> Report positive evidence for **resolved modes** without converting unresolved modes into evidence of absence.

That is the claim retained in v1.0.
