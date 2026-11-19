# Methods patch (2026-10-19) — not yet merged into Methods v0.3

Four items that must be fixed regardless of the GEARS outcome (round-2 review items 4.4, 4.6, 5.2, 6.3). Each block below is written to be dropped into Methods verbatim; the anchor names the subsection it replaces or follows. Every formula was read off the frozen code (file and function named in each block); nothing here changes a computation.

---

## P-1. Replaces the sentence on cluster-robust inference (review 4.4: GEE must be independently executable)

**Interval estimation.** Let the analysis set be one (method, model, k, matching mode) combination. Each row is one *test unit*: one target gene in one split, so a gene contributes up to 21 rows (3 rounds × 7 folds, each gene a test unit once per round) and rows of the same gene are not independent. Write y_{gs} ∈ {0,1} for the coverage indicator of gene g in split s, and c(g) for the TSS group of g.

*Coverage per tier* (`analysis_v06.gee_tiers`). For each effect-strength quintile separately we fit an intercept-only GEE

  g(E[y_{gs}]) = β₀,  g = logit, marginal variance from the binomial family, working correlation = independence, clusters = TSS group c(g),

with `statsmodels.api.GEE(..., family=Binomial(), cov_struct=Independence())`. All rows of the same TSS group — including the same gene's rows across splits and different genes sharing a TSS group — are one cluster, so the sandwich (Liang–Zeger) covariance absorbs both the within-gene repetition and the within-TSS-group dependence. The reported coverage is expit(β̂₀) and the interval expit(β̂₀ ± 1.96·SE(β̂₀)) with the sandwich SE. The independence working correlation makes the point estimate numerically equal to the pooled mean, so only the interval depends on the GEE. When a tier is all-covered or none-covered (β̂₀ not finite), we report the pooled proportion with a cluster bootstrap interval instead (below) and flag the row as `method = "bootstrap (separation)"`.

*Paired coverage difference* (`analysis_v06.paired_T5`). For a method m versus the reference M0 in one tier, the two arms contain the *same* (split, gene) rows (asserted in code). We stack them, add an arm indicator m_i ∈ {0,1}, and fit the same logistic GEE with design [1, m] and the same clusters. The difference and its interval are obtained by the delta method on the two fitted probabilities p₁ = expit(β̂₀+β̂₁), p₀ = expit(β̂₀):

  Δ̂ = p₁ − p₀,  ∇ = [p₁(1−p₁) − p₀(1−p₀),  p₁(1−p₁)],  SE(Δ̂) = √(∇ᵀ V̂ ∇),

with V̂ the sandwich covariance of (β̂₀, β̂₁); the interval is Δ̂ ± 1.96·SE(Δ̂) and the two-sided p-value 2Φ(−|Δ̂|/SE). If either arm is all-covered or none-covered we refit with an identity link (`family=Gaussian()`), in which case Δ̂ = β̂₁ directly; such rows are flagged `how = "identity-link GEE (separation)"`. E1 is this quantity in the top quintile; E2 is the same quantity on the decision subset (P-2 below).

*Cluster bootstrap* (`analysis_v06.cluster_boot`). TSS groups are resampled with replacement (2,000 draws unless stated), the statistic is recomputed on each resampled frame, and the interval is the 2.5/97.5 percentiles. This is used for separated tiers, for the confounding analysis, and wherever a statistic is not a single coverage probability.

*What the interval covers.* These intervals quantify sampling variability of the coverage estimate over genes/TSS groups and splits. They are not conformal guarantees, and they are not valid for comparisons selected after inspecting the same rows; the pre-registered endpoints were fixed before the confirmation data were scored.

---

## P-2. New subsection "Calibration methods compared" (review 4.6)

Notation for one split: calibration scores S^cal = (S_1..S_{n_cal}), test units with proxies, and q̂_α(·) the split-conformal quantile, i.e. the ⌈(n+1)(1−α)⌉-th smallest value of its argument (+∞ if that index exceeds n; `methods_v06.qhat`). Prediction sets are {y : S(x,y) ≤ t_i}; the reported width is the capped radius (Methods v0.3). Every method below uses only calibration scores and proxies; no test response enters any threshold.

| name | family | proxy | threshold for test unit i |
|---|---|---|---|
| M0 | global | — | t_i = q̂_α(S^cal) |
| M1 | normalised | ĥ (predicted effect strength) | t_i = q̂_α(S^cal_j /(ĥ_j+c)) · (ĥ_i+c), c = γ·median(ĥ^oof) |
| M2 | Mondrian | ĥ | t_i = q̂_α({S^cal_j : b(ĥ_j)=b(ĥ_i)}), bins b from quantiles of ĥ^oof |
| M3 | conditional | ĥ | randomised Gibbs–Cherian–Candès set for Φ(ĥ) = span{1, z, z²} or natural cubic splines (df 3, 5) |
| Mdist-N | normalised | d̂ (similarity / training-quality difficulty) | as M1 with ĥ → d̂ |
| Mdist-M | Mondrian | d̂ | as M2 with ĥ → d̂ |
| MD-GO | Mondrian | D_GO (raw GO distance) | as M2 with ĥ → D_GO |
| MD-STR | Mondrian | D_STRING (raw STRING distance) | as M2 with ĥ → D_STRING |
| Munc-N | normalised | σ̂_boot (model's own prediction dispersion) | as M1 with ĥ → σ̂_boot |
| Munc-M | Mondrian | σ̂_boot | as M2 with ĥ → σ̂_boot |
| Mscore-N | normalised | ŝ (all ĥ features refitted to the out-of-fold score) | as M1 with ĥ → ŝ |
| Mcomb-N | normalised | σ̂_comb = exp(OLS of log out-of-fold score on [1, log ĥ, log d̂]) | as M1 with ĥ → σ̂_comb |
| Mcomb-3 | conditional | (ĥ, d̂) | M3 with Φ = span{1, z(log ĥ), z(log d̂)} |
| reference method | normalised | measured effect strength E | as M1 with ĥ → E (not deployable; see Results) |
| M_perm | negative control | permuted ĥ | any of the above with ĥ permuted independently within the calibration set and within the test set |

Implementation details that a reader needs:

1. **Bin edges and normalisation constants never use calibration or test units.** `mondrian_edges` takes quantiles of the *training out-of-fold* proxy, and c = γ·median of the same out-of-fold vector, with γ ∈ {0, 0.25, 0.5, 1} (`methods_v06.mondrian_edges`, `thresholds`).
2. **M3** (`methods_v06.m3_threshold`) implements the randomised conditional set of Gibbs, Cherian & Candès (2025) exactly, as a single linear program: with g(e) = max{Σᵢ ηᵢ Sᵢ : ηᵢ ∈ [−α, 1−α], Σᵢ ηᵢ φᵢ = −e φ_i}, the test dual satisfies {S : η_{n+1}(S) < U} = {S ≤ −g′(U)}, and −g′(U) is read from the equality-constraint marginals of that LP. U_i = u_i − α with u_i ~ Uniform(0,1) fixed per split by seed. Infeasible programs mean U lies outside the attainable range of η_{n+1} and give ±∞ (all / empty); two documented numerical fallbacks (interior-point, then bisection on η_{n+1}) compute the same quantity. The implementation was cross-checked against `conditionalconformal` 0.0.5 (protocol v0.6.1 §11c).
3. **σ̂_boot** (`ext_v07dev.compute_ext`, `_unc_from_moments`): B = 50 refits of linear-PCA on TSS-group bootstrap resamples of the training pool; per unit, the σ-scaled prediction dispersion (variance with ddof = 1) averaged over the 50 genes with the largest |mean prediction|, target column(s) excluded. Out-of-fold values for training units use only refits in which the unit was out of bag (minimum 3).
4. **D_GO, D_STRING** (`pipeline_v05.distances`): prior-only distances from the query unit to the training pool, identical to the covariates used elsewhere (asserted against the stage-1 cache at run time); out-of-fold versions for training units use the inner folds.
5. **M_perm** (`methods_v06.permute_within`, `eval_v06.N_PERM = 20`): the proxy is permuted within calibration and within test separately, 20 permutations with fixed seeds; the decision subset is always selected with the real ĥ, so the control changes only the calibration, not the subset.
6. **Width matching** (`confirm_v07.rand_rows`): for each method the level α′ is chosen on the grid 0.010–0.400 (step 0.005) so that the capped mean width equals M0's; the two bracketing grid points are mixed with a per-unit uniform draw (probability π) to hit the target in expectation. Statuses `unmatched_too_narrow` / `unmatched_too_wide` record methods whose width cannot reach M0's anywhere on the grid (this is the case for the reference method; Supplementary Table S1). After matching, the nominal level is α′, not α; matching is an **empirical equal-width comparison** and does not by itself confer a finite-sample guarantee at any level (see Methods v0.3, "Width matching").

---

## P-3. Replaces the "computed once" sentence (review 5.2)

**What was computed when.** Each confirmatory endpoint was computed once from the frozen protocol, on stage-1 caches produced before any analysis choice could be revisited, and no confirmatory configuration, criterion or classification was changed afterwards. For the non-library target set the confirmatory block was computed first (the descriptive summary over 9.2 million coverage rows exceeded the available memory at the time) and the full descriptive run followed; when it completed, every confirmatory block was **recomputed numerically and matched the first computation exactly**, with randomised-matching rows byte-identical. The same ordering and the same exact agreement hold for RPE1. "Computed once" therefore refers to the decision: a confirmatory number was never re-derived under a different configuration, although it was re-executed for verification. The full chronology, including this ordering, is Supplementary Table S3.

---

## P-4. Addition to "Evaluation units" and to Limitations (review 6.3)

**Cell-count floor applies to the reference, not to the fitted models.** Evaluation units require ≥30 cells. Inside a split, the reference quantities μ_tr and σ_g — and therefore the mean-response model, whose prediction is exactly μ_tr — are computed only from training-pool units with ≥30 cells (`pipeline_v05.MIN_CELLS_REF = 30`, asserted at run time). The linear models' fitting columns and the neighbour features of the deployable proxies (ĥ's GO- and STRING-weighted neighbour effect strength, and the neighbour cell-count/noise features) carry **no** cell-count floor, so responses measured in fewer than 30 cells do enter the gene embeddings, the fitted coefficients and the neighbour averages. The asymmetry is deliberate for the reference (which must be stable) but was not imposed on the predictors; it means that part of the proxy's input is noisier than the reference it is compared against. A sensitivity analysis weighting neighbour features by cell number, or applying the same floor to them, was not run and is stated as future work.

*Limitations addition.* Together with the exclusion of perturbations measured in fewer than 30 cells (Supplementary Table S4), this is a two-sided caveat on cell number: the excluded units have systematically larger measured effects but the least precise measurements, and among the units that are kept, low-cell-count responses still enter the predictors and the proxies while being excluded from the reference.

---

## Pointers for the four items

| item | file / function |
|---|---|
| GEE, paired difference, bootstrap | `analysis_v06.py`: `gee_tiers`, `paired_T5`, `cluster_boot`, `decision_paired` |
| method registry and thresholds | `methods_v06.py`: `METHODS`, `qhat`, `m1`, `m2`, `mondrian_edges`, `phi_basis`, `phi_lin2`, `m3_threshold`, `thresholds`, `match_width`, `capped_width`, `permute_within` |
| extension proxies | `ext_v07dev.py`: `compute_ext`, `_unc_from_moments`, `oof_distances` |
| proxy wiring per method | `eval_v06.py`: `proxies`, `apply_methods`, `N_PERM` |
| randomised matching rows | `confirm_v07.py`: `rand_rows`, `e2`, `DEC_FRAC = 0.10` |
| cell-count floor | `pipeline_v05.py`: `MIN_CELLS_REF`, `train_objects`; `proxy_v06.py`: `assemble` |
