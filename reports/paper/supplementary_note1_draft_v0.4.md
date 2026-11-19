# Supplementary Note 1 (draft v0.4): post-confirmation exploratory and illustrative analyses

## 1. Exploratory analysis (protocol explore_modularity_v1.0, frozen before running)

*This analysis was specified after both confirmations were unblinded (protocol frozen before running; reports/protocol/explore_modularity_v1.0.md) and does not alter any confirmatory conclusion.*

*Measurement noise does not explain the failure.* Split-half reliability of effect strength (Spearman–Brown, cells split into halves A and B) was ≈0.99 in all three datasets (non-library target set 0.992, 95% CI 0.992–0.993), and correcting ρ(ĥ, E) for attenuation left it essentially unchanged (non-library target set 0.110). By the pre-specified rule, low predictability on non-library target set is attributed to a lack of information in the prior neighbours rather than to noise in E. Split-half reliability measures repeatability across cell halves, so components shared by both halves (for example a cell-number-dependent noise floor) count as signal; it is not a measure of validity. The non-library target set distribution of log effect strength was also narrower (SD 0.34 vs 0.44 on RPE1), which may contribute to lower ρ. Removing noise shared between width allocation and coverage assessment (allocating width with half A, assessing coverage with half B) lowered reference-method T5 coverage only slightly (nominal α; non-library target set 0.898 → 0.849; RPE1 0.990 → 0.994), so the attainable repair remains large.

*Effect strength is network-modular among essential but not non-essential genes.* The correlation between a gene's effect strength and the STRING-weighted mean of its neighbours' effect strengths was 0.41 on RPE1 and 0.11 on non-library target set (difference 0.29, 95% CI 0.23–0.35). Within every dataset, ĥ predicted effect strength better for CORUM complex members and for high-STRING-degree targets (high minus low degree tertile: RPE1 Δρ = +0.40, non-library target set +0.14; 7 of 8 pre-specified contrasts significant after Holm correction; Supplementary Fig. S3a). Most non-library target set targets lie in sparse network regions (66% below the RPE1 lower degree tertile), but a gap remained within matched degree strata (RPE1 minus non-library target set Δρ 0.11–0.30), so network composition explains only part of the difference. Descriptively, the repair concentrated in densely connected genes (non-library target set high-degree, fixed cut-points: +10.2 pp) and the T5 coverage change in the lowest-degree tertile was negative in all three datasets (CIs excluding zero on RPE1 and non-library target set) (e.g. RPE1 -14.9 pp). Two caveats apply. Because ĥ is built from neighbours' effect strengths, network autocorrelation and ĥ predictability overlap conceptually. And STRING degree partly reflects how well a gene has been studied.


## 2. Illustrative analysis (outside any protocol; requested after the exploratory results)

*Purpose.* Split-half reliability counts components shared by both halves as signal; a cell-number-dependent noise floor in the selection of the top-50 genes is the most obvious such component. Two checks were added: (a) a cross-split effect strength, in which the top-50 specific genes are selected on half A and the magnitude is evaluated on half B (and vice versa; mean of both), so that selection noise cannot inflate the magnitude; (b) the reliability of E_A vs E_B after partialling out log cell number (rank-based partial correlation, Spearman–Brown step-up). Pre-written reading: if ρ(ĥ, cross-split E) on non-library target set is close to 0.110 the "prior lacks information" interpretation stands; if clearly higher (e.g. > 0.25) the framing must be re-discussed. 95% CIs: TSS-group cluster bootstrap, 2,000 resamples.

| | development | RPE1 | non-library target set |
|---|---|---|---|
| ρ(ĥ, E_true) | 0.494 | 0.427 | 0.110 |
| ρ(ĥ, cross-split E) | 0.502 [0.464, 0.538] | 0.368 [0.327, 0.409] | 0.119 [0.093, 0.144] |
| ρ(E_true, cross-split E) | 0.936 | 0.927 | 0.900 |
| reliability (Spearman–Brown) | 0.989 | 0.987 | 0.992 |
| reliability of cross-split E | 0.991 | 0.991 | 0.972 |
| partial reliability given log cells | 0.987 [0.984, 0.988] | 0.976 [0.972, 0.979] | 0.960 [0.955, 0.964] |
| ρ(E_true, log cells) | -0.256 | -0.512 | -0.799 [-0.813, -0.784] |
| ρ(cross-split E, log cells) | -0.002 | -0.229 | -0.542 |

*Reading (pre-written rule).* non-library target set ρ(ĥ, cross-split E) = 0.119, close to 0.110; the stop rule (> 0.25) was not triggered, and the "prior lacks information" interpretation stands. Partial reliability given cell number remains high (non-library target set 0.960).

*An observation that must be reported.* On non-library target set, measured effect strength is strongly anti-correlated with the number of cells per perturbation (ρ = -0.80; development -0.26, RPE1 -0.51). The cross-split definition removes selection inflation and weakens, but does not remove, this coupling (non-library target set -0.54; development ≈ 0). Two explanations cannot be separated with these data: (i) a noise floor — with fewer cells the RMS of residual sampling noise on any gene set is larger, which the cross-split construction does not remove; (ii) biology — perturbations with larger transcriptional effects also impair proliferation and therefore yield fewer cells. Under (i), part of the "strong-effect" tier on non-library target set reflects measurement precision rather than biological effect size; the coverage failure would still be real for those perturbations, but its interpretation as an effect-strength phenomenon would be weaker on non-library target set than on the essential-gene screens. Because cell number is an experimental outcome, it is excluded from ĥ by design, which may contribute to low predictability on non-library target set. A noise-floor-corrected effect strength (e.g. subtracting the effect strength of cell-number-matched pseudo-perturbations sampled from non-targeting control cells) would address this; it has not been computed.


## 3. Illustrative analysis: noise-floor-corrected effect strength (outside any protocol; requested 2026-10-16)

*Construction.* For each dataset, non-targeting control cells (development 10,691; RPE1 11,485; non-library target set: random 20,000 of 75,328) were transformed exactly as in the pseudobulk pipeline. Pseudo-perturbations of n cells were drawn without replacement (geometric grid of 20 values of n from 30 to the largest cell number among E units; 200 draws per n). Their effect strength was computed with the E_true definition and each split's own training reference; control cells are not part of that reference. For each perturbation, E_corr = E_true − m(n) and E_z = (E_true − m(n))/s(n), where m and s are the null mean and SD at its cell number (linear interpolation in log n). Perturbations were re-tiered by E_corr (and E_z) with the same quintile rule, and M0 coverage by the new tiers was estimated with the same GEE as the confirmatory B0. Recomputed E_true and original tiers matched the stored caches exactly, and the original tiers reproduced the confirmatory B0 estimate exactly. Pre-written reading: non-library target set M0 T5 coverage under corrected tiers > 0.75 → the non-library target set failure is mainly attributable to measurement precision.

*Null.* The null effect strength of a no-effect perturbation falls steeply with cell number (split 0: development 5.30 at n = 30 to 0.99 at n = 1996; RPE1 4.63 to 1.11). This floor is far from negligible relative to measured effect strengths, and 9% (development) and 9% (RPE1) of perturbations do not exceed it.

| | development | RPE1 | non-library target set |
|---|---|---|---|
| ρ(E_true, log cells) | -0.256 | -0.512 | -0.799 |
| ρ(E_corr, log cells) | 0.221 | -0.044 | 0.003 |
| ρ(E_true, E_corr) | 0.808 | 0.844 | 0.468 |
| T5 overlap (original vs corrected tiers) | 0.81 | 0.87 | 0.59 |
| **M0 T5 coverage, original tiers** | 0.516 [0.468, 0.564] | 0.480 [0.431, 0.528] | 0.518 [0.490, 0.546] |
| **M0 T5 coverage, E_corr tiers** | 0.552 [0.503, 0.600] | 0.492 [0.443, 0.542] | 0.664 [0.637, 0.690] |
| M0 T5 coverage, E_z tiers | 0.621 [0.572, 0.668] | 0.551 [0.501, 0.600] | 0.711 [0.685, 0.736] |
| ρ(ĥ, E_true) / ρ(ĥ, E_corr) | 0.494 / 0.513 | 0.427 / 0.367 | 0.110 / 0.127 |

*Reading (essential-gene screens).* On both essential-gene screens the coverage failure is essentially unchanged after noise-floor correction (T5 0.55 and 0.49), while the correction removes the cell-number association (development ρ −0.26 → +0.22; RPE1 −0.51 → −0.04). Standardising by the null SD (E_z) moves some coverage loss into T3–T4 and raises T5 modestly (0.62, 0.55), as expected for a score that ranks perturbations partly by precision; the failure remains. ĥ's predictability is similar for corrected effect strength. *Reading (non-library target set, pre-written rule).* M0 T5 coverage under noise-corrected tiers is 0.664 [0.637, 0.690], below the 0.75 threshold, so the "effect-strength" framing is retained for non-library target set and the cell-number association is reported as a limitation. The correction nonetheless matters much more on non-library target set than on the essential-gene screens: the null effect strength of a no-effect pseudo-perturbation falls from 7.34 at n = 30 to 1.13 at n = 1,170 (split 0); 21% of non-library target set perturbations do not exceed the null mean at their cell number (essential screens ≈ 9%); corrected and original effect strength correlate only moderately (ρ = 0.47); only 59% of the original T5 remains in the corrected T5; and T5 coverage rises by ~15 pp (to 0.71 after SD standardisation), with the coverage loss partly redistributed to weaker tiers (T1 0.956). The proxy's predictability is not rescued (ρ(ĥ, E_corr) = 0.127). We therefore read non-library target set as follows: a coverage failure for the perturbations with the strongest effects beyond the noise floor persists, but a sizeable share of the uncorrected non-library target set failure reflects measurement precision of low-cell-number perturbations.


## 4. Exploratory: extended K/d grids on RPE1 (protocol explore_gridK_rpe1_v1.0, frozen before running)

Candidate grids were extended to K ∈ {5, 10, 20, 30, 50, 100, 200} and d ∈ {2, 5, 10, 20, 30, 50, 100, 200}, with everything else frozen. A pre-check with the original grids reproduced the confirmatory stage-1 cache of split 0 byte for byte. Selected K: 50 (8/21 splits), 100 (8/21), 200 (5/21); selected d: 200 in 21/21 splits (still at the grid ceiling; reported as is, no further extension). Validation S1 error improved only marginally relative to the confirmatory choice (mean −0.35% for K, −1.5% for d).

| RPE1 | confirmatory (K ≤ 50, d ≤ 100) | extended grids |
|---|---|---|
| M0 T5 coverage | 0.480 [0.431, 0.528] | 0.491 [0.443, 0.539] |
| step 1 E1 (M1 − M0) | +9.7 [5.8, 13.6] pp | +9.5 [5.6, 13.4] pp |
| step 1 E2 (subset width ratio) | +17.4 [12.6, 22.3] pp (1.38) | +17.8 [12.9, 22.7] pp (1.38) |
| step 2 E1 (M1 − Mdist-M) | +10.9 [7.8, 14.0] pp | +10.7 [7.4, 13.9] pp |
| reference-method E1 | +52.0 [47.1, 56.8] pp | +50.7 [45.9, 55.6] pp |
| nominal marginal coverage | 0.900 | 0.899 |

*Reading (pre-written rule).* M0 T5 coverage stays < 0.85 and M1 stays above M0: the failure and the partial repair persist under more extensive tuning of the linear models. The pre-specified classification wording would be the same as in the confirmation; this run is exploratory and does not replace it.

## 5. The conformity score is close to a multiple of effect strength (descriptive; existing outputs only)

*Purpose.* The coverage failure is heterogeneity in the scale of the conformity score. This section quantifies how much of the score is the perturbation's measured effect strength and how much is model-specific residual, using the stored coverage rows of the confirmatory and development runs (nominal α, k = 50, M0, no new model fit).

*Quantities.* For every test row, the ratio score / effect strength, its coefficient of variation (CV) across rows, the Spearman correlation between score and effect strength, the coefficient of determination of log score on log effect strength, and the relative error reduction 1 − mean(score)/mean(effect strength), which is the share of the mean-model error that the fitted model removes.

| dataset | model | CV of ratio | median ratio [5–95%] | Spearman | R² (log) | relative error reduction |
|---|---|---|---|---|---|---|
| development | mean response | 0 (by construction) | 1.000 | 1.000 | 1.000 | 0 |
| development | linear-PCA | 0.107 | 0.994 [0.81–1.04] | 0.941 | 0.888 | 4.3% |
| development | linear-GO | 0.148 | 0.983 [0.61–1.04] | 0.888 | 0.807 | 8.5% |
| RPE1 | linear-PCA | 0.097 | 0.991 [0.81–1.07] | 0.963 | 0.935 | 4.3% |
| RPE1 | linear-GO | 0.116 | 0.991 [0.71–1.06] | 0.943 | 0.906 | 6.5% |
| non-library target set | linear-PCA | 0.008 | 1.000 [0.99–1.01] | 1.000 | 0.999 | 0.1% |
| non-library target set | linear-GO | 0.014 | 1.000 [0.99–1.01] | 0.999 | 0.998 | 0.2% |

Within the top tier the reduction is larger than overall (7.6%, 7.9% and 0.3% for linear-PCA), and in the weakest tier it is around zero or negative, i.e. for weakly perturbed genes the fitted model is no better than predicting the training mean.

*Reading.* The predictors used here capture little perturbation-specific signal: they remove 4–8% of the mean-model error on the development set, 4–6% on RPE1 and 0.1–0.2% on the non-library target set, so the conformity score is largely determined by the actual effect scale of each perturbation (Supplementary Fig. S2). Two consequences must be stated plainly. First, the argument that heterogeneous coverage is heterogeneous scale is strongest where a residual remains (CV ≈ 0.10–0.15 on the two essential-gene screens) and is **close to definitional on the non-library target set**, where score and effect strength are nearly identical (CV 0.008, Spearman 1.000); conclusions about that dataset therefore rest on the pre-registered coverage comparison, not on this relation. Second, the relation is a property of these predictors, not of perturbation prediction in general: a model that captured more perturbation-specific signal would lower the ratio and could change the picture. The planned GEARS robustness run uses the same splits and score and is the test of whether a stronger, non-linear model changes this proportionality.

## 6. Coverage on prospectively identifiable subsets, and a hypothesis about network position (post-confirmation, descriptive and exploratory)

*Why this section exists.* Tiers of measured effect strength are defined by the outcome, and under a score that is a monotone function of effect strength their coverage is fixed by the quantile construction (Results; Supplementary Table S8). Subsets defined by information available **before** the experiment are not fixed in that way, so they are reported separately. Split conformal (M0), linear-PCA, k = 50, nominal α = 0.1, 21 grouped splits; 95% intervals from a cluster bootstrap over TSS groups. "Expected from tier composition" is the subset's fraction of top-tier units times the observed top-tier coverage plus the remainder times one.

| dataset | subset | coverage [95% CI] | top-tier fraction | expected from composition | difference |
|---|---|---|---|---|---|
| development | top 10% by predicted effect strength | 0.740 [0.684, 0.795] | 0.60 | 0.709 | +0.031 |
| RPE1 | top 10% by predicted effect strength | 0.708 [0.643, 0.773] | 0.51 | 0.734 | -0.026 |
| non-library target set | top 10% by predicted effect strength | 0.822 [0.795, 0.849] | 0.34 | 0.835 | -0.012 |
| development | top 20% by STRING degree | 0.871 [0.842, 0.904] | 0.39 | 0.812 | +0.059 |
| RPE1 | top 20% by STRING degree | 0.869 [0.838, 0.903] | 0.31 | 0.841 | +0.028 |
| non-library target set | top 20% by STRING degree | 0.861 [0.841, 0.881] | 0.25 | 0.880 | -0.019 |
| development | CORUM complex member | 0.886 [0.867, 0.902] | 0.23 | 0.888 | -0.002 |
| RPE1 | CORUM complex member | 0.890 [0.872, 0.909] | 0.22 | 0.883 | +0.007 |
| non-library target set | CORUM complex member | 0.892 [0.877, 0.905] | 0.23 | 0.891 | +0.001 |

Full tables: Supplementary Tables S7 (prospective subsets) and S9 (CORUM; exploratory).

*Reading.* Subsets a practitioner can name in advance do under-cover — a 10% selection by predicted effect strength is covered 0.74, 0.71 and 0.82 against a nominal 0.9, with intervals excluding 0.9 — and this is the operationally relevant statement. It is, however, not independent evidence of a calibration defect beyond the construction: the deficits agree with what the subsets' top-tier composition alone predicts, to within 3.1 percentage points for the proxy-defined subsets and 0.7 for CORUM membership. Complex membership therefore adds nothing once effect strength is accounted for.

*Hypothesis (not tested here).* One subset behaves differently and in the opposite direction: genes with many STRING partners are covered **better** than their effect-strength composition predicts — by +5.9 pp (development) and +2.8 pp (RPE1) in the top degree quintile, and +11.7 pp and +8.6 pp in the top decile, while the non-library target set shows no such excess (-1.9 pp). A possible reading is that for densely connected essential genes the training set contains many functionally similar perturbations, so the fitted model removes relatively more of the error at a given effect strength — consistent with the higher network autocorrelation of effect strength on the essential-gene screens (§1) and with the proxy's higher predictability there. This is a hypothesis generated after the confirmations were unblinded, on subsets chosen by us, with degree defined within each screen's own STRING sub-network (score ≥ 0.7); it is not a pre-registered endpoint, it was not corrected for the number of subsets examined, and it would need its own protocol and an independent screen to test.
