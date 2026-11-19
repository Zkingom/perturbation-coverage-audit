# Supplementary Note 1 (draft v0.2): post-confirmation exploratory and illustrative analyses

## 1. Exploratory analysis (protocol explore_modularity_v1.0, frozen before running)

*This analysis was specified after both confirmations were unblinded (protocol frozen before running; reports/protocol/explore_modularity_v1.0.md) and does not alter any confirmatory conclusion.*

*Measurement noise does not explain the failure.* Split-half reliability of effect strength (Spearman–Brown, cells split into halves A and B) was ≈0.99 in all three datasets (GWPS 0.992, 95% CI 0.992–0.993), and correcting ρ(ĥ, E) for attenuation left it essentially unchanged (GWPS 0.110). By the pre-specified rule, low predictability on GWPS is attributed to a lack of information in the prior neighbours rather than to noise in E. Split-half reliability measures repeatability across cell halves, so components shared by both halves (for example a cell-number-dependent noise floor) count as signal; it is not a measure of validity. The GWPS distribution of log effect strength was also narrower (SD 0.34 vs 0.44 on RPE1), which may contribute to lower ρ. Removing noise shared between width allocation and coverage assessment (allocating width with half A, assessing coverage with half B) lowered oracle T5 coverage only slightly (nominal α; GWPS 0.898 → 0.849; RPE1 0.990 → 0.994), so the attainable repair remains large.

*Effect strength is network-modular among essential but not non-essential genes.* The correlation between a gene's effect strength and the STRING-weighted mean of its neighbours' effect strengths was 0.41 on RPE1 and 0.11 on GWPS (difference 0.29, 95% CI 0.23–0.35). Within every dataset, ĥ predicted effect strength better for CORUM complex members and for high-STRING-degree targets (high minus low degree tertile: RPE1 Δρ = +0.40, GWPS +0.14; 7 of 8 pre-specified contrasts significant after Holm correction; Fig. 1d). Most GWPS targets lie in sparse network regions (66% below the RPE1 lower degree tertile), but a gap remained within matched degree strata (RPE1 minus GWPS Δρ 0.11–0.30), so network composition explains only part of the difference. Descriptively, the repair concentrated in densely connected genes (GWPS high-degree, fixed cut-points: +10.2 pp) and the T5 coverage change in the lowest-degree tertile was negative in all three datasets (CIs excluding zero on RPE1 and GWPS) (e.g. RPE1 -14.9 pp). Two caveats apply. Because ĥ is built from neighbours' effect strengths, network autocorrelation and ĥ predictability overlap conceptually. And STRING degree partly reflects how well a gene has been studied.


## 2. Illustrative analysis (outside any protocol; requested after the exploratory results)

*Purpose.* Split-half reliability counts components shared by both halves as signal; a cell-number-dependent noise floor in the selection of the top-50 genes is the most obvious such component. Two checks were added: (a) a cross-split effect strength, in which the top-50 specific genes are selected on half A and the magnitude is evaluated on half B (and vice versa; mean of both), so that selection noise cannot inflate the magnitude; (b) the reliability of E_A vs E_B after partialling out log cell number (rank-based partial correlation, Spearman–Brown step-up). Pre-written reading: if ρ(ĥ, cross-split E) on GWPS is close to 0.110 the "prior lacks information" interpretation stands; if clearly higher (e.g. > 0.25) the framing must be re-discussed. 95% CIs: TSS-group cluster bootstrap, 2,000 resamples.

| | development | RPE1 | GWPS |
|---|---|---|---|
| ρ(ĥ, E_true) | 0.494 | 0.427 | 0.110 |
| ρ(ĥ, cross-split E) | 0.502 [0.464, 0.538] | 0.368 [0.327, 0.409] | 0.119 [0.093, 0.144] |
| ρ(E_true, cross-split E) | 0.936 | 0.927 | 0.900 |
| reliability (Spearman–Brown) | 0.989 | 0.987 | 0.992 |
| reliability of cross-split E | 0.991 | 0.991 | 0.972 |
| partial reliability given log cells | 0.987 [0.984, 0.988] | 0.976 [0.972, 0.979] | 0.960 [0.955, 0.964] |
| ρ(E_true, log cells) | -0.256 | -0.512 | -0.799 [-0.813, -0.784] |
| ρ(cross-split E, log cells) | -0.002 | -0.229 | -0.542 |

*Reading (pre-written rule).* GWPS ρ(ĥ, cross-split E) = 0.119, close to 0.110; the stop rule (> 0.25) was not triggered, and the "prior lacks information" interpretation stands. Partial reliability given cell number remains high (GWPS 0.960).

*An observation that must be reported.* On GWPS, measured effect strength is strongly anti-correlated with the number of cells per perturbation (ρ = -0.80; development -0.26, RPE1 -0.51). The cross-split definition removes selection inflation and weakens, but does not remove, this coupling (GWPS -0.54; development ≈ 0). Two explanations cannot be separated with these data: (i) a noise floor — with fewer cells the RMS of residual sampling noise on any gene set is larger, which the cross-split construction does not remove; (ii) biology — perturbations with larger transcriptional effects also impair proliferation and therefore yield fewer cells. Under (i), part of the "strong-effect" tier on GWPS reflects measurement precision rather than biological effect size; the coverage failure would still be real for those perturbations, but its interpretation as an effect-strength phenomenon would be weaker on GWPS than on the essential-gene screens. Because cell number is an experimental outcome, it is excluded from ĥ by design, which may contribute to low predictability on GWPS. A noise-floor-corrected effect strength (e.g. subtracting the effect strength of cell-number-matched pseudo-perturbations sampled from non-targeting control cells) would address this; it has not been computed.


## 3. Illustrative analysis: noise-floor-corrected effect strength (outside any protocol; requested 2026-10-16)

*Construction.* For each dataset, non-targeting control cells (development 10,691; RPE1 11,485; GWPS: random 20,000 of 75,328) were transformed exactly as in the pseudobulk pipeline. Pseudo-perturbations of n cells were drawn without replacement (geometric grid of 20 values of n from 30 to the largest cell number among E units; 200 draws per n). Their effect strength was computed with the E_true definition and each split's own training reference; control cells are not part of that reference. For each perturbation, E_corr = E_true − m(n) and E_z = (E_true − m(n))/s(n), where m and s are the null mean and SD at its cell number (linear interpolation in log n). Perturbations were re-tiered by E_corr (and E_z) with the same quintile rule, and M0 coverage by the new tiers was estimated with the same GEE as the confirmatory B0. Recomputed E_true and original tiers matched the stored caches exactly, and the original tiers reproduced the confirmatory B0 estimate exactly. Pre-written reading: GWPS M0 T5 coverage under corrected tiers > 0.75 → the GWPS failure is mainly attributable to measurement precision.

*Null.* The null effect strength of a no-effect perturbation falls steeply with cell number (split 0: development 5.30 at n = 30 to 0.99 at n = 1996; RPE1 4.63 to 1.11). This floor is far from negligible relative to measured effect strengths, and 9% (development) and 9% (RPE1) of perturbations do not exceed it.

| | development | RPE1 | GWPS |
|---|---|---|---|
| ρ(E_true, log cells) | -0.256 | -0.512 | -0.799 |
| ρ(E_corr, log cells) | 0.221 | -0.044 | 0.003 |
| ρ(E_true, E_corr) | 0.808 | 0.844 | 0.468 |
| T5 overlap (original vs corrected tiers) | 0.81 | 0.87 | 0.59 |
| **M0 T5 coverage, original tiers** | 0.516 [0.468, 0.564] | 0.480 [0.431, 0.528] | 0.518 [0.490, 0.546] |
| **M0 T5 coverage, E_corr tiers** | 0.552 [0.503, 0.600] | 0.492 [0.443, 0.542] | 0.664 [0.637, 0.690] |
| M0 T5 coverage, E_z tiers | 0.621 [0.572, 0.668] | 0.551 [0.501, 0.600] | 0.711 [0.685, 0.736] |
| ρ(ĥ, E_true) / ρ(ĥ, E_corr) | 0.494 / 0.513 | 0.427 / 0.367 | 0.110 / 0.127 |

*Reading (essential-gene screens).* On both essential-gene screens the coverage failure is essentially unchanged after noise-floor correction (T5 0.55 and 0.49), while the correction removes the cell-number association (development ρ −0.26 → +0.22; RPE1 −0.51 → −0.04). Standardising by the null SD (E_z) moves some coverage loss into T3–T4 and raises T5 modestly (0.62, 0.55), as expected for a score that ranks perturbations partly by precision; the failure remains. ĥ's predictability is similar for corrected effect strength. *Reading (GWPS, pre-written rule).* M0 T5 coverage under noise-corrected tiers is 0.664 [0.637, 0.690], below the 0.75 threshold, so the "effect-strength" framing is retained for GWPS and the cell-number association is reported as a limitation. The correction nonetheless matters much more on GWPS than on the essential-gene screens: the null effect strength of a no-effect pseudo-perturbation falls from 7.34 at n = 30 to 1.13 at n = 1,170 (split 0); 21% of GWPS perturbations do not exceed the null mean at their cell number (essential screens ≈ 9%); corrected and original effect strength correlate only moderately (ρ = 0.47); only 59% of the original T5 remains in the corrected T5; and T5 coverage rises by ~15 pp (to 0.71 after SD standardisation), with the coverage loss partly redistributed to weaker tiers (T1 0.956). The proxy's predictability is not rescued (ρ(ĥ, E_corr) = 0.127). We therefore read GWPS as follows: a coverage failure for the perturbations with the strongest effects beyond the noise floor persists, but a sizeable share of the uncorrected GWPS failure reflects measurement precision of low-cell-number perturbations.


## 4. Exploratory: extended K/d grids on RPE1 (protocol explore_gridK_rpe1_v1.0, frozen before running)

Candidate grids were extended to K ∈ {5, 10, 20, 30, 50, 100, 200} and d ∈ {2, 5, 10, 20, 30, 50, 100, 200}, with everything else frozen. A pre-check with the original grids reproduced the confirmatory stage-1 cache of split 0 byte for byte. Selected K: 50 (8/21 splits), 100 (8/21), 200 (5/21); selected d: 200 in 21/21 splits (still at the grid ceiling; reported as is, no further extension). Validation S1 error improved only marginally relative to the confirmatory choice (mean −0.35% for K, −1.5% for d).

| RPE1 | confirmatory (K ≤ 50, d ≤ 100) | extended grids |
|---|---|---|
| M0 T5 coverage | 0.480 [0.431, 0.528] | 0.491 [0.443, 0.539] |
| step 1 E1 (M1 − M0) | +9.7 [5.8, 13.6] pp | +9.5 [5.6, 13.4] pp |
| step 1 E2 (subset width ratio) | +17.4 [12.6, 22.3] pp (1.38) | +17.8 [12.9, 22.7] pp (1.38) |
| step 2 E1 (M1 − Mdist-M) | +10.9 [7.8, 14.0] pp | +10.7 [7.4, 13.9] pp |
| oracle E1 | +52.0 [47.1, 56.8] pp | +50.7 [45.9, 55.6] pp |
| nominal marginal coverage | 0.900 | 0.899 |

*Reading (pre-written rule).* M0 T5 coverage stays < 0.85 and M1 stays above M0: the failure and the partial repair persist under more extensive tuning of the linear models. The pre-specified classification wording would be the same as in the confirmation; this run is exploratory and does not replace it.
