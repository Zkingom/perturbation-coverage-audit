# Supplementary Note 1 (draft v0.1): post-confirmation exploratory and illustrative analyses

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
