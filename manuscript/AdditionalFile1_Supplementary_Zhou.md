---
title: "Supplementary information"
subtitle: "Additional file 1 to: A diagnostic for when effect-stratified coverage of single-cell perturbation predictions is informative"
author: "Xu Zhou, Independent Researcher (aurel1501@outlook.com)"
date: "8 November 2026"
---

# Supplementary Note 1

## 1. Exploratory analysis (protocol explore_modularity_v1.0, frozen before running)

*This analysis was specified after both confirmations were unblinded (protocol frozen before running; reports/protocol/explore_modularity_v1.0.md) and does not alter any confirmatory conclusion.*

*Measurement noise does not explain the failure.* Split-half reliability of effect strength (Spearman–Brown, cells split into halves A and B) was ≈0.99 in all three datasets (non-library target set 0.992, 95% CI 0.99 to 0.99), and correcting ρ(ĥ, E) for attenuation left it essentially unchanged (non-library target set 0.110). By the pre-specified rule, low predictability on the non-library target set is attributed to a lack of information in the prior neighbours rather than to noise in E. Split-half reliability measures repeatability across cell halves, so components shared by both halves (for example a cell-number-dependent noise floor) count as signal; it is not a measure of validity. The non-library target set distribution of log effect strength was also narrower (SD 0.34 vs 0.44 on RPE1), which may contribute to lower ρ. Removing noise shared between width allocation and coverage assessment (allocating width with half A, assessing coverage with half B) lowered reference-method T5 coverage only slightly (nominal α; non-library target set 0.898 → 0.849; RPE1 0.990 → 0.994), so the attainable repair remains large.

*Effect strength is network-modular among essential-library targets but not among non-library targets.* The correlation between a gene's effect strength and the STRING-weighted mean of its neighbours' effect strengths was 0.41 on RPE1 and 0.11 on the non-library target set (difference 0.29, 95% CI 0.23 to 0.35). Within every dataset, ĥ predicted effect strength better for CORUM complex members and for high-STRING-degree targets (high minus low degree tertile: RPE1 Δρ = +0.40, non-library target set +0.14; 7 of 8 pre-specified contrasts significant after Holm correction; Supplementary Fig. S3a). Most non-library target set targets lie in sparse network regions (66% below the RPE1 lower degree tertile), but a gap remained within matched degree strata (RPE1 minus non-library target set Δρ 0.11–0.30), so network composition explains only part of the difference. Descriptively, the repair concentrated in densely connected genes (non-library target set high-degree, fixed cut-points: +10.2 pp) and the T5 coverage change in the lowest-degree tertile was negative in all three datasets (CIs excluding zero on RPE1 and non-library target set) (e.g. RPE1 −14.9 pp). Two caveats apply. Because ĥ is built from neighbours' effect strengths, network autocorrelation and ĥ predictability overlap conceptually. And STRING degree partly reflects how well a gene has been studied.


## 2. Illustrative analysis (outside any protocol; requested after the exploratory results)

*Purpose.* Split-half reliability counts components shared by both halves as signal; a cell-number-dependent noise floor in the selection of the top-50 genes is the most obvious such component. Two checks were added: (a) a cross-split effect strength, in which the top-50 specific genes are selected on half A and the magnitude is evaluated on half B (and vice versa; mean of both), so that selection noise cannot inflate the magnitude; (b) the reliability of E_A vs E_B after partialling out log cell number (rank-based partial correlation, Spearman–Brown step-up). Pre-written reading: if ρ(ĥ, cross-split E) on the non-library target set is close to 0.110 the "prior lacks information" interpretation stands; if clearly higher (e.g. > 0.25) the framing must be re-discussed. 95% CIs: TSS-group cluster bootstrap, 2,000 resamples.

| | development | RPE1 | non-library target set |
|---|---|---|---|
| ρ(ĥ, E_true) | 0.494 | 0.427 | 0.110 |
| ρ(ĥ, cross-split E) | 0.502 [0.46, 0.54] | 0.368 [0.33, 0.41] | 0.119 [0.09, 0.14] |
| ρ(E_true, cross-split E) | 0.936 | 0.927 | 0.900 |
| reliability (Spearman–Brown) | 0.989 | 0.987 | 0.992 |
| reliability of cross-split E | 0.991 | 0.991 | 0.972 |
| partial reliability given log cells | 0.987 [0.98, 0.99] | 0.976 [0.97, 0.98] | 0.960 [0.95, 0.96] |
| ρ(E_true, log cells) | −0.256 | −0.512 | −0.799 [−0.81, −0.78] |
| ρ(cross-split E, log cells) | −0.002 | −0.229 | −0.542 |

*Reading (pre-written rule).* non-library target set ρ(ĥ, cross-split E) = 0.119, close to 0.110; the stop rule (> 0.25) was not triggered, and the "prior lacks information" interpretation stands. Partial reliability given cell number remains high (non-library target set 0.960).

*An observation that must be reported.* On non-library target set, measured effect strength is strongly anti-correlated with the number of cells per perturbation (ρ = −0.80; development −0.26, RPE1 −0.51). The cross-split definition removes selection inflation and weakens, but does not remove, this coupling (non-library target set −0.54; development ≈ 0). Two explanations cannot be separated with these data: (i) a noise floor — with fewer cells the RMS of residual sampling noise on any gene set is larger, which the cross-split construction does not remove; (ii) biology — perturbations with larger transcriptional effects also impair proliferation and therefore yield fewer cells. Under (i), part of the "strong-effect" tier on the non-library target set reflects measurement precision rather than biological effect size; the coverage failure would still be real for those perturbations, but its interpretation as an effect-strength phenomenon would be weaker on the non-library target set than on the essential-gene screens. Because cell number is an experimental outcome, it is excluded from ĥ by design, which may contribute to low predictability on the non-library target set. A noise-floor-corrected effect strength (e.g. subtracting the effect strength of cell-number-matched pseudo-perturbations sampled from non-targeting control cells) would address this; this correction was computed subsequently and is reported in S3.


## 3. Illustrative analysis: noise-floor-corrected effect strength (outside any protocol; requested 2026-10-16)

*Construction.* For each dataset, non-targeting control cells (development 10,691; RPE1 11,485; non-library target set: random 20,000 of 75,328) were transformed exactly as in the pseudobulk pipeline. Pseudo-perturbations of n cells were drawn without replacement (geometric grid of 20 values of n from 30 to the largest cell number among E units; 200 draws per n). Their effect strength was computed with the E_true definition and each split's own training reference; control cells are not part of that reference. For each perturbation, E_corr = E_true − m(n) and E_z = (E_true − m(n))/s(n), where m and s are the null mean and SD at its cell number (linear interpolation in log n). Perturbations were re-tiered by E_corr (and E_z) with the same quintile rule, and M0 coverage by the new tiers was estimated with the same GEE as the confirmatory B0. Recomputed E_true and original tiers matched the stored caches exactly, and the original tiers reproduced the confirmatory B0 estimate exactly. Pre-written reading: non-library target set M0 T5 coverage under corrected tiers > 0.75 → the non-library target set failure is mainly attributable to measurement precision.

*Null.* The null effect strength of a no-effect perturbation falls steeply with cell number (split 0: development 5.30 at n = 30 to 0.99 at n = 1996; RPE1 4.63 to 1.11). This floor is far from negligible relative to measured effect strengths, and 9% (development) and 9% (RPE1) of perturbations do not exceed it.

| | development | RPE1 | non-library target set |
|---|---|---|---|
| ρ(E_true, log cells) | −0.256 | −0.512 | −0.799 |
| ρ(E_corr, log cells) | 0.221 | −0.044 | 0.003 |
| ρ(E_true, E_corr) | 0.808 | 0.844 | 0.468 |
| T5 overlap (original vs corrected tiers) | 0.810 | 0.866 | 0.586 |
| **M0 T5 coverage, original tiers** | 0.516 [0.47, 0.56] | 0.480 [0.43, 0.53] | 0.518 [0.49, 0.55] |
| **M0 T5 coverage, E_corr tiers** | 0.552 [0.50, 0.60] | 0.492 [0.44, 0.54] | 0.664 [0.64, 0.69] |
| M0 T5 coverage, E_z tiers | 0.621 [0.57, 0.67] | 0.551 [0.50, 0.60] | 0.711 [0.68, 0.74] |
| ρ(ĥ, E_true) / ρ(ĥ, E_corr) | 0.494 / 0.513 | 0.427 / 0.367 | 0.110 / 0.127 |

*Reading (essential-gene screens).* On both essential-gene screens the coverage failure is essentially unchanged after noise-floor correction (T5 0.552 and 0.492), while the correction removes the cell-number association (development ρ −0.26 → +0.22; RPE1 −0.51 → −0.04). Standardising by the null SD (E_z) moves some coverage loss into T3–T4 and raises T5 modestly (0.621, 0.551), as expected for a score that ranks perturbations partly by precision; the failure remains. ĥ's predictability is similar for corrected effect strength. *Reading (non-library target set, pre-written rule).* M0 T5 coverage under noise-corrected tiers is 0.664 (95% CI 0.64 to 0.69), below the 0.75 threshold, so the "effect-strength" framing is retained for non-library target set and the cell-number association is reported as a limitation. The correction nonetheless matters much more on the non-library target set than on the essential-gene screens: the null effect strength of a no-effect pseudo-perturbation falls from 7.34 at n = 30 to 1.13 at n = 1,170 (split 0); 21% of non-library target set perturbations do not exceed the null mean at their cell number (essential screens ≈ 9%); corrected and original effect strength correlate only moderately (ρ = 0.47); only 59% of the original T5 remains in the corrected T5; and T5 coverage rises by 14.6 pp (to 0.711 after SD standardisation), with the coverage loss partly redistributed to weaker tiers (T1 0.956). The proxy's predictability is not rescued (ρ(ĥ, E_corr) = 0.127). We therefore read non-library target set as follows: a coverage failure for the perturbations with the strongest effects beyond the noise floor persists, but a sizeable share of the uncorrected non-library target set failure reflects measurement precision of low-cell-number perturbations.


## 4. Exploratory: extended K/d grids on RPE1 (protocol explore_gridK_rpe1_v1.0, frozen before running)

Candidate grids were extended to K ∈ {5, 10, 20, 30, 50, 100, 200} and d ∈ {2, 5, 10, 20, 30, 50, 100, 200}, with everything else frozen. A pre-check with the original grids reproduced the confirmatory stage-1 cache of split 0 byte for byte. Selected K: 50 (8/21 splits), 100 (8/21), 200 (5/21); selected d: 200 in 21/21 splits (still at the grid ceiling; reported as is, no further extension). Validation S1 error improved only marginally relative to the confirmatory choice (mean −0.35% for K, −1.5% for d).

| RPE1 | confirmatory (K ≤ 50, d ≤ 100) | extended grids |
|---|---|---|
| M0 T5 coverage | 0.480 [0.43, 0.53] | 0.491 [0.44, 0.54] |
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

*Reading.* The predictors used here capture little perturbation-specific signal: they remove 4–8% of the mean-model error on the development set, 4–6% on RPE1 and 0.1–0.2% on the non-library target set, so the conformity score is largely determined by the actual effect scale of each perturbation (Supplementary Fig. S1). Two consequences must be stated plainly. (i) The argument that heterogeneous coverage is heterogeneous scale is strongest where a residual remains (CV ≈ 0.10–0.15 on the two essential-gene screens) and is **close to definitional on the non-library target set**, where score and effect strength are nearly identical (CV 0.008, Spearman 1.000); conclusions about that dataset therefore rest on the pre-registered coverage comparison, not on this relation. (ii) The relation is a property of these predictors, not of perturbation prediction in general: a model that captured more perturbation-specific signal would lower the ratio and could change the picture. The post-confirmation GEARS coupling test used the same splits and score to ask whether a non-linear model changes this proportionality; it detected no reduction (Supplementary Table S8).

## 6. Coverage on prospectively identifiable subsets, and a hypothesis about network position (post-confirmation, descriptive and exploratory)

*Why this section exists.* Tiers of measured effect strength are defined by the outcome, and under a score that is a monotone function of effect strength their coverage is fixed by the quantile construction (Results; Supplementary Table S4). Subsets defined by information available **before** the experiment are not fixed in that way, so they are reported separately. Split conformal (M0), linear-PCA, k = 50, nominal α = 0.1, 21 grouped splits; 95% intervals from a cluster bootstrap over TSS groups. "Expected from tier composition" is the subset's fraction of T5 units times the observed T5 coverage plus the remainder times one.

| dataset | subset | coverage [95% CI] | T5 fraction | expected from composition | difference |
|---|---|---|---|---|---|
| development | top 10% by predicted effect strength | 0.740 [0.68, 0.80] | 0.60 | 0.709 | +0.031 |
| RPE1 | top 10% by predicted effect strength | 0.708 [0.64, 0.77] | 0.51 | 0.734 | −0.026 |
| non-library target set | top 10% by predicted effect strength | 0.822 [0.80, 0.85] | 0.34 | 0.835 | −0.012 |
| development | top 20% by STRING degree | 0.871 [0.84, 0.90] | 0.39 | 0.812 | +0.059 |
| RPE1 | top 20% by STRING degree | 0.869 [0.84, 0.90] | 0.31 | 0.841 | +0.028 |
| non-library target set | top 20% by STRING degree | 0.861 [0.84, 0.88] | 0.25 | 0.880 | −0.019 |
| development | CORUM complex member | 0.886 [0.87, 0.90] | 0.23 | 0.888 | −0.002 |
| RPE1 | CORUM complex member | 0.890 [0.87, 0.91] | 0.22 | 0.883 | +0.007 |
| non-library target set | CORUM complex member | 0.892 [0.88, 0.90] | 0.23 | 0.891 | +0.001 |

Full tables: Supplementary Tables S12 (prospective subsets) and S13 (CORUM; exploratory).

*Reading.* Subsets a practitioner can name in advance do under-cover — a 10% selection by predicted effect strength is covered 0.740, 0.708 and 0.822 against a nominal 0.9, with intervals excluding 0.9 — and this is the operationally relevant statement. It is, however, not independent evidence of a calibration defect beyond the construction: the deficits agree with what the subsets' T5 composition alone predicts, to within 3.1 percentage points for the proxy-defined subsets and 0.7 for CORUM membership. Complex membership therefore adds nothing once effect strength is accounted for.

*Hypothesis (not tested here).* One subset behaves differently and in the opposite direction: genes with many STRING partners are covered **better** than their effect-strength composition predicts — by +5.9 pp (development) and +2.8 pp (RPE1) in the top degree quintile, and +11.7 pp and +8.6 pp in the top decile, while the non-library target set shows no such excess (−1.9 pp). A possible reading is that for densely connected essential genes the training set contains many functionally similar perturbations, so the fitted model removes relatively more of the error at a given effect strength — consistent with the higher network autocorrelation of effect strength on the essential-gene screens (S1) and with the proxy's higher predictability there. This is a hypothesis generated after the confirmations were unblinded, on subsets chosen by us, with degree defined within each screen's own STRING sub-network (score ≥ 0.7); it is not a pre-registered endpoint, it was not corrected for the number of subsets examined, and it would need its own protocol and an independent screen to test.

# Supplementary Note 2. Execution record of the GEARS run

*Moved from Methods (manuscript v1.0) without change of content.* The server-side steps were carried out by a separate browser-based agent on a rented GPU host, outside this analysis environment; the record below is assembled from its logs and audit files (`gears_audit_v07.tar.gz`, sha256 e92e9bb9…; results `gears_results_v07.tar.gz`, sha256 05a957c4…, both verified on download).

1. *Hardware.* The protocol named an A800 GPU; the run that produced the results used an **NVIDIA GeForce RTX 3060** (recorded in every split's completion record) in a container with a 15 GB memory limit.
2. *Violation of the stop-on-error rule.* The hand-off rules required the executing agent to stop at any error and return it unchanged. The server timing logs (`logs/run_all_timing_failed_attempts_*.txt`) record the sequence. On 2026-10-25 between 23:37 and 2026-10-26 00:07, all 21 splits were killed for memory within about one minute of starting (exit code 137); this is where the rule should have stopped the run. Instead, the agent continued. Three further attempts on split 00 between 11:17 and 11:39 on 2026-10-26 ended after about 4 min each (exit codes 137, 137 and 1); according to the agent's own report the frozen code had been modified at that stage [not confirmable from the records: this rests on the executing agent's own account; the logs do not show what was modified]. After the rules reached it, the agent restored `code/` to the frozen state (bundle checksums verified) and moved all changes into runtime patches loaded from outside the bundle; it also deleted the automatically generated `__pycache__` directory and the residual directory of an unfinished split-00 attempt. A run started at 11:49 under an earlier version of these patches ended each of its four attempted splits with an out-of-memory kill after 79–120 min (13:09–18:32), consistent with the post-training memory failure that patch (c) below addresses; a further attempt on split 00 ended the same way at 20:11 (78 min). How the 11:49 run itself was stopped after its fourth split is not recorded [not confirmable from the records; presumably stopped manually and restarted]. A debug run of split 00 with the epoch count overridden to 1 (2026-10-26 20:22–20:31) was executed in a separate directory with its own results folder (`debug_bundle/`; its script began by stopping any running job); the final run started at 20:36 in the main bundle. None of these attempts produced a completion record in the final results directory: completion records are written only after prediction, every failed attempt ended with a non-zero exit code before that point, the debug run wrote to a separate directory, and the final run executed all 21 splits from the start (each started once, none skipped as already done; split 00 training 4,779 s versus 269 s in the one-epoch debug run). We therefore find no indication that output of any earlier attempt entered the reported results, or that results were selected among attempts. Whether the `debug_bundle/` directory was later deleted is not recorded [not confirmable from the records; it does not affect the final results directory].
3. *Runtime patches* (loaded via `sitecustomize` from `runtime_patches/`; frozen `code/` unchanged): (a) the co-expression and GO similarity graphs were computed in a separate lightweight process with GEARS's own function and its default parameters (co-expression threshold 0.4, k = 20) and read from file during training, which is GEARS's own cache-hit path; (b) the per-epoch evaluation stores predictions in pre-allocated arrays, on the GPU for the large training-set evaluation, with an on-server bitwise comparison against the original in four cases (`test_mem_patch.out`: TEST PASS); (c) the two post-training analysis functions used only for logging were rewritten to avoid densifying the full expression matrix (`test_analysis_patch.out`: PASS); (d) memory was returned to the operating system after evaluation and prediction; (e) at prediction time the control expression matrix was densified once instead of once per perturbation, with the same random-number call per perturbation (on-server check: inputs identical to the original for 8 test perturbations). We read the patch source and found no change to the training procedure, the epoch-selection rule or its inputs. These checks are the executing agent's; we did not rerun them.
4. *Leakage check* (performed locally). The 21 co-expression graph files are distinct; each file's name encodes a hash of its split's sorted training condition list, and recomputing these hashes from the frozen split files matched all 21; the SHA-256 of every file matched the audit manifest; training conditions are disjoint from validation and calibration+test conditions in every split.
5. *Co-expression graphs across machines.* The graphs were computed on a different machine from an earlier attempt on another host. We did not compare the graphs edge by edge across machines. On the run host, the graph for split 00 computed in float64 and the cached graph had the same 13,728 edges; selected correlations lie as close as 2.9 × 10⁻⁵ to the 0.4 threshold, so floating-point differences between machines could move individual edges.
6. *Environment.* `pip freeze` on the server matched the frozen pins (cell-gears 0.1.2, torch 2.3.0 + cu121, torch_geometric 2.5.3, numpy 1.26.4, pandas 2.2.2, scanpy 1.10.1, anndata 0.10.7).
7. *Non-determinism.* GPU training of GEARS is not bitwise reproducible; it was not quantified (no repeated split).

# Supplementary figures

![(descriptive; computed from the stored coverage rows of the confirmatory and development runs, no new model fit). (a) Conformity score S1 of the linear-PCA model against measured effect strength (S1 of the mean-response model) for every test row at nominal α (k = 50, 21 grouped splits; 5,358 / 5,253 / 17,049 rows for the development set, RPE1 and the non-library target set); small points, individual rows; connected points, medians of 12 equal-count bins of effect strength; dotted line, identity. (b) Coefficient of variation of the ratio score / effect strength, by model and dataset; the mean-response model is 0 by construction because its score equals effect strength. Values are in Supplementary Table S6.](figures/figS_scale.png){width=95%}

![(post hoc illustrative analysis, outside the pre-registered protocol). (a) T5 coverage of split-conformal prediction sets (M0; linear-PCA, k = 50, α = 0.1, randomised width matching, 21 grouped splits) when perturbations are tiered by the pre-registered effect strength (filled circles) or by noise-corrected effect strength (open diamonds: effect strength minus the mean effect strength of cell-number-matched pseudo-perturbations drawn from non-targeting control cells); 95% CI from GEE with TSS-group clusters; 1,786 / 1,751 / 5,683 genes. (b) Null effect strength of pseudo-perturbations (200 draws per cell number, 20 cell numbers per dataset, averaged over the 21 split-specific references) against cell number; control cells used: 10,691 (development), 11,485 (RPE1), 20,000 of 75,328 (non-library target set).](figures/figS_noise_floor.png){width=95%}

![(post-confirmation exploratory; protocol explore_modularity_v1.0). (a) Spearman ρ between the prior-neighbour proxy ĥ and measured effect strength within per-dataset tertiles of STRING degree (genome-wide partners, experimental + database score ≥ 0.7); 95% CI from TSS-group cluster bootstrap (2,000 resamples); open symbols, development set. (b) T5 coverage at nominal α = 0.1 (no width matching) for split conformal, the reference (normalisation by measured effect strength) and a split-half reference (width allocated with effect strength from half A of the cells, coverage judged with half B).](figures/figS_explore_modularity.png){width=95%}

# Supplementary tables

| time (UTC unless noted) | event | kind | note |
|---|---|---|---|
| 2026-10-04T06:23:27Z | protocol v0.5 frozen (manifests/freeze_v0.5.json) | freeze | G1 on the development screen |
| 2026-10-04T07:55:54Z | protocol v0.5.1 frozen | freeze | C5 demoted to secondary before any G1 coverage result was read |
| 2026-10-06T10:38:40Z | protocol v0.6 frozen | freeze | design-level: dev = K562 essential, RPE1 primary / non-library target set secondary confirmation, isolation rule |
| 2026-10-06T11:52:14Z | protocol v0.6.1 frozen | freeze | narrowed claims, Mdist/Mscore baselines, width matching; before any dev analysis of v0.6 |
| 2026-10-08T10:36:01Z | post-freeze record (freeze_v0.6.1.json): post_freeze_code_fix | pure code fix | affects conclusions: no |
| 2026-10-08T11:35:41Z | INITIAL READ of development-set responses under v0.6.x (D1 stage 1) | initial data read | d1_stage1.log line 1; v0.6.2 was frozen the same day, before the rerun |
| 2026-10-08T11:53:59Z | protocol v0.6.2 frozen | freeze | learner-unit rule + single E_true definition, after the stage-1 crash and before any dev result existed |
| 2026-10-09T04:36:24Z | post-freeze record (freeze_v0.6.2.json): post_freeze_additions | added analysis module (v0.7 sensitivity), declared | affects conclusions: no (sensitivity only) |
| 2026-10-10 20:58 local | RPE1 confirmatory result written (primary helper) | confirmation result | confirm_v07_result_primary.json |
| 2026-10-10T04:01:30Z | post-freeze record (freeze_v0.6.2.json): post_freeze_code_fix_speedup | pure speed-up, byte-identical output | affects conclusions: no |
| 2026-10-10T04:28:02Z | v0.7-development rerun with the added baselines (ext arrays) | dev analysis | results/d1_v0.7dev |
| 2026-10-10T09:05:55Z | protocol v0.7 frozen (manifests/freeze_v0.7.json) | freeze | confirmation design; user reviewed the full protocol before the freeze |
| 2026-10-10T09:07:37Z | post-freeze record: post_freeze_code_fix | pure code fix | affects conclusions: no |
| 2026-10-10T09:07:49Z | RPE1 C0 started (file hashes; pseudobulk; unit table; splits; response-free power) | confirmation C0 | c0.log; C0 builds expression but computes no effect strength, score or coverage |
| 2026-10-10T09:11:54Z | post-freeze record: post_freeze_code_fix_2 | pure code fix | affects conclusions: no |
| 2026-10-10T11:47:21Z | INITIAL READ of RPE1 responses in model fitting / scores (C1 stage 1, split 00) | initial data read | c1_stage1.log line 1; i.e. after the v0.7 freeze and after C0 |
| 2026-10-10T11:56:20Z | post-freeze record: post_freeze_execution_helper | execution helper calling frozen code | affects conclusions: no |
| 2026-10-10T12:59:19Z | post-freeze record: post_freeze_execution_helper_2 | execution order (primary before descriptive) | affects conclusions: no (primary blocks reproduced exactly later) |
| 2026-10-11 11:02 local | RPE1 full descriptive run; all confirmatory blocks reproduced exactly | confirmation result | confirm_v07_result.json |
| 2026-10-12T05:34:52Z | post-freeze record: post_freeze_execution_resume | engineering: resumable writes | affects conclusions: no (validated byte-identical on dev) |
| 2026-10-12T06:10:52Z | non-library target set C0 started | confirmation C0 | c0.log |
| 2026-10-12T10:54:05Z | post-freeze record: post_freeze_execution_lowmem | engineering: low-memory loader | affects conclusions: no (validated byte-identical on dev) |
| 2026-10-12T10:54:35Z | non-library target set C1 stage 1 started (INITIAL READ of non-library target set responses in scores) | initial data read | c1_resume.log line 1 |
| 2026-10-12T11:41:47Z | post-freeze record: post_freeze_execution_lowmem_2 | engineering: chunked row/column computation | affects conclusions: no (validated byte-identical on dev) |
| 2026-10-14T00:07:51Z | post-freeze record: post_freeze_execution_f1 | engineering: string interning in finish | affects conclusions: no (validated byte-identical on dev; exceeded memory cap on non-library target set, not used) |
| 2026-10-14T12:15:51Z | post-freeze record: post_freeze_execution_f2 | engineering: column store in finish | affects conclusions: no (validated on dev; exceeded cap on non-library target set, not used) |
| 2026-10-15 16:49 local | non-library target set confirmatory result written (primary helper) | confirmation result | confirm_v07_result_primary.json |
| 2026-10-15 17:59 local | non-library target set full descriptive finish; confirmatory blocks matched exactly | confirmation result | confirm_v07_result.json |
| 2026-10-15T08:11:25Z | post-freeze record: post_freeze_execution_helper_gwps_primary | execution helper (non-library target set primary before finish) | affects conclusions: no (all confirmatory blocks matched the later finish exactly) |
| 2026-10-15T09:13:27Z | post-freeze record: post_freeze_execution_f2c | engineering: integer-coded column store in finish | affects conclusions: no (validated byte-identical on dev; used for the non-library target set finish) |
| 2026-10-15T12:23:55Z | post-freeze record: post_freeze_execution_helper_confounding_lowmem | filtered reader for the confounding analysis | affects conclusions: no (validated byte-identical on RPE1) |
| 2026-10-15T12:23:55Z | post-freeze record: post_confirmation_exploratory_modularity | post-confirmation exploratory (own frozen protocol) | affects conclusions: no (exploratory tier) |
| 2026-10-15T13:15:46Z | post-freeze record: post_confirmation_exploratory_gridK_rpe1 | post-confirmation exploratory (own frozen protocol) | affects conclusions: no (exploratory tier) |
| 2026-10-15T13:15:46Z | post-freeze record: illustrative_reliability_analysis | illustrative, outside every protocol | affects conclusions: no (illustrative tier) |
| 2026-10-16T00:16:55Z | post-freeze record: illustrative_noise_floor_analysis | illustrative, outside every protocol | affects conclusions: no (illustrative tier) |
| 2026-10-16T05:25:31Z | post-freeze record: gears_robustness_preparation | robustness check prepared, not run | affects conclusions: not applicable yet |
| 2026-10-19T11:50:36Z | protocol gears_coupling v1.0 frozen (manifests/freeze_gears_coupling_v1.0.json) | freeze | post-confirmation coupling test with GEARS on RPE1; before any GEARS output was read locally |
| 2026-10-28T00:17:18Z | protocol gears_coupling v1.1 frozen (manifests/freeze_gears_coupling_v1.1.json) | freeze | blank, negative and split-half positive controls; before any GEARS metric was computed locally |
| 2026-11-01T08:39:44Z | protocol combo_crosscell v1.2 frozen (manifests/freeze_v1.2.json) | freeze | Norman combinations and Nadig HepG2/Jurkat; no file of these datasets present on disk at freeze |

Table: Freeze timeline and post-freeze deviations (all protocols).

| dataset | split | alpha_lo | alpha_hi | pi | status | W0 | expected_width | expected_over_W0 | realized_width | realized_over_W0 | T5_cov | n_test |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| development | 0 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.44 | 4.61 | 0.72 | 4.61 | 0.72 | 1.000 | 255 |
| development | 1 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.62 | 4.91 | 0.74 | 4.91 | 0.74 | 0.981 | 256 |
| development | 2 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.97 | 4.96 | 0.71 | 4.96 | 0.71 | 1.000 | 255 |
| development | 3 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.79 | 4.91 | 0.72 | 4.91 | 0.72 | 1.000 | 255 |
| development | 4 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.19 | 4.93 | 0.80 | 4.93 | 0.80 | 1.000 | 255 |
| development | 5 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.69 | 4.90 | 0.73 | 4.90 | 0.73 | 1.000 | 255 |
| development | 6 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.30 | 5.32 | 0.84 | 5.32 | 0.84 | 1.000 | 255 |
| development | 7 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.67 | 4.49 | 0.67 | 4.49 | 0.67 | 1.000 | 256 |
| development | 8 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.59 | 4.94 | 0.75 | 4.94 | 0.75 | 0.982 | 255 |
| development | 9 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.68 | 4.87 | 0.73 | 4.87 | 0.73 | 1.000 | 255 |
| development | 10 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.44 | 4.67 | 0.73 | 4.67 | 0.73 | 1.000 | 255 |
| development | 11 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.42 | 4.80 | 0.75 | 4.80 | 0.75 | 1.000 | 255 |
| development | 12 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.15 | 4.75 | 0.77 | 4.75 | 0.77 | 1.000 | 255 |
| development | 13 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.58 | 4.89 | 0.74 | 4.89 | 0.74 | 1.000 | 255 |
| development | 14 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.52 | 4.76 | 0.73 | 4.76 | 0.73 | 1.000 | 256 |
| development | 15 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.44 | 4.74 | 0.74 | 4.74 | 0.74 | 1.000 | 255 |
| development | 16 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.21 | 4.99 | 0.80 | 4.99 | 0.80 | 1.000 | 255 |
| development | 17 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.70 | 4.84 | 0.72 | 4.84 | 0.72 | 1.000 | 255 |
| development | 18 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.12 | 4.92 | 0.80 | 4.92 | 0.80 | 1.000 | 255 |
| development | 19 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.78 | 4.46 | 0.66 | 4.46 | 0.66 | 1.000 | 255 |
| development | 20 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.39 | 4.72 | 0.64 | 4.72 | 0.64 | 1.000 | 255 |
| RPE1 | 0 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.38 | 5.47 | 0.74 | 5.47 | 0.74 | 1.000 | 251 |
| RPE1 | 1 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.58 | 5.83 | 0.89 | 5.83 | 0.89 | 1.000 | 250 |
| RPE1 | 2 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.88 | 6.10 | 0.89 | 6.10 | 0.89 | 1.000 | 250 |
| RPE1 | 3 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.16 | 5.65 | 0.79 | 5.65 | 0.79 | 1.000 | 250 |
| RPE1 | 4 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.02 | 5.47 | 0.78 | 5.47 | 0.78 | 0.978 | 250 |
| RPE1 | 5 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.91 | 5.26 | 0.76 | 5.26 | 0.76 | 1.000 | 250 |
| RPE1 | 6 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.96 | 5.32 | 0.67 | 5.32 | 0.67 | 1.000 | 250 |
| RPE1 | 7 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.73 | 5.41 | 0.70 | 5.41 | 0.70 | 1.000 | 251 |
| RPE1 | 8 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.33 | 5.81 | 0.79 | 5.81 | 0.79 | 1.000 | 250 |
| RPE1 | 9 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.89 | 5.91 | 0.86 | 5.91 | 0.86 | 1.000 | 250 |
| RPE1 | 10 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.42 | 6.29 | 0.98 | 6.29 | 0.98 | 1.000 | 250 |
| RPE1 | 11 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.80 | 5.25 | 0.77 | 5.25 | 0.77 | 1.000 | 250 |
| RPE1 | 12 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.85 | 5.72 | 0.83 | 5.72 | 0.83 | 1.000 | 250 |
| RPE1 | 13 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.23 | 5.55 | 0.77 | 5.55 | 0.77 | 1.000 | 250 |
| RPE1 | 14 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.56 | 5.89 | 0.90 | 5.89 | 0.90 | 1.000 | 251 |
| RPE1 | 15 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 7.32 | 5.50 | 0.75 | 5.50 | 0.75 | 1.000 | 250 |
| RPE1 | 16 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.82 | 5.53 | 0.81 | 5.53 | 0.81 | 1.000 | 250 |
| RPE1 | 17 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.79 | 6.61 | 0.97 | 6.61 | 0.97 | 1.000 | 250 |
| RPE1 | 18 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.98 | 5.40 | 0.77 | 5.40 | 0.77 | 1.000 | 250 |
| RPE1 | 19 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 6.93 | 5.72 | 0.83 | 5.72 | 0.83 | 1.000 | 250 |
| RPE1 | 20 | 0.010 | 0.015 | 0.927 | matched | 7.54 | 7.54 | 1.00 | 7.55 | 1.00 | 1.000 | 250 |
| non-library target set | 0 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.72 | 3.88 | 0.68 | 3.88 | 0.68 | 0.984 | 812 |
| non-library target set | 1 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.73 | 3.80 | 0.66 | 3.80 | 0.66 | 0.976 | 812 |
| non-library target set | 2 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.39 | 3.68 | 0.68 | 3.68 | 0.68 | 0.979 | 812 |
| non-library target set | 3 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.51 | 3.79 | 0.69 | 3.79 | 0.69 | 0.982 | 812 |
| non-library target set | 4 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.56 | 3.68 | 0.66 | 3.68 | 0.66 | 1.000 | 812 |
| non-library target set | 5 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.63 | 3.64 | 0.65 | 3.64 | 0.65 | 1.000 | 811 |
| non-library target set | 6 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.17 | 3.90 | 0.75 | 3.90 | 0.75 | 0.994 | 812 |
| non-library target set | 7 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.45 | 3.64 | 0.67 | 3.64 | 0.67 | 0.993 | 812 |
| non-library target set | 8 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.43 | 3.79 | 0.70 | 3.79 | 0.70 | 0.977 | 812 |
| non-library target set | 9 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.58 | 3.78 | 0.68 | 3.78 | 0.68 | 0.980 | 812 |
| non-library target set | 10 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.41 | 3.77 | 0.70 | 3.77 | 0.70 | 0.988 | 812 |
| non-library target set | 11 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.58 | 3.79 | 0.68 | 3.79 | 0.68 | 0.994 | 813 |
| non-library target set | 12 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.60 | 3.77 | 0.67 | 3.77 | 0.67 | 0.994 | 811 |
| non-library target set | 13 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.64 | 3.78 | 0.67 | 3.78 | 0.67 | 0.994 | 811 |
| non-library target set | 14 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.39 | 3.75 | 0.70 | 3.75 | 0.70 | 0.988 | 812 |
| non-library target set | 15 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.65 | 3.79 | 0.67 | 3.79 | 0.67 | 0.988 | 812 |
| non-library target set | 16 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.40 | 3.74 | 0.69 | 3.74 | 0.69 | 1.000 | 812 |
| non-library target set | 17 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.57 | 3.83 | 0.69 | 3.83 | 0.69 | 0.994 | 812 |
| non-library target set | 18 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.42 | 3.79 | 0.70 | 3.79 | 0.70 | 0.982 | 812 |
| non-library target set | 19 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.45 | 3.71 | 0.68 | 3.71 | 0.68 | 0.993 | 812 |
| non-library target set | 20 | 0.010 | 0.010 | 1.000 | unmatched_too_narrow | 5.57 | 3.78 | 0.68 | 3.78 | 0.68 | 0.976 | 811 |

Table: Reference method (width allocated by measured effect strength): per-split width matching. W0 = mean M0 width; α′ = matched nominal level range; status unmatched_too_narrow means no α′ on the grid reached M0's width. Where π = 1 (all development splits), the matched level is a single α′ and the expected width equals the realised width by definition; the development-set rows were filled accordingly (the development pipeline did not store the expected width separately).

| dataset | item | value | verdict | note |
|---|---|---|---|---|
| RPE1 | B0: M0 T5 coverage | 0.480 (95% CI 0.43 to 0.53) (n=1051, gee) | pass | rule: point estimate < 0.85 and CI upper < 0.90 |
| RPE1 | P: rho(h, E_true) | 0.43 (95% CI 0.39 to 0.47) | P+ | n_genes=1751 |
| RPE1 | P: AUROC(T5) | 0.74 (95% CI 0.71 to 0.77) | P+ | P+ needs rho_lo>=0.30 and auroc_lo>=0.65; P- if rho_hi<0.20 or auroc_hi<0.60 |
| RPE1 | step1 E1 (T5, M1 - M0) | +9.7 (95% CI 5.8 to 13.6) pp | pass | margin 5 pp; SE=0.0201; logistic GEE + delta; n=1051 |
| RPE1 | step1 E2 (top-10% by h, M1 - M0) | +17.4 (95% CI 12.6 to 22.3) pp; subset width ratio 1.38 | pass | coverage 0.883 vs 0.708; n=528; T5 share 0.51 |
| RPE1 | step1 overall | passed | pass | E1 and E2 both at/above margin with CI lower bound > 0 |
| RPE1 | step1 E1 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=0.974; 0.50x_dev=0.686 (SE from C0 scaling=0.0206) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| RPE1 | step1 E1 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=0.977; 0.50x_dev=0.691 (observed SE=0.0201) | post hoc | same dev estimates, SE from this dataset |
| RPE1 | step1 E2 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=0.997; 0.50x_dev=0.873 (SE from C0 scaling=0.0236) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| RPE1 | step1 E2 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=0.996; 0.50x_dev=0.861 (observed SE=0.0248) | post hoc | same dev estimates, SE from this dataset |
| RPE1 | step2 E1 (T5, M1 - Mdist-M) | +10.9 (95% CI 7.8 to 14.0) pp | pass | margin 3 pp; SE=0.0158; logistic GEE + delta; n=1051 |
| RPE1 | step2 E2 (top-10% by h, M1 - Mdist-M) | +14.6 (95% CI 10.4 to 18.7) pp; subset width ratio 1.25 | pass | coverage 0.883 vs 0.737; n=528; T5 share 0.51 |
| RPE1 | step2 overall | passed | pass | E1 and E2 both at/above margin with CI lower bound > 0 |
| RPE1 | step2 E1 power - C0 pre-specified | 1.00x_dev=0.976; 0.75x_dev=0.839; 0.50x_dev=0.502 (SE from C0 scaling=0.0161) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| RPE1 | step2 E1 power - post hoc at observed SE | 1.00x_dev=0.980; 0.75x_dev=0.852; 0.50x_dev=0.518 (observed SE=0.0158) | post hoc | same dev estimates, SE from this dataset |
| RPE1 | step2 E2 power - C0 pre-specified | 1.00x_dev=0.993; 0.75x_dev=0.905; 0.50x_dev=0.559 (SE from C0 scaling=0.0138) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| RPE1 | step2 E2 power - post hoc at observed SE | 1.00x_dev=0.857; 0.75x_dev=0.622; 0.50x_dev=0.328 (observed SE=0.0212) | post hoc | same dev estimates, SE from this dataset |
| RPE1 | guard: nominal_marginal_in_band | True | pass |  |
| RPE1 | guard: M_perm | True | pass |  |
| RPE1 | guard: matching_primary_comparators | True | pass |  |
| RPE1 | guard: P_not_minus | True | pass |  |
| RPE1 | nominal marginal coverage (primary, unmatched) | 0.900 | pass | band [0.88, 0.92] |
| RPE1 | M_perm mean E1 / permutations passing step 1 | −5.2 pp / 0 of 20 | pass | rule: mean < 2 pp and at most 1 passing |
| RPE1 | final classification | steps 1 and 2 passed: the deployable prior-neighbour proxy repairs part of the failure and adds to the strongest distance / difficulty baseline |  |  |
| non-library target set | B0: M0 T5 coverage | 0.518 (95% CI 0.49 to 0.55) (n=3417, gee) | pass | rule: point estimate < 0.85 and CI upper < 0.90 |
| non-library target set | P: rho(h, E_true) | 0.11 (95% CI 0.08 to 0.14) | P- | n_genes=5683 |
| non-library target set | P: AUROC(T5) | 0.57 (95% CI 0.55 to 0.59) | P- | P+ needs rho_lo>=0.30 and auroc_lo>=0.65; P- if rho_hi<0.20 or auroc_hi<0.60 |
| non-library target set | step1 E1 (T5, M1 - M0) | +1.0 (95% CI −0.1 to 2.1) pp | fail | margin 5 pp; SE=0.0057; logistic GEE + delta; n=3417 |
| non-library target set | step1 E2 (top-10% by h, M1 - M0) | +4.4 (95% CI 3.1 to 5.7) pp; subset width ratio 1.09 | fail | coverage 0.866 vs 0.822; n=1722; T5 share 0.34 |
| non-library target set | step1 overall | not passed | fail | E1 and E2 both at/above margin with CI lower bound > 0 |
| non-library target set | step1 E1 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=0.809 (SE from C0 scaling=0.0115) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| non-library target set | step1 E1 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=0.960 (observed SE=0.0057) | post hoc | same dev estimates, SE from this dataset |
| non-library target set | step1 E2 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=0.980 (SE from C0 scaling=0.0131) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| non-library target set | step1 E2 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=1.000 (observed SE=0.0066) | post hoc | same dev estimates, SE from this dataset |
| non-library target set | step2 E1 (T5, M1 - Mdist-M) | +1.8 (95% CI 0.8 to 2.9) pp | fail | margin 3 pp; SE=0.0053; logistic GEE + delta; n=3417 |
| non-library target set | step2 E2 (top-10% by h, M1 - Mdist-M) | +2.6 (95% CI 1.6 to 3.5) pp; subset width ratio 1.04 | fail | coverage 0.866 vs 0.841; n=1722; T5 share 0.34 |
| non-library target set | step2 overall | not passed | fail | E1 and E2 both at/above margin with CI lower bound > 0 |
| non-library target set | step2 E1 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=0.975; 0.50x_dev=0.572 (SE from C0 scaling=0.0089) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| non-library target set | step2 E1 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=0.999; 0.50x_dev=0.620 (observed SE=0.0053) | post hoc | same dev estimates, SE from this dataset |
| non-library target set | step2 E2 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=0.991; 0.50x_dev=0.606 (SE from C0 scaling=0.0077) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| non-library target set | step2 E2 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=0.667 (observed SE=0.0048) | post hoc | same dev estimates, SE from this dataset |
| non-library target set | guard: nominal_marginal_in_band | True | pass |  |
| non-library target set | guard: M_perm | True | pass |  |
| non-library target set | guard: matching_primary_comparators | True | pass |  |
| non-library target set | guard: P_not_minus | False | FAIL |  |
| non-library target set | nominal marginal coverage (primary, unmatched) | 0.903 | pass | band [0.88, 0.92] |
| non-library target set | M_perm mean E1 / permutations passing step 1 | −0.9 pp / 0 of 20 | pass | rule: mean < 2 pp and at most 1 passing |
| non-library target set | final classification | step 1 not passed: repair not confirmed |  |  |

Table: Confirmatory decision tables (protocol v0.7) for RPE1 and the non-library target set; values from results/conf_v0.7/<dataset>/confirm_v07_result.json and c0_report.json.

| dataset | model | n_test | n_T5 | T5_cov | marg_cov | cv_ratio | resid_var_share_log | corr_logratio_logE | constructive_ref_meanmodel | mc_reference | excess_vs_mc_pp | excess_vs_constructive_pp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| development | linear_go | 5,358 | 1,075 | 0.526 | 0.904 | 0.148 | 0.193 | −0.30 | 0.511 | 0.517 | +0.9 | +1.5 |
| development | linear_pca | 5,358 | 1,075 | 0.516 | 0.903 | 0.107 | 0.112 | −0.24 | 0.511 | 0.517 | −0.0 | +0.6 |
| development | mean | 5,358 | 1,075 | 0.511 | 0.902 | 0.000 | 0.000 |  | 0.511 | 0.517 | −0.6 | +0.0 |
| non-library target set | linear_go | 17,049 | 3,417 | 0.517 | 0.903 | 0.014 | 0.002 | −0.15 | 0.520 | 0.504 | +1.3 | −0.3 |
| non-library target set | linear_pca | 17,049 | 3,417 | 0.518 | 0.903 | 0.008 | 0.001 | −0.14 | 0.520 | 0.504 | +1.4 | −0.1 |
| non-library target set | mean | 17,049 | 3,417 | 0.520 | 0.904 | 0.000 | 0.000 |  | 0.520 | 0.504 | +1.6 | +0.0 |
| RPE1 | linear_go | 5,253 | 1,051 | 0.504 | 0.901 | 0.116 | 0.094 | −0.38 | 0.493 | 0.508 | −0.4 | +1.1 |
| RPE1 | linear_pca | 5,253 | 1,051 | 0.480 | 0.896 | 0.097 | 0.065 | −0.39 | 0.493 | 0.508 | −2.8 | −1.3 |
| RPE1 | mean | 5,253 | 1,051 | 0.493 | 0.899 | 0.000 | 0.000 |  | 0.493 | 0.508 | −1.5 | +0.0 |

Table: Constructive reference and observed T5 coverage by model (M0, k = 50). constructive_ref_meanmodel = T5 coverage of the mean model, whose score equals E_true; mc_reference = Monte Carlo reference from scripts/constructive_reference_mc.py (seed 20261001, 4,000 replicates); resid_var_share_log = R.

| dataset | n_units_total | n_eligible_in_R_and_GO | n_included_E | n_excluded_cells_lt30 | frac_excluded | ncells_excluded_median | ncells_included_median | E_excluded | E_included | mannwhitney_p | frac_excluded_above_T5_cut | spearman_E_logncells_eligible |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| development | 2,057 | 1,863 | 1,786 | 77 | 0.041 | 22.00 | 122.00 | {'median': 7.854, 'q25': 6.195, 'q75': 9.77} | {'median': 3.557, 'q25': 2.779, 'q75': 4.917} | 9e-37 | 0.961 | −0.35 |
| RPE1 | 2,393 | 2,101 | 1,751 | 350 | 0.167 | 19.00 | 79.00 | {'median': 6.603, 'q25': 5.55, 'q75': 8.371} | {'median': 3.892, 'q25': 2.92, 'q75': 5.19} | 6e-104 | 0.731 | −0.64 |
| non-library target set | 7,813 | 5,772 | 5,683 | 89 | 0.015 | 21.00 | 177.00 | {'median': 9.58, 'q25': 8.456, 'q75': 11.946} | {'median': 3.222, 'q25': 2.714, 'q75': 4.066} | 2e-56 | 1.000 | −0.81 |

Table: Selection of evaluation units: perturbations excluded for fewer than 30 cells versus included units (E_excluded/E_included = median measured effect strength).

| dataset | model | n_rows | cv_ratio | median_ratio | q05 | q95 | spearman_S1_E | rel_error_reduction_vs_mean | rel_reduction_T5 | rel_reduction_T1 | r2_log |
|---|---|---|---|---|---|---|---|---|---|---|---|
| development | linear_go | 5,358 | 0.148 | 0.98 | 0.61 | 1.04 | 0.89 | 0.085 | 0.137 | 0.010 | 0.807 |
| development | linear_pca | 5,358 | 0.107 | 0.99 | 0.81 | 1.04 | 0.94 | 0.043 | 0.076 | −0.004 | 0.888 |
| development | mean | 5,358 | 0.000 | 1 | 1.00 | 1.00 | 1 | 0 | 0 | 0 | 1.000 |
| non-library target set | linear_go | 17,049 | 0.014 | 1.00 | 0.99 | 1.01 | 1.00 | 0.002 | 0.006 | −0.000 | 0.998 |
| non-library target set | linear_pca | 17,049 | 0.008 | 1.00 | 0.99 | 1.01 | 1.00 | 0.001 | 0.003 | −0.000 | 1.000 |
| non-library target set | mean | 17,049 | 0.000 | 1 | 1.00 | 1.00 | 1 | 0 | 0 | 0 | 1.000 |
| RPE1 | linear_go | 5,253 | 0.116 | 0.99 | 0.71 | 1.06 | 0.94 | 0.065 | 0.109 | −0.016 | 0.906 |
| RPE1 | linear_pca | 5,253 | 0.097 | 0.99 | 0.81 | 1.07 | 0.96 | 0.043 | 0.079 | −0.034 | 0.935 |
| RPE1 | mean | 5,253 | 0.000 | 1 | 1.00 | 1.00 | 1 | 0 | 0 | 0 | 1.000 |

Table: Relation between the conformity score S1 and measured effect strength by model (descriptive). cv_ratio = coefficient of variation of S1/E_true; rel_error_reduction_vs_mean = error reduction relative to the mean model; r2_log = R² of log S1 on log E_true.

| trace level | entries |
|---|---|
| value-level: value present in the designated file; up to four or five candidate fields listed | 370 |
| field-level: file and field recorded when the number was written | 320 |
| field-level: unique field in the designated file | 101 |
| field-level: several fields share the value; the field matching the context is given | 62 |
| protocol parameter or pre-specified threshold | 29 |
| repeated from manuscript v1.0 (traced in the v1.0 rows) | 6 |
| protocol constant (1.4826) | 1 |

Table: Number-to-source map (summary). The full map (889 entries: every decimal number in the main text of manuscript v1.0, Supplementary Note 1 and the figure captions, plus every decimal number in the sections added in version 1.1) gives for each number its context, source field(s) and the archived result file(s) in which the value is present at the stated precision; provided as supplementary_table_S7_number_source_map.csv. No number is untraced; for the value-level entries the value is confirmed in the designated file but the specific field is ambiguous because several fields share it.

| Quantity | Value |
|---|---|
| R (share of log S1 not explained by log E_true) | linear-PCA 0.065; GEARS 0.059; ratio 0.91 (0.76 to 1.11) |
| CV of S1/E_true | linear-PCA 0.097; GEARS 0.112; ratio 1.16 (1.06 to 1.27) |
| Error reduction vs mean model, all | linear-PCA 4.3%; GEARS 6.3%; +1.9 pp (1.5 to 2.4) |
| Error reduction vs mean model, T5 | linear-PCA 7.9%; GEARS 9.7%; +1.8 pp (0.8 to 2.8) |
| Test rows | linear-PCA 5,253; GEARS 5,253; same (split, gene) rows |
| M0 T5 coverage (GEARS scores) | GEARS 0.498; analytic reference 0.5 |
| M0 marginal coverage (GEARS scores) | GEARS 0.899; nominal 0.9 |
| Pre-specified verdict (v1.0) | not reduced |
| Control (v1.1): R, split-half predictor (95% CI) | RPE1 0.830 (0.79 to 0.87); development 0.980 (0.96 to 0.99) |
| Control (v1.1): R, linear-PCA | RPE1 0.039; development 0.104 |
| Control (v1.1): R, GEARS | RPE1 0.036 |
| Control (v1.1): Error reduction, control − linear-PCA, all | RPE1 +28.9 pp (27.5 to 30.3); development +26.2 pp (24.6 to 27.9) |
| Control (v1.1): Error reduction, control − linear-PCA, T5 | RPE1 +44.9 pp (42.3 to 47.6); development +45.5 pp (42.3 to 48.5) |
| Control (v1.1): Positive-control verdict | RPE1 threshold reachable; development threshold reachable |
| Control (v1.1): Negative control (20 permutations): linear-PCA verdicts | RPE1 not reduced 20 |
| Control (v1.1): Negative control (20 permutations): GEARS verdicts | RPE1 not reduced 19, coupling metrics up with accuracy down 1 |
| Control (v1.1): Mean error reduction of permuted predictions (outside protocol) | RPE1 linear-PCA 0.6%; GEARS 3.4% |

Table: Post-confirmation coupling test with GEARS on RPE1 (protocols gears_coupling v1.0 and v1.1). Values from results/gears_coupling_v1.1/gears_coupling_result.json and controls_summary.json; 95% CIs from TSS-group cluster bootstrap.

| Quantity | mean-response model | additive baseline | linear-PCA (pair = sum of gene embeddings) |
|---|---|---|---|
| R (95% CI) | 0 (by construction) | 0.662 (0.54 to 0.78) | 0.486 (0.38 to 0.60) |
| CV of S1/E_true (95% CI) | 0 | 0.396 (0.35 to 0.44) | 0.285 (0.25 to 0.32) |
| Error reduction vs mean model, all (95% CI) | 0 | 0.498 (0.46 to 0.54) | 0.339 (0.30 to 0.38) |
| Error reduction vs mean model, T5 (95% CI) | 0 | 0.605 (0.54 to 0.66) | 0.440 (0.35 to 0.51) |
| Pre-specified verdict | reference | R clearly elevated | R clearly elevated |
| Marginal coverage (95% CI) | 0.897 (0.84 to 0.93) | 0.904 (0.85 to 0.94) | 0.868 (0.81 to 0.91) |
| T1 coverage (95% CI) | 1.000 (all covered) | 1.000 (all covered) | 1.000 (all covered) |
| T2 coverage (95% CI) | 1.000 (all covered) | 1.000 (all covered) | 1.000 (all covered) |
| T3 coverage (95% CI) | 1.000 (all covered) | 0.987 (0.92 to 1.00) | 0.974 (0.91 to 0.99) |
| T4 coverage (95% CI) | 0.962 (0.89 to 0.99) | 0.823 (0.66 to 0.92) | 0.949 (0.88 to 0.98) |
| T5 coverage (95% CI) | 0.532 (0.39 to 0.67) | 0.722 (0.54 to 0.85) | 0.430 (0.29 to 0.59) |
| T5 coverage minus mean model, pp (95% CI) | – | +19.0 (0.0 to 38.0) | −10.1 (−28.6 to 8.3) |
| Mean conformal radius q̂ over splits | 11.86 | 5.32 | 6.77 |
| Median relative width q̂/E_true, T1 / T5 | 3.37 / 1.01 | 1.50 / 0.46 | 1.93 / 0.58 |
| Negative control 1 (E_true permuted, 20×): mean R; mean T5 − T1 coverage (pp) | 0.997; −0.7 | 0.996; −0.2 | 0.996; −0.7 |
| Negative control 2 (predictions permuted, 20×): mean R; mean error reduction; verdicts | – | 0.465; −0.145; 20 × 'R elevated without accuracy gain' | 0.314; −0.082; 20 × 'R elevated without accuracy gain' |
| Half-data basis: R (95% CI) | 0 | 0.603 (0.47 to 0.73) | 0.386 (0.28 to 0.51) |
| Split-half control (half A predicts half B): R (95% CI); error reduction (95% CI) | 0.818 (0.66 to 0.94); 0.451 (0.40 to 0.50) |  |  |
| Two-way gene resampling (post hoc, v1.3): R 95% CI | 0 | 0.44 to 0.86 | 0.27 to 0.69 |
| Two-way gene resampling: error reduction 95% CI | 0 | 0.41 to 0.57 | 0.23 to 0.44 |
| Two-way gene resampling: T5 coverage 95% CI | 0.26 to 0.83 | 0.44 to 0.97 | 0.16 to 0.71 |
| Two-way gene resampling: T5 coverage minus mean model, pp, 95% CI | – | −18.8 to 51.9 | −47.2 to 26.7 |
| One-way gene resampling (sensitivity): R 95% CI | 0 | 0.57 to 0.76 | 0.39 to 0.58 |
| Post hoc verdict with two-way intervals | – | R clearly elevated | R clearly elevated |

Table: Combinatorial perturbations (Norman et al.; protocol combo_crosscell v1.2): 129 gene pairs, 21 unseen-combination splits, 387 test predictions (79 in T5, from 27 pairs). Coverage intervals from an intercept-only binomial GEE clustered by pair (undefined when every row is covered); all other intervals from a pair-level cluster bootstrap with 2,000 replicates. Values from results/v1.2/norman/norman_v12_result.json.

| Setting | predictor | pairs / test predictions | R (95% CI) | error reduction vs mean model (95% CI) | CV of S1/E_true (95% CI) |
|---|---|---|---|---|---|
| one gene held out | linear-PCA, K = 50 | 129 / 618 | 0.442 (0.24 to 0.70) | 0.188 (0.11 to 0.28) | 0.325 (0.21 to 0.43) |
| one gene held out | linear-PCA, K = 10 | 129 / 618 | 0.279 (0.05 to 0.59) | 0.110 (0.03 to 0.22) | 0.212 (0.09 to 0.33) |
| one gene held out | measured single-gene component only | 129 / 618 | 0.468 (0.29 to 0.62) | 0.267 (0.20 to 0.33) | 0.290 (0.21 to 0.37) |
| both genes held out | linear-PCA, K = 50 | 63 / 78 | 0.302 (0.03 to 0.83) | −0.057 (−0.24 to 0.04) | 0.315 (0.08 to 0.53) |
| both genes held out | linear-PCA, K = 10 | 63 / 78 | 0.273 (0.01 to 0.83) | −0.007 (−0.16 to 0.11) | 0.291 (0.03 to 0.54) |

Table: Unseen-component ablation on the Norman screen (post hoc; protocol v1.3, frozen before computation). The 71 genes forming the evaluated pairs were split into 5 folds, 3 times; all perturbations involving a held-out gene were removed from training. Intervals from two-way gene-level resampling (2,000 replicates). Values from results/v1.3/norman/ablation_result.json; two-way resampling of the v1.2 predictions from results/v1.3/norman/bootstrap_result.json. Pre-written reading for linear-PCA with both genes held out: uncertain.

| Quantity | HepG2 | Jurkat |
|---|---|---|
| Evaluated target genes (TSS groups); test predictions | 1,632 (1,569); 4,896 | 1,908 (1,835); 5,724 |
| linear-PCA: R (95% CI) | 0.086 (0.07 to 0.11) | 0.038 (0.03 to 0.05) |
| linear-PCA: CV of S1/E_true (95% CI) | 0.097 (0.09 to 0.11) | 0.063 (0.06 to 0.07) |
| linear-PCA: error reduction vs mean model (95% CI) | 0.051 (0.04 to 0.06) | 0.023 (0.02 to 0.03) |
| linear-PCA: error reduction in T5 (95% CI) | 0.117 (0.10 to 0.13) | 0.038 (0.03 to 0.05) |
| Pre-specified verdict (prediction: R not elevated) | R not elevated | R not elevated |
| Marginal coverage, mean model / linear-PCA | 0.906 / 0.906 | 0.904 / 0.903 |
| T5 coverage, mean model (95% CI) | 0.531 (0.48 to 0.58) | 0.521 (0.47 to 0.57) |
| T5 coverage, linear-PCA (95% CI) | 0.532 (0.48 to 0.58) | 0.515 (0.47 to 0.56) |
| T5 coverage, linear-PCA minus mean model, pp (95% CI) | 0.1 (−3.9 to 4.2) | −0.6 (−3.2 to 2.0) |
| Negative control 1 (E_true permuted, 20×): mean R; mean T5 − T1 coverage, pp | 1.000; −0.0 | 1.000; −0.6 |
| Negative control 2 (predictions permuted, 20×): verdicts | R not elevated 20 | R not elevated 20 |
| Split-half control: R (95% CI) | 0.943 (0.91 to 0.97) | 0.702 (0.65 to 0.75) |
| K selected = 50 (of 21 splits) | 21 | 18 |

Table: Two further cell lines (Nadig et al.; protocol combo_crosscell v1.2 S3.4 with amendment v1.2.1; post-confirmation, frozen before the data were read). Coverage intervals from an intercept-only binomial GEE clustered by transcription-start-site group (undefined when every row is covered); other intervals from a TSS-group cluster bootstrap with 2,000 replicates. Values from results/v1.2/hepg2/nadig_v12_result.json and results/v1.2/jurkat/nadig_v12_result.json; data checks in results/v1.2/{hepg2,jurkat}/qc_v121.json.

| dataset | subset | prospective | n_rows | coverage | mean_E_true | frac_of_T5 | ci_lo | ci_hi | expected_from_tier_mix | deficit_beyond_mix |
|---|---|---|---|---|---|---|---|---|---|---|
| development | T5 (outcome-defined, top 20% by measured effect) | 0 | 1,075 | 0.516 | 7.76 | 1.000 | 0.47 | 0.56 | 0.516 | 0.000 |
| development | T1-T4 (outcome-defined) | 0 | 4,283 | 1.000 | 3.50 | 0.000 | 1.00 | 1.00 | 1.000 | 0.000 |
| development | top 10% by h (prospective) | 1 | 546 | 0.740 | 6.60 | 0.603 | 0.68 | 0.80 | 0.709 | 0.031 |
| development | top 20% by h (prospective) | 1 | 1,074 | 0.775 | 6.03 | 0.510 | 0.74 | 0.81 | 0.753 | 0.022 |
| development | top 20% by STRING degree (prospective) | 1 | 1,086 | 0.871 | 5.23 | 0.389 | 0.84 | 0.90 | 0.812 | 0.059 |
| development | top 10% by STRING degree (prospective) | 1 | 552 | 0.913 | 5.25 | 0.422 | 0.88 | 0.94 | 0.796 | 0.117 |
| development | bottom 20% by STRING degree (prospective) | 1 | 1,175 | 0.946 | 3.73 | 0.104 | 0.92 | 0.97 | 0.950 | −0.004 |
| non-library target set | T5 (outcome-defined, top 20% by measured effect) | 0 | 3,417 | 0.518 | 6.12 | 1.000 | 0.49 | 0.54 | 0.518 | 0.000 |
| non-library target set | T1-T4 (outcome-defined) | 0 | 13,632 | 1.000 | 3.10 | 0.000 | 1.00 | 1.00 | 1.000 | 0.000 |
| non-library target set | top 10% by h (prospective) | 1 | 1,722 | 0.822 | 4.35 | 0.343 | 0.80 | 0.85 | 0.835 | −0.012 |
| non-library target set | top 20% by h (prospective) | 1 | 3,423 | 0.859 | 4.06 | 0.286 | 0.84 | 0.88 | 0.862 | −0.003 |
| non-library target set | top 20% by STRING degree (prospective) | 1 | 3,621 | 0.861 | 3.93 | 0.250 | 0.84 | 0.88 | 0.880 | −0.019 |
| non-library target set | top 10% by STRING degree (prospective) | 1 | 1,796 | 0.837 | 4.05 | 0.276 | 0.81 | 0.87 | 0.867 | −0.030 |
| non-library target set | bottom 20% by STRING degree (prospective) | 1 | 6,222 | 0.918 | 3.60 | 0.180 | 0.91 | 0.93 | 0.914 | 0.005 |
| RPE1 | T5 (outcome-defined, top 20% by measured effect) | 0 | 1,051 | 0.480 | 8.43 | 1.000 | 0.43 | 0.53 | 0.480 | 0.000 |
| RPE1 | T1-T4 (outcome-defined) | 0 | 4,202 | 1.000 | 3.78 | 0.000 | 1.00 | 1.00 | 1.000 | 0.000 |
| RPE1 | top 10% by h (prospective) | 1 | 528 | 0.708 | 7.12 | 0.511 | 0.64 | 0.77 | 0.734 | −0.026 |
| RPE1 | top 20% by h (prospective) | 1 | 1,053 | 0.769 | 6.35 | 0.435 | 0.73 | 0.81 | 0.774 | −0.004 |
| RPE1 | top 20% by STRING degree (prospective) | 1 | 1,063 | 0.869 | 5.40 | 0.306 | 0.84 | 0.90 | 0.841 | 0.028 |
| RPE1 | top 10% by STRING degree (prospective) | 1 | 536 | 0.927 | 5.10 | 0.304 | 0.89 | 0.96 | 0.842 | 0.086 |
| RPE1 | bottom 20% by STRING degree (prospective) | 1 | 1,194 | 0.931 | 4.05 | 0.137 | 0.90 | 0.95 | 0.929 | 0.002 |

Table: M0 coverage on prospectively identifiable subsets versus the coverage expected from their T5 composition (post-confirmation, descriptive).

| dataset | subset | prospective | n_rows | n_genes | coverage | ci_lo | ci_hi | mean_E_true | frac_of_T5 | expected_from_tier_mix | deficit_beyond_mix |
|---|---|---|---|---|---|---|---|---|---|---|---|
| development | CORUM complex member (prospective) | 1 | 3,552 | 1,184 | 0.886 | 0.87 | 0.90 | 4.58 | 0.232 | 0.888 | −0.002 |
| development | not a CORUM member (prospective) | 1 | 1,806 | 602 | 0.936 | 0.92 | 0.95 | 3.90 | 0.138 | 0.933 | 0.003 |
| RPE1 | CORUM complex member (prospective) | 1 | 3,462 | 1,154 | 0.890 | 0.87 | 0.91 | 4.92 | 0.224 | 0.884 | 0.007 |
| RPE1 | not a CORUM member (prospective) | 1 | 1,791 | 597 | 0.907 | 0.88 | 0.93 | 4.30 | 0.154 | 0.920 | −0.013 |
| non-library target set | CORUM complex member (prospective) | 1 | 5,718 | 1,906 | 0.892 | 0.88 | 0.90 | 3.81 | 0.227 | 0.891 | 0.001 |
| non-library target set | not a CORUM member (prospective) | 1 | 11,331 | 3,777 | 0.909 | 0.90 | 0.92 | 3.65 | 0.187 | 0.910 | −0.001 |

Table: M0 coverage of CORUM complex members and non-members versus composition-expected coverage (post-confirmation, exploratory).

| Gate | Joung et al. 2023 (combinatorial TF overexpression, hESC) | Wessels et al. 2023 (CaRPool-seq, Cas13d knockdown, THP-1) |
|---|---|---|
| G1: ≥ 60 evaluable pairs | **fails**: 44 doubles (and 3 triples) designed | passes by design: 158 pairs; cells per pair to be verified from the file |
| G2: ≥ 90% of pairs with both singles measured | passes by design (10 singles) | passes by design (26 regulators measured alone) |
| G3: ≥ 500 non-targeting control cells (new gate) | not assessable from the paper | not assessable from the paper |
| G4: raw integer counts with per-cell perturbation labels | not assessable from the paper | not assessable from the paper |
| G5: ≥ 10 distinct genes in pairs (new gate) | borderline: 10 | passes: 26 |
| Data-level check (protocol v1.4 S3; metadata and counts only) | not read | file identity confirmed (30,707 cells; paper 31,308); G1 136 evaluable pairs (pass); **G2 136/157 = 86.6% (fails)**; **G3 424 control cells (fails)**; G4 integer counts (pass); G5 24 genes (pass) |
| Outcome | not included (G1) | **not analysed** (fails G2 and G3 on the data; conclusions restricted to Norman as pre-specified) |
| Outcome without the two new gates (G3, G5) | not included (G1) | not analysed (G2 alone fails) |


Table: Applicability gates for a second combinatorial screen (protocol v1.4): first assessed from the original papers only, then, for the screen that passed on paper, verified on the data file (metadata and counts only; no score computed). Data-level values from results/v1.4/wessels/gates_v14.json. Sources: Joung et al., Cell 2023, doi:10.1016/j.cell.2022.11.026; Wessels et al., Nat Methods 2023, doi:10.1038/s41592-022-01705-x. Supplementary note: Joung et al. report that averaging single-TF profiles predicted combination profiles better than regression, so a summation baseline may be ill-suited to overexpression; this is noted but was not the reason for exclusion. Gates G3 and G5 were set in the same session in which the papers were read; removing them leaves both outcomes unchanged. The rule for the case that only one dataset passes was added after the gate outcomes were known and before any data were read.
