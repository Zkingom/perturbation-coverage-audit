# 改稿对照摘要（diff）

日期：2026-10-17。原稿保留，新稿为 Results v0.7、Methods v0.2、Discussion v0.4、主图 fig1_main_v4（含新图注）。确认性数值一处未改。

## Results v0.6 -> v0.7

| 操作 | 原文（节选） | 新文（节选） |
|---|---|---|
| replace | # Results (draft v0.6, 2026-10-15; development set, RPE1 and GWPS; exploratory analysis summarised, details in Supplementary Note 1) | # Results (draft v0.7, 2026-10-17; development set, RPE1 and the non-library target set; exploratory analysis summarised, details in Supplementary Note 1) |
| replace | *Draft for internal review. Development-set numbers carry selection optimism (methods and thresholds were chosen on the same units). For GWPS the confirmatory a | *Evidence tiers used throughout, as defined in Methods: **confirmatory** (pre-registered endpoint of a frozen protocol, computed once on a confirmation set), ** / *Draft for internal review. Development-set numbers carry selection optimism (methods and thresholds were chosen on the same units). For the non-library targ |
| replace | We calibrated perturbation-level prediction sets for a linear low-rank perturbation model (linear-PCA; score S1, the σ-scaled RMS error on the 50 most perturbat | We calibrated perturbation-level prediction sets for a linear low-rank perturbation model (linear-PCA; score S1, the σ-scaled RMS error on the 50 most perturbat |
| replace | Under randomised width matching to M0 (same expected mean interval width), none of six deployable distance- or uncertainty-based proxies improved T5 coverage in | Under randomised width matching to M0 (same expected mean interval width), none of six deployable distance- or uncertainty-based proxies improved T5 coverage in |
| replace | An oracle that normalises the score by the measured effect strength restores T5 coverage almost completely (width-matched T5 gain +52.0 pp (95% CI 47.1 to 56.8) | **Heterogeneous coverage is heterogeneous scale (descriptive).** The conformity score is close to proportional to measured effect strength: the ratio S1/E_true  |
| replace | *Essential genes (RPE1; cross-cell-line).* ĥ predicted effect strength out of sample (Spearman ρ = 0.43, 95% CI 0.39–0.47; AUROC for T5 0.74). Normalising the c | *Essential-gene screen (RPE1; cross-cell-line).* ĥ predicted effect strength out of sample (Spearman ρ = 0.43, 95% CI 0.39–0.47; AUROC for T5 0.74). Normalising |
| replace | *Non-essential genes (GWPS; cross-gene-set).* ĥ barely predicted effect strength (ρ = 0.11, 95% CI 0.08–0.14; AUROC 0.57, 95% CI 0.55–0.59; class P− under the p | *Non-library target set (cross-gene-set).* ĥ barely predicted effect strength (ρ = 0.11, 95% CI 0.08–0.14; AUROC 0.57, 95% CI 0.55–0.59; class P− under the pre- |
| replace | Taken together, the failure axis (effect strength) and its oracle repair generalise across datasets, whereas the prior-neighbour proxy transfers across cell lin | Taken together, the failure along effect strength replicated in both confirmation sets, whereas the prior-neighbour proxy transferred across cell lines within e |
| replace | In analyses specified after unblinding (protocol frozen before running; Supplementary Note 1, Supplementary Fig. S_explore), effect strength was highly repeatab | In analyses specified after unblinding (protocol frozen before running; Supplementary Note 1, Supplementary Fig. S_explore), effect strength was highly repeatab |
| replace | In a logistic model of M0 miscoverage, log effect strength retained a strong positive association after adjustment for split-half noise and reliability, cell nu | In a logistic model of M0 miscoverage, log effect strength retained a strong positive association after adjustment for split-half noise and reliability, cell nu |
| replace | - RPE1 targets largely overlap the development set (cross-cell-line test); GWPS tests a disjoint gene set but the same cell line as the development set. / - For GWPS, the confirmatory analysis was computed before the descriptive run because of local memory limits (recorded deviation of execution order); the full r | - RPE1 targets largely overlap the development set (cross-cell-line test); the non-library target set tests a disjoint gene set but the same cell line as the de / - For the non-library target set, the confirmatory analysis was computed before the descriptive run because of local memory limits (recorded deviation of exe |
| replace | - K and d reached the upper grid bound in all 21 RPE1 splits, but only in 11 of 21 GWPS splits each; a post-confirmation exploratory analysis with an extended g | - K and d reached the upper grid bound in all 21 RPE1 splits, but only in 11 of 21 non-library-target-set splits each. A post-confirmation exploratory analysis  / - Evaluation units require at least 30 cells. The perturbations excluded by this rule have systematically higher measured effect strength than the included o |
| replace | Main Figure 1 = fig1_main_v3 (panels a–c: confirmatory and descriptive results only). The exploratory panels (proxy predictability by STRING-degree tertile; spl | Main Figure 1 = fig1_main_v4 (panels a–c: confirmatory and descriptive results only; the measured-effect-strength reference method was removed from panel c beca |

## Methods v0.1 -> v0.2

| 操作 | 原文（节选） | 新文（节选） |
|---|---|---|
| replace | # Methods (draft v0.1, 2026-10-15) | # Methods (draft v0.2, 2026-10-17) / **Evidence tiers.** Results are labelled as **confirmatory** (a pre-registered endpoint of a protocol frozen with SHA-256 manifests before any response of that  |
| replace | We used the genome-scale CRISPRi Perturb-seq screens of Replogle et al.¹ (raw single-cell and bulk files, Figshare+ doi:10.25452/figshare.plus.20029387.v1): K56 | We used the genome-scale CRISPRi Perturb-seq screens of Replogle et al.¹ (raw single-cell and bulk files, Figshare+ doi:10.25452/figshare.plus.20029387.v1): K56 |
| replace | **Units and leakage control.** The exchangeable unit is the target gene. Evaluation units ("E units") were targets with at least 30 cells that met the pre-speci | **Units and leakage control.** The exchangeable unit is the target gene. Evaluation units ("E units") were targets with at least 30 cells that met the pre-speci |
| replace | - **Secondary and reference methods**: Mscore-N (all ĥ features fitted directly to the out-of-fold score), M2 (Mondrian on ĥ), M3 (conditional conformal of Gibb | - **Secondary and reference methods**: Mscore-N (all ĥ features fitted directly to the out-of-fold score), M2 (Mondrian on ĥ), M3 (conditional conformal of Gibb |
| replace | **Width matching.** Because conditional-coverage gains can be bought with wider sets, methods were compared at equal mean width: for each split, α of each metho | **Width matching.** Because conditional-coverage gains can be bought with wider sets, methods were compared at equal mean width: for each split, α of each metho |
| insert | （无） | **Selection of evaluation units and its bias.** Evaluation units require at least 30 cells; perturbations failing only this rule remain in training. Because thi / **Functional enrichment (development set only; descriptive).** Target sets were compared against the background of all evaluation units of the screen with on |
| replace | After both confirmations were unblinded we ran an exploratory analysis whose protocol was frozen before running (explore_modularity_v1.0): split-half reliabilit | After both confirmations were unblinded we ran an exploratory analysis whose protocol was frozen before running (explore_modularity_v1.0): split-half reliabilit |
| replace | Master seed 20261001; all seeds derived deterministically per split. Every post-freeze change was either (i) pure engineering verified byte-identical (or bitwis | Master seed 20261001; all seeds derived deterministically per split. Every post-freeze change was either (i) pure engineering verified byte-identical (or bitwis |

## Discussion v0.3 -> v0.4

| 操作 | 原文（节选） | 新文（节选） |
|---|---|---|
| replace | # Discussion (draft v0.3, 2026-10-15) | # Discussion (draft v0.4, 2026-10-15) |
| replace | **Marginal validity hides a systematic blind spot.** Across three Perturb-seq screens, two cell lines and both essential and non-essential genes, split-conforma | **Marginal validity hides a systematic blind spot.** Across three Perturb-seq screens, two cell lines and both an essential-gene library and a set of targets ou |
| replace | **Distance is the wrong axis.** Out-of-distribution work on perturbation prediction emphasises how far a test perturbation is from the training data. Six deploy | **Distance-based proxies did not capture this axis.** Out-of-distribution work on perturbation prediction emphasises how far a test perturbation is from the tra |
| replace | **The bottleneck is predicting effect strength before the experiment.** An oracle that knows effect strength restores T5 coverage (≈ +47–52 pp under width match | **Heterogeneous coverage is heterogeneous scale.** The conformity score is nearly proportional to measured effect strength (coefficient of variation of S1/E_tru |
| replace | **The repair has a cost.** Normalising by ĥ at equal total width moves width from perturbations predicted to be weak to those predicted to be strong. Perturbati | **The repair has a cost.** Normalising by ĥ at equal total width moves width from perturbations predicted to be weak to those predicted to be strong. Perturbati |
| replace | 1. *Circularity.* ĥ is built from neighbours' effect strengths, so "low network autocorrelation" and "low ĥ predictability" overlap conceptually and the stateme / 2. *Reliability is not validity.* Split-half reliability measures repeatability across cell halves; components shared by both halves (for example a cell-numb | 1. *Circularity.* ĥ is built from neighbours' effect strengths, so "low network autocorrelation" and "low ĥ predictability" overlap conceptually and the stateme / 2. *Reliability is not validity.* Split-half reliability measures repeatability across cell halves; components shared by both halves (for example a cell-numb |
| replace | 4. *Scope of the confirmation.* All three screens come from one laboratory and platform, and two are in K562; most RPE1 targets overlap the development targets, | 4. *Selection of evaluation units.* Perturbations with fewer than 30 cells are excluded from evaluation, and these have systematically higher measured effect st / 5. *Scope of the confirmation.* All three screens come from one laboratory and platform, and two are in K562; most RPE1 targets overlap the development targe |
| replace | **Conclusion.** Uncertainty statements for perturbation-response predictions should be audited along effect strength, not only distance. Marginal conformal guar | **Conclusion.** Uncertainty statements for perturbation-response predictions should be audited along effect strength, not only distance. Marginal conformal guar |

## Fig1 caption v3 -> v4

| 操作 | 原文（节选） | 新文（节选） |
|---|---|---|
| replace | Figure 1 (draft v3). (a) Coverage of split-conformal prediction sets (randomised width matching, M0) by quintile of measured effect strength; linear-PCA, α = 0. | Figure 1 (draft v4). (a) Coverage of split-conformal prediction sets (randomised width matching, M0) by quintile of measured effect strength; linear-PCA, α = 0. |

## 本轮改稿要点

1. **oracle 重新定位**：Results 新增"覆盖异质性即尺度异质性"一段，给出 S1/E_true 的变异系数（0.107 / 0.097 / 0.008）与 Spearman（0.94 / 0.96 / 1.00），说明对均值模型是定义性的、在非文库靶基因集上近乎定义性；删除全部"宽度匹配后 oracle 增益 +52.0 / +47.0 pp"的表述，改为给出实际宽度比（0.67–1.00 / 0.65–0.75）并说明未完成宽度匹配。Discussion 小标题由"瓶颈在于预测效应强度"改为"覆盖异质性即尺度异质性"。主图 c 删除 oracle 标注，并在图注中说明原因。
2. **宽度匹配的保证范围**：Methods 写明选定 α′ 后理论保证为 1 − α′（RPE1 的 M1 为 0.060–0.115），宽度匹配是比较方法如何分配固定宽度预算的评估手段，不构成 1 − 0.1 的覆盖保证；同时写明"宽度"是封顶半径。
3. **弱化的表述**：`distance is the wrong axis` → `distance-based proxies did not capture this axis`，并限定在六个代理、一个分数、两个线性模型；`rules out measurement noise` → `an unlikely sole explanation`，并补上"可靠性把两半共享的成分算作信号"；`bounds any repair to about 2 pp` → 改为"置信区间与至多约 2 pp 的真实效应相容（限本代理、本分数、本模型）"；`The proxy's predictability explains the difference` → `differs in the same direction`；结论段删除"the fix requires predicting effect strength"，改为尺度信息的表述。
4. **Mscore-N 的归因**改为引用开发集的特征配置比较（描述性），不再由 M1 ≈ Mscore-N 直接推断"增益来自先验邻居特征"。
5. **命名统一**：全文将 GWPS 统一为 "non-library target set"（K562 全基因组筛选剔除 essential 文库全部 2,057 个靶基因后的集合），并在 Methods 写明它按文库归属定义，而非生物学意义上的非必需。
6. **三稿状态统一**：GWPS 混杂分析由"pending"改为已完成；K/d 网格扩展由"planned"改为已完成并给出结果。
7. **新增内容**：Methods 补入富集分析方法（背景集、检验、多重校正、注释版本、三个基因集，以及度数与复合体共线性的说明）；Methods 与 Discussion 补入细胞数 <30 的入选偏倚（描述性、方向不可定号）；Methods 开头与 Results 开头加入四档证据等级定义。
8. **新增补充项**：Table S_oracle（逐划分匹配表）、Table S_timeline（冻结时间线与偏离清单）。
