# 第二阶段：确认性决策表（RPE1 与 GWPS）

日期：2026-10-17。全部数值取自冻结代码的输出（`results/conf_v0.7/<dataset>/confirm_v07_result.json`），没有重算。主模型 linear-PCA，k = 50，α = 0.1，随机化宽度匹配，21 次划分。功效分两栏：**C0 事先功效**取自 `c0_report.json` 的 `power` 字段（在看到任何覆盖结果之前按开发集效应量与按单元数缩放的 SE 计算）；**事后情景功效**取自 `confirm_v07_result.json` 的 `power_at_observed_se`（同样的开发集效应量，但用本数据集实测的 SE）。两者都不是"观测到的效应的功效"。

## RPE1（主确认集）

| 项目 | 数值 | 判定 | 说明 |
|---|---|---|---|
| B0: M0 T5 coverage | 0.480 [0.431, 0.528] (n=1051, gee) | pass | rule: point estimate < 0.85 and CI upper < 0.90 |
| P: rho(h, E_true) | 0.427 [0.387, 0.467] | P+ | n_genes=1751 |
| P: AUROC(T5) | 0.742 [0.714, 0.772] | P+ | P+ needs rho_lo>=0.30 and auroc_lo>=0.65; P- if rho_hi<0.20 or auroc_hi<0.60 |
| step1 E1 (T5, M1 - M0) | +9.7 [5.8, 13.6] pp | pass | margin 5 pp; SE=0.0201; logistic GEE + delta; n=1051 |
| step1 E2 (top-10% by h, M1 - M0) | +17.4 [12.6, 22.3] pp; subset width ratio 1.38 | pass | coverage 0.883 vs 0.708; n=528; T5 share 0.51 |
| step1 overall | passed | pass | E1 and E2 both at/above margin with CI lower bound > 0 |
| step1 E1 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=0.974; 0.50x_dev=0.686 (SE from C0 scaling=0.0206) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| step1 E1 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=0.977; 0.50x_dev=0.691 (observed SE=0.0201) | post hoc | same dev estimates, SE from this dataset |
| step1 E2 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=0.997; 0.50x_dev=0.873 (SE from C0 scaling=0.0236) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| step1 E2 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=0.996; 0.50x_dev=0.861 (observed SE=0.0248) | post hoc | same dev estimates, SE from this dataset |
| step2 E1 (T5, M1 - Mdist-M) | +10.9 [7.8, 14.0] pp | pass | margin 3 pp; SE=0.0158; logistic GEE + delta; n=1051 |
| step2 E2 (top-10% by h, M1 - Mdist-M) | +14.6 [10.4, 18.7] pp; subset width ratio 1.25 | pass | coverage 0.883 vs 0.737; n=528; T5 share 0.51 |
| step2 overall | passed | pass | E1 and E2 both at/above margin with CI lower bound > 0 |
| step2 E1 power - C0 pre-specified | 1.00x_dev=0.976; 0.75x_dev=0.839; 0.50x_dev=0.502 (SE from C0 scaling=0.0161) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| step2 E1 power - post hoc at observed SE | 1.00x_dev=0.980; 0.75x_dev=0.852; 0.50x_dev=0.518 (observed SE=0.0158) | post hoc | same dev estimates, SE from this dataset |
| step2 E2 power - C0 pre-specified | 1.00x_dev=0.993; 0.75x_dev=0.905; 0.50x_dev=0.559 (SE from C0 scaling=0.0138) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| step2 E2 power - post hoc at observed SE | 1.00x_dev=0.857; 0.75x_dev=0.622; 0.50x_dev=0.328 (observed SE=0.0212) | post hoc | same dev estimates, SE from this dataset |
| guard: nominal_marginal_in_band | True | pass |  |
| guard: M_perm | True | pass |  |
| guard: matching_primary_comparators | True | pass |  |
| guard: P_not_minus | True | pass |  |
| nominal marginal coverage (primary, unmatched) | 0.900 | pass | band [0.88, 0.92] |
| M_perm mean E1 / permutations passing step 1 | -5.2 pp / 0 of 20 | pass | rule: mean < 2 pp and at most 1 passing |
| final classification | steps 1 and 2 passed: the deployable prior-neighbour proxy repairs part of the failure and adds to the strongest distance / difficulty baseline |  |  |

## GWPS（次确认集）

| 项目 | 数值 | 判定 | 说明 |
|---|---|---|---|
| B0: M0 T5 coverage | 0.518 [0.490, 0.546] (n=3417, gee) | pass | rule: point estimate < 0.85 and CI upper < 0.90 |
| P: rho(h, E_true) | 0.110 [0.085, 0.138] | P- | n_genes=5683 |
| P: AUROC(T5) | 0.569 [0.550, 0.589] | P- | P+ needs rho_lo>=0.30 and auroc_lo>=0.65; P- if rho_hi<0.20 or auroc_hi<0.60 |
| step1 E1 (T5, M1 - M0) | +1.0 [-0.1, 2.1] pp | fail | margin 5 pp; SE=0.0057; logistic GEE + delta; n=3417 |
| step1 E2 (top-10% by h, M1 - M0) | +4.4 [3.1, 5.7] pp; subset width ratio 1.09 | fail | coverage 0.866 vs 0.822; n=1722; T5 share 0.34 |
| step1 overall | not passed | fail | E1 and E2 both at/above margin with CI lower bound > 0 |
| step1 E1 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=0.809 (SE from C0 scaling=0.0115) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| step1 E1 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=0.960 (observed SE=0.0057) | post hoc | same dev estimates, SE from this dataset |
| step1 E2 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=0.980 (SE from C0 scaling=0.0131) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| step1 E2 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=1.000 (observed SE=0.0066) | post hoc | same dev estimates, SE from this dataset |
| step2 E1 (T5, M1 - Mdist-M) | +1.8 [0.8, 2.9] pp | fail | margin 3 pp; SE=0.0053; logistic GEE + delta; n=3417 |
| step2 E2 (top-10% by h, M1 - Mdist-M) | +2.6 [1.6, 3.5] pp; subset width ratio 1.04 | fail | coverage 0.866 vs 0.841; n=1722; T5 share 0.34 |
| step2 overall | not passed | fail | E1 and E2 both at/above margin with CI lower bound > 0 |
| step2 E1 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=0.975; 0.50x_dev=0.572 (SE from C0 scaling=0.0089) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| step2 E1 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=0.999; 0.50x_dev=0.620 (observed SE=0.0053) | post hoc | same dev estimates, SE from this dataset |
| step2 E2 power - C0 pre-specified | 1.00x_dev=1.000; 0.75x_dev=0.991; 0.50x_dev=0.606 (SE from C0 scaling=0.0077) | pre-specified | dev estimate x {1.00, 0.75, 0.50} |
| step2 E2 power - post hoc at observed SE | 1.00x_dev=1.000; 0.75x_dev=1.000; 0.50x_dev=0.667 (observed SE=0.0048) | post hoc | same dev estimates, SE from this dataset |
| guard: nominal_marginal_in_band | True | pass |  |
| guard: M_perm | True | pass |  |
| guard: matching_primary_comparators | True | pass |  |
| guard: P_not_minus | False | FAIL |  |
| nominal marginal coverage (primary, unmatched) | 0.903 | pass | band [0.88, 0.92] |
| M_perm mean E1 / permutations passing step 1 | -0.9 pp / 0 of 20 | pass | rule: mean < 2 pp and at most 1 passing |
| final classification | step 1 not passed: repair not confirmed |  |  |

## 两点说明

1. **GWPS 的第 2 步是在第 1 步未通过的情况下计算的**，`step2.tested = false`。按方案的层级检验，第 1 步不通过时第 2 步不作判定，表中数值仅供参考。
2. **GWPS 的保护条件 `P_not_minus` 未满足**（ĥ 的预测力落入预设的 P− 类）。按 `classify()` 的规则，最终分类为"第 1 步未通过：修复未得到确认"，并附加降级标注。这与方案一致，不作任何回溯修改。
