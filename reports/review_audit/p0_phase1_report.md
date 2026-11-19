# 外审意见 P0 统计实现核查：第一阶段报告（只读，未重跑、未改稿）

日期：2026-10-17。核查的代码哈希与冻结清单 `manifests/freeze_v0.7.json` 一致：`code_sha256` 下各文件一致；`data_v07.py` 与 `post_freeze_code_fix_2.sha256_after` 一致；`confirm_run_v07.py` 与 `post_freeze_execution_f2c` 一致。

## 停止原因：发现一处实现与稿件描述的实质不一致（P0-7）

- **稿件写法**：Results 写的是 "An oracle … restores T5 coverage almost completely (**width-matched** T5 gain +52.0 pp … RPE1 and +47.0 pp … GWPS; T5 coverage 0.999 and 0.988)"；Discussion 写的是 "≈ +47–52 pp **under width matching**"；补充说明 1 的表中也有 "oracle E1"。
- **实际实现**：oracle 并没有做到宽度匹配。在 `confirm_v07_result.json → matching.oracle.status` 中：
  - RPE1：21 次划分里有 20 次为 `unmatched_too_narrow`，1 次为 matched；
  - GWPS：21/21 次为 `unmatched_too_narrow`；
  - 开发集：21/21 次为 `unmatched_too_narrow`。

  `unmatched_too_narrow` 的意思是：即使取网格上最小的 α（0.01，也就是最宽的集合），oracle 的平均宽度仍小于 M0 的宽度。实际宽度与 M0 宽度之比（`realized_width_ratio_range`）为：RPE1 0.67–1.00，GWPS 0.65–0.75，开发集 0.64–0.84。oracle 的边际覆盖为 RPE1 0.994、GWPS 0.992。所以，报告的 oracle 增益是在 α′ = 0.01、集合比 M0 更窄的条件下得到的，不是在等宽条件下得到的。
- **原因：oracle 归一化后的分数几乎是常数。** 计算 linear-PCA 的 S1 / E_true（M0 行，`confirm_rand_rows.parquet`）：

| | 中位数 | 5%–95% 分位 | 变异系数 | Spearman(S1, E_true) |
|---|---|---|---|---|
| 开发集 | 0.994 | 0.813–1.043 | 0.107 | 0.941 |
| RPE1 | 0.991 | 0.807–1.071 | 0.097 | 0.963 |
| GWPS | 1.000 | 0.990–1.007 | **0.008** | **1.000** |

  对均值模型来说，按定义 S1 = E_true，比值恒为 1，oracle 完全退化。在 GWPS 上，linear-PCA 的分数与 E_true 几乎完全相同，说明它在 GWPS 上的预测几乎等于训练质心，接近均值模型。

### 对已有结论的影响

1. **确认性分类不受影响。** oracle 不参与 B0、P、第 1 步、第 2 步或任何保护条件（见 `confirm_v07.analyze`），其余 15 种方法在全部划分上都是 `matched`。
2. **描述性的 oracle 数字（+52.0 / +47.0 / +48.2 pp）必须改称**为"oracle 在 α′ = 0.01 下、宽度 ≤ M0 时的 T5 覆盖"，不能再称"宽度匹配"。它们是"已知真实效应尺度时"的上界，而且对 S1 而言**基本是由定义决定的**：S1 ≈ E_true，所以 S1/E_true ≈ 1，用 E_true 分配宽度几乎必然覆盖所有单元。这与你记录的改稿原则一致，但"width-matched"这一具体错误必须改正。"oracle repair generalises across datasets"应删除；它反映的是分数的构造，不是可推广的修复。
3. **一个更深的问题（不属于实现不一致，但审稿人一定会追问）。** 在 GWPS 上 S1 与 E_true 的 Spearman 为 1.000，因此 GWPS 的"T5 欠覆盖"几乎就是"均值模型分数按效应强度排序"的直接结果。"The failure therefore replicated … across a disjoint set of non-essential target genes"这句话在 GWPS 上的信息量很低，需要收窄措辞，同时报告 linear-PCA 在 GWPS 上相对均值模型几乎没有预测能力。
4. **补充说明 1 中另有一处表述错误**："lowered oracle T5 coverage only slightly (… RPE1 0.990 → 0.994)"——RPE1 的数值是上升的，不是下降。

建议：改稿时只修改措辞，并把 oracle 的匹配状态和宽度比写进 Methods 与补充表；不重跑，不改冻结结果。是否另外补一个"oracle 在等宽条件下"的结果，例如在网格之外放宽 α 的下界，留到最后决定是否补实验时再定。这属于方案之外的描述性分析。

## 逐项核查

### P0-1：S1 集合与"宽度"
- **公式**：s_p = (y_p − μ_tr)/σ，ŝ_p = (ŷ_p − μ_tr)/σ。T_k(p) 是 |s_pg| 最大的 k = 50 个读出基因，排除靶基因所在列。S1(p) = sqrt(mean_{g∈T_k(p)} (s_pg − ŝ_pg)²)。集合定义为 C(p) = {y : S1(p; y) ≤ t_p}。**T_k(p) 由观测到的响应 y_p 选出**，因此 C(p) 无法在实验前写成表达空间中的显式集合；它只能用来审计覆盖与否。"宽度"是半径 t_p 在测试单元上的平均值（S1 单位）；t = +∞ 时截断为 2 × max(校准分数)，t < 0 时取 0。
- **代码**：`models_v04.score_S1`（L69）；`models_v04.train_reference`（L12；μ 与 σ = 1.4826·MAD 只用当次训练池中细胞数 ≥ 30 的扰动，σ 的下限为其 5% 分位数）；`proxy_v06.e_true`（L54；E_true 即均值模型的 S1）；`methods_v06.capped_width`（L218）。
- **输出**：`confirm_rand_rows.parquet` 中的列 `score`、`thr`、`w`、`covered`、`W0`、`expected_width`。
- **与稿件是否一致**：公式一致。但 Results 和 Discussion 把 C(p) 称为 "prediction sets"，没有说明 T_k 依赖测试响应、因而只用于审计。**需要补写说明（不是实现错误）。**

### P0-2：宽度匹配
- **流程**（`confirm_v07.rand_rows`，L34；网格为 `methods_v06.ALPHA_GRID`，0.01–0.40，步长 0.005，共 79 点）：
  - W0 = capped_width(M0 在 α = 0.1 时的阈值)，只依赖校准分数；
  - 对每个候选方法和每个 α，计算其阈值的截断平均宽度 W(α)。阈值由校准分数、校准单元的 ĥ 和**测试单元的 ĥ** 决定（M1：t_i = q̂_α(S_cal/ĥ_cal)·ĥ_i）；
  - 用二分法找到使 W(lo) ≥ W0 ≥ W(hi) 的相邻两个网格点，令 π = (W0 − W(hi))/(W(lo) − W(hi))；
  - 每个测试单元按独立均匀数 v_i < π 取 lo 阈值，否则取 hi 阈值，使期望平均宽度等于 W0。
- **使用的信息**：校准分数、校准与测试单元的代理值（ĥ、d̂ 等）、固定种子的均匀随机数。**不使用测试响应和测试分数**；代码中 `s_test` 只出现在最后的 `covered = s_test <= t` 一行。匹配用到了测试单元的 ĥ，属于协变量层面的转导，Methods 中应写明。
- **保证的适用范围**：固定 α 时，M1 满足边际覆盖 ≥ 1 − α，因为 ĥ 只依赖先验和训练数据，(ĥ, S) 在校准与测试单元之间可交换。但 α′ 是依据校准分数（以及测试单元的 ĥ）选出的，所以**匹配后的集合在任何预先给定的水平上都没有有限样本保证**，只能作为等宽比较的工具。匹配后边际覆盖的观测值为：RPE1 中 M1 0.906、M0 0.896；GWPS 中 M1 0.906、M0 0.903。选出的 α′ 范围为：RPE1 0.06–0.115，GWPS 0.085–0.11。
- **匹配质量**：`matching.*.max_rel_dev_expected_width` 约 1e−16；实际宽度与 W0 之比为 0.996–1.002。唯一的例外是 oracle（见上文）。
- **与稿件是否一致**：Methods 描述的网格和随机化与实现一致。缺少两点：匹配使用了测试单元的 ĥ；匹配后不再有名义水平的保证。**需补写。**

### P0-3：TSS 整组分配与 GEE 聚类
- **随机对象**：`splits_v06.repeated_group_kfold`（L32）每轮以种子 base + r 打乱 TSS 组，按单元数均衡地分到 7 个测试折；其余单元再按 TSS 组随机分为 val、cal、train（`_fill`）。覆盖保证依赖的是**整组随机分配所带来的 TSS 组之间的可交换性**，严格来说可交换单位是 TSS 组，不是基因。测试中每个 TSS 组含多个基因的情况为：开发集 65 组含 2 个基因（另有 1,656 组为单基因）；RPE1 53 组；GWPS 232 组含 2 个、1 组含 3 个。
- **GEE**：`analysis_v06.gee_tiers`（L59）、`paired_T5`（L81）、`decision_paired`（L181）都使用 logistic GEE，独立工作相关矩阵，`groups = tss_group`。每个基因在 3 轮中各被测试一次，3 行都落在同一个 TSS 组簇内（核对结果：`rows_per_gene` 三个数据集均为 3），所以同一基因跨划分的重复已被聚类覆盖。
- **未覆盖的依赖**：同一划分内的测试单元共用同一个 q̂ 和同一个训练模型；不同划分之间训练集和校准集有重叠。GEE 没有把划分当作聚类单元，CI 可能偏窄。v0.5 的 GLMM 有 (1|split)，从 v0.6 起不再包含。这属于方法局限，可作为补充敏感性分析的候选（例如按 TSS 组与划分做双向聚类，或按划分做块 bootstrap），由你决定。
- **与稿件是否一致**：Methods 写 "The exchangeable unit is the target gene"，**不准确**，应改为"可交换单位是 TSS 组（96–97% 为单基因）"。GEE 的描述一致。

### P0-4：ĥ
- **公式**：ĥ = exp(f(x))，f 是以 log E_true 为目标的 ridge 回归（`proxy_v06.fit_proxy`，L208；冻结配置 `sim_prior|ridge`）。因此 **ĥ 由构造保证为正**，尺度与 E_true（S1 单位）相同，S1/ĥ 无量纲。
- **流程**：学习器只用当次划分的 E 训练单元训练（v0.6.2 规则），always-train 单元只作为邻居。训练单元按 TSS 组分为 5 个内折（`inner_folds`，L135，种子来自 `_seed(MASTER_SEED, 6000, split)`）。每折的特征只用折外邻居计算；学习器在 OOF 特征上拟合；val/cal/test 单元的特征用全部训练池计算，再由全量学习器预测（`eval_v06.compute_split` 的 proxies 段，以及 `eval_v06.proxies`，L173）。ĥ_OOF 来自在其余 4 折上重新拟合的学习器，只用于确定 c 和 Mondrian 的分箱边界。主方法 γ = 0，c = 0，因此 **M1 不使用 ĥ_OOF**。
- **观测范围**（测试行）：RPE1 2.50–9.27，GWPS 2.82–4.90，开发集 1.98–8.58。
- **与稿件是否一致**：一致。Methods 应补一句"ĥ = exp(ridge prediction)，因而为正"。

### P0-7：oracle
- **分母**：`METHODS["oracle"] = ("norm", "E")`（`methods_v06.py` L172），P["E"] = 校准与测试单元的 `cache["E_calt"]`，也就是用全部细胞测得的 E_true。冻结配置为 γ = 0（`FROZEN.hypers.oracle = "0.0"`），所以分数是 S1/E_true。
- **分布**：见上文表格。均值模型完全退化（比值恒为 1）；linear-PCA 在开发集和 RPE1 上近似退化（变异系数约 0.1），在 GWPS 上几乎完全退化（变异系数 0.008）。
- **结论**：见"停止原因"。

### P0-8：GWPS 集合
- **定义**：`confirm_run_v07.do_c0`（L247–249）把 `p0v3_units_all_targets.csv` 中 K562 essential 筛选的**全部 2,057 个靶基因**作为 `exclude_targets`。`data_v07.build_units`（L115）先删除这些靶基因（包括 K562 essential 中的非 E 单元），再对 GWPS 其余靶基因应用与开发集相同的 E 规则：靶基因在读出基因中、在 GEARS GO 图中、细胞数 ≥ 30。结果为 5,683 个 E 单元、2,001 个 always-train 单元。
- **与稿件是否一致**：Methods 的写法（"from which all targets of the K562 essential screen were removed"）与实现一致。Results 和 Discussion 中的 "non-essential targets / non-essential genes" **命名不准确**：实现只保证靶基因"不在 K562 essential 筛选文库中"，并没有对必需性重新判定。按你的原则，统一改名为例如 "GWPS targets outside the K562 essential-screen library"。

### P0-5、P0-6
你的消息里写了"综合意见附后"，但我只看到 P0-1、2、3、4、7、8 的重点说明，没有看到审稿意见原文；在本项目的文件中也没有找到。**P0-5 和 P0-6 的内容我不知道，所以没有核查**，请把审稿意见原文发给我。

## 尚未进行
第二阶段（决策表）和第三阶段（数值溯源表）按你的要求暂停，等你对上述不一致作出决定后再做。
