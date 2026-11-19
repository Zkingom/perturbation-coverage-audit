# 确认后探索性分析方案草案：效应强度在先验网络上的模块化与 ĥ 的可预测性

版本：**v1.0（2026-10-15 冻结，运行前）**。v0.1 草案经用户确认，并按用户要求补充第 8 节（a–d）和第 9 节（实现细节）。冻结时，本分析的任何结果都尚未计算。

性质：**确认后的探索性分析**，在 v0.7 的 RPE1 主确认和 GWPS 次确认都揭盲之后才提出。结果不修改、也不重新解释任何确认性结论。论文中只能作为探索性结果呈现，并注明提出时间晚于确认结果。

## 1. 问题

GWPS（非 essential 基因）上 ĥ 的 ρ = 0.11（P−），RPE1（essential 基因）上 ρ = 0.43。这里检验一个假设：

> 效应强度在先验网络（STRING、CORUM 复合体）上是否呈模块化分布，即相邻基因的效应强度相近。ĥ 的可预测性是否取决于这种模块化的程度。

**可证伪的预测**：
- P-a：RPE1 上效应强度的网络自相关 ρ_net 高于 GWPS。
- P-b：在两套数据中，复合体成员和 STRING 高度数基因上 ĥ 的 ρ 和 AUROC 都更高。
- P-c：如果按度数分层后，两套数据在相同度数层内的 ĥ 可预测性差距明显缩小，就说明差异有一部分可以用网络位置的组成差异来解释；如果差距不变，这个解释就不成立。

## 2. 数据与已有量（不新增任何模型拟合）

- RPE1 用 `results/conf_v0.7/rpe1`，GWPS 用 `results/conf_v0.7/gwps`。两者都直接读 21 个 stage-1 缓存里已有的 ĥ（冻结配置 h = `sim_prior|ridge`，主模型 linear-PCA）和 E_true，读取方式与 `confirm_v07.analyze` 中 `ev.proxy_table` 相同。
- 开发集（K562 essential）只作参照，并标注"选择性乐观"。
- **不重新拟合模型，不改动 ĥ，不加入新特征。**

## 3. 分组（预先固定，看到结果后不改）

1. **G-CORUM**：靶基因是否属于 CORUM 人类复合体。使用 enrichment v0.7 已登记的 `corum_human_current.txt`，其 sha256 已记录在 `results/enrich_v07` 的 annotation_manifest 中。分为两组：成员 / 非成员。
2. **G-deg**：STRING 度数。使用与先验相同的 STRING v12.0 文件和同一通道，取 max(experimental, database) ≥ 0.7 的全基因组伙伴数。这与富集分析 `string_degree` 的定义相同。
   - 主分组：在**每套数据内部**按三分位分为低 / 中 / 高三组。
   - 敏感性分组：以 RPE1 的三分位切点作为固定阈值，套用到 GWPS 上，使两套数据的分组在绝对度数上可比。这一分组只用于 P-c。
3. 不做其他分组。符号匹配沿用富集分析的 `symbol_mapping` 规则。匹配不上的基因单列计数，不归入任何组。

## 4. 要报告的量（全部报告，不挑选）

对每套数据 × 每种分组 × 每个组：
- n（基因数），以及 T5 占比；
- ρ(ĥ, E_true)，95% CI；
- AUROC(T5)，95% CI；
- 次要量：M1 相对 M0 的 T5 覆盖增量 E1（随机宽度匹配，`an.paired_T5`）及 95% CI。

另报告每套数据的效应强度网络自相关 **ρ_net**：Spearman(log E_true(g), 该基因 STRING 邻居 log E_true 的加权平均)。邻居取 E ∪ always-train 中 score ≥ 0.7 的基因，权重为 STRING score。95% CI 由 TSS 组聚类 bootstrap 给出。

**组间对比**：每种分组取两个预先指定的对比，分别是复合体成员 − 非成员，以及度数高三分位 − 低三分位。对每个对比报告 Δρ 和 ΔAUROC 及 bootstrap CI。

**CI 的算法**：ρ 和 AUROC 的 CI 沿用 `an.predictability` 的做法，即按基因聚合、按 TSS 组聚类 bootstrap，重复次数与之相同；Δ 的 CI 用同一套 bootstrap 重抽样。随机种子用 `master_seed + 9100`。

## 5. 多重比较

主对比共 8 个：2 套数据 × 2 种分组 × (Δρ, ΔAUROC)。对它们做 Holm 校正，报告校正后 p 值。ρ_net 是 2 个独立估计，另报告两者的差及 CI。E1 分组结果只作描述，不做检验。

## 6. 解读规则（预先写定）

- P-a 和 P-b 都成立时，写作"与模块化假设相符（探索性）"。
- 任一不成立时，如实写"不支持"。
- 无论结果如何，都不改变 GWPS 上"第 1 步未通过"的确认性结论，也不据此调整方法再回到确认集上检验。

## 7. 资源

只读缓存和注释文件，不需要大内存（预计 < 1.5 GB），用时 < 30 分钟。代码会写成一个独立脚本，与冻结代码分开，并在冻结清单中登记为确认后探索。


## 8. 补充（用户要求，2026-10-15，运行前写入）

### 8a. E_true 的分半可靠性与衰减校正
- 对每个数据集（开发集、RPE1、GWPS）的每个测试单元，用 A 半细胞计算 E_A，用 B 半细胞计算 E_B。参照量 μ/σ 用该划分的参照，靶基因列的排除规则与 E_true 相同。E_A 直接取缓存中的 `EA_calt`；E_B 用同一函数（`proxy_v06.e_true`）在 B 半 pseudobulk 上计算。运行时校验：用同一方式重算 E_A 和全量 E，必须与缓存逐位相等。
- 每个基因取其各次测试出现的平均值，聚合方式与 `gene_proxy` 相同。
- 报告以下量：
  - r_AB = Spearman(E_A, E_B)，同时报告 log 尺度上的 Pearson 相关作为补充；
  - Spearman–Brown 可靠性 rel = 2r/(1+r)；
  - 衰减校正值 ρ_corr = ρ(ĥ, E_true) / √rel。只校正 E_true 一侧，ĥ 视为确定的预测值。
- 上述量在全集和第 3 节的每个分组内都报告，95% CI 来自同一套 TSS 组聚类 bootstrap。
- **局限（预先说明）**：每一半只有约一半细胞，而且 top-50 基因是在各半内部分别选出的，所以单半的 E 存在噪声底抬高。Spearman 的衰减校正只是近似。rel < 0.2 时 ρ_corr 不稳定，这种情况下只报告，不解读。

### 8b. 分半 oracle
- 沿用冻结代码中已经算好的 `variant = "half"` 结果：oracle 用 E_A 归一化分数并分配宽度，覆盖用基于 B 半响应的分数 S1 判定，T5 按 E_A 分层。把它与全量 oracle（宽度和覆盖判定用的是同一批细胞）对比。两者都在名义 α 下，模型为 linear-PCA，k = 50。
- 对 M0 和 oracle、full 和 half 的组合，报告 T5 覆盖及 95% CI（`analysis_v06.gee_tiers`，按划分和靶基因排序）。三个数据集都报告。这一步不重新计算任何覆盖。

### 8c. 范围限制
- 每个分组同时报告 E_true 的 SD 和 IQR（原尺度和 log 尺度），以及 T5 占比。
- 解读组间 ρ 差异时，必须结合离散程度说明：离散度小的组，ρ 天然偏低。

### 8d. 预先写定的解读规则（主线表述）
GWPS 全集上：
- **规则 N（噪声解释）**：rel 的 95% CI 上界 < 0.6，并且 ρ_corr 的点估计 ≥ 0.30（P+ 的 ρ 门槛）。满足时，"先验缺乏信息"的解释不成立，主线改为"非 essential 基因的效应强度测量噪声限制了可预测性"。
- **规则 I（信息解释）**：rel 的 95% CI 下界 ≥ 0.6，或者 ρ_corr 的 95% CI 上界 < 0.20。满足时，保留"先验缺乏信息"的解释，同时说明噪声对可预测性的限制程度（以 ρ_corr 为准）。
- 两条规则都不满足时，写作"测量噪声与先验信息不足的贡献无法区分"。
- 无论落在哪条规则下，都不修改任何确认性结论。

## 9. 实现细节（运行前固定）
- **脚本**：`code/explore_modularity_v07.py`，执行前登记 sha256。
- **ĥ 与 E_true**：取自各数据集的 `proxy_rows.parquet`（model = linear_pca，kind = h，proxy = sim_prior|ridge）。运行时校验：用 `an.predictability` 算出的 ρ 必须与确认结果中记录的 P.rho 相等。
- **CORUM 成员身份**：使用 enrichment v0.7 的注释（GO 2026-05-19 的 GAF 同义词，经 `map_symbols` 映射；CORUM 文件 `corum_human_current.txt`）。
- **STRING 度数**：全基因组 v12.0，max(experimental, database) ≥ 0.7。
- **基因映射**：不在 STRING info 文件中的靶基因记为未映射；CORUM 映射不上的靶基因同样记为未映射。两类都单独计数，不归入任何组。
- **三分位规则**：q1、q2 取度数分布的 1/3 和 2/3 分位数，低组为 deg ≤ q1，中组为 q1 < deg ≤ q2，高组为 deg > q2。固定切点版本使用 RPE1 的 q1 和 q2。
- **bootstrap**：B = 2000，按 TSS 组聚类重抽样，抽样方式与 `cluster_boot` 相同。种子为 master_seed + 9100 + {开发集 0，RPE1 10，GWPS 20}；ρ_net 的种子在此基础上再 +1。各分组与组间对比共用同一次重抽样。跨数据集差值取两个数据集各自第 i 次重抽样之差（两数据集的 bootstrap 相互独立）。bootstrap p 值 = 2 × min(P(Δ* ≤ 0), P(Δ* ≥ 0))。
- **ρ_net 的邻居**：E 单元与 always-train 单元中，在该数据集 STRING 宇宙（`string_expdb_universe`）内 score ≥ 0.7 的基因。E 单元取基因层面的平均 E_true；always-train 单元取各划分 E_at 的平均。
- **分组 E1**：在随机宽度匹配行中（RPE1、GWPS 用 `confirm_rand_rows.parquet`，开发集用 `sens_rand_match_rows_linear_pca.parquet`），取 M1（γ = 0）与 M0，用 `an.paired_T5` 计算，只作描述。
