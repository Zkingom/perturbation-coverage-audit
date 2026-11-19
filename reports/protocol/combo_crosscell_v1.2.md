# 方案 v1.2：组合扰动（Norman 2019）与新细胞类型（Nadig 2025 HepG2 / Jurkat）上的 R 诊断量（冻结版）

**冻结时间**：见 `manifests/freeze_v1.2.json`（UTC 与北京时间）。冻结时本机没有任何 Norman、HepG2、Jurkat 数据文件（`G:\biology\scperturb_13350497`、`G:\biology\geo_GSE264667` 均不存在，manifest 中有记录），因此本方案写定于读取这些数据的任何内容（包括 obs 元数据）之前。

**证据等级**：确认结果观察之后的新增分析（post-confirmation），但所用数据集此前从未被读取过。按本方案事先写定的预测，可以称为"在未接触数据上事先登记的检验"，不能称为 v0.7 意义上的确认性结果。论文中统一表述为 *pre-registered post-confirmation analysis on previously unseen data*。

**与已冻结方案的关系**：v0.7（审计）、gears_coupling v1.0/v1.1 的全部内容不变。本方案复用其冻结的函数（`models_v04.train_reference`、`gene_pca_embedding`、`fit_linear`、`score_S1`；`pseudobulk_v06` 的单细胞变换；`lowmem_v07` 的两遍读取；`splits_v06.repeated_group_kfold`），只新增数据接入、加性基线、组合扰动的线性模型嵌入（第 3.3 节，**需用户确认的算法扩展**）和本方案的判定规则。

---

## 1. 问题与预测

R = 1 − R²(log S1 ~ log E_true)，即分数中不能由效应强度解释的份额（定义同 Methods 与 gears_coupling v1.0）。在三个 CRISPRi 筛选上，linear‑PCA 的 R 为 0.112（开发集）、0.065（RPE1）、0.0006（非文库靶基因集）；R 很小时，按效应强度分层的覆盖率主要由构造决定（构造性参照约 0.5）。

**问题 1（Norman，组合扰动）**：在"未见组合"任务上，加性基线（两个单基因效应相加）掌握了组合效应的主体，分数主要剩下非加性的相互作用残差。我们预期这会让分数与效应尺度的耦合变弱，即 R 明显升高。若成立，这就是 R 在真实数据中离开 0 的第一个例子，此时分层覆盖率开始携带模型相关的信息。

**问题 2（HepG2 / Jurkat，新细胞类型）**：在同一平台、未接触过的两个细胞系上，linear‑PCA 的 R 是否仍然很小（预测：不升高）。

两个问题都按第 6 节的规则判定，**无论结果方向如何都如实报告，不改规则**。

## 2. 数据

| 数据集 | 文件 | 来源 | 校验 |
|---|---|---|---|
| Norman 2019 (K562, CRISPRa) | `NormanWeissman2019_filtered.h5ad`，698,680,199 字节 | scPerturb v1.4，Zenodo 记录 13350497 | MD5 `c870e6967d91c017d9da827bab183cd6` |
| Nadig 2025 HepG2 (CRISPRi, essential) | `GSE264667_hepg2_raw_singlecell_01.h5ad`，5,614,460,941 字节 | GEO GSE264667 | GEO 未公布 MD5：核对字节数与 h5ad 结构，记录 SHA‑256 |
| Nadig 2025 Jurkat (CRISPRi, essential) | `GSE264667_jurkat_raw_singlecell_01.h5ad`，9,366,490,264 字节 | GEO GSE264667 | 同上 |

文件由用户手动下载（速度低于 3 MB/s 或文件大于 2 GB 时的约定）。读取一律采用 backed 模式，分块计算 pseudobulk（`lowmem_v07` 的两遍读取，块大小按可用内存设定，进程内存上限 4 GB）。全程断点续跑，每步写完成标记（`resume_v07` 的原子写 + sha256 标记）。

### 2.1 C0 检查（只读元数据和计数，不计算任何分数）

1. 文件核对（第 2 节）。确认 X 为非负整数计数：抽查前 2,000 个细胞，所有非零值都是整数。若不是原始计数，**停下报告**（不能用现有变换）。
2. 单细胞变换与 v0.6 相同：x = log1p(1e4 · counts / UMI_count)，float64 计算，按扰动取均值得到 pseudobulk；同时按固定种子把每个扰动的细胞随机分成 A、B 两半，分别得到半数 pseudobulk（用于第 5.3 节阳性对照）。分半种子 = MASTER_SEED + 15000 + 数据集编号（Norman 0、HepG2 1、Jurkat 2）。
3. 读出基因：
   - Norman：对照细胞中平均原始计数 ≥ 0.1 的基因，**再并入所有靶基因**（文件中存在的）。阈值只用对照细胞计算，不用任何扰动细胞。并入靶基因是因为 CRISPRa 的靶基因在 K562 对照中往往低表达，若不并入，嵌入就无法覆盖它们；而靶基因列本身在分数中始终被排除（3.5），所以并入不会改变任何 S1 的计算基因。
   - Nadig：沿用文件自带的基因集，与 Replogle 文件处理方式一致。
4. 单元表、细胞数，以及每个单元的靶基因列（靶基因在读出基因中时，排除对应的列）。

## 3. 单元、划分与模型

### 3.1 Norman 的单元

- 扰动标签按 scPerturb 的 `perturbation` 列解析，组合统一写成按字母排序的无序基因对，若同一基因对以两种顺序出现，合并其细胞。
- **单基因单元（single）**：细胞数 ≥ 30。
- **E 单元（评估单元）= 双基因单元**，需同时满足：细胞数 ≥ 30；两个组成基因都有单基因单元；两个组成基因都在读出基因中（按 2.1 第 3 条，文件中存在的靶基因都已并入读出基因；linear‑PCA 嵌入需要，见 3.3）。
- **不满足的双基因单元**：只要细胞数 ≥ 30，仍可进入训练池（作为拟合和参照的一部分），但不参与评估。
- 对照细胞：标签为 `control` 的细胞，用于 y_ctrl。
- **事先规定的样本量下限**：若 E 单元少于 60 个，Norman 分析照常执行并报告，但 R 的判定标注为"样本量不足（underpowered）"，不进入摘要。

### 3.2 Norman 的划分："未见组合"

- 组 = 无序基因对。在 E 单元上做 3 轮重复分组 7 折（`splits_v06.repeated_group_kfold`，种子基数 MASTER_SEED + 14000），共 21 次划分。
- 每次划分中：测试 = 当折的 E 单元；其余 E 单元按组分配到验证（10%）、校准（15%）和训练。
- **所有单基因单元和非 E 的双基因单元始终在训练池中**，测试组合的两个组成基因都作为单基因扰动出现在训练中。这是"未见组合、已见单基因"（GEARS 中的 seen 2/2）任务的标准定义，不属于泄漏。
- **泄漏断言**（运行时）：
  - 测试、校准、验证组合不出现在训练池中；
  - 参照量 μ/σ 只用训练池计算；
  - 加性基线只读训练池中的单基因 pseudobulk 和对照细胞。

### 3.3 模型（Norman）

1. **均值模型**（空白对照）：预测 μ_tr。其分数按定义等于 E_true，所以 R ≡ 0。
2. **加性基线（主模型）**：ŷ_AB = y_A + y_B − y_ctrl（对数归一化 pseudobulk 空间，即两个单基因效应之和；与 Norman 2019 和 Ahlmann‑Eltze 2025 的加性基线一致）。
3. **linear‑PCA（次模型；需用户确认的算法扩展）**：
   - G = 训练表达矩阵（基因 × 训练扰动，含对照列）的前 K 个主成分（`gene_pca_embedding`，method 与 v0.7 相同）。
   - 单基因扰动 A 的嵌入为 P_A = G[A]；**双基因扰动 AB 的嵌入定义为 P_AB = G[A] + G[B]**。
   - 拟合 W（`fit_linear`，λ = 0.1），预测 ŷ_AB = G W P_ABᵀ + b。
   - K 在验证折上从 {5, 10, 20, 30, 50} 中选取（准则与 v0.7 相同：验证集 S1 均值最小）。
   - 原论文的线性模型只针对单基因扰动定义，"嵌入相加"是本方案为组合扰动所做的扩展。若用户不同意，在读取任何数据之前以 v1.2.1 修订（删除该模型或替换定义）。

### 3.4 HepG2 / Jurkat

- 单元、E 单元和划分完全沿用 v0.7 第 3 节（C0）：
  - TSS 分组规则不变（合并距离 1,000 bp）；
  - 每单元细胞数 ≥ 30；
  - 3 轮 × 7 折，共 21 次划分；
  - 种子基数 MASTER_SEED + 16000 + 数据集编号。
- GEO 没有提供 bulk 文件，单元元数据由单细胞 obs 聚合得到。若 TSS 分组所需的注释无法复现，改为按靶基因分组，并作为偏离记录。
- 模型：均值模型与 linear‑PCA（与 v0.7 相同，K 在验证折上选取）。不运行 linear‑GO、GEARS、scGPT 或代理学习器。

### 3.5 分数、效应强度、分层与区间

- 参照：μ_tr、σ_g 用训练池中细胞数 ≥ 30 的扰动计算（`train_reference`）。
- 分数：S1 = `score_S1(y, ŷ, ref, tcols, k=50)`。
- E_true = 均值模型的 S1。靶基因列：所有在读出基因中的靶基因都排除，双基因单元排除两列。
- 校准：M0 全局 split‑conformal，α = 0.1，阈值 q̂ = 校准分数的第 ⌈(n+1)(1−α)⌉ 小值。
- 分层：每次划分内，按 E_true 在该划分全部 E 单元上的五分位分为 T1–T5（规则同 `eval_v06.tiers_of`）。
- 区间宽度：M0 下所有单元的半径都是 q̂，按定义在各层之间相同。因此另外报告"相对宽度" q̂ / E_true 的各层中位数，它反映区间相对于效应大小有多宽，属描述性量。

## 4. 计算的量（每个数据集、每个模型；21 次划分的测试行合并计算）

1. R = 1 − R²(log S1 ~ log E_true)，普通最小二乘。
2. CV = CV(S1 / E_true)。
3. RER = 1 − mean(S1) / mean(E_true)，即相对均值模型的误差下降，整体和 T5 内各算一次。
4. M0 的边际覆盖率，以及 T1–T5 各层的覆盖率：截距模型 GEE，logit 连接，独立工作相关，聚类见第 4 节末尾。
5. 各层的 q̂ 均值与相对宽度中位数。
6. ΔT5 = T5 覆盖率(模型) − T5 覆盖率(均值模型)：在同一测试行上配对计算。均值模型的 T5 覆盖率就是构造性参照在本数据上的实现。

95% 置信区间：

- **聚类**：Norman 以无序基因对为聚类（每对 3 行），HepG2/Jurkat 以 TSS 组为聚类。
- **自助法**：聚类重抽样 2,000 次，种子 MASTER_SEED + 17000。所有模型在同一次重抽样中一起计算，配对差也在同一次重抽样中算出。
- **已知局限**：Norman 中共享同一基因的组合彼此相关，以基因对为聚类不能完全刻画这种相关，CI 可能偏窄。

## 5. 对照

### 5.1 阴性对照 1：打乱效应强度（用户 A3）

- 每次划分内，把 E_true 在该划分的 E 单元之间随机错排（derangement），然后重新分层。
- 用打乱后的 E 重算 R、CV、各层覆盖率，S1 不变。
- 共 20 次置换，种子 (MASTER_SEED + 13000, split_id, r)。
- **通过条件（事先写定）**：对每个模型，20 次置换的平均值同时满足：
  - (a) 覆盖率层间差 T5 − T1 落在 [−0.10, 0.10] 内；
  - (b) |R_m − R_mean| < 0.05，即各模型之间的 R 差异消失（打乱后两者都接近 1）。

  任一条件不满足，说明分层或 R 的计算流程本身产生了伪信号，**所有 R 判定不予采信**，如实报告。

### 5.2 阴性对照 2：打乱预测（同 v1.1 §2c）

- 每次划分内，把模型对测试单元的预测在测试单元之间随机错排；top‑50 基因仍按真值选取。
- 对加性基线和 linear‑PCA 各做 20 次置换，种子 (MASTER_SEED + 12500, split_id, r)。
- **通过条件**：20 次中被第 6 节规则判为"R 明显升高"的次数为 0。若不为 0，说明规则不能拒绝"由噪声冒充的解耦"，**"R 明显升高"的判定不予采信**。

### 5.3 阳性对照：分半预测器（同 v1.1 §2，半数据口径）

- 用同一单元 A 半的 pseudobulk 作为"预测"，B 半作为真值，E_B = 均值模型在 B 半上的 S1；所有模型同时在半数据口径下重算。
- **作用**：检验第 6 节的阈值在该数据上能否达到。若达不到，在结果中如实说明；它不参与主判定。

## 6. 判定规则（运行前写定，不得修改）

常数：

- R_lo = 0.112：三个筛选上 linear‑PCA R 的最大值，取自 `reports/review/constructive_reference_by_model_v2.csv`，开发集 0.11221，取三位小数。
- R_hi = 2 × R_lo = 0.224。

对每个数据集的每个非均值模型 m：

- **准确度前提**：RER_m（整体）的 95% CI 下界 > 0。
- **"R 明显升高"**：准确度前提成立，且 R_m 点估计 ≥ R_hi，且 R_m 的 95% CI 下界 > R_lo。
- **"R 升高但准确度未改善"**：R_m 点估计 ≥ R_hi，但准确度前提不成立。这种情况不解释为耦合变弱。
- **"R 未升高"**：R_m 的 95% CI 上界 < R_hi。
- **其余情况**："不确定"。

**主检验**：

- Norman 上的加性基线。linear‑PCA 为次要结果。
- HepG2 与 Jurkat 上的 linear‑PCA，事先预测为"R 未升高"。

**附带读数（描述性，不作判定）**：

- 当判定为"R 明显升高"时，ΔT5 的点估计与 CI，用于说明分层覆盖率是否偏离了构造性参照。
- 当判定为"R 未升高"时，T5 覆盖率是否与均值模型的 T5 覆盖率一致。

## 7. 预写表述

- "明显升高"：*On combinatorial perturbations, the additive baseline left a share R = ⟨⟩ (95% CI ⟨⟩) of the score unexplained by effect strength, above the pre-specified threshold of 0.224 (twice the largest value on the three single-perturbation screens), with an error reduction of ⟨⟩ relative to the mean model. Stratified coverage on this dataset therefore carries model-dependent information: T5 coverage differed from the constructive reference by ⟨⟩ pp (95% CI ⟨⟩).*
- "未升高"：*R = ⟨⟩ (95% CI ⟨⟩), below the pre-specified threshold; as on the single-perturbation screens, stratified coverage is close to its constructive value.*
- "不确定"或"样本量不足"：如实写出 CI，不作方向性陈述。
- HepG2 / Jurkat 的预测成立："In two cell lines not previously analysed, R stayed below the pre-specified threshold (HepG2 ⟨⟩, Jurkat ⟨⟩)." 若不成立，照实写，并把它作为对"R ≈ 0 是筛选普遍性质"这一说法的反例写入正文。

## 8. 不做的事

- 不运行 GEARS、scGPT 或任何深度模型。
- 不调整 k、α、阈值或分层方法。
- 不在看到结果后增加模型或子集。
- 任何事后追加的分析都单独标注为 post hoc。

## 9. 执行顺序与断点

C0（2.1）→ 划分写出（只读）→ 模型与分数（每次划分写缓存和完成标记）→ 对照 5.1–5.3 → 汇总与判定。

断点记录写在 `results/v1.2/<dataset>/resume_state.json`。距离断电不足 3 小时时，不启动 pseudobulk 的完整一遍读取。

## 10. 修订记录

- v1.2（本版）：初版冻结。
