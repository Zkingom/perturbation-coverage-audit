# 审计实验方案：以扰动为交换单元的 split conformal 在结构化 OOD 下的覆盖率

版本 v0.2（2026-10-01）。v0.1 经用户批准，v0.2 的修改见“v0.2 修订”一节；正文中与修订冲突之处，以修订为准。所有参数在下文“预注册参数”中固定，运行前如需改动，必须以新版本号记录改动及理由。

## v0.2 修订（相对 v0.1）

1. **新增决策关口 G1（位于 P1 之后、任何 GPU 工作之前）。** 只运行均值和线性模型：审计 A 100 次划分，审计 B 20 次划分。按第 8 节的判据给出初步结论，去掉涉及 GEARS 的条件，即支持条件改为“线性模型在最远分箱欠覆盖、斜率为负、两种距离方向一致、线性减均值的覆盖差为负、审计 B 边际欠覆盖”。根据 G1 的结论再决定是否进入 GPU 阶段。
2. **P0/P1 改用 Replogle K562 essential 数据集**，即 GEARS 0.1.2 `PertData` 提供的版本（Harvard Dataverse datafile 7458695，记录 SHA256）。GWPS 推迟到主线确立之后。P0 报告该数据集中 E、E_noGO 的规模，并据此重新计算功效。
3. **推迟到第二阶段：** Norman、GEARS-unc、E_strong。第一阶段的敏感性分析只做 k ∈ {20, 200} 和 D_STRING。
4. **S1 的定位：** S1 的 top-k 基因集由真实响应决定，在部署时无法计算，所以 S1 只用于审计，用来量化现有预测器在 OOD 下的误差分布。到方法阶段，需要的是预测时可计算的集合，例如由训练数据或预测本身决定基因集的分数（如 S2，或用预测 top-k 定义的分数）。
5. **运行前的资源核查：** 确认本机可用的硬盘和内存，并估算租用单张 A10/A100 跑 GEARS 试运行（P2）的费用，见第 10 节附注。
6. **K562 essential 版本的数据说明。** GEARS 版本是已经 log-normalize 的单细胞矩阵，且基因数已被 GEARS 预先筛选，所以第 1 节的表达表示改为“直接用 GEARS 矩阵按靶基因取细胞均值”，读出基因集 R = GEARS 矩阵中的全部基因。z-normalized 版本的敏感性分析暂缓。

## 0. 对原始设定的修改

我对原设定做了五处修改或补充，每一处都附有理由；不同意的话，可以在运行前改回。

1. **非一致性分数的参照点从“对照细胞”改为“训练扰动的质心”。** 以对照为参照选 top-k 差异基因，结果会被所有扰动共有的系统性偏移主导。Systema 把每个扰动的偏移与平均扰动效应之间的余弦相似度定义为系统性变异，并建议把参照点换成所有扰动质心的平均（perturbed centroid）（[Viñas Torné et al. 2025](https://doi.org/10.1038/s41587-025-02777-8)）。本方案中，这个质心只用训练扰动计算。
2. **主分数给出的是半径，不是逐基因区间；另设一个逐基因区间分数。** 如果 top-k 基因由真实响应决定，conformal 集合是“在该扰动自身 top-k 特异基因上 RMSE ≤ q̂”的球，宽度就是 q̂。要得到逐基因区间，需要一个由训练数据预先固定的基因面板，见第 4 节的 S2。
3. **随机划分只能检验条件覆盖，所以增设结构化划分。** 随机划分下校准集和测试集可交换，边际覆盖率由定理保证；“远 OOD 欠覆盖”在随机划分里只能表现为按距离分箱的条件覆盖失败，而且必然由近端分箱的过覆盖来抵消。真实的部署场景是测试扰动系统性地比校准扰动更远，需要另做 GO 簇整体留出的划分（审计 B）来检验。建议审计 A（随机划分）作为主分析，审计 B 作为必要的补充。
4. **所有模型使用同一评估总体。** GEARS v0.1.2 会静默删除不在其 GO 图中的扰动（源码 `pertdata.py` 第 184–194 行）。而这些扰动恰恰是先验信息最少、最远的 OOD。线性 baseline 也只能预测靶基因属于读出基因的扰动（[Ahlmann-Eltze et al. 2025](https://doi.org/10.1038/s41592-025-02772-6)）。因此主分析只在两者的交集上进行；不在 GO 中的扰动单独成为一个“无先验”分箱，只评估均值和线性模型。
5. **把均值 baseline 作为距离趋势的阴性对照。** 均值 baseline 的预测不依赖距离，它在各分箱上的覆盖差异只反映效应大小随距离的分布变化。例如注释少的基因可能效应较弱，这样远端分箱反而会过覆盖。所以模型本身外推造成的欠覆盖，要看“模型覆盖率减去均值覆盖率”的趋势。

## 1. 数据与版本

### 主数据集：Replogle 2022 K562 全基因组 CRISPRi（GWPS）

数据来源用原始 Figshare+ 存档 [10.25452/figshare.plus.20029387](https://doi.org/10.25452/figshare.plus.20029387)，不用 scPerturb 的二次处理版本。理由是 GWPS 的标准分析单元和 z-normalization 都由原作者定义，二次处理会多出一层不透明步骤。该存档每个群体提供四个文件：原始单细胞、原始 pseudobulk、按 gemgroup z-normalize 的单细胞、按 gemgroup z-normalize 的 pseudobulk，均只保留 >0.01 UMI/cell 的基因。下载时记录文件名、Figshare 文件 ID、字节数、MD5/SHA256 和下载日期。scPerturb（[Peidli et al. 2024](https://doi.org/10.1038/s41592-023-02144-y)）只用于交叉核对扰动数和细胞数；如果使用，记录其 Zenodo 记录号和版本号（本次我没能在线核实 scPerturb 当前的 Zenodo 版本号，下载时必须记录）。

表达表示（主分析）：从 `K562_gwps_raw_singlecell_01.h5ad` 以 backed 模式分块读取，每个细胞做 library-size 归一化到 1e4 后取 log1p，再按扰动对细胞求均值得到 pseudobulk y_p。对照 y_ctrl 是全部非靶向对照细胞的均值。选这个表示，是因为 GEARS 的输入和输出都在这个尺度上，三个模型必须预测同一个量。z-normalized pseudobulk 只用于均值和线性模型的敏感性分析。

文件中扰动按 “gene_transcript” 标注，同一基因可能有多个 TSS 行。分析单元是**靶基因**。主分析中每个基因保留细胞数最多的一行；敏感性分析保留所有行，并对靶基因做聚类 bootstrap。

### 次要数据集：Norman 2019 CRISPRa

使用 GEARS v0.1.2 `PertData` 提供的 `norman` 版本（Harvard Dataverse），以便和 GEARS 文献可比，并记录文件 MD5。单基因与双基因扰动的确切数量以下载文件为准，不引用二手数字。原始论文为 [Norman et al. 2019](https://doi.org/10.1126/science.aax4438)。

### 先验信息与软件

- GO：使用 GEARS 自带的 `gene2go_all.pkl`（Dataverse datafile 6153417），确保距离和 GEARS 使用同一份先验，并记录 MD5。该文件的 GO 版本日期未知，这是一个需要处理的问题（见 L9）。
- STRING（第二种距离）：v12.0，只使用 experiments 和 databases 两个通道，排除 co-expression 和 text-mining。
- GEARS：`cell-gears==0.1.2`（PyPI 当前最新版；Ahlmann-Eltze et al. 2025 使用的也是 0.1.2）。我已下载源码并检查了与泄漏相关的代码路径，见第 6 节。
- 环境：`conda env export` 和 `pip freeze` 存档；记录 torch、torch_geometric、scanpy、anndata、numpy、scipy 的版本，分析代码的 git commit hash，以及 GPU 型号、驱动和 CUDA 版本。

## 2. 评估总体与划分

### 总体定义（只用实验前可知的信息或实验设计层面的 QC）

- **E（主总体）**：靶基因 ∈ 读出基因集 R（线性模型的要求，也意味着靶基因在 K562 中有表达），且 ∈ GEARS 的 `pert_names`（有 GO 注释），且细胞数 ≥ 30。细胞数是测序设计层面的 QC，不是响应结局。
- **E_noGO**：靶基因 ∈ R，但不在 GO 图中。只评估均值和线性模型，作为“无先验”分箱。
- **E_strong（条件分析）**：E 中被 Replogle 2022 判定为有显著转录表型的扰动（使用原文补充表中的能量距离检验结果）。这是按结局定义的子集，校准集和测试集必须同样限制在此子集内。由此得到的保证不能用于部署，因为未见扰动是否“有表型”事先未知。论文中必须明确标注这一点。

R 取文件中的全部基因（>0.01 UMI/cell）。这个阈值由原作者在所有细胞上计算，不随我们的划分变化，属于轻微的全局信息，需要在方法部分写明。

### 划分单元与比例

单元是靶基因。为防止 CRISPRi 在双向启动子上对相邻基因的抑制造成“隐性重复”，先按 GENCODE（与原作者参考基因组注释一致的版本，下载时记录）把 TSS 相距 ≤1 kb 的基因合并为一个划分组；敏感性分析中不做合并。

审计 A（随机）：每次划分以划分组为单位随机分配，训练 60% / 验证 10% / 校准 15% / 测试 15%。验证集只用于 GEARS 的 early stopping，线性和均值模型不用它；验证集永远不能与校准集共用。

审计 B（结构化）：在 GO Jaccard 图上做 Leiden 聚类（只用先验，分辨率预先固定），每次随机抽取若干整簇作为测试集，使测试集约占 15%；校准集则从其余簇中随机抽 15%。这模拟了“在已知基因上校准、部署到功能上陌生的基因”的场景。

重复次数：三个模型在同一组划分上配对比较，审计 A 和审计 B 各 20 次。均值和线性模型代价很低，审计 A 额外跑到 100 次，用来刻画划分本身带来的变异。

### Norman 的划分

单基因按靶基因随机划分（同样是 60/10/15/15）。双基因按两个组分中有几个属于训练基因，标记为 seen2、seen1、seen0，校准按 seen 类别分层（Mondrian），因为只有在类别内部，校准和测试才可交换。Norman 的结果只做描述性报告：每个分箱只有个位数到十几个扰动，置信区间会很宽。

## 3. 模型

1. **均值 baseline**：ŷ_p = b，即训练扰动 pseudobulk 的平均。在扰动特异空间中，它的预测恒为 0。
2. **线性 baseline（严格照 Ahlmann-Eltze et al. 2025 的设定）**：Y_train 是“基因 × 训练扰动”矩阵，b 为其行均值。对 Y_train 做 PCA，取前 K=10 个主成分作为 G。P 取 G 中对应训练扰动靶基因的那些行。W = (GᵀG+λI)⁻¹Gᵀ(Y_train−b)P(PᵀP+λI)⁻¹，λ=0.1。预测为 Ŷ = G W P̃ᵀ + b，其中 P̃ 是 G 中对应测试扰动靶基因的行。每次划分重新拟合 PCA，且只用训练扰动。Norman 双基因扰动使用原文的 additive baseline。
3. **GEARS v0.1.2**：使用默认超参数（hidden_size=64，num_similar_genes_go_graph=20，num_similar_genes_co_express_graph=20，coexpress_threshold=0.4，epochs=20），uncertainty=False。不在校准集或测试集上调参；如需调参，只在 split 0 的训练/验证集上进行，然后冻结。GWPS 训练时，每个训练扰动最多抽 32 个细胞（抽样种子固定），再加固定数量的对照细胞；评估时使用 GEARS 对每个扰动的平均预测。
   - 可选分支 GEARS-unc：uncertainty=True，把 GEARS 自带的 log-variance 用作自适应分数的归一化因子，检验模型自带的不确定性能否修复远端欠覆盖。这个分支会使 GPU 开销翻倍，放在第二阶段。

## 4. 非一致性分数

记号：Δ_p = y_p − y_ctrl。训练质心 μ_tr = mean_{q∈train} Δ_q。扰动特异信号 s_p = Δ_p − μ_tr，模型预测 ŝ_p = Δ̂_p − μ_tr。基因尺度 σ_g = 1.4826·MAD_{q∈train}(s_qg)，下限截断在训练集 σ 的 5% 分位数。μ_tr 和 σ_g 每次划分只用训练扰动计算。

**S1（主分数，top-k 特异基因 RMSE）**：T_k(p) 是按 |s_pg| 排序的前 k 个基因，排除靶基因本身，k=50。

S1(p) = sqrt( (1/k) Σ_{g∈T_k(p)} ((s_pg − ŝ_pg)/σ_g)² )

集合为 C(p) = { y : S1(p; y) ≤ q̂ }，宽度报告为半径 q̂，单位是归一化的。选这个分数有三个理由：(a) 误差只在该扰动真正特异响应的基因上计算，回应了 Systema 的批评；(b) 所有模型共用同一个分数函数，q̂ 可以直接比较，q̂ 越小，在同等保证下信息越多；(c) 对均值 baseline 来说 S1 就是特异效应的大小，它构成一个可解释的零模型。已知的缺点是：弱效应扰动的 top-k 会选到噪声最大值（赢家诅咒）。但这对所有模型相同，并且会在 E_strong 中做敏感性分析。

**S2（次分数，固定面板上的逐基因区间）**：面板 H 是 s 在训练扰动间方差最大的 500 个基因（只用训练数据）。S2(p) = max_{g∈H} |s_pg − ŝ_pg|/σ_g，对应逐基因区间 ŝ_pg ± q̂·σ_g，并且在整个面板上同时覆盖。宽度报告为面板平均的 2q̂σ_g。

**敏感性**：k ∈ {20, 200}；把靶基因纳入 T_k；以对照而非 μ_tr 为参照（用来展示参照点选择的影响）。

conformal 分位数取校准分数的第 ⌈(n+1)(1−α)⌉ 个顺序统计量（[Lei et al. 2018](https://doi.org/10.1080/01621459.2017.1307116)；[Angelopoulos & Bates 2021](https://arxiv.org/abs/2107.07511)）。主分析 α=0.1，次要分析 α ∈ {0.05, 0.2}。

## 5. OOD 距离与分箱（只用先验）

- **D_GO（主距离）**：使用 GEARS 同一份 gene2go 文件和同一个 Jaccard 定义。D_GO(p) = 1 − 测试扰动 p 与训练集基因之间最大 20 个 Jaccard 相似度的均值。20 与 GEARS 默认的 GO 图邻居数一致，因此该距离直接衡量 GEARS 能从训练集“借到”多少信息。敏感性分析用 1 − 最大 Jaccard。
- **D_STRING（第二距离）**：1 − 与训练集基因的最大 STRING 分数（v12.0，仅 experiments 和 databases 通道）。
- **分箱**：每次划分，对 E 中所有校准和测试单元计算距离，按五分位数分为 Q1–Q5。如果 D=1（与任何训练基因都无 GO 重叠）的单元超过 10%，就把它们单列为 Q_iso，其余单元按四分位数分箱。另有 E_noGO 分箱。分箱边界只由距离决定，不看任何表达数据。

需要指出一点循环性：D_GO 用的正是 GEARS 的输入先验。所以“GEARS 在 D_GO 远端欠覆盖”的解释是“GEARS 缺少可借用的邻居”，而不是一般意义上的生物学距离。D_STRING 用来检查结论是否依赖于具体选用的先验。

## 6. 泄漏检查点

每个检查点都要写成代码里的 assert，并输出日志，不能只做人工确认。

| 编号 | 步骤 | 风险 | 检查 / 措施 |
|---|---|---|---|
| L1 | 划分 | 同一靶基因的不同 guide 或 TSS 行出现在两侧 | assert：每个靶基因的所有行只出现在一个 fold |
| L2 | 划分 | 双向启动子邻近基因被 CRISPRi 同时抑制 | TSS ≤1 kb 的基因合并为一个划分组；敏感性分析不合并 |
| L3 | 预处理 | 用全部扰动选高变基因 | R 由原作者固定；面板 H、μ_tr、σ_g、PCA 只用训练扰动；assert 计算时传入的索引 ⊆ train |
| L4 | 线性模型 | PCA 包含校准或测试列 | assert：Y_train 的列 ⊆ train |
| L5 | GEARS 共表达图 | 图由训练条件的细胞计算（源码 `utils.py::get_coexpression_network_from_train`），但会以 `{split}_{seed}_{train_gene_set_size}_{threshold}_{k}` 为文件名**缓存**；自定义划分如果复用同一名字，会读入其他划分的图 | 每次划分使用唯一的 split 名，并且每次运行前清空 `data_path` 缓存；记录读入图文件的哈希 |
| L6 | GEARS GO 图 | `make_GO` 缓存到 `./data/go_essential_{data_name}.csv` | 按划分命名 data_name，或确认 GO 图与划分无关（它只依赖先验，所以可以共享，但必须核对哈希） |
| L7 | GEARS 模型选择 | 按验证集 top-20 DE MSE 选最佳 epoch（源码 `gears.py` 第 570–572 行），其 DE 基因由 `rank_genes_groups_cov_all` 预先计算 | 验证集只来自训练侧基因；不使用 GEARS 内部的 DE 注释做任何评估；我们的 T_k 独立计算 |
| L8 | GEARS 过滤 | 不在 GO 中的扰动被静默删除，导致评估总体不一致 | 在进入 GEARS 之前显式构造 E；assert：GEARS 输出的扰动集合 = E ∩ test |
| L9 | 先验的时间泄漏 | 2022 年 6 月后的 GO 注释可能部分来自 Replogle 2022 本身的结果 | 敏感性分析：用 2022-06 之前的 GO 归档版本重算 D_GO 和分箱（需先确认能获取归档版本）；记录 gene2go_all.pkl 的哈希 |
| L10 | 校准集 | 用校准集调参或做 early stopping | 校准集只在训练结束后读入一次；在代码层面把 cal/test 索引写入只读文件 |
| L11 | 结局定义的总体 | E_strong 用到了测试结局 | 只作为条件分析报告，校准和测试同样限定，不作为主结论 |
| L12 | Norman | 双基因扰动的组分已出现在训练中 | 按 seen 类别分层校准并分开报告 |

## 7. 统计分析

每次划分 s、每个模型 m、每个分箱 b，计算覆盖率 c_smb 和半径 q̂_sm。对非自适应分数，q̂ 对所有分箱相同，这是结构性的；分箱间的差异体现在覆盖率上。另外报告“达到名义覆盖所需的半径” r_smb，即该分箱内测试分数的 (1−α) 分位数，用来量化远端区间需要放宽多少。

报告三种区间：

1. **跨划分的分布**：20 次划分的均值，以及 2.5%–97.5% 范围。这刻画单次部署时可能遇到的覆盖率变异。理论参照是：给定 n_cal 时，split conformal 的覆盖率服从 Beta(n+1−l, l)，其中 l = ⌊(n+1)α⌋。观测到的边际覆盖分布应落在这个参照带内，否则提示流程有 bug。
2. **均值的推断区间**：同一基因会在多次划分中反复进入测试集，各次划分不独立，直接对跨划分均值做 t 区间会偏乐观。主推断采用二项 GLMM：miscover ~ model × bin + log(n_cells) + (1|target gene) + (1|split)，从中得到各 (model, bin) 的覆盖率及 95% CI，以及距离趋势的斜率。
3. **配对差异**：模型 − 均值 baseline 在同一分箱上的覆盖差，用同样的 GLMM 对比来估计。

多重比较：3 个模型 × 2 种距离 × 主 α，趋势检验采用 Holm 校正。

功效粗估：E 中约有数千个基因（确切数字在 P0 阶段统计）。若取 n_cal ≈ n_test ≈ 1,000，则每个分箱约 200 个测试单元，单次划分的二项 SE ≈ √(0.9·0.1/200) ≈ 2.1 个百分点。预设的最小有意义偏差 δ = 5 个百分点，GLMM 汇总 20 次划分后，应有足够功效检测出来。真实的功效以 P0 统计出的实际数量重新计算。

## 8. 支持或否定假设的判据（预注册）

假设 H：在 D 的远端，perturbation-level split conformal 欠覆盖。

**支持 H**，需要同时满足：

1. 审计 A 中，GEARS 和线性模型至少一个在最远分箱（Q5 或 Q_iso）的覆盖率 < 1−α−δ，且 GLMM 95% CI 的上界 < 1−α；
2. 覆盖率随 D 的斜率为负，Holm 校正后显著；
3. 在 D_GO 和 D_STRING 两种距离下方向一致；
4. 在同一分箱上，“模型 − 均值 baseline”的覆盖差为负，表明欠覆盖来自模型外推，而不只是效应大小的分布；
5. 审计 B 的边际覆盖率低于 1−α−δ。

**否定 H**：所有分箱的覆盖率都在 1−α ± δ 以内，CI 半宽 < δ（说明功效足够），斜率不显著，审计 B 的边际覆盖也接近名义值。

**需要重新解释的结果**：

- 均值 baseline 也呈现同样的距离趋势，而模型相对均值没有差异：说明远端扰动本身的效应分布不同，这是异方差或条件覆盖问题，不是模型外推问题。论文的叙事要相应改变。
- 远端**过覆盖**、近端欠覆盖：这同样是条件校准失败，但方向与 H 相反，H 被否定，结论需要如实改写。
- 趋势在加入 log(n_cells) 后消失：说明是测序深度效应，不是 OOD 效应。
- 只在 E_strong 中出现：保证仅在结局定义的子集上失败，部署意义有限。
- 只在 Norman 中出现、GWPS 中不出现：样本量太小，不作为结论。

## 9. 对照与流程自检

- **阳性对照（植入欠覆盖）**：在测试集预测上人为加入与 D 成正比的特异噪声，确认分析流程能检测出预设大小的覆盖下降。
- **阴性对照（打乱距离）**：在测试单元之间随机置换 D，趋势应当消失。
- **边际覆盖自检**：审计 A 中每个模型的边际覆盖率应落在 Beta 参照带内。

## 10. 算力与运行时间估计

当前本机为 Windows，20 核，16 GB 内存（目前可用约 1 GB）；沙箱无法访问 GPU，也没有配置远程计算资源。

| 任务 | 资源 | 估计耗时 | 说明 |
|---|---|---|---|
| 下载 GWPS 原始单细胞 h5ad 和 pseudobulk | 磁盘：数十 GB（以实际文件大小为准，下载前确认） | 取决于网速 | 必要时需申请 figshare 域名的网络权限 |
| 构建 log-normalized pseudobulk | 本机，backed 模式分块读取，内存 <4 GB | 1–3 小时（I/O 受限） | 一次性完成，结果作为 checkpoint 保存 |
| GO / STRING 距离 | 本机 | 数分钟 | 每次划分重新计算到训练集的距离 |
| 均值 + 线性，100 次划分 | 本机，多核 | <1 小时 | K=10 的截断 SVD，每次划分只需数秒 |
| conformal + GLMM + bootstrap | 本机 | <1 小时 | GLMM 约 10⁵ 行 |
| GEARS，GWPS，20 次划分 | GPU（≥16 GB 显存），主机内存 32–64 GB | 每次约 2–6 GPU 小时，合计约 40–120 GPU 小时 | **估计不确定，可能偏差 2–3 倍**。依据是一个非同行评审的二手来源：Norman（约 11 万细胞）在 16 GB GPU 上训练 20 epoch 约 1–2 小时。实际数字必须在 split 0 的试运行中测量 |
| GEARS，Norman，20 次划分 | 同上 | 合计约 20–40 GPU 小时 | 同上 |
| GEARS-unc（可选） | 同上 | 再翻倍 | 放在第二阶段 |

结论：均值和线性部分可以在本机完成；GEARS 部分在本机不可行，需要在 Customize → Compute 中添加带 GPU 的远程主机（单张 A10、A100 或同级别即可）。如果 GPU 预算不足，可以把 GEARS 的划分数减到 10 次，但这是对“至少 20 次”的偏离，必须在论文中声明，不能静默缩减。

### 第 10 节附注（v0.2）：运行前资源核查与 GPU 费用估算

- 本机（2026-10-01 核查）：内存 15.7 GB，P1 运行时可用 1.5–2.2 GB；可用硬盘约 33.7 GB（已扣除 K562 essential、STRING 等约 6 GB 数据）。K562 essential 的 P0/P1 在本机分块处理即可完成。GWPS 原始单细胞文件预计远大于剩余空间，主线确立后要换到外部存储或远程主机上处理。
- GPU 租用报价（2026 年公开价目，各平台随时变动，下单前需重新确认）：A100 80GB 按需价约 1.19–1.39 美元/小时（RunPod community/secure cloud），Lambda 约 2.79 美元/小时；Lambda 的 A10 约 1.29 美元/小时。
- P2（GEARS，K562 essential，split 0，两个训练种子）：每次运行估计 1–3 GPU 小时，另加 1–2 小时环境配置与调试，合计约 3–8 小时，即约 **4–22 美元**。实例的主机内存应 ≥32 GB，因为 GEARS 会为每个细胞构建 PyG 图对象。
- 若进入 P3（GEARS，20 次划分）：约 20–60 GPU 小时，即约 25–170 美元。以上均是外推估计，以 P2 实测为准。

## 11. 分阶段执行

1. **P0 数据审计（本机；v0.2 改为 K562 essential）**：下载数据、记录哈希；统计 E、E_noGO、E_strong 的规模，GO 覆盖率，每个扰动的细胞数分布，以及双向启动子分组数；据此重新计算功效。
2. **P1 流程验证（本机）**：均值和线性模型跑 3 次划分，完成三项自检（第 9 节）；所有泄漏 assert 必须通过。
   - **G1 决策关口（v0.2 新增）**：均值和线性模型在审计 A 上跑 100 次划分、在审计 B 上跑 20 次划分，按修订后的第 8 节判据给出初步结论，再决定是否进入 P2。
3. **P2 GEARS 试运行（GPU）**：只跑 split 0，测量时间和内存，并用两个训练种子重复，量化训练随机性。
4. **P3 全量运行**：审计 A 和审计 B 各 20 次划分，三个模型。
5. **P4 分析与报告**：GLMM、敏感性分析，并按第 8 节的判据给出结论。

## 12. 预注册参数（机器可读）

```yaml
protocol_version: "0.2"
master_seed: 20261001
split_seeds: "master_seed + i, i = 0..99"   # GEARS uses i = 0..19
gears_train_seeds: [0, 1]                   # second seed only on split 0 (P2)
data:
  gwps_source: "Figshare+ 10.25452/figshare.plus.20029387"
  gwps_files: ["K562_gwps_raw_singlecell_01.h5ad", "K562_gwps_normalized_bulk_01.h5ad"]
  record: [file_id, bytes, sha256, download_date]
  norman_source: "GEARS 0.1.2 PertData('norman') via Harvard Dataverse"
  expression: "per-cell normalize_total(1e4) + log1p, mean per target gene"
population:
  min_cells: 30
  require_target_in_readout: true
  require_target_in_gears_go: true
  bidirectional_promoter_merge_bp: 1000
split:
  unit: "target-gene group"
  fractions: {train: 0.60, val: 0.10, cal: 0.15, test: 0.15}
  n_splits_all_models: 20
  n_splits_mean_linear_extra: 100
  structured_split: "Leiden clusters on GO Jaccard graph, whole-cluster test holdout"
models:
  mean: {}
  linear: {K: 10, ridge_lambda: 0.1, source: "Ahlmann-Eltze et al. 2025 Methods"}
  gears: {version: "0.1.2", hidden_size: 64, num_similar_genes_go_graph: 20,
          num_similar_genes_co_express_graph: 20, coexpress_threshold: 0.4,
          epochs: 20, uncertainty: false, max_cells_per_pert_train: 32}
score:
  reference: "train perturbed centroid (Systema)"
  S1: {k: 50, exclude_target: true, gene_scale: "1.4826*MAD over train perts, floor q05"}
  S2: {panel_size: 500, panel_rule: "top variance of s over train perts", aggregate: "max"}
  sensitivity_k: [20, 200]
conformal:
  alpha_primary: 0.10
  alpha_secondary: [0.05, 0.20]
  quantile: "ceil((n+1)(1-alpha))-th order statistic"
distance:
  primary: "D_GO = 1 - mean(top-20 Jaccard to train genes), GEARS gene2go_all.pkl"
  secondary: "D_STRING = 1 - max STRING v12.0 score (experiments+databases only)"
  bins: "quintiles over cal+test; separate Q_iso if >=10% have D=1; separate noGO bin"
analysis:
  delta_pp: 5
  glmm: "miscover ~ model*bin + log(n_cells) + (1|gene) + (1|split)"
  multiplicity: "Holm over models x distances (primary alpha)"
```

## 13. 预期产出

- 每次划分 × 模型 × 分箱的结果表（parquet）：覆盖率、q̂、r_b、n。
- 主图：覆盖率随距离分箱的变化（每个模型一条线，带 GLMM 95% CI，名义值用虚线标出）；所需半径 r_b 随距离分箱的变化。
- 补充材料：审计 B 的边际覆盖率、全部敏感性分析、三项自检结果、完整的版本与哈希清单。
