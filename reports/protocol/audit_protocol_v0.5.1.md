# 审计实验方案：以扰动为交换单元的 split conformal 在结构化 OOD 下的覆盖率

版本 v0.5.1（2026-10-04）。修订历史依次见“v0.5.1 修订”“v0.5 修订”“v0.4.2 修订”“v0.4.1 修订”“v0.4 修订”“v0.3 修订”“v0.2 修订”各节；与正文冲突时，以最新一次修订为准。第 12 节的 YAML 是 v0.5 的完整机器可读设定。所有参数在下文“预注册参数”中固定，运行前如需改动，必须以新版本号记录改动及理由。

## v0.5.1 修订（相对 v0.5；2026-10-04T07:55:54Z，G1 计算进行中、查看任何覆盖结果之前）

**唯一改动：C5 从“支持 H”的必要条件降为次要结果。**

- **新的“支持 H”规则**：至少一个线性模型（linear-PCA 或 linear-GO，选出的版本）同时满足 C1、C2、C3，并且 C4 成立。
- **C5 的处理**：C5（审计 B 的 GEE 边际覆盖率 < 0.85）照常计算并报告。满足时作为补充证据；不满足时，按 v0.5 第 1 节的预先声明解读，即“在这个偏移幅度下未检测到影响”，不能据此否定 H。
- **其他一律不变**：“否定 H”的条件、C1–C4 的定义、Holm 族、R1–R5、全部模型、数据和计算代码都保持原样。

**理由。** v0.5 第 1 节已经预先声明：审计 B 相对审计 A 的偏移很小（D_GO 中位数约高 0.04），检验力有限。在这种情况下还把 C5 设为“支持 H”的必要条件，与这条预先声明自相矛盾。这处修改由用户提出。

**时间顺序（如实记录）。**
1. v0.5 于 2026-10-04T06:23:27Z 冻结。
2. G1 计算从约 06:30Z 开始，由子任务执行冻结代码。
3. 用户在 G1 计算期间提出本次修改。此时审计 A 的 100 次划分已经算完，审计 B 正在计算，GLMM、GEE 和分类都还没有运行。
4. 到写下本节时，主分析者只看过 G1 进度日志里的测试集比例、选出的 K 和 d、每次划分的耗时，没有看过任何 G1 的覆盖率、斜率或分类结果。

**代码。** 新增 `code/g1_analysis_v051.py`，与 `g1_analysis_v05.py` 只差分类规则那几行，输出写到 `g1_results_v0.5.1.json`。v0.5 的冻结文件都不修改。子任务按 v0.5 原计划运行冻结脚本得到的 `g1_results_v0.5.json` 原样保留，报告中同时列出两种分类，主结论以 v0.5.1 为准。

## v0.5 修订（相对 v0.4.2，2026-10-04；G1 运行版本，冻结）

本节在看到任何“覆盖率随距离变化”的结果之前写成。冻结后不再修改任何参数；如需改动，另立新版本号。

### 1. 审计 B 的设计（用户 2026-10-04 确认）

- **图**：节点为 E 的 1,786 个单元，基于原始 GEARS gene2go，不删任何 GO 项。每个基因保留 Jaccard 最大的 20 个邻居，只保留 J > 0 的边；对称化时取较大值，边权就是 Jaccard。
- **聚类**：Leiden，`RBConfigurationVertexPartition`，加权，分辨率 8，种子 = master_seed（20261001）。代码在 `code/auditB_clusters_v05.py`，输出冻结为 `data/processed/v0.5/auditB_clusters_original_r8.csv`。
- **测试集停止规则 C**（`splits_v04.cluster_holdout_split_C`）：按随机顺序考察每个簇；把它连同它牵涉的 TSS 组一起加入后，若测试集离 15%·n 更近或一样近，才加入；测试集达到 15%·n 时停止。余下单元按 TSS 组随机分到验证集（10%·n）、校准集（15%·n），其余进训练集。
- **测试集比例**：容许范围 13.5%–16.5%，可行率门槛 95%。G1 的 20 次划分都照常使用，不会因为测试集比例而重抽或丢弃。
- **冻结版本与设计扫描之间的一处差别。** 冻结代码选邻居时用 numpy 的稳定排序（`kind="stable"`），以保证并列值的处理可复现。设计扫描用的是默认的 quicksort，所以并列值的取舍不同：边数是 28,186 对 28,160，两次划分的 ARI 为 0.858。冻结划分的特征如下：
  - 69 个簇，最大簇 69 个单元（占 E 的 3.9%）；
  - 规则 C 下，200 个种子中 98.5% 的测试集比例落在容许范围内；
  - G1 实际使用的 20 个种子，测试集比例在 0.149–0.165 之间，全部落在容许范围内。

  这些特征仍满足用户确认的两条标准。
- **对偏移幅度的预先声明（在看任何结果之前写下）。** 设计扫描显示，审计 B 测试集到训练池的 D_GO 中位数只比审计 A 高约 0.04（0.807 对 0.770），而 E 的 D_GO 四分位距约为 0.66–0.81。因此：
  - 若审计 B 与审计 A 的覆盖率相近，只能解读为“在这个偏移幅度下未检测到影响”，**不能**解读为“结构化 OOD 不影响覆盖”；
  - G1 报告中会同时给出审计 B 实际的偏移量（D_GO、D_GO_filt、D_STRING）。

### 2. G1 统计分析（全部预先设定）

- **主数据**：S1(k = 50)，α = 0.1，3 个主模型：均值、linear-PCA（选出的 K）、linear-GO（选出的 d）。rank(D) 是该次划分 cal ∪ test 内的百分位秩；只有测试单元进入模型。
- **主检验（GLMM，`code/g1_glmm_v05.R`，lme4）**：miscover ~ model × rank(D) + log(n_cells)（中心化）+ (1|靶基因) + (1|划分)，二项分布，bobyqa 优化器。分别对审计 A、B 和三种距离（D_GO、D_STRING、D_GO_filt）各拟合一次；另对审计 A、B 在 D_GO 上各拟合一次不含 log(n_cells) 的版本。
  - **注意尺度。** GLMM 斜率是条件（基因特异）尺度。同一基因在不同模型、不同划分中的 miscover 高度相关，所以基因随机效应的方差大，GLMM 斜率的数值会明显大于群体平均（GEE）的斜率。合成数据检验证实了这一点。因此判据只用 GLMM 斜率的**符号和 p 值**，覆盖率的**幅度**一律以 GEE 为准。
- **Holm 族**：审计 A，3 个主模型 × {D_GO, D_STRING}，共 6 个双侧 Wald p 值。审计 B 的趋势、D_GO_filt 都是次要结果，不进入 Holm 族。
- **GEE（`code/g1_analysis_v05.py`）**：logistic，独立工作相关，按 TSS 组聚类，稳健标准误。它报告各模型的斜率、在 rank(D) = 0.1 和 0.9 处的群体平均覆盖率（log(n_cells) 取均值）及 95% CI，以及审计 A、B 的边际覆盖率。
- **敏感性分析（GEE，审计 A，rank(D_GO)）**：k = 20、k = 200、S2、α = 0.05、α = 0.2；复现版模型（K = 10、d = 2）分别在 D_GO、D_STRING、D_GO_filt 上。
- **审计 A 的边际覆盖自检**：规则同 P1 自检 (i)（二项 P(≥ n_out) > 0.01）。若不通过，G1 结果标记为“流程可能有误，不作解读”。

### 3. 第 8 节判据的 G1 版本（不含 GEARS；δ = 5 个百分点）

对每个线性模型 m（linear-PCA、linear-GO，均为选出的版本）：

- **C1（远端欠覆盖）**：审计 A，D_GO；GEE 在 rank = 0.9 处的覆盖率 < 0.85，且 95% CI 上界 < 0.90。
- **C2（趋势）**：审计 A，D_GO；GLMM 的 miscover 斜率 > 0，且 Holm 校正后 p < 0.05。
- **C3（方向一致）**：审计 A；GLMM 斜率在 D_GO 和 D_STRING 下符号相同。
- **C4（模型外推，暂定）**：审计 A，D_GO；GLMM 中 linear-GO 与均值的斜率差 > 0，且 p < 0.05。不做 Holm 校正，GEARS 阶段再确认。另外报告 linear-PCA 与均值的斜率差。
- **C5（审计 B 边际欠覆盖）**：GEE 估计的审计 B 边际覆盖率 < 0.85。

**分类规则：**

- **支持 H**：至少一个线性模型同时满足 C1、C2、C3、C5，并且 C4 成立。
- **否定 H**：同时满足以下三条：
  - 3 个主模型在 rank = 0.1 和 0.9 处的 GEE 覆盖率都在 [0.85, 0.95] 内，且 CI 半宽 < 0.05；
  - 3 个主模型的审计 B 边际覆盖率都在 [0.85, 0.95] 内；
  - Holm 族中没有任何一个显著的斜率（无论正负）。
- **其他情况**：部分支持或无法判定，逐项报告。

**需要重新解释的情形：**

- **R1**：均值模型也有 Holm 显著的正斜率（D_GO），而 C4 不成立。这说明问题在于远端扰动本身的效应分布（异方差或条件覆盖），而不是模型外推。
- **R2**：Holm 族中出现显著的**负**斜率，即远端过覆盖。
- **R3**：某个模型在不含 log(n_cells) 时斜率为正且 p < 0.05，加入 log(n_cells) 后 p ≥ 0.05。这说明趋势来自测序深度。
- **R4**（只在 E_strong 中出现）、**R5**（只在 Norman 中出现）：本阶段不评估，因为这两部分已推迟。

### 4. 代码与验证

- **v0.5 的代码**：
  - `code/pipeline_v05.py`：P1 与 G1 的运行脚本。由 `p1_v04.py` 派生；模型、分数和参照量的代码不变。新增审计 B 划分、rank(D) 列和逐划分输出；P1 不再生成按距离分箱的覆盖表。
  - `code/splits_v04.py`：规则 C 划分函数。
  - `code/auditB_clusters_v05.py`：审计 B 聚类。
  - `code/g1_analysis_v05.py`：GEE、自检和判据。
  - `code/g1_glmm_v05.R`：GLMM。
- **合成数据端到端检验**（表达用随机数，不含任何真实数据；输出保存在 `results/_synthetic_tests_v0.5/`）：
  - **零情形**（审计 A 3 次划分、审计 B 2 次划分）：全部 assert、canary 和 Gram/SVD 检查通过；边际自检通过；Holm 族没有显著斜率；分类为“否定 H”。
  - **植入情形**：只把 linear-GO 的 miscover 改为随 rank(D_GO) 上升（真实的 logistic 斜率为 3）。分类正确判为“支持 H”，均值和 linear-PCA 都没有被误判。普通 logistic 回归估计的斜率为 3.27；GLMM 的条件斜率为 5.6，原因见上面的“注意尺度”。审计 A 的边际自检也正确地报告了失败。

### 5. 冻结

- 冻结时间、方案文件和全部代码与输入文件的 SHA256，记录在 `manifests/freeze_v0.5.json`。
- 方案存档为 `reports/protocol/audit_protocol_v0.5.md`，设为只读。

## v0.4.2 修订（相对 v0.4.1，2026-10-03）

1. **bulk 原始文件恢复。** `K562_essential_raw_bulk_01.h5ad` 曾从 `<DL>` 中缺失，用户已手动重新下载。我独立核对了文件大小、MD5 和 SHA256，三者都与 2026-10-01/02 的记录一致，已重新登记并设为只读。冒烟测试第 1 项再次通过。

2. **调参对称化。**
   - linear-GO 的 d 网格扩为 {2, 5, 10, 20, 30, 50, 100}。GO 嵌入改为一次算出 100 维（`go_ase100.npz`），前 20 列与 `go_ase20.npz` 一致，最大差 5.1e-14。
   - linear-PCA 的 K 也在验证集上选，网格为 {5, 10, 20, 30, 50}，规则与 d 相同：比较验证集上 S1(k=50) 的平均误差，并列时取较小值，选定后不在验证集上重新拟合。
   - PCA 只算一次 K = 50 的解，各候选 K 取其前 K 列。主成分得分是嵌套的，所以这与分别按各个 K 计算完全相同。
   - linear-GO 的基因嵌入 G 仍固定为 PCA 的前 10 列。
   - K = 10 和 d = 2 作为复现版一并报告。选中的值若再次落在网格上界，如实报告，不再扩大网格。
   - **时间顺序**：这项改动是在看到 P1（v0.4.1）的验证误差之后、看到任何“覆盖率随距离变化”的结果之前做出的。P1 中已经看过的覆盖相关输出只有两项：自检 (i) 的边际覆盖，和各模型的 q̂。`descriptive_coverage_by_DGO_bin.csv` 生成后没有打开，也没有解读。
   - 代码中，Gram 与 SVD 的一致性比较扩展到 K = 50 的 G、全部候选模型的预测，以及选中的 (K, d)；canary A 要求选中的 (K, d) 不变。
   - 评估模型改为 5 个：均值、linear-PCA（选中的 K）、linear-PCA K = 10 复现版、linear-GO（选中的 d）、linear-GO d = 2 复现版。主检验与 Holm 校正的范围仍是均值、linear-PCA（选中）、linear-GO（选中）3 个模型 × 2 种距离；复现版只作敏感性分析。
   - 代码已在合成表达数据（随机数，不含任何真实表达）上完整跑通：全部 assert、canary 和自检通过。真实数据上尚未重跑。

3. **审计 B 的聚类设计（待用户确认后固定）。** 设计扫描的结果见 `results/auditB_design_v0.4.2/` 和 `reports/auditB_design_v0.4.2.md`，全部只用先验信息。另外发现，v0.2 的 `audit_lib.cluster_holdout_split` 从未运行过，而且有缺陷：它传给 `random_group_split` 的比例里没有 test，各比例也不是相对全部单元定义的。现已改用 `code/splits_v04.py` 中的新实现。

4. **冗余文件清理。** 用户批准的冗余文件已移入回收站。有 4 个空目录无法通过回收站删除，暂时保留。

## v0.4.1 修订（相对 v0.4，2026-10-03，在 P1 运行之前）

1. **第 12 节 YAML 已改写**，与 v0.4 的实际设定逐项一致：数据、表达表示与 L 的规则、三个模型（含 linear-GO 的 d 网格和选择规则）、参照量的计算样本、D_GO_filt、连续趋势 GLMM、Holm 范围。未启用的部分（GWPS、Norman、GEARS、GEARS-unc、E_strong、审计 B 的具体参数）标为 deferred。
2. **均值模型的预测改为 μ_tr**：
   - v0.4 中，均值模型的预测是全部训练扰动（训练池，含细胞数 < 30 的始终训练扰动）的平均，而 μ_tr 只用细胞数 ≥ 30 的训练扰动计算。两者不相等，于是均值模型在特异空间 s = y − μ_tr 中的预测不再恒为 0，与第 3 节第 1 条和第 4 节对 S1 的解释（“对均值 baseline，S1 就是特异效应的大小”）不一致。
   - v0.4.1 起，均值模型的预测与 μ_tr 使用同一组扰动：训练池中细胞数 ≥ 30 的那些。它的特异空间预测因此恒为 0，第 3、4 节的解释成立。代码中用一条 assert 检查这一点。
   - 线性模型照旧使用全部训练池扰动拟合（linear-PCA 要求靶基因在读出基因中，linear-GO 要求靶基因在 GO 图中）。

3. **PCA 计算方法的真实数据验证（P1，2026-10-03）**：3 次划分中，numpy SVD（gesdd）在真实数据上均收敛，不需要换用 gesvd。Gram 与 SVD 的比较结果：
   - G 按列对齐符号后的相对最大差：7.3e-14（阈值 1e-8）；
   - 全部线性模型候选在 val、cal、test 上预测的最大绝对差：1.2e-12（阈值 1e-8）；
   - 每次划分选出的 d 两种方法相同（均为 20）。
   三条标准都满足，因此采用 Gram，P1 的结果都用 Gram 计算。
   说明：当初改用 Gram，是因为 gesdd 在调试用的替代数据（bulk 文件先平均再取 log）上不收敛；在真实的“先取 log 再平均”数据上，这个问题没有出现。每次划分都会重复这项比较，结果记录在 `results/p1_v0.4/p1_v0.4_summary.json` 的 `pca_check` 字段。

## v0.4 修订（相对 v0.3）

### 0. 存储与范围

- 项目根目录由 `<ROOT>/config.toml` 定义，下载目录为 `download_root`（当前是 G:\biology）。原始数据只在 `<DL>` 保存一份，`data/raw/raw_manifest.json` 记录每个文件相对 `<DL>` 的路径和完整哈希。
- 凡是需要保留的中间文件和划分，一律直接写到 `<ROOT>`。
- v0.2 已归档，不再复现。其报告、结果和 manifest 只作为历史记录保留。

### (a) 评估总体

接受 E = 1,786（分属 1,721 个 TSS 组）。其余定义与 v0.3 相同：始终训练的扰动 257 个；与 E 同处一个 TSS 组的非 E 扰动 14 个，予以剔除。

### (b) 表达表示

- 单细胞文件为 `K562_essential_raw_singlecell_01.h5ad`（Figshare 文件 35773219，10,661,879,995 字节），下载到 `<DL>\figshare_20029387_v1\`，核对 Figshare 公布的 MD5 `4f1122ce1c7f13299a68df6459a266d3` 和文件大小后才可使用。
- 每个细胞的表达 = log1p(10⁴ · counts / L)。L 是文库大小，规则现在固定：若细胞元数据中有每细胞总 UMI 数（列名 `UMI_count`），就用它作为 L；没有则用这 8,563 个基因的计数之和。实际采用哪一种，写入运行 manifest。两种 L 的区别在于：前者是该细胞的全转录组文库大小，后者只覆盖原作者保留的基因，是一个已知的近似。
- 先对每个细胞取 log，再在每个扰动的**全部细胞**上求平均，得到 y_p；对照 y_ctrl 是全部 non-targeting 细胞的平均。整个计算分块进行。
- 评估用的真实响应，以及均值和线性模型的训练 pseudobulk，都使用每个扰动的全部细胞。
- **已知限制：** GEARS 训练时每个扰动最多只用 32 个细胞，这是出于内存考虑（见 v0.3 的 P0 资源估计）。因此 GEARS 的训练信号比均值和线性模型更嘈杂，三者比较时对 GEARS 略为不利。这一点要写进论文的局限部分。

### (c) linear-GO

- 版本核实（2026-10-02，Zenodo API）：10.5281/zenodo.14832393 是概念 DOI，即论文 Code availability 部分引用的那个 DOI（文献 33），它下面共有 5 个版本。其中 16092690 是第 5 个、也是最新的版本，创建于 2025-07-18，晚于接受日期 2025-06-24，早于在线发表日期 2025-08-04；在它之后没有更新的版本。它的代码包 MD5 为 e352192c…，与同日的 GitHub “Publication” 版本 16092611 字节一致。因此，16092690 就是论文发表时概念 DOI 解析到的版本。
- 复现依据（已发布代码，Zenodo 记录 16092690（论文所引 DOI 10.5281/zenodo.14832393，2026-10-01 访问时解析到此记录；下载并校验的正是该记录的文件，见 data/raw/raw_manifest.json），commit 8ea205c）：
  - `benchmark/src/extract_pert_embedding_from_gears.R` 第 31 行：`emb <- igraph::embed_adjacency_matrix(gr, no = 2, which = "lm", scaled = TRUE)`，嵌入维数固定为 2；
  - 同一文件第 4 行定义了 `--pca_dim`（默认 10），但脚本没有使用这个参数；
  - `benchmark/submission/run_perturbation_benchmark.R` 第 13 行传入 `gears_extract_pert_embedding(list(pca_dim = 10L))`，这个值没有起作用。
- 主版本 linear-GO 的嵌入维数 d 从 {2, 5, 10, 20} 中选，选择规则现在固定：
  1. 每次划分只用训练扰动拟合 4 个候选模型（G 为训练 PCA 的前 K = 10 维，λ = 0.1，W 只用靶基因在 GO 图中的训练扰动拟合）；
  2. 在验证集上计算每个候选的验证误差。误差定义为：对每个验证扰动，在它自身真实特异效应最大的 50 个基因上（排除靶基因），按 σ_g 缩放后计算预测误差的 RMS，即 S1(k=50)，μ_tr 和 σ_g 只用训练集计算；再对所有验证扰动取平均；
  3. 选平均误差最小的 d。若有并列（相对差 < 1e-12），取较小的 d；
  4. 选定后不在训练集加验证集上重新拟合。验证集只用于选 d，不进入训练，也不进入校准。
- 嵌入用 GO 图的 adjacency spectral embedding，一次算出 d = 20，取前 d 列作为各候选。前 2 列和前 10 列要与已存档的 `go_ase_embeddings.npz` 一致（容差 1e-8）。
- 复现版本固定 d = 2，与主版本一并报告。
- 我需要指出一处不对称：linear-GO 的 d 用验证集调过，而 linear-PCA 的 K 固定为 10，没有调。如果 linear-GO 优于 linear-PCA，其中有一部分可能来自调参。所以 P1 的验证误差中会同时报告“d = 2 复现版”，作为未调参的对照。

### (d) 参照量的计算样本

- μ_tr、σ_g 和面板 H，只用细胞数 ≥ 30 的训练扰动计算，包括始终训练集中细胞数 ≥ 30 的那部分。
- 模型拟合仍使用全部训练扰动。这意味着均值模型的预测 b（全部训练扰动的平均）和参照 μ_tr 不再相同。这样做是有意的：参照量要求低噪声，而模型应当利用全部可用的训练数据。

### (e) v0.4 的 P1 范围（待用户确认后执行）

- 模型：均值、linear-PCA、linear-GO（主版本和 d = 2 复现版），3 次划分（种子为 master_seed + 0、1、2）。
- 第 6 节的全部泄漏 assert，加上三项自检（边际覆盖、植入欠覆盖的阳性对照、置换距离的阴性对照），必须全部通过。
- 报告各模型在验证集上的预测误差，以及相对均值模型的改进：(均值误差 − 模型误差) / 均值误差，按划分分别列出。
- 划分文件写到 `splits/v0.4/`，写入后设为只读。

### (f) P1 的实现细节与通过标准（2026-10-02 固定，看到真实数据结果之前）

代码在 `code/p1_v04.py` 和 `code/models_v04.py`。在接入真实数据之前，用 bulk 文件构造了一份“先平均再取 log”的替代数据，只用来检查代码能否跑通，**结果不作任何报告**。

- **训练池** = E 中当次划分的训练集 + 257 个始终训练扰动。三个模型和三种先验距离（D_GO、D_GO_filt、D_STRING）都以这个训练池为参照。
- **线性模型** 照已发布代码 `solve_y_axb` 实现：ctrl 作为一列训练数据，其嵌入为零向量；先按拟合列的行均值把 Y 中心化，两侧 ridge 都取 λ = 0.1。
  - G 是“训练池 + ctrl”伪 bulk 矩阵（基因 × 列，按列中心化）的前 10 个主成分得分。
  - **PCA 计算方法（算法改动，暂定）**：v0.2 用的是 numpy SVD（LAPACK gesdd）。在替代数据上 gesdd 报“不收敛”，所以暂时改为对列 × 列 Gram 矩阵做特征分解，取 G = Xc·V。理论上它与 SVD 的 U·S 只差每列的正负号，而两种线性模型的预测都不受 G 的列符号影响。
  - **真实数据上的验证（看任何覆盖结果之前完成）**：每次划分都用 Gram 和 SVD 各算一次 G 及全部模型的预测。SVD 先用 numpy gesdd；若 gesdd 不收敛，再用 scipy 的 gesvd 驱动，并通知用户。判定标准如下，现在固定：
    1. G 按列对齐符号后，max|G_gram − G_svd| ≤ 1e-8 × max|G_svd|；
    2. linear-PCA 和 linear-GO 全部候选在 val、cal、test 上的预测，max 绝对差 ≤ 1e-8（单位是 log1p(CP10k)）；
    3. 每次划分选出的 d 完全相同。
  - 三条全部满足，才采用 Gram，并在报告中列出两种方法的最大差值；任何一条不满足，就改回 SVD（gesdd；不收敛时用 gesvd）。替代数据只用于调试，不计入这项验证。
- **canary A**：把 cal 和 test 的表达替换成随机数，要求参照量、所有候选模型的预测、三种距离以及选出的 d 全部不变。
- **canary B**：把 val 的表达替换成随机数，要求参照量、所有候选模型的预测和三种距离不变（d 可以变，因为 d 本来就由验证集决定）。
- **自检 (i) 边际覆盖**：在 3 次划分 × 4 个模型 × 4 种分数 × 3 个 α，共 144 个组合上，计算覆盖数落在 Beta-二项 95% 参照带之外的个数，记为 n_out。
  - 通过标准：在“每个组合以 5% 的概率落带外”的二项零假设下，单侧 P(≥ n_out) > 0.01。
  - 这些组合彼此相关，所以这只是一个近似标准。
  - 主组合（S1_k50，α = 0.1）是否落在带外，逐个列出。
- **自检 (ii) 阳性对照**：只在 linear-PCA 的测试集预测上植入噪声，强度为 τ·u·σ_g，其中 u 是测试单元 D_GO 的百分位，τ = 6。
  - τ = 6 沿用 v0.2 的设定；当时 τ = 2 不够大，τ = 6 是事后选的，这里如实标注。
  - 通过标准：在主组合上，覆盖失败对 rank(D_GO) 的趋势斜率 > 0，且 p < 0.05。检验用二项 GLM，按靶基因聚类的稳健标准误，3 次划分合并。
- **自检 (iii) 阴性对照**：在测试单元之间随机打乱 D_GO。通过标准：每个模型在主组合上的趋势检验都满足 p ≥ 0.05。
- **报告内容**：各模型在验证集上的 S1(k=50) 平均误差，以及相对均值模型的改进 1 − err_model/err_mean；每次划分选出的 d，以及 4 个候选的验证误差。若选出的 d 落在网格边界（d = 20），要注明。


## v0.3 修订（相对 v0.2）

### 1. 数据：改用 Replogle 原始 Figshare+ 的 K562 essential 文件

- 来源：Figshare+ 条目 20029387 的第 1 版（DOI 10.25452/figshare.plus.20029387.v1，最后修改 2023-05-30）。
- P0 使用 `K562_essential_raw_bulk_01.h5ad`（文件 ID 35773070）。MD5 `8321d5d3ffc99db2a5c71edca4189735`，与 Figshare 公布值一致；SHA256 记录在 manifest 中。
- 读出基因集 R：文件中的 8,563 个基因，即原作者筛选的 >0.01 UMI/cell 的基因。其中 2 个基因名重复，后续以 gene_id 作为主键。
- 每个靶基因保留细胞数最多的一行 TSS（共丢弃 119 行次要 TSS）。这些次要行**一律不进入任何 fold**，否则同一靶基因可能同时出现在训练集和测试集。
- 表达表示（P1 之前待定，见 P0 报告）：GEARS 的预测量是“每个细胞 log1p(CP10k) 后再取均值”。若所有模型都要预测同一个量，就需要下载 10.7 GB 的单细胞文件（ID 35773219，MD5 `4f1122ce1c7f13299a68df6459a266d3`），在本机分块计算。

### 2. 新增第三个模型 linear-GO（与 linear-PCA 并列）

- 依据：Ahlmann-Eltze et al. 2025 的 “LM with P from GEARS”。原文方法写的是：对 GEARS 的 GO 通路成员矩阵做 spectral embedding 得到 P。我查了其公开代码（Zenodo 记录 16092690（论文所引 DOI 10.5281/zenodo.14832393，2026-10-01 访问时解析到此记录；下载并校验的正是该记录的文件，见 data/raw/raw_manifest.json），`extract_pert_embedding_from_gears.R` 与 `run_linear_pretrained_model.R`，commit 8ea205c），实际做法如下：
  - 用 GEARS 的 `go_essential_all.csv`（Dataverse datafile 6934319；即两两 GO Jaccard > 0.1 的边表，共 9,853 个基因，边表对称）构造无向加权图；
  - 用 `igraph::embed_adjacency_matrix(no = 2, which = "lm", scaled = TRUE)` 做 adjacency spectral embedding（Sussman et al. 2012, JASA, 10.1080/01621459.2012.699795）。
- **代码与文字不一致：** 上游脚本通过参数传入 `pca_dim = 10`，但嵌入脚本没有使用这个参数，固定取了 no = 2。我只能确认已发布代码的行为，无法确认论文实际运行的是不是这个版本。
- 预注册：主版本 linear-GO 照已发布代码取 d = 2；敏感性版本取 d = 10。两者都在看任何覆盖结果之前固定。
- Python 复现：图取对称 max（等价于 igraph 的 `mode = "undirected"`）；去掉自环后，在对角线放 strength/(n−1)，对应 igraph `cvec` 的默认值；取模最大的 d 个特征值，X = U·|λ|^{1/2}。igraph 对自环的具体处理我没有逐行核实，与 R 版数值是否一致，计划在 P1 用 R igraph 交叉检查。
- 只用先验：嵌入只依赖 GO 边表，与划分无关，全数据只算一次（`go_ase_embeddings.npz`）。测试基因作为图中的节点参与嵌入，但这只用到先验注释，不涉及表达。
- 只用训练扰动：
  - G 是 Y_train（基因 × 训练扰动，含 ctrl 列）做列中心化 PCA 后的前 K = 10 个得分；
  - W = (GᵀG+λI)⁻¹ Gᵀ (Y_train − b) Pᵀ_train (P_train Pᵀ_train + λI)⁻¹，λ = 0.1，b 是 Y_train 的行均值；
  - 预测为 G W P_test + b。
  - W 的拟合只用靶基因在 GO 图中的训练扰动。
- 与已发布代码对齐的一处细节：已发布代码把 ctrl 作为一个训练列，其嵌入为零向量、change 为零。v0.3 中 linear-PCA 和 linear-GO 都照此处理（v0.2 的 P1 实现没有包含 ctrl 列，这是一处偏差，P1 重跑时修正）。

### 3. 非 E 扰动始终放入训练集

- 所有不属于 E、且不与 E 基因同处一个 TSS 合并组的扰动，固定放入每次划分的训练集，供所有模型使用。
- 进入各模型 W 拟合的条件：linear-PCA 要求靶基因 ∈ R；linear-GO 要求靶基因 ∈ GO 图。均值模型和以后的 GEARS 使用全部这些扰动。
- 与 E 基因同组的非 E 扰动直接剔除，不进入训练集，防止双向启动子造成的隐性泄漏。
- **需要声明的影响：**
  - μ_tr、σ_g、面板 H 和 PCA 都会包含这些扰动，其中有细胞数 <30、pseudobulk 噪声较大的扰动；
  - 所有 E 单元到训练集的 D_GO 和 D_STRING 会系统性变小，分箱边界也随之改变；
  - 结论只适用于“部署时训练集包含这批辅助扰动”的场景。

### 4. G1 主结论的范围

- G1 的主结论只回答一个问题：已部署预测器（均值、linear-PCA、linear-GO）的条件覆盖是否随先验距离下降。
- 判据 4（欠覆盖是否来自模型外推）改为用 linear-GO 与均值模型之间的覆盖差随距离的趋势做初步判断；到 GEARS 阶段再确认。

### 5. D_GO 的敏感性版本（在看到任何覆盖结果之前固定）

- **规则：** 计算 Jaccard 之前，删除在 GEARS GO 图基因全集（9,853 个基因）中注释基因数 >500 的 GO 项。
- **理由：**
  1. 500 是 GSEA 软件对基因集大小上限的默认值，常被用来排除过于宽泛、信息量低的基因集。这个数值来自软件默认设置，方法本身见 Subramanian et al. 2005, PNAS, 10.1073/pnas.0506580102；
  2. 在 9,853 个基因中，500 约占 5%。
- 词频要在这个去冗余的全集上统计，因为 `gene2go_all.pkl` 的 67,832 个键里包含大量别名（例如 A1BG、A1B、ABG、GAB 都映射到同一组注释），按全部键计数会把词频放大数倍。
- **只看先验的统计（不涉及任何表达或覆盖结果）：**
  - 被删除的 GO 项共 29 个，包括 GO:0005515 protein binding（7,628 个基因）、nucleus、cytosol 等；
  - E 基因的注释平均减少 32.4%，1.2% 的 E 基因删除后不剩任何 GO 项；
  - E 基因两两之间 Jaccard > 0.1 的比例从 19.9% 降到 1.4%；
  - 以全部 E 基因为池（留一）计算，D_GO 的 IQR 从 0.67–0.81 变为 0.77–0.90，与原始 D_GO 的 Spearman 相关为 0.77。
- **解读上的注意：** GEARS 自身的 GO 图使用的是未过滤的 Jaccard，所以过滤后的 D_GO 衡量的是更具体的功能距离，不再是“GEARS 能借到多少邻居信息”。

### 6. 主检验改为连续趋势检验

- 主检验：miscover ~ model × rank(D) + log(n_cells) + (1|target gene) + (1|split)。
  - 二项 GLMM；按基因的随机截距承担聚类，同一基因在不同划分中的重复观测因此不被当作独立样本。
  - rank(D) 是 D 在每次划分的 cal ∪ test 内的百分位秩。
  - 主要参数是各模型的 rank(D) 斜率，以及 linear-GO 与均值的斜率差。
- 多重比较：Holm 校正，范围是 3 个模型 × 2 种距离 = 6 个检验。
- 分箱覆盖率只作描述性报告。
- 稳健性：另用 GEE（按 TSS 组聚类，独立工作相关）加 cluster-robust SE 重复估计一次。

### 7. 已知的分析约束

E 单元按 TSS 组划分。在 v0.3 的数据下，E 有 1,786 个单元、1,721 个 TSS 组。每次划分的 train/val/cal/test 约为 1,074 / 178 / 267 / 267，另加 257 个始终在训练集的扰动。

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
protocol_version: "0.5.1"
master_seed: 20261001
split_seeds: "master_seed + i"          # P1: i = 0..2; G1 audit A: i = 0..99; audit B: i = 0..19
data:
  dataset: "Replogle et al. 2022, K562 essential CRISPRi (Figshare+ 10.25452/figshare.plus.20029387.v1)"
  bulk_file: {name: "K562_essential_raw_bulk_01.h5ad", figshare_file_id: 35773070, bytes: 79766954,
              md5: "8321d5d3ffc99db2a5c71edca4189735", use: "unit table (targets, num_cells_filtered, readout genes)"}
  singlecell_file: {name: "K562_essential_raw_singlecell_01.h5ad", figshare_file_id: 35773219, bytes: 10661879995,
                    md5: "4f1122ce1c7f13299a68df6459a266d3", use: "expression representation"}
  readout_genes: 8563                    # original authors' filter (> 0.01 UMI/cell)
  hashes: "data/raw/raw_manifest.json (sha256 of every raw file)"
expression:
  per_cell: "log1p(1e4 * counts / L)"
  library_size_L: "obs.UMI_count if present (whole-transcriptome), else sum over the 8,563 retained genes"
  per_perturbation: "mean over ALL cells of the perturbation (truth and mean/linear training)"
  control: "mean over all non-targeting cells"
  gears_training_cells: "<= 32 per perturbation (known limitation; GEARS deferred)"
population:
  unit: "one TSS row per target gene (max num_cells_filtered); other TSS rows never used"
  E: {require_target_in_readout: true, require_target_in_gears_go: true, min_cells: 30, n: 1786, tss_groups: 1721}
  always_train: {rule: "non-E and not sharing a TSS group with an E gene", n: 257}
  dropped_tss_conflict: 14
  tss_merge_bp: 1000
  coordinates: "Ensembl release 116, GRCh38.p14 (frozen in p0v3_units_all_targets.csv)"
split:
  unit: "TSS group"
  fractions: {train: 0.60, val: 0.10, cal: 0.15, test: 0.15}
  train_pool: "E-train + always_train"
  folds_written: "splits/v0.4/*.json, read-only"
reference:
  rows: "training-pool perturbations with >= 30 cells"
  mu_tr: "mean of y over reference rows (Systema perturbed centroid)"
  sigma_g: "1.4826 * MAD of s over reference rows, floor at 5% quantile"
  panel_H: "500 genes with largest variance of s over reference rows"
models:
  mean: {prediction: "mu_tr (reference rows only); specific-space prediction == 0"}
  linear_pca: {source: "Ahlmann-Eltze et al. 2025, released code Zenodo 16092690 (concept DOI 10.5281/zenodo.14832393), solve_y_axb",
               G: "top-K PCA scores of (training pool + ctrl) pseudobulk, genes x columns, column-centred",
               K_grid: [5, 10, 20, 30, 50],
               K_selection: "per split, min mean S1(k=50) on val (train-only mu/sigma); ties -> smaller K; no refit on val",
               reproduction_version: {K: 10},
               P: "rows of G for target genes; ctrl column with zero embedding", ridge_lambda: 0.1,
               fit_rows: "training-pool perturbations whose target is a readout gene, plus ctrl"}
  linear_go: {G: "top-10 PCA scores (fixed K = 10)", P: "adjacency spectral embedding of GEARS go_essential_all.csv (first d of 100 columns)",
              d_grid: [2, 5, 10, 20, 30, 50, 100],
              d_selection: "per split, min mean S1(k=50) on val (train-only mu/sigma); ties -> smaller d; no refit on val",
              reproduction_version: {d: 2, source: "extract_pert_embedding_from_gears.R line 31 (no = 2)"},
              ridge_lambda: 0.1, fit_rows: "training-pool perturbations whose target is in the GO graph, plus ctrl"}
  pca_method: "Gram eigendecomposition; verified against numpy SVD on real data in every P1 split (G rel. diff <= 7.3e-14, predictions <= 1.2e-12, same d); the check is repeated in every split"
  gears: deferred
  gears_unc: deferred
score:
  S1: {k: 50, exclude_target: true, units: "sigma-scaled RMSE on top-k |s| genes"}
  S2: {panel: "H", aggregate: "max", exclude_target: true}
  sensitivity_k: [20, 200]
conformal:
  alpha_primary: 0.10
  alpha_secondary: [0.05, 0.20]
  quantile: "ceil((n+1)(1-alpha))-th order statistic of calibration scores"
distance:                                # pool = training pool (E-train + always_train)
  D_GO: "1 - mean top-20 Jaccard to pool, GEARS gene2go_all.pkl"
  D_GO_filt: "as D_GO after dropping GO terms annotating > 500 of the 9,853 GEARS GO-universe genes (29 terms)"
  D_STRING: "1 - max STRING v12.0 combined experiments+databases score to pool"
analysis:
  primary_test: "binomial GLMM: miscover ~ model * rank(D) + log(n_cells) + (1|target gene) + (1|split); rank(D) = within-split percentile over cal+test"
  primary_parameters: "rank(D) slope per model; linear_go minus mean slope difference (criterion 4, provisional)"
  multiplicity: "Holm over 3 models (mean, linear_pca selected K, linear_go selected d) x 2 distances (D_GO, D_STRING), primary score and alpha"
  grid_edge_rule: "if the selected K or d is the grid maximum, report it; the grid is not extended again"
  robustness: "GEE clustered by TSS group, independence working correlation, cluster-robust SE"
  bins: "quintiles of D over cal+test; descriptive only"
  sensitivity: ["k in {20, 200}", "D_STRING", "D_GO_filt", "linear_pca K = 10 (reproduction)", "linear_go d = 2 (reproduction)"]
p1_selfchecks:
  marginal: "count of (split x 5 models x 4 scores x 3 alpha) combos outside beta-binomial 95% band; pass if P(>= n_out | Bin(n_combos, 0.05)) > 0.01"
  planted: {model: linear_pca, noise: "tau * u * sigma_g, u = test D_GO percentile", tau: 6, pass: "slope > 0 and p < 0.05"}
  permuted: {pass: "p >= 0.05 for every model"}
g1_gate: {audit_A_splits: 100, audit_B_splits: 20, models: [mean, linear_pca, linear_go], reproduction_models: [linear_pca_K10_repro, linear_go_d2_repro]}
audit_B:
  nodes: "E units (1,786)"
  graph: "original GEARS gene2go Jaccard; kNN = 20 (numpy stable argsort), J > 0, symmetrised by max, weight = J"
  clustering: {method: "Leiden RBConfigurationVertexPartition (leidenalg 0.12.0, igraph 1.0.0)", resolution: 8, seed: 20261001,
               frozen_file: "data/processed/v0.5/auditB_clusters_original_r8.csv", n_clusters: 69, max_cluster_frac_E: 0.039}
  stopping_rule: "C (closest-to-15% acceptance; whole clusters plus touched TSS groups)"
  test_fraction_tolerance: [0.135, 0.165]
  feasibility_threshold: 0.95
  resampling: none
  shift_caveat: "design-scan median test D_GO shift vs audit A ~ +0.04; similar coverage = no effect detected at this shift, not absence of structured-OOD effect"
g1_analysis:
  glmm: "lme4 glmer, binomial, bobyqa; miscover ~ model * rank(D) + log_ncells_c + (1|target) + (1|split); fits: {A,B} x {D_GO, D_STRING, D_GO_filt}, plus {A,B} x D_GO without log_ncells"
  glmm_scale_note: "conditional (gene-specific) scale; criteria use sign and p only; magnitudes from GEE"
  holm_family: "audit A, 3 primary models x {D_GO, D_STRING}, two-sided Wald"
  gee: "logistic, independence, clusters = TSS group, robust SE; coverage at rank 0.1 and 0.9; marginal coverage A and B"
  criteria: {C1: "GEE coverage at rank 0.9 < 0.85 and upper CI < 0.90 (A, D_GO)", C2: "GLMM slope > 0, Holm p < 0.05 (A, D_GO)",
             C3: "same GLMM slope sign D_GO vs D_STRING (A)", C4: "GLMM linear_go - mean slope > 0, p < 0.05 (provisional)",
             C5: "GEE audit-B marginal coverage < 0.85",
             support: "v0.5.1: some linear model with C1&C2&C3, and C4 (C5 secondary; reported)",
             refute: "all primary models: GEE coverage at rank 0.1 and 0.9 in [0.85, 0.95] with CI half-width < 0.05; audit-B marginal in [0.85, 0.95]; no Holm-significant slope"}
  reinterpretation: {R1: "mean has Holm-significant positive slope and C4 fails", R2: "Holm-significant negative slope",
                     R3: "positive p < 0.05 without log_ncells but p >= 0.05 with it", R4: "deferred", R5: "deferred"}
  selfcheck: "audit A marginal beta-binomial band, pass if P(>= n_out | Bin(n, 0.05)) > 0.01; otherwise results not interpreted"
code_frozen: [pipeline_v05.py, models_v04.py, audit_lib.py, splits_v04.py, auditB_clusters_v05.py, g1_analysis_v05.py, g1_glmm_v05.R, paths.py]
deferred:
  - "Replogle GWPS (K562 genome-wide)"
  - "Norman 2019 CRISPRa"
  - "GEARS and GEARS-unc (GPU source not decided)"
  - "E_strong subset"
  - "z-normalised expression sensitivity"
  - "L9 GO-annotation temporal leakage check"
```

## 13. 预期产出

- 每次划分 × 模型 × 分箱的结果表（parquet）：覆盖率、q̂、r_b、n。
- 主图：覆盖率随距离分箱的变化（每个模型一条线，带 GLMM 95% CI，名义值用虚线标出）；所需半径 r_b 随距离分箱的变化。
- 补充材料：审计 B 的边际覆盖率、全部敏感性分析、三项自检结果、完整的版本与哈希清单。
