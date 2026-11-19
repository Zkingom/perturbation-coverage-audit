# 新颖性核查：单细胞扰动预测在结构化 OOD 下的不确定性量化 / conformal prediction

检索日期：2026-10-01。检索来源：Europe PMC（含 bioRxiv 预印本）、arXiv API、CrossRef、通用网页检索。Semantic Scholar 在本次检索中被限流（HTTP 429），OpenReview 和 Google Scholar 未覆盖，因此 NeurIPS/ICLR 的 workshop 论文可能有遗漏，投稿前需要再查一次。下表所有 DOI 和 arXiv 编号均已通过 CrossRef 或 arXiv 解析。

## 结论先行

“扰动预测 + 不确定性”这件事已经有人做了，“conformal prediction + Perturb-seq 数据”也已经有人碰过，但“在结构化 OOD 划分（未见扰动、未见细胞类型）下，为扰动响应预测模型提供有限样本的覆盖率/风险控制保证，并系统审计现有模型在这些划分下的校准情况”，在本次检索范围内没有找到直接做过的工作。最接近的三篇是 PRESCRIBE（[Cheng et al. 2025](https://arxiv.org/abs/2510.07964)，证据回归，没有覆盖率保证）、dbDiffusion（[Shang et al. 2025](https://doi.org/10.1073/pnas.2525268122)，用 prediction-powered inference 给出平均效应的 CI，但只用“CI 重叠比例”评估，没有报告覆盖率）、以及 [Asiaee et al. 2026](https://arxiv.org/abs/2603.02204)（在 Replogle K562 上做 selective conformal，但问题是用因果结构挑选可交换的校准干预，不是预测未见扰动的响应）。这意味着：只做“把 split conformal 套到 GEARS 上”不够新；只有把交换性单元、结构化偏移和生物学决策这三件事一起讲清楚，才有方法贡献。

## 1. 扰动预测模型及其自带的不确定性

现有模型中，自带不确定性的大致分三类：异方差方差头（GEARS）、潜空间距离或贝叶斯后验（CPA、GPerturb、PRESCRIBE），以及基于统计推断的置信区间（dbDiffusion）。它们有一个共同点：要么只评估“不确定性和误差是否相关”（排序意义上的校准），要么评估 ECE 一类的平均校准，都不提供分布无关的覆盖率保证，也没有按 OOD 程度分层检验覆盖率。

| 工作 | 年份 | DOI / arXiv | 做了什么 | 与本课题的区别 |
|---|---|---|---|---|
| scGen: predicts single-cell perturbation responses ([Lotfollahi et al. 2019](https://doi.org/10.1038/s41592-019-0494-8)) | 2019 | 10.1038/s41592-019-0494-8 | VAE 潜空间向量运算，预测留出细胞类型的刺激响应 | 最早的“未见细胞类型”设定之一；只有点预测 |
| CPA: Predicting cellular responses to complex perturbations in high-throughput screens ([Lotfollahi et al. 2023](https://doi.org/10.15252/msb.202211517)) | 2023 | 10.15252/msb.202211517 | 组合扰动/剂量/协变量的解耦自编码器，按潜空间与训练条件的距离给出 uncertainty 值 | 启发式距离分数，无覆盖率含义 |
| chemCPA: Predicting cellular responses to novel drug perturbations ([Hetzel et al. 2022](https://arxiv.org/abs/2204.13545)) | 2022 | arXiv:2204.13545 | 用分子表示把 CPA 推广到未见药物，借助 bulk 数据迁移 | 未见扰动（药物）设定；无 UQ |
| GEARS: Predicting transcriptional outcomes of novel multigene perturbations ([Roohani et al. 2023](https://doi.org/10.1038/s41587-023-01905-6)) | 2023 | 10.1038/s41587-023-01905-6 | GO 知识图谱 + GNN 预测未见单/双基因扰动；高斯似然学习 log-variance 作为不确定性分数，正文称其与模型表现负相关 | 最常用的 baseline；其不确定性是异方差回归分数，未经校准，不保证覆盖率 |
| CellOT: Learning single-cell perturbation responses using neural optimal transport ([Bunne et al. 2023](https://doi.org/10.1038/s41592-023-01969-x)) | 2023 | 10.1038/s41592-023-01969-x | 神经最优传输学习未配对的扰动前后分布映射，评估了分布外样本泛化 | 分布级预测，无 UQ |
| biolord: Disentanglement of single-cell data ([Piran et al. 2024](https://doi.org/10.1038/s41587-023-02079-x)) | 2024 | 10.1038/s41587-023-02079-x | 解耦已知属性，用于反事实预测（含未见扰动） | 无 UQ |
| SAMS-VAE ([Bereket & Karaletsos 2023](https://arxiv.org/abs/2311.02794)) | 2023 | arXiv:2311.02794 | 稀疏加性机制偏移 VAE，变分后验 | 贝叶斯后验，无频率学覆盖保证 |
| Sparse mechanism shift modeling ([Lopez et al. 2022](https://arxiv.org/abs/2211.03553)) | 2022 | arXiv:2211.03553 | 稀疏机制偏移的因果表示学习，用于未见扰动组合 | 无 UQ |
| scGPT ([Cui et al. 2024](https://doi.org/10.1038/s41592-024-02201-0)) | 2024 | 10.1038/s41592-024-02201-0 | 单细胞基础模型，微调后做扰动预测 | 点预测；预训练语料可能与测试细胞系重叠（需自查） |
| scFoundation ([Hao et al. 2024](https://doi.org/10.1038/s41592-024-02305-7)) | 2024 | 10.1038/s41592-024-02305-7 | 大规模基础模型，嵌入接入 GEARS 做扰动预测 | 同上 |
| Geneformer ([Theodoris et al. 2023](https://doi.org/10.1038/s41586-023-06139-9)) | 2023 | 10.1038/s41586-023-06139-9 | rank-based 基础模型，in silico 删除 | 不是定量响应预测器；无 UQ |
| GenePert ([Chen & Zou 2024](https://doi.org/10.1101/2024.10.27.620513)) | 2024 | 10.1101/2024.10.27.620513 | 用 GenePT（文本 LLM）基因嵌入做正则化回归 | 预印本；文本嵌入可能编码了关于测试基因的文献知识（潜在泄漏源） |
| scLAMBDA ([Wang et al. 2024](https://doi.org/10.1101/2024.12.04.626878)) | 2024 | 10.1101/2024.12.04.626878 | LLM 嵌入 + 解耦生成模型，预测未见多基因扰动 | 预印本；无覆盖率 |
| PerturbNet ([Yu et al. 2025](https://doi.org/10.1038/s44320-025-00131-3)) | 2025 | 10.1038/s44320-025-00131-3 | 条件生成模型预测未见药物/基因扰动的细胞状态分布 | 生成分布 ≠ 校准的预测区间 |
| GPerturb ([Xing & Yau 2025](https://doi.org/10.1038/s41467-025-61165-7)) | 2025 | 10.1038/s41467-025-61165-7 | 高斯过程稀疏回归，给出基因层面效应存在与强度的不确定性 | 模型内后验，依赖模型正确设定；非 OOD 校准 |
| dbDiffusion ([Shang et al. 2025](https://doi.org/10.1073/pnas.2525268122)) | 2025 | 10.1073/pnas.2525268122 | 潜扩散生成 + 借鉴 prediction-powered inference 的去偏步骤，对未见扰动的每基因平均效应给 CI，并用于 DE 推断 | 最接近的“有统计推断”的竞品。CI 依赖“相似扰动簇的偏差可迁移”这一假设；评估指标是 CI 与观测 CI 的重叠比例，未报告覆盖率，未做细胞类型 OOD |
| PRESCRIBE ([Cheng et al. 2025](https://arxiv.org/abs/2510.07964)) | 2025 | arXiv:2510.07964（NeurIPS 2025） | 多元深度证据回归，联合估计偶然/认知不确定性（pseudo E-distance），报告 ECE，并用置信度过滤预测以提高准确率；baseline 含 GEARS-Dropout/Ensemble | 最直接的竞品。属于启发式校准（ECE、排序），无有限样本覆盖保证；只针对未见基因，不涉及未见细胞类型 |
| State ([Adduri et al. 2025](https://doi.org/10.1101/2025.06.26.661135)) | 2025 | 10.1101/2025.06.26.661135 | 跨上下文（细胞类型）的扰动响应模型，附 Cell-Eval 评估框架 | 预印本；跨上下文泛化的主要大模型竞品，无 UQ |
| CellFlow ([Klein et al. 2025](https://doi.org/10.1101/2025.04.11.648220)) | 2025 | 10.1101/2025.04.11.648220 | flow matching 生成单细胞表型 | 预印本；无覆盖率 |
| CFM-GP ([Abir et al. 2026](https://doi.org/10.1093/nargab/lqag087)) | 2026 | 10.1093/nargab/lqag087 | 条件 flow matching，跨细胞类型学习基因扰动 | 跨细胞类型设定；无 UQ |

## 2. 基准测试与对现有模型的批评

2024–2025 年的一批基准研究形成了一个共识：在常用指标上，深度模型和基础模型通常打不过均值或线性 baseline（[Ahlmann-Eltze et al. 2025](https://doi.org/10.1038/s41592-025-02772-6)；[Csendes et al. 2025](https://doi.org/10.1186/s12864-025-11600-2)；[Kernfeld et al. 2023](https://doi.org/10.1101/2023.07.28.551039)）。随后的工作把问题部分归因于指标本身：Systema 指出常用指标被扰动无关的“系统性变异”主导（[Viñas Torné et al. 2025](https://doi.org/10.1038/s41587-025-02777-8)），Miller 等人则认为 MSE、Pearson(Δctrl) 这类指标校准很差，换用校准良好的指标后深度模型确实优于无信息 baseline（[Miller et al. 2025](https://doi.org/10.1101/2025.10.20.683304)，预印本）。对本课题最关键的经验事实是：误差随 OOD 程度增大，两个大规模基准都报告了扰动效应大小和源/目标细胞类型的转录组距离会显著影响精度（[Li et al. 2024](https://doi.org/10.1101/2024.12.23.630036)；[Wei et al. 2025](https://doi.org/10.1038/s41592-025-02980-0)）。这正是“随机划分上校准的 conformal 区间在远 OOD 单元上会欠覆盖”这一假设的依据，但据我检索，还没有人直接测过。

| 工作 | 年份 | DOI / arXiv | 做了什么 | 与本课题的关系 |
|---|---|---|---|---|
| Norman et al., genetic interaction manifolds ([Norman et al. 2019](https://doi.org/10.1126/science.aax4438)) | 2019 | 10.1126/science.aax4438 | K562 CRISPRa 单/双基因 Perturb-seq | 核心数据；Miller 等人指出加性 baseline 在该数据上已饱和，组合覆盖有限 |
| Replogle et al., genome-scale Perturb-seq ([Replogle et al. 2022](https://doi.org/10.1016/j.cell.2022.05.013)) | 2022 | 10.1016/j.cell.2022.05.013 | K562 全基因组 + RPE1 必需基因 CRISPRi | 核心数据；两个细胞系共享大量扰动，是“未见细胞类型”划分中的泄漏风险点 |
| Nadig et al., perturbation atlases DE ([Nadig et al. 2025](https://doi.org/10.1038/s41588-025-02169-3)) | 2025 | 10.1038/s41588-025-02169-3 | 新增 HepG2、Jurkat 等 Perturb-seq，提出 TRADE 估计 | 可作为更多细胞类型/较晚发布的伪前瞻测试集 |
| sci-Plex ([Srivatsan et al. 2020](https://doi.org/10.1126/science.aax6234)) | 2020 | 10.1126/science.aax6234 | 3 个细胞系 × 188 化合物 × 剂量 | 化学扰动的标准数据；细胞类型只有 3 个 |
| Tahoe-100M ([Zhang et al. 2025](https://doi.org/10.1101/2025.02.20.639398)) | 2025 | 10.1101/2025.02.20.639398 | 约 50 个癌细胞系 × 上千药物条件的单细胞图谱 | 预印本；细胞系数量足以做 leave-one-cell-line-out 的群组层面覆盖率评估 |
| scPerturb ([Peidli et al. 2024](https://doi.org/10.1038/s41592-023-02144-y)) | 2024 | 10.1038/s41592-023-02144-y | 44 个数据集的统一处理 + E-distance | 数据版本来源，需记录版本号 |
| Ahlmann-Eltze et al. ([2025](https://doi.org/10.1038/s41592-025-02772-6)) | 2025 | 10.1038/s41592-025-02772-6 | 5 个基础模型 + 2 个深度模型均未超过简单 baseline | 你的 conformal 实验必须包含 “均值/线性 baseline + conformal” |
| Systema ([Viñas Torné et al. 2025](https://doi.org/10.1038/s41587-025-02777-8)) | 2025 | 10.1038/s41587-025-02777-8 | 区分系统性变异与扰动特异效应的评估框架 | 非一致性分数若基于全基因 MSE，会被系统性变异主导 |
| Miller et al. ([2025](https://doi.org/10.1101/2025.10.20.683304)) | 2025 | 10.1101/2025.10.20.683304 | 用阳性/阴性对照衡量指标校准（DRF） | 预印本；可借来定义有生物学意义的非一致性分数 |
| Diversity by design ([Mejia et al. 2025](https://arxiv.org/abs/2506.22641)) | 2025 | arXiv:2506.22641 | 指出模式坍塌，在校准良好指标上重新评估 | 预印本 |
| Csendes et al. ([2025](https://doi.org/10.1186/s12864-025-11600-2)) | 2025 | 10.1186/s12864-025-11600-2 | 训练均值 baseline 优于 scGPT、scFoundation | 同上 |
| Kernfeld et al., expression forecasting ([2023](https://doi.org/10.1101/2023.07.28.551039)) | 2023 | 10.1101/2023.07.28.551039 | PEREGGRN 平台，结论高度依赖指标 | 预印本 |
| PertEval-scFM ([Wenteler et al. 2024](https://doi.org/10.1101/2024.10.02.616248)) | 2024 | 10.1101/2024.10.02.616248 | 基础模型嵌入用于扰动效应预测的基准，关注分布偏移 | 预印本 |
| PerturBench ([Wu et al. 2024](https://arxiv.org/abs/2408.10609)) | 2024 | arXiv:2408.10609 | 多任务扰动预测基准（含协变量迁移） | 可复用其划分 |
| Li et al., systematic comparison ([2024](https://doi.org/10.1101/2024.12.23.630036)) | 2024 | 10.1101/2024.12.23.630036 | 13 方法 × 25 数据集；效应大小与细胞类型距离影响精度 | 预印本；“误差随 OOD 程度增大”的直接证据 |
| Wei et al., generalizable prediction benchmark ([2025](https://doi.org/10.1038/s41592-025-02980-0)) | 2025 | 10.1038/s41592-025-02980-0 | 在未见细胞上下文和扰动场景下评估泛化 | 同上 |
| Virtual Cell Challenge ([Roohani et al. 2025](https://doi.org/10.1016/j.cell.2025.06.008)) | 2025 | 10.1016/j.cell.2025.06.008 | 虚拟细胞挑战赛设定 | 我没有核实其评分是否包含 UQ，不作断言 |
| How to build the virtual cell with AI ([Bunne et al. 2024](https://doi.org/10.1016/j.cell.2024.11.015)) | 2024 | 10.1016/j.cell.2024.11.015 | 虚拟细胞路线图 | 背景引用 |
| Toward trustworthy virtual cells ([Li et al. 2026](https://doi.org/10.3389/fcell.2026.1900624)) | 2026 | 10.3389/fcell.2026.1900624 | 综述，明确把“未见细胞类型/扰动的泛化”和“不确定性校准”列为应有的评估标准 | 可作为“领域认为这是空白”的引用，但它是观点性综述 |
| Explicit mechanistic operators in virtual cells ([Stańczak & Niazi 2026](https://doi.org/10.3390/cells15181711)) | 2026 | 10.3390/cells15181711 | 观点文章，总结排名依赖划分、统计单元、指标 | 观点性；支持“划分和统计单元要严格”的论点 |
| Can We Trust In-Distribution Success? ([Shoeibi & Yousefi 2026](https://arxiv.org/abs/2608.00152)) | 2026 | arXiv:2608.00152 | 锁定评估下 CRISPRi 预测的迁移失败，按靶基因做 cluster bootstrap | 预印本；其按靶基因聚类的 bootstrap 方案值得借鉴 |

## 3. 分布偏移下的 conformal 方法

### 核心理论工具

conformal 文献已经为本课题的每一种偏移提供了现成工具，难点在于把生物学结构映射到正确的交换性假设上。标准 split conformal（[Lei et al. 2018](https://doi.org/10.1080/01621459.2017.1307116)）和 CQR（[Romano et al. 2019](https://arxiv.org/abs/1905.03222)）要求校准与测试单元可交换。协变量偏移可用加权 conformal 处理（[Tibshirani et al. 2019](https://arxiv.org/abs/1904.06019)），非交换数据的覆盖损失可以被界定（[Barber et al. 2023](https://doi.org/10.1214/23-AOS2276)），在线漂移可用自适应 conformal（[Gibbs & Candès 2021](https://arxiv.org/abs/2106.00170)）。对“未见细胞类型”最贴切的是分层/两层交换性：新群组（新细胞类型）的预测集合（[Dunn et al. 2022](https://doi.org/10.1080/01621459.2022.2060112)；[Lee et al. 2026](https://doi.org/10.1145/3786352)），以及已有少量新群组样本时的 GHCP（[Mallick et al. 2026](https://arxiv.org/abs/2608.15500)，预印本），后者正好对应“新细胞类型先做一个小规模 pilot 筛选”的场景。按 OOD 程度分层的覆盖可以用条件保证框架（[Gibbs et al. 2025](https://doi.org/10.1093/jrsssb/qkaf008)）。超越覆盖率的损失（如 DEG 集合的假阴性率）可用 conformal risk control（[Angelopoulos et al. 2022](https://arxiv.org/abs/2208.02814)，ICLR 2024）和 Learn-then-Test（[Angelopoulos et al. 2021](https://arxiv.org/abs/2110.01052)）。

| 工作 | 年份 | DOI / arXiv | 做了什么 | 对本课题的用处 / 区别 |
|---|---|---|---|---|
| Distribution-free predictive inference for regression ([Lei et al. 2018](https://doi.org/10.1080/01621459.2017.1307116)) | 2018 | 10.1080/01621459.2017.1307116 | split conformal 回归 | 基础 |
| Conformalized quantile regression ([Romano et al. 2019](https://arxiv.org/abs/1905.03222)) | 2019 | arXiv:1905.03222 | 自适应宽度区间 | 适合异方差的扰动效应 |
| Conformal prediction under covariate shift ([Tibshirani et al. 2019](https://arxiv.org/abs/1904.06019)) | 2019 | arXiv:1904.06019 | 似然比加权 conformal | 可在扰动嵌入空间估计权重；高维权重估计不稳定是已知问题 |
| Conformal prediction beyond exchangeability ([Barber et al. 2023](https://doi.org/10.1214/23-AOS2276)) | 2023 | 10.1214/23-AOS2276 | 非交换下的覆盖差界 | 理论上界定 OOD 下的覆盖损失 |
| Adaptive conformal inference under distribution shift ([Gibbs & Candès 2021](https://arxiv.org/abs/2106.00170)) | 2021 | arXiv:2106.00170 | 在线调整分位数 | 适合序贯筛选/主动学习场景 |
| Conformal prediction with conditional guarantees ([Gibbs et al. 2025](https://doi.org/10.1093/jrsssb/qkaf008)) | 2025 | 10.1093/jrsssb/qkaf008 | 对函数类（如分组、协变量偏移族）的条件覆盖 | 可实现“按与训练扰动的距离分层”的覆盖 |
| Robust validation ([Cauchois et al. 2024](https://doi.org/10.1080/01621459.2023.2298037)) | 2024 | 10.1080/01621459.2023.2298037 | f-散度球内最坏情况覆盖 | 未见细胞类型的保守备选 |
| PAC prediction sets under covariate shift ([Park et al. 2021](https://arxiv.org/abs/2106.09848)) | 2021 | arXiv:2106.09848 | PAC 型集合 + 重要性权重 | 训练-条件保证，适合小校准集 |
| Label shift CP ([Podkopaev & Ramdas 2021](https://arxiv.org/abs/2103.03323)) | 2021 | arXiv:2103.03323 | 标签偏移下的分类 CP | 细胞类型比例变化时相关性有限 |
| Two-layer hierarchical models ([Dunn et al. 2022](https://doi.org/10.1080/01621459.2022.2060112)) | 2022 | 10.1080/01621459.2022.2060112 | 新群组的分布无关预测集合 | “未见细胞类型”的标准框架；群组数少时保证很弱 |
| Hierarchical data ([Lee et al. 2026](https://doi.org/10.1145/3786352)) | 2026 | 10.1145/3786352 | 分层数据的分布无关推断 | 同上 |
| Generalized hierarchical CP ([Mallick et al. 2026](https://arxiv.org/abs/2608.15500)) | 2026 | arXiv:2608.15500 | 新群组已有少量观测时的 HCP | 预印本；直接对应 pilot 筛选 |
| Few-shot CP with auxiliary tasks ([Fisch et al. 2021](https://arxiv.org/abs/2102.08898)) | 2021 | arXiv:2102.08898 | 元学习 + 任务层交换性 | 每个扰动/细胞类型可视为一个任务 |
| Conformal risk control ([Angelopoulos et al. 2022](https://arxiv.org/abs/2208.02814)) | 2022 | arXiv:2208.02814 | 控制单调损失的期望 | 控制预测 DEG 集合的 FNR 等生物学损失 |
| Learn then Test ([Angelopoulos et al. 2021](https://arxiv.org/abs/2110.01052)) | 2021 | arXiv:2110.01052 | 多重检验框架下的风险控制 | 非单调损失或多指标同时控制 |
| Copula-based multi-target CP ([Messoudi et al. 2021](https://doi.org/10.1016/j.patcog.2021.108101)) | 2021 | 10.1016/j.patcog.2021.108101 | 多输出联合覆盖 | 数千基因的联合覆盖 |
| Testing for outliers with conformal p-values ([Bates et al. 2023](https://doi.org/10.1214/22-AOS2244)) | 2023 | 10.1214/22-AOS2244 | conformal p 值与 FDR | 识别“模型不可信”的 OOD 扰动 |
| Selection by prediction with conformal p-values ([Jin & Candès 2023](https://arxiv.org/abs/2210.01408)) | 2023 | arXiv:2210.01408（JMLR） | 带 FDR 控制的选择 | 筛选“预测强效应”的候选扰动 |
| Selection from hierarchical data with conformal e-values ([Lee & Ren 2025](https://arxiv.org/abs/2501.02514)) | 2025 | arXiv:2501.02514 | 分层数据下的 FDR 控制选择 | 跨细胞类型的候选选择 |
| Conformal inference of counterfactuals and ITEs ([Lei & Candès 2021](https://doi.org/10.1111/rssb.12445)) | 2021 | 10.1111/rssb.12445 | 反事实/个体处理效应的 conformal 区间 | 扰动效应本质上是反事实量，概念上最近 |
| Robust conformal sensitivity analysis of ITEs ([Jin et al. 2023](https://doi.org/10.1073/pnas.2214889120)) | 2023 | 10.1073/pnas.2214889120 | 混杂偏移下的稳健 conformal | 偏移有界时的稳健版本 |

### 已用于生物学的 conformal

conformal 已经进入生物学，但集中在分类、bulk 药物敏感性和分子设计，还没有进入“单细胞扰动响应预测 + 结构化 OOD”。

| 工作 | 年份 | DOI / arXiv | 做了什么 | 区别 |
|---|---|---|---|---|
| **Partial causal structure learning for valid selective conformal inference under interventions ([Asiaee et al. 2026](https://arxiv.org/abs/2603.02204))** | 2026 | arXiv:2603.02204 | 在干预实验中，只用“不影响目标变量”的干预作校准集以收紧区间；估计干预-靶点关系，给出误纳入时的覆盖退化与校正；在 Replogle K562 上实验 | **目前找到的唯一一篇 conformal + Perturb-seq 方法文。** 它预测的是某变量在已知干预下的取值，校准依据是因果不变性；不涉及用预测模型外推到未见扰动或未见细胞类型。论文中必须引用并区分 |
| Conformal prediction under feedback covariate shift ([Fannjiang et al. 2022](https://doi.org/10.1073/pnas.2204569119)) | 2022 | 10.1073/pnas.2204569119 | 蛋白设计中测试点由模型选出导致的偏移 | 若做“模型推荐下一批扰动”，就是这种反馈偏移，必须考虑 |
| CoDrug ([Laghuvarapu et al. 2023](https://arxiv.org/abs/2310.12033)) | 2023 | arXiv:2310.12033 | 分子性质预测中用密度估计做加权 conformal | 未见化合物的类比；单输出 |
| Reliable anti-cancer drug sensitivity prediction ([Lenhof et al. 2024](https://doi.org/10.1038/s41598-024-62956-6)) | 2024 | 10.1038/s41598-024-62956-6 | bulk 细胞系药敏的 conformal 回归/分类 | bulk、单输出，非单细胞转录组响应 |
| Conformal inference for reliable scRNA-seq annotation ([López-De-Castro et al. 2025](https://doi.org/10.1093/bioinformatics/btaf521)) | 2025 | 10.1093/bioinformatics/btaf521 | 细胞注释的 conformal 集合 + OOD 细胞类型检测 | 分类任务，非扰动预测 |
| Conformal cell type annotation with graph constraints ([Corbetta et al. 2024](https://arxiv.org/abs/2410.23786)) | 2024 | arXiv:2410.23786 | 带层级约束的注释 conformal | 同上 |
| CILANTRO-SL ([Hua et al. 2026](https://doi.org/10.64898/2026.02.25.708096)) | 2026 | 10.64898/2026.02.25.708096 | 基础模型 in silico 敲除嵌入 + conformal 给合成致死对校准置信度 | 预印本；二分类，不是表达响应 |
| ReliST ([Zhang et al. 2026](https://doi.org/10.1016/j.isci.2026.117206)) | 2026 | 10.1016/j.isci.2026.117206 | 空间转录组去卷积的模型无关风险层 | 不同任务；但“模型无关可靠性层”的定位与本课题相似 |

## 真正的空白

### 什么已经被做过

“为未见扰动给出不确定性分数并验证其与误差相关”已经被 GEARS、PRESCRIBE 和 dbDiffusion 覆盖，单独做这件事没有新颖性。“conformal 用在 Perturb-seq 数据上”也已经出现（Asiaee et al. 2026），所以不能声称“首次把 conformal 引入扰动数据”。

### 本次检索中没有找到的

第一，**以扰动为交换单元、按 OOD 程度分层的覆盖率审计**。现有工作都在细胞层或平均意义上评估校准，没人报告过：在随机扰动划分上校准的区间，对“与训练扰动功能距离远”的扰动和“未见细胞类型”是否欠覆盖、欠多少。基准文献（Li et al. 2024；Wei et al. 2025）已证明误差随距离增大，所以这个结果大概率是“欠覆盖且有系统规律”，这本身就是可发表的经验发现。

第二，**未见细胞类型的群组层面保证**。分层 conformal 的理论（Dunn et al. 2022；Lee et al. 2026；GHCP 2026）已经存在，但没有被用在扰动预测上；“用新细胞类型的小规模 pilot 筛选做少样本校准”这一场景尤其契合 GHCP，并且对应真实的实验决策。

第三，**对生物学可操作输出做风险控制，而不是对每个基因做区间**。生物学家需要的是“这个扰动会让哪些基因显著上调/下调”或“会不会使某个程序（如红系分化）改变”，可以用 conformal risk control 控制预测 DEG 集合的假阴性率或方向错误率，或用 conformal selection 在候选扰动中控制 FDR。这在扰动预测中没有找到先例。

第四，**把不确定性作为解决“深度模型是否胜过线性 baseline”争论的新维度**。在同样的覆盖率下比较区间/集合的宽度（效率），可以回答“深度模型是否在分布外提供了更多可信信息”，这比再报一个 Pearson 相关更有说服力，也直接回应 Ahlmann-Eltze、Systema 和 Miller 之间的分歧。

## 审稿人视角的主要风险

### 新颖性风险

如果方法部分只是把现成的加权 conformal 或 CRC 套到 GEARS/scGPT 上，生物学顶刊会把它视为应用，CS 审稿人会视为增量。方法贡献至少要有一个具体点，例如：扰动相似性加权且有保证的局部化校准；扰动×细胞类型的双层交换性下的覆盖；或基因间相关结构下的 DEG 集合风险控制。

### 统计效力风险

“未见细胞类型”的群组数量太少。Replogle + Nadig 只有 K562、RPE1、HepG2、Jurkat 四个细胞系，sci-Plex 只有三个；在 K=4 时，群组层面的覆盖率保证几乎无法经验验证（每个 leave-one-out 只给一个覆盖率观测）。需要 Tahoe-100M 这类约 50 个细胞系的数据（预印本，需核实许可与版本）。另外，Norman 只有约百级的单/双扰动，以扰动为单元的校准集很小，覆盖率本身的抽样波动可达数个百分点，必须报告多次随机划分的分布和置信区间，而不是单次结果。

### 泄漏风险

需要逐项排查：(a) Norman 双扰动中，“两个单扰动都已在训练集”与“完全未见”必须分开（GEARS 的 seen 0/1/2 划分）；(b) 同一基因扰动在 K562 与 RPE1 中都出现，“未见细胞类型”划分要说清该扰动是否在其他细胞系中被见过；(c) GO 图、GenePT 文本嵌入可能编码了关于测试基因扰动效应的文献知识；(d) 基础模型预训练语料是否包含测试细胞系（本次未核实，需自查）；(e) 同一靶基因的多个 guide 必须划在同一侧。

### 指标风险

如果非一致性分数用全基因 MSE，平凡的“训练均值 + conformal”也能达到名义覆盖，而且区间宽度可能和深度模型差不多（Systema 和 Miller 的论点在这里同样适用）。主要比较对象必须包括均值/线性 baseline + 同样的 conformal 流程，评价指标应是同覆盖下的效率和在扰动特异信号上的信息量。

## 最可能的生物学应用场景

### 带 FDR 控制的 in silico 筛选优先级

这是最有说服力的场景：从全基因组候选中，选出“预测会显著改变目标程序”的扰动送去实验，并保证被选集合中假阳性比例不超过设定值。方法上对应 conformal selection（Jin & Candès 2023；分层版本 Lee & Ren 2025），可以用 Replogle 全基因组 K562 回顾性模拟：只用部分扰动训练，看在其余扰动中选出的集合实际 FDR 是否受控。若采用“模型推荐 → 实验 → 再训练”的循环，需要处理反馈偏移（Fannjiang et al. 2022）。

### 迁移到新细胞类型时的“信任或补实验”决策

在新细胞系或原代细胞中，先做小规模 pilot 扰动，再用 GHCP 一类的少样本校准判断模型预测在新环境中是否足够可信，从而决定剩余扰动是否需要实测。这对应真实的实验预算决策，适合在 Tahoe-100M（药物、约 50 细胞系）或 Replogle/Nadig（四细胞系，仅作案例）上展示。

### 组合扰动的遗传相互作用

可以用来判断哪些双基因组合的预测“足够不确定、值得测”。但 Miller 等人指出 Norman 数据上加性 baseline 已饱和，组合覆盖有限，所以这个场景作为主线风险较高，建议作为次要案例。

## 下一步建议

1. 投稿前再用 Google Scholar/OpenReview 检索 “conformal” + “perturbation”/“Perturb-seq”/“virtual cell”，并追踪 PRESCRIBE、dbDiffusion、Asiaee et al. 2026 的施引文献。
2. 先做一个低成本的“审计实验”：Norman + Replogle K562，GEARS、线性、均值三个模型，split conformal（以扰动为单元），按与训练扰动的距离分箱报告覆盖率，重复 ≥20 次随机划分。如果远 OOD 分箱确实欠覆盖，就确立了问题；如果不欠覆盖，这个方向的动机需要重新考虑。
3. 早期确定是否能使用 Tahoe-100M 或其他多细胞类型数据集，否则“未见细胞类型”这一半的结论在统计上站不住。
