# 编号方案（提案，2026-10-28，尚未修改正文）

规则：按 Abstract → Introduction → Methods → Results → Discussion 的顺序，图、表各自按首次引用依次编号。

## 主图
| 新编号 | 内容 | 现编号 |
|---|---|---|
| Fig. 1 | (a) 三个数据集的分层覆盖；(b) 距离/难度代理；(c) 先验邻居代理的修复 | Fig. 1（fig1_main_v4，不变） |
| Fig. 2（可选，需新画） | 耦合标尺：空白 / linear-PCA / linear-GO / GEARS / 正对照 / 阴性对照的 R，附 GEARS 的分层覆盖与解析参照 0.5 | 无 |

## 补充图
| 新编号 | 内容 | 现编号 |
|---|---|---|
| Fig. S1 | 分数与实测效应强度的关系，以及各模型的 CV | Fig. S2（figS_scale） |
| Fig. S2 | 噪声底校正（方案外说明性分析） | Fig. S1（figS_noise_floor） |
| Fig. S3 | 确认后探索分析：度数分层与分半参照 | Fig. S3（figS_explore_modularity；正文中残留的写法"S_explore"统一改掉） |

## 补充表
| 新编号 | 内容 | 现编号 |
|---|---|---|
| Table S1 | 冻结时间线与全部偏离 | S3 |
| Table S2 | 参照方法逐划分的宽度匹配 | S1 |
| Table S3 | 构造性参照与实测 T5（含存档蒙特卡洛列） | S8 |
| Table S4 | 评估单元的入选偏倚（<30 细胞） | S4 |
| Table S5 | 分数与效应强度的关系（各模型） | S6 |
| Table S6 | 确认性决策表（RPE1、非文库集） | S2（**正文目前未引用，需补引用**） |
| Table S7 | 数字溯源表 | S5（**正文目前未引用，需补引用**） |
| Table S8 | 事先可识别子集的覆盖 | S7 |
| Table S9 | CORUM 子集覆盖（探索性） | S9（不变） |
| Table S10 | GEARS 耦合检验（v1.0 主要量、v1.1 三项对照） | S10（已定） |

**Supplementary Note 1** 保持唯一一份。另建议把 GEARS 的偏离与审计细节单独整理成 **Supplementary Note 2**（可选），Methods 中只保留摘要。

## 需要你决定
1. 是否接受上面的重新编号。按首次引用顺序编号会改动大部分表号；另一种做法是保留现有编号，只给新增内容补号（S10）。
2. Fig. 2 是否要新画。不画的话，GEARS 结果只放在 Table S10。
3. 章节顺序：多数期刊（例如 Genome Biology）把 Methods 放在 Discussion 之后。如果改成 Abstract → Introduction → Results → Discussion → Methods，编号会按新顺序重排。
4. 是否设 Supplementary Note 2。
