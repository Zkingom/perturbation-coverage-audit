# 迁移计划：工作文件迁至 `F:\codex_in_this\claude_science`（待确认，尚未执行）

日期：2026-10-02。本文件只是计划，没有复制、下载或改动任何文件。

## 0. 先说明一个重要情况：原工作区的文件已被清理

本会话恢复时，系统提示工作区在闲置期间被清理了，我已核实：工作区里只剩环境自身的文件。下面这些都**已经不在磁盘上**：
- 全部原始下载（h5ad、GO 先验、STRING、Zenodo 代码包等）；
- 中间数据（pseudobulk、索引、STRING 子网等）；
- `splits\` 下的只读划分文件；
- GEARS 源码。

仍然完好的，只有应用管理的 artifact 存储中的 19 个文件版本（含方案的 3 个历史版本）。

因此，“先复制、后验证”需要分成两条来源执行：

1. **artifact（19 个文件）**：从 artifact 存储**只读**复制到新路径，然后逐个比对 SHA256。只读取，不修改 `~/.claude-science` 内的任何文件。
2. **原始数据（11 个文件）**：无法复制，只能重新下载到 `G:\biology`，再与 `data_manifest_v0.3.json` 里记录的 SHA256，以及发布方公布的 MD5（Figshare、Dataverse、Zenodo）或 SHA256（PyPI）核对；全部一致后，再复制到 `data\raw\`。
3. **中间数据和划分文件（6 项）**：只能用原代码和固定种子重新生成。它们以前没有记录过哈希，所以“与迁移前一致”要用内容比对来证明（见第 5 节的冒烟测试）。

第 3 条“全部一致之前不删除原文件”自然满足。工作区里已经没有原文件可删；artifact 存储里的副本由应用管理，**本计划不删除任何东西**，所以也不需要你再单独确认删除。

## 1. 需要你授权的两个路径

| 路径 | 用途 | 权限 |
|---|---|---|
| `F:\codex_in_this\claude_science` | 项目根目录（下文写作 `<ROOT>`） | 读写 |
| `G:\biology` | 下载专用目录（下文写作 `<DL>`） | 读写 |

执行前我会通过授权弹窗分别申请这两个路径。另外需要你告诉我 F: 和 G: 的可用空间。本次迁移估计 F: 和 G: 各需约 4.5 GB（下载 0.99 GB，加解压后的 GEARS 版 h5ad 2.90 GB，加解压后的 GO 边表约 0.4 GB）；以后的单细胞文件（10.7 GB）需另算。

## 2. 目录结构

```
<ROOT>/                                F:\codex_in_this\claude_science
├── config.toml                        project_root, download_root, master_seed (NOT the app config)
├── code/                              audit_lib.py, p1_run.py, p0_build.py, paths.py, smoke_test.py
├── data/
│   ├── raw/                           verified copies of downloads (read-only after verification)
│   │   ├── replogle2022_figshare/
│   │   ├── gears_priors/              gene2go_all.pkl, essential_all_data_pert_genes.pkl, go_essential_all.tar.gz
│   │   ├── gears_k562_essential/      replogle_k562_essential.zip, perturb_processed.h5ad
│   │   ├── string_v12.0/
│   │   ├── ahlmann_eltze_2025/
│   │   └── software/                  cell-gears-0.1.2.tar.gz (+ extracted source)
│   └── processed/
│       ├── priors/                    go_ase_embeddings.npz, string_expdb_universe_*.parquet
│       ├── v0.2_gears_k562ess/        p0_units_E.csv, pseudobulk_mean.npy, conditions.csv, readout_genes.csv
│       └── v0.3_figshare_k562ess/     p0v3_units_all_targets.csv
├── splits/
│   └── v0.2/                          A_000.json ... A_002.json (read-only, L10)
├── results/
│   ├── p1_v0.2/
│   └── p0_v0.3/
├── reports/
│   ├── protocol/                      audit_protocol_v0.1.md ... v0.4.md
│   └── *.md                           novelty review, P0/P1 reports, project_status.md
└── manifests/                         data_manifest_v0.3.json, p1_run_manifest_v0.2.json,
                                       migration_checksums.csv, env_pert-audit.yml
```

这里的 `config.toml` 是**项目自己的配置文件**，放在 `<ROOT>` 下，与应用的 config.toml 无关；应用的 config.toml 不会被改动。

## 3. 完整清单

### 3a. artifact 副本（从 artifact 存储只读复制）

| 文件 | 当前位置（artifact version） | 大小 | 目标路径（相对 `<ROOT>`） |
|---|---|---|---|
| audit_protocol.md (v0.1) | artifact `4bbd6664-d9a5-48dd-94a2-4dfdbad0154d` | 0.02 MB | `reports\protocol\audit_protocol_v0.1.md` |
| audit_protocol.md (v0.2) | artifact `eb243319-9ed6-4b10-b801-268ce519def2` | 0.03 MB | `reports\protocol\audit_protocol_v0.2.md` |
| audit_protocol.md (v0.3) | artifact `2e5d6518-4b6f-45bd-8144-02a0b7ec3bc9` | 0.03 MB | `reports\protocol\audit_protocol_v0.3.md` |
| novelty_review.md | artifact `5ce88d31-80b5-4333-911e-ca46d3603527` | 0.03 MB | `reports\novelty_review.md` |
| p0_p1_report.md | artifact `e6a80d5d-c417-4d83-8d84-96b605798854` | 0.01 MB | `reports\p0_p1_report_v0.2.md` |
| p0_v03_report.md | artifact `9b7e8c9b-623a-4e20-a1bf-1260a0f54b86` | 0.01 MB | `reports\p0_report_v0.3.md` |
| audit_lib.py | artifact `0daac09c-d3e1-4354-989e-27320d9f6bfd` | 0.01 MB | `code\audit_lib.py` |
| p1_run.py | artifact `29682c39-c6ab-4592-a233-8194d1cb4033` | 0.00 MB | `code\p1_run.py` |
| p0_units_E.csv | artifact `03245c8d-1188-44a1-9fee-1b94ca38b30f` | 0.04 MB | `data\processed\v0.2_gears_k562ess\p0_units_E.csv` |
| p0v3_units_all_targets.csv | artifact `5942576d-89ca-4b79-b0b9-b5809d63167e` | 0.41 MB | `data\processed\v0.3_figshare_k562ess\p0v3_units_all_targets.csv` |
| go_ase_embeddings.npz | artifact `4047de68-e5f3-4576-8aa4-9770579b2123` | 1.34 MB | `data\processed\priors\go_ase_embeddings.npz` |
| p1_results.parquet | artifact `518b28b7-a70a-487c-9988-ee92bcd6569d` | 0.04 MB | `results\p1_v0.2\p1_results.parquet` |
| p1_planted.parquet | artifact `7eced9d3-78d7-4921-8e3d-30253e669115` | 0.04 MB | `results\p1_v0.2\p1_planted_tau2.parquet` |
| p1_planted_tau6.parquet | artifact `59f8e58d-1d3a-4076-8007-5dee572e5f8f` | 0.04 MB | `results\p1_v0.2\p1_planted_tau6.parquet` |
| p1_permuted.parquet | artifact `0083927e-176f-4b1b-908a-c4300f9d3d5e` | 0.04 MB | `results\p1_v0.2\p1_permuted.parquet` |
| p1_selfcheck_marginal.csv | artifact `ce6b8a5e-b90a-455e-966b-a384af6d38d4` | 0.00 MB | `results\p1_v0.2\p1_selfcheck_marginal.csv` |
| p0v3_power_trend.csv | artifact `7736bf61-b92d-4c4d-a631-367b426ce967` | 0.00 MB | `results\p0_v0.3\p0v3_power_trend.csv` |
| p1_run_manifest.json | artifact `70d71958-5397-4715-bd55-389fd72cf94e` | 0.00 MB | `manifests\p1_run_manifest_v0.2.json` |
| data_manifest_v03.json | artifact `a63c2f0f-7391-47e6-9940-a2015af4c1e3` | 0.00 MB | `manifests\data_manifest_v0.3.json` |

### 3b. 原始数据（原工作区 `data_raw\` 中的文件已丢失，需要重新下载）

| 文件 | 来源 | 大小 | 已记录 SHA256（前 8 位） | 发布方校验 | 下载到 `<DL>\` | 目标路径（相对 `<ROOT>`） |
|---|---|---|---|---|---|---|
| K562_essential_raw_bulk_01.h5ad | Figshare+ 20029387 v1 / file 35773070 | 79.77 MB | `80de95e5…` | Figshare MD5 8321d5d3ffc99db2a5c71edca4189735 | `figshare_20029387_v1\\` | `data\raw\replogle2022_figshare\\` |
| gene2go_all.pkl | Dataverse datafile 6153417 | 9.46 MB | `f145c5e8…` | Dataverse API MD5（执行时获取） | `dataverse_gears\\` | `data\raw\gears_priors\\` |
| essential_all_data_pert_genes.pkl | Dataverse datafile 6934320 | 0.56 MB | `46c3dfe3…` | Dataverse API MD5（执行时获取） | `dataverse_gears\\` | `data\raw\gears_priors\\` |
| go_essential_all.tar.gz | Dataverse datafile 6934319 | 60.65 MB | `98a14a60…` | Dataverse API MD5（执行时获取） | `dataverse_gears\\` | `data\raw\gears_priors\\` |
| replogle_k562_essential.zip | Dataverse datafile 7458695（GEARS 版，v0.2 用） | 669.70 MB | `e769c89a…` | Dataverse API MD5（执行时获取） | `dataverse_gears\\` | `data\raw\gears_k562_essential\\` |
| perturb_processed.h5ad（由上面的 zip 解压） | — | 2.90 GB | `550fde89…` | 无（只能比对已记录的 SHA256） | `（解压到 G:）` | `data\raw\gears_k562_essential\\` |
| 9606.protein.info.v12.0.txt.gz | STRING v12.0 | 1.97 MB | `144de4b0…` | 无公布校验和 | `string_v12.0\\` | `data\raw\string_v12.0\\` |
| 9606.protein.links.detailed.v12.0.txt.gz | STRING v12.0 | 139.63 MB | `28a2ea3b…` | 无公布校验和 | `string_v12.0\\` | `data\raw\string_v12.0\\` |
| linear_perturbation_prediction-Paper-publication.zip | Zenodo 16092690（Ahlmann-Eltze 代码） | 25.88 MB | `2cdd871c…` | Zenodo API MD5（执行时获取） | `zenodo_16092690\\` | `data\raw\ahlmann_eltze_2025\\` |
| single_perturbation_jobs_stats.tsv | Zenodo 16092690 | 0.07 MB | `395aa602…` | Zenodo API MD5（执行时获取） | `zenodo_16092690\\` | `data\raw\ahlmann_eltze_2025\\` |
| cell-gears-0.1.2.tar.gz | PyPI sdist | — | `未记录` | PyPI 公布的 SHA256 | `pypi\\` | `data\raw\software\\` |

上表中的“已记录 SHA256”都取自 `data_manifest_v0.3.json`，比对时用完整的 64 位哈希。

### 3c. 中间数据与划分（已丢失，需按原逻辑重新生成）

| 内容 | 说明 | 重建方式 | 目标路径（相对 `<ROOT>`） |
|---|---|---|---|
| pseudobulk_mean.npy、conditions.csv、readout_genes.csv | GEARS 版 K562 的 pseudobulk 与索引（v0.2 P1 的输入） | 按原代码从 perturb_processed.h5ad 重建 | `data\processed\v0.2_gears_k562ess\` |
| string_expdb_universe.parquet | STRING 的 experiments+databases 子网（v0.2 宇宙） | 按原代码重建；v0.3 宇宙另生成一份 | `data\processed\priors\` |
| go_essential_all.csv | tar 解压出的 GO 边表（约 1,210 万行） | 解压 | `data\processed\priors\（或只保留在 G:）` |
| splits\A_000–A_002.json | v0.2 P1 的只读划分 | 用固定种子 20261001+i 重新生成，并与 p1_results.parquet 里各次划分的测试基因集逐一比对 | `splits\v0.2\` |
| gears_src\ | GEARS 0.1.2 源码（用于核对行号） | 由上面的 sdist 解压 | `data\raw\software\cell-gears-0.1.2\` |
| P0 代码（原来写在 notebook 单元里） | v0.2 与 v0.3 的 P0 构建逻辑 | 整理成 code\p0_build.py，逻辑逐行照搬，不改算法 | `code\` |

Ensembl 坐标（release 116）不重新查询。`p0v3_units_all_targets.csv` 已经保存了每个靶基因的 chr、strand、start、end 和 tss_group；若重新查询，Ensembl 一旦发布新版本，结果就可能变化，反而破坏可复现性。

## 4. 执行步骤

1. **申请授权**：申请 `<ROOT>` 和 `<DL>` 的读写权限，检查两个盘的剩余空间，建立第 2 节的目录结构。
2. **复制 artifact**：把 3a 中的 19 个文件从 artifact 存储只读复制到目标路径。对每个文件，计算源文件（`host.artifact_path`）和目标文件的 SHA256，两者必须一致。
3. **重新下载原始数据**：下载到 `<DL>`，逐个计算 SHA256 和 MD5。
   - SHA256 与 `data_manifest_v0.3.json` 的记录比对，**必须一致**；
   - Figshare 文件还要比对 Figshare 公布的 MD5；Dataverse 和 Zenodo 文件比对其 API 返回的 MD5；PyPI 文件比对其公布的 SHA256。
   - 任何一项不一致，就停止，并报告是哪个文件、哪个哈希不一致，**不用不一致的文件继续往下做**。
   - 全部一致后，复制到 `<ROOT>\data\raw\...`，再比对一次 SHA256，并把 `data\raw\` 下的文件设为只读。
4. **改代码路径**：
   - 新增 `code/paths.py`：从 `<ROOT>/config.toml` 读取 `project_root` 和 `download_root`；若设置了环境变量 `PERTAUDIT_ROOT` 或 `PERTAUDIT_DL`，以环境变量为准；
   - 所有脚本里的路径改为 `ROOT / "data" / "raw" / ...` 这样的形式，**代码里不出现任何盘符**；
   - 两个 manifest 中的文件路径改为相对 `<ROOT>` 的路径。原文件保留一份 `.orig` 备份，改动后的文件另行计算哈希；
   - `p1_run.py` 的划分输出目录改为 `ROOT/"splits"/<protocol_version>`。
5. **整理 P0 代码**：把原先写在 notebook 单元里的 P0 逻辑原样整理成 `code/p0_build.py`。这一步只改组织方式，不改任何算法和参数。
6. **重新生成中间数据和划分**：按 3c 重建。`splits\v0.2\` 中的文件生成后设为只读，然后核验：只读属性确实生效；尝试写入时会报错。
7. **校验表**：生成 `manifests/migration_checksums.csv`，每个文件一行，列包括：
   - 文件、来源类型（artifact / 下载 / 重建）、大小；
   - 源 SHA256、目标 SHA256、是否一致；
   - 已记录 SHA256、是否一致；
   - 发布方校验和、是否一致；
   - 是否只读。
8. **冒烟测试**：见第 5 节。
9. **生成 `reports/project_status.md`**，同时保存为 artifact，然后停下来报告，不执行 P1。

## 5. 冒烟测试：证明迁移后结果与迁移前完全一致

`code/smoke_test.py` 只从 `<ROOT>` 读取数据，依次检查下列各项，全部通过才算迁移成功：

1. **v0.3 的 P0**：用新路径下的原始文件重新跑 P0 构建，要求：
   - E = 1,786，TSS 组 = 1,721，始终训练集 = 257，因 TSS 冲突剔除的 = 14；
   - 每次划分的 fold 大小与已记录值一致（种子 20261001–20261003 下为 1,074/178/267/267、1,075/177/267/267、1,074/178/267/267）；
   - 重新生成的单元表与 `p0v3_units_all_targets.csv` 的副本逐列完全相同（`pandas.testing.assert_frame_equal`）。
2. **v0.2 的 P0**：E = 406、TSS 组 = 402，重新生成的单元表与 `p0_units_E.csv` 的副本完全相同。
3. **v0.2 的划分**：重新生成的 `A_000–A_002`，其测试基因集与 `p1_results.parquet` 中对应划分的测试基因集完全相同。这是唯一能证明“划分与迁移前一致”的方法，因为原 JSON 文件已经丢失。
4. **v0.2 的 P1 数值复现**：在新路径下重跑 v0.2 的 split 0，得到的 score_value、qhat、covered 与 `p1_results.parquet` 中 split 0 的对应行一致（数值容差 1e-10；覆盖指示要求完全相同）。
5. **GO 嵌入**：重新计算的 ASE 嵌入与 `go_ase_embeddings.npz` 的副本在逐列符号对齐后一致（容差 1e-8）。特征向量的符号本身不确定，所以必须先对齐符号再比较。

## 6. 迁移之后的写入规则

- 所有新的输出直接写到 `<ROOT>` 下的相应目录，所有新的下载先落到 `<DL>`，经校验后再复制到 `<ROOT>\data\raw\`。
- 需要在应用里展示的结果（报告、图、表），会额外复制一份到会话工作区再保存为 artifact。原因是 artifact 只能从工作区，或从已注册的本地代码库保存。
- 如果你希望 artifact 直接从 `<ROOT>` 保存，可以把 `<ROOT>` 注册为本地代码库根目录。但注册时应用会在 `<ROOT>` 下建一个 `.venv`，所以默认不这样做。

## 7. 需要你确认的事项

1. 是否同意按本计划执行，包括重新下载 11 个原始文件。共约 0.99 GB；其中 GEARS 版的 zip 有 0.67 GB，它只用于复现 v0.2，如果你不需要复现 v0.2，可以省掉。
2. F: 和 G: 各有多少可用空间。
3. 原始数据是否**同时**存放在 `<DL>` 和 `<ROOT>\data\raw\` 两处（本计划的默认做法）。另一种做法是 `data\raw\` 只放指向 `<DL>` 的清单文件。以后有 10.7 GB 的单细胞文件，这个选择会影响空间。
4. 你提到“当前方案版本 v0.4”，但 v0.4 的具体内容（例如上次 P0 报告末尾提出的四个待定事项）你还没有决定。我计划在 `project_status.md` 中写“v0.4：待你确定上述四项后定稿”，不擅自替你做决定。
