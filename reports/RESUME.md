> **2026-10-10 更新：v0.7 草案（未冻结）在 reports/protocol/audit_protocol_v0.7_draft.md；开发集补充结果在 reports/d1/d1_v0.7dev_report.md；下一步是用户审阅草案。**

> **2026-10-09 更新：D1 第 2–4 阶段已完成，见 reports/d1/d1_report_v0.6.2.md；下一步是用户就报告第 6 节的事项做决定，然后起草 v0.7。**

# 快速恢复说明（写于 2026-10-08T14:32:06Z，会话暂停时）

## 当前状态
- **方案**：v0.6.2 已冻结。文件为 `reports/protocol/audit_protocol_v0.6.2.md`，冻结清单为 `manifests/freeze_v0.6.2.json`，两者均只读。
- **D1a 已完成**：开发集分半 pseudobulk 在 `data/processed/v0.6/k562ess/`，只读。
- **D1 第 1 阶段已完成**：21 个缓存在 `results/d1_v0.6.1/cache/split_00..20.pkl`；划分文件在 `splits/v0.6.1_dev/`，只读。
- **尚未运行**：第 2–4 阶段（D1b–D1d），即覆盖率、选择规则、混杂分析和增量分析。开发集上**还没有任何覆盖结果**。
- **确认集**（RPE1、GWPS）只做过大小和 MD5 核对，没有读取任何数值。
- **没有任何后台任务在运行。**

## 下次启动的步骤
1. **检查环境**
   - 本机可用内存：第 1 阶段实测进程峰值为 2.14 GB，第 2–4 阶段尚未实测，建议可用内存 ≥ 3.5 GB。
   - 核对冻结清单：`freeze_v0.6.2.json` 中 code_sha256 列出的各文件，哈希都应一致。
2. **写第 2–4 阶段的运行入口**
   - `d1_run_v061.py` 目前只有 `pseudobulk` 和 `stage1` 两种模式。需要新增一个 `stages234` 模式，调用 `pipeline_v06.run_dataset(load_dev_v06(PB_DIR), DEV_SEED_BASE, results/d1_v0.6.1)`。它会跳过已缓存的划分，并带内存监控。
   - 这会改变 `d1_run_v061.py` 的哈希，所以要在冻结清单中追加改动记录，原哈希保留。这只是新增运行入口，不改变任何算法。
   - 先计时：只在 1 次划分上运行一遍第 2 阶段，以估算总耗时。开发集的校准集约有 267 个单元，比合成数据大，M3 和宽度匹配的 LP 会更慢。
3. **执行**：由 Sonnet 5.5 子任务执行并监控，规则不变（报错就停，不做变通处理）；结果解读由 Opus 5.5 完成。
4. **D1 的必报项**：
   - 效应强度与 miscoverage 的关系，以及混杂分析（含分半版本）；
   - ĥ 的可预测性和增量价值，即 (a)(b)(c) 三项；
   - 宽度匹配下各方法的比较：oracle、M1/M2、M_dist、M0、M3、Mcomb、Mscore、M_perm；
   - **Mdist-N 和 Mdist-M 相对 M0 在 D_GO 和 D_STRING 分层上的覆盖改进**。如果几乎没有改进，在 v0.7 之前提出更强的类 PRESCRIBE 变体，由用户决定。
   - 照实报告：d 在 21/21 次划分中都取到网格上界 100，K 在 16/21 次取到上界 50。
5. **之后**：起草 v0.7，包含 `reports/protocol/v0.7_draft_notes.md` 中已确定的各项。v0.7 必须经用户审阅后才能冻结。

## 环境与路径
- Python 环境 `pert-audit`，位于 `%USERPROFILE%\.claude-science\conda\envs\pert-audit`。用子进程运行时，PATH 需加上该环境的根目录、`Library\bin` 和 `Scripts`。
- R 环境 `pert-audit-r`（GLMM）；`condconf-ref`（M3 交叉核对用，numpy 1.26）。
- 原始数据：G:\biology（K562 essential、RPE1）；H:\cs+bio_file（GWPS 单细胞文件；沙箱无法访问，需要由用户在本机运行）。
