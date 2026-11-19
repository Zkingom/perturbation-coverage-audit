"""postchecks.py — 全部 split 跑完后（GPU 进程已退出、内存空出来）再做的两项核查，只读，输出到 gears/POSTCHECK/。
(A) fast_predict_patch：在真实 ctrl_adata 上，对 split 00 的若干测试扰动比较原版/补丁版生成的预测输入（同一随机种子）是否逐位相同。
(B) 共表达网络对浮点舍入的敏感性（split 00）：用与 GEARS 完全相同的步骤，但全程 float64 重算一遍，
    与缓存里的 float32 结果比较边集；并统计“选中边离阈值 0.4 的最小距离”和“每个基因第 k+1 与第 k+2 名的最小间隔”。
    若这两个间隔远大于 float32 舍入误差（~1e-6），说明 A800/3060 之间的舍入差异不可能改变边集。
"""
import sys, os, json, shutil, pickle as _real_pickle
from pathlib import Path
import numpy as np, pandas as pd

B = Path("/root/autodl-tmp/gears/gears_bundle_v07"); WORK = B / "work" / "gears_data"; DS = "rpe1_cap16"
G = Path("/root/autodl-tmp/gears"); OUT = G / "POSTCHECK"; OUT.mkdir(exist_ok=True)
import gears.pertdata as gp


class _Skip:
    def __getattr__(self, n): return getattr(_real_pickle, n)
    @staticmethod
    def load(f, *a, **k):
        return {} if str(getattr(f, "name", "")).endswith("cell_graphs.pkl") else _real_pickle.load(f, *a, **k)


gp.pickle = _Skip()
from gears import PertData
import fast_predict_patch as fpp
from mem_patch import coexp_cache_name

name = "post_s00"; sdir = WORK / name
shutil.rmtree(sdir, ignore_errors=True); (sdir / "data_pyg").mkdir(parents=True)
os.symlink(WORK / DS / "perturb_processed.h5ad", sdir / "perturb_processed.h5ad")
os.symlink(WORK / DS / "data_pyg" / "cell_graphs.pkl", sdir / "data_pyg" / "cell_graphs.pkl")
pdat = PertData(str(WORK)); pdat.load(data_path=str(sdir))
adata = pdat.adata
s2c = json.loads((B / "splits" / "set2conditions_00.json").read_text())
rep = []

# (A)
ctrl = adata[adata.obs["condition"] == "ctrl"]
perts = [[c.split("+")[0]] for c in s2c["test"][:8]]
try:
    ok = fpp.verify_equivalence(ctrl, list(pdat.pert_names), perts, device="cpu", seed=20270701)
    rep.append(f"(A) fast_predict inputs identical to original for {len(perts)} real test perturbations: {ok}")
except AssertionError as e:
    rep.append(f"(A) fast_predict MISMATCH: {e!r}")

# (B)
k, thr = 20, 0.4
gene_list = list(adata.var.gene_name.values)
X = adata.X
X_tr = X[np.isin(adata.obs.condition, [i for i in s2c["train"] if "ctrl" in i])]
X_tr = X_tr.toarray().astype(np.float64)
xv = X_tr - X_tr.mean(axis=0); del X_tr
ss = (xv * xv).sum(axis=0)
out = np.matmul(xv.T, xv) / np.sqrt(np.outer(ss, ss)); del xv
out = np.maximum(np.minimum(out, 1.0), -1.0); out[np.isnan(out)] = 0; out = np.abs(out)
srt = np.sort(out, axis=1)
gap = srt[:, -(k + 1)] - srt[:, -(k + 2)]          # 第 k+1 名与第 k+2 名的间隔（决定 top-(k+1) 选谁）
idx = np.argsort(out)[:, -(k + 1):]
ref = set()
for i in range(idx.shape[0]):
    for j in range(idx.shape[1]):
        v = out[i, idx[i, j]]
        if v > thr:
            ref.add((gene_list[idx[i, j]], gene_list[i]))
sel = srt[:, -(k + 1):]
dist_thr = np.abs(sel - thr).min()
cache = pd.read_csv(G / "netcache" / coexp_cache_name(s2c, "custom", 1, 0.75, thr, k))
got = set(zip(cache.source, cache.target))
rep.append(f"(B) split00 float64 edges={len(ref)} cached(float32 on this machine) edges={len(got)} "
           f"identical_edge_set={ref == got} only_in_f64={len(ref - got)} only_in_cache={len(got - ref)}")
rep.append(f"(B) min gap between rank k+1 and k+2 over genes = {gap.min():.3e}; "
           f"genes with gap<1e-5: {int((gap < 1e-5).sum())}; min |selected-|r| - 0.4| = {dist_thr:.3e}")
(OUT / "postchecks.txt").write_text("\n".join(rep) + "\n")
shutil.rmtree(sdir, ignore_errors=True)
print("\n".join(rep))
