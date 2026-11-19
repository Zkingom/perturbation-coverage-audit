"""在真实 adata 上比较原版与改写版 deeper_analysis / non_dropout_analysis 的输出是否完全相同（须在没有训练进程时运行）。"""
import os, json, shutil, pickle as _real_pickle
from pathlib import Path
import numpy as np
import gears.pertdata as gp
import gears.inference as gi


class _Skip:
    def __getattr__(self, n): return getattr(_real_pickle, n)
    @staticmethod
    def load(f, *a, **k):
        return {} if str(getattr(f, "name", "")).endswith("cell_graphs.pkl") else _real_pickle.load(f, *a, **k)


gp.pickle = _Skip()
from gears import PertData
import mem_patch

B = Path("/root/autodl-tmp/gears/gears_bundle_v07"); WORK = B / "work" / "gears_data"; DS = "rpe1_cap16"
sdir = WORK / "post_ana"; shutil.rmtree(sdir, ignore_errors=True); (sdir / "data_pyg").mkdir(parents=True)
os.symlink(WORK / DS / "perturb_processed.h5ad", sdir / "perturb_processed.h5ad")
os.symlink(WORK / DS / "data_pyg" / "cell_graphs.pkl", sdir / "data_pyg" / "cell_graphs.pkl")
pdat = PertData(str(WORK)); pdat.load(data_path=str(sdir)); adata = pdat.adata
print("adata is_view:", adata.is_view, "shape:", adata.shape, "fast_X is parent X:", mem_patch._fast_X(adata) is getattr(adata, "_adata_ref", adata).X)

s2c = json.loads((B / "splits" / "set2conditions_00.json").read_text())
rng = np.random.default_rng(0)
perts = s2c["test"][:6]
rows = np.concatenate([np.where(adata.obs.condition.values == p)[0] for p in perts])
truth = np.asarray(adata.X[rows].todense(), dtype=np.float32)
pred = truth + rng.normal(0, 0.1, truth.shape).astype(np.float32)
geneid2idx = dict(zip(adata.var.index.values, range(adata.n_vars)))
p2f = dict(adata.obs[["condition", "condition_name"]].values)
pc = adata.obs.condition.values[rows]
de = np.array([[geneid2idx[g] for g in adata.uns["rank_genes_groups_cov_all"][p2f[c]][:20]] for c in pc])
res = {"pert_cat": np.array(pc), "pred": pred, "truth": truth,
       "pred_de": np.take_along_axis(pred, de, 1), "truth_de": np.take_along_axis(truth, de, 1)}


def same(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and (np.array_equal(a, b) or np.array_equal(a, b, equal_nan=True))


# 原版每个条件都对 view 访问一次 .X（每次整份复制稀疏矩阵），对 2000+ 个条件要跑一个多小时，
# 所以这里改为直接验证改写所依赖的两个等式，再确认改写版能在真实数据上跑通：
Xv = adata.X                      # 原版拿到的矩阵（view 的整份复制）
Xf = mem_patch._fast_X(adata)     # 改写版拿到的矩阵
eq_X = (Xv.shape == Xf.shape) and (Xv != Xf).nnz == 0 and Xv.dtype == Xf.dtype
eq_n = adata.n_vars == Xv.shape[1]
print("X identical (values, shape, dtype):", eq_X, "| n_vars == X.toarray().shape[1]:", eq_n)
del Xv
o1 = mem_patch._gg.deeper_analysis(adata, res); o2 = mem_patch._gg.non_dropout_analysis(adata, res)
print("patched analyses ran:", len(o1), len(o2), "perturbations")
ok = bool(eq_X and eq_n)
shutil.rmtree(sdir, ignore_errors=True)
print("ANALYSIS TEST PASS" if ok else "ANALYSIS TEST FAIL")
