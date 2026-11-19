"""
precompute_nets.py SPLIT — 在轻量进程里为某个 split 预先算好 GEARS 的共表达网络（以及一次性的 GO 网络），
存到 /root/autodl-tmp/gears/netcache，供 mem_patch.py 读取。
- 加载 adata 的方式与正式运行完全相同（PertData.load），只是不把 11GB 的 cell_graphs.pkl 读进内存
  （这一步只影响图数据集，不影响 adata 和网络计算）。
- 网络用 GEARS 原函数 gears.utils.get_similarity_network、与 GEARS.model_initialize 默认值相同的参数计算：
  co-express threshold=0.4, k=20；GO k=20；split='custom', seed=1, train_gene_set_size=0.75（与 run_split 一致）。
"""
import sys, os, json, shutil, pickle as _real_pickle
from pathlib import Path

i = int(sys.argv[1])
B = Path("/root/autodl-tmp/gears/gears_bundle_v07")
WORK = B / "work" / "gears_data"
DS = "rpe1_cap16"
CACHE = Path("/root/autodl-tmp/gears/netcache"); CACHE.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(B / "code"))

import gears.pertdata as gp


class _SkipCellGraphs:
    """只拦截 cell_graphs.pkl 的读取（返回空 dict），其余 pickle 调用原样转发。"""
    def __getattr__(self, name):
        return getattr(_real_pickle, name)

    @staticmethod
    def load(f, *a, **kw):
        if str(getattr(f, "name", "")).endswith("cell_graphs.pkl"):
            return {}
        return _real_pickle.load(f, *a, **kw)


gp.pickle = _SkipCellGraphs()

from gears import PertData
from gears.utils import get_similarity_network as orig
from mem_patch import coexp_cache_name, go_cache_name  # 同一套缓存文件名

s2c = json.loads((B / "splits" / f"set2conditions_{i:02d}.json").read_text())
split, seed, tgs, thr, k = "custom", 1, 0.75, 0.4, 20
target = CACHE / coexp_cache_name(s2c, split, seed, tgs, thr, k)

if not target.exists():
    name = f"pre_s{i:02d}"
    sdir = WORK / name
    shutil.rmtree(sdir, ignore_errors=True)
    (sdir / "data_pyg").mkdir(parents=True)
    os.symlink(WORK / DS / "perturb_processed.h5ad", sdir / "perturb_processed.h5ad")
    os.symlink(WORK / DS / "data_pyg" / "cell_graphs.pkl", sdir / "data_pyg" / "cell_graphs.pkl")
    pdat = PertData(str(WORK))
    pdat.load(data_path=str(sdir))
    orig("co-express", adata=pdat.adata, threshold=thr, k=k, data_path=str(WORK), data_name=name,
         split=split, seed=seed, train_gene_set_size=tgs, set2conditions=s2c)
    produced = sdir / f"{split}_{seed}_{tgs}_{thr}_{k}_co_expression_network.csv"
    tmp = target.with_suffix(".tmp")
    shutil.copyfile(produced, tmp); os.replace(tmp, target)
    shutil.rmtree(sdir, ignore_errors=True)
    print(f"[precompute] split {i:02d} co-express -> {target.name}", flush=True)
    go_adata = pdat.adata
else:
    print(f"[precompute] split {i:02d} co-express cache exists", flush=True)
    go_adata = None

gof = CACHE / go_cache_name(k)
if not gof.exists():
    df = orig("go", adata=go_adata, threshold=thr, k=k, data_path=str(WORK), data_name="go",
              split=split, seed=seed, train_gene_set_size=tgs, set2conditions=s2c, default_pert_graph=True)
    tmp = gof.with_suffix(".tmp"); df.to_pickle(tmp); os.replace(tmp, gof)
    print(f"[precompute] GO -> {gof.name}", flush=True)
