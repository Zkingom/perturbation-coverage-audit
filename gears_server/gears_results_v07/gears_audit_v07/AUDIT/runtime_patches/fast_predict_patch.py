"""
fast_predict_patch.py  —  GEARS 0.1.2 预测加速补丁（结果与原版逐位一致）

原版 create_cell_graph_dataset_for_prediction 对每个扰动都会：
  1) 在 2.6GB 的 AnnData 视图上做一次 ctrl_adata[idx, :] 切片（很慢，是主要瓶颈）
  2) 把 300 个小图逐个 .to(device)（30 万次小拷贝）
本补丁：ctrl 表达矩阵只稠密化一次并缓存；每个扰动仍调用同一个
np.random.randint(0, n_ctrl, 300)，随机数序列与原版完全相同，所以采样到的细胞、
预测结果都不变。图先留在 CPU，由 GEARS.predict 里原有的 batch.to(device) 一次搬上 GPU。

用法：在 gears_server_run.py 里 `from gears import GEARS` 那一行之后加一行
    import fast_predict_patch  # noqa
（把本文件放在 code/ 目录下，和 gears_server_run.py 同级）

自检：python fast_predict_patch.py --selftest  （不需要，可选；见 verify_equivalence）
"""
import numpy as np
import torch
from torch_geometric.data import Data

import gears.gears as _gg
import gears.utils as _gu

_orig = _gu.create_cell_graph_dataset_for_prediction
_cache = {"key": None, "X": None}


def _dense_ctrl(ctrl_adata):
    # predict() 每次调用都会重新生成 ctrl_adata 视图；用 (id(父对象), n_obs) 作缓存键
    key = (id(getattr(ctrl_adata, "_adata_ref", ctrl_adata)), ctrl_adata.n_obs, ctrl_adata.n_vars)
    if _cache["key"] != key:
        X = ctrl_adata.X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        _cache["X"] = np.ascontiguousarray(X)
        _cache["key"] = key
    return _cache["X"]


def fast_create_cell_graph_dataset_for_prediction(pert_gene, ctrl_adata, gene_names,
                                                  device, num_samples=300):
    pert_idx = [np.where(p == np.array(gene_names))[0][0] for p in pert_gene]
    Xc = _dense_ctrl(ctrl_adata)
    idx = np.random.randint(0, len(ctrl_adata), num_samples)  # 与原版同一次 RNG 调用
    Xs = Xc[idx]
    return [Data(x=torch.Tensor(X).T, pert_idx=pert_idx, pert=pert_gene) for X in Xs]


def verify_equivalence(ctrl_adata, gene_names, perts, device="cpu", seed=123):
    """对给定的几个扰动，比较原版和补丁版生成的输入是否完全一致。"""
    for p in perts:
        np.random.seed(seed); a = _orig(p, ctrl_adata, gene_names, device)
        np.random.seed(seed); b = fast_create_cell_graph_dataset_for_prediction(p, ctrl_adata, gene_names, device)
        assert len(a) == len(b)
        for ga, gb in zip(a, b):
            assert torch.equal(ga.x.cpu(), gb.x.cpu()) and ga.pert_idx == gb.pert_idx
    return True


# GEARS.predict 在 gears/gears.py 的命名空间里查找这个函数，所以要替换那里的引用
_gg.create_cell_graph_dataset_for_prediction = fast_create_cell_graph_dataset_for_prediction
print("[fast_predict_patch] GEARS predict 加速补丁已启用", flush=True)
