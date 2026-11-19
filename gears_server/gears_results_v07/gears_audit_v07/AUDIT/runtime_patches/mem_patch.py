"""
mem_patch.py — 让 GEARS 在 15GB 内存的容器里跑得动（结果不变）

问题：一个 split 里，cell_graphs 读进来就占约 11GB；GEARS 初始化时还要现算
  (1) 共表达网络：把全部训练细胞转成稠密矩阵算 Pearson 相关，额外占约 4–7GB
  (2) GO 网络：读 354MB 的 go_essential_all.csv 再 groupby，额外占数 GB
峰值超过容器 15GB 上限 -> 被系统杀掉 (rc=137)。

做法：这两张网络由 precompute_nets.py 在一个只加载 adata 的轻量进程里，用 GEARS
自己的同一个函数、同一组参数提前算好存盘；正式运行时这里直接读盘：
  - 共表达网络：返回 pd.read_csv(缓存文件)。这和 GEARS 自身“缓存命中”分支
    (utils.get_coexpression_network_from_train 里的 `return pd.read_csv(fname)`) 完全相同。
  - GO 网络：返回 pickle 的 DataFrame，逐位一致。
两张网络的计算都不消耗随机数，所以训练的随机数序列也不受影响。
缓存不存在时自动退回 GEARS 原始计算（行为与原版完全一样）。
"""
import os
import hashlib
import pandas as pd
import gears.gears as _gg
import gears.utils as _gu

CACHE = "/root/autodl-tmp/gears/netcache"
_orig = _gu.get_similarity_network


def coexp_cache_name(set2conditions, split, seed, train_gene_set_size, threshold, k):
    h = hashlib.sha1("\n".join(sorted(set2conditions["train"])).encode()).hexdigest()[:16]
    return f"coexp_{h}_{split}_{seed}_{train_gene_set_size}_{threshold}_{k}.csv"


def go_cache_name(k):
    return f"go_default_k{k}.pkl"


def get_similarity_network(network_type, adata, threshold, k, data_path, data_name, split, seed,
                           train_gene_set_size, set2conditions, default_pert_graph=True, pert_list=None):
    if network_type == "co-express":
        f = os.path.join(CACHE, coexp_cache_name(set2conditions, split, seed, train_gene_set_size, threshold, k))
        if os.path.exists(f):
            print(f"[mem_patch] co-express network from cache {os.path.basename(f)}", flush=True)
            return pd.read_csv(f)
    elif network_type == "go" and default_pert_graph:
        f = os.path.join(CACHE, go_cache_name(k))
        if os.path.exists(f):
            print(f"[mem_patch] GO network from cache {os.path.basename(f)}", flush=True)
            return pd.read_pickle(f)
    print(f"[mem_patch] no cache for {network_type}, computing with original GEARS code", flush=True)
    return _orig(network_type, adata, threshold, k, data_path, data_name, split, seed,
                 train_gene_set_size, set2conditions, default_pert_graph=default_pert_graph, pert_list=pert_list)


_gg.get_similarity_network = get_similarity_network
print("[mem_patch] GEARS similarity-network cache enabled", flush=True)

# ---- (3) evaluate(): 结果逐位相同，但内存占用更小 ----
# 原版把每个 batch 的行张量都存进 list，最后 torch.stack 再复制一份；对整个训练集评估时，
# pred+truth 要占 3GB 以上的 CPU 内存，在 15GB 容器里会爆。这里：
#  - 小的评估集(val/test)：预分配 numpy 数组逐 batch 填入，结果与原版完全相同；
#  - 大的评估集(训练集，每个 epoch 只用于打印 train MSE)：把 pred/truth 留在 GPU（3060 有 12GB 显存空着），
#    compute_metrics 里按扰动取行 results['pred'][p_idx] 时再拷回 CPU。取回的行数值与原版相同，
#    np.mean 作用在同样形状的连续数组上，所以 train 指标也是逐位相同。
import numpy as np
import torch


class _GpuRows:
    """把 GPU 上的 (N, G) 张量包装成只支持“按行号取行”的数组：x[idx] -> numpy。"""
    def __init__(self, t):
        self.t = t

    def __getitem__(self, idx):
        if isinstance(idx, np.ndarray):
            idx = torch.as_tensor(idx, device=self.t.device)
        return self.t[idx].cpu().numpy()

    @property
    def shape(self):
        return tuple(self.t.shape)

    def __len__(self):
        return self.t.shape[0]


GPU_BYTES_THRESHOLD = 1_000_000_000  # pred+truth 超过约 1GB 才放 GPU


def evaluate(loader, model, uncertainty, device):
    model.eval()
    model.to(device)
    n = len(loader.dataset)
    pert_cat, pred_de, truth_de, logvar = [], [], [], []
    pred = truth = None
    on_gpu = False
    k = 0
    for itr, batch in enumerate(loader):
        batch.to(device)
        pert_cat.extend(batch.pert)
        with torch.no_grad():
            if uncertainty:
                p, unc = model(batch)
                logvar.extend(unc.cpu())
            else:
                p = model(batch)
            t = batch.y
            if pred is None:
                need = 2 * n * int(np.prod(p.shape[1:])) * p.element_size()
                on_gpu = (p.is_cuda and need > GPU_BYTES_THRESHOLD
                          and torch.cuda.mem_get_info(p.device)[0] > need + 1_000_000_000)
                if on_gpu:
                    pred = torch.empty((n,) + tuple(p.shape[1:]), dtype=p.dtype, device=p.device)
                    truth = torch.empty((n,) + tuple(t.shape[1:]), dtype=t.dtype, device=t.device)
                else:
                    pred = np.empty((n,) + tuple(p.shape[1:]), dtype=p.cpu().numpy().dtype)
                    truth = np.empty((n,) + tuple(t.shape[1:]), dtype=t.cpu().numpy().dtype)
            b = p.shape[0]
            if on_gpu:
                pred[k:k + b] = p
                truth[k:k + b] = t
            else:
                pred[k:k + b] = p.cpu().numpy()
                truth[k:k + b] = t.cpu().numpy()
            k += b
            for j, de_idx in enumerate(batch.de_idx):
                pred_de.append(p[j, de_idx])
                truth_de.append(t[j, de_idx])
    # 训练集 loader 是 drop_last=True，实际行数 k 可能小于 len(dataset)；与原版一样只保留实际评估过的 k 行
    pred, truth = pred[:k], truth[:k]
    if on_gpu:
        pred, truth = _GpuRows(pred), _GpuRows(truth)
    results = {"pert_cat": np.array(pert_cat), "pred": pred, "truth": truth}
    results["pred_de"] = torch.stack(pred_de).detach().cpu().numpy()
    results["truth_de"] = torch.stack(truth_de).detach().cpu().numpy()
    if uncertainty:
        results["logvar"] = torch.stack(logvar).detach().cpu().numpy()
    return results


_gg.evaluate = evaluate
print("[mem_patch] lean evaluate enabled", flush=True)

# ---- (4) deeper_analysis / non_dropout_analysis：去掉两处“把整个 adata.X 转成稠密矩阵/整份复制”的写法 ----
# GEARS 在测试结束后调用这两个分析函数（结果只用于 wandb 日志，本冻结流程不保存它们）。原代码里：
#   mean_expression = ....reshape(..., adata.X.toarray().shape[1])   # 为了取列数把整个 X 稠密化(数 GB)
#   for i, j in conditions2index.items(): np.mean(adata.X[j], axis=0)  # adata 是 view，每次访问 .X 都整份复制稀疏矩阵
# 在 15GB 容器里，训练 20 个 epoch 全部跑完后正是在这里被 OOM 杀掉（rc=137）。
# 改法是对函数源码做两处字面替换（替换次数有断言）：
#   adata.X.toarray().shape[1]  ->  adata.n_vars                  （同一个数）
#   np.mean(adata.X[j], axis = 0) -> np.mean(_X[j], axis = 0)，_X 在循环前取一次：
#     若 adata 是“覆盖父对象全部行列、顺序不变”的 view（本流程 GO 过滤没有删掉任何扰动，就是这种情况），
#     _X 直接用父对象的 X（数据完全相同，不复制）；否则退回 adata.X（与原版一致）。
import inspect
import gears.inference as _gi


def _fast_X(adata):
    if getattr(adata, "is_view", False):
        ref = adata._adata_ref
        oidx, vidx = adata._oidx, adata._vidx

        def _full(ix, n):
            if isinstance(ix, slice):
                return ix == slice(None, None, None) or (ix.start in (None, 0) and ix.stop in (None, n) and ix.step in (None, 1))
            ix = np.asarray(ix)
            if ix.dtype == bool:
                return ix.shape[0] == n and bool(ix.all())
            return ix.shape[0] == n and bool((ix == np.arange(n)).all())

        if _full(oidx, ref.n_obs) and _full(vidx, ref.n_vars):
            return ref.X
    return adata.X


def _rewrite(fn):
    src = inspect.getsource(fn)
    a, b = "adata.X.toarray().shape[1]", "np.mean(adata.X[j], axis = 0)"
    c = "    condition2mean_expression = {}\n"
    assert src.count(a) == 1 and src.count(b) == 1 and src.count(c) == 1, f"unexpected source for {fn.__name__}"
    src = src.replace(a, "adata.n_vars").replace(b, "np.mean(_X[j], axis = 0)").replace(c, c + "    _X = _fast_X(adata)\n")
    ns = dict(_gi.__dict__); ns["_fast_X"] = _fast_X
    exec(compile(src, f"<mem_patch:{fn.__name__}>", "exec"), ns)
    return ns[fn.__name__]


_gg.deeper_analysis = _rewrite(_gi.deeper_analysis)
_gg.non_dropout_analysis = _rewrite(_gi.non_dropout_analysis)
print("[mem_patch] deeper_analysis / non_dropout_analysis without full densify", flush=True)

# ---- (5) 释放 glibc 堆碎片：在每次 evaluate / compute_metrics / 两个分析函数 / predict 之后
#          gc.collect() + malloc_trim(0)，把已经释放、但 glibc 还占着不还给系统的内存交还。
#          只影响进程占用的物理内存，不影响任何计算结果。
import ctypes as _ct, gc as _gc
try:
    _libc = _ct.CDLL("libc.so.6")
except OSError:
    _libc = None


def _trim():
    _gc.collect()
    if _libc is not None:
        _libc.malloc_trim(0)


def _with_trim(f):
    def w(*a, **k):
        r = f(*a, **k)
        _trim()
        return r
    w.__wrapped__ = f
    return w


for _n in ["evaluate", "compute_metrics", "deeper_analysis", "non_dropout_analysis"]:
    setattr(_gg, _n, _with_trim(getattr(_gg, _n)))
_gg.GEARS.predict = _with_trim(_gg.GEARS.predict)
print("[mem_patch] malloc_trim after evaluate/metrics/analyses/predict enabled", flush=True)
