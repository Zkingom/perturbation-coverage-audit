"""对比 GEARS 原版 evaluate 与 mem_patch.evaluate（GPU 路径和 numpy 路径）：compute_metrics 结果必须逐位相同。"""
import numpy as np, torch
from torch_geometric.data import Data, DataLoader
import gears.inference as gi
import mem_patch

G = 600
rng = np.random.default_rng(0)
graphs = []
for pi in range(50):
    for c in range(16):
        graphs.append(Data(x=torch.randn(G, 1), y=torch.randn(1, G),
                           de_idx=rng.choice(G, 20, replace=False), pert=f"P{pi}+ctrl"))


class M(torch.nn.Module):
    def __init__(self):
        super().__init__(); self.l = torch.nn.Linear(1, 1)
    def forward(self, b):
        return self.l(b.x).view(-1, G) * 1.7 + 0.3


m = M().cuda()
ok = True
for name, thr, dl in [("gpu", 0, False), ("numpy", 10**18, False), ("gpu-droplast", 0, True), ("numpy-droplast", 10**18, True)]:
    torch.manual_seed(0); r1 = gi.evaluate(DataLoader(graphs, batch_size=33, shuffle=True, drop_last=dl), m, False, "cuda")
    mem_patch.GPU_BYTES_THRESHOLD = thr
    torch.manual_seed(0); r2 = mem_patch.evaluate(DataLoader(graphs, batch_size=33, shuffle=True, drop_last=dl), m, False, "cuda")
    a, ap = gi.compute_metrics(r1); b, bp = gi.compute_metrics(r2)
    same = (list(r1["pert_cat"]) == list(r2["pert_cat"])
            and all(np.array_equal(np.asarray(a[k]), np.asarray(b[k])) for k in a)
            and all(ap[p][k] == bp[p][k] for p in ap for k in ap[p])
            and np.array_equal(r1["pred_de"], r2["pred_de"]) and np.array_equal(r1["truth_de"], r2["truth_de"]))
    idx = np.arange(0, len(r1["pert_cat"]), 7)
    same = same and len(r1["pert_cat"]) == len(r2["pred"])
    same = same and np.array_equal(r1["pred"][idx], r2["pred"][idx]) and np.array_equal(r1["truth"][idx], r2["truth"][idx])
    print(name, "IDENTICAL" if same else "DIFFERENT", "on_gpu" if isinstance(r2["pred"], mem_patch._GpuRows) else "numpy")
    ok &= same
print("TEST PASS" if ok else "TEST FAIL")
