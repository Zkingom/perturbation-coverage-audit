# qc_nadig_wessels_v121_v14.py -- data QC before any analysis (no score, no model, no result).
#   python qc_nadig_wessels_v121_v14.py nadig    protocol v1.2.1 section 4, checks Q1-Q5 + library-size records
#   python qc_nadig_wessels_v121_v14.py wessels  protocol v1.4 section 3, gate verification G1-G5 + file identity
# Reads h5py in row blocks (memory bounded); writes results/v1.2/<cell line>/qc_v121.json and
# results/v1.4/wessels/gates_v14.json (create-only: refuses to overwrite).
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import RESULTS  # noqa: E402
from safe_io import create_new  # noqa: E402

DL = Path(r"G:\biology\scperturb_13350497")
GEO = {"hepg2": (145473, 9624), "jurkat": (262956, 8882)}  # GEO header (read 2026-11-04)
N_TARGETS_DESIGN = 2393  # Nadig et al. 2025, main text
BLOCK = 2000


def col(o, k):
    x = o[k]
    if isinstance(x, h5py.Group):
        c = np.array([v.decode() if isinstance(v, bytes) else v for v in x["categories"][:]], dtype=object)
        return c[x["codes"][:]]
    v = x[:]
    return np.array([a.decode() if isinstance(a, bytes) else a for a in v], dtype=object) if v.dtype.kind in "SO" else v


def nadig(cl):
    t0 = time.time()
    with h5py.File(DL / f"NadigOConner2024_{cl}.h5ad", "r") as h:
        X = h["X"]
        n, g = X.shape
        nonint = neg = nnz = 0
        rowsum = np.empty(n)
        for s in range(0, n, BLOCK):
            b = X[s:s + BLOCK].astype(np.float64)
            nzm = b != 0
            nnz += int(nzm.sum())
            neg += int((b < 0).sum())
            nonint += int((np.abs(b - np.round(b)) > 1e-6).sum())
            rowsum[s:s + BLOCK] = b.sum(1)
        o = h["obs"]
        pert = col(o, "perturbation")
        ncounts = col(o, "ncounts").astype(float)
        umi = col(o, "UMI_count").astype(float)
    targets = sorted(set(pert) - {"control"})
    n_ctrl = int((pert == "control").sum())
    gn, gg = GEO[cl]
    q = {
        "Q1_integer_counts": {"nonzero_values": nnz, "non_integer": nonint, "negative": neg, "pass": nonint == 0 and neg == 0},
        "Q2_cells": {"n": int(n), "geo": gn, "ratio": n / gn, "pass": 0.80 <= n / gn <= 1.00},
        "Q3_genes": {"n": int(g), "geo": gg, "ratio": g / gg, "pass": g >= 0.90 * gg},
        "Q4_targets": {"n": len(targets), "design": N_TARGETS_DESIGN, "ratio": len(targets) / N_TARGETS_DESIGN,
                       "pass": len(targets) >= 0.90 * N_TARGETS_DESIGN},
        "Q5_control_cells": {"label": "control", "n": n_ctrl, "pass": n_ctrl >= 1000},
    }
    q["all_pass"] = all(v["pass"] for v in q.values())
    q["library_size_record"] = {
        "ncounts_minus_rowsum_max_abs": float(np.abs(ncounts - rowsum).max()),
        "ncounts_over_rowsum_median": float(np.median(ncounts / rowsum)),
        "corr_UMI_count_rowsum": float(np.corrcoef(umi, rowsum)[0, 1]),
        "UMI_count_minus_rowsum_max_abs": float(np.abs(umi - rowsum).max()),
        "rule": "v1.2.1 section 5: library size = ncounts if present"}
    q["seconds"] = round(time.time() - t0, 1)
    od = RESULTS / "v1.2" / cl
    create_new(od / "qc_v121.json", json.dumps(q, indent=1, default=float))
    print(json.dumps({cl: q}, default=float))


def wessels():
    with h5py.File(DL / "WesselsSatija2023.h5ad", "r") as h:
        X = h["X"]
        shape = tuple(int(x) for x in X.attrs["shape"])
        d = X["data"]
        neg = nonint = 0
        for s in range(0, d.shape[0], 20_000_000):
            b = d[s:s + 20_000_000]
            neg += int((b < 0).sum())
            if b.dtype.kind == "f":
                nonint += int((np.abs(b - np.round(b)) > 1e-6).sum())
        o = h["obs"]
        pert = col(o, "perturbation")
        npert = col(o, "nperts").astype(int)
        cls = col(o, "Guide.Class")
        dtype = str(d.dtype)
    s = pd.Series(pert)
    cnt = s.value_counts()
    singles = {k: int(v) for k, v in cnt.items() if "_" not in k and k not in ("control",) and npert[s == k][0] == 1}
    doubles = {k: int(v) for k, v in cnt.items() if "_" in k}
    n_ctrl = int((npert == 0).sum())
    d30 = [k for k, v in doubles.items() if v >= 30]
    both30 = [k for k in d30 if all(singles.get(gn, 0) >= 30 for gn in k.split("_"))]
    genes = sorted({gn for k in both30 for gn in k.split("_")})
    gates = {
        "identity": {"n_cells": shape[0], "paper_cells": 31308, "ratio": shape[0] / 31308,
                     "pass": 0.8 <= shape[0] / 31308 <= 1.2 and len(doubles) > 0,
                     "guide_classes": pd.Series(cls).value_counts().to_dict()},
        "G1_evaluable_pairs": {"n": len(both30), "threshold": 60, "pass": len(both30) >= 60},
        "G2_single_coverage": {"pairs_ge30": len(d30), "with_both_singles_ge30": len(both30),
                               "frac": len(both30) / max(len(d30), 1), "pass": len(both30) >= 0.9 * len(d30)},
        "G3_control_cells": {"label": "nperts == 0 (Guide.Class NT)", "n": n_ctrl, "threshold": 500, "pass": n_ctrl >= 500},
        "G4_integer_counts": {"dtype": dtype, "negative": neg, "non_integer": nonint, "pass": neg == 0 and nonint == 0},
        "G5_genes": {"n": len(genes), "threshold": 10, "pass": len(genes) >= 10},
        "n_singles": len(singles), "singles_ge30": int(sum(v >= 30 for v in singles.values())),
        "n_doubles": len(doubles), "doubles_ge30": len(d30),
        "doubles_cells_median": float(np.median(list(doubles.values()))) if doubles else None,
        "max_pairs_per_gene": int(pd.Series([gn for k in both30 for gn in k.split("_")]).value_counts().max()) if both30 else 0}
    gates["all_pass"] = all(gates[k]["pass"] for k in ("identity", "G1_evaluable_pairs", "G2_single_coverage",
                                                      "G3_control_cells", "G4_integer_counts", "G5_genes"))
    create_new(RESULTS / "v1.4" / "wessels" / "gates_v14.json", json.dumps(gates, indent=1, default=float))
    print(json.dumps(gates, default=float))


if __name__ == "__main__":
    if sys.argv[1] == "nadig":
        for c in ("hepg2", "jurkat"):
            nadig(c)
    else:
        wessels()
