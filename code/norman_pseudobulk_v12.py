# norman_pseudobulk_v12.py -- protocol combo_crosscell_v1.2, section 2.1 (C0) for Norman 2019 (scPerturb file).
# Backed, chunked reading of the CSC count matrix (cells x genes, column-compressed): genes are processed in
# column blocks, so memory is bounded by one block. Two passes:
#   pass 1: per-cell total counts (library size) and per-gene mean raw count in control cells;
#   pass 2: per-cell transform x = log1p(1e4 * counts / total_counts) (float64), as in pseudobulk_v06, averaged per
#           unit (full, half A, half B).
# Library size: the scPerturb obs column UMI_count holds guide UMIs, not transcript UMIs, so the total transcript
# count of each cell is recomputed from X (and compared with obs 'ncounts'). This matches the definition used for
# the Replogle files, where UMI_count is the total transcript UMI count.
# Units: unordered gene sets ('A' or 'A_B' with genes sorted), 'control'. Half split: per unit (units in sorted
# order), cells in file order are permuted with rng = default_rng(MASTER_SEED + 15000 + 0); A = first floor(n/2).
# Outputs (atomic): data/processed/v1.2/norman/{pb_full,pb_A,pb_B}.npy (units x all genes, float64),
# units.csv, genes.csv, ctrl_mean_counts.npy, c0_build.json.
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED  # noqa: E402

SRC = Path(r"G:\biology\geo_GSE264667\NormanWeissman2019_filtered.h5ad")
OUT = PROCESSED / "v1.2" / "norman"
BLOCK = 1500


def _cat(o, k):
    x = o[k]
    c = np.array([v.decode() if isinstance(v, bytes) else v for v in x["categories"][:]], dtype=object)
    return c[x["codes"][:]]


def unit_label(s):
    return "control" if s == "control" else "_".join(sorted(s.split("_")))


def blocks(n):
    for s in range(0, n, BLOCK):
        yield s, min(n, s + BLOCK)


def read_block(X, ip, s, e, n_cells):
    a, b = ip[s], ip[e]
    data = X["data"][a:b].astype(np.float64)
    ind = X["indices"][a:b]
    return sps.csc_matrix((data, ind, ip[s:e + 1] - a), shape=(n_cells, e - s))


def save_npy(arr, path):
    tmp = path.with_suffix(".npy.part")
    with open(tmp, "wb") as fh:
        np.save(fh, arr)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    h = h5py.File(SRC, "r")
    o, X = h["obs"], h["X"]
    n_cells, n_genes = (int(v) for v in X.attrs["shape"])
    ip = X["indptr"][:].astype(np.int64)
    genes = np.array([g.decode() for g in h["var/_index"][:]], dtype=object)
    lab = pd.Series(_cat(o, "perturbation")).map(unit_label).values
    units = np.array(sorted(set(lab)), dtype=object)
    uidx = {u: i for i, u in enumerate(units)}
    cu = np.array([uidx[v] for v in lab])
    ctrl = lab == "control"
    # pass 1
    lib = np.zeros(n_cells)
    ctrl_sum = np.zeros(n_genes)
    nonint = 0
    for s, e in blocks(n_genes):
        B = read_block(X, ip, s, e, n_cells)
        nonint += int((B.data != np.round(B.data)).sum()) + int((B.data < 0).sum())
        lib += np.asarray(B.sum(1)).ravel()
        ctrl_sum[s:e] = np.asarray(B[ctrl].sum(0)).ravel()
    assert nonint == 0, f"{nonint} non-integer or negative values: X is not raw counts"
    ncounts = o["ncounts"][:]
    assert (lib > 0).all()
    # half split
    rng = np.random.default_rng(MASTER_SEED + 15000 + 0)
    half = np.full(n_cells, -1)
    for i, u in enumerate(units):
        idx = np.where(cu == i)[0]
        p = rng.permutation(idx)
        half[p[: len(idx) // 2]] = 0
        half[p[len(idx) // 2:]] = 1
    assert (half >= 0).all()
    ncell = np.bincount(cu, minlength=len(units))
    nA = np.bincount(cu[half == 0], minlength=len(units))
    nB = ncell - nA
    ind = {k: sps.csr_matrix((np.ones(m.sum()), (cu[m], np.where(m)[0])), shape=(len(units), n_cells))
           for k, m in (("full", np.ones(n_cells, bool)), ("A", half == 0), ("B", half == 1))}
    den = {"full": ncell, "A": nA, "B": nB}
    pb = {k: np.zeros((len(units), n_genes)) for k in ind}
    scale = 1e4 / lib
    # pass 2
    for s, e in blocks(n_genes):
        B = read_block(X, ip, s, e, n_cells).tocsr()
        B.data = np.log1p(B.data * np.repeat(scale, np.diff(B.indptr)))
        for k in ind:
            pb[k][:, s:e] = (ind[k] @ B).toarray() / den[k][:, None]
    for k in pb:
        save_npy(pb[k], OUT / f"pb_{k}.npy")
    save_npy(ctrl_sum / ctrl.sum(), OUT / "ctrl_mean_counts.npy")
    U = pd.DataFrame({"unit": units, "n_cells": ncell, "n_A": nA, "n_B": nB,
                      "n_genes_in_unit": [0 if u == "control" else len(u.split("_")) for u in units]})
    U.to_csv(OUT / "units.csv", index=False)
    pd.DataFrame({"gene": genes}).to_csv(OUT / "genes.csv", index=False)
    info = {"source": str(SRC), "n_cells": n_cells, "n_genes": n_genes, "n_units": len(units),
            "lib_vs_obs_ncounts_max_abs_diff": float(np.abs(lib - ncounts).max()),
            "lib_vs_obs_ncounts_corr": float(np.corrcoef(lib, ncounts)[0, 1]),
            "nonint_values": nonint, "half_seed": MASTER_SEED + 15000, "block": BLOCK,
            "seconds": round(time.time() - t0, 1)}
    (OUT / "c0_build.json.tmp").write_text(json.dumps(info, indent=1))
    os.replace(OUT / "c0_build.json.tmp", OUT / "c0_build.json")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
