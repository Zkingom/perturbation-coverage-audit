# pseudobulk_v04.py -- protocol v0.4 expression representation (section (b)).
# Per cell: x = log1p(1e4 * counts / L), L = obs 'UMI_count' if present (whole-transcriptome
# library size), else the sum over the 8,563 retained genes. Per perturbation row
# (gene_transcript): mean of x over ALL its cells. Controls: all non-targeting cells.
# Each perturbation's value depends only on its own cells, so no cross-split information flows.
# Chunked over cells; memory ~ chunk_rows x n_genes x 8 bytes.
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PROCESSED, raw  # noqa: E402

OUT = PROCESSED / "v0.4"
CHUNK_ROWS = 4000


def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def read_obs_column(f, name):
    """Handle both anndata encodings: legacy (obs/__categories/<name>) and >=0.8 (categorical group)."""
    node = f["obs"][name]
    if isinstance(node, h5py.Group):                                  # anndata >= 0.8 categorical
        return _dec(node["categories"][:])[node["codes"][:]]
    vals = node[:]
    if "__categories" in f["obs"] and name in f["obs/__categories"]:    # legacy categorical
        return _dec(f["obs/__categories"][name][:])[vals]
    return _dec(vals) if vals.dtype.kind in "SO" else vals


def var_gene_names(f):
    node = f["var"]["gene_name"]
    if isinstance(node, h5py.Group):
        return _dec(node["categories"][:])[node["codes"][:]]
    if "__categories" in f["var"] and "gene_name" in f["var/__categories"]:
        return _dec(f["var/__categories/gene_name"][:])[node[:]]
    return _dec(node[:])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    path = raw("K562_essential_raw_singlecell_01.h5ad")
    t0 = time.time()
    with h5py.File(path, "r") as f:
        obs_keys = list(f["obs"].keys())
        gt = read_obs_column(f, "gene_transcript")
        lib_source = "obs.UMI_count" if "UMI_count" in obs_keys else "sum over retained genes"
        umi = read_obs_column(f, "UMI_count").astype(np.float64) if "UMI_count" in obs_keys else None
        genes = var_gene_names(f)
        X = f["X"]
        assert isinstance(X, h5py.Dataset), "expected dense X"
        n, G = X.shape
        rows, inv = np.unique(gt, return_inverse=True)
        sums = np.zeros((len(rows), G))
        ncell = np.bincount(inv, minlength=len(rows))
        sum_sum_genes = 0.0
        for s in range(0, n, CHUNK_ROWS):
            e = min(n, s + CHUNK_ROWS)
            c = X[s:e].astype(np.float64)
            tot = c.sum(1)
            L = umi[s:e] if umi is not None else tot
            sum_sum_genes += tot.sum()
            x = np.log1p(c * (1e4 / L)[:, None])
            ind = sp.csr_matrix((np.ones(e - s), (inv[s:e], np.arange(e - s))), shape=(len(rows), e - s))
            sums += ind @ x
            if (s // CHUNK_ROWS) % 10 == 0:
                print(f"{e}/{n} cells {time.time() - t0:.0f}s", flush=True)
    pb = sums / ncell[:, None]
    np.save(OUT / "pseudobulk_logmean.npy", pb.astype(np.float32))
    pd.DataFrame({"gene_transcript": rows, "n_cells": ncell}).to_csv(OUT / "pseudobulk_rows.csv", index=False)
    pd.Series(genes, name="gene_name").to_csv(OUT / "pseudobulk_genes.csv", index=False)
    info = {"cells": int(n), "genes": int(G), "rows": int(len(rows)), "library_size": lib_source,
            "mean_fraction_of_UMI_in_retained_genes": (float(sum_sum_genes / umi.sum()) if umi is not None else None),
            "obs_keys": obs_keys, "seconds": round(time.time() - t0)}
    (OUT / "pseudobulk_v0.4_build.json").write_text(json.dumps(info, indent=1), encoding="utf8")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
