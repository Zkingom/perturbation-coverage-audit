# pseudobulk_v06.py -- protocol v0.6 expression representation for any Replogle single-cell h5ad
# (dense X: cells x genes, raw UMI counts). Same per-cell transform as pseudobulk_v04.py:
#   x = log1p(1e4 * counts / L), L = obs 'UMI_count'; per perturbation row (gene_transcript): mean over cells.
# Additionally (v0.6 section 2): cells of every row are split into halves A / B (seed master_seed + 3000,
# fixed once, independent of the gene pass), giving pseudobulk_A / pseudobulk_B; full = (sum_A + sum_B) / n.
# Quality check against the publisher's bulk file (v0.6 section 1.0): per-row mean raw UMI from the single
# cells vs bulk X; only the maximum absolute / relative difference and a pass flag are written. The
# per-row QC means themselves are discarded (never saved).
# Memory: --passes P splits the genes into P groups and reads the file P times (low-memory mode for GWPS).
# Before starting, and after every chunk, available physical memory is checked: < --min-free-gb at
# start -> stop with a message; < --pause-free-gb during the run -> pause until >= --min-free-gb again.
#
# Usage: python pseudobulk_v06.py --input X.h5ad --bulk X_bulk.h5ad --outdir DIR [--passes 2]
#        [--chunk-rows 4000] [--status-file S.json] [--min-free-gb 2.5] [--pause-free-gb 1.0]
import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED  # noqa: E402
from pseudobulk_v04 import read_obs_column, var_gene_names  # noqa: E402

HALF_SEED = MASTER_SEED + 3000
QC_ATOL, QC_RTOL = 1e-3, 1e-4


def avail_gb():
    import psutil
    return psutil.virtual_memory().available / 1024 ** 3


def utc():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def write_status(path, d):
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(d, default=str), encoding="utf8")
    os.replace(tmp, path)


def half_labels(inv, n_rows, seed=HALF_SEED):
    """For every row, a random half of its cells (floor(n/2)) goes to A, the rest to B. Deterministic."""
    rng = np.random.default_rng(seed)
    half = np.zeros(len(inv), np.int8)          # 0 = A, 1 = B
    order = np.argsort(inv, kind="stable")
    bounds = np.r_[0, np.cumsum(np.bincount(inv, minlength=n_rows))]
    for r in range(n_rows):
        idx = order[bounds[r]:bounds[r + 1]]
        perm = rng.permutation(len(idx))
        half[idx[perm[len(idx) // 2:]]] = 1
    return half


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def build(inp, outdir, bulk=None, passes=1, chunk_rows=4000, status=None, min_free=2.5, pause_free=1.0,
          log=print):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    st = {"input": str(inp), "outdir": str(outdir), "passes": passes, "started_utc": utc(), "pid": os.getpid()}
    free = avail_gb()
    if free < min_free:
        msg = (f"available memory {free:.2f} GB < {min_free} GB: close other programs and rerun "
               f"(nothing was read)")
        write_status(status, {**st, "phase": "stopped_low_memory", "available_gb": free, "updated_utc": utc()})
        log(msg)
        return {"status": "stopped_low_memory", "available_gb": free}
    t0 = time.time()
    with h5py.File(inp, "r") as f:
        X = f["X"]
        assert isinstance(X, h5py.Dataset), "expected dense X"
        n, G = X.shape
        obs_keys = list(f["obs"].keys())
        assert "UMI_count" in obs_keys, "obs.UMI_count required (library size)"
        gt = read_obs_column(f, "gene_transcript")
        umi = read_obs_column(f, "UMI_count").astype(np.float64)
        genes = var_gene_names(f)
        rows, inv = np.unique(gt, return_inverse=True)
        R = len(rows)
        half = half_labels(inv, R)
        ncell = np.bincount(inv, minlength=R)
        nA = np.bincount(inv[half == 0], minlength=R)
        nB = ncell - nA
        st.update({"cells": int(n), "genes": int(G), "rows": int(R)})
        log(f"{inp}: {n} cells x {G} genes, {R} rows, passes={passes}, free={free:.2f} GB")
        # outputs are written column-block by column-block into .npy memmaps (never fully in RAM)
        mm = {nm: np.lib.format.open_memmap(outdir / f"{nm}.npy.part", mode="w+", dtype=np.float32, shape=(R, G))
              for nm in ("pseudobulk_logmean", "pseudobulk_logmean_A", "pseudobulk_logmean_B")}
        qc_acc = None
        if bulk is not None:
            fb = h5py.File(bulk, "r")
            bgt = read_obs_column(fb, "gene_transcript")
            bgenes = var_gene_names(fb)
            bnc = read_obs_column(fb, "num_cells_filtered") if "num_cells_filtered" in fb["obs"] else None
            bix = {g: i for i, g in enumerate(bgt)}
            common = np.array([i for i, r in enumerate(rows) if r in bix], dtype=int)
            bi = np.array([bix[rows[i]] for i in common], dtype=int)
            same_genes = len(bgenes) == len(genes) and bool((np.asarray(bgenes) == np.asarray(genes)).all())
            qc_acc = {"max_abs": 0.0, "max_rel": 0.0, "n_out": 0}
        groups = np.array_split(np.arange(G), passes)
        done_chunks, total_chunks = 0, passes * int(np.ceil(n / chunk_rows))
        for gi, gcols in enumerate(groups):
            g0, g1 = int(gcols[0]), int(gcols[-1]) + 1
            sA = np.zeros((R, g1 - g0))
            sB = np.zeros((R, g1 - g0))
            sR = np.zeros((R, g1 - g0)) if bulk is not None else None
            for s in range(0, n, chunk_rows):
                e = min(n, s + chunk_rows)
                while avail_gb() < pause_free:
                    write_status(status, {**st, "phase": "paused_low_memory", "available_gb": avail_gb(),
                                          "updated_utc": utc()})
                    log(f"paused: available memory below {pause_free} GB; close other programs")
                    time.sleep(30)
                    if avail_gb() >= min_free:
                        break
                c = X[s:e, g0:g1].astype(np.float64)
                x = np.log1p(c * (1e4 / umi[s:e])[:, None])
                hA = half[s:e] == 0
                ind = lambda m: sp.csr_matrix((np.ones(int(m.sum())), (inv[s:e][m], np.where(m)[0])),   # noqa: E731
                                              shape=(R, e - s))
                sA += ind(hA) @ x
                sB += ind(~hA) @ x
                if sR is not None:
                    sR += ind(np.ones(e - s, bool)) @ c
                del c, x
                done_chunks += 1
                if done_chunks % 25 == 0 or done_chunks == total_chunks:
                    el = time.time() - t0
                    write_status(status, {**st, "phase": f"pass {gi + 1}/{passes}", "chunks": done_chunks,
                                          "total_chunks": total_chunks, "elapsed_s": round(el),
                                          "eta_s": round(el / done_chunks * (total_chunks - done_chunks)),
                                          "available_gb": round(avail_gb(), 2), "updated_utc": utc()})
            with np.errstate(invalid="ignore", divide="ignore"):
                mm["pseudobulk_logmean"][:, g0:g1] = ((sA + sB) / ncell[:, None]).astype(np.float32)
                mm["pseudobulk_logmean_A"][:, g0:g1] = (sA / nA[:, None]).astype(np.float32)
                mm["pseudobulk_logmean_B"][:, g0:g1] = (sB / nB[:, None]).astype(np.float32)
            del sA, sB
            if sR is not None and same_genes and len(common):
                rm = sR[common] / ncell[common][:, None]                 # per-row QC means: transient only
                bx = fb["X"][:, g0:g1][bi].astype(np.float64)
                dif = np.abs(rm - bx)
                qc_acc["max_abs"] = max(qc_acc["max_abs"], float(dif.max()))
                qc_acc["max_rel"] = max(qc_acc["max_rel"], float((dif / np.maximum(np.abs(bx), 1e-12)).max()))
                qc_acc["n_out"] += int((dif > QC_ATOL + QC_RTOL * np.abs(bx)).sum())
                del rm, bx, dif
            del sR
        for nm_ in list(mm):
            mm[nm_].flush()
            mm[nm_]._mmap.close()               # release the Windows file handle before renaming
        mm.clear()
        for nm in ("pseudobulk_logmean", "pseudobulk_logmean_A", "pseudobulk_logmean_B"):
            os.replace(outdir / f"{nm}.npy.part", outdir / f"{nm}.npy")
    qc = None
    if bulk is not None:
        if same_genes and len(common):
            ncell_ok = None if bnc is None else bool((np.asarray(bnc)[bi].astype(int) == ncell[common]).all())
            qc = {"rows_single_cell": int(R), "rows_bulk": int(len(bgt)), "rows_compared": int(len(common)),
                  "genes_identical": True, "max_abs_diff": qc_acc["max_abs"], "max_rel_diff": qc_acc["max_rel"],
                  "n_entries_outside_tol": qc_acc["n_out"], "tolerance": f"{QC_ATOL} + {QC_RTOL} x |bulk|",
                  "cell_counts_match_bulk": ncell_ok, "pass": bool(qc_acc["n_out"] == 0 and ncell_ok is not False)}
        else:
            qc = {"genes_identical": same_genes, "rows_compared": int(len(common)), "pass": False}
        fb.close()
    pd.DataFrame({"gene_transcript": rows, "n_cells": ncell, "n_A": nA, "n_B": nB}).to_csv(outdir / "pseudobulk_rows.csv",
                                                                                          index=False)
    pd.Series(genes, name="gene_name").to_csv(outdir / "pseudobulk_genes.csv", index=False)
    info = {"input": str(inp), "cells": int(n), "genes": int(G), "rows": int(R), "passes": passes,
            "chunk_rows": chunk_rows, "half_seed": HALF_SEED, "library_size": "obs.UMI_count",
            "qc_vs_bulk": qc, "seconds": round(time.time() - t0), "obs_keys": obs_keys,
            "outputs_sha256": {p.name: sha256(p) for p in sorted(outdir.glob("pseudobulk_*"))}}
    (outdir / "pseudobulk_build.json").write_text(json.dumps(info, indent=1), encoding="utf8")
    write_status(status, {**st, "phase": "done", "qc_pass": None if qc is None else qc["pass"], "updated_utc": utc()})
    log(json.dumps({k: info[k] for k in ("cells", "genes", "rows", "passes", "seconds", "qc_vs_bulk")}))
    return {"status": "done", **info}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bulk", default=None)
    ap.add_argument("--passes", type=int, default=1)
    ap.add_argument("--chunk-rows", type=int, default=4000)
    ap.add_argument("--status-file", default=None)
    ap.add_argument("--min-free-gb", type=float, default=2.5)
    ap.add_argument("--pause-free-gb", type=float, default=1.0)
    a = ap.parse_args()
    r = build(a.input, a.outdir, a.bulk, a.passes, a.chunk_rows, Path(a.status_file) if a.status_file else None,
              a.min_free_gb, a.pause_free_gb)
    return 0 if r["status"] == "done" else 2


if __name__ == "__main__":
    raise SystemExit(main())
