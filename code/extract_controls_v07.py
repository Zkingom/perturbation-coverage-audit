# extract_controls_v07.py -- ILLUSTRATIVE noise-floor analysis (outside every protocol; user request 2026-10-16),
# step 1: per-cell expression of non-targeting control cells, transformed exactly as in pseudobulk_v06
# (x = log1p(1e4 * counts / UMI_count), computed in float64), stored as float32 rows (cells x readout genes).
# At most --max-cells control cells are kept (uniform random subset, fixed seed) to bound memory.
#   python extract_controls_v07.py --dataset {dev,rpe1,gwps} [--input SC.h5ad] [--max-cells 20000]
# Output: data/processed/noise_floor_v07/<dataset>/ctrl_x.npy, ctrl_cells.json (atomic writes).
import argparse
import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED  # noqa: E402
from pseudobulk_v04 import read_obs_column, var_gene_names  # noqa: E402

PB = {"dev": PROCESSED / "v0.6" / "k562ess", "rpe1": PROCESSED / "v0.7" / "rpe1", "gwps": PROCESSED / "v0.7" / "gwps"}
SEED = MASTER_SEED + 9500


def default_input(ds):
    import confirm_run_v07 as cr
    return {"dev": lambda: cr.resolve("K562_essential_raw_singlecell_01.h5ad"),
            "rpe1": lambda: cr.resolve("rpe1_raw_singlecell_01.h5ad"),
            "gwps": lambda: Path(r"H:\cs+bio_file\K562_gwps_raw_singlecell_01.h5ad")}[ds]()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=list(PB))
    ap.add_argument("--input", default=None)
    ap.add_argument("--max-cells", type=int, default=20000)
    ap.add_argument("--chunk-rows", type=int, default=4000)
    a = ap.parse_args()
    inp = Path(a.input) if a.input else default_input(a.dataset)
    od = PROCESSED / "noise_floor_v07" / a.dataset
    od.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    pb_genes = pd.read_csv(PB[a.dataset] / "pseudobulk_genes.csv").gene_name.values
    with h5py.File(inp, "r") as f:
        X = f["X"]
        n, G = X.shape
        gt = read_obs_column(f, "gene_transcript")
        umi = read_obs_column(f, "UMI_count").astype(np.float64)
        genes = np.asarray(var_gene_names(f))
        assert len(genes) == len(pb_genes) and (genes == pb_genes).all(), "gene order differs from pseudobulk"
        ctrl = np.flatnonzero(pd.Series(gt).str.contains("non-targeting", case=False).values)
        keep = np.sort(np.random.default_rng(SEED).choice(ctrl, min(a.max_cells, len(ctrl)), replace=False))
        print(f"{inp}: {n} cells, {len(ctrl)} non-targeting, keeping {len(keep)}", flush=True)
        tmp = od / "ctrl_x.npy.part"
        mm = np.lib.format.open_memmap(tmp, mode="w+", dtype=np.float32, shape=(len(keep), G))
        pos = 0
        for s in range(0, n, a.chunk_rows):
            e = min(n, s + a.chunk_rows)
            sel = keep[(keep >= s) & (keep < e)]
            if len(sel) == 0:
                continue
            c = X[s:e][sel - s].astype(np.float64)
            mm[pos:pos + len(sel)] = np.log1p(c * (1e4 / umi[sel])[:, None]).astype(np.float32)
            pos += len(sel)
            if (s // a.chunk_rows) % 20 == 0:
                print(f"{time.strftime('%H:%M:%S')} rows {e}/{n} kept {pos}", flush=True)
        assert pos == len(keep)
        mm.flush()
        mm._mmap.close()
        del mm
    os.replace(tmp, od / "ctrl_x.npy")
    info = {"dataset": a.dataset, "input": str(inp), "n_cells_file": int(n), "n_ctrl_total": int(len(ctrl)),
            "n_ctrl_kept": int(len(keep)), "seed": SEED, "genes": int(G), "per_cell_transform": "log1p(1e4*counts/UMI_count), float64 -> stored float32",
            "control_guides": int(pd.Series(gt[ctrl]).nunique()), "seconds": round(time.time() - t0)}
    (od / "ctrl_cells.json.tmp").write_text(json.dumps(info, indent=1))
    os.replace(od / "ctrl_cells.json.tmp", od / "ctrl_cells.json")
    print(json.dumps(info), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
