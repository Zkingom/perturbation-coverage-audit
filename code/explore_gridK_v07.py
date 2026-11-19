# explore_gridK_v07.py -- POST-CONFIRMATION EXPLORATORY run on RPE1 with extended K / d grids
# (protocol reports/protocol/explore_gridK_rpe1_v1.0.md). Frozen code is not edited: candidate grids, the GO-embedding
# matrix (for d = 200) and the output directory are overridden at runtime, then the validated resumable pipeline
# (confirm_run_v07.do_run, low-memory loader, F2c finish) is called.
#   python explore_gridK_v07.py prep                       GO ASE with 200 columns (prior only) + nesting check
#   python explore_gridK_v07.py check  [--out DIR]         ORIGINAL grids, split 0, own output dir (byte check)
#   python explore_gridK_v07.py run    [--only ...]        EXTENDED grids, all splits + F2c finish
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PROCESSED, RESULTS  # noqa: E402

K_EXT = (5, 10, 20, 30, 50, 100, 200)
D_EXT = (2, 5, 10, 20, 30, 50, 100, 200)
OUT = RESULTS / "explore_v0.7" / "rpe1_gridK"
ASE200 = PROCESSED / "priors_explore" / "go_ase200_hybrid.npz"


def prep():
    """ASE with 200 columns of the same GO graph (same construction as build_priors_v04.go_ase). Columns 1-100 are
    kept EXACTLY as in go_ase100.npz (so every d <= 100 is unchanged); columns 101-200 come from the k = 200
    solution. Check: the first 100 columns of the k = 200 solution equal go_ase100 up to sign."""
    import tarfile
    import tempfile
    import pandas as pd
    import scipy.sparse as sp
    from scipy.sparse.linalg import eigsh
    from paths import DL, raw
    with tempfile.TemporaryDirectory(dir=str(DL)) as td:
        with tarfile.open(raw("go_essential_all.tar.gz")) as tf:
            tf.extractall(td)
        goe = pd.read_csv(Path(td) / "go_essential_all" / "go_essential_all.csv")
    gl = np.array(sorted(set(goe.source) | set(goe.target)))
    gi = {g: i for i, g in enumerate(gl)}
    A = sp.csr_matrix((goe.importance.values.astype(np.float64),
                       (goe.source.map(gi).values, goe.target.map(gi).values)), shape=(len(gl), len(gl)))
    A = A.maximum(A.T)
    off = A.copy()
    off.setdiag(0)
    off.eliminate_zeros()
    strength = np.asarray(off.sum(1)).ravel()
    Aug = off + sp.diags(strength / (len(gl) - 1))
    vals, vecs = eigsh(Aug, k=200, which="LM")
    o = np.argsort(-np.abs(vals))
    X = vecs[:, o] * np.sqrt(np.abs(vals[o]))
    z = np.load(PROCESSED / "priors" / "go_ase100.npz", allow_pickle=True)
    assert np.array_equal(z["genes"], gl)
    X100 = z["ase"]
    sign = np.sign((X[:, :100] * X100).sum(0))
    dif = np.abs(X[:, :100] * sign - X100).max(0)
    ev_dif = float(np.abs(np.abs(vals[o][:100]) - np.abs(z["eigenvalues"])).max())
    H = np.hstack([X100, X[:, 100:]])
    ASE200.parent.mkdir(parents=True, exist_ok=True)
    np.savez(ASE200, genes=gl, ase=H, eigenvalues=vals[o],
             note="columns 1-100 = go_ase100.npz exactly; columns 101-200 from k=200 eigsh of the same matrix")
    chk = {"max_abs_diff_first100_cols": float(dif.max()), "n_cols_diff_gt_1e-6": int((dif > 1e-6).sum()),
           "max_abs_eigenvalue_diff_first100": ev_dif, "shape": list(H.shape)}
    (ASE200.parent / "go_ase200_hybrid_check.json").write_text(json.dumps(chk, indent=1))
    print(json.dumps(chk), flush=True)
    return 0


def _setup(out_dir, extended):
    import confirm_run_v07 as cr
    import pipeline_v05 as pl5
    pdir, sdir, _ = cr.dirs("rpe1")
    cr.dirs = lambda ds: (pdir, sdir, Path(out_dir))           # results go to the exploratory directory only
    if extended:
        pl5.K_CANDIDATES, pl5.K_MAX, pl5.D_CANDIDATES = K_EXT, max(K_EXT), D_EXT
        z = np.load(ASE200, allow_pickle=True)
        orig = cr.load_lowmem

        def load_ext(ds):
            D = orig(ds)
            g0 = D["pri"]["ase_ix"]
            assert {g: i for i, g in enumerate(z["genes"])} == g0, "GO gene order differs"
            assert np.array_equal(D["pri"]["ase"], z["ase"][:, :100]), "first 100 columns differ from go_ase100"
            D["pri"]["ase"] = z["ase"]
            return D
        cr.load_lowmem = load_ext
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    (Path(out_dir) / "explore_gridK_setup.json").write_text(json.dumps(
        {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "extended": extended,
         "K_CANDIDATES": list(pl5.K_CANDIDATES), "D_CANDIDATES": list(pl5.D_CANDIDATES), "out": str(out_dir)}, indent=1))
    return cr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["prep", "check", "run"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--only", type=int, nargs="*", default=None)
    a = ap.parse_args()
    if a.mode == "prep":
        return prep()
    if a.mode == "check":
        cr = _setup(a.out or (RESULTS / "explore_v0.7" / "_gridK_original_check"), extended=False)
        return cr.do_run("rpe1", [0])
    cr = _setup(a.out or OUT, extended=True)
    return cr.do_run("rpe1", a.only, finish_f2c=a.only is None)


if __name__ == "__main__":
    raise SystemExit(main())
