# smoke_test.py -- post-migration checks (protocol v0.4). Reads only from ROOT / DL via paths.py.
# 1) P0 v0.3 counts and unit table identical to the archived table
# 2) fold sizes for seeds master_seed+0..2 identical to the recorded values
# 3) GO adjacency spectral embedding identical (up to column sign) to the archived npz
import json
import sys
import time
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import DL, MASTER_SEED, PROCESSED, RESULTS, raw  # noqa: E402
from audit_lib import random_group_split  # noqa: E402
from p0_build import UNITS_ARCHIVE, build_units  # noqa: E402

EXPECTED = {"E": 1786, "tss_groups_E": 1721, "always_train": 257, "dropped": 14,
            "secondary_rows": 119, "readout_genes": 8563}
EXPECTED_FOLDS = [{"train": 1074, "val": 178, "cal": 267, "test": 267},
                  {"train": 1075, "val": 177, "cal": 267, "test": 267},
                  {"train": 1074, "val": 178, "cal": 267, "test": 267}]


def check_units():
    U, readout, n_sec = build_units()
    got = {"E": int(U.E.sum()), "tss_groups_E": int(U[U.E].tss_group.nunique()),
           "always_train": int(U.always_train.sum()), "dropped": int(U.dropped_tss_conflict.sum()),
           "secondary_rows": n_sec, "readout_genes": len(readout)}
    import filecmp
    import tempfile
    with tempfile.TemporaryDirectory(dir=str(DL)) as td:   # rebuilt table is not kept (v0.4.2)
        out = Path(td) / "p0_units_rebuilt.csv"
        U.to_csv(out, index=False)
        a = pd.read_csv(UNITS_ARCHIVE)
        b = pd.read_csv(out)
        pd.testing.assert_frame_equal(a, b, check_exact=True)
        assert filecmp.cmp(out, UNITS_ARCHIVE, shallow=False), "rebuilt unit CSV not byte-identical"
        out = "byte-identical (temporary file, removed)"
    folds = []
    Eu = U[U.E].reset_index(drop=True)
    for i in range(3):
        f = random_group_split(Eu.tss_group.values, MASTER_SEED + i)
        folds.append({k: len(v) for k, v in f.items()})
    return got, folds, out


def check_ase(tmpdir):
    with tarfile.open(raw("go_essential_all.tar.gz")) as tf:
        tf.extractall(tmpdir)
    goe = pd.read_csv(Path(tmpdir) / "go_essential_all" / "go_essential_all.csv")
    gnodes = set(goe.source) | set(goe.target)
    gl = np.array(sorted(gnodes))
    gi = {g: i for i, g in enumerate(gl)}
    A = sp.csr_matrix((goe.importance.values.astype(np.float64),
                       (goe.source.map(gi).values, goe.target.map(gi).values)), shape=(len(gl), len(gl)))
    del goe
    A = A.maximum(A.T)
    off = A.copy()
    off.setdiag(0)
    off.eliminate_zeros()
    strength = np.asarray(off.sum(1)).ravel()
    Aug = off + sp.diags(strength / (len(gl) - 1))
    ref = np.load(PROCESSED / "priors" / "go_ase_embeddings.npz", allow_pickle=True)
    assert np.array_equal(ref["genes"], gl), "gene order differs"
    res = {}
    for d in (2, 10):
        vals, vecs = eigsh(Aug, k=d, which="LM")
        o = np.argsort(-np.abs(vals))
        X = vecs[:, o] * np.sqrt(np.abs(vals[o]))
        R = ref[f"ase{d}"]
        sign = np.sign((X * R).sum(0))
        res[d] = float(np.abs(X * sign - R).max())
    return res


if __name__ == "__main__":
    import tempfile
    got, folds, out = check_units()
    with tempfile.TemporaryDirectory(dir=str(DL)) as td:
        ase = check_ase(td)
    report = {"counts": got, "counts_expected": EXPECTED, "counts_ok": got == EXPECTED,
              "unit_table_identical": True, "folds": folds, "folds_ok": folds == EXPECTED_FOLDS,
              "ase_max_abs_diff": ase, "ase_ok": all(v < 1e-8 for v in ase.values())}
    (RESULTS / f"smoke_test_{time.strftime('%Y-%m-%d')}.json").write_text(json.dumps(report, indent=1), encoding="utf8")
    print(json.dumps(report))
    assert report["counts_ok"] and report["folds_ok"] and report["ase_ok"]
