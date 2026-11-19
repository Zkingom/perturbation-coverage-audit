# build_priors_v04.py -- prior-only objects for protocol v0.4 (no expression data is read here).
#  * GO adjacency spectral embedding (ASE) with d = 100 on GEARS go_essential_all.csv; candidate
#    embeddings d in {2, 5, 10, 20, 30, 50, 100} are the first d columns (v0.4.2). Checked against the archived d=2/d=10.
#  * STRING v12.0 experiments+databases combined score restricted to all screen target genes.
import json
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import DL, PROCESSED, raw  # noqa: E402

OUT = PROCESSED / "priors"
D_MAX = 100   # v0.4.2: d grid {2,5,10,20,30,50,100}


def go_ase():
    with tempfile.TemporaryDirectory(dir=str(DL)) as td:
        with tarfile.open(raw("go_essential_all.tar.gz")) as tf:
            tf.extractall(td)
        goe = pd.read_csv(Path(td) / "go_essential_all" / "go_essential_all.csv")
    gl = np.array(sorted(set(goe.source) | set(goe.target)))
    gi = {g: i for i, g in enumerate(gl)}
    A = sp.csr_matrix((goe.importance.values.astype(np.float64),
                       (goe.source.map(gi).values, goe.target.map(gi).values)), shape=(len(gl), len(gl)))
    del goe
    A = A.maximum(A.T)                       # igraph mode="undirected" == "max"
    off = A.copy()
    off.setdiag(0)
    off.eliminate_zeros()
    strength = np.asarray(off.sum(1)).ravel()
    Aug = off + sp.diags(strength / (len(gl) - 1))   # igraph default cvec
    vals, vecs = eigsh(Aug, k=D_MAX, which="LM")
    o = np.argsort(-np.abs(vals))
    X = vecs[:, o] * np.sqrt(np.abs(vals[o]))
    ref = np.load(OUT / "go_ase_embeddings.npz", allow_pickle=True)
    assert np.array_equal(ref["genes"], gl)
    diffs = {}
    for d in (2, 10):
        R = ref[f"ase{d}"]
        sign = np.sign((X[:, :d] * R).sum(0))
        diffs[d] = float(np.abs(X[:, :d] * sign - R).max())
        X[:, :d] = X[:, :d] * sign            # align signs with archive (sign is arbitrary)
    np.savez(OUT / f"go_ase{D_MAX}.npz", genes=gl, ase=X, eigenvalues=vals[o],
             note=f"ASE of GEARS go_essential_all.csv, d={D_MAX}, columns ordered by |eigenvalue|")
    return diffs, vals[o]


def string_universe(targets):
    info = pd.read_csv(raw("9606.protein.info.v12.0.txt.gz"), sep="\t")
    info.columns = [c.lstrip("#") for c in info.columns]
    p2g = dict(zip(info.string_protein_id, info.preferred_name))
    univ = set(targets)
    keep = []
    for ch in pd.read_csv(raw("9606.protein.links.detailed.v12.0.txt.gz"), sep=" ", chunksize=2_000_000,
                          usecols=["protein1", "protein2", "experimental", "database"]):
        g1 = ch.protein1.map(p2g)
        g2 = ch.protein2.map(p2g)
        m = g1.isin(univ) & g2.isin(univ) & ((ch.experimental > 0) | (ch.database > 0))
        keep.append(pd.DataFrame({"g1": g1[m], "g2": g2[m], "exp": ch.experimental[m] / 1000, "db": ch.database[m] / 1000}))
    st = pd.concat(keep, ignore_index=True)
    prior = 0.041

    def nopr(s):
        return np.clip((s - prior) / (1 - prior), 0, None)
    comb = 1 - (1 - nopr(st.exp)) * (1 - nopr(st.db))
    st["score"] = comb + prior * (1 - comb)
    st = st.groupby(["g1", "g2"], as_index=False).score.max()
    st.to_parquet(OUT / "string_expdb_v0.4_universe.parquet")
    return len(st)


if __name__ == "__main__":
    units = pd.read_csv(PROCESSED / "v0.3_figshare_k562ess" / "p0v3_units_all_targets.csv")
    diffs, ev = go_ase()
    n_edges = int(len(pd.read_parquet(OUT / "string_expdb_v0.4_universe.parquet")))  # unchanged; not rebuilt
    rep = {"ase_check_max_abs_diff_vs_archive": diffs, "ase_top_eigenvalues": [round(float(x), 3) for x in ev],
           "string_edges_in_universe": n_edges}
    (OUT / f"priors_v0.4.2_ase{D_MAX}_build.json").write_text(json.dumps(rep, indent=1), encoding="utf8")
    print(json.dumps(rep))
