# auditB_clusters_v05.py -- frozen audit-B clusters (protocol v0.5). Prior-only: no expression is read.
# Graph: E units; each gene keeps its 20 largest GO Jaccard neighbours (original GEARS gene2go, no term
# filter), J > 0 only; symmetrised by max; edge weight = Jaccard. Leiden RBConfigurationVertexPartition,
# resolution 8, seed = master_seed. Output: data/processed/v0.5/auditB_clusters_original_r8.csv
import json
import pickle
import sys
from pathlib import Path

import igraph as ig
import leidenalg as la
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED, raw  # noqa: E402
from models_v04 import go_incidence, jaccard_all  # noqa: E402

KNN, RESOLUTION = 20, 8.0


def build():
    U = pd.read_csv(PROCESSED / "v0.3_figshare_k562ess" / "p0v3_units_all_targets.csv")
    U = U[U.E | U.always_train].reset_index(drop=True)
    E_idx = np.where(U.E.values)[0]
    with open(raw("gene2go_all.pkl"), "rb") as fh:
        gene2go = pickle.load(fh)
    J = jaccard_all(go_incidence(U.target.values[E_idx], gene2go))
    n = len(J)
    W = np.zeros_like(J)
    for i in range(n):
        nb = np.argsort(-J[i], kind="stable")[:KNN]
        nb = nb[J[i, nb] > 0]
        W[i, nb] = J[i, nb]
    W = np.maximum(W, W.T)
    iu = np.triu_indices(n, 1)
    m = W[iu] > 0
    g = ig.Graph(n=n, edges=list(zip(iu[0][m], iu[1][m])))
    g.es["weight"] = W[iu][m]
    part = la.find_partition(g, la.RBConfigurationVertexPartition, weights="weight",
                             resolution_parameter=RESOLUTION, seed=MASTER_SEED)
    memb = np.array(part.membership)
    out = pd.DataFrame({"target": U.target.values[E_idx], "tss_group": U.tss_group.values[E_idx], "cluster": memb})
    d = PROCESSED / "v0.5"
    d.mkdir(parents=True, exist_ok=True)
    out.to_csv(d / "auditB_clusters_original_r8.csv", index=False)
    sz = np.bincount(memb)
    info = {"n_units": n, "n_edges": g.ecount(), "n_clusters": int(len(sz)), "max_cluster": int(sz.max()),
            "max_frac_E": float(sz.max() / n), "knn": KNN, "resolution": RESOLUTION, "seed": MASTER_SEED,
            "igraph": ig.__version__, "leidenalg": la.version}
    (d / "auditB_clusters_build.json").write_text(json.dumps(info, indent=1), encoding="utf8")
    return out, info


if __name__ == "__main__":
    print(json.dumps(build()[1]))
