# gears_coupling_eval_v11.py -- protocol gears_coupling_v1.0 (+ v1.1 controls) on RPE1: scores of the downloaded GEARS
# predictions with the frozen S1, frozen per-split references and frozen cal/test units.
#   python gears_coupling_eval_v11.py <gears_results_dir> <out_dir>
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.dont_write_bytecode = True
import data_v07 as dv  # noqa: E402
import models_v04 as m4  # noqa: E402
from paths import MASTER_SEED, PROCESSED, RESULTS  # noqa: E402

K, N_PERM = 50, 20


def derangement(rng, n):
    while True:
        p = rng.permutation(n)
        if not (p == np.arange(n)).any():
            return p


def main(gdir, out):
    gdir, out = Path(gdir), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    pdir = PROCESSED / "v0.7" / "rpe1"
    U = pd.read_csv(pdir / "units_all_targets.csv")
    D = dv.load_dataset(U, pdir, pdir / "string_expdb_universe.parquet")
    U, Y, YB = D["U"], D["Y"], D["YB"]
    genes = pd.read_csv(pdir / "pseudobulk_genes.csv").gene_name.astype(str).tolist()
    rows, negs = [], []
    for f in sorted((RESULTS / "conf_v0.7" / "rpe1" / "cache").glob("split_*.pkl")):
        c = pickle.load(open(f, "rb"))
        i = int(c["split_id"])
        sd = gdir / f"split_{i:02d}"
        pt = json.loads((sd / "pred_targets.json").read_text())
        gj = json.loads((sd / "genes.json").read_text())
        assert [str(x) for x in gj["gene_name"]] == genes, f"gene order split {i}"
        assert pt["test_cal_and_test"] == list(c["target_calt"]), f"unit order split {i}"
        P = np.load(sd / "pred.npy").astype(np.float64)
        nv = len(pt["val"])
        Pc = P[nv:]
        assert Pc.shape[0] == len(c["target_calt"]) and np.isfinite(Pc).all()
        calt = np.concatenate([c["cal"], c["test"]])
        nc = len(c["cal"])
        ref = {"mu": c["ref_mu"], "sigma": c["ref_sigma"]}
        tc = [U.tcols.values[u] for u in calt]
        sg = m4.score_S1(Y[calt], Pc, ref, tc, K)
        sgh = m4.score_S1(YB[calt], Pc, ref, tc, K)
        cal = np.sort(sg[:nc])
        kq = int(np.ceil((nc + 1) * 0.9))
        qh = np.inf if kq > nc else cal[kq - 1]
        te = slice(nc, None)
        rows.append(pd.DataFrame({"split": i, "target": c["target_calt"][te], "tss_group": c["tss_calt"][te],
                                  "tier": c["tier_calt"][te], "tierA": c["tierA_calt"][te],
                                  "gears_full": sg[te], "gears_half": sgh[te], "gears_covered": sg[te] <= qh}))
        tidx, pt_ = calt[nc:], Pc[nc:]
        tct = tc[nc:]
        for r in range(N_PERM):
            p = derangement(np.random.default_rng([MASTER_SEED + 12000, i, r]), len(tidx))
            negs.append(pd.DataFrame({"perm": r, "split": i, "target": c["target_calt"][te],
                                      "neg_full": m4.score_S1(Y[tidx], pt_[p], ref, tct, K),
                                      "neg_half": m4.score_S1(YB[tidx], pt_[p], ref, tct, K)}))
        print(f"split {i:02d}", flush=True)
    pd.concat(rows, ignore_index=True).to_parquet(out / "gears_rows.parquet")
    pd.concat(negs, ignore_index=True).to_parquet(out / "gears_negative_rows.parquet")
    print("done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
