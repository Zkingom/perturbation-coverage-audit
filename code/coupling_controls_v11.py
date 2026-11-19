# coupling_controls_v11.py -- gears_coupling_v1.1 controls (blank = mean model, positive = split-half, negative =
# permuted linear-PCA predictions). Reads only frozen caches and processed data; never reads any GEARS output.
#   python coupling_controls_v11.py <dataset: rpe1|dev> <out_dir>
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.dont_write_bytecode = True
import eval_v06 as ev  # noqa: E402
import models_v04 as m4  # noqa: E402
from paths import MASTER_SEED, PROCESSED, RESULTS  # noqa: E402

N_PERM, N_BOOT, N_BOOT_PERM, K = 20, 2000, 500, 50


def load(ds):
    if ds == "rpe1":
        import data_v07 as dv
        pdir = PROCESSED / "v0.7" / "rpe1"
        U = pd.read_csv(pdir / "units_all_targets.csv")
        return dv.load_dataset(U, pdir, pdir / "string_expdb_universe.parquet"), RESULTS / "conf_v0.7" / "rpe1" / "cache"
    import pipeline_v06 as p6
    return p6.load_dev_v06(PROCESSED / "v0.6" / "k562ess"), RESULTS / "d1_v0.6.1" / "cache"


def derangement(rng, n):
    while True:
        p = rng.permutation(n)
        if not (p == np.arange(n)).any():
            return p


def main(ds, out):
    import pickle
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    D, cdir = load(ds)
    U, Y, YA, YB, y_ctrl, pri = D["U"], D["Y"], D["YA"], D["YB"], D["y_ctrl"], D["pri"]
    rows, negs, checks = [], [], []
    t0 = time.time()
    for f in sorted(cdir.glob("split_*.pkl")):
        c = pickle.load(open(f, "rb"))
        i = c["split_id"]
        calt = np.concatenate([c["cal"], c["test"]])
        assert (U.target.values[calt] == c["target_calt"]).all(), "unit order"
        nc = len(c["cal"])
        ref = {"mu": c["ref_mu"], "sigma": c["ref_sigma"]}
        tc = [U.tcols.values[u] for u in calt]
        mu = np.broadcast_to(c["ref_mu"], (len(calt), Y.shape[1]))
        EB = m4.score_S1(YB[calt], mu, ref, tc, K)
        ok_blank = np.array_equal(EB, c["scores"][("mean", K, "half")])
        Ef = m4.score_S1(Y[calt], mu, ref, tc, K)
        ok_full = np.array_equal(Ef, c["E_calt"])
        s_ctrl = m4.score_S1(YB[calt], YA[calt], ref, tc, K)
        Kc = int(c["chosen"]["K"])
        pred = ev.refit_predict(U, Y, y_ctrl, pri, np.asarray(c["pool"]), calt, "pca", Kc)
        s_pca_re = m4.score_S1(Y[calt], pred, ref, tc, K)
        dpca = float(np.max(np.abs(s_pca_re - c["scores"][("linear_pca", K, "full")])))
        checks.append({"split": int(i), "blank_half_equal": bool(ok_blank), "E_full_equal": bool(ok_full),
                       "pca_refit_max_abs_diff": dpca, "K": Kc})
        assert ok_blank and ok_full and dpca < 1e-8, checks[-1]
        te = slice(nc, None)
        q = m4  # noqa: F841
        # positive-control M0 coverage on half truth (cal -> test)
        cal_s = np.sort(s_ctrl[:nc])
        kq = int(np.ceil((nc + 1) * 0.9))
        qh = np.inf if kq > nc else cal_s[kq - 1]
        rows.append(pd.DataFrame({
            "split": i, "target": c["target_calt"][te], "tss_group": c["tss_calt"][te],
            "tier": c["tier_calt"][te], "tierA": c["tierA_calt"][te],
            "E_full": Ef[te], "E_B": EB[te],
            "pca_full": c["scores"][("linear_pca", K, "full")][te], "go_full": c["scores"][("linear_go", K, "full")][te],
            "pca_half": c["scores"][("linear_pca", K, "half")][te], "go_half": c["scores"][("linear_go", K, "half")][te],
            "ctrl_half": s_ctrl[te], "ctrl_covered": s_ctrl[te] <= qh}))
        tidx = calt[nc:]
        pt = pred[nc:]
        tct = tc[nc:]
        for r in range(N_PERM):
            p = derangement(np.random.default_rng([MASTER_SEED + 12000, int(i), r]), len(tidx))
            negs.append(pd.DataFrame({"perm": r, "split": i, "target": c["target_calt"][te], "tss_group": c["tss_calt"][te],
                                      "tier": c["tier_calt"][te], "tierA": c["tierA_calt"][te],
                                      "neg_full": m4.score_S1(Y[tidx], pt[p], ref, tct, K),
                                      "neg_half": m4.score_S1(YB[tidx], pt[p], ref, tct, K)}))
        print(f"split {i:02d} done {time.time() - t0:.0f}s", flush=True)
    A = pd.concat(rows, ignore_index=True)
    N = pd.concat(negs, ignore_index=True)
    A.to_parquet(out / "controls_rows.parquet")
    N.to_parquet(out / "negative_rows.parquet")
    json.dump(checks, open(out / "checks.json", "w"), indent=1)
    print("done", len(A), len(N), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
