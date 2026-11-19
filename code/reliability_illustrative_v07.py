# reliability_illustrative_v07.py -- ILLUSTRATIVE analysis (outside every protocol; user request 2026-10-15) of the
# split-half reliability of effect strength E_true. Two questions:
#   (a) cross-split effect strength: top-50 specific genes selected on half A, magnitude evaluated on half B (and
#       B->A; mean of the two). Selection noise no longer inflates the magnitude, so a cell-number-dependent noise
#       floor shared by both halves cannot masquerade as signal. Report rho(h, E_cross) vs the original rho(h, E_true).
#   (b) reliability of E_A vs E_B after partialling out log(number of cells) (rank-based partial correlation),
#       with the Spearman-Brown step-up.
# Pre-written reading (user): if rho(h, E_cross) on GWPS is close to 0.110 the 'prior lacks information' conclusion
# stands; if it is clearly higher (e.g. > 0.25) STOP and report -- the framing must be re-discussed.
#   python reliability_illustrative_v07.py [--boot 2000]
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, RESULTS  # noqa: E402
import explore_modularity_v07 as xm  # noqa: E402

OUT = RESULTS / "reliability_illustrative_v07"
SEED_OFF = {"dev": 0, "rpe1": 10, "gwps": 20}


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


def _sel_eval(Ysel, Yev, ref, tcols, k=50):
    """Top-k genes by |s| of Ysel (target columns excluded, as score_S1), RMS of s of Yev on those genes."""
    ss = (Ysel - ref["mu"]) / ref["sigma"]
    se = (Yev - ref["mu"]) / ref["sigma"]
    a = np.abs(ss)
    for i, cols in enumerate(tcols):
        a[i, cols] = -np.inf
    top = np.argsort(-a, axis=1)[:, :k]
    return np.sqrt((np.take_along_axis(se, top, axis=1) ** 2).mean(1))


def per_dataset(ds):
    import proxy_v06 as px
    cfg = xm.dataset_cfg(ds)
    pt = pd.read_parquet(cfg["res"] / "proxy_rows.parquet",
                         filters=[("kind", "==", "h"), ("model", "==", xm.MODEL), ("proxy", "==", xm.HCFG)])
    pt = pt.rename(columns={"value": "h"})
    D = xm.load_D(ds)
    Ud = D["U"]
    pos = {t: i for i, t in enumerate(Ud.target.values)}
    rows = []
    for f in sorted(Path(cfg["cache"]).glob("split_*.pkl")):
        c = pickle.load(open(f, "rb"))
        nt = len(c["test"])
        tg = c["target_calt"][-nt:]
        idx = np.array([pos[t] for t in tg])
        ref = {"mu": c["ref_mu"], "sigma": c["ref_sigma"]}
        tc = [Ud.tcols.values[i] for i in idx]
        YA, YB = np.asarray(D["YA"][idx]), np.asarray(D["YB"][idx])
        EA = _sel_eval(YA, YA, ref, tc)
        assert np.array_equal(EA, px.e_true(YA, ref, tc)), f"{f.name}: own selection/evaluation differs from e_true"
        assert np.array_equal(EA, c["EA_calt"][-nt:]), f"{f.name}: E_A differs from cache"
        EB = _sel_eval(YB, YB, ref, tc)
        xAB = _sel_eval(YA, YB, ref, tc)
        xBA = _sel_eval(YB, YA, ref, tc)
        rows.append(pd.DataFrame({"split": c["split_id"], "target": tg, "EA": EA, "EB": EB, "xAB": xAB, "xBA": xBA}))
    del D
    r = pd.concat(rows, ignore_index=True)
    r["Ecross"] = (r.xAB + r.xBA) / 2
    m = pt[["split", "target"]].merge(r, on=["split", "target"], how="left", validate="one_to_one")
    assert m.EA.notna().all(), f"{ds}: join failed"
    g = pt.groupby(["target", "tss_group"]).agg(h=("h", "mean"), E=("E_true", "mean")).reset_index()
    g = g.join(r.groupby("target")[["EA", "EB", "xAB", "xBA", "Ecross"]].mean(), on="target")
    nc = dict(zip(Ud.target.values, Ud.ncells.values))
    g["logn"] = np.log(g.target.map(nc).astype(float))
    assert g.notna().all().all()
    return g


def partial_rank_corr(a, b, z):
    ra, rb, rz = rankdata(a), rankdata(b), rankdata(z)
    X = np.column_stack([np.ones(len(rz)), rz])
    res = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]   # noqa: E731
    ea, eb = res(ra), res(rb)
    return float(np.corrcoef(ea, eb)[0, 1])


def stats(g):
    rho = spearmanr(g.h, g.E).correlation
    rx = spearmanr(g.h, g.Ecross).correlation
    rAB = spearmanr(g.EA, g.EB).correlation
    rxx = spearmanr(g.xAB, g.xBA).correlation
    pr = partial_rank_corr(g.EA, g.EB, g.logn)
    sb = lambda r: 2 * r / (1 + r)   # noqa: E731
    return {"rho_h_E": float(rho), "rho_h_Ecross": float(rx), "rho_h_Ecross_minus_rho": float(rx - rho),
            "rho_E_Ecross": float(spearmanr(g.E, g.Ecross).correlation),
            "r_EA_EB": float(rAB), "rel_SB": float(sb(rAB)),
            "r_xAB_xBA": float(rxx), "rel_cross_SB": float(sb(rxx)),
            "partial_r_EA_EB_given_logn": pr, "partial_rel_SB": float(sb(pr)),
            "rho_E_logn": float(spearmanr(g.E, g.logn).correlation),
            "rho_Ecross_logn": float(spearmanr(g.Ecross, g.logn).correlation)}


def boot(g, n, seed):
    rng = np.random.default_rng(seed)
    gs = g.tss_group.unique()
    posd = g.groupby("tss_group", sort=False).indices
    out = {}
    for b in range(n):
        ix = np.concatenate([posd[k] for k in rng.choice(gs, len(gs), replace=True)])
        for k, v in stats(g.take(ix).reset_index(drop=True)).items():
            out.setdefault(k, []).append(v)
    return {k: [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))] for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--datasets", nargs="*", default=["dev", "rpe1", "gwps"])
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    od = Path(a.out)
    od.mkdir(parents=True, exist_ok=True)
    res = {"label": "ILLUSTRATIVE (outside protocol; user request 2026-10-15)", "n_boot": a.boot, "datasets": {}}
    for ds in a.datasets:
        log(f"{ds}: effect strengths")
        g = per_dataset(ds)
        ref_rho = json.loads(Path(xm.dataset_cfg(ds)["ref_json"]).read_text(encoding="utf8"))["P"]["rho"]
        st = stats(g)
        assert st["rho_h_E"] == ref_rho, f"{ds}: rho mismatch {st['rho_h_E']} vs {ref_rho}"
        log(f"{ds}: bootstrap")
        st["ci"] = boot(g, a.boot, MASTER_SEED + 9200 + SEED_OFF[ds])
        st["n_genes"] = int(len(g))
        res["datasets"][ds] = st
        g.to_csv(od / f"gene_level_{ds}.csv", index=False)
    if "gwps" in res["datasets"]:
        v = res["datasets"]["gwps"]["rho_h_Ecross"]
        res["reading_c"] = {"gwps_rho_h_Ecross": v, "stop_rule_gt_0.25": bool(v > 0.25)}
    (od / "reliability_illustrative_results.json").write_text(json.dumps(res, indent=1), encoding="utf8")
    log("done" + (" -- STOP RULE TRIGGERED (GWPS rho_h_Ecross > 0.25)" if res.get("reading_c", {}).get("stop_rule_gt_0.25") else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
