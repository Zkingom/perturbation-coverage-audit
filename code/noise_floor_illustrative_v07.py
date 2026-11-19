# noise_floor_illustrative_v07.py -- ILLUSTRATIVE analysis (outside every protocol; user request 2026-10-16).
# Does the coverage failure along effect strength survive correction for the cell-number-dependent noise floor?
#  (a) Null: pseudo-perturbations of n cells drawn without replacement from non-targeting control cells
#      (extract_controls_v07.py; control cells are never part of the training-perturbation reference mu / sigma,
#      which is built from perturbation pseudobulks only). Their effect strength is computed with exactly the
#      E_true definition (proxy_v06.e_true: sigma-scaled RMS on the top-50 |y - mu| / sigma genes) using each
#      split's own reference; no target column exists for controls, so none is excluded. Grid of n (geometric,
#      30 .. max cells of E units), R draws per n; mean m(n) and SD s(n) per split, interpolated linearly in log n.
#  (b) E_corr = E_true - m(n_cells);  E_z = (E_true - m(n)) / s(n).
#  (c) rho(E_corr, log cells); tiers of E_corr (same quintile rule as eval_v06.tiers_of, over all E units of the
#      split) and M0 coverage by those tiers (same GEE as the confirmatory B0); rho(h, E_corr) at gene level.
#  (d) Pre-written reading (user): GWPS M0 T5 coverage under E_corr tiers > 0.75 -> failure mainly attributable to
#      measurement precision on GWPS -> stop and report; otherwise keep the effect-strength framing and report the
#      cell-number association as a limitation.
#   python noise_floor_illustrative_v07.py --dataset {dev,rpe1,gwps} [--reps 200] [--grid 20]
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED, RESULTS  # noqa: E402
import explore_modularity_v07 as xm  # noqa: E402

OUT = RESULTS / "noise_floor_illustrative_v07"
SEED = MASTER_SEED + 9600
CH = 256


def log(m):
    print(time.strftime("%H:%M:%S"), m, flush=True)


def e_rows(Ym, ref, tcols):
    import proxy_v06 as px
    out = np.empty(len(Ym))
    for s in range(0, len(Ym), CH):
        out[s:s + CH] = px.e_true(np.asarray(Ym[s:s + CH], dtype=np.float64), ref, tcols[s:s + CH])
    return out


def main():
    import analysis_v06 as an
    import eval_v06 as ev
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=["dev", "rpe1", "gwps"])
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--grid", type=int, default=20)
    a = ap.parse_args()
    ds = a.dataset
    od = OUT / ds
    od.mkdir(parents=True, exist_ok=True)
    cfg = xm.dataset_cfg(ds)
    D = xm.load_D(ds)
    U, Y, E_idx = D["U"], D["Y"], D["E_idx"]
    nE = U.ncells.values[E_idx].astype(float)
    Xc = np.load(PROCESSED / "noise_floor_v07" / ds / "ctrl_x.npy")          # float32, cells x genes
    assert Xc.shape[1] == Y.shape[1]
    nmax = int(min(nE.max(), len(Xc) - 1))
    grid = np.unique(np.round(np.geomspace(30, nmax, a.grid)).astype(int))
    log(f"{ds}: {len(Xc)} control cells, grid {grid.tolist()}, reps {a.reps}")
    rng = np.random.default_rng(SEED)
    M = np.empty((len(grid) * a.reps, Y.shape[1]))
    k = 0
    for n in grid:
        for _ in range(a.reps):
            idx = rng.choice(len(Xc), n, replace=False)
            M[k] = Xc[idx].astype(np.float64).mean(0).astype(np.float32)   # pseudobulks are stored as float32
            k += 1
    del Xc
    gcol = np.repeat(grid, a.reps)
    empty = [np.array([], dtype=int)] * len(M)
    pos = {t: i for i, t in enumerate(U.target.values)}
    tcE = [U.tcols.values[i] for i in E_idx]
    rows, curves, tier_rows = [], {}, []
    files = sorted(Path(cfg["cache"]).glob("split_*.pkl"))
    for f in files:
        c = pickle.load(open(f, "rb"))
        sid = int(c["split_id"])
        ref = {"mu": c["ref_mu"], "sigma": c["ref_sigma"]}
        En = e_rows(M, ref, empty)
        mu_n = np.array([En[gcol == n].mean() for n in grid])
        sd_n = np.array([En[gcol == n].std(ddof=1) for n in grid])
        curves[sid] = {"mean": mu_n.tolist(), "sd": sd_n.tolist()}
        Eall = e_rows(Y[E_idx], ref, tcE)
        ln = np.log(np.clip(nE, grid[0], grid[-1]))
        m_ = np.interp(ln, np.log(grid), mu_n)
        s_ = np.interp(ln, np.log(grid), sd_n)
        Ec, Ez = Eall - m_, (Eall - m_) / s_
        tg = U.target.values[E_idx]
        t_orig = ev.tiers_of(Eall, 5)
        t_c = ev.tiers_of(Ec, 5)
        t_z = ev.tiers_of(Ez, 5)
        calt = c["target_calt"]
        ix = np.array([np.flatnonzero(tg == t)[0] for t in calt])
        assert np.array_equal(Eall[ix], c["E_calt"]), f"split {sid}: E_true differs from cache"
        assert np.array_equal(t_orig[ix], c["tier_calt"]), f"split {sid}: tiers differ from cache"
        nt = len(c["test"])
        tier_rows.append(pd.DataFrame({"split": sid, "target": calt[-nt:], "tier_orig": t_orig[ix][-nt:],
                                       "tier_corr": t_c[ix][-nt:], "tier_z": t_z[ix][-nt:]}))
        rows.append(pd.DataFrame({"split": sid, "target": tg, "E": Eall, "E_corr": Ec, "E_z": Ez, "ncells": nE}))
        log(f"split {sid:02d}: null mean at n={grid[0]}: {mu_n[0]:.3f}, at n={grid[-1]}: {mu_n[-1]:.3f}")
    del M, D
    allr = pd.concat(rows, ignore_index=True)
    tr = pd.concat(tier_rows, ignore_index=True)
    # M0 coverage rows (primary model, k = 50, randomised matching, unpermuted) exactly as the confirmatory B0
    rm = pd.read_parquet(cfg["rand"])
    d0 = rm[(rm.method == "M0") & (rm.model == ev.PRIMARY_MODEL) & (rm.k == ev.PRIMARY_K)]
    if "variant" in d0:
        d0 = d0[d0.variant == "full"]
    if "perm" in d0:
        d0 = d0[(d0.perm.isna()) | (d0.perm == -1)]
    d0 = d0.merge(tr, on=["split", "target"], how="left", validate="one_to_one")
    assert d0.tier_orig.notna().all() and (d0.tier_orig.astype(int) == d0.tier.astype(int)).all(), "tier mismatch"
    ref_json = json.loads(Path(cfg["ref_json"]).read_text(encoding="utf8"))
    g_orig = an.gee_tiers(d0.assign(tier=d0.tier_orig.astype(int)))
    b0_ref = ref_json["B0"]["M0_T5"]
    repro = abs(g_orig[5]["cov"] - b0_ref["cov"]) < 1e-12 and abs(g_orig[5]["hi"] - b0_ref["hi"]) < 1e-12
    g_corr = an.gee_tiers(d0.assign(tier=d0.tier_corr.astype(int)))
    g_z = an.gee_tiers(d0.assign(tier=d0.tier_z.astype(int)))
    # gene-level correlations (mean over splits)
    gl = allr.groupby("target").agg(E=("E", "mean"), E_corr=("E_corr", "mean"), E_z=("E_z", "mean"),
                                    ncells=("ncells", "first")).reset_index()
    pt = pd.read_parquet(cfg["res"] / "proxy_rows.parquet",
                         filters=[("kind", "==", "h"), ("model", "==", xm.MODEL), ("proxy", "==", xm.HCFG)])
    hg = pt.groupby("target").value.mean().rename("h")
    gl = gl.join(hg, on="target")
    gt = gl.dropna(subset=["h"])
    rho = lambda x, y: float(spearmanr(x, y).correlation)   # noqa: E731
    t5o = set(map(tuple, tr[tr.tier_orig == 5][["split", "target"]].values))
    t5c = set(map(tuple, tr[tr.tier_corr == 5][["split", "target"]].values))
    res = {"label": "ILLUSTRATIVE (outside protocol; user request 2026-10-16)", "dataset": ds,
           "n_ctrl_cells_used": None, "grid": grid.tolist(), "reps": a.reps, "seed": SEED,
           "null_curves_by_split": curves,
           "rho_logn": {"E": rho(gl.E, np.log(gl.ncells)), "E_corr": rho(gl.E_corr, np.log(gl.ncells)),
                        "E_z": rho(gl.E_z, np.log(gl.ncells))},
           "rho_h": {"E": rho(gt.h, gt.E), "E_corr": rho(gt.h, gt.E_corr), "E_z": rho(gt.h, gt.E_z), "n_genes": int(len(gt))},
           "rho_E_Ecorr": rho(gl.E, gl.E_corr),
           "frac_E_corr_le_0": float((gl.E_corr <= 0).mean()),
           "T5_overlap_orig_corr": len(t5o & t5c) / max(1, len(t5o)),
           "M0_tiers_orig": g_orig, "M0_tiers_orig_reproduces_confirmatory_B0": bool(repro),
           "M0_tiers_corr": g_corr, "M0_tiers_z": g_z}
    if ds == "gwps":
        res["reading_d"] = {"M0_T5_corr": g_corr[5]["cov"], "stop_rule_gt_0.75": bool(g_corr[5]["cov"] > 0.75)}
    info = json.loads((PROCESSED / "noise_floor_v07" / ds / "ctrl_cells.json").read_text())
    res["n_ctrl_cells_used"] = info["n_ctrl_kept"]
    (od / "noise_floor_results.json").write_text(json.dumps(res, indent=1), encoding="utf8")
    gl.to_csv(od / "gene_level.csv", index=False)
    log(f"done T5 orig {g_orig[5]['cov']:.3f} corr {g_corr[5]['cov']:.3f} z {g_z[5]['cov']:.3f}"
        + (" -- STOP RULE TRIGGERED (GWPS M0 T5 under corrected tiers > 0.75)" if res.get("reading_d", {}).get("stop_rule_gt_0.75") else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
