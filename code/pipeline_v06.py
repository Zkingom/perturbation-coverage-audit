# pipeline_v06.py -- dataset-level driver for protocol v0.6.1.
#   stage 1  compute_split for all 21 splits (repeated grouped 7-fold), cached per split
#   stage 2  rule 7.1 (h configuration) and 7.1b (d / s learner) on the primary model -> nominal-alpha rows
#            for every model x k (+ split-half variant at k = 50); width-matched rows (test-unit and
#            calibration-unit matching) for the primary and secondary model at k = 50
#   stage 3  rules 7.2 / 7.3 on matched rows (dev) or the frozen choice (confirmation); best M_dist;
#            permuted-h control (20 runs, nominal and matched) for the primary method
#   stage 4  summaries, confounding, incremental value, criteria (on dev only as a code check)
# load_dev_v06() builds the development dataset; it is NOT run in D0' (synthetic data only).
import json
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd

import analysis_v06 as an
import eval_v06 as ev
from splits_v06 import repeated_group_kfold


def stage1(D, seed_base, cache_dir, n_rounds=3, log=print):
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    splits = repeated_group_kfold(D["U"].tss_group.values[D["E_idx"]], seed_base, n_rounds=n_rounds)
    for i, sp in enumerate(splits):
        f = cache_dir / f"split_{i:02d}.pkl"
        if f.exists():
            continue
        t0 = time.time()
        c = ev.compute_split(i, sp, D)
        tmp = f.with_suffix(".tmp")
        with open(tmp, "wb") as fh:
            pickle.dump(c, fh)
        tmp.replace(f)
        log(f"stage1 split {i:02d} round {sp['round']} fold {sp['fold']} K={c['chosen']['K']} d={c['chosen']['d']} "
            f"n_test={len(c['test'])} {time.time() - t0:.1f}s")
    return splits


def load_caches(cache_dir):
    out = []
    for f in sorted(Path(cache_dir).glob("split_*.pkl")):
        with open(f, "rb") as fh:
            out.append(pickle.load(fh))
    return out


def _hyper_value(meth, hy):
    from methods_v06 import METHODS
    fam = METHODS[meth][0]
    if fam == "mondrian":
        return int(hy)
    if fam in ("gcc", "gcc2"):
        return hy
    return float(hy)


def run_dataset(D, seed_base, out_dir, n_rounds=3, frozen=None, log=print, ext_dir=None):
    """frozen = None (development: rules 7.1-7.3) or {h, d, s, hypers, primary, best_dist} (confirmation)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stage1(D, seed_base, out / "cache", n_rounds, log)
    caches = load_caches(out / "cache")
    if ext_dir is not None:                  # v0.7 development additions (ext_v07dev.py)
        import ext_v07dev as ex
        for c in caches:
            z = np.load(Path(ext_dir) / f"ext_split_{c['split_id']:02d}.npz")
            ex.attach(c, {k: z[k] for k in z.files})
        log(f"attached ext arrays from {ext_dir}")
    pt = pd.concat([ev.proxy_table(c) for c in caches], ignore_index=True)
    pt.to_parquet(out / "proxy_rows.parquet")
    P = ev.PRIMARY_MODEL
    if frozen is None:
        hcfg, rho_h = an.select_proxy(pt[(pt.model == P) & (pt.kind == "h")])
        dln, rho_d = an.select_learner(pt[(pt.model == P) & (pt.kind == "d") & (pt.k == ev.PRIMARY_K)])
        sln, rho_s = an.select_learner(pt[(pt.model == P) & (pt.kind == "s") & (pt.k == ev.PRIMARY_K)])
    else:
        hcfg, dln, sln, rho_h, rho_d, rho_s = frozen["h"], frozen["d"], frozen["s"], None, None, None
    sel = {"h": tuple(hcfg.split("|")), "d": dln, "s": sln}
    log(f"selected proxies: h={hcfg} d={dln} s={sln}")
    parts = []
    t0 = time.time()
    for c in caches:
        for m in ev.MODELS:
            for k in ev.KS:
                parts.append(ev.apply_methods(c, m, k, sel))
            parts.append(ev.apply_methods(c, m, ev.PRIMARY_K, sel, variant="half"))
        for m in (P, ev.SECONDARY_MODEL):
            parts.append(ev.apply_methods(c, m, ev.PRIMARY_K, sel, mode="matched"))
        parts.append(ev.apply_methods(c, P, ev.PRIMARY_K, sel, mode="matched_cal"))
    cov = pd.concat(parts, ignore_index=True)
    log(f"stage 2 rows {len(cov)} in {time.time() - t0:.0f}s")
    pm = lambda d, mode, model=P: d[(d.model == model) & (d.k == ev.PRIMARY_K) & (d.variant == "full")   # noqa: E731
                                    & (d["mode"] == mode) & (d.perm == -1)]
    if frozen is None:
        hyp = an.select_hyper(pm(cov, "matched"))
        prim = an.select_primary(hyp)
        best_dist = an.select_primary(hyp, [b for b in an.BASELINES if b in hyp], tie=0.0)
    else:
        hyp = {m: {"chosen": h} for m, h in frozen["hypers"].items()}
        prim, best_dist = frozen["primary"], frozen["best_dist"]
    prim_h = hyp[prim]["chosen"]
    log(f"primary {prim} {prim_h}; best M_dist {best_dist} {hyp[best_dist]['chosen']}")
    perm = []
    for c in caches:
        for r in range(ev.N_PERM):
            for mode in ("nominal", "matched"):
                perm.append(ev.apply_methods(c, P, ev.PRIMARY_K, sel, grid=[(prim, _hyper_value(prim, prim_h))],
                                             mode=mode, perm_r=r))
    cov = pd.concat([cov] + perm, ignore_index=True)
    cov.to_parquet(out / "coverage_rows.parquet")
    covs = pd.concat([c["covariates"].assign(split=c["split_id"]) for c in caches], ignore_index=True)
    covs.to_parquet(out / "covariates_analysis_only.parquet")
    summ, _, _ = summarise(out, cov, pt, covs, hcfg, rho_h, dln, rho_d, sln, rho_s, hyp, prim, prim_h, best_dist,
                           [c["chosen"]["K"] for c in caches], [c["chosen"]["d"] for c in caches])
    return summ, cov, pt


def summarise(out, cov, pt, covs, hcfg, rho_h, dln, rho_d, sln, rho_s, hyp, prim, prim_h, best_dist, chosen_K, chosen_d):
    """Stage 4 of run_dataset (summary statistics, criteria). Split out on 2026-10-10 so that it can be rerun
    from saved coverage rows; the computation is unchanged."""
    out = Path(out)
    P = ev.PRIMARY_MODEL
    pm = lambda d, mode, model=P: d[(d.model == model) & (d.k == ev.PRIMARY_K) & (d.variant == "full")   # noqa: E731
                                    & (d["mode"] == mode) & (d.perm == -1)]
    # stage 4
    chosen = lambda d, meth: d[(d.method == meth) & (d.hyper == hyp[meth]["chosen"])]        # noqa: E731
    summ = {"selection": {"h": hcfg, "rho_h": rho_h, "d": dln, "rho_d": rho_d, "s": sln, "rho_s": rho_s,
                          "hypers": {m: v["chosen"] for m, v in hyp.items()}, "hyper_detail": hyp,
                          "primary": [prim, prim_h], "best_dist": [best_dist, hyp[best_dist]["chosen"]]},
            "chosen_K": chosen_K, "chosen_d": chosen_d}
    tabs = {}
    for (m, k, var, mode, meth, hy), d in cov[cov.perm == -1].groupby(["model", "k", "variant", "mode", "method", "hyper"]):
        tabs[f"{m}|{k}|{var}|{mode}|{meth}|{hy}"] = an.tier_table(d)
    summ["tier_tables"] = tabs
    summ["monotone_all"] = bool(cov.monotone.all())
    import methods_v06 as mv
    summ["m3_solver_fallbacks"] = dict(mv.M3_FALLBACKS)
    m_rows = pm(cov, "matched")
    m0n = pm(cov, "nominal")
    summ["dgo_strata"] = {f"{mode}|{meth}": an.strata_table(chosen(pm(cov, mode), meth) if meth != "M0"
                                                              else pm(cov, mode)[pm(cov, mode).method == "M0"], "D_GO")
                          for mode in ("nominal", "matched")
                          for meth in ["M0", "M1", prim, "Mscore-N"] + [b for b in an.BASELINES if b in hyp]}
    # matched paired comparisons
    d0m = m_rows[m_rows.method == "M0"].sort_values(["split", "target"])
    comp = {}
    for meth in [prim, best_dist, "Mcomb-N", "Mcomb-3", "Mscore-N", "oracle"]:
        dm = chosen(m_rows, meth).sort_values(["split", "target"])
        comp[f"{meth}_vs_M0"] = an.paired_T5(dm, d0m)
    dbd = chosen(m_rows, best_dist).sort_values(["split", "target"])
    for meth in [prim, "Mcomb-N", "Mcomb-3", "Mscore-N"]:
        comp[f"{meth}_vs_{best_dist}"] = an.paired_T5(chosen(m_rows, meth).sort_values(["split", "target"]), dbd)
    bases = [b for b in an.BASELINES if b in hyp]
    for b in bases:                          # v0.7 development: every baseline vs M0, and prim / Mscore-N vs it
        db = chosen(m_rows, b).sort_values(["split", "target"])
        comp[f"{b}_vs_M0"] = an.paired_T5(db, d0m)
        for meth in sorted({prim, "Mscore-N"}):
            if meth != b:
                comp[f"{meth}_vs_{b}"] = an.paired_T5(chosen(m_rows, meth).sort_values(["split", "target"]), db)
    summ["matched_paired_T5"] = comp
    dec = {}
    for meth in [prim, "Mscore-N"] + bases:
        dec[f"{meth}_vs_M0"] = an.decision_paired(chosen(m_rows, meth), d0m)
    for meth in sorted({prim, "Mscore-N"}):
        dec[f"{meth}_vs_{best_dist}"] = an.decision_paired(chosen(m_rows, meth), chosen(m_rows, best_dist))
    summ["decision_paired_top10_matched"] = dec
    # confounding and incremental (M0 nominal, k = 50, linear models)
    rows_lin = m0n.iloc[:0]
    for m in ev.LINEAR:
        rows_lin = pd.concat([rows_lin, pm(cov, "nominal", m)[pm(cov, "nominal", m).method == "M0"]])
    summ["confounding"] = an.confounding(an.confounding_frame(rows_lin, covs))
    hr = cov[(cov.k == ev.PRIMARY_K) & (cov.variant == "half") & (cov["mode"] == "nominal") & (cov.method == "M0")
             & cov.model.isin(ev.LINEAR) & (cov.perm == -1)]
    summ["confounding_half"] = an.confounding(an.confounding_frame(hr, covs, half=True))
    summ["incremental"] = an.incremental(rows_lin, covs)
    # criteria (nominal as in v0.6, and width-matched)
    pth = pt[(pt.kind == "h") & (pt.proxy == hcfg)].rename(columns={"value": "h"})
    for mode in ("nominal", "matched"):
        cm = cov[(cov.variant == "full") & (cov["mode"] == mode)]
        summ[f"criteria_{mode}"] = an.criteria(cm, pth, prim, prim_h, P, ev.PRIMARY_K)
    summ["decision_metrics_matched"] = {"M0": an.decision_metrics(d0m), prim: an.decision_metrics(chosen(m_rows, prim))}
    summ["gene_level_matched"] = {"M0": an.gene_level(d0m), prim: an.gene_level(chosen(m_rows, prim))}
    summ["tier_stability_frac_genes_changing"] = an.tier_stability(m0n[m0n.method == "M0"])
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=_js), encoding="utf8")
    return summ, cov, pt


def load_dev_v06(pb_dir):
    """Development dataset (K562 essential) for D1: v0.5 loader + split halves from pseudobulk_v06 output
    in pb_dir + number of GO annotations per target. Written in D0'; first run in D1a."""
    import pipeline_v05 as pl
    from paths import raw
    U, Y, y_ctrl, genes, pri = pl.load()
    pb_dir = Path(pb_dir)
    rows = pd.read_csv(pb_dir / "pseudobulk_rows.csv")
    gv = pd.read_csv(pb_dir / "pseudobulk_genes.csv").gene_name.values
    assert (gv == genes).all(), "gene order differs between v0.4 and v0.6 pseudobulk"
    rix = {g: i for i, g in enumerate(rows.gene_transcript)}
    r = U.gene_transcript.map(rix).values
    YA = np.load(pb_dir / "pseudobulk_logmean_A.npy").astype(np.float64)[r]
    YB = np.load(pb_dir / "pseudobulk_logmean_B.npy").astype(np.float64)[r]
    with open(raw("gene2go_all.pkl"), "rb") as fh:
        gene2go = pickle.load(fh)
    pri["n_go"] = np.array([len(gene2go.get(t, ())) for t in U.target.values], float)
    return {"U": U, "Y": Y, "YA": YA, "YB": YB, "y_ctrl": y_ctrl, "pri": pri,
            "E_idx": np.where(U.E.values)[0], "AT_idx": np.where(U.always_train.values)[0]}


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, tuple):
        return list(o)
    return str(o)
