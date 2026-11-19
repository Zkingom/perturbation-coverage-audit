# tests_v06.py -- D0 checks for protocol v0.6 on synthetic data only (no confirmation data, no dev data).
#   python tests_v06.py unit       fast unit tests (splits, methods, M3, analysis, proxy canary, pseudobulk)
#   python tests_v06.py e2e <pos|neg> <outdir>   end-to-end positive / negative control (21 splits)
#   python tests_v06.py e2e_v07 <pos|neg> <outdir>   v0.7 confirmation path (frozen selection, new synthetic set)
#   (t_v07_data reads the K562 dev bulk and pseudobulk to check the generic loader; no confirmation data)
# Results are written as JSON next to the outputs; every check is an assert with a message.
import json
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis_v06 as an  # noqa: E402
import eval_v06 as ev  # noqa: E402
import methods_v06 as mv  # noqa: E402
import pipeline_v06 as p6  # noqa: E402
import pseudobulk_v06 as pb  # noqa: E402
import synth_v06  # noqa: E402
from splits_v06 import repeated_group_kfold  # noqa: E402

RES = {}


def check(name, cond, detail=None):
    RES[name] = {"pass": bool(cond), "detail": detail}
    print(("PASS " if cond else "FAIL ") + name + ("" if detail is None else f"  {detail}"), flush=True)
    return bool(cond)


# ------------------------------------------------------------------ splits

def t_splits():
    rng = np.random.default_rng(0)
    gid = np.array([f"g{x}" for x in rng.integers(0, 1700, 1786)])
    sps = repeated_group_kfold(gid, ev.DEV_SEED_BASE)
    n = len(gid)
    cnt = np.zeros(n, int)
    ok_part, ok_group, sizes = True, True, []
    for s in sps:
        allu = np.concatenate([s[k] for k in ("train", "val", "cal", "test")])
        ok_part &= len(allu) == n and len(np.unique(allu)) == n
        fold_of = np.empty(n, object)
        for k in ("train", "val", "cal", "test"):
            fold_of[s[k]] = k
        ok_group &= bool(pd.DataFrame({"g": gid, "f": fold_of}).groupby("g").f.nunique().max() == 1)
        cnt[s["test"]] += 1
        sizes.append((len(s["test"]), len(s["val"]), len(s["cal"])))
    sz = np.array(sizes)
    check("splits: 21 splits", len(sps) == 21)
    check("splits: each split is a partition", ok_part)
    check("splits: TSS groups never span folds", ok_group)
    check("splits: every unit is a test unit exactly 3 times", (cnt == 3).all(), {"min": int(cnt.min()), "max": int(cnt.max())})
    # greedy whole-group filling stops a fold once the next group does not fit, so val / cal can fall a
    # few units short of their quota (same rule as v0.5's random_group_split); allowed shortfall <= 3%.
    check("splits: test ~ n/7; val, cal within 3% below floor(0.10 n), floor(0.15 n)",
          abs(sz[:, 0].mean() - n / 7) < 2 and (sz[:, 1] <= int(0.10 * n)).all() and (sz[:, 1] >= 0.97 * int(0.10 * n)).all()
          and (sz[:, 2] <= int(0.15 * n)).all() and (sz[:, 2] >= 0.97 * int(0.15 * n)).all(),
          {"test": [int(sz[:, 0].min()), int(sz[:, 0].max())], "val": [int(sz[:, 1].min()), int(sz[:, 1].max())],
           "cal": [int(sz[:, 2].min()), int(sz[:, 2].max())]})
    sps2 = repeated_group_kfold(gid, ev.DEV_SEED_BASE)
    check("splits: deterministic given seed", all((a["test"] == b["test"]).all() and (a["cal"] == b["cal"]).all()
                                                  for a, b in zip(sps, sps2)))


# ------------------------------------------------------------------ methods

def t_methods():
    rng = np.random.default_rng(1)
    n = 267
    s = rng.exponential(1, n)
    h = rng.uniform(0.5, 3, n)
    one = lambda x: np.ones((len(np.ravel(x)), 1))                       # noqa: E731
    t = mv.m3_threshold(one(h), s, np.ones(1), 0.9 - 1e-9)
    check("M3: Phi = {1} at U -> 1 - alpha equals split-conformal qhat", abs(t - mv.qhat(s)) < 1e-8,
          {"m3": t, "qhat": mv.qhat(s)})
    # one-LP threshold equals bisection on the dual (quad and spline bases)
    worst = 0.0
    for kind in ("quad", "ns3", "ns5"):
        f = mv.phi_basis(kind, rng.uniform(0.5, 3, 500))
        Pc = f(h)
        for ht, U in [(0.7, 0.3), (1.5, 0.85), (2.8, 0.5), (1.0, -0.05)]:
            pt = f(np.array([ht]))[0]
            a = mv.m3_threshold(Pc, s, pt, U)
            lo, hi = 0.0, 4 * s.max() + 1
            for _ in range(60):
                mid = (lo + hi) / 2
                if mv._eta_test(Pc, s, pt, mid, 0.1) < U:
                    lo = mid
                else:
                    hi = mid
            if np.isfinite(a):
                worst = max(worst, abs(a - lo) / max(1, lo))
    check("M3: one-LP threshold equals dual bisection", worst < 1e-6, {"max_rel_diff": worst})
    # coverage guarantees over replicates; heteroscedastic scores, true scale observed
    reps = 60
    acc = {m: [] for m in ("M0", "M1", "M2", "M3")}
    wsum = []
    for r in range(reps):
        g = np.random.default_rng(100 + r)
        ho = g.uniform(0.5, 3, 800)
        hc, ht = g.uniform(0.5, 3, 267), g.uniform(0.5, 3, 255)
        sc, st = hc * g.exponential(1, 267), ht * g.exponential(1, 255)
        for nm, hy in (("M0", None), ("M1", 0.0), ("M2", 5), ("M3", "quad")):
            _, cv = mv.run_method(nm, hy, sc, st, hc, ht, ho, seed=r)
            top = ht >= np.quantile(ht, 0.8)
            acc[nm].append((cv.mean(), cv[top].mean()))
            if nm == "M3":
                z = (ht - ho.mean()) / ho.std()
                wsum.append(((cv - 0.9) * z).mean())
    tab = {m: [round(float(np.mean([x[0] for x in v])), 3), round(float(np.mean([x[1] for x in v])), 3)] for m, v in acc.items()}
    check("methods: marginal coverage within 0.90 +/- 0.015 (60 replicates)",
          all(abs(v[0] - 0.9) < 0.015 for v in tab.values()), tab)
    check("methods: M0 under-covers the top-20% h; M1/M2/M3 near 0.90 there",
          tab["M0"][1] < 0.85 and all(abs(tab[m][1] - 0.9) < 0.03 for m in ("M1", "M2", "M3")), tab)
    m, se = float(np.mean(wsum)), float(np.std(wsum) / np.sqrt(reps))
    check("M3: E[(cov - 0.9) * h] = 0 (exact guarantee for f in Phi)", abs(m) < 3 * se + 1e-3, {"mean": m, "se": se})
    # far-extrapolated test points: no exception; the infeasible-LP branch returns +/-inf. Recorded, not
    # asserted: an extrapolated basis can also give a finite negative threshold (empty set).
    f = mv.phi_basis("ns5", h)
    tt, cv = mv.m3(s, h, np.array([1.0, 1.0, 1.0]), np.array([50.0, -40.0, 1.0]), f, seed=3)
    one = lambda x: np.ones((len(np.ravel(x)), 1))                        # noqa: E731
    Pc = np.column_stack([np.ones(n), h])
    t_inf = mv.m3_threshold(Pc, s, np.array([1.0, 1e6]), 0.5)
    check("M3: extrapolated test points run without error; infeasible branch returns +/-inf",
          np.isfinite(tt[2]) and np.isinf(t_inf), {"thr_h=50,-40,1": tt.tolist(), "thr_h=1e6_linear": t_inf})
    # D0' regression case: dual simplex 'Not Set' failure (ns5, alpha = 0.01) -> interior-point fallback;
    # result must agree with bisection on the dual and with the attainable-range rule
    fx = Path(__file__).resolve().parents[1] / "results" / "d0p_v0.6.1" / "m3_fallback_case_split02.npz"
    if fx.exists():
        z = np.load(fx)
        t_fb = mv.m3_threshold(z["Phi_cal"], z["s_cal"], z["phi_t"], float(z["U"]), float(z["alpha"]))
        emax = mv._eta_range_max(z["Phi_cal"], z["phi_t"], float(z["alpha"]))
        check("M3: solver-failure fallback gives the bisection/range answer (+inf here)",
              np.isposinf(t_fb) and float(z["U"]) > emax, {"thr": t_fb, "U": float(z["U"]), "eta_max": emax})
    # stress: small alpha, spline bases, extrapolated points -> every threshold equals dual bisection
    worst, n_cmp = 0.0, 0
    g = np.random.default_rng(77)
    for kind in ("ns5", "ns3"):
        for a_ in (0.01, 0.02):
            hh = g.uniform(0.5, 3, 106)
            ss = hh * g.exponential(1, 106)
            f = mv.phi_basis(kind, g.uniform(0.5, 3, 500))
            Pc_ = f(hh)
            for ht in (0.2, 0.5, 1.5, 3.0, 3.6):
                Uu = g.uniform(-a_, 1 - a_)
                tt = mv.m3_threshold(Pc_, ss, f(np.array([ht]))[0], Uu, a_)
                pt_s = f(np.array([ht]))[0]
                lo, hi = -4 * ss.max() - 1, 4 * ss.max() + 1
                if Uu > mv._eta_range_max(Pc_, pt_s, a_):
                    ref_ = np.inf
                else:
                    while mv._eta_test(Pc_, ss, pt_s, hi, a_) < Uu:
                        hi *= 2
                    for _ in range(60):
                        mid = (lo + hi) / 2
                        if mv._eta_test(Pc_, ss, f(np.array([ht]))[0], mid, a_) < Uu:
                            lo = mid
                        else:
                            hi = mid
                    ref_ = lo
                n_cmp += 1
                if np.isinf(ref_) or np.isinf(tt):
                    worst = max(worst, 0.0 if (np.isinf(ref_) and np.isinf(tt)) else np.inf)
                elif ref_ > -4 * ss.max():
                    worst = max(worst, abs(tt - ref_) / max(1, abs(ref_)))
    check("M3: small-alpha / extrapolation stress test agrees with dual bisection", worst < 1e-6,
          {"n": n_cmp, "max_rel": worst, "fallbacks": dict(mv.M3_FALLBACKS)})
    # M2: bins with too few calibration units give +inf thresholds
    q = mv.qhat(np.arange(5.0))
    check("qhat: +inf when ceil((n+1)(1-alpha)) > n", np.isinf(q))


# ------------------------------------------------------------------ analysis

def t_analysis():
    rng = np.random.default_rng(5)
    rows = []
    for sp in range(21):
        for gi in range(400):
            tier = gi % 5 + 1
            p0 = [0.97, 0.95, 0.92, 0.88, 0.70][tier - 1]
            rows.append({"split": sp, "target": f"g{gi}", "tss_group": f"t{gi // 2}", "tier": tier,
                         "covered": rng.uniform() < p0, "thr": 1.0})
    d0 = pd.DataFrame(rows)
    g = an.gee_tiers(d0)
    check("GEE: tier coverage estimates and CIs cover the truth",
          all(g[t]["lo"] <= p <= g[t]["hi"] for t, p in zip(range(1, 6), [0.97, 0.95, 0.92, 0.88, 0.70])),
          {t: [round(g[t]["cov"], 3), round(g[t]["lo"], 3), round(g[t]["hi"], 3)] for t in g})
    d1 = d0.copy()
    m5 = d1.tier == 5
    d1.loc[m5, "covered"] = d1.loc[m5, "covered"] | (rng.uniform(size=m5.sum()) < 0.5)   # T5 -> ~0.85
    pr = an.paired_T5(d1, d0)
    check("paired T5: difference ~ +0.15 with CI excluding 0", 0.10 < pr["diff"] < 0.20 and pr["lo"] > 0, pr)
    d2 = d0.copy()
    d2.loc[d2.tier == 1, "covered"] = True
    g2 = an.gee_tiers(d2)
    check("GEE: separation in a tier handled by bootstrap", g2[1]["method"].startswith("bootstrap"), g2[1])


# ------------------------------------------------------------------ proxy leakage canary

def _all_proxies(c):
    out = {}
    for kind in ("h", "d", "s"):
        for key, v in c[kind].items():
            out[(kind,) + key] = (v["q"], v["oof"])
    return out


def _same_proxies(c0, c1, val_too=True):
    a, b = _all_proxies(c0), _all_proxies(c1)
    nv = c0["nv"]
    ok = True
    for key in a:
        q0, q1 = a[key][0], b[key][0]
        ok &= np.allclose(a[key][1], b[key][1]) and np.allclose(q0[nv:], q1[nv:]) and (not val_too or np.allclose(q0, q1))
    return bool(ok)


def t_canary():
    D = synth_v06.make_dataset(21)
    sps = repeated_group_kfold(D["U"].tss_group.values[D["E_idx"]], ev.DEV_SEED_BASE, n_rounds=1)
    sp = sps[0]
    cfg = [("sim_prior_base_model+nbq", "ridge"), ("sim_prior", "gbm")]
    c0 = ev.compute_split(0, sp, D, cfg)
    rng = np.random.default_rng(9)
    ct = D["E_idx"][np.concatenate([sp["cal"], sp["test"]])]
    # (1) corrupt cal + test expression (full and both halves)
    D2 = dict(D)
    for key in ("Y", "YA", "YB"):
        Z = D[key].copy()
        Z[ct] = rng.normal(5, 3, Z[ct].shape)
        D2[key] = Z
    c1 = ev.compute_split(0, sp, D2, cfg)
    check("canary: corrupting cal+test expression leaves every h, d, s unchanged", _same_proxies(c0, c1))
    check("canary: train-derived quantities unchanged (E_pool, OOF scores, mu, sigma, chosen K/d)",
          np.allclose(c0["E_pool"], c1["E_pool"]) and np.allclose(c0["ref_mu"], c1["ref_mu"])
          and np.allclose(c0["ref_sigma"], c1["ref_sigma"]) and c0["chosen"] == c1["chosen"]
          and all(np.allclose(c0["oof_score"][k], c1["oof_score"][k]) for k in c0["oof_score"]))
    check("canary: scores of cal+test do change (corruption took effect)",
          not np.allclose(c0["scores"][("linear_pca", 50, "full")], c1["scores"][("linear_pca", 50, "full")]))
    # (2) v0.6.1: replace cal + test cell counts by random numbers
    D3 = dict(D)
    U3 = D["U"].copy()
    U3.loc[ct, "ncells"] = rng.integers(5, 5000, len(ct))
    D3["U"] = U3
    c3 = ev.compute_split(0, sp, D3, cfg)
    check("canary (v0.6.1): random cal+test cell counts leave every h, d, s unchanged", _same_proxies(c0, c3))
    # (3) corrupting val may change chosen K/d but never the mean-model proxies or E_pool
    D4 = dict(D)
    vv = D["E_idx"][sp["val"]]
    Z = D["Y"].copy()
    Z[vv] = rng.normal(5, 3, Z[vv].shape)
    D4["Y"] = Z
    c4 = ev.compute_split(0, sp, D4, cfg)
    k = ("mean", "sim_prior", "gbm")
    nv = c0["nv"]
    check("canary B: corrupting val leaves mean-model h (cal+test) and E_pool unchanged",
          np.allclose(c0["h"][k]["q"][nv:], c4["h"][k]["q"][nv:]) and np.allclose(c0["E_pool"], c4["E_pool"]))
    # (4) v0.6.1: width matching never uses test scores
    sel = {"h": ("sim_prior_base_model+nbq", "ridge"), "d": "ridge", "s": "ridge"}
    grid = [("M1", 0.0), ("M2", 3), ("M3", "quad"), ("Mdist-N", 0.0), ("Mcomb-3", "lin2")]
    r0 = ev.apply_methods(c0, "linear_pca", 50, sel, grid=grid, mode="matched")
    c5 = dict(c0)
    c5["scores"] = dict(c0["scores"])
    nc = len(c0["cal"])
    sc = c0["scores"][("linear_pca", 50, "full")].copy()
    sc[nc:] = rng.permutation(sc[nc:]) * rng.uniform(0.2, 5, len(sc) - nc)
    c5["scores"][("linear_pca", 50, "full")] = sc
    r5 = ev.apply_methods(c5, "linear_pca", 50, sel, grid=grid, mode="matched")
    check("canary (v0.6.1): matched alpha' and thresholds invariant to test-unit scores",
          np.array_equal(r0.alpha_used.values, r5.alpha_used.values) and np.allclose(r0.thr.values, r5.thr.values),
          {"alpha_used": r0.groupby("method").alpha_used.first().to_dict()})
    check("width matching: widths monotone in alpha on the evaluated grid points", bool(r0.monotone.all()))


def e_def_independent(Y_rows, targets, genes, mu, sigma, k=50):
    """E_true by the written definition, independent of U.tcols: mean-model S1(k), excluding every readout
    column whose gene name equals the target (nothing excluded if the target is not a readout gene)."""
    s = (Y_rows - mu) / sigma
    out = np.empty(len(Y_rows))
    for i, tg in enumerate(targets):
        a = np.abs(s[i]).copy()
        a[genes == tg] = -np.inf
        top = np.argsort(-a)[:k]
        out[i] = np.sqrt((s[i, top] ** 2).mean())
    return out


def t_v062():
    """v0.6.2: learner units = E training units; AT units (incl. unpredictable ones) are model-fitting units
    and neighbours only; one E_true definition for all units."""
    D = synth_v06.make_dataset(41, unpred_at_frac=0.7)
    U = D["U"]
    unpred = D["truth"]["unpredictable_at"]
    sp = repeated_group_kfold(U.tss_group.values[D["E_idx"]], ev.DEV_SEED_BASE, n_rounds=1)[0]
    cfg = [("sim_prior_base_model+nbq", "ridge"), ("sim_prior", "gbm")]
    c = ev.compute_split(0, sp, D, cfg)
    lrn_true = D["E_idx"][sp["train"]]
    check("v0.6.2: compute_split runs with unpredictable AT units in the pool",
          bool(unpred[c["at"]].sum() > 0), {"n_unpredictable_at": int(unpred.sum()), "n_at": int(len(c["at"]))})
    lens = {len(v["oof"]) for kind in ("h", "d", "s") for v in c[kind].values()}
    check("v0.6.2: learner units are exactly the E training units (no AT unit), all OOF arrays of that length",
          np.array_equal(np.sort(c["lrn"]), np.sort(lrn_true)) and not np.isin(c["lrn"], D["AT_idx"]).any()
          and lens == {len(lrn_true)} and len(c["E_pool"]) == len(lrn_true) and
          all(len(v) == len(lrn_true) for v in c["oof_score"].values()), {"oof_lengths": sorted(lens)})
    G = D["genes"]
    e_at = e_def_independent(D["Y"][c["at"]], U.target.values[c["at"]], G, c["ref_mu"], c["ref_sigma"])
    e_l = e_def_independent(D["Y"][c["lrn"]], U.target.values[c["lrn"]], G, c["ref_mu"], c["ref_sigma"])
    check("v0.6.2: E_true of AT and learner units equals the written definition (independent recompute)",
          np.allclose(e_at, c["E_at"], rtol=1e-10, atol=1e-12) and np.allclose(e_l, c["E_pool"], rtol=1e-10, atol=1e-12),
          {"max_abs_diff_at": float(np.abs(e_at - c["E_at"]).max()), "n_at_target_not_readout": int(unpred[c["at"]].sum())})
    # AT units really are neighbours: changing AT responses changes query h (mean model, prior features)
    D2 = dict(D)
    rng = np.random.default_rng(5)
    for key in ("Y", "YA", "YB"):
        Z = D[key].copy()
        Z[D["AT_idx"]] = Z[D["AT_idx"]] + rng.normal(0, 2, Z[D["AT_idx"]].shape)
        D2[key] = Z
    c2 = ev.compute_split(0, sp, D2, cfg)
    k = ("mean", "sim_prior", "gbm")
    check("v0.6.2: AT responses enter the neighbour features (perturbing them changes query h)",
          not np.allclose(c["h"][k]["q"], c2["h"][k]["q"]))
    sel = {"h": ("sim_prior_base_model+nbq", "ridge"), "d": "ridge", "s": "ridge"}
    r = ev.apply_methods(c, "linear_pca", 50, sel, grid=[("M0", None), ("M1", 0.0), ("Mdist-N", 0.0), ("oracle", 0.0),
                                                        ("M3", "quad")], mode="matched")
    check("v0.6.2: stage 2 runs on this split (matched mode)", len(r) > 0 and r.thr.notna().all())


# ------------------------------------------------------------------ pseudobulk

def _write_h5ad(path, X, gt, umi, genes, ncf=None, legacy=True):
    with h5py.File(path, "w") as f:
        f.create_dataset("X", data=X.astype(np.float32))
        o = f.create_group("obs")
        cats, codes = np.unique(gt, return_inverse=True)
        if legacy:
            c = o.create_group("__categories")
            c.create_dataset("gene_transcript", data=np.array(cats, dtype="S"))
            o.create_dataset("gene_transcript", data=codes.astype(np.int16))
        o.create_dataset("UMI_count", data=np.asarray(umi, np.float64))
        if ncf is not None:
            o.create_dataset("num_cells_filtered", data=np.asarray(ncf))
        v = f.create_group("var")
        gc, gcode = np.unique(genes, return_inverse=True)
        vc = v.create_group("__categories")
        vc.create_dataset("gene_name", data=np.array(gc, dtype="S"))
        v.create_dataset("gene_name", data=gcode.astype(np.int32))


def t_pseudobulk():
    rng = np.random.default_rng(4)
    n, G, R = 3000, 120, 40
    rows = np.array([f"P{r:03d}_x_y" for r in range(R - 2)] + ["non-targeting_a", "non-targeting_b"])
    gt = rows[rng.integers(0, R, n)]
    lam = rng.gamma(2, 1, (R, G))
    rid = np.searchsorted(np.sort(rows), gt)
    srt = np.sort(rows)
    X = rng.poisson(lam[rid]).astype(np.float32)
    umi = X.sum(1) + rng.integers(50, 200, n)
    genes = np.array([f"gene{g:03d}" for g in range(G)])
    tmp = Path(tempfile.mkdtemp(prefix="pb_v06_"))
    sc = tmp / "sc.h5ad"
    _write_h5ad(sc, X, gt, umi, genes)
    # bulk file: per-row mean raw counts (what the publisher's raw bulk X contains)
    bmean = np.vstack([X[gt == r].mean(0) for r in srt])
    ncf = np.array([(gt == r).sum() for r in srt])
    bk = tmp / "bulk.h5ad"
    _write_h5ad(bk, bmean, srt, np.ones(R), genes, ncf)
    bad = bmean.copy()
    bad[3, 7] += 0.5
    bkbad = tmp / "bulk_bad.h5ad"
    _write_h5ad(bkbad, bad, srt, np.ones(R), genes, ncf)
    r1 = pb.build(sc, tmp / "p1", bk, passes=1, chunk_rows=700, min_free=0.0, pause_free=0.0, log=lambda *a: None)
    r2 = pb.build(sc, tmp / "p2", bk, passes=2, chunk_rows=700, min_free=0.0, pause_free=0.0, log=lambda *a: None)
    F1 = np.load(tmp / "p1/pseudobulk_logmean.npy")
    F2 = np.load(tmp / "p2/pseudobulk_logmean.npy")
    A, B = np.load(tmp / "p2/pseudobulk_logmean_A.npy"), np.load(tmp / "p2/pseudobulk_logmean_B.npy")
    rr = pd.read_csv(tmp / "p2/pseudobulk_rows.csv")
    x = np.log1p(X.astype(np.float64) * (1e4 / umi)[:, None])
    ref = np.vstack([x[gt == r].mean(0) for r in srt])
    check("pseudobulk: 1-pass equals direct computation", np.allclose(F1, ref, atol=1e-5))
    check("pseudobulk: 2-pass (low memory) identical to 1-pass", np.array_equal(F1, F2))
    comb = (A * rr.n_A.values[:, None] + B * rr.n_B.values[:, None]) / rr.n_cells.values[:, None]
    check("pseudobulk: halves recombine to the full mean", np.allclose(comb, F2, atol=1e-5))
    check("pseudobulk: halves balanced (|n_A - n_B| <= 1) and deterministic",
          (np.abs(rr.n_A - rr.n_B) <= 1).all() and np.array_equal(
              pd.read_csv(tmp / "p1/pseudobulk_rows.csv").n_A.values, rr.n_A.values))
    check("pseudobulk QC: consistent bulk file passes", r1["qc_vs_bulk"]["pass"] and r2["qc_vs_bulk"]["pass"],
          r2["qc_vs_bulk"])
    r3 = pb.build(sc, tmp / "p3", bkbad, passes=2, chunk_rows=700, min_free=0.0, pause_free=0.0, log=lambda *a: None)
    check("pseudobulk QC: a single perturbed bulk entry is detected", (not r3["qc_vs_bulk"]["pass"])
          and r3["qc_vs_bulk"]["n_entries_outside_tol"] == 1, r3["qc_vs_bulk"])
    r4 = pb.build(sc, tmp / "p4", bk, passes=2, min_free=1e6, log=lambda *a: None)
    check("pseudobulk: stops before reading when available memory < threshold",
          r4["status"] == "stopped_low_memory" and not (tmp / "p4/pseudobulk_logmean.npy").exists())
    check("pseudobulk: no per-row QC means saved", not any("raw" in p.name.lower() for p in (tmp / "p2").iterdir()),
          sorted(p.name for p in (tmp / "p2").iterdir()))


# ------------------------------------------------------------------ end-to-end

SCENARIOS = {"pos": dict(seed=31, b_mod=0.8, b_base=0.3, sd_e=0.4),
             "neg": dict(seed=32, b_mod=0.0, b_base=0.0, sd_e=0.6),
             "dist": dict(seed=33, b_mod=0.0, b_base=0.0, sd_e=0.6, orphan_frac=0.3)}


def expected_11d(kind, summ, pt):
    """Protocol v0.6.1 section 11d expectations (code checks on synthetic data)."""
    cm = summ["criteria_matched"]
    comp = summ["matched_paired_T5"]
    prim, bd = summ["selection"]["primary"][0], summ["selection"]["best_dist"][0]
    inc = summ["incremental"]["linear_pca"]["b2_miscover"]
    out = {}
    if kind == "pos":
        out["P+"] = cm["P"] == "P+"
        out["matched T5 primary - M0: CI lower > 0"] = comp[f"{prim}_vs_M0"]["lo"] > 0
        out["(b2) h coefficient CI lower > 1 (OR)"] = inc is not None and inc["zh"]["lo"] > 1
        out["matched T5 Mcomb-N - best M_dist > 0"] = comp[f"Mcomb-N_vs_{bd}"]["diff"] > 0
    elif kind == "neg":
        out["P-"] = cm["P"] == "P-"
        p = comp[f"{prim}_vs_M0"]
        out["matched T5 primary - M0: CI includes 0 or diff < 2 pp"] = (p["lo"] <= 0 <= p["hi"]) or p["diff"] < 0.02
        ps = cm["perm_summary"]
        out["M_perm: 0 of 20 meet W+ (matched)"] = ps is not None and ps["n_meeting_W+"] == 0
    else:
        x = pt[(pt.model == "linear_pca") & (pt.k == 50)]
        dln, sl = summ["selection"]["d"], summ["selection"]["h"]
        rd = an._gene_rho(x[(x.kind == "d") & (x.proxy == dln)], "value", "score")
        hrows = pt[(pt.model == "linear_pca") & (pt.kind == "h") & (pt.proxy == sl)]
        sc = x[(x.kind == "d") & (x.proxy == dln)][["split", "target", "score"]]
        rh = an._gene_rho(hrows.merge(sc, on=["split", "target"], suffixes=("_h", "")), "value", "score")
        out["rho(d, S1) > rho(h, S1)"] = rd > rh
        out["_rho_d_S1"], out["_rho_h_S1"] = rd, rh
        s0 = summ["dgo_strata"]["matched|M0"]["max_dev"]
        sd = summ["dgo_strata"]["matched|Mdist-N"]["max_dev"]
        out["matched D_GO-quintile max dev: Mdist-N <= 0.5 x M0"] = sd <= 0.5 * s0
        out["_dgo_maxdev_M0_MdistN"] = [s0, sd]
        out["(b2) h coefficient CI includes 1 (OR)"] = inc is not None and inc["zh"]["lo"] <= 1 <= inc["zh"]["hi"]
        c = comp[f"Mcomb-N_vs_{bd}"]
        out["matched T5 Mcomb-N - best M_dist: CI includes 0"] = c["lo"] <= 0 <= c["hi"]
    return out


def e2e(kind, outdir):
    D = synth_v06.make_dataset(**SCENARIOS[kind])
    t0 = time.time()
    summ, cov, pt = p6.run_dataset(D, ev.DEV_SEED_BASE, outdir, n_rounds=3)
    exp = expected_11d(kind, summ, pt)
    cm, cn = summ["criteria_matched"], summ["criteria_nominal"]
    out = {"kind": kind, "seconds": round(time.time() - t0), "selection": summ["selection"],
           "expected_11d": exp, "all_expected_pass": bool(all(v for k, v in exp.items() if not k.startswith("_"))),
           "matched": {k: cm.get(k) for k in ("B0", "P", "W", "W_perm_mean", "perm_summary", "paired_T5", "W_detail", "overall")},
           "nominal": {k: cn.get(k) for k in ("B0", "P", "W", "W_perm_mean", "perm_summary", "paired_T5", "overall")},
           "predictability": cm.get("predictability"), "matched_paired_T5": summ["matched_paired_T5"],
           "incremental_linear_pca": summ["incremental"]["linear_pca"], "monotone_all": summ["monotone_all"],
           "confounding_linear_pca_logE": {"unadj": summ["confounding"].get("linear_pca", {}).get("unadjusted"),
                                           "adj": summ["confounding"].get("linear_pca", {}).get("adjusted", {}).get("zlogE")},
           "tier_stability": summ["tier_stability_frac_genes_changing"],
           "genes_tested_3x": bool((cov[(cov.method == "M0") & (cov.model == "linear_pca") & (cov.k == 50)
                                        & (cov.variant == "full") & (cov["mode"] == "nominal")].groupby("target").size() == 3).all())}
    Path(outdir, "e2e_result.json").write_text(json.dumps(out, indent=1, default=p6._js), encoding="utf8")
    print(json.dumps({"kind": kind, "expected_11d": exp, "all_expected_pass": out["all_expected_pass"]}, default=p6._js),
          flush=True)


def t_v07():
    """v0.7: randomised width matching in the confirmatory module, E2 reporting, pre-specified wording, power,
    frozen run path. Synthetic data only (D0' e2e_pos caches)."""
    import pickle
    import confirm_v07 as cf
    import sens_v07 as sv
    root = Path(__file__).resolve().parents[1]
    src = root / "results" / "d0p_v0.6.2" / "e2e_pos"
    c = pickle.load(open(src / "cache" / "split_00.pkl", "rb"))
    sl = json.loads((src / "summary.json").read_text(encoding="utf8"))["selection"]
    sel = {"h": tuple(sl["h"].split("|")), "d": sl["d"], "s": sl["s"]}
    grid = [(m, sv._hv(m, h)) for m, h in sl["hypers"].items() if m in ("M1", "M2", "Mdist-N", "Mdist-M", "Mscore-N")]
    a = sv.rand_match_rows(c, "linear_pca", 50, sel, grid).sort_values(["method", "target"]).reset_index(drop=True)
    b = cf.rand_rows(c, "linear_pca", 50, sel, grid).sort_values(["method", "target"]).reset_index(drop=True)
    check("v0.7: confirm rand_rows identical to sens_v07 randomised matching",
          len(a) == len(b) and bool((a.thr.values == b.thr.values).all() and (a.covered.values == b.covered.values).all()))
    m = b[b.match_status == "matched"].groupby("method")[["W0", "expected_width"]].first()
    check("v0.7: expected width equals M0 width for matched methods",
          len(m) > 0 and bool(((m.expected_width - m.W0).abs() / m.W0).max() < 1e-9), {"n_methods": int(len(m))})
    p1 = cf.rand_rows(c, "linear_pca", 50, sel, grid[:1], perm_r=3)
    p2 = cf.rand_rows(c, "linear_pca", 50, sel, grid[:1], perm_r=3)
    check("v0.7: permuted rows reproducible; decision subset uses the real h",
          bool((p1.thr.values == p2.thr.values).all() and np.allclose(
              p1[p1.method == "M0"].sort_values("target").h.values,
              b[b.method == "M0"].sort_values("target").h.values)))
    d0 = b[b.method == "M0"].sort_values(["split", "target"])
    dm = b[b.method == grid[0][0]].sort_values(["split", "target"])
    r = cf.e2(dm, d0)
    check("v0.7: E2 always carries the subset width ratio",
          "width_ratio_subset" in r and "subset width ratio" in r["report"], r["report"])
    E = lambda d, lo, hi: {"diff": d, "lo": lo, "hi": hi}  # noqa: E731
    g_ok = {"nominal_marginal_in_band": True, "M_perm": True}
    s1p = {"pass": True}
    w = {"b0": cf.classify(s1p, {"pass": False}, False, g_ok),
         "s1fail": cf.classify({"pass": False}, {"pass": False}, True, g_ok),
         "both": cf.classify(s1p, {"pass": True}, True, g_ok),
         "lowpow": cf.classify(s1p, {"pass": False, "E1": E(0.02, -0.01, 0.05), "E2": E(0.03, 0.001, 0.06)}, True, g_ok),
         "small": cf.classify(s1p, {"pass": False, "E1": E(0.02, 0.005, 0.035), "E2": E(0.025, 0.01, 0.04)}, True, g_ok),
         "guard": cf.classify(s1p, {"pass": True}, True, {"nominal_marginal_in_band": False, "M_perm": True})}
    check("v0.7: pre-specified wording (limited power; never 'no increment'; guard downgrade)",
          w["b0"].startswith("B0 not met") and w["s1fail"].startswith("step 1 not passed")
          and w["both"].startswith("steps 1 and 2 passed") and "limited power" in w["lowpow"]
          and "below the pre-specified 3 pp margin" in w["small"] and "downgraded" in w["guard"]
          and not any("no increment" in x for x in w.values()), w)
    pw = [cf.power_at_se(x, 0.016, 0.03) for x in (0.0, 0.03, 0.06, 0.09)]
    check("v0.7: power increases with effect and is ~alpha-free at 0", bool(np.all(np.diff(pw) > 0) and pw[0] < 0.05),
          [round(x, 3) for x in pw])
    # frozen path: no selection function may be called
    orig = (an.select_hyper, an.select_primary, an.select_proxy, an.select_learner)

    def boom(*a_, **k_):
        raise AssertionError("selection called in frozen mode")
    an.select_hyper = an.select_primary = an.select_proxy = an.select_learner = boom
    try:
        with tempfile.TemporaryDirectory() as td:
            import shutil
            (Path(td) / "cache").mkdir()
            for i in range(2):
                shutil.copy2(src / "cache" / f"split_{i:02d}.pkl", Path(td) / "cache" / f"split_{i:02d}.pkl")
            Dm = {"U": None}
            p6.stage1 = (lambda *a_, **k_: None)
            fz = {"h": sl["h"], "d": sl["d"], "s": sl["s"], "hypers": dict(sl["hypers"]), "primary": sl["primary"][0],
                  "best_dist": sl["best_dist"][0]}
            summ, _, _ = p6.run_dataset(Dm, ev.CONF_SEED_BASE, td, n_rounds=3, frozen=fz, log=lambda *_: None)
            ok = summ["selection"]["primary"][0] == fz["primary"]
    except AssertionError as e:
        ok = False
        print(e)
    finally:
        an.select_hyper, an.select_primary, an.select_proxy, an.select_learner = orig
        import importlib
        importlib.reload(p6)
    check("v0.7: run_dataset(frozen) runs without any selection rule", ok)


def t_v07_data():
    """v0.7: dataset-generic builder/loader reproduce the frozen K562 dev objects exactly (dev data, no
    confirmation data)."""
    import data_v07 as dv
    from paths import PROCESSED, raw
    with tempfile.TemporaryDirectory() as td:
        U, _, info = dv.build_units(raw("K562_essential_raw_bulk_01.h5ad"), Path(td) / "coords.csv")
        no_query = not (Path(td) / "coords.csv").exists()
    A = pd.read_csv(PROCESSED / "v0.3_figshare_k562ess" / "p0v3_units_all_targets.csv").sort_values("target")
    B = U.sort_values("target")
    cols = ["gene_transcript", "ncells", "E", "always_train", "dropped_tss_conflict", "tss_group", "chr", "start"]
    check("v0.7 data: unit table identical to the frozen K562 table (no Ensembl query)",
          no_query and all((A[c].astype(str).values == B[c].astype(str).values).all() for c in cols), info)
    Dg = dv.load_dataset(U, PROCESSED / "v0.6" / "k562ess", PROCESSED / "priors" / "string_expdb_v0.4_universe.parquet")
    Dd = p6.load_dev_v06(PROCESSED / "v0.6" / "k562ess")
    same = all(np.array_equal(Dg[k], Dd[k]) for k in ("Y", "YA", "YB", "y_ctrl", "E_idx", "AT_idx")) and \
        all(np.array_equal(Dg["pri"][k], Dd["pri"][k]) for k in ("J", "Jf", "S", "ase", "n_go")) and \
        bool((Dg["U"].target.values == Dd["U"].target.values).all())
    check("v0.7 data: load_dataset identical to load_dev_v06", same)


def e2e_v07(kind, outdir):
    """v0.7 confirmation-path e2e on synthetic data: selection frozen from the D0' (v0.6.2) e2e run of the same
    scenario; a NEW synthetic dataset (seed + 1000) plays the confirmation set; stage 1 -> ext -> frozen
    run_dataset -> confirm_v07.analyze, exactly as confirm_run_v07 does."""
    import confirm_v07 as cf
    import ext_v07dev as ex
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    sl = json.loads((root / "results" / "d0p_v0.6.2" / f"e2e_{kind}" / "summary.json").read_text(encoding="utf8"))["selection"]
    fz = {"h": sl["h"], "d": sl["d"], "s": sl["s"],
          "hypers": {**sl["hypers"], "MD-GO": "5", "MD-STR": "5", "Munc-N": "1.0", "Munc-M": "3"},
          "primary": sl["primary"][0], "best_dist": sl["best_dist"][0]}
    sc = dict(SCENARIOS[kind])
    sc["seed"] = sc["seed"] + 1000
    D = synth_v06.make_dataset(**sc)
    t0 = time.time()
    p6.stage1(D, ev.CONF_SEED_BASE, out / "cache", 3, print)
    (out / "ext_cache").mkdir(exist_ok=True)
    for c in p6.load_caches(out / "cache"):
        f = out / "ext_cache" / f"ext_split_{c['split_id']:02d}.npz"
        if not f.exists():
            np.savez(f, **ex.compute_ext(c, D, log=lambda *_: None))
    summ, cov, _ = p6.run_dataset(D, ev.CONF_SEED_BASE, out, n_rounds=3, frozen=fz, ext_dir=out / "ext_cache")
    caches = p6.load_caches(out / "cache")
    for c in caches:
        z = np.load(out / "ext_cache" / f"ext_split_{c['split_id']:02d}.npz")
        ex.attach(c, {k: z[k] for k in z.files})
    res, _ = cf.analyze(caches, fz, cov_grid=cov)
    cf.write(res, out / "confirm_v07_result.json")
    s1 = res["step1"]
    if kind == "pos":
        exp = {"B0": res["B0"]["pass"], "step1_pass": s1["pass"], "guard_M_perm": res["guards"]["M_perm"],
               "classification_step1_passed": res["classification"].startswith(("step 1 passed", "steps 1 and 2 passed"))}
    else:
        exp = {"step1_not_pass": not s1["pass"], "guard_M_perm": res["guards"]["M_perm"],
               "classification_not_success": not res["classification"].startswith("steps 1 and 2 passed")}
    exp["E2_reports_width_ratio"] = "subset width ratio" in s1["E2"]["report"]
    exp["no_selection_in_confirmation"] = summ["selection"]["primary"][0] == fz["primary"]
    o = {"kind": kind, "seconds": round(time.time() - t0), "frozen": fz, "expected_v07": exp,
         "all_expected_pass": bool(all(exp.values())), "classification": res["classification"],
         "step1": {e: {k: s1[e].get(k) for k in ("diff", "lo", "hi", "report")} for e in ("E1", "E2")},
         "step2": {e: {k: res["step2"][e].get(k) for k in ("diff", "lo", "hi", "report")} for e in ("E1", "E2")},
         "guards": res["guards"], "B0": res["B0"], "matching": res["matching"]}
    (out / "e2e_v07_result.json").write_text(json.dumps(o, indent=1, default=cf._js), encoding="utf8")
    print(json.dumps({"kind": kind, "expected_v07": exp, "all_expected_pass": o["all_expected_pass"]}, default=cf._js),
          flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "unit":
        for fn in (t_splits, t_methods, t_analysis, t_canary, t_v062, t_v07, t_v07_data, t_pseudobulk):
            fn()
        out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
        out.mkdir(parents=True, exist_ok=True)
        (out / "unit_results.json").write_text(json.dumps(RES, indent=1, default=p6._js), encoding="utf8")
        print(f"{sum(v['pass'] for v in RES.values())}/{len(RES)} passed")
    elif sys.argv[1] == "e2e_v07":
        e2e_v07(sys.argv[2], sys.argv[3])
    else:
        e2e(sys.argv[2], sys.argv[3])
