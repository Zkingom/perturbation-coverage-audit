# confirm_v07.py -- protocol v0.7 section 4: confirmatory endpoints under randomised width matching.
# Self-contained given the per-split stage-1 caches (optionally carrying ext arrays) and the frozen selection.
#   E1  paired T5 coverage difference (method - comparator) under randomised width matching
#   E2  decision scenario: per split the ceil(10%) test units with the largest h (same units in both arms);
#       paired coverage difference, ALWAYS reported together with the subset width ratio (method / comparator)
#   Step 1: M1 - M0; Step 2 (only if step 1 passes): M1 - best distance / difficulty baseline (frozen).
#   Pass at a step = for both E1 and E2: 95% CI lower bound > 0 AND point estimate >= margin (5 pp, 3 pp).
#   Guards: B0, P (reported), nominal marginal coverage of M1 in [0.88, 0.92], M_perm (20 permutations of h
#   under randomised matching: mean E1 < 2 pp and <= 1 permutation individually meeting the step-1 E1 rule),
#   matching status.
# Randomised matching is the v0.7 primary width matching (v0.7_draft_notes item 1). For perm_r = None the
# computation is identical to sens_v07.rand_match_rows (checked in tests_v06.t_v07).
import json
import math

import numpy as np
import pandas as pd
from scipy.stats import norm

import analysis_v06 as an
import eval_v06 as ev
from methods_v06 import ALPHA_GRID, METHODS, capped_width, permute_within, thresholds
from paths import MASTER_SEED
from sens_v07 import MIX_SEED, _hv

MARGIN_STEP1, MARGIN_STEP2 = 0.05, 0.03
PERM_MEAN_MAX, PERM_MAX_PASSING = 0.02, 1
NOMINAL_BAND = (0.88, 0.92)
DEC_FRAC = 0.10
# dev-set point estimates under randomised width matching (results/d1_v0.7dev/confirm_v07_on_dev/confirm_v07_dev.json)
DEV_ESTIMATES = {"E1_step1": 0.120000, "E2_step1": 0.153846, "E1_step2": 0.063256, "E2_step2": 0.064103}


def rand_rows(cache, model, k, sel, grid_methods, perm_r=None):
    """Randomised width matching rows (mode 'matched_rand'). perm_r: None, or 0..N_PERM-1 for permuted h."""
    nc = len(cache["cal"])
    s = cache["scores"][(model, k, "full")]
    s_cal, s_test = s[:nc], s[nc:]
    P_cal, P_test, P_oof, _ = ev.proxies(cache, model, k, sel)
    h_rank = np.asarray(P_test["h"]).copy()          # decision subset always uses the real h
    if perm_r is not None:
        P_cal["h"], P_test["h"] = permute_within(P_cal["h"], P_test["h"],
                                                 ev._seed(ev.PERM_SEED_BASE + perm_r, cache["split_id"]))
    wi = ev.MODELS.index(model)
    u_test = np.random.default_rng(ev._seed(MASTER_SEED, 5000, cache["split_id"], wi, k, 0,
                                            0 if perm_r is None else perm_r + 1)).uniform(size=len(s_test))
    vseed = ev._seed(MIX_SEED, cache["split_id"], wi, k) if perm_r is None else \
        ev._seed(MIX_SEED, cache["split_id"], wi, k, perm_r + 1)
    v = np.random.default_rng(vseed).uniform(size=len(s_test))
    t0 = thresholds("M0", None, s_cal, P_cal, P_test, P_oof, u_test, ev.ALPHA)
    W0 = capped_width(t0, s_cal)
    cap = 2.0 * float(np.max(s_cal))
    rows = []
    avail = lambda nm: METHODS[nm][1] is None or isinstance(METHODS[nm][1], tuple) or METHODS[nm][1] in P_cal  # noqa: E731
    for name, hy in [("M0", None)] + [g_ for g_ in grid_methods if g_[0] != "M0" and avail(g_[0])]:
        memo = {}

        def W(i):
            if i not in memo:
                t = thresholds(name, _hv(name, hy), s_cal, P_cal, P_test, P_oof, u_test, ALPHA_GRID[i])
                memo[i] = (capped_width(t, s_cal), t)
            return memo[i][0]
        exp_w = W0
        if name == "M0":
            t, a_lo, a_hi, pi, status = t0, ev.ALPHA, ev.ALPHA, 1.0, "reference"
        else:
            n = len(ALPHA_GRID)
            if W(0) < W0:
                t, a_lo, a_hi, pi, status = memo[0][1], ALPHA_GRID[0], ALPHA_GRID[0], 1.0, "unmatched_too_narrow"
                exp_w = W(0)
            elif W(n - 1) > W0:
                t, a_lo, a_hi, pi, status = memo[n - 1][1], ALPHA_GRID[-1], ALPHA_GRID[-1], 1.0, "unmatched_too_wide"
                exp_w = W(n - 1)
            else:
                lo, hi = 0, n - 1
                if W(hi) == W0:
                    lo = hi
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    if W(mid) >= W0:
                        lo = mid
                    else:
                        hi = mid
                wl, wh = W(lo), W(hi)
                pi = 1.0 if wl == wh else (W0 - wh) / (wl - wh)
                t = np.where(v < pi, memo[lo][1], memo[hi][1])
                a_lo, a_hi = ALPHA_GRID[lo], ALPHA_GRID[hi]
                status = "matched" if wl >= W0 >= wh else "bracket_nonmonotone"
                exp_w = pi * wl + (1 - pi) * wh
        rows.append(pd.DataFrame({
            "split": cache["split_id"], "model": model, "k": k, "variant": "full", "mode": "matched_rand",
            "method": name, "hyper": str(hy), "perm": -1 if perm_r is None else perm_r, "alpha_lo": a_lo,
            "alpha_hi": a_hi, "pi": pi, "match_status": status, "W0": W0, "expected_width": exp_w,
            "target": cache["target_calt"][nc:], "tss_group": cache["tss_calt"][nc:],
            "tier": cache["tier_calt"][nc:], "tier3": cache["tier3_calt"][nc:], "E_true": cache["E_calt"][nc:],
            "h": h_rank, "d": P_test["d"], "score": s_test, "thr": t,
            "w": np.clip(np.where(np.isposinf(t), cap, t), 0.0, cap), "covered": s_test <= t}))
    return pd.concat(rows, ignore_index=True)


def e2(df_m, df_0):
    r = an.decision_paired(df_m, df_0, frac=DEC_FRAC)
    r["width_ratio_subset"] = r["width_method"] / r["width_M0"]
    r["report"] = (f"{100 * r['diff']:+.1f} pp [{100 * r['lo']:.1f}, {100 * r['hi']:.1f}]; "
                   f"subset width ratio {r['width_ratio_subset']:.2f}")
    return r


def _ok(r, margin):
    return bool(np.isfinite(r["lo"]) and r["lo"] > 0 and r["diff"] >= margin)


def power_at_se(est, se, margin):
    crit = max(margin, 1.96 * se)
    return float(1 - norm.cdf((crit - est) / se)) if se > 0 else float("nan")


def classify(step1, step2, b0, guards):
    """Pre-specified wording (v0.7 section 4). Never states 'no increment'."""
    if not b0:
        return "B0 not met: the failure was not reproduced on this dataset; steps 1-2 are not interpreted"
    if not step1["pass"]:
        return "step 1 not passed: repair not confirmed"
    if step2["pass"]:
        s = ("steps 1 and 2 passed: the deployable prior-neighbour proxy repairs part of the failure and adds to "
             "the strongest distance / difficulty baseline")
    else:
        e = [step2["E1"], step2["E2"]]
        if all(x["diff"] > 0 for x in e) and any(x["lo"] <= 0 <= x["hi"] for x in e):
            s = ("step 1 passed; step 2 not passed: increment over the distance / difficulty baseline not confirmed "
                 "(limited power; see power at the observed SE)")
        elif all(x["lo"] > 0 for x in e):
            s = ("step 1 passed; step 2 not passed: statistical evidence of an increment, but below the "
                 "pre-specified 3 pp margin")
        else:
            s = "step 1 passed; step 2 not passed: increment over the distance / difficulty baseline not confirmed"
    bad = [k for k, v in guards.items() if v is False]
    if bad:
        s += " [downgraded: guard(s) not met: " + ", ".join(bad) + "]"
    return s


def analyze(caches, frozen, log=print, model=ev.PRIMARY_MODEL, k=ev.PRIMARY_K, n_perm=ev.N_PERM, cov_grid=None):
    sel = {"h": tuple(frozen["h"].split("|")), "d": frozen["d"], "s": frozen["s"]}
    hyp = frozen["hypers"]
    prim, best = frozen["primary"], frozen["best_dist"]
    grid = [(m, _hv(m, h)) for m, h in hyp.items()]
    rm = pd.concat([rand_rows(c, model, k, sel, grid) for c in caches], ignore_index=True)
    log(f"confirm_v07: randomised matching rows {len(rm)}")
    pick = lambda d, m: d[(d.method == m)].sort_values(["split", "target"])  # noqa: E731
    d0 = pick(rm, "M0")
    out = {"frozen": frozen, "n_splits": len(caches), "n_test_rows": int(len(d0)),
           "n_genes": int(d0.target.nunique()), "n_T5_rows": int((d0.tier == 5).sum())}
    # matching quality
    mq = {}
    for m in rm.method.unique():
        dm = pick(rm, m)
        sp = dm.groupby("split").agg(W0=("W0", "first"), ew=("expected_width", "first"), w=("w", "mean"))
        mq[m] = {"status": dm.groupby("split").match_status.first().value_counts().to_dict(),
                 "max_rel_dev_expected_width": float(((sp.ew - sp.W0).abs() / sp.W0).max()),
                 "realized_width_ratio_range": [float((sp.w / sp.W0).min()), float((sp.w / sp.W0).max())]}
    out["matching"] = mq
    out["tier_tables_rand"] = {m: an.tier_table(pick(rm, m)) for m in rm.method.unique()}
    # B0 and P
    g0 = an.gee_tiers(d0)
    b0 = bool(g0[5]["cov"] < 0.85 and g0[5]["hi"] < 0.90)
    out["B0"] = {"pass": b0, "M0_T5": g0[5]}
    pt = pd.concat([ev.proxy_table(c) for c in caches], ignore_index=True)
    pth = pt[(pt.model == model) & (pt.kind == "h") & (pt.proxy == frozen["h"])].rename(columns={"value": "h"})
    out["P"] = an.predictability(pth)
    (rlo, rhi), (alo, ahi) = out["P"]["rho_ci"], out["P"]["auroc_ci"]
    out["P"]["class"] = "P+" if (rlo >= 0.30 and alo >= 0.65) else ("P-" if (rhi < 0.20 or ahi < 0.60) else "P0")
    # steps
    dm, db = pick(rm, prim), pick(rm, best)
    s1 = {"E1": an.paired_T5(dm, d0), "E2": e2(dm, d0)}
    s1["pass"] = _ok(s1["E1"], MARGIN_STEP1) and _ok(s1["E2"], MARGIN_STEP1)
    s2 = {"E1": an.paired_T5(dm, db), "E2": e2(dm, db), "tested": bool(s1["pass"])}
    s2["pass"] = bool(s1["pass"] and _ok(s2["E1"], MARGIN_STEP2) and _ok(s2["E2"], MARGIN_STEP2))
    for st, nm in ((s1, "step1"), (s2, "step2")):
        mg = MARGIN_STEP1 if nm == "step1" else MARGIN_STEP2
        for e in ("E1", "E2"):
            se = (st[e]["hi"] - st[e]["lo"]) / 3.92
            dev = DEV_ESTIMATES[f"{e}_{nm}"]
            st[e]["se"] = se
            st[e]["power_at_observed_se"] = {f"{f:.2f}x_dev": power_at_se(f * dev, se, mg) for f in (1.0, 0.75, 0.5)}
    out["step1"], out["step2"] = s1, s2
    # guards
    nom = pd.concat([ev.apply_methods(c, model, k, sel, grid=[(prim, _hv(prim, hyp[prim]))], mode="nominal")
                     for c in caches], ignore_index=True)
    nom_marg = float(nom.covered.mean())
    perm = []
    for r in range(n_perm):
        rp = pd.concat([rand_rows(c, model, k, sel, [(prim, _hv(prim, hyp[prim]))], perm_r=r) for c in caches],
                       ignore_index=True)
        x = an.paired_T5(pick(rp, prim), pick(rp[rp.method == "M0"], "M0"))
        perm.append({"perm": r, **{kk: x[kk] for kk in ("diff", "lo", "hi")}, "meets_step1_E1": _ok(x, MARGIN_STEP1)})
    pdf = pd.DataFrame(perm)
    out["M_perm"] = {"rows": perm, "mean_E1": float(pdf["diff"].mean()), "n_meeting_step1_E1": int(pdf.meets_step1_E1.sum())}
    guards = {"nominal_marginal_in_band": bool(NOMINAL_BAND[0] <= nom_marg <= NOMINAL_BAND[1]),
              "M_perm": bool(pdf["diff"].mean() < PERM_MEAN_MAX and pdf.meets_step1_E1.sum() <= PERM_MAX_PASSING),
              "matching_primary_comparators": bool(all(set(mq[m]["status"]) <= {"matched", "reference"}
                                                       for m in ("M0", prim, best) if m in mq)),
              "P_not_minus": out["P"]["class"] != "P-"}
    out["guards"] = guards
    out["nominal_marginal_primary"] = nom_marg
    out["classification"] = classify(s1, s2, b0, guards)
    # descriptive: secondary method and every baseline under randomised matching
    desc = {}
    for m in [x for x in hyp if x not in (prim,)]:
        dx = pick(rm, m)
        if not len(dx):
            continue
        desc[m] = {"E1_vs_M0": an.paired_T5(dx, d0), "E2_vs_M0": e2(dx, d0)}
        if m != best:
            desc[m][f"E1_vs_{best}"] = an.paired_T5(dx, db)
    out["descriptive_rand"] = desc
    if cov_grid is not None:                 # sensitivity: v0.6.2 grid matching
        g = cov_grid[(cov_grid.model == model) & (cov_grid.k == k) & (cov_grid.variant == "full")
                     & (cov_grid["mode"] == "matched") & (cov_grid.perm == -1)]
        gp = lambda m: g[(g.method == m) & (g.hyper == ("None" if m == "M0" else str(hyp[m])))].sort_values(  # noqa: E731
            ["split", "target"])
        out["sensitivity_grid_matching"] = {
            "step1": {"E1": an.paired_T5(gp(prim), gp("M0")), "E2": e2(gp(prim), gp("M0"))},
            "step2": {"E1": an.paired_T5(gp(prim), gp(best)), "E2": e2(gp(prim), gp(best))}}
    return out, rm


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def write(out, path):
    path.write_text(json.dumps(out, indent=1, default=_js), encoding="utf8")
