# sens_v07.py -- analyses added for the D1 report under protocol v0.6.2 (they do not change any v0.6.2 result):
#   (1) v0.7-style randomised width matching (sensitivity; reports/protocol/v0.7_draft_notes.md item 1):
#       for each method at its v0.6.2-chosen hyperparameter, the two adjacent grid alphas whose capped mean
#       widths bracket M0's width W0 are mixed per test unit with probability pi so that the expected mean
#       width equals W0. Never uses test-unit scores. M3-type methods use the same fixed u as eval_v06.
#   (2) W status with criterion (c) judged at nominal alpha and (a), (b), (d) under (1) (draft note item 2).
#   (3) Confounding: logistic coefficients with sign + gene-level (TSS-group) cluster-bootstrap percentile CIs;
#       Firth-penalised logistic as sensitivity (draft note item 3).
#   (4) User-required D1 item: M_dist (Mdist-N, Mdist-M) vs M0 coverage by D_GO and D_STRING strata.
import json
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm

import analysis_v06 as an
import eval_v06 as ev
from methods_v06 import ALPHA_GRID, METHODS, capped_width, thresholds
from paths import MASTER_SEED

MIX_SEED = MASTER_SEED + 8000
BOOT_SEED = MASTER_SEED + 9000
N_BOOT_CONF = 500
N_BOOT_FIRTH = 200
SEP_ABS_COEF = 15.0


def _hv(meth, hy):
    from methods_v06 import METHODS
    fam = METHODS[meth][0]
    if fam == "mondrian":
        return int(hy)
    if fam in ("gcc", "gcc2"):
        return hy
    if hy in (None, "None"):
        return None
    return float(hy)


# ------------------------------------------------------------------ (1) randomised width matching

def rand_match_rows(cache, model, k, sel, grid_methods):
    nc = len(cache["cal"])
    s = cache["scores"][(model, k, "full")]
    s_cal, s_test = s[:nc], s[nc:]
    P_cal, P_test, P_oof, _ = ev.proxies(cache, model, k, sel)
    wi = ev.MODELS.index(model)
    u_test = np.random.default_rng(ev._seed(MASTER_SEED, 5000, cache["split_id"], wi, k, 0, 0)).uniform(size=len(s_test))
    v = np.random.default_rng(ev._seed(MIX_SEED, cache["split_id"], wi, k)).uniform(size=len(s_test))
    t0 = thresholds("M0", None, s_cal, P_cal, P_test, P_oof, u_test, ev.ALPHA)
    W0 = capped_width(t0, s_cal)
    cap = 2.0 * float(np.max(s_cal))
    rows = []
    avail = lambda nm: METHODS[nm][1] is None or isinstance(METHODS[nm][1], tuple) or METHODS[nm][1] in P_cal  # noqa: E731
    for name, hy in [("M0", None)] + [g_ for g_ in grid_methods if avail(g_[0])]:
        memo = {}

        def W(i):
            if i not in memo:
                t = thresholds(name, _hv(name, hy), s_cal, P_cal, P_test, P_oof, u_test, ALPHA_GRID[i])
                memo[i] = (capped_width(t, s_cal), t)
            return memo[i][0]
        if name == "M0":
            t, a_lo, a_hi, pi, status = t0, ev.ALPHA, ev.ALPHA, 1.0, "reference"
        else:
            n = len(ALPHA_GRID)
            if W(0) < W0:
                t, a_lo, a_hi, pi, status = memo[0][1], ALPHA_GRID[0], ALPHA_GRID[0], 1.0, "unmatched_too_narrow"
            elif W(n - 1) > W0:
                t, a_lo, a_hi, pi, status = memo[n - 1][1], ALPHA_GRID[-1], ALPHA_GRID[-1], 1.0, "unmatched_too_wide"
            else:
                lo, hi = 0, n - 1                       # invariant W(lo) >= W0 > W(hi)
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
                a_lo, a_hi, status = ALPHA_GRID[lo], ALPHA_GRID[hi], ("matched" if wl >= W0 >= wh else "bracket_nonmonotone")
        rows.append(pd.DataFrame({
            "split": cache["split_id"], "model": model, "k": k, "variant": "full", "mode": "matched_rand",
            "method": name, "hyper": str(hy), "perm": -1, "alpha_lo": a_lo, "alpha_hi": a_hi, "pi": pi,
            "match_status": status, "W0": W0, "target": cache["target_calt"][nc:], "tss_group": cache["tss_calt"][nc:],
            "tier": cache["tier_calt"][nc:], "tier3": cache["tier3_calt"][nc:], "E_true": cache["E_calt"][nc:],
            "h": P_test["h"], "d": P_test["d"], "score": s_test, "thr": t,
            "w": np.clip(np.where(np.isposinf(t), cap, t), 0.0, cap), "covered": s_test <= t}))
    return pd.concat(rows, ignore_index=True)


def rand_match_summary(rm, cov_nominal, summ, model):
    """Tier tables, paired T5 vs M0 and vs the best M_dist, and W status with (c) at nominal alpha."""
    sel = summ["selection"]
    hyp = sel["hypers"]
    prim, bd = sel["primary"][0], sel["best_dist"][0]
    pick = lambda d, m: d[(d.method == m) & (d.hyper == ("None" if m == "M0" else str(hyp[m])))].sort_values(  # noqa: E731
        ["split", "target"])
    d0 = pick(rm, "M0")
    out = {"model": model, "W0_mean": float(rm.W0.mean()), "methods": {}}
    for m in [x for x in an.TUNED if x in hyp]:
        dm = pick(rm, m)
        if not len(dm):
            continue
        tt = an.tier_table(dm)
        realized = dm.groupby("split").w.mean() / d0.groupby("split").w.mean()
        out["methods"][m] = {"hyper": hyp[m], "tier_table": tt, "match_status": dm.match_status.value_counts().to_dict(),
                             "realized_width_ratio_range": [float(realized.min()), float(realized.max())],
                             "paired_T5_vs_M0": an.paired_T5(dm, d0)}
        if m != bd and bd in hyp:
            out["methods"][m][f"paired_T5_vs_{bd}"] = an.paired_T5(dm, pick(rm, bd))
    # W status, v0.7 rule: (a), (b), (d) under randomised matching; (c) under nominal alpha
    dm = pick(rm, prim)
    st, det = an.w_status(dm, d0, an.paired_T5(dm, d0))
    nom = cov_nominal[(cov_nominal.method == prim) & (cov_nominal.hyper == str(hyp[prim]))]
    c_nom = 0.88 <= float(nom.covered.mean()) <= 0.92
    a, b, d = det["a"], det["b"], det["d"]
    if a and b and c_nom and d:
        w7 = "W+"
    elif a and b and c_nom:
        w7 = "W±"
    elif det.get("T5_cov") is not None and out["methods"][prim]["paired_T5_vs_M0"]["lo"] <= 0 <= \
            out["methods"][prim]["paired_T5_vs_M0"]["hi"]:
        w7 = "W-"
    elif out["methods"][prim]["paired_T5_vs_M0"]["diff"] < 0.02:
        w7 = "W-"
    else:
        w7 = "partial"
    out["W_v07_rule"] = {"status": w7, "a": a, "b": b, "c_nominal": bool(c_nom), "nominal_marginal": float(nom.covered.mean()),
                         "d": d, "primary": [prim, hyp[prim]], "W_v062_rule_under_rand_matching": st}
    return out


# ------------------------------------------------------------------ (3) confounding: bootstrap + Firth

def firth_logit(X, y, maxiter=100, tol=1e-8):
    """Firth (1993) bias-reduced logistic regression (Heinze & Schemper 2002 iteration with step halving)."""
    n, p = X.shape
    beta = np.zeros(p)

    def pen_ll(b):
        eta = X @ b
        pr = 1 / (1 + np.exp(-eta))
        W = pr * (1 - pr)
        info = X.T @ (X * W[:, None])
        sign, logdet = np.linalg.slogdet(info)
        return float(np.sum(y * eta - np.logaddexp(0, eta)) + 0.5 * logdet), pr, W, info
    ll, pr, W, info = pen_ll(beta)
    for it in range(maxiter):
        inv = np.linalg.inv(info)
        h = np.einsum("ij,jk,ik->i", X * np.sqrt(W)[:, None], inv, X * np.sqrt(W)[:, None])
        score = X.T @ (y - pr + h * (0.5 - pr))
        step = inv @ score
        for _ in range(30):
            nb = beta + step
            nll, npr, nW, ninfo = pen_ll(nb)
            if nll >= ll - 1e-12:
                break
            step /= 2
        beta, ll, pr, W, info = nb, nll, npr, nW, ninfo
        if np.max(np.abs(step)) < tol:
            break
    se = np.sqrt(np.diag(np.linalg.inv(info)))
    return beta, se, it + 1


def conf_bootstrap(d, cols, seed, B=N_BOOT_CONF):
    X = sm.add_constant(d[cols].astype(float).values, has_constant="add")
    y = d.miscover.values.astype(float)
    names = ["const"] + list(cols)

    def fit(Xs, ys):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                r = sm.Logit(ys, Xs).fit(disp=0, maxiter=200, method="newton")
                ok = bool(r.mle_retvals.get("converged", False)) and np.all(np.isfinite(r.params)) and \
                    np.max(np.abs(r.params)) < SEP_ABS_COEF
                return r.params, ok
            except Exception:                                           # noqa: BLE001
                return np.full(Xs.shape[1], np.nan), False
    point, point_ok = fit(X, y)
    groups = d.tss_group.values
    gidx = {g: np.where(groups == g)[0] for g in np.unique(groups)}
    keys = np.array(list(gidx))
    rng = np.random.default_rng(seed)
    boots, nfail = [], 0
    for _ in range(B):
        pick = rng.choice(keys, len(keys), replace=True)
        ix = np.concatenate([gidx[g] for g in pick])
        b, ok = fit(X[ix], y[ix])
        if ok:
            boots.append(b)
        else:
            nfail += 1
    boots = np.array(boots)
    fb, fse, fit_iter = firth_logit(X, y)
    # if more than 5% of plain-logistic bootstrap fits fail or separate, the percentile CI is unreliable;
    # then a Firth bootstrap (B = N_BOOT_FIRTH, same resampling scheme) is added and reported alongside.
    fboot = None
    if nfail > 0.05 * B:
        rng2 = np.random.default_rng(seed + 1)
        fb_list = []
        for _ in range(N_BOOT_FIRTH):
            pick = rng2.choice(keys, len(keys), replace=True)
            ix = np.concatenate([gidx[g] for g in pick])
            fb_list.append(firth_logit(X[ix], y[ix], maxiter=50)[0])
        fboot = np.array(fb_list)
    out = {}
    for j, nm in enumerate(names):
        col = boots[:, j] if len(boots) else np.array([np.nan])
        out[nm] = {"coef": float(point[j]), "point_fit_ok": bool(point_ok),
                   "boot_lo": float(np.nanquantile(col, 0.025)), "boot_hi": float(np.nanquantile(col, 0.975)),
                   "frac_boot_positive": float(np.nanmean(col > 0)),
                   "sign_supported": bool(np.nanquantile(col, 0.025) > 0 or np.nanquantile(col, 0.975) < 0),
                   "firth_coef": float(fb[j]), "firth_lo": float(fb[j] - 1.96 * fse[j]), "firth_hi": float(fb[j] + 1.96 * fse[j])}
        if fboot is not None:
            out[nm].update({"firth_boot_lo": float(np.quantile(fboot[:, j], 0.025)),
                            "firth_boot_hi": float(np.quantile(fboot[:, j], 0.975)),
                            "firth_boot_frac_positive": float(np.mean(fboot[:, j] > 0))})
    return {"coefs_per_SD_logodds": out, "n_rows": int(len(d)), "n_genes": int(d.target.nunique()),
            "n_tss_groups": int(len(keys)), "B": B, "n_boot_failed_or_separated": int(nfail), "firth_iterations": int(fit_iter),
            "plain_bootstrap_reliable": bool(nfail <= 0.05 * B), "firth_bootstrap_B": None if fboot is None else int(len(fboot)),
            "miscoverage_rate": float(y.mean())}


def confounding_v07(cov, covs, half=False):
    if half:
        rows = cov[(cov.k == ev.PRIMARY_K) & (cov.variant == "half") & (cov["mode"] == "nominal") & (cov.method == "M0")
                   & cov.model.isin(ev.LINEAR) & (cov.perm == -1)]
    else:
        rows = cov[(cov.k == ev.PRIMARY_K) & (cov.variant == "full") & (cov["mode"] == "nominal") & (cov.method == "M0")
                   & cov.model.isin(ev.LINEAR) & (cov.perm == -1)]
    fr = an.confounding_frame(rows, covs, half=half)
    out = {}
    for i, (m, x) in enumerate(fr.groupby("model")):
        out[m] = {"unadjusted": conf_bootstrap(x, ["zlogE"], BOOT_SEED + 10 * i + (5 if half else 0)),
                  "adjusted": conf_bootstrap(x, an.CONF_COLS, BOOT_SEED + 10 * i + 1 + (5 if half else 0))}
    return out


# ------------------------------------------------------------------ (4) M_dist by prior-distance strata

def dist_strata(cov, covs, summ, model=ev.PRIMARY_MODEL, n=5):
    hyp = summ["selection"]["hypers"]
    prim = summ["selection"]["primary"][0]
    base = cov[(cov.model == model) & (cov.k == ev.PRIMARY_K) & (cov.variant == "full") & (cov.perm == -1)
               & cov["mode"].isin(["nominal", "matched"])]
    base = base.merge(covs[["split", "target", "D_STRING"]], on=["split", "target"], how="left", validate="many_to_one")
    meths = [("M0", "None"), (prim, str(hyp[prim]))] + [(m, str(hyp[m])) for m in ["Mscore-N"] + an.BASELINES
                                                          if m in hyp and m != prim]
    out = {}
    for dist in ("D_GO", "D_STRING"):
        x = base.copy()
        x["q"] = x.groupby(["split", "mode", "method", "hyper"])[dist].transform(
            lambda s: pd.qcut(s.rank(method="first"), n, labels=False) + 1)
        iso = np.isclose(x[dist], 1.0)
        res = {"n_isolated_unit_rows": int(iso[(x.method == "M0") & (x["mode"] == "nominal")].sum()),
               "ties_note": "quintiles by within-split rank (ties broken by order); D = 1 (no neighbour) also reported"}
        for mode in ("nominal", "matched"):
            d0 = x[(x["mode"] == mode) & (x.method == "M0")].sort_values(["split", "target"])
            for meth, hy in meths:
                dm = x[(x["mode"] == mode) & (x.method == meth) & (x.hyper == hy)].sort_values(["split", "target"])
                if not len(dm):
                    continue
                g = dm.groupby("q").covered.mean()
                r = {"cov_by_quintile": {int(a): float(b) for a, b in g.items()},
                     "max_dev": float((g - an.NOMINAL).abs().max()),
                     "cov_isolated": float(dm.covered[np.isclose(dm[dist], 1.0)].mean()) if np.isclose(dm[dist], 1.0).any() else None,
                     "width_ratio_q5": float(dm.w[dm.q == n].mean() / d0.w[d0.q == n].mean())}
                if meth != "M0":
                    r["paired_q5_vs_M0"] = an.paired_T5(dm.assign(tier=dm.q), d0.assign(tier=d0.q), tier=n)
                    r["paired_q1_vs_M0"] = an.paired_T5(dm.assign(tier=dm.q), d0.assign(tier=d0.q), tier=1)
                res[f"{mode}|{meth}|{hy}"] = r
        out[dist] = res
    return out


def run_all(out_dir, caches, log=print):
    from pathlib import Path
    out_dir = Path(out_dir)
    summ = json.loads((out_dir / "summary.json").read_text(encoding="utf8"))
    cov = pd.read_parquet(out_dir / "coverage_rows.parquet")
    covs = pd.read_parquet(out_dir / "covariates_analysis_only.parquet")
    sel = {"h": tuple(summ["selection"]["h"].split("|")), "d": summ["selection"]["d"], "s": summ["selection"]["s"]}
    hyp = summ["selection"]["hypers"]
    grid = [(m, hyp[m]) for m in an.TUNED if m in hyp]
    res = {"protocol": "0.6.2 results unchanged; these are v0.7-style sensitivity analyses and the user-required M_dist strata"}
    for model in (ev.PRIMARY_MODEL, ev.SECONDARY_MODEL):
        rm = pd.concat([rand_match_rows(c, model, ev.PRIMARY_K, sel, grid) for c in caches], ignore_index=True)
        rm.to_parquet(out_dir / f"sens_rand_match_rows_{model}.parquet")
        nom = cov[(cov.model == model) & (cov.k == ev.PRIMARY_K) & (cov.variant == "full") & (cov["mode"] == "nominal")
                  & (cov.perm == -1)]
        res[f"rand_match|{model}"] = rand_match_summary(rm, nom, summ, model)
        log(f"sens: randomised matching done for {model}")
    res["confounding_bootstrap"] = confounding_v07(cov, covs)
    res["confounding_bootstrap_half"] = confounding_v07(cov, covs, half=True)
    log("sens: confounding bootstrap done")
    res["dist_strata"] = dist_strata(cov, covs, summ)
    log("sens: distance strata done")
    (out_dir / "sens_v07.json").write_text(json.dumps(res, indent=1, default=_js), encoding="utf8")
    return res


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)
