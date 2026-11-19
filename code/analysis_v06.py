# analysis_v06.py -- protocol v0.6.1 sections 5-10 (v0.6 sections 6-8): gene-level analysis, width-matched
# comparisons, selection rules 7.1 / 7.1b / 7.2 / 7.3, confounding and incremental-value analyses, criteria.
# Input: coverage rows from eval_v06.apply_methods (one row per split x test unit x method/hyper) and
# proxy rows from eval_v06.proxy_table. Units of inference are genes; clustering is by TSS group.
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit
from scipy.stats import norm, spearmanr
from sklearn.metrics import roc_auc_score

NOMINAL = 0.90
N_BOOT = 2000
WIDTH_CAP = 1.25
TIE_HYPER = 0.005          # 0.5 pp: "tie" tolerance for hyperparameter selection (rule 7.2)
TIE_METHOD = 0.01          # 1 pp: rule 7.3
TIE_RHO = 0.02             # rule 7.1
_G = ["1.0", "0.5", "0.25", "0.0"]
SIMPLER = {"M1": _G, "M2": ["3", "5"], "M3": ["quad", "ns3", "ns5"], "Mdist-N": _G, "Mdist-M": ["3", "5"],
           "Mscore-N": _G, "Mcomb-N": _G, "Mcomb-3": ["lin2"], "oracle": _G,     # simplest first
           "MD-GO": ["5"], "MD-STR": ["5"], "Munc-N": _G, "Munc-M": ["3", "5"]}   # v0.7 development additions
METHOD_ORDER = ["M1", "M2", "M3"]
TUNED = ["M1", "M2", "M3", "Mdist-N", "Mdist-M", "Mscore-N", "Mcomb-N", "Mcomb-3", "oracle",
         "MD-GO", "MD-STR", "Munc-N", "Munc-M"]
# v0.7 development: distance / difficulty baselines; the strongest (highest width-matched T5, no tie rule)
# becomes the comparator of step 2 of the hierarchical test
BASELINES = ["Mdist-N", "Mdist-M", "MD-GO", "MD-STR", "Munc-N", "Munc-M"]


# ------------------------------------------------------------------ descriptive summaries

def tier_table(df):
    """Pooled coverage per tier (= independence-GEE point estimate); strong (T5) / mid (T2-T4) / weak (T1);
    worst group and gap over T1..T5; capped widths (column w); fraction of infinite thresholds."""
    g = df.groupby("tier").covered.mean()
    grp = df.tier.map(lambda t: "strong" if t == 5 else ("weak" if t == 1 else "mid"))
    g3 = df.groupby("tier3").covered.mean() if "tier3" in df else pd.Series(dtype=float)
    return {"marginal": float(df.covered.mean()), "tier_cov": {int(k): float(v) for k, v in g.items()},
            "group_cov": df.groupby(grp).covered.mean().to_dict(),
            "tertile_cov": {int(k): float(v) for k, v in g3.items()},
            "worst": float(g.min()), "gap": float(g.max() - g.min()), "max_dev": float((g - NOMINAL).abs().max()),
            "width": float(df.w.mean()) if "w" in df else float(df.thr.mean()),
            "width_median": float(df.w.median()) if "w" in df else float(df.thr.median()),
            "width_T5": float(df.w[df.tier == 5].mean()) if "w" in df else float(df.thr[df.tier == 5].mean()),
            "frac_inf": float(np.isinf(df.thr).mean()),
            "alpha_used_median": float(df.alpha_used.median()) if "alpha_used" in df else None,
            "width_resid_max": float(df.width_resid.max()) if "width_resid" in df else None}


def strata_table(df, col, n=5):
    """Coverage by per-split quantile strata of a covariate (e.g. true D_GO; used for the dist scenario)."""
    q = df.groupby("split")[col].transform(lambda x: pd.qcut(x.rank(method="first"), n, labels=False))
    g = df.groupby(q).covered.mean()
    return {"cov": {int(k): float(v) for k, v in g.items()}, "max_dev": float((g - NOMINAL).abs().max())}


def gee_tiers(df):
    """Logistic GEE, independence working correlation, clusters = TSS groups, one intercept per tier.
    Returns coverage and 95% CI per tier. Tiers with all-covered or none-covered rows (separation) are
    reported with a cluster-bootstrap CI instead and flagged."""
    out = {}
    for t, d in df.groupby("tier"):
        p = d.covered.mean()
        if 0 < p < 1:
            X = np.ones((len(d), 1))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = sm.GEE(d.covered.astype(float).values, X, groups=d.tss_group.values,
                           family=sm.families.Binomial(), cov_struct=sm.cov_struct.Independence()).fit()
            b, se = float(r.params[0]), float(r.bse[0])
            out[int(t)] = {"cov": float(expit(b)), "lo": float(expit(b - 1.96 * se)), "hi": float(expit(b + 1.96 * se)),
                           "n": int(len(d)), "method": "gee"}
        else:
            lo, hi = cluster_boot(d, lambda x: x.covered.mean(), seed=int(t))
            out[int(t)] = {"cov": float(p), "lo": lo, "hi": hi, "n": int(len(d)), "method": "bootstrap (separation)"}
    return out


def paired_T5(df_m, df_0, tier=5):
    """Coverage difference (method - M0) in one tier on the same (split, gene) rows.
    Stacked logistic GEE with a method indicator, clusters = TSS groups; difference and 95% CI by the
    delta method. Falls back to an identity-link GEE if either arm is all-covered or none-covered."""
    a = df_m[df_m.tier == tier][["split", "target", "tss_group", "covered"]].assign(m=1.0)
    b = df_0[df_0.tier == tier][["split", "target", "tss_group", "covered"]].assign(m=0.0)
    assert len(a) == len(b), "paired arms differ in size"
    d = pd.concat([a, b], ignore_index=True)
    pa, pb = a.covered.mean(), b.covered.mean()
    X = np.column_stack([np.ones(len(d)), d.m.values])
    y = d.covered.astype(float).values
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if 0 < pa < 1 and 0 < pb < 1:
            r = sm.GEE(y, X, groups=d.tss_group.values, family=sm.families.Binomial(),
                       cov_struct=sm.cov_struct.Independence()).fit()
            b0, b1 = r.params
            V = r.cov_params()
            p1, p0 = expit(b0 + b1), expit(b0)
            grad = np.array([p1 * (1 - p1) - p0 * (1 - p0), p1 * (1 - p1)])
            diff, se, how = p1 - p0, float(np.sqrt(grad @ V @ grad)), "logistic GEE + delta"
        else:
            r = sm.GEE(y, X, groups=d.tss_group.values, family=sm.families.Gaussian(),
                       cov_struct=sm.cov_struct.Independence()).fit()
            diff, se, how = float(r.params[1]), float(r.bse[1]), "identity-link GEE (separation)"
    return {"diff": float(diff), "lo": float(diff - 1.96 * se), "hi": float(diff + 1.96 * se),
            "p": float(2 * norm.sf(abs(diff / se))) if se > 0 else np.nan, "n_rows": int(len(a)), "how": how}


def cluster_boot(df, stat, n=N_BOOT, seed=0, level=0.95):
    """Percentile CI of stat(df) under resampling of TSS groups.
    2026-10-10 speed-up (pure code change, user-approved): the resampled frame is built with one positional
    take instead of concatenating one small frame per group. Same RNG calls in the same order and the same
    rows in the same order (rows of a group keep their original order; groups in drawn order), so every
    resampled frame and hence every statistic is identical to the previous implementation."""
    rng = np.random.default_rng(seed)
    gs = df.tss_group.unique()
    pos = df.groupby("tss_group", sort=False).indices
    vals = []
    for _ in range(n):
        pick = rng.choice(gs, len(gs), replace=True)
        ix = np.concatenate([pos[g] for g in pick])
        vals.append(stat(df.take(ix).reset_index(drop=True)))
    a = (1 - level) / 2
    return float(np.nanquantile(vals, a)), float(np.nanquantile(vals, 1 - a))


def gene_level(df):
    """Gene-level coverage frequency: mean coverage per gene over its test appearances; gene tier =
    median of its per-split tiers (rounded down); then mean per tier, CI by TSS-group bootstrap."""
    g = df.groupby(["target", "tss_group"]).agg(cov=("covered", "mean"), tier=("tier", "median")).reset_index()
    g["tier"] = np.floor(g.tier).astype(int)
    out = {}
    for t, d in g.groupby("tier"):
        lo, hi = cluster_boot(d, lambda x: x["cov"].mean(), seed=100 + int(t))
        out[int(t)] = {"cov": float(d["cov"].mean()), "lo": lo, "hi": hi, "n_genes": int(len(d))}
    return out


def tier_stability(df):
    """Fraction of genes whose tier differs between their test appearances (per-split reference)."""
    n = df.groupby("target").tier.nunique()
    return float((n > 1).mean())


# ------------------------------------------------------------------ predictability (question i)

def gene_proxy(pt):
    """Per gene: mean h and mean E_true over its test appearances; tier = median tier."""
    g = pt.groupby(["target", "tss_group"]).agg(h=("h", "mean"), E=("E_true", "mean"), tier=("tier", "median")).reset_index()
    g["T5"] = (np.floor(g.tier) == 5).astype(int)
    return g


def predictability(pt, boot=True, seed=0):
    g = gene_proxy(pt)
    rho = float(spearmanr(g.h, g.E).correlation)
    auc = float(roc_auc_score(g.T5, g.h)) if g.T5.nunique() == 2 else np.nan
    out = {"rho": rho, "auroc_T5": auc, "n_genes": int(len(g))}
    if boot:
        out["rho_ci"] = cluster_boot(g, lambda x: spearmanr(x.h, x.E).correlation, seed=seed)
        out["auroc_ci"] = cluster_boot(g, lambda x: roc_auc_score(x.T5, x.h) if x.T5.nunique() == 2 else np.nan,
                                       seed=seed + 1)
    return out


def decision_metrics(df, fracs=(0.10, 0.20)):
    """Per split: test units ranked by h; top 10% / 20%: coverage, mean width, precision for T5."""
    rows = []
    for (sp,), d in df.groupby(["split"]):
        d = d.sort_values("h", ascending=False)
        for f in fracs:
            top = d.iloc[:max(1, int(np.ceil(f * len(d))))]
            rows.append({"split": sp, "frac": f, "cov": top.covered.mean(), "width": top.thr.mean(),
                         "precision_T5": (top.tier == 5).mean()})
    r = pd.DataFrame(rows)
    return {f: {"cov": float(x["cov"].mean()), "width": float(x.width.mean()), "precision_T5": float(x.precision_T5.mean()),
                "cov_range": (float(x["cov"].min()), float(x["cov"].max()))} for f, x in r.groupby("frac")}


def decision_paired(df_m, df_0, frac=0.10):
    """v0.7 development: per split, the ceil(frac x n) test units with the largest h (ties: first in
    (split, target) order; the same units for both arms, since the ranking uses h only) -> paired coverage
    difference method - M0 with the logistic GEE of paired_T5 (clusters = TSS groups). Also returns both
    coverages and mean capped widths. If both arms are constant, the difference is reported without a CI."""
    def top(d):
        d = d.sort_values(["split", "target"]).reset_index(drop=True)
        r = d.groupby("split").h.rank(ascending=False, method="first")
        n = d.groupby("split").h.transform("size")
        return d[r <= np.maximum(1, np.ceil(frac * n))]
    a, b = top(df_m), top(df_0)
    assert (a[["split", "target"]].values == b[["split", "target"]].values).all(), "arms select different units"
    out = {"cov_method": float(a.covered.mean()), "cov_M0": float(b.covered.mean()),
           "width_method": float(a.w.mean()), "width_M0": float(b.w.mean()), "frac": frac,
           "n_rows": int(len(a)), "frac_T5_in_selection": float((b.tier == 5).mean())}
    if a.covered.nunique() == 1 and b.covered.nunique() == 1:
        out.update({"diff": out["cov_method"] - out["cov_M0"], "lo": np.nan, "hi": np.nan, "p": np.nan,
                    "how": "degenerate (both arms constant)"})
        return out
    r = paired_T5(a.assign(tier=5), b.assign(tier=5))
    r.update(out)
    return r

# ------------------------------------------------------------------ selection rules (v0.6.1 section 9; dev only)

def _gene_rho(d, x, y):
    g = d.groupby("target").agg(x=(x, "mean"), y=(y, "mean"))
    return float(spearmanr(g.x, g.y).correlation)


def select_proxy(pt):
    """Rule 7.1 on the primary model: highest gene-level Spearman rho(h, E_true); within TIE_RHO of the best,
    the configuration with the fewest substantive features, then ridge before gbm."""
    from proxy_v06 import n_substantive
    res = {c: _gene_rho(d, "value", "E_true") for c, d in pt.groupby("proxy")}
    best = max(res.values())
    cand = sorted([c for c, r in res.items() if r >= best - TIE_RHO],
                  key=lambda c: (n_substantive(c.split("|")[0]), 0 if c.endswith("ridge") else 1, -res[c]))
    return cand[0], res


def select_learner(pt):
    """Rule 7.1b (d_hat, and s_hat by the same rule): highest gene-level rho(proxy, S1); ties -> ridge."""
    res = {ln: _gene_rho(d, "value", "score") for ln, d in pt.groupby("proxy")}
    best = max(res.values())
    return ("ridge" if res.get("ridge", -9) >= best - TIE_RHO else max(res, key=res.get)), res


def _t5(d):
    return float(d.covered[d.tier == 5].mean())


def select_hyper(cov_matched):
    """Rule 7.2 on width-matched rows (primary model, k = 50, full data): per method, the hyperparameter with
    the highest T5 coverage; within TIE_HYPER of the best, the simpler value."""
    out = {}
    for meth in TUNED:
        d = cov_matched[cov_matched.method == meth]
        if not len(d):
            continue
        t5 = {hy: _t5(x) for hy, x in d.groupby("hyper")}
        best = max(t5.values())
        chosen = next(h for h in SIMPLER[meth] if h in t5 and t5[h] >= best - TIE_HYPER)
        out[meth] = {"chosen": chosen, "T5_matched": t5}
    return out


def select_primary(hyp, methods=METHOD_ORDER, tie=TIE_METHOD):
    """Rule 7.3: among `methods` (chosen hypers), highest matched T5 coverage; within `tie`, earlier in order."""
    cand = {m: hyp[m]["T5_matched"][hyp[m]["chosen"]] for m in methods if m in hyp}
    best = max(cand.values())
    return next(m for m in methods if m in cand and cand[m] >= best - tie)


# ------------------------------------------------------------------ confounding (v0.6.1 section 6)

def _z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / (x.std() + 1e-12)


CONF_COLS = ["zlogE", "rel", "zlognoise", "zlogcells", "zDGO", "zDSTR", "zbase", "zsim", "zshared"]


def confounding_frame(rows_M0, covs, half=False):
    """rows_M0: M0 nominal rows (linear models); covs: per-split covariate frames (with split column)."""
    d = rows_M0.merge(covs, on=["split", "target", "tss_group"], how="left", validate="many_to_one")
    E = d.E_true_A if half else d.E_true_y
    d = d.assign(miscover=(~d.covered).astype(float), zlogE=_z(np.log(E)), zlognoise=_z(np.log(d.noise)),
                 zlogcells=_z(np.log(d.ncells.astype(float))), zDGO=_z(d.D_GO_y), zDSTR=_z(d.D_STRING),
                 zbase=_z(d.f_base), zsim=_z(d.sim_maxJ), zshared=_z(d.shared), rel=_z(d.rel))
    return d


def _gee_logit(d, cols):
    X = sm.add_constant(d[cols].astype(float).values)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = sm.GEE(d.miscover.values, X, groups=d.tss_group.values, family=sm.families.Binomial(),
                   cov_struct=sm.cov_struct.Independence()).fit()
    names = ["const"] + list(cols)
    ci = r.conf_int()
    return {n: {"OR_perSD": float(np.exp(r.params[i])), "lo": float(np.exp(ci[i, 0])), "hi": float(np.exp(ci[i, 1])),
                "p": float(r.pvalues[i])} for i, n in enumerate(names)}


def vif(d, cols):
    X = d[cols].astype(float).values
    out = {}
    for j, c in enumerate(cols):
        others = np.delete(X, j, axis=1)
        A = np.column_stack([np.ones(len(X)), others])
        res = X[:, j] - A @ np.linalg.lstsq(A, X[:, j], rcond=None)[0]
        r2 = 1 - res.var() / X[:, j].var()
        out[c] = float(1 / max(1e-12, 1 - r2))
    return out


def confounding(d):
    """Per linear model: unadjusted and adjusted log-E coefficients; pooled model with model x logE."""
    out = {}
    for m, x in d.groupby("model"):
        if x.miscover.nunique() < 2:
            out[m] = {"note": "no variation in miscoverage"}
            continue
        out[m] = {"unadjusted": _gee_logit(x, ["zlogE"])["zlogE"], "adjusted": _gee_logit(x, CONF_COLS),
                  "vif": vif(x, CONF_COLS), "n_rows": int(len(x)), "n_genes": int(x.target.nunique())}
    d2 = d.assign(is_go=(d.model == "linear_go").astype(float))
    d2 = d2.assign(go_x_logE=d2.is_go * d2.zlogE)
    if d2.miscover.nunique() == 2 and d2.model.nunique() == 2:
        out["pooled"] = _gee_logit(d2, CONF_COLS + ["is_go", "go_x_logE"])
    return out


# ------------------------------------------------------------------ incremental value (v0.6.1 section 5)

def incremental(rows_M0, covs, boot=True):
    """(a) gene-level Spearman of h with d, D_GO, D_STRING, sim_maxJ, predicted size; (b1) GEE of log E on
    z log h + z log d; (b2) logistic GEE of M0 miscoverage on z log h + z log d per model, with the gene-level
    AUROC gain of adding h to d (cluster bootstrap CI)."""
    d = rows_M0.merge(covs, on=["split", "target", "tss_group"], how="left", validate="many_to_one")
    out = {}
    for m, x in d.groupby("model"):
        g = x.groupby(["target", "tss_group"]).agg(h=("h", "mean"), dd=("d", "mean"), E=("E_true_x", "mean"),
                                                   DGO=("D_GO_y", "mean"), DSTR=("D_STRING", "mean"),
                                                   sim=("sim_maxJ", "mean"), mis=("covered", lambda c: 1 - c.mean()),
                                                   size=(f"pred_size_{m}" if f"pred_size_{m}" in x else "h", "mean")
                                                   ).reset_index()
        a = {c: float(spearmanr(g.h, g[c]).correlation) for c in ("dd", "DGO", "DSTR", "sim", "size")}
        xx = x.assign(zh=_z(np.log(x.h)), zd=_z(np.log(x.d)), logE=np.log(x.E_true_x),
                      miscover=(~x.covered).astype(float))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r1 = sm.GEE(xx.logE.values, sm.add_constant(xx[["zh", "zd"]].values), groups=xx.tss_group.values,
                        family=sm.families.Gaussian(), cov_struct=sm.cov_struct.Independence()).fit()
        ci1 = r1.conf_int()
        b1 = {"h": [float(r1.params[1]), float(ci1[1, 0]), float(ci1[1, 1])],
              "d": [float(r1.params[2]), float(ci1[2, 0]), float(ci1[2, 1])]}
        b2 = _gee_logit(xx, ["zh", "zd"]) if xx.miscover.nunique() == 2 else None
        g = g.assign(misbin=(g.mis > 0.5).astype(int), lh=np.log(g.h), ld=np.log(g.dd))

        def dauc(gg):
            if gg.misbin.nunique() < 2:
                return np.nan
            from sklearn.linear_model import LogisticRegression
            y = gg.misbin.values
            p1 = LogisticRegression().fit(gg[["ld"]], y).predict_proba(gg[["ld"]])[:, 1]
            p2 = LogisticRegression().fit(gg[["ld", "lh"]], y).predict_proba(gg[["ld", "lh"]])[:, 1]
            return roc_auc_score(y, p2) - roc_auc_score(y, p1)
        out[m] = {"a_spearman_h_vs": a, "b1_logE": b1, "b2_miscover": b2, "dAUROC_add_h": float(dauc(g)),
                  "dAUROC_ci": cluster_boot(g, dauc, n=500, seed=3) if boot else None,
                  "note": "dAUROC is in-sample (gene-level logistic, d vs d + h); descriptive"}
    return out


# ------------------------------------------------------------------ criteria (section 8; confirmation)

def w_status(df_m, df_0, paired):
    tm, t0 = tier_table(df_m), tier_table(df_0)
    a = (tm["tier_cov"].get(5, np.nan) >= 0.85) and (paired["diff"] >= 0.05) and (paired["lo"] > 0)
    b = tm["max_dev"] <= 0.5 * t0["max_dev"]
    c = 0.88 <= tm["marginal"] <= 0.92
    d = (tm["width"] / t0["width"]) <= WIDTH_CAP
    if a and b and c and d:
        st = "W+"
    elif a and b and c:
        st = "W±"
    elif paired["lo"] <= 0 <= paired["hi"] or paired["diff"] < 0.02:
        st = "W-"
    else:
        st = "partial"
    return st, {"a": bool(a), "b": bool(b), "c": bool(c), "d": bool(d), "T5_cov": tm["tier_cov"].get(5),
                "max_dev": tm["max_dev"], "max_dev_M0": t0["max_dev"], "marginal": tm["marginal"],
                "width_ratio": tm["width"] / t0["width"]}


def criteria(cov, pt, primary_method, primary_hyper, model, k):
    """Section 8 classification for one dataset. cov: coverage rows incl. M0, the primary method (perm = -1)
    and its 20 permuted-h runs (perm = 0..19). pt: proxy rows of the selected proxy config."""
    sel = lambda m, hy, perm: cov[(cov.model == model) & (cov.k == k) & (cov.method == m)          # noqa: E731
                                  & (cov.hyper == str(hy)) & (cov.perm == perm)]
    d0 = sel("M0", None, -1).sort_values(["split", "target"])
    dm = sel(primary_method, primary_hyper, -1).sort_values(["split", "target"])
    g0 = gee_tiers(d0)
    B0 = g0[5]["cov"] < 0.85 and g0[5]["hi"] < 0.90
    pr = predictability(pt[pt.model == model])
    rlo, rhi = pr["rho_ci"]
    alo, ahi = pr["auroc_ci"]
    P = "P+" if (rlo >= 0.30 and alo >= 0.65) else ("P-" if (rhi < 0.20 or ahi < 0.60) else "P0")
    paired = paired_T5(dm, d0)
    W, wdet = w_status(dm, d0, paired)
    # permuted-h negative control: average over the 20 permutations, plus their distribution
    perm_rows = []
    for r in sorted(cov.perm.unique()):
        if r < 0:
            continue
        dp = sel(primary_method, primary_hyper, r).sort_values(["split", "target"])
        if len(dp) == 0:
            continue
        pp = paired_T5(dp, d0)
        st, det = w_status(dp, d0, pp)
        perm_rows.append({"perm": int(r), "T5_diff": pp["diff"], "max_dev": det["max_dev"],
                          "width_ratio": det["width_ratio"], "W_status": st})
    perm_df = pd.DataFrame(perm_rows)
    if len(perm_df):
        # averaged permuted control: mean coverage indicator per (split, gene) over the 20 runs
        allp = cov[(cov.model == model) & (cov.k == k) & (cov.method == primary_method)
                   & (cov.hyper == str(primary_hyper)) & (cov.perm >= 0)]
        avg = allp.groupby(["split", "target", "tss_group", "tier"], as_index=False).agg(covered=("covered", "mean"),
                                                                                    thr=("thr", "mean"))
        avg = avg.sort_values(["split", "target"])
        pa = paired_T5_mean(avg, d0)
        W_perm, wdet_perm = w_status(avg, d0, pa)
    else:
        W_perm, wdet_perm, pa = None, None, None
    attributable = W_perm != "W+"
    if not B0:
        overall = "effect-strength-driven conditional failure not reproduced (B0 fails); (ii) descriptive only"
    elif P == "P-":
        overall = "effect strength not predictable; conditional coverage hard to achieve"
    elif P == "P+" and W in ("W+", "W±") and attributable:
        overall = "method effective" + (" (inefficient: width > 1.25 x M0)" if W == "W±" else "")
    elif P == "P+" and W in ("W+", "W±") and not attributable:
        overall = "improvement not attributable to h (permuted-h control also meets W+)"
    elif P == "P+" and W == "W-":
        overall = "effect strength predictable but insufficient to repair conditional coverage"
    else:
        overall = f"partial: B0 holds, {P}, {W}"
    return {"B0": bool(B0), "M0_T5": g0[5], "P": P, "predictability": pr, "W": W, "W_detail": wdet,
            "paired_T5": paired, "W_perm_mean": W_perm, "W_perm_detail": wdet_perm, "paired_T5_perm_mean": pa,
            "perm_distribution": perm_df.to_dict("records"),
            "perm_summary": ({"n": int(len(perm_df)), "n_meeting_W+": int((perm_df.W_status == "W+").sum()),
                              **{f"{c}_median": float(perm_df[c].median()) for c in ("T5_diff", "max_dev", "width_ratio")},
                              **{f"{c}_range": (float(perm_df[c].min()), float(perm_df[c].max()))
                                 for c in ("T5_diff", "max_dev", "width_ratio")}} if len(perm_df) else None),
            "overall": overall}


def paired_T5_mean(df_avg, df_0, tier=5):
    """Paired T5 difference for the permutation-averaged control (coverage is a fraction in [0, 1]);
    identity-link GEE, clusters = TSS groups."""
    a = df_avg[df_avg.tier == tier][["tss_group", "covered"]].assign(m=1.0)
    b = df_0[df_0.tier == tier][["tss_group", "covered"]].assign(m=0.0)
    d = pd.concat([a, b], ignore_index=True)
    X = np.column_stack([np.ones(len(d)), d.m.values])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = sm.GEE(d.covered.astype(float).values, X, groups=d.tss_group.values, family=sm.families.Gaussian(),
                   cov_struct=sm.cov_struct.Independence()).fit()
    diff, se = float(r.params[1]), float(r.bse[1])
    return {"diff": diff, "lo": diff - 1.96 * se, "hi": diff + 1.96 * se, "how": "identity-link GEE on averaged coverage"}
