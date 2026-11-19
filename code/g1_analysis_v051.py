# g1_analysis_v051.py -- G1 statistics for protocol v0.5.1. Identical to g1_analysis_v05.py except the
# 'support H' rule: C5 (audit-B marginal coverage) is a secondary outcome, not a requirement.
# Writes g1_results_v0.5.1.json (the v0.5 output g1_results_v0.5.json is left untouched).
#   python g1_analysis_v05.py prep     -> results/g1_v0.5/analysis/glmm_input.parquet (+ self-check, shift)
#   Rscript g1_glmm_v05.R              -> results/g1_v0.5/analysis/glmm_results.json
#   python g1_analysis_v05.py analyze  -> GEE, audit-B marginal coverage, sensitivity, criteria (section 8)
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.special import expit
from scipy.stats import binom, norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import RESULTS  # noqa: E402
from audit_lib import betabinom_band  # noqa: E402

G1 = RESULTS / "g1_v0.5"
AN = G1 / "analysis"
ALPHA, DELTA = 0.10, 0.05
N_A, N_B = 100, 20
PRIMARY_SCORE = "S1_k50"
PRIMARY_MODELS = ("mean", "linear_pca", "linear_go")
REPRO_MODELS = ("linear_pca_K10_repro", "linear_go_d2_repro")
DISTS = ("D_GO", "D_STRING", "D_GO_filt")
COLS = ["audit", "split", "model", "score", "alpha", "target", "tss_group", "D_GO", "D_STRING", "D_GO_filt",
        "rank_D_GO", "rank_D_STRING", "rank_D_GO_filt", "covered", "qhat", "n_cal", "n_cells"]


def read_all():
    parts = []
    for audit, n in (("A", N_A), ("B", N_B)):
        files = sorted((G1 / audit).glob("split_*.parquet"))
        assert len(files) == n, f"audit {audit}: {len(files)} split files, expected {n}"
        parts += [pd.read_parquet(f, columns=COLS) for f in files]
    R = pd.concat(parts, ignore_index=True)
    R["miscover"] = (~R.covered).astype(int)
    R["log_ncells"] = np.log(R.n_cells)
    return R


def prep():
    AN.mkdir(parents=True, exist_ok=True)
    R = read_all()
    # ---- self-check (i) on audit A: marginal coverage vs beta-binomial band (all models/scores/alphas)
    A = R[R.audit == "A"]
    band = []
    for (sp, m, s, a), g in A.groupby(["split", "model", "score", "alpha"]):
        lo, hi = betabinom_band(int(g.n_cal.iloc[0]), len(g), a)
        k = int(g.covered.sum())
        band.append((sp, m, s, a, k, len(g), lo, hi, lo <= k <= hi))
    B = pd.DataFrame(band, columns=["split", "model", "score", "alpha", "n_covered", "n_test", "band_lo", "band_hi", "in_band"])
    B.to_csv(AN / "selfcheck_marginal_auditA.csv", index=False)
    n_out = int((~B.in_band).sum())
    sc = {"n_combos": len(B), "n_out_of_band": n_out, "p_excess_vs_5pct": float(binom.sf(n_out - 1, len(B), 0.05))}
    sc["pass"] = sc["p_excess_vs_5pct"] > 0.01
    # ---- GLMM input: primary score/alpha, primary + reproduction models
    P = R[(R.score == PRIMARY_SCORE) & (R.alpha == ALPHA) & R.model.isin(PRIMARY_MODELS + REPRO_MODELS)].copy()
    P["log_ncells_c"] = P.log_ncells - P.log_ncells.mean()
    P[["audit", "split", "model", "target", "tss_group", "miscover", "log_ncells_c",
       "rank_D_GO", "rank_D_STRING", "rank_D_GO_filt"]].to_parquet(AN / "glmm_input.parquet")
    # ---- actual audit-B shift: per split median test distance, then median over splits
    one = R[(R.model == "mean") & (R.score == PRIMARY_SCORE) & (R.alpha == ALPHA)]
    med = one.groupby(["audit", "split"])[list(DISTS)].median().groupby("audit").median()
    shift = {d: {"A": float(med.loc["A", d]), "B": float(med.loc["B", d]), "B_minus_A": float(med.loc["B", d] - med.loc["A", d])}
             for d in DISTS}
    qh = R[R.alpha == ALPHA].groupby(["audit", "model", "score"]).qhat.agg(["mean", "min", "max"]).reset_index()
    qh.to_csv(AN / "qhat_by_model.csv", index=False)
    out = {"selfcheck_marginal_auditA": sc, "auditB_shift": shift, "n_rows": int(len(R)),
           "n_rows_glmm_input": int(len(P))}
    (AN / "prep_summary.json").write_text(json.dumps(out, indent=1), encoding="utf8")
    print(json.dumps(out))


# ------------------------------------------------------------------ GEE helpers
def gee_trend(df, rank_col, with_ncells=True):
    """Population-averaged logistic GEE (independence, clusters = TSS group, robust SE)."""
    d = df.copy()
    d["rk"] = d[rank_col]
    f = "miscover ~ C(model, Treatment('mean')) * rk" + (" + log_ncells_c" if with_ncells else "")
    fit = smf.gee(f, groups="tss_group", data=d, family=sm.families.Binomial(),
                  cov_struct=sm.cov_struct.Independence()).fit()
    names = list(fit.params.index)
    V = fit.cov_params().values
    out = {}
    for m in sorted(d.model.unique()):
        main = "Intercept" if m == "mean" else f"C(model, Treatment('mean'))[T.{m}]"
        inter = None if m == "mean" else f"C(model, Treatment('mean'))[T.{m}]:rk"
        cs = np.zeros(len(names))
        cs[names.index("rk")] = 1
        if inter:
            cs[names.index(inter)] = 1
        b, se = float(cs @ fit.params.values), float(np.sqrt(cs @ V @ cs))
        res = {"slope": b, "se": se, "p": float(2 * norm.sf(abs(b / se)))}
        for r in (0.1, 0.9):
            c = np.zeros(len(names))
            c[names.index("Intercept")] = 1
            if m != "mean":
                c[names.index(main)] = 1
            c += r * cs
            eta, s_eta = float(c @ fit.params.values), float(np.sqrt(c @ V @ c))
            res[f"coverage_at_{r}"] = float(1 - expit(eta))
            res[f"coverage_at_{r}_ci"] = [float(1 - expit(eta + 1.96 * s_eta)), float(1 - expit(eta - 1.96 * s_eta))]
        out[m] = res
    return out


def gee_marginal(df):
    out = {}
    for m, g in df.groupby("model"):
        fit = smf.gee("miscover ~ 1", groups="tss_group", data=g, family=sm.families.Binomial(),
                      cov_struct=sm.cov_struct.Independence()).fit()
        eta, se = float(fit.params.iloc[0]), float(fit.bse.iloc[0])
        per_split = 1 - g.groupby("split").miscover.mean()
        out[m] = {"coverage": float(1 - expit(eta)), "ci": [float(1 - expit(eta + 1.96 * se)), float(1 - expit(eta - 1.96 * se))],
                  "per_split_mean": float(per_split.mean()),
                  "per_split_2.5_97.5": [float(per_split.quantile(0.025)), float(per_split.quantile(0.975))]}
    return out


def holm(pvals):
    keys = sorted(pvals, key=pvals.get)
    m, adj, run = len(keys), {}, 0.0
    for i, k in enumerate(keys):
        run = max(run, min(1.0, (m - i) * pvals[k]))
        adj[k] = run
    return adj


def analyze():
    R = read_all()
    R["log_ncells_c"] = R.log_ncells - R.log_ncells.mean()
    glmm = json.loads((AN / "glmm_results.json").read_text(encoding="utf8"))
    prep_s = json.loads((AN / "prep_summary.json").read_text(encoding="utf8"))
    prim = R[(R.score == PRIMARY_SCORE) & (R.alpha == ALPHA)]
    res = {"gee": {}, "auditB_marginal": None, "sensitivity_gee": {}}
    for audit in ("A", "B"):
        d = prim[(prim.audit == audit) & prim.model.isin(PRIMARY_MODELS)]
        for dist in DISTS:
            res["gee"][f"{audit}_{dist}"] = gee_trend(d, f"rank_{dist}")
        res["gee"][f"{audit}_D_GO_no_ncells"] = gee_trend(d, "rank_D_GO", with_ncells=False)
    res["auditA_marginal"] = gee_marginal(prim[(prim.audit == "A") & prim.model.isin(PRIMARY_MODELS + REPRO_MODELS)])
    res["auditB_marginal"] = gee_marginal(prim[(prim.audit == "B") & prim.model.isin(PRIMARY_MODELS + REPRO_MODELS)])
    # sensitivity (GEE, audit A, rank D_GO)
    A = R[R.audit == "A"]
    sens = {"S1_k20_a0.1": (A.score == "S1_k20") & (A.alpha == 0.10), "S1_k200_a0.1": (A.score == "S1_k200") & (A.alpha == 0.10),
            "S2_a0.1": (A.score == "S2") & (A.alpha == 0.10), "S1_k50_a0.05": (A.score == "S1_k50") & (A.alpha == 0.05),
            "S1_k50_a0.2": (A.score == "S1_k50") & (A.alpha == 0.20)}
    for lab, msk in sens.items():
        res["sensitivity_gee"][lab] = gee_trend(A[msk & A.model.isin(PRIMARY_MODELS)], "rank_D_GO")
    res["sensitivity_gee"]["repro_models_S1_k50_a0.1"] = gee_trend(
        prim[(prim.audit == "A") & prim.model.isin(("mean",) + REPRO_MODELS)], "rank_D_GO")
    for dist in ("D_STRING", "D_GO_filt"):
        res["sensitivity_gee"][f"repro_models_{dist}"] = gee_trend(
            prim[(prim.audit == "A") & prim.model.isin(("mean",) + REPRO_MODELS)], f"rank_{dist}")

    # ------------------------------------------------ criteria (protocol v0.5 section 8, G1 version)
    g = glmm["fits"]
    fam = {f"{m}|{dist}": g[f"A_{dist}_ncells"]["slopes"][m]["p"] for m in PRIMARY_MODELS for dist in ("D_GO", "D_STRING")}
    adj = holm(fam)
    slope = lambda fitkey, m: g[fitkey]["slopes"][m]["est"]
    crit = {"holm_family": {k: {"p": fam[k], "p_holm": adj[k],
                                "slope_miscover": slope(f"A_{k.split('|')[1]}_ncells", k.split('|')[0])} for k in fam}}
    per_model = {}
    for m in ("linear_pca", "linear_go"):
        ge = res["gee"]["A_D_GO"][m]
        c1 = ge["coverage_at_0.9"] < 1 - ALPHA - DELTA and ge["coverage_at_0.9_ci"][1] < 1 - ALPHA
        c2 = slope("A_D_GO_ncells", m) > 0 and adj[f"{m}|D_GO"] < 0.05
        c3 = np.sign(slope("A_D_GO_ncells", m)) == np.sign(slope("A_D_STRING_ncells", m))
        c5 = res["auditB_marginal"][m]["coverage"] < 1 - ALPHA - DELTA
        per_model[m] = {"C1_far_end_undercoverage": bool(c1), "C2_positive_miscover_slope_holm": bool(c2),
                        "C3_same_direction_GO_STRING": bool(c3), "C5_auditB_marginal_below_0.85": bool(c5),
                        "coverage_at_rank0.9_GEE": ge["coverage_at_0.9"], "ci": ge["coverage_at_0.9_ci"],
                        "auditB_coverage": res["auditB_marginal"][m]["coverage"]}
    d4 = g["A_D_GO_ncells"]["diffs"]
    c4 = d4["linear_go-mean"]["est"] > 0 and d4["linear_go-mean"]["p"] < 0.05
    crit["per_linear_model"] = per_model
    crit["C4_linear_go_minus_mean_slope"] = {"est": d4["linear_go-mean"]["est"], "p": d4["linear_go-mean"]["p"], "pass": bool(c4),
                                             "linear_pca-mean": d4["linear_pca-mean"]}
    # v0.5.1: support requires C1 & C2 & C3 for at least one linear model, and C4; C5 is secondary.
    support = any(all(v[k] for k in ("C1_far_end_undercoverage", "C2_positive_miscover_slope_holm",
                                     "C3_same_direction_GO_STRING")) for v in per_model.values()) and c4
    crit["C5_secondary"] = {m: per_model[m]["C5_auditB_marginal_below_0.85"] for m in per_model}
    crit["support_H_v0.5_rule_including_C5"] = bool(any(all(v[k] for k in (
        "C1_far_end_undercoverage", "C2_positive_miscover_slope_holm", "C3_same_direction_GO_STRING",
        "C5_auditB_marginal_below_0.85")) for v in per_model.values()) and c4)
    within = lambda x: 1 - ALPHA - DELTA <= x <= 1 - ALPHA + DELTA
    refute_bits = []
    for m in PRIMARY_MODELS:
        ge = res["gee"]["A_D_GO"][m]
        for r in (0.1, 0.9):
            lo, hi = ge[f"coverage_at_{r}_ci"]
            refute_bits.append(within(ge[f"coverage_at_{r}"]) and (hi - lo) / 2 < DELTA)
        refute_bits.append(within(res["auditB_marginal"][m]["coverage"]))
    refute_bits.append(all(v >= 0.05 for v in adj.values()))
    refute = all(refute_bits)
    crit["support_H"] = bool(support)
    crit["refute_H"] = bool(refute)
    crit["classification"] = "support" if support else ("refute" if refute else "neither (partial / inconclusive)")
    nn = g["A_D_GO_no_ncells"]["slopes"]
    wn = g["A_D_GO_ncells"]["slopes"]
    crit["reinterpretation"] = {
        "R1_mean_same_trend_no_extrapolation": bool(slope("A_D_GO_ncells", "mean") > 0 and adj["mean|D_GO"] < 0.05 and not c4),
        "R2_far_end_overcoverage": bool(any(crit["holm_family"][k]["slope_miscover"] < 0 and crit["holm_family"][k]["p_holm"] < 0.05 for k in fam)),
        "R3_trend_vanishes_with_log_ncells": {m: bool(nn[m]["est"] > 0 and nn[m]["p"] < 0.05 and wn[m]["p"] >= 0.05) for m in PRIMARY_MODELS},
        "R4_E_strong_only": "not assessed (E_strong deferred)",
        "R5_Norman_only": "not assessed (Norman deferred)"}
    out = {"protocol_version": "0.5.1", "prep": prep_s, "glmm_convergence": {k: v.get("messages") for k, v in g.items()},
           **res, "criteria": crit}
    (AN / "g1_results_v0.5.1.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf8")
    print(json.dumps({"classification": crit["classification"], "holm": crit["holm_family"], "C4": crit["C4_linear_go_minus_mean_slope"],
                      "per_model": per_model, "reinterpretation": crit["reinterpretation"]}, default=float))


if __name__ == "__main__":
    {"prep": prep, "analyze": analyze}[sys.argv[1]]()
