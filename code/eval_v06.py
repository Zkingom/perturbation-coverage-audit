# eval_v06.py -- protocol v0.6.1 per-split computation.
# Stage 1, compute_split() (v0.6.2 learner-unit rule): train models on the training pool (frozen v0.5 code: pipeline_v05.train_objects /
#   predict / select_params / leakage_asserts); E_true and tiers; out-of-fold refits of the linear models;
#   deployable proxies h (all configurations), d and s (per model x k, both learners); calibration / test
#   scores; analysis-only covariates. No calibration yet.
# Stage 2, apply_methods(): conformal methods for one model x k x variant, nominal alpha or width-matched.
# A dataset D: U, Y, YA, YB (split-half pseudobulk), y_ctrl, pri, E_idx, AT_idx.
import numpy as np
import pandas as pd

import pipeline_v05 as pl
import proxy_v06 as px
from methods_v06 import (METHODS, M1_GAMMAS, M2_BINS, M3_PHI, capped_width, match_width, permute_within,
                         thresholds)
from models_v04 import fit_linear, gene_pca_embedding, score_S1
from paths import MASTER_SEED

MODELS = ("mean", "linear_pca", "linear_go")
LINEAR = ("linear_pca", "linear_go")
PRIMARY_MODEL = "linear_pca"
SECONDARY_MODEL = "linear_go"
KS = (20, 50, 200)
PRIMARY_K = 50
ALPHA = 0.10
N_TIERS = 5
DEV_SEED_BASE = MASTER_SEED + 1000
CONF_SEED_BASE = MASTER_SEED + 2000
HALF_SEED = MASTER_SEED + 3000          # used by pseudobulk_v06 when the halves are generated
PERM_SEED_BASE = MASTER_SEED + 4000
N_PERM = 20


def _seed(*words):
    return int(np.random.SeedSequence([int(w) for w in words]).generate_state(1)[0])


def refit_predict(U, Y, y_ctrl, pri, train_units, query, kind, param):
    """Refit one linear model (same recipe as pipeline_v05.train_objects) on train_units; predict query."""
    G = gene_pca_embedding(np.vstack([Y[train_units], y_ctrl[None, :]]), K=param if kind == "pca" else pl.K_GO,
                           method="gram")
    if kind == "pca":
        GK = G[:, :param]
        fit = train_units[np.array([len(U.tcols.values[i]) > 0 for i in train_units])]
        P = np.vstack([GK[U.tcols.values[i][0]] for i in fit] + [np.zeros(param)])
        f = fit_linear(np.vstack([Y[fit], y_ctrl]), P, GK, pl.LAM)
        return f(np.vstack([GK[U.tcols.values[i][0]] for i in query]))
    Gg = G[:, :pl.K_GO]
    fit = train_units[np.array([U.target.values[i] in pri["ase_ix"] for i in train_units])]
    P = np.vstack([pri["ase"][pri["ase_ix"][U.target.values[i]], :param] for i in fit] + [np.zeros(param)])
    f = fit_linear(np.vstack([Y[fit], y_ctrl]), P, Gg, pl.LAM)
    return f(np.vstack([pri["ase"][pri["ase_ix"][U.target.values[i]], :param] for i in query]))


def tiers_of(E_values, n):
    """Quantile tiers (1..n) of E over all E units of the dataset in this split."""
    return pd.qcut(pd.Series(E_values).rank(method="first"), n, labels=False).values + 1


def compute_split(split_id, sp, D, h_configs=None):
    U, Y, y_ctrl, pri, E_idx, AT_idx = (D[k] for k in ("U", "Y", "y_ctrl", "pri", "E_idx", "AT_idx"))
    assert "YA" in D and "YB" in D, "v0.6.1 needs split-half pseudobulk (nbq_noise, rel, noise, half variant)"
    YA, YB = D["YA"], D["YB"]
    folds = {k: E_idx[sp[k]] for k in ("train", "val", "cal", "test")}
    pool = np.concatenate([folds["train"], AT_idx])
    ref, G, models, fit_cols, ref_rows = pl.train_objects(U, Y, y_ctrl, pri, pool, "gram")
    ev = np.concatenate([folds["val"], folds["cal"], folds["test"]])
    nv, nc = len(folds["val"]), len(folds["cal"])
    pr = pl.predict(U, models, pri, G, ev)
    chosen, val_grid = pl.select_params(U, Y, ref, {m: p[:nv] for m, p in pr.items()}, folds["val"])
    pl.leakage_asserts(U, E_idx, folds, pool, fit_cols, ref_rows)
    key = {"mean": "mean", "linear_pca": f"linear_pca_K{chosen['K']}", "linear_go": f"linear_go_d{chosen['d']}"}
    spec = {"linear_pca": ("pca", chosen["K"]), "linear_go": ("go", chosen["d"])}
    tc = lambda idx: [U.tcols.values[i] for i in idx]                     # noqa: E731
    allE = np.concatenate([folds[k] for k in ("train", "val", "cal", "test")])
    E_all = pd.Series(px.e_true(Y[allE], ref, tc(allE)), index=allE)
    EA_all = pd.Series(px.e_true(YA[allE], ref, tc(allE)), index=allE)
    # v0.6.2: proxy learners (h, d, s) are trained on the E training units of this split only ("lrn").
    # Always-train (AT) units stay in model fitting and in the reference (train_objects above) and serve as
    # neighbours for the neighbour features. E_true has one definition for all units (mean-model S1, k = 50,
    # target-gene column(s) excluded whenever the target is a readout gene).
    lrn = folds["train"]
    at = np.asarray(AT_idx)
    assert np.array_equal(pool, np.concatenate([lrn, at])), "pool must be E-train followed by AT units"
    E_lrn = px.e_true(Y[lrn], ref, tc(lrn))
    E_at = px.e_true(Y[at], ref, tc(at))
    assert np.isfinite(E_lrn).all() and np.isfinite(E_at).all() and (E_at > 0).all(), "E_true not computable"
    calt = ev[nv:]
    stat = px.static_features(U, pri, y_ctrl)
    groups = U.tss_group.values[lrn]
    fold = px.inner_folds(groups, _seed(MASTER_SEED, 6000, split_id))
    _, noise_pool = px.split_half_stats(YA[pool], YB[pool], ref, tc(pool))
    lc = np.log(U.ncells.values[pool].astype(float))
    lnz = np.log(noise_pool)
    nl = len(lrn)
    nbq = {"lc_lrn": lc[:nl], "lc_at": lc[nl:], "ln_lrn": lnz[:nl], "ln_at": lnz[nl:]}
    # out-of-fold refits of the linear models on (other learner folds + all AT units) -> f_model OOF and
    # OOF scores of learner units (all E units, hence predictable by both linear models)
    oof_score, size_oof, size_full = {}, {"mean": None}, {"mean": None}
    mu_tile = np.tile(ref["mu"], (nl, 1))
    for k in KS:
        oof_score[("mean", k)] = score_S1(Y[lrn], mu_tile, ref, tc(lrn), k)
    for m in LINEAR:
        kind, par = spec[m]
        Yo = np.empty((nl, Y.shape[1]))
        for f in range(px.N_INNER):
            te = fold == f
            Yo[te] = refit_predict(U, Y, y_ctrl, pri, np.concatenate([lrn[~te], at]), lrn[te], kind, par)
        size_oof[m] = px.predicted_size(Yo, ref, tc(lrn))
        size_full[m] = px.predicted_size(pr[key[m]], ref, tc(ev))
        for k in KS:
            oof_score[(m, k)] = score_S1(Y[lrn], Yo, ref, tc(lrn), k)
        del Yo
    cache = {"split_id": split_id, "round": sp.get("round"), "fold": sp.get("fold"), "chosen": chosen,
             "val_grid": val_grid, "pool": pool, "cal": folds["cal"], "test": folds["test"], "val": folds["val"],
             "nv": nv, "lrn": lrn, "at": at, "E_pool": E_lrn, "E_at": E_at, "oof_score": oof_score,
             "E_calt": E_all.loc[calt].values, "EA_calt": EA_all.loc[calt].values,
             "tier_calt": pd.Series(tiers_of(E_all.values, 5), index=allE).loc[calt].values,
             "tier3_calt": pd.Series(tiers_of(E_all.values, 3), index=allE).loc[calt].values,
             "tierA_calt": pd.Series(tiers_of(EA_all.values, 5), index=allE).loc[calt].values,
             "target_calt": U.target.values[calt], "tss_calt": U.tss_group.values[calt],
             "ref_mu": ref["mu"], "ref_sigma": ref["sigma"], "scores": {}, "h": {}, "d": {}, "s": {}}
    for m in MODELS:
        Yp = pr[key[m]][nv:]
        for k in KS:
            cache["scores"][(m, k, "full")] = score_S1(Y[calt], Yp, ref, tc(calt), k)
        cache["scores"][(m, PRIMARY_K, "half")] = score_S1(YB[calt], Yp, ref, tc(calt), PRIMARY_K)
    # deployable proxies
    sd = _seed(MASTER_SEED, 7000, split_id)
    h_cfgs = h_configs or [(fs, ln) for fs in px.FEATURE_SETS for ln in px.LEARNERS]
    for m in MODELS:
        F_oof, F_q = px.feature_tables(lrn, fold, at, ev, E_lrn, E_at, pri, stat, size_oof[m], size_full[m], nbq)
        done = {}
        for fs, ln in h_cfgs:
            fs_m = fs if m != "mean" else px.without_model(fs)
            if (fs_m, ln) not in done:
                q, oof, lam = px.fit_proxy(F_oof, F_q, px.FEATURE_SETS[fs_m], np.log(E_lrn), fold, groups, ln, sd)
                done[(fs_m, ln)] = {"q": q, "oof": oof, "lambda": lam}
            cache["h"][(m, fs, ln)] = done[(fs_m, ln)]
        full_set = px.FULL_SET if m != "mean" else px.without_model(px.FULL_SET)
        for k in KS:
            for ln in px.LEARNERS:
                t = np.log(oof_score[(m, k)])
                q, oof, lam = px.fit_proxy(F_oof, F_q, px.FEATURE_SETS[px.D_SET], t, fold, groups, ln, sd + 1)
                cache["d"][(m, k, ln)] = {"q": q, "oof": oof, "lambda": lam}
                q, oof, lam = px.fit_proxy(F_oof, F_q, px.FEATURE_SETS[full_set], t, fold, groups, ln, sd + 2)
                cache["s"][(m, k, ln)] = {"q": q, "oof": oof, "lambda": lam}
        if m == "mean":
            cache["sim_calt"] = {c: F_q[c][nv:] for c in ("D_GO", "D_STRING", "sim_maxJ", "f_base")}
    # analysis-only covariates of cal + test units (never used by any deployable method)
    rel, noise = px.split_half_stats(YA[calt], YB[calt], ref, tc(calt))
    cache["covariates"] = pd.DataFrame({
        "target": U.target.values[calt], "tss_group": U.tss_group.values[calt],
        "ncells": U.ncells.values[calt], "rel": rel, "noise": noise,
        "shared": px.shared_strength(Y[calt], y_ctrl, ref), "E_true": cache["E_calt"], "E_true_A": cache["EA_calt"],
        **cache["sim_calt"],
        **{f"pred_size_{m}": size_full[m][nv:] for m in LINEAR}})
    return cache


# ------------------------------------------------------------------ stage 2

def method_grid():
    g = [("M0", None)]
    g += [("M1", x) for x in M1_GAMMAS] + [("M2", b) for b in M2_BINS] + [("M3", p) for p in M3_PHI]
    g += [("Mdist-N", x) for x in M1_GAMMAS] + [("Mdist-M", b) for b in M2_BINS]
    g += [("Mscore-N", x) for x in M1_GAMMAS] + [("Mcomb-N", x) for x in M1_GAMMAS] + [("Mcomb-3", "lin2")]
    g += [("oracle", x) for x in M1_GAMMAS]
    # v0.7 development additions (used only when the split carries ext arrays, see proxies / apply_methods)
    g += [("MD-GO", 5), ("MD-STR", 5)] + [("Munc-N", x) for x in M1_GAMMAS] + [("Munc-M", b) for b in M2_BINS]
    return g


def proxies(cache, model, k, sel):
    """Proxy dicts (cal, test, oof) for one model x k. sel = {"h": (fs, ln), "d": ln, "s": ln}."""
    nv, nc = cache["nv"], len(cache["cal"])
    H = cache["h"][(model, *sel["h"])]
    Dd = cache["d"][(model, k, sel["d"])]
    Ss = cache["s"][(model, k, sel["s"])]
    y = np.log(cache["oof_score"][(model, k)])
    Xo = np.column_stack([np.ones(len(y)), np.log(H["oof"]), np.log(Dd["oof"])])
    beta = np.linalg.lstsq(Xo, y, rcond=None)[0]
    comb_q = np.exp(np.column_stack([np.ones(len(H["q"])), np.log(H["q"]), np.log(Dd["q"])]) @ beta)
    calt = {"h": H["q"][nv:], "d": Dd["q"][nv:], "s": Ss["q"][nv:], "comb": comb_q[nv:], "E": cache["E_calt"]}
    P_cal = {v: a[:nc] for v, a in calt.items()}
    P_test = {v: a[nc:] for v, a in calt.items()}
    P_oof = {"h": H["oof"], "d": Dd["oof"], "s": Ss["oof"], "comb": np.exp(Xo @ beta), "E": cache["E_pool"]}
    ext = cache.get("ext")
    if ext is not None:                      # v0.7 development additions (ext_v07dev.py)
        keys = [("dgo", "dgo"), ("dstr", "dstr")] + ([("unc", "unc")] if model == "linear_pca" else [])
        for v, src in keys:
            P_cal[v], P_test[v] = ext[f"{src}_calt"][:nc], ext[f"{src}_calt"][nc:]
            oo = np.asarray(ext[f"{src}_oof"], float)
            P_oof[v] = oo[np.isfinite(oo)]
    return P_cal, P_test, P_oof, beta


def apply_methods(cache, model, k, sel, grid=None, variant="full", mode="nominal", perm_r=None, alpha=ALPHA):
    """Coverage rows for one model x k x variant. mode: "nominal" (alpha = 0.1), "matched" (alpha' chosen on
    test-unit proxies to match M0's capped mean width) or "matched_cal" (alpha' chosen on calibration-unit
    proxies; sensitivity). perm_r: None for the real h, r = 0..N_PERM-1 for the permuted-h control."""
    nc = len(cache["cal"])
    s = cache["scores"][(model, k, variant)]
    s_cal, s_test = s[:nc], s[nc:]
    P_cal, P_test, P_oof, _ = proxies(cache, model, k, sel)
    if perm_r is not None:
        P_cal["h"], P_test["h"] = permute_within(P_cal["h"], P_test["h"],
                                                 _seed(PERM_SEED_BASE + perm_r, cache["split_id"]))
    wi = MODELS.index(model)
    rng = np.random.default_rng(_seed(MASTER_SEED, 5000, cache["split_id"], wi, k, 0 if variant == "full" else 1,
                                      0 if perm_r is None else perm_r + 1))
    u_test, u_cal = rng.uniform(size=len(s_test)), rng.uniform(size=nc)
    tier = cache["tierA_calt"] if variant == "half" else cache["tier_calt"]
    E = cache["EA_calt"] if variant == "half" else cache["E_calt"]
    W0 = capped_width(thresholds("M0", None, s_cal, P_cal, P_test, P_oof, u_test, alpha), s_cal)
    rows = []
    avail = lambda nm: METHODS[nm][1] is None or isinstance(METHODS[nm][1], tuple) or METHODS[nm][1] in P_cal  # noqa: E731
    for name, hy in [g_ for g_ in (grid or method_grid()) if avail(g_[0])]:
        if mode == "nominal":
            a_used, resid, mono = alpha, np.nan, True
            t = thresholds(name, hy, s_cal, P_cal, P_test, P_oof, u_test, alpha)
        elif mode == "matched":
            a_used, t, resid, mono = match_width(name, hy, s_cal, P_cal, P_test, P_oof, u_test, W0)
        else:
            a_used, _, resid, mono = match_width(name, hy, s_cal, P_cal, P_cal, P_oof, u_cal, W0)
            t = thresholds(name, hy, s_cal, P_cal, P_test, P_oof, u_test, a_used)
        cap = 2.0 * float(np.max(s_cal))
        rows.append(pd.DataFrame({
            "split": cache["split_id"], "model": model, "k": k, "variant": variant, "mode": mode, "method": name,
            "hyper": str(hy), "perm": -1 if perm_r is None else perm_r, "alpha_used": a_used, "width_resid": resid,
            "monotone": mono, "target": cache["target_calt"][nc:], "tss_group": cache["tss_calt"][nc:],
            "tier": tier[nc:], "tier3": cache["tier3_calt"][nc:], "E_true": E[nc:], "h": P_test["h"],
            "d": P_test["d"], "D_GO": cache["sim_calt"]["D_GO"][nc:], "score": s_test, "thr": t,
            "w": np.clip(np.where(np.isposinf(t), cap, t), 0.0, cap), "covered": s_test <= t}))
    return pd.concat(rows, ignore_index=True)


def proxy_table(cache):
    """Test-unit proxy values (h for every configuration; d and s for every learner) for rules 7.1 / 7.1b."""
    nv, nc = cache["nv"], len(cache["cal"])
    base = {"split": cache["split_id"], "target": cache["target_calt"][nc:], "tss_group": cache["tss_calt"][nc:],
            "tier": cache["tier_calt"][nc:], "E_true": cache["E_calt"][nc:]}
    out = []
    for (m, fs, ln), v in cache["h"].items():
        out.append(pd.DataFrame({**base, "model": m, "kind": "h", "k": 0, "proxy": f"{fs}|{ln}",
                                 "value": v["q"][nv + nc:], "score": np.nan, "lambda": v["lambda"]}))
    for kind in ("d", "s"):
        for (m, k, ln), v in cache[kind].items():
            out.append(pd.DataFrame({**base, "model": m, "kind": kind, "k": k, "proxy": ln,
                                     "value": v["q"][nv + nc:], "score": cache["scores"][(m, k, "full")][nc:],
                                     "lambda": v["lambda"]}))
    return pd.concat(out, ignore_index=True)
