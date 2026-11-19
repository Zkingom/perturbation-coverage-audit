# p1_v04.py -- protocol v0.4.2, phase P1/G1 runner: mean, linear-PCA (selected K and K=10),
# linear-GO (selected d and d=2),
# 3 random TSS-group splits (seeds master_seed + 0..2), all leakage asserts, canaries and the
# three self-checks. Writes read-only folds to splits/v0.4 and results to results/p1_v0.4.
import json
import os
import pickle
import stat
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED, RESULTS, SPLITS, raw  # noqa: E402
from audit_lib import betabinom_band, make_bins, random_group_split  # noqa: E402
from models_v04 import (d_go, d_string, fit_linear, fit_mean, gene_pca_embedding, go_incidence,  # noqa: E402
                        jaccard_all, score_S1, score_S2, train_reference)

OUT = RESULTS / "p1_v0.4"
SPLIT_DIR = SPLITS / "v0.4"
N_SPLITS = 3
LAM = 0.1
K_CANDIDATES = (5, 10, 20, 30, 50)        # v0.4.2: linear-PCA K selected on val
K_GO = 10                                   # gene embedding G of linear-GO (unchanged)
K_MAX = max(K_CANDIDATES)
D_CANDIDATES = (2, 5, 10, 20, 30, 50, 100)  # v0.4.2
K_REPRO, D_REPRO = 10, 2
KS = (20, 50, 200)
ALPHAS = (0.05, 0.10, 0.20)
MIN_CELLS_REF = 30
PLANT_TAU = 6.0          # fixed from the v0.2 P1 positive control (chosen post hoc there; disclosed)
GO_TERM_MAX = 500        # v0.3 item 5: drop GO terms annotating > 500 genes of the GEARS universe
PRIMARY = ("S1_k50", 0.10)
PCA_TOL_G = 1e-8          # relative to max|G_svd| (protocol v0.4 (f))
PCA_TOL_PRED = 1e-8       # absolute, log1p(CP10k) units


# ------------------------------------------------------------------ data
def load():
    U = pd.read_csv(PROCESSED / "v0.3_figshare_k562ess" / "p0v3_units_all_targets.csv")
    pb = np.load(PROCESSED / "v0.4" / "pseudobulk_logmean.npy").astype(np.float64)
    rows = pd.read_csv(PROCESSED / "v0.4" / "pseudobulk_rows.csv")
    genes = pd.read_csv(PROCESSED / "v0.4" / "pseudobulk_genes.csv").gene_name.values
    rix = {g: i for i, g in enumerate(rows.gene_transcript)}
    assert U.gene_transcript.isin(rix).all(), "unit rows missing from pseudobulk"
    ctrl_rows = np.where(rows.gene_transcript.str.contains("non-targeting", case=False))[0]
    w = rows.n_cells.values[ctrl_rows].astype(np.float64)
    assert len(ctrl_rows) > 0 and np.isfinite(w).all() and (w > 0).all(), "control cell counts invalid"
    y_ctrl = (pb[ctrl_rows] * w[:, None]).sum(0) / w.sum()          # mean over all control cells
    U = U[U.E | U.always_train].reset_index(drop=True)                # dropped / secondary rows never used
    U["row"] = U.gene_transcript.map(rix)
    U["n_cells_sc"] = rows.n_cells.values[U.row.values]
    gcols = {}
    for j, g in enumerate(genes):
        gcols.setdefault(g, []).append(j)
    U["tcols"] = [np.array(gcols.get(t, []), dtype=int) for t in U.target]
    Y = pb[U.row.values]
    with open(raw("gene2go_all.pkl"), "rb") as fh:
        gene2go = pickle.load(fh)
    with open(raw("essential_all_data_pert_genes.pkl"), "rb") as fh:
        ess = pickle.load(fh)
    universe = [g for g in ess if g in gene2go]
    cnt = Counter(t for g in set(universe) for t in gene2go[g])
    broad = frozenset(t for t, c in cnt.items() if c > GO_TERM_MAX)
    ase = np.load(PROCESSED / "priors" / "go_ase100.npz", allow_pickle=True)
    ase_ix = {g: i for i, g in enumerate(ase["genes"])}
    st = pd.read_parquet(PROCESSED / "priors" / "string_expdb_v0.4_universe.parquet")
    pos = {g: i for i, g in enumerate(U.target)}
    S = np.zeros((len(U), len(U)))
    st = st[st.g1.isin(pos) & st.g2.isin(pos)]
    for a, b, s in st.itertuples(index=False):
        i, j = pos[a], pos[b]
        S[i, j] = S[j, i] = max(S[i, j], s)
    np.fill_diagonal(S, 0.0)
    pri = {"J": jaccard_all(go_incidence(U.target.values, gene2go)),
           "Jf": jaccard_all(go_incidence(U.target.values, gene2go, broad)),
           "S": S, "ase": ase["ase"], "ase_ix": ase_ix,
           "n_universe": len(set(universe)), "n_broad": len(broad)}
    return U, Y, y_ctrl, genes, pri


# ------------------------------------------------------------------ one split: train-derived objects
def gene_index(U, idx):
    return [U.tcols.values[i] for i in idx]


def train_objects(U, Y, y_ctrl, pri, pool, pca_method="gram"):
    """Everything here may use expression ONLY of rows in `pool` (training pool) plus controls."""
    ref_rows = pool[U.ncells.values[pool] >= MIN_CELLS_REF]
    ref = train_reference(Y[ref_rows])
    # PCA scores are nested: the first K columns of the K_MAX solution are the K-dimensional solution.
    G = gene_pca_embedding(np.vstack([Y[pool], y_ctrl[None, :]]), K=K_MAX, method=pca_method)
    # v0.4.1: the mean model predicts exactly mu_tr (same >=30-cell training rows as the reference),
    # so its prediction in the specific space s = y - mu_tr is identically 0.
    models = {"mean": ("mean", fit_mean(Y[ref_rows]))}
    inR = pool[np.array([len(U.tcols.values[i]) > 0 for i in pool])]
    for K in K_CANDIDATES:
        GK = G[:, :K]
        Ppca = np.vstack([GK[U.tcols.values[i][0]] for i in inR] + [np.zeros(K)])
        models[f"linear_pca_K{K}"] = (("pca", K), fit_linear(np.vstack([Y[inR], y_ctrl]), Ppca, GK, LAM))
    inGO = pool[np.array([U.target.values[i] in pri["ase_ix"] for i in pool])]
    Ggo = G[:, :K_GO]
    for d in D_CANDIDATES:
        Pgo = np.vstack([pri["ase"][pri["ase_ix"][U.target.values[i]], :d] for i in inGO] + [np.zeros(d)])
        models[f"linear_go_d{d}"] = (("go", d), fit_linear(np.vstack([Y[inGO], y_ctrl]), Pgo, Ggo, LAM))
    fit_cols = {"linear_pca": inR, "linear_go": inGO}
    return ref, G, models, fit_cols, ref_rows


def predict(U, models, pri, G, idx):
    out = {}
    for m, (kind, f) in models.items():
        if kind == "mean":
            out[m] = f(len(idx))
        elif kind[0] == "pca":
            out[m] = f(np.vstack([G[U.tcols.values[i][0], :kind[1]] for i in idx]))
        else:
            d = kind[1]
            out[m] = f(np.vstack([pri["ase"][pri["ase_ix"][U.target.values[i]], :d] for i in idx]))
    return out


def _pick(err, grid):
    best = min(err.values())
    return min(c for c in grid if err[c] - best <= 1e-12 * max(1.0, abs(best)))   # ties -> smaller value


def select_params(U, Y, ref, preds_val, val):
    """v0.4.2: K (linear-PCA) and d (linear-GO) chosen by mean val S1(k=50); no refit on val."""
    tc = gene_index(U, val)
    errK = {K: float(score_S1(Y[val], preds_val[f"linear_pca_K{K}"], ref, tc, 50).mean()) for K in K_CANDIDATES}
    errd = {d: float(score_S1(Y[val], preds_val[f"linear_go_d{d}"], ref, tc, 50).mean()) for d in D_CANDIDATES}
    return {"K": _pick(errK, K_CANDIDATES), "d": _pick(errd, D_CANDIDATES)}, {"K": errK, "d": errd}


def distances(pri, ev, pool):
    return {"D_GO": d_go(pri["J"][np.ix_(ev, pool)]), "D_GO_filt": d_go(pri["Jf"][np.ix_(ev, pool)]),
            "D_STRING": d_string(pri["S"][np.ix_(ev, pool)])}


# ------------------------------------------------------------------ leakage checks
def leakage_asserts(U, E_idx, folds, pool, fit_cols, ref_rows):
    held = np.concatenate([folds["val"], folds["cal"], folds["test"]])
    allE = np.concatenate(list(folds.values()))
    assert len(allE) == len(E_idx) and len(np.unique(allE)) == len(E_idx), "L1: folds are not a partition of E"
    assert U.target.is_unique, "L1: one row per target gene"
    fold_of = pd.Series(index=allE, data=np.concatenate([[k] * len(v) for k, v in folds.items()]))
    gfold = pd.DataFrame({"g": U.tss_group.values[allE], "f": fold_of.loc[allE].values}).groupby("g").f.nunique()
    assert (gfold == 1).all(), "L2: a TSS group spans several folds"
    assert not U.E.values[U.always_train.values].any(), "always-train overlaps E"
    assert set(pool).isdisjoint(held), "train pool contains held-out units"
    assert set(U.tss_group.values[pool]).isdisjoint(U.tss_group.values[held]), "L2: held-out TSS group in train pool"
    for k, cols in fit_cols.items():
        assert set(cols) <= set(pool), f"{k}: W fitted on non-training columns"
    assert set(ref_rows) <= set(pool) and (U.ncells.values[ref_rows] >= MIN_CELLS_REF).all(), "(d) reference rows"


def canary(U, Y, y_ctrl, pri, folds, pool, corrupt, base, seed, pca_method):
    """Replace expression of `corrupt` folds with noise and recompute train-derived objects."""
    rng = np.random.default_rng(seed + 999)
    Yc = Y.copy()
    rows = np.concatenate([folds[k] for k in corrupt])
    Yc[rows] = rng.normal(5.0, 3.0, size=Yc[rows].shape)
    ref2, G2, models2, _, _ = train_objects(U, Yc, y_ctrl, pri, pool, pca_method)
    ev = np.concatenate([folds["val"], folds["cal"], folds["test"]])
    p2 = predict(U, models2, pri, G2, ev)
    same = (all(np.allclose(base["ref"][k], ref2[k]) for k in ("mu", "sigma"))
            and np.array_equal(base["ref"]["panel"], ref2["panel"])
            and all(np.allclose(base["preds_all"][m], p2[m]) for m in p2))
    d2 = distances(pri, ev, pool)
    same = same and all(np.allclose(base["dist_all"][k], d2[k]) for k in d2)
    chosen2 = None
    if "val" not in corrupt:
        nv = len(folds["val"])
        chosen2, _ = select_params(U, Yc, ref2, {m: v[:nv] for m, v in p2.items()}, folds["val"])
    return same, chosen2


def write_readonly(folds, U, name):
    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    p = SPLIT_DIR / f"{name}.json"
    if p.exists():
        os.chmod(p, stat.S_IWRITE | stat.S_IREAD)
    p.write_text(json.dumps({k: U.target.values[v].tolist() for k, v in folds.items()}), encoding="utf8")
    os.chmod(p, stat.S_IREAD)
    try:
        open(p, "a").close()
        blocked = False
    except PermissionError:
        blocked = True
    assert blocked, "L10: fold file is writable"
    return str(p.relative_to(SPLITS.parent))


# ------------------------------------------------------------------ run
def run_split(i, U, Y, y_ctrl, pri, E_idx, AT_idx):
    seed = MASTER_SEED + i
    sub = random_group_split(U.tss_group.values[E_idx], seed)
    folds = {k: E_idx[v] for k, v in sub.items()}
    pool = np.concatenate([folds["train"], AT_idx])
    ev_all = np.concatenate([folds["val"], folds["cal"], folds["test"]])
    nv = len(folds["val"])
    # ---- PCA method check (protocol v0.4 (f)): Gram eigendecomposition vs SVD on real data
    obj = {}
    for meth in ("gram", "svd"):
        try:
            o = train_objects(U, Y, y_ctrl, pri, pool, meth)
        except np.linalg.LinAlgError:
            meth = "gesvd"
            o = train_objects(U, Y, y_ctrl, pri, pool, meth)
        pr = predict(U, o[2], pri, o[1], ev_all)
        cd, _ = select_params(U, Y, o[0], {m: p[:nv] for m, p in pr.items()}, folds["val"])
        obj[meth] = (o, pr, cd)
    svd_name = [k for k in obj if k != "gram"][0]
    (o_g, pr_g, cd_g), (o_s, pr_s, cd_s) = obj["gram"], obj[svd_name]
    Gg, Gs = o_g[1], o_s[1]
    sign = np.sign((Gg * Gs).sum(0))
    g_diff = float(np.abs(Gg * sign - Gs).max() / np.abs(Gs).max())
    p_diff = float(max(np.abs(pr_g[m] - pr_s[m]).max() for m in pr_g if m != "mean"))
    pca_ok = (g_diff <= PCA_TOL_G) and (p_diff <= PCA_TOL_PRED) and (cd_g == cd_s)
    use = "gram" if pca_ok else svd_name
    pca_check = {"svd_method": svd_name, "G_rel_maxdiff": g_diff, "pred_maxabsdiff": p_diff,
                 "chosen_gram": cd_g, "chosen_svd": cd_s, "pass": pca_ok, "used": use}
    (ref, G, models, fit_cols, ref_rows), preds_all, chosen = obj[use]
    leakage_asserts(U, E_idx, folds, pool, fit_cols, ref_rows)
    _, val_grid = select_params(U, Y, ref, {m: p[:nv] for m, p in preds_all.items()}, folds["val"])
    dist_all = distances(pri, ev_all, pool)
    assert np.allclose(preds_all["mean"][0], ref["mu"]), "v0.4.1: mean prediction must equal mu_tr"
    base = {"ref": ref, "preds_all": preds_all, "dist_all": dist_all}
    okA, chosenA = canary(U, Y, y_ctrl, pri, folds, pool, ("cal", "test"), base, seed, use)
    assert okA and chosenA == chosen, "canary A (cal+test corrupted) changed a train/val-derived quantity"
    okB, _ = canary(U, Y, y_ctrl, pri, folds, pool, ("val",), base, seed + 1, use)
    assert okB, "canary B (val corrupted) changed a train-derived quantity"
    fold_file = write_readonly(folds, U, f"A_{i:03d}")

    # validation error of every model (primary score S1_k50, val units only)
    tc_val = gene_index(U, folds["val"])
    val_err = {m: float(score_S1(Y[folds["val"]], p[:nv], ref, tc_val, 50).mean()) for m, p in preds_all.items()}

    # conformal evaluation on cal+test
    ev = ev_all[nv:]
    ncal = len(folds["cal"])
    is_cal = np.r_[np.ones(ncal, bool), np.zeros(len(folds["test"]), bool)]
    dist = {k: v[nv:] for k, v in dist_all.items()}
    eval_models = {"mean": "mean", "linear_pca": f"linear_pca_K{chosen['K']}", f"linear_pca_K{K_REPRO}_repro": f"linear_pca_K{K_REPRO}",
                   "linear_go": f"linear_go_d{chosen['d']}", f"linear_go_d{D_REPRO}_repro": f"linear_go_d{D_REPRO}"}
    Yt = Y[ev]
    tc = gene_index(U, ev)
    rng = np.random.default_rng(seed + 12345)
    u_test = pd.Series(dist["D_GO"][~is_cal]).rank(pct=True).values
    perm = rng.permutation(np.where(~is_cal)[0])
    rows = []
    for variant in ("real", "planted", "permuted"):
        for lab, key in eval_models.items():
            Yp = preds_all[key][nv:].copy()
            if variant == "planted":
                if lab != "linear_pca":
                    continue
                Yp[~is_cal] += PLANT_TAU * u_test[:, None] * ref["sigma"][None, :] * rng.standard_normal(Yp[~is_cal].shape)
            sc = {f"S1_k{k}": score_S1(Yt, Yp, ref, tc, k) for k in KS}
            sc["S2"] = score_S2(Yt, Yp, ref, tc)
            dG = dist["D_GO"].copy()
            if variant == "permuted":
                t = np.where(~is_cal)[0]
                dG[t] = dist["D_GO"][perm]
            for sname, s in sc.items():
                for a in ALPHAS:
                    n = len(s[is_cal])
                    r = int(np.ceil((n + 1) * (1 - a)))
                    q = np.inf if r > n else np.sort(s[is_cal])[r - 1]
                    for j in np.where(~is_cal)[0]:
                        rows.append((variant, i, lab, sname, a, U.target.values[ev[j]], dG[j], dist["D_GO_filt"][j],
                                     dist["D_STRING"][j], s[j], q, bool(s[j] <= q), ncal, int(U.ncells.values[ev[j]])))
    res = pd.DataFrame(rows, columns=["variant", "split", "model", "score", "alpha", "target", "D_GO", "D_GO_filt",
                                      "D_STRING", "score_value", "qhat", "covered", "n_cal", "n_cells"])
    meta = {"split": i, "seed": seed, "fold_sizes": {k: len(v) for k, v in folds.items()},
            "train_pool": int(len(pool)), "ref_rows": int(len(ref_rows)), "fit_cols_pca": int(len(fit_cols["linear_pca"])),
            "fit_cols_go": int(len(fit_cols["linear_go"])), "chosen": chosen, "chosen_at_grid_max": {"K": chosen["K"] == max(K_CANDIDATES), "d": chosen["d"] == max(D_CANDIDATES)}, "val_err_candidates": val_grid,
            "val_err": val_err, "pca_check": pca_check, "canary_A_ok": okA, "canary_B_ok": okB, "fold_file": fold_file}
    return res, meta


def trend(df, cov="D_GO"):
    d = df.copy()
    d["u"] = d.groupby("split")[cov].rank(pct=True)
    fit = sm.GLM((~d.covered).astype(float), sm.add_constant(d[["u"]]), family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": pd.factorize(d.target)[0]})
    return float(fit.params["u"]), float(fit.bse["u"]), float(fit.pvalues["u"])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    U, Y, y_ctrl, genes, pri = load()
    assert pri["n_universe"] == 9853 and pri["n_broad"] == 29, (pri["n_universe"], pri["n_broad"])
    E_idx = np.where(U.E.values)[0]
    AT_idx = np.where(U.always_train.values)[0]
    assert len(E_idx) == 1786 and U.tss_group.values[E_idx].size and len(AT_idx) == 257
    assert all(t in pri["ase_ix"] for t in U.target.values[E_idx]), "E target missing from GO graph"
    res, metas = [], []
    for i in range(N_SPLITS):
        r, m = run_split(i, U, Y, y_ctrl, pri, E_idx, AT_idx)
        res.append(r)
        metas.append(m)
        print(json.dumps({k: m[k] for k in ("split", "fold_sizes", "chosen", "pca_check")}), flush=True)
    R = pd.concat(res, ignore_index=True)
    R.to_parquet(OUT / "p1_v0.4_results.parquet")

    # ---- self-check (i): marginal coverage vs beta-binomial band, real variant
    real = R[R.variant == "real"]
    band = []
    for (sp_, m, s, a), g in real.groupby(["split", "model", "score", "alpha"]):
        lo, hi = betabinom_band(int(g.n_cal.iloc[0]), len(g), a)
        k = int(g.covered.sum())
        band.append((sp_, m, s, a, k, len(g), lo, hi, lo <= k <= hi))
    B = pd.DataFrame(band, columns=["split", "model", "score", "alpha", "n_covered", "n_test", "band_lo", "band_hi", "in_band"])
    B.to_csv(OUT / "selfcheck_marginal.csv", index=False)
    from scipy.stats import binom
    n_out = int((~B.in_band).sum())
    p_excess = float(binom.sf(n_out - 1, len(B), 0.05))
    prim = B[(B.score == PRIMARY[0]) & (B.alpha == PRIMARY[1])]
    check_i = p_excess > 0.01

    # ---- self-checks (ii) planted and (iii) permuted, primary score
    sel = lambda D: D[(D.score == PRIMARY[0]) & (D.alpha == PRIMARY[1])]
    pl = sel(R[R.variant == "planted"])
    b_pl, se_pl, p_pl = trend(pl)
    check_ii = (b_pl > 0) and (p_pl < 0.05)
    perm_tr = {m: trend(sel(R[(R.variant == "permuted") & (R.model == m)])) for m in sel(R).model.unique()}
    check_iii = all(v[2] >= 0.05 for v in perm_tr.values())

    # ---- descriptive (3 splits; no inference)
    rp = sel(real).copy()
    rp["bin_GO"] = rp.groupby("split").D_GO.transform(lambda x: pd.Series(make_bins(x.values), index=x.index))
    desc = rp.groupby(["model", "bin_GO"]).covered.mean().unstack().round(3)
    desc.to_csv(OUT / "descriptive_coverage_by_DGO_bin.csv")
    qh = real.groupby(["model", "score", "alpha"]).qhat.agg(["mean", "min", "max"]).reset_index()
    qh.to_csv(OUT / "qhat_by_model.csv", index=False)
    marg = real.groupby(["model", "score", "alpha"]).covered.mean().reset_index()
    marg.to_csv(OUT / "marginal_coverage.csv", index=False)

    val = pd.DataFrame([{"split": m["split"], "model": k, "val_err_S1k50": v} for m in metas for k, v in m["val_err"].items()])
    base_err = val[val.model == "mean"].set_index("split").val_err_S1k50
    val["rel_improvement_vs_mean"] = 1 - val.val_err_S1k50 / val.split.map(base_err)
    val.to_csv(OUT / "validation_error.csv", index=False)

    summary = {"protocol_version": "0.4.2", "phase": "P1", "splits": metas,
               "selfcheck_i": {"n_combos": len(B), "n_out_of_band": n_out, "p_excess_vs_5pct": p_excess,
                               "primary_out_of_band": prim[~prim.in_band][["split", "model"]].values.tolist(), "pass": check_i},
               "selfcheck_ii_planted": {"tau": PLANT_TAU, "model": "linear_pca", "slope": b_pl, "se": se_pl, "p": p_pl, "pass": check_ii},
               "selfcheck_iii_permuted": {m: {"slope": v[0], "se": v[1], "p": v[2]} for m, v in perm_tr.items()} | {"pass": check_iii},
               "n_universe_go": pri["n_universe"], "n_broad_terms_removed": pri["n_broad"],
               "ncells_agreement_bulk_vs_singlecell": float((U.ncells == U.n_cells_sc).mean())}
    (OUT / "p1_v0.4_summary.json").write_text(json.dumps(summary, indent=1, default=float), encoding="utf8")
    print(json.dumps({"selfcheck": [check_i, check_ii, check_iii], "planted": [b_pl, p_pl],
                      "permuted_p": {m: v[2] for m, v in perm_tr.items()}, "n_out_of_band": n_out, "n_combos": len(B)}, default=float))
    assert check_i and check_ii and check_iii, "a self-check failed"


if __name__ == "__main__":
    main()
