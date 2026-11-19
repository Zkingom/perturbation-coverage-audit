# ext_v07dev.py -- additional difficulty baselines for the v0.7 development rerun (added 2026-10-10, AFTER the
# D1 results of v0.6.2 were seen; user decision 2026-10-10).
#   (a) MD-GO / MD-STR : Mondrian conformal on quintiles of the raw prior distance D_GO / D_STRING of a unit to the
#       model-fitting pool. Edges = quintiles of the out-of-fold distances of the learner units (learner unit in
#       inner fold f -> distance to the learner units outside f plus all always-train units), i.e. the same
#       "edges from training-pool OOF values" rule as M2 / Mdist-M. Prior-only, deployable.
#   (b) Munc-N / Munc-M : "PRESCRIBE-like" model-uncertainty difficulty proxy (we do NOT reproduce PRESCRIBE).
#       B bootstrap refits of linear-PCA (TSS groups of the model-fitting pool resampled with replacement; the
#       split's selected K; same recipe as eval_v06.refit_predict). For a query unit q, with standardised
#       predictions z_b(q) = (yhat_b(q) - mu) / sigma:
#         unc(q) = sqrt( mean over the 50 genes with largest |mean_b z_b(q)| (target column(s) excluded)
#                        of the across-bootstrap variance (ddof = 1) of z_b(q, g) ).
#       Cal / test units use all B refits; learner units use their out-of-bag refits only (OOB value; used
#       only for the normalisation constant c and the Mondrian edges). Never uses any response of a
#       cal / test unit.
# Output per split: results/<run>/ext_cache/ext_split_NN.npz (calt = cal followed by test, as in the cache).
import numpy as np

import eval_v06 as ev
import pipeline_v05 as pl5
import proxy_v06 as px
from paths import MASTER_SEED

B_BOOT = 50
TOP_K = 50
EXT_SEED = MASTER_SEED + 10000
MIN_OOB = 3


def oof_distances(pri, lrn, at, fold):
    out = {"D_GO": np.empty(len(lrn)), "D_STRING": np.empty(len(lrn))}
    for f in np.unique(fold):
        te = fold == f
        d = pl5.distances(pri, lrn[te], np.concatenate([lrn[~te], at]))
        out["D_GO"][te] = d["D_GO"]
        out["D_STRING"][te] = d["D_STRING"]
    return out


def _unc_from_moments(s1, s2, n, tcols):
    """Running sums s1 = sum_b z_b, s2 = sum_b z_b^2 over n refits (per unit x gene) -> dispersion on top-K genes
    by |mean| (target column(s) excluded); variance with ddof = 1."""
    m = s1 / n[:, None]
    v = (s2 - s1 * m) / (n[:, None] - 1)
    v = np.maximum(v, 0.0)
    a = np.abs(m)
    for i, c in enumerate(tcols):
        a[i, c] = -np.inf
    top = np.argsort(-a, axis=1)[:, :TOP_K]
    return np.sqrt(np.take_along_axis(v, top, axis=1).mean(1))


def compute_ext(cache, D, log=print):
    U, Y, y_ctrl, pri = D["U"], D["Y"], D["y_ctrl"], D["pri"]
    i = cache["split_id"]
    lrn, at, pool = np.asarray(cache["lrn"]), np.asarray(cache["at"]), np.asarray(cache["pool"])
    calt = np.concatenate([cache["cal"], cache["test"]])
    assert (U.target.values[calt] == cache["target_calt"]).all()
    fold = px.inner_folds(U.tss_group.values[lrn], ev._seed(MASTER_SEED, 6000, i))
    dc = pl5.distances(pri, calt, pool)
    assert np.allclose(dc["D_GO"], cache["covariates"].D_GO.values) and \
        np.allclose(dc["D_STRING"], cache["covariates"].D_STRING.values), "distance recomputation differs from cache"
    do = oof_distances(pri, lrn, at, fold)
    # (b) bootstrap refits of linear-PCA
    K = cache["chosen"]["K"]
    ref = {"mu": cache["ref_mu"], "sigma": cache["ref_sigma"]}
    groups = U.tss_group.values[pool]
    ug = np.unique(groups)
    members = {g: np.where(groups == g)[0] for g in ug}
    rng = np.random.default_rng([EXT_SEED, i])
    query = np.concatenate([calt, lrn])
    nq = len(calt)
    G = Y.shape[1]
    c1, c2 = np.zeros((nq, G)), np.zeros((nq, G))
    l1, l2 = np.zeros((len(lrn), G)), np.zeros((len(lrn), G))
    n_oob = np.zeros(len(lrn), int)
    for b in range(B_BOOT):
        pick = rng.choice(ug, len(ug), replace=True)
        tr = pool[np.concatenate([members[g] for g in pick])]
        oob = ~np.isin(lrn, tr)
        z = (ev.refit_predict(U, Y, y_ctrl, pri, tr, query, "pca", K) - ref["mu"]) / ref["sigma"]
        c1 += z[:nq]
        c2 += z[:nq] ** 2
        zl = z[nq:][oob]
        l1[oob] += zl
        l2[oob] += zl ** 2
        n_oob += oob
    tc_calt = [U.tcols.values[u] for u in calt]
    tc_lrn = [U.tcols.values[u] for u in lrn]
    unc_calt = _unc_from_moments(c1, c2, np.full(nq, B_BOOT), tc_calt)
    ok = n_oob >= MIN_OOB
    unc_oof = np.full(len(lrn), np.nan)
    unc_oof[ok] = _unc_from_moments(l1[ok], l2[ok], n_oob[ok], [t for t, o in zip(tc_lrn, ok) if o])
    out = {"dgo_calt": dc["D_GO"], "dstr_calt": dc["D_STRING"], "dgo_oof": do["D_GO"], "dstr_oof": do["D_STRING"],
           "unc_calt": unc_calt, "unc_oof": unc_oof, "n_oob": n_oob, "B": np.array(B_BOOT), "K": np.array(K)}
    log(f"ext split {i:02d}: K={K} unc_calt median {np.median(unc_calt):.4f} oof median {np.nanmedian(unc_oof):.4f} "
        f"n_oob min {n_oob.min()} (<{MIN_OOB}: {(n_oob < MIN_OOB).sum()})")
    return out


def attach(cache, ext):
    """Adds the extension arrays to a stage-1 cache in memory (the cache file itself is never modified)."""
    cache["ext"] = {k: np.asarray(v) for k, v in ext.items()}
    return cache
