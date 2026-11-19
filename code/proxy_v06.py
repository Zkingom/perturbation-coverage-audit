# proxy_v06.py -- deployable proxies of protocol v0.6.1 (sections 2-4).
#   h_hat  predicted effect strength (target log E_true of training-pool units)
#   d_hat  similarity / training-data-quality difficulty proxy (target log out-of-fold score; features sim+nbq)
#   s_hat  all h-features fitted to the log out-of-fold score (M_score)
# Deployable inputs only: GO / STRING priors, training-pool units' measured responses, cell counts and
# split-half noise (as neighbour statistics), target-gene baseline expression in control cells, and the
# trained model's own prediction. Never: responses, cell counts or noise of the query unit itself,
# of any val / cal / test unit, knockdown efficiency, E_true of non-pool units.
#
# Cross-fitting (v0.6.2): the learner units (E training units of the split; always-train units are
# neighbours only) are split into N_INNER folds by TSS group (fold vector supplied by
# the caller). Out-of-fold (OOF) features of fold f use neighbours outside f and models refitted without
# f. Learners are fitted on the OOF features of the whole pool and predict query units from features
# computed against the full pool. OOF predictions of pool units come from learners fitted on 4 folds.
import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits

from models_v04 import score_S1
from pipeline_v05 import d_go, d_string

N_INNER = 5
TOP_GO = 20
RIDGE_LAMBDAS = (0.01, 0.1, 1.0, 10.0, 100.0)
GBM_PARAMS = dict(max_iter=200, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=20, l2_regularization=1.0)

SIM = ["D_GO", "D_GO_filt", "D_STRING", "sim_maxJ", "sim_top20J", "str_n", "str_max"]
PRIOR = ["f_go", "f_go_miss", "f_str", "f_str_miss", "f_ann", "f_deg"]
BASE = ["f_base", "f_base_miss"]
MODEL = ["f_model"]
NBQ = ["nbq_cells", "nbq_noise"]
FEATURE_SETS = {
    "sim": SIM, "sim_prior": SIM + PRIOR, "sim_prior_base": SIM + PRIOR + BASE,
    "sim_prior_base_model": SIM + PRIOR + BASE + MODEL,
    "sim+nbq": SIM + NBQ, "sim_prior+nbq": SIM + PRIOR + NBQ, "sim_prior_base+nbq": SIM + PRIOR + BASE + NBQ,
    "sim_prior_base_model+nbq": SIM + PRIOR + BASE + MODEL + NBQ,
}
FULL_SET = "sim_prior_base_model+nbq"
D_SET = "sim+nbq"
LEARNERS = ("ridge", "gbm")


def n_substantive(fs):
    return sum(1 for c in FEATURE_SETS[fs] if not c.endswith("_miss"))


def without_model(fs):
    return fs.replace("_model", "")


# ------------------------------------------------------------------ targets / unit statistics

def e_true(Y_rows, ref, tcols):
    """Mean-model S1(k=50); the mean-model prediction is exactly mu_tr (protocol v0.4.1)."""
    return score_S1(Y_rows, np.tile(ref["mu"], (len(Y_rows), 1)), ref, tcols, 50)


def predicted_size(Y_pred, ref, tcols, k=50):
    """sigma-scaled RMS of the predicted specific signal on its own predicted top-k genes."""
    sh = (Y_pred - ref["mu"]) / ref["sigma"]
    a = np.abs(sh).copy()
    for i, c in enumerate(tcols):
        a[i, c] = -np.inf
    top = np.argsort(-a, axis=1)[:, :k]
    return np.sqrt((np.take_along_axis(sh, top, axis=1) ** 2).mean(1))


def split_half_stats(YA, YB, ref, tcols):
    """rel = Pearson(s_A, s_B) over readout genes; noise = RMS((s_A - s_B) / 2); target columns excluded."""
    sA = (YA - ref["mu"]) / ref["sigma"]
    sB = (YB - ref["mu"]) / ref["sigma"]
    rel, noise = np.empty(len(YA)), np.empty(len(YA))
    for i, c in enumerate(tcols):
        m = np.ones(YA.shape[1], bool)
        m[c] = False
        a, b = sA[i, m], sB[i, m]
        rel[i] = np.corrcoef(a, b)[0, 1]
        noise[i] = np.sqrt((((a - b) / 2) ** 2).mean())
    return rel, noise


def shared_strength(Y_rows, y_ctrl, ref):
    """Projection of sigma-scaled (y - y_ctrl) on the unit vector of sigma-scaled (mu_tr - y_ctrl)."""
    v = (ref["mu"] - y_ctrl) / ref["sigma"]
    v = v / (np.linalg.norm(v) + 1e-12)
    return ((Y_rows - y_ctrl) / ref["sigma"]) @ v


# ------------------------------------------------------------------ features

def static_features(U, pri, y_ctrl):
    """Split-independent prior features: f_ann, f_deg, f_base (control cells only)."""
    base = np.array([y_ctrl[c[0]] if len(c) else np.nan for c in U.tcols.values], dtype=np.float64)
    return {"f_ann": np.log1p(np.asarray(pri["n_go"], dtype=np.float64)),
            "f_deg": np.log1p(pri["S"].sum(1)), "f_base_raw": base}


def assemble(q, nb, E_nb, pri, stat, size_pred, logcells_nb, lognoise_nb):
    """Feature dict for query units q from neighbour units nb (index arrays into U). A unit is never its
    own neighbour. Missing values are imputed with neighbour-set medians."""
    q, nb = np.asarray(q), np.asarray(nb)
    J = pri["J"][np.ix_(q, nb)].copy()
    Jf = pri["Jf"][np.ix_(q, nb)].copy()
    S = pri["S"][np.ix_(q, nb)].copy()
    same = q[:, None] == nb[None, :]
    J[same] = Jf[same] = S[same] = 0.0
    top = np.argsort(-J, axis=1)[:, :TOP_GO]
    Jt = np.take_along_axis(J, top, axis=1)
    wsum = Jt.sum(1)
    pos = wsum > 0
    wdiv = np.where(pos, wsum, 1.0)
    f_go = np.where(pos, (Jt * E_nb[top]).sum(1) / wdiv, np.nan)
    nbq_c = np.where(pos, (Jt * logcells_nb[top]).sum(1) / wdiv, np.nan)
    nbq_n = np.where(pos, (Jt * lognoise_nb[top]).sum(1) / wdiv, np.nan)
    ssum = S.sum(1)
    f_str = np.where(ssum > 0, (S * E_nb[None, :]).sum(1) / np.where(ssum > 0, ssum, 1), np.nan)
    f_go, f_str = np.log(f_go), np.log(f_str)
    med_logE = float(np.median(np.log(E_nb)))
    base_q = stat["f_base_raw"][q]
    base_med = float(np.nanmedian(stat["f_base_raw"][nb]))
    F = {"D_GO": d_go(J), "D_GO_filt": d_go(Jf), "D_STRING": d_string(S),
         "sim_maxJ": J.max(1), "sim_top20J": Jt.mean(1), "str_n": np.log1p((S > 0).sum(1)), "str_max": S.max(1),
         "f_go_miss": np.isnan(f_go).astype(float), "f_go": np.where(np.isnan(f_go), med_logE, f_go),
         "f_str_miss": np.isnan(f_str).astype(float), "f_str": np.where(np.isnan(f_str), med_logE, f_str),
         "f_ann": stat["f_ann"][q], "f_deg": stat["f_deg"][q],
         "f_base_miss": np.isnan(base_q).astype(float), "f_base": np.where(np.isnan(base_q), base_med, base_q),
         "nbq_cells": np.where(np.isnan(nbq_c), float(np.median(logcells_nb)), nbq_c),
         "nbq_noise": np.where(np.isnan(nbq_n), float(np.median(lognoise_nb)), nbq_n)}
    if size_pred is not None:
        F["f_model"] = np.log(size_pred + 1e-6)
    return F


def inner_folds(group_id, seed, n=N_INNER):
    """Assign TSS groups of the pool to n folds (balanced by unit count)."""
    rng = np.random.default_rng(seed)
    gid = np.asarray(group_id)
    size = np.zeros(n, int)
    fold = np.empty(len(gid), int)
    for g in rng.permutation(np.unique(gid)):
        m = gid == g
        f = int(np.argmin(size))
        fold[m] = f
        size[f] += int(m.sum())
    return fold


def feature_tables(lrn, fold, at, query, E_lrn, E_at, pri, stat, size_oof, size_full, nbq):
    """v0.6.2. OOF feature table of the learner units (E training units) and feature table of the query units.
    Neighbours of learner fold f: learner units outside f plus all always-train (AT) units. Neighbours of
    query units: all learner units plus all AT units. nbq: log cell counts / log split-half noise of
    learner (lc_lrn, ln_lrn) and AT units (lc_at, ln_at)."""
    lrn, at = np.asarray(lrn), np.asarray(at)
    F_oof = None
    for f in range(N_INNER):
        te, tr = fold == f, fold != f
        Ff = assemble(lrn[te], np.concatenate([lrn[tr], at]), np.concatenate([E_lrn[tr], E_at]), pri, stat,
                      None if size_oof is None else size_oof[te],
                      np.concatenate([nbq["lc_lrn"][tr], nbq["lc_at"]]), np.concatenate([nbq["ln_lrn"][tr], nbq["ln_at"]]))
        if F_oof is None:
            F_oof = {k: np.empty(len(lrn)) for k in Ff}
        for k, v in Ff.items():
            F_oof[k][te] = v
    F_q = assemble(query, np.concatenate([lrn, at]), np.concatenate([E_lrn, E_at]), pri, stat, size_full,
                   np.concatenate([nbq["lc_lrn"], nbq["lc_at"]]), np.concatenate([nbq["ln_lrn"], nbq["ln_at"]]))
    return F_oof, F_q


# ------------------------------------------------------------------ learners

class Learner:
    """ridge: standardised features, lambda by grouped inner CV on the fitting units (ties -> stronger
    penalty); gbm: HistGradientBoostingRegressor with fixed parameters, single thread."""

    def __init__(self, kind, seed):
        self.kind, self.seed = kind, seed

    def fit(self, X, y, groups):
        X = np.asarray(X, dtype=np.float64)
        if self.kind == "ridge":
            self.m_, self.s_ = X.mean(0), X.std(0) + 1e-12
            Z = (X - self.m_) / self.s_
            fold = inner_folds(groups, self.seed + 7)
            cv = {}
            for lam in RIDGE_LAMBDAS:
                err = 0.0
                for f in range(N_INNER):
                    tr, te = fold != f, fold == f
                    err += ((Ridge(alpha=lam).fit(Z[tr], y[tr]).predict(Z[te]) - y[te]) ** 2).sum()
                cv[lam] = err
            best = min(cv.values())
            self.lam_ = max(l for l in RIDGE_LAMBDAS if cv[l] <= best * (1 + 1e-9))
            self.model_ = Ridge(alpha=self.lam_).fit(Z, y)
        else:
            with threadpool_limits(limits=1):
                self.model_ = HistGradientBoostingRegressor(random_state=self.seed, **GBM_PARAMS).fit(X, y)
        return self

    def predict(self, X):
        X = np.asarray(X, dtype=np.float64)
        if self.kind == "ridge":
            return self.model_.predict((X - self.m_) / self.s_)
        with threadpool_limits(limits=1):
            return self.model_.predict(X)


def fit_proxy(F_oof, F_q, cols, target_log, fold, groups, kind, seed):
    """exp(learner) for query units and OOF values for the pool."""
    Xo = np.column_stack([F_oof[c] for c in cols])
    Xq = np.column_stack([F_q[c] for c in cols])
    full = Learner(kind, seed).fit(Xo, target_log, groups)
    q = np.exp(full.predict(Xq))
    oof = np.empty(len(target_log))
    for f in range(N_INNER):
        m = fold == f
        oof[m] = np.exp(Learner(kind, seed + 100 + f).fit(Xo[~m], target_log[~m], groups[~m]).predict(Xo[m]))
    return q, oof, getattr(full, "lam_", None)
