# methods_v06.py -- conformal calibration methods of protocol v0.6.1 (sections 4, 7, 8).
# All functions are pure: they receive calibration scores, test scores and per-unit proxies and return,
# for every test unit, the score threshold t_i (interval = {score <= t_i}). Nothing here sees expression.
#
#   M0        global split conformal                          t = qhat
#   M1        S / (h + c), h = predicted effect strength      t = qhat_norm * (h + c), c = gamma * median(h_oof)
#   M2        Mondrian on h (edges = quantiles of h_oof)
#   M3        Gibbs, Cherian & Candes (2025) conditional conformal, Phi(h); exact one-LP implementation
#   Mdist-N   as M1 with d = similarity / training-quality difficulty proxy (PRESCRIBE-like baseline)
#   Mdist-M   as M2 with d
#   Mscore-N  as M1 with s_hat = all h-features fitted to the out-of-fold score
#   Mcomb-N   as M1 with sigma_comb = exp(OLS of log OOF score on [log h, log d])
#   Mcomb-3   M3 with Phi = span{1, z(log h), z(log d)}
#   oracle    as M1 with the true effect strength E (mechanistic upper bound, not deployable)
# M3-type randomisation: each test unit has a fixed u ~ Unif(0, 1); U = u - alpha, so thresholds are
# monotone in alpha (needed for width matching).
import numpy as np
from scipy.optimize import linprog

ALPHA = 0.10


def qhat(scores, alpha=ALPHA):
    """Split-conformal quantile: the ceil((n+1)(1-alpha))-th smallest score; +inf if that exceeds n."""
    s = np.sort(np.asarray(scores, dtype=np.float64))
    n = len(s)
    k = int(np.ceil((n + 1) * (1 - alpha)))
    return np.inf if k > n else float(s[k - 1])


def m0(s_cal, s_test, alpha=ALPHA):
    t = np.full(len(s_test), qhat(s_cal, alpha))
    return t


def m1(s_cal, h_cal, s_test, h_test, c, alpha=ALPHA):
    q = qhat(np.asarray(s_cal) / (np.asarray(h_cal) + c), alpha)
    return q * (np.asarray(h_test) + c)


def mondrian_edges(h_oof, n_bins):
    """Interior bin edges from quantiles of the training out-of-fold proxy (never cal/test)."""
    return np.quantile(np.asarray(h_oof, dtype=np.float64), np.linspace(0, 1, n_bins + 1)[1:-1])


def m2(s_cal, h_cal, s_test, h_test, edges, alpha=ALPHA):
    b_cal = np.searchsorted(edges, h_cal, side="right")
    b_test = np.searchsorted(edges, h_test, side="right")
    q = {b: qhat(np.asarray(s_cal)[b_cal == b], alpha) for b in range(len(edges) + 1)}
    return np.array([q[b] for b in b_test]), b_test


# ------------------------------------------------------------------ M3: Gibbs-Cherian-Candes

def phi_basis(kind, h_ref):
    """Returns a function h -> design matrix. kind: "quad" = span{1, h, h^2} (h standardised by the
    training out-of-fold proxy); "ns3"/"ns5" = natural cubic spline of h with df 3 / 5 (+ intercept),
    knots at quantiles of h_oof."""
    h_ref = np.asarray(h_ref, dtype=np.float64)
    mu, sd = h_ref.mean(), h_ref.std() + 1e-12
    if kind == "quad":
        def f(h):
            z = (np.asarray(h, dtype=np.float64) - mu) / sd
            return np.column_stack([np.ones_like(z), z, z ** 2])
        return f
    df = {"ns3": 3, "ns5": 5}[kind]
    # natural cubic spline, truncated-power basis (Hastie et al., ESL eq. 5.4-5.5); df basis functions
    knots = np.quantile(h_ref, np.linspace(0, 1, df + 1))     # df+1 knots incl. boundary -> df columns
    K = len(knots)

    def d(z, k):
        return (np.maximum(z - knots[k], 0) ** 3 - np.maximum(z - knots[K - 1], 0) ** 3) / (knots[K - 1] - knots[k])

    def f(h):
        z = np.asarray(h, dtype=np.float64)
        cols = [np.ones_like(z), (z - mu) / sd]
        for k in range(K - 2):
            cols.append((d(z, k) - d(z, K - 2)) / sd ** 3)
        return np.column_stack(cols)
    return f


def _eta_test(Phi_cal, s_cal, phi_t, s_t, alpha):
    """Dual of the augmented quantile regression (level tau = 1 - alpha):
         max  sum_i eta_i S_i   s.t.  eta_i in [-alpha, 1 - alpha],  sum_i eta_i phi_i = 0
       over calibration points plus the test point with imputed score s_t. Returns eta_{n+1}."""
    A = np.vstack([Phi_cal, phi_t[None, :]]).T
    c = -np.r_[s_cal, s_t]
    r = linprog(c, A_eq=A, b_eq=np.zeros(A.shape[0]), bounds=[(-alpha, 1 - alpha)] * A.shape[1],
                method="highs")
    if r.status != 0:
        raise RuntimeError(f"M3 LP failed: {r.message}")
    return float(r.x[-1])


M3_FALLBACKS = {"ipm": 0, "bisect": 0}


def _eta_range_max(Phi_cal, phi_t, alpha):
    A = np.vstack([Phi_cal, phi_t[None, :]]).T
    c = np.zeros(A.shape[1])
    c[-1] = -1.0
    for meth in ("highs", "highs-ipm"):
        r = linprog(c, A_eq=A, b_eq=np.zeros(A.shape[0]), bounds=[(-alpha, 1 - alpha)] * A.shape[1], method=meth)
        if r.status == 0:
            return -r.fun
    raise RuntimeError(f"M3 range LP failed: {r.message}")


def m3_threshold(Phi_cal, s_cal, phi_t, U, alpha=ALPHA):
    """Exact randomised threshold of the Gibbs-Cherian-Candes set for one test point.
    With g(e) = max{sum_i eta_i S_i : eta_i in [-alpha, 1-alpha], sum_i eta_i phi_i = -e phi_t}, the
    test dual eta_{n+1}(S) = argmax_e g(e) + e S, so {S : eta_{n+1}(S) < U} = {S <= -g'(U)}.
    -g'(U) is read off the equality-constraint marginals y of one LP: threshold = -(y . phi_t) with
    scipy's sign convention (cross-checked against conditionalconformal 0.0.5, protocol v0.6.1 11c).
    Numerical fallback (D0', 2026-10-06; same mathematical quantity): if the dual simplex returns a
    solver failure (HiGHS 'Not Set', seen near the feasibility boundary at small alpha with ns5), retry
    with the interior-point solver; if that also fails, bisect on eta_{n+1}(S) (_eta_test)."""
    s_cal = np.asarray(s_cal, dtype=np.float64)
    kw = dict(A_eq=Phi_cal.T, b_eq=-U * phi_t, bounds=[(-alpha, 1 - alpha)] * len(s_cal))
    r = linprog(-s_cal, method="highs", **kw)
    if r.status not in (0, 2):
        M3_FALLBACKS["ipm"] += 1
        r = linprog(-s_cal, method="highs-ipm", **kw)
    if r.status == 2:
        # eta_{n+1} = U infeasible: U outside the attainable range of eta_{n+1}.
        # U > e_max -> eta_{n+1}(S) < U for every S -> +inf (everything); U < e_min -> -inf (empty set).
        return np.inf if U > _eta_range_max(Phi_cal, phi_t, alpha) else -np.inf
    if r.status == 0:
        return -float(r.eqlin.marginals @ phi_t)
    M3_FALLBACKS["bisect"] += 1
    if U > _eta_range_max(Phi_cal, phi_t, alpha):
        return np.inf
    hi = 4 * float(np.abs(s_cal).max()) + 1
    for _ in range(60):                      # eta_{n+1}(S) -> e_max as S -> inf, so hi exists
        if _eta_test(Phi_cal, s_cal, phi_t, hi, alpha) >= U:
            break
        hi *= 2
    lo = -4 * float(np.abs(s_cal).max()) - 1
    for _ in range(60):
        if _eta_test(Phi_cal, s_cal, phi_t, lo, alpha) < U:
            break
        lo *= 2
    else:
        return -np.inf
    for _ in range(60):
        mid = (lo + hi) / 2
        if _eta_test(Phi_cal, s_cal, phi_t, mid, alpha) < U:
            lo = mid
        else:
            hi = mid
    return lo


def m3(s_cal, h_cal, s_test, h_test, phi_fn, alpha=ALPHA, seed=0):
    """Randomised Gibbs-Cherian-Candes conditional conformal; U_i = u_i - alpha, u_i ~ Unif(0, 1)."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(size=len(s_test))
    thr = _gcc_thresholds(np.asarray(s_cal, float), phi_fn(h_cal), phi_fn(h_test), u, alpha)
    return thr, np.asarray(s_test, dtype=np.float64) <= thr


# ------------------------------------------------------------------ dispatcher

M1_GAMMAS = (0.0, 0.25, 0.5, 1.0)
M2_BINS = (3, 5)
M3_PHI = ("quad", "ns3", "ns5")
ALPHA_GRID = tuple(np.round(np.arange(0.010, 0.4001, 0.005), 3))
# method -> (family, proxy variable)
METHODS = {"M0": ("global", None), "M1": ("norm", "h"), "M2": ("mondrian", "h"), "M3": ("gcc", "h"),
           "Mdist-N": ("norm", "d"), "Mdist-M": ("mondrian", "d"), "Mscore-N": ("norm", "s"),
           "Mcomb-N": ("norm", "comb"), "Mcomb-3": ("gcc2", ("h", "d")), "oracle": ("norm", "E"),
           # v0.7 development additions (2026-10-10, after D1 of v0.6.2; see ext_v07dev.py)
           "MD-GO": ("mondrian", "dgo"), "MD-STR": ("mondrian", "dstr"),
           "Munc-N": ("norm", "unc"), "Munc-M": ("mondrian", "unc")}


def phi_lin2(h_ref, d_ref):
    """Phi = span{1, z(log h), z(log d)}, standardised by the training out-of-fold proxies."""
    lh, ld = np.log(np.asarray(h_ref, float)), np.log(np.asarray(d_ref, float))
    mh, sh, md, sd = lh.mean(), lh.std() + 1e-12, ld.mean(), ld.std() + 1e-12

    def f(hd):
        hd = np.asarray(hd, float).reshape(-1, 2)
        return np.column_stack([np.ones(len(hd)), (np.log(hd[:, 0]) - mh) / sh, (np.log(hd[:, 1]) - md) / sd])
    return f


def _gcc_thresholds(s_cal, Phi_cal, Phi_q, u, alpha):
    return np.array([m3_threshold(Phi_cal, s_cal, Phi_q[i], u[i] - alpha, alpha) for i in range(len(Phi_q))])


def thresholds(name, hyper, s_cal, P_cal, P_q, P_oof, u_q, alpha=ALPHA):
    """Thresholds for query units q (test units, or calibration units for the matching sensitivity).
    P_cal / P_q / P_oof: dicts of proxy arrays ("h", "d", "s", "comb", "E") for calibration units, query
    units and training-pool out-of-fold values. u_q: fixed uniforms for M3-type methods."""
    fam, var = METHODS[name]
    s_cal = np.asarray(s_cal, float)
    nq = len(u_q)
    if fam == "global":
        return np.full(nq, qhat(s_cal, alpha))
    if fam == "norm":
        c = hyper * float(np.median(P_oof[var]))
        return m1(s_cal, P_cal[var], np.zeros(nq), P_q[var], c, alpha)
    if fam == "mondrian":
        t, _ = m2(s_cal, P_cal[var], np.zeros(nq), P_q[var], mondrian_edges(P_oof[var], hyper), alpha)
        return t
    if fam == "gcc":
        f = phi_basis(hyper, P_oof[var])
        return _gcc_thresholds(s_cal, f(P_cal[var]), f(P_q[var]), u_q, alpha)
    if fam == "gcc2":
        f = phi_lin2(P_oof["h"], P_oof["d"])
        stack = lambda P: np.column_stack([P["h"], P["d"]])           # noqa: E731
        return _gcc_thresholds(s_cal, f(stack(P_cal)), f(stack(P_q)), u_q, alpha)
    raise ValueError(name)


def capped_width(t, s_cal):
    """Mean interval radius with +inf -> 2 x max calibration score and empty sets (t < 0) -> 0."""
    cap = 2.0 * float(np.max(s_cal))
    return float(np.clip(np.where(np.isposinf(t), cap, t), 0.0, cap).mean())


def match_width(name, hyper, s_cal, P_cal, P_q, P_oof, u_q, target, grid=ALPHA_GRID):
    """alpha' on the grid whose capped mean width over the query units is closest to `target`.
    Never uses query-unit scores. Cheap families: full grid. M3-type: bisection on the grid index
    (width is non-increasing in alpha), then the closer of the two bracketing points.
    Returns (alpha', thresholds at alpha', relative width residual, monotone_ok)."""
    fam = METHODS[name][0]
    cache = {}

    def W(i):
        if i not in cache:
            t = thresholds(name, hyper, s_cal, P_cal, P_q, P_oof, u_q, grid[i])
            cache[i] = (capped_width(t, s_cal), t)
        return cache[i][0]
    if fam in ("gcc", "gcc2"):
        lo, hi = 0, len(grid) - 1                    # W(lo) >= W(hi)
        if W(lo) <= target:
            best = lo
        elif W(hi) >= target:
            best = hi
        else:
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if W(mid) >= target:
                    lo = mid
                else:
                    hi = mid
            best = lo if abs(W(lo) - target) <= abs(W(hi) - target) else hi
    else:
        for i in range(len(grid)):
            W(i)
        best = min(cache, key=lambda i: (abs(cache[i][0] - target), abs(grid[i] - ALPHA)))   # ties -> nearest 0.1
    ks = sorted(cache)
    ws = [cache[i][0] for i in ks]
    mono = bool(all(ws[j] >= ws[j + 1] - 1e-9 for j in range(len(ws) - 1)))
    return float(grid[best]), cache[best][1], abs(cache[best][0] - target) / target, mono


def permute_within(h_cal, h_test, seed):
    rng = np.random.default_rng(seed)
    return rng.permutation(np.asarray(h_cal)), rng.permutation(np.asarray(h_test))


def run_method(name, hyper, s_cal, s_test, h_cal=None, h_test=None, h_oof=None, alpha=ALPHA, seed=0):
    """v0.6 interface (kept for the unit tests): proxy h only; returns (thresholds, covered)."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(size=len(s_test))
    nm = "oracle" if name == "oracle" else name
    var = METHODS[nm][1]
    P_cal, P_q, P_oof = {}, {}, {}
    if var is not None:
        key = var if isinstance(var, str) else "h"
        P_cal[key], P_q[key], P_oof[key] = h_cal, h_test, h_oof
    t = thresholds(nm, hyper, s_cal, P_cal, P_q, P_oof, u, alpha)
    return t, np.asarray(s_test, float) <= t
