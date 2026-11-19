# audit_lib.py -- Perturbation-level split-conformal audit (protocol v0.2, phases P0/P1).
# Models: training-mean baseline and the Ahlmann-Eltze et al. (2025) linear baseline.
# All training-derived quantities are computed from the training fold only; leakage
# checks are implemented as assertions and as a "canary" invariance test.
import numpy as np
import pandas as pd
from scipy.stats import betabinom

MASTER_SEED = 20261001
FRACTIONS = {"train": 0.60, "val": 0.10, "cal": 0.15, "test": 0.15}


# ---------------------------------------------------------------- grouping / splits
def tss_groups(genes, chrom, strand, start, end, max_bp=1000):
    """Union-find grouping of genes whose TSSs lie within max_bp on the same chromosome."""
    tss = np.where(np.asarray(strand) == "+", np.asarray(start), np.asarray(end)).astype(np.int64)
    parent = list(range(len(genes)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    df = pd.DataFrame({"i": np.arange(len(genes)), "chr": chrom, "tss": tss}).sort_values(["chr", "tss"])
    for _, sub in df.groupby("chr"):
        idx, pos = sub.i.values, sub.tss.values
        for a in range(len(idx)):
            b = a + 1
            while b < len(idx) and pos[b] - pos[a] <= max_bp:
                ra, rb = find(idx[a]), find(idx[b])
                if ra != rb:
                    parent[rb] = ra
                b += 1
    roots = np.array([find(i) for i in range(len(genes))])
    _, gid = np.unique(roots, return_inverse=True)
    return gid


def random_group_split(group_id, seed, fractions=FRACTIONS):
    """Assign whole groups to folds in random order until each fold reaches its unit quota."""
    rng = np.random.default_rng(seed)
    groups = np.unique(group_id)
    order = rng.permutation(groups)
    n = len(group_id)
    names = list(fractions)
    quota = np.floor(np.array([fractions[k] for k in names]) * n).astype(int)
    quota[0] = n - quota[1:].sum()  # remainder goes to train
    fold = np.empty(n, dtype=object)
    counts = np.zeros(len(names), int)
    # fill the small folds first (cal, test, val) so that their sizes are exact-ish
    fill_order = [names.index(k) for k in ("test", "cal", "val", "train")]
    k = 0
    for g in order:
        members = np.where(group_id == g)[0]
        while k < len(fill_order) - 1 and counts[fill_order[k]] + len(members) > quota[fill_order[k]]:
            k += 1
        fold[members] = names[fill_order[k]]
        counts[fill_order[k]] += len(members)
    return {nm: np.where(fold == nm)[0] for nm in names}


def cluster_holdout_split(group_id, cluster_id, seed, fractions=FRACTIONS):
    """Audit B: whole prior-defined clusters go to test; cal/val drawn at random from the rest."""
    rng = np.random.default_rng(seed)
    n = len(group_id)
    n_test = int(np.floor(fractions["test"] * n))
    test = np.zeros(n, bool)
    for c in rng.permutation(np.unique(cluster_id)):
        if test.sum() >= n_test:
            break
        test[cluster_id == c] = True
    # keep TSS groups intact: a group touching test goes entirely to test
    for g in np.unique(group_id[test]):
        test[group_id == g] = True
    rest = np.where(~test)[0]
    sub = random_group_split(group_id[rest], seed + 7,
                             {"train": fractions["train"], "val": fractions["val"], "cal": fractions["cal"]})
    out = {k: rest[v] for k, v in sub.items()}
    out["test"] = np.where(test)[0]
    return out


# ---------------------------------------------------------------- models
def fit_mean(Y_train):
    """Training-mean baseline. Y_train: (n_train_perts, n_genes)."""
    b = Y_train.mean(0)
    return lambda target_gene_idx: np.tile(b, (len(target_gene_idx), 1))


def fit_linear(Y_train, train_target_gene_idx, K=10, lam=0.1):
    """Linear baseline of Ahlmann-Eltze et al. 2025 (Nat Methods, Eq. 1-3).
    Y (genes x perts) ~ G W P^T + b; G = top-K PC scores of genes from Y_train,
    P = rows of G for the perturbed genes; ridge lambda for numerical stability."""
    Y = Y_train.T.astype(np.float64)                     # genes x n_train
    b = Y.mean(1, keepdims=True)                          # row means
    Xc = Y - Y.mean(0, keepdims=True)                     # prcomp-style column centring
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    G = U[:, :K] * S[:K]                                  # genes x K
    P = G[train_target_gene_idx]                          # n_train x K
    I = np.eye(K)
    W = np.linalg.solve(G.T @ G + lam * I, G.T @ (Y - b) @ P) @ np.linalg.inv(P.T @ P + lam * I)

    def predict(target_gene_idx):
        Pt = G[np.asarray(target_gene_idx)]
        return (G @ W @ Pt.T + b).T                       # n_test x genes
    return predict


# ---------------------------------------------------------------- scores
def train_reference(Y_train, panel_size=500):
    """Train-only Systema-style reference (perturbed centroid), robust gene scale, fixed panel."""
    mu = Y_train.mean(0)
    s = Y_train - mu
    sigma = 1.4826 * np.median(np.abs(s - np.median(s, 0)), 0)
    sigma = np.maximum(sigma, np.quantile(sigma, 0.05))
    panel = np.argsort(-s.var(0))[:panel_size]
    return {"mu": mu, "sigma": sigma, "panel": panel}


def score_S1(Y_true, Y_pred, ref, target_gene_idx, k=50):
    """RMSE (sigma-scaled) on the k genes with the largest |true specific effect|, target gene excluded."""
    s = (Y_true - ref["mu"]) / ref["sigma"]
    sh = (Y_pred - ref["mu"]) / ref["sigma"]
    a = np.abs(s).copy()
    a[np.arange(len(a)), target_gene_idx] = -np.inf
    top = np.argsort(-a, axis=1)[:, :k]
    r = np.take_along_axis(s - sh, top, axis=1)
    return np.sqrt((r ** 2).mean(1))


def score_S2(Y_true, Y_pred, ref, target_gene_idx):
    """Max sigma-scaled absolute error over the fixed train panel, own target gene excluded."""
    H = ref["panel"]
    e = np.abs((Y_true[:, H] - Y_pred[:, H]) / ref["sigma"][H])
    own = (H[None, :] == np.asarray(target_gene_idx)[:, None])
    e[own] = -np.inf
    return e.max(1)


def conformal_qhat(cal_scores, alpha):
    n = len(cal_scores)
    r = int(np.ceil((n + 1) * (1 - alpha)))
    return np.inf if r > n else np.sort(cal_scores)[r - 1]


def betabinom_band(n_cal, n_test, alpha, level=0.95):
    """Predictive interval for the number of covered test points under exchangeability."""
    l = int(np.floor((n_cal + 1) * alpha))
    a, b = n_cal + 1 - l, l
    lo, hi = betabinom.ppf([(1 - level) / 2, 1 - (1 - level) / 2], n_test, a, b)
    return int(lo), int(hi)


# ---------------------------------------------------------------- prior distances
def jaccard_matrix(genes_a, genes_b, gene2go):
    A = [gene2go.get(g, set()) for g in genes_a]
    B = [gene2go.get(g, set()) for g in genes_b]
    J = np.zeros((len(A), len(B)))
    for i, x in enumerate(A):
        for j, y in enumerate(B):
            u = len(x | y)
            J[i, j] = len(x & y) / u if u else 0.0
    return J


def d_go(J_rows_to_train, top=20):
    t = min(top, J_rows_to_train.shape[1])
    return 1.0 - np.sort(J_rows_to_train, axis=1)[:, -t:].mean(1)


def d_string(S_rows_to_train):
    return 1.0 - (S_rows_to_train.max(1) if S_rows_to_train.shape[1] else 0.0)


def make_bins(d_eval, iso_frac=0.10):
    """Quintiles over cal+test distances; separate 'Qiso' bin if >=10% have D==1."""
    d_eval = np.asarray(d_eval)
    iso = np.isclose(d_eval, 1.0)
    labels = np.empty(len(d_eval), dtype=object)
    if iso.mean() >= iso_frac:
        labels[iso] = "Qiso"
        rest = ~iso
        q = np.quantile(d_eval[rest], [0.25, 0.5, 0.75])
        labels[rest] = np.array(["Q1", "Q2", "Q3", "Q4"])[np.searchsorted(q, d_eval[rest], side="right")]
    else:
        q = np.quantile(d_eval, [0.2, 0.4, 0.6, 0.8])
        labels[:] = np.array(["Q1", "Q2", "Q3", "Q4", "Q5"])[np.searchsorted(q, d_eval, side="right")]
    return labels
