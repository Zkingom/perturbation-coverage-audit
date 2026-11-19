# models_v04.py -- models, reference and prior distances for protocol v0.4.
# Linear models follow the released code of Ahlmann-Eltze et al. 2025 (Zenodo record 16092690,
# benchmark/src/run_linear_pretrained_model.R, solve_y_axb): ctrl is a training column with a
# zero perturbation embedding; Y is centred by its row means over the fitted columns; the gene
# embedding G is the top-K PCA scores of the training pseudobulk (genes x training columns,
# column-centred as in prcomp); ridge penalty lambda = 0.1 on both sides.
import numpy as np
import scipy.sparse as sp


# ------------------------------------------------------------------ reference (protocol (d))
def train_reference(Y_ref, panel_size=500):
    """Systema-style train-only reference; Y_ref = training perturbations with >= 30 cells."""
    mu = Y_ref.mean(0)
    s = Y_ref - mu
    sigma = 1.4826 * np.median(np.abs(s - np.median(s, 0)), 0)
    sigma = np.maximum(sigma, np.quantile(sigma, 0.05))
    panel = np.argsort(-s.var(0))[:panel_size]
    return {"mu": mu, "sigma": sigma, "panel": panel}


# ------------------------------------------------------------------ models
def fit_mean(Y_train):
    b = Y_train.mean(0)
    return lambda n: np.tile(b, (n, 1))


def gene_pca_embedding(Y_cols, K=10, method="gram"):
    """Y_cols: (n_columns, n_genes) training pseudobulk incl. ctrl. Returns G (n_genes x K).
    method: "svd"   numpy SVD (LAPACK gesdd), as in v0.2;
            "gesvd" scipy SVD with the LAPACK gesvd driver (same decomposition, other algorithm);
            "gram"  eigendecomposition of the columns x columns Gram matrix; G = Xc V equals U S
                    up to column sign. PROVISIONAL: adopted only if it matches SVD on real data
                    (protocol v0.4 section (f)); otherwise SVD is restored."""
    X = Y_cols.T.astype(np.float64)                 # genes x columns
    Xc = X - X.mean(0, keepdims=True)               # prcomp centres each column (variable)
    if method == "gram":
        w, V = np.linalg.eigh(Xc.T @ Xc)
        order = np.argsort(w)[::-1][:K]
        return Xc @ V[:, order]
    if method == "svd":
        U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    elif method == "gesvd":
        import scipy.linalg
        U, S, _ = scipy.linalg.svd(Xc, full_matrices=False, lapack_driver="gesvd")
    else:
        raise ValueError(method)
    return U[:, :K] * S[:K]


def solve_y_axb(Y, A, B, lam=0.1):
    """Y (genes x n), A (genes x K), B (d x n). Returns K_mat (K x d), centre (genes,)."""
    centre = Y.mean(1)
    Yc = Y - centre[:, None]
    K = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ Yc @ B.T) @ np.linalg.inv(B @ B.T + lam * np.eye(B.shape[0]))
    return K, centre


def fit_linear(Y_fit, P_fit, G, lam=0.1):
    """Y_fit: (n_fit, genes) incl. ctrl row; P_fit: (n_fit, d) with zeros for ctrl."""
    Km, centre = solve_y_axb(Y_fit.T.astype(np.float64), G, P_fit.T.astype(np.float64), lam)

    def predict(P_new):
        return (G @ Km @ np.asarray(P_new, dtype=np.float64).T + centre[:, None]).T
    return predict


# ------------------------------------------------------------------ scores
def score_S1(Y_true, Y_pred, ref, target_cols, k=50):
    """sigma-scaled RMSE on the k genes with largest |true specific effect|; own target gene(s) excluded.
    target_cols: list of arrays of readout-column indices to exclude for each row."""
    s = (Y_true - ref["mu"]) / ref["sigma"]
    sh = (Y_pred - ref["mu"]) / ref["sigma"]
    a = np.abs(s).copy()
    for i, cols in enumerate(target_cols):
        a[i, cols] = -np.inf
    top = np.argsort(-a, axis=1)[:, :k]
    r = np.take_along_axis(s - sh, top, axis=1)
    return np.sqrt((r ** 2).mean(1))


def score_S2(Y_true, Y_pred, ref, target_cols):
    H = ref["panel"]
    e = np.abs((Y_true[:, H] - Y_pred[:, H]) / ref["sigma"][H])
    for i, cols in enumerate(target_cols):
        e[i, np.isin(H, cols)] = -np.inf
    return e.max(1)


# ------------------------------------------------------------------ prior distances
def go_incidence(genes, gene2go, drop_terms=frozenset()):
    terms = sorted({t for g in genes for t in gene2go.get(g, ()) if t not in drop_terms})
    ti = {t: j for j, t in enumerate(terms)}
    r, c = [], []
    for i, g in enumerate(genes):
        for t in gene2go.get(g, ()):
            if t in ti:
                r.append(i)
                c.append(ti[t])
    return sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(genes), len(terms)))


def jaccard_all(M):
    inter = (M @ M.T).toarray()
    sz = np.asarray(M.sum(1)).ravel()
    union = sz[:, None] + sz[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        J = np.where(union > 0, inter / union, 0.0)
    np.fill_diagonal(J, 0.0)
    return J


def d_go(J_rows_to_pool, top=20):
    t = min(top, J_rows_to_pool.shape[1])
    return 1.0 - np.sort(J_rows_to_pool, axis=1)[:, -t:].mean(1)


def d_string(S_rows_to_pool):
    return 1.0 - S_rows_to_pool.max(1)
