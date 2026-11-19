# lowmem_v07.py -- low-memory dataset loader for the v0.7 confirmation runs (engineering only, 2026-10-12, user
# request for GWPS under a 4 GB cap). Returns the same D dict as data_v07.load_dataset, but with the three largest
# kinds of objects represented without changing a single value:
#   Y, YA, YB  kept in RAM as the float32 rows read from the .npy files (the on-disk precision); every access
#              Y[idx] returns float64 = exactly what the frozen loader's np.load(...).astype(np.float64)[rows]
#              held (float32 -> float64 conversion is exact).
#   J, Jf      not materialised (n x n float64 each); a block J[np.ix_(q, nb)] is computed from the sparse GO
#              incidence matrix with the same arithmetic as models_v04.jaccard_all (integer intersection counts,
#              same union / division / where / zero diagonal), in row chunks.
#   S          stored sparse (same max-over-duplicates construction, zero diagonal); a block is densified on
#              access; S.sum(1) is computed on dense row chunks.
# All code that uses these objects does so only via Y[idx] / Y.shape / X[np.ix_(q, nb)] / S.sum(1)
# (checked by grep on 2026-10-12). Equality with the frozen loader is tested in test_lowmem_v07.py.
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from data_v07 import GO_TERM_MAX, gears_universe
from models_v04 import go_incidence
from paths import PROCESSED

CHUNK = 1024


class RowsF64:
    def __init__(self, a32):
        self.a = a32
        self.shape = a32.shape
        self.ndim = 2
        self.dtype = np.dtype(np.float64)

    def __len__(self):
        return self.shape[0]

    def __getitem__(self, k):
        return np.asarray(self.a[k], dtype=np.float32).astype(np.float64)


def _ix(k):
    assert isinstance(k, tuple) and len(k) == 2, "only X[np.ix_(rows, cols)] access is supported"
    return np.asarray(k[0]).ravel(), np.asarray(k[1]).ravel()


class LazyJaccard:
    def __init__(self, M):
        self.M = sp.csr_matrix(M)
        self.sz = np.asarray(self.M.sum(1)).ravel()
        self.shape = (self.M.shape[0], self.M.shape[0])

    def __getitem__(self, k):
        q, nb = _ix(k)
        out = np.empty((len(q), len(nb)))
        Mn = self.M[nb].T.tocsc()
        for s in range(0, len(q), CHUNK):
            qq = q[s:s + CHUNK]
            inter = (self.M[qq] @ Mn).toarray()
            union = self.sz[qq][:, None] + self.sz[nb][None, :] - inter
            with np.errstate(invalid="ignore", divide="ignore"):
                J = np.where(union > 0, inter / union, 0.0)
            J[qq[:, None] == nb[None, :]] = 0.0
            out[s:s + CHUNK] = J
        return out


class SparseSym:
    def __init__(self, csr):
        self.A = sp.csr_matrix(csr)
        self.shape = self.A.shape

    def __getitem__(self, k):
        q, nb = _ix(k)
        return self.A[q][:, nb].toarray()

    def sum(self, axis):
        assert axis == 1
        n = self.shape[0]
        out = np.empty(n)
        for s in range(0, n, CHUNK):
            out[s:s + CHUNK] = self.A[s:s + CHUNK].toarray().sum(1)
        return out


def load_dataset_lowmem(U_all, pb_dir, string_parquet):
    pb_dir = Path(pb_dir)
    pb = np.load(pb_dir / "pseudobulk_logmean.npy", mmap_mode="r")
    rows = pd.read_csv(pb_dir / "pseudobulk_rows.csv")
    genes = pd.read_csv(pb_dir / "pseudobulk_genes.csv").gene_name.values
    rix = {g: i for i, g in enumerate(rows.gene_transcript)}
    U = U_all.copy()
    assert U.gene_transcript.isin(rix).all(), "unit rows missing from pseudobulk"
    ctrl_rows = np.where(rows.gene_transcript.str.contains("non-targeting", case=False))[0]
    w = rows.n_cells.values[ctrl_rows].astype(np.float64)
    assert len(ctrl_rows) > 0 and np.isfinite(w).all() and (w > 0).all(), "control cell counts invalid"
    y_ctrl = (np.asarray(pb[ctrl_rows]).astype(np.float64) * w[:, None]).sum(0) / w.sum()
    U = U[U.E | U.always_train].reset_index(drop=True)
    U["row"] = U.gene_transcript.map(rix)
    U["n_cells_sc"] = rows.n_cells.values[U.row.values]
    gcols = {}
    for j, g in enumerate(genes):
        gcols.setdefault(g, []).append(j)
    U["tcols"] = [np.array(gcols.get(t, []), dtype=int) for t in U.target]
    r = U.row.values
    Y = RowsF64(np.ascontiguousarray(pb[r]))
    del pb
    YA = RowsF64(np.ascontiguousarray(np.load(pb_dir / "pseudobulk_logmean_A.npy", mmap_mode="r")[r]))
    YB = RowsF64(np.ascontiguousarray(np.load(pb_dir / "pseudobulk_logmean_B.npy", mmap_mode="r")[r]))
    gene2go, universe = gears_universe()
    cnt = Counter(t for g in set(universe) for t in gene2go[g])
    broad = frozenset(t for t, c in cnt.items() if c > GO_TERM_MAX)
    ase = np.load(PROCESSED / "priors" / "go_ase100.npz", allow_pickle=True)
    ase_ix = {g: i for i, g in enumerate(ase["genes"])}
    st = pd.read_parquet(string_parquet)
    pos = {g: i for i, g in enumerate(U.target)}
    st = st[st.g1.isin(pos) & st.g2.isin(pos)]
    d = {}
    for a, b, s in st.itertuples(index=False):          # same max-over-duplicates rule as data_v07.load_dataset
        i, j = pos[a], pos[b]
        v = max(d.get((i, j), 0.0), s)
        d[(i, j)] = v
        d[(j, i)] = v
    keys = [k for k in d if k[0] != k[1]]                 # fill_diagonal(S, 0)
    n = len(U)
    S = sp.csr_matrix(([d[k] for k in keys], ([k[0] for k in keys], [k[1] for k in keys])), shape=(n, n))
    pri = {"J": LazyJaccard(go_incidence(U.target.values, gene2go)),
           "Jf": LazyJaccard(go_incidence(U.target.values, gene2go, broad)),
           "S": SparseSym(S), "ase": ase["ase"], "ase_ix": ase_ix,
           "n_universe": len(set(universe)), "n_broad": len(broad)}
    pri["n_go"] = np.array([len(gene2go.get(t, ())) for t in U.target.values], float)
    return {"U": U, "Y": Y, "YA": YA, "YB": YB, "y_ctrl": y_ctrl, "pri": pri, "genes": genes,
            "E_idx": np.where(U.E.values)[0], "AT_idx": np.where(U.always_train.values)[0]}


# ------------------------------------------------------------------ round 2 (2026-10-12): computation-side memory
# Installed by install(); the frozen files are not edited. Rules (user): chunk ONLY along independent directions
# (rows for per-perturbation computations; gene columns for per-gene reductions across perturbations); BLAS calls
# (matrix products, dot products) and library/thread settings unchanged.
ROW_CHUNK = 256
COL_CHUNK = 512
_FROZEN = {}


def _rows(n):
    for s in range(0, n, ROW_CHUNK):
        yield s, min(s + ROW_CHUNK, n)


def score_S1_chunked(Y_true, Y_pred, ref, target_cols, k=50):
    """models_v04.score_S1, rows processed in chunks (every operation is per row: elementwise arithmetic,
    per-row argsort, per-row mean)."""
    f = _FROZEN["score_S1"]
    tc = list(target_cols)
    n = len(tc)
    out = np.empty(n)
    for a, b in _rows(n):
        out[a:b] = f(np.asarray(Y_true[a:b]), np.asarray(Y_pred[a:b]), ref, tc[a:b], k)
    return out


def e_true_chunked(Y_rows, ref, tcols):
    """proxy_v06.e_true = score_S1(Y, tile(mu), ...); the tile is built per row chunk (same values)."""
    f = _FROZEN["score_S1"]
    tc = list(tcols)
    out = np.empty(len(tc))
    for a, b in _rows(len(tc)):
        out[a:b] = f(np.asarray(Y_rows[a:b]), np.tile(ref["mu"], (b - a, 1)), ref, tc[a:b], 50)
    return out


def predicted_size_chunked(Y_pred, ref, tcols, k=50):
    f = _FROZEN["predicted_size"]
    tc = list(tcols)
    out = np.empty(len(tc))
    for a, b in _rows(len(tc)):
        out[a:b] = f(np.asarray(Y_pred[a:b]), ref, tc[a:b], k)
    return out


def split_half_stats_chunked(YA, YB, ref, tcols):
    f = _FROZEN["split_half_stats"]
    tc = list(tcols)
    rel, noise = np.empty(len(tc)), np.empty(len(tc))
    for a, b in _rows(len(tc)):
        rel[a:b], noise[a:b] = f(np.asarray(YA[a:b]), np.asarray(YB[a:b]), ref, tc[a:b])
    return rel, noise


def train_reference_chunked(Y_ref, panel_size=500):
    """models_v04.train_reference with the per-gene statistics (mean, median, MAD, var over perturbations)
    computed on contiguous gene-column blocks; the across-gene steps (5% quantile floor, panel ranking) are
    then applied to the full vectors exactly as in the frozen function."""
    n, G = Y_ref.shape
    mu, sigma, var = np.empty(G), np.empty(G), np.empty(G)
    for a in range(0, G, COL_CHUNK):
        b = min(a + COL_CHUNK, G)
        X = np.ascontiguousarray(Y_ref[:, a:b])
        m = X.mean(0)
        s = X - m
        mu[a:b] = m
        sigma[a:b] = 1.4826 * np.median(np.abs(s - np.median(s, 0)), 0)
        var[a:b] = s.var(0)
    sigma = np.maximum(sigma, np.quantile(sigma, 0.05))
    panel = np.argsort(-var)[:panel_size]
    return {"mu": mu, "sigma": sigma, "panel": panel}


class _PredStore:
    """Holds pipeline_v05.predict outputs in float64 memmaps on F: (two alternating file sets, overwritten)."""
    def __init__(self, d):
        self.d = Path(d)
        self.d.mkdir(parents=True, exist_ok=True)
        self.call = 0


_STORE = None


def predict_memmap(U, models, pri, G, idx):
    """pipeline_v05.predict, one model at a time (identical calls); each result is copied into a float64
    memmap on disk and the in-RAM array released."""
    import gc
    f = _FROZEN["predict"]
    slot = _STORE.call % 2
    _STORE.call += 1
    gc.collect()
    out = {}
    for m in models:
        y = f(U, {m: models[m]}, pri, G, idx)[m]
        mm = np.lib.format.open_memmap(_STORE.d / f"pred_{slot}_{m}.npy", mode="w+", dtype=y.dtype, shape=y.shape)
        mm[:] = y
        mm.flush()
        del y
        out[m] = np.lib.format.open_memmap(_STORE.d / f"pred_{slot}_{m}.npy", mode="r")
        del mm
    return out


def install(tmp_dir, chunk_train_reference=True, memmap_predictions=True):
    """Route the frozen modules to the chunked implementations (same values). Call once per process."""
    global _STORE
    import eval_v06
    import models_v04
    import pipeline_v05
    import proxy_v06
    if _FROZEN:
        return
    _FROZEN.update({"score_S1": models_v04.score_S1, "predicted_size": proxy_v06.predicted_size,
                    "split_half_stats": proxy_v06.split_half_stats, "train_reference": models_v04.train_reference,
                    "predict": pipeline_v05.predict})
    for mod in (eval_v06, pipeline_v05, proxy_v06):
        if hasattr(mod, "score_S1"):
            mod.score_S1 = score_S1_chunked
    proxy_v06.e_true = e_true_chunked
    proxy_v06.predicted_size = predicted_size_chunked
    proxy_v06.split_half_stats = split_half_stats_chunked
    if chunk_train_reference:
        pipeline_v05.train_reference = train_reference_chunked
    if memmap_predictions:
        _STORE = _PredStore(tmp_dir)
        pipeline_v05.predict = predict_memmap
