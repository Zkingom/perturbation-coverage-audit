# synth_v06.py -- synthetic Perturb-seq-like datasets for D0 end-to-end tests (protocol v0.6 section 11).
# Effect strength a_t = exp(mu + b_mod * z_module(t) + b_base * zbase_t + e_t). The response direction of
# target t is a fixed linear map of a gene latent vector, so linear models can learn the direction but not
# the magnitude; model error therefore scales with a_t (the planted "score proportional to h" mechanism).
#   positive control: b_mod = 0.8, b_base = 0.3 -> a_t predictable from GO/STRING neighbours and f_base
#   negative control: b_mod = b_base = 0      -> a_t unpredictable
#   distance-driven (v0.6.1): b_mod = b_base = 0 and orphan_frac > 0 -> a_t unpredictable; "orphan" targets
#     have no module GO terms / STRING edges and a latent response direction unrelated to their expression
#     loadings, so their direction cannot be learned: difficulty is driven by distance to training units.
# v0.6.1 change: the readout loading of a target gene equals its latent vector plus noise (co-regulated
# genes respond alike), so linear-PCA can learn response directions of non-orphan targets.
# Cell counts are independent of a_t. Optional split halves YA / YB (independent noise, half the cells).
import numpy as np
import pandas as pd


def make_dataset(seed, b_mod=0.8, b_base=0.3, sd_e=0.4, n_E=700, n_AT=120, G=1000, n_mod=40, halves=True,
                 orphan_frac=0.0, unpred_at_frac=0.0):
    rng = np.random.default_rng(seed)
    n = n_E + n_AT
    mod = rng.integers(0, n_mod, n)
    zmod = rng.normal(0, 1, n_mod)
    # control baseline (log1p CP10k scale); target genes are readout columns 0..n-1
    y_ctrl = rng.normal(1.0, 0.4, G).clip(0.05)
    zbase = (y_ctrl[:n] - y_ctrl[:n].mean()) / y_ctrl[:n].std()
    loga = np.log(0.6) + b_mod * zmod[mod] + b_base * zbase + rng.normal(0, sd_e, n)
    a = np.exp(loga)
    # gene latent vectors (module-structured for target genes) and response loadings
    L = 6
    cent = rng.normal(0, 1, (n_mod, L))
    v = rng.normal(0, 1, (G, L))
    v[:n] = 0.8 * cent[mod] + 0.6 * rng.normal(0, 1, (n, L))
    orphan = rng.uniform(size=n) < orphan_frac
    Lam = rng.normal(0, 1, (G, L))
    Lam[:n] = v[:n] + 0.3 * rng.normal(0, 1, (n, L))                 # target genes: loading ~ latent
    v[:n][orphan] = rng.normal(0, 1, (int(orphan.sum()), L))          # orphans: direction not learnable
    tau = rng.uniform(0.2, 0.5, G)                                    # per-gene cell-level noise sd
    ncells = np.maximum(10, rng.lognormal(np.log(80), 0.5, n)).astype(int)
    U_dir = (Lam @ v[:n].T).T
    U_dir /= np.linalg.norm(U_dir, axis=1, keepdims=True) / np.sqrt(G)
    sig_unit = tau[None, :] / np.sqrt(ncells[:, None])
    resp = a[:, None] * 0.25 * U_dir * tau[None, :]
    resp[np.arange(n), np.arange(n)] = -1.0                            # knockdown of the target itself
    Y = y_ctrl[None, :] + resp + rng.normal(0, 1, (n, G)) * sig_unit
    YA = YB = None
    if halves:
        sh = sig_unit * np.sqrt(2)
        YA = y_ctrl[None, :] + resp + rng.normal(0, 1, (n, G)) * sh
        YB = y_ctrl[None, :] + resp + rng.normal(0, 1, (n, G)) * sh
        Y = (YA + YB) / 2                                               # full = average of the halves
    # GO annotations: 3 of the module's 5 terms + 2 random global terms
    n_terms = n_mod * 5 + 300
    go = []
    for t in range(n):
        own = rng.choice(5, 0 if orphan[t] else 3, replace=False) + 5 * mod[t]
        glob = rng.choice(300, 3 if orphan[t] else 2, replace=False) + 5 * n_mod
        go.append(set(own.tolist()) | set(glob.tolist()))
    M = np.zeros((n, n_terms))
    for t, s in enumerate(go):
        M[t, list(s)] = 1
    inter = M @ M.T
    cnt = M.sum(1)
    J = inter / (cnt[:, None] + cnt[None, :] - inter)
    np.fill_diagonal(J, 0)
    # STRING-like weights
    same = mod[:, None] == mod[None, :]
    p = np.where(same, 0.15, 0.005)
    p[orphan, :] = 0.005
    p[:, orphan] = 0.005
    E = np.triu(rng.uniform(size=(n, n)) < p, 1)
    S = np.where(E, rng.uniform(0.4, 0.99, (n, n)), 0.0)
    S = S + S.T
    targets = np.array([f"G{t:04d}" for t in range(n)])
    tss = np.arange(n).astype(object)
    pair = rng.choice(n, int(0.05 * n), replace=False)
    for i in range(0, len(pair) - 1, 2):
        tss[pair[i + 1]] = tss[pair[i]]
    U = pd.DataFrame({"target": targets, "tss_group": [f"tss{x}" for x in tss], "ncells": ncells,
                      "E": np.r_[np.ones(n_E, bool), np.zeros(n_AT, bool)]})
    U["always_train"] = ~U.E
    # a TSS group must not mix E and always-train units: move such AT units into E
    mixed = U.groupby("tss_group").E.transform(lambda x: x.any() and not x.all())
    U.loc[mixed, "E"] = True
    U["always_train"] = ~U.E
    # readout gene names: column t is the target gene of unit t. v0.6.2: a fraction of always-train units
    # get a target that is NOT a readout gene (column renamed) and that is missing from the GO embedding,
    # mimicking the real data (targets not measured / not embedded); their responses are unchanged.
    genes = np.array([f"G{t:04d}" for t in range(n)] + [f"R{g:04d}" for g in range(n, G)], dtype=object)
    at_ix = np.where(U.always_train.values)[0]
    unpred = np.zeros(n, bool)
    if unpred_at_frac > 0:
        unpred[rng.choice(at_ix, int(round(unpred_at_frac * len(at_ix))), replace=False)] = True
        genes[np.where(unpred)[0]] = [f"X{t:04d}" for t in np.where(unpred)[0]]
    U["tcols"] = [np.array([t]) if not unpred[t] else np.array([], dtype=int) for t in range(n)]
    ase = 0.7 * np.repeat(rng.normal(0, 1, (n_mod, 100)), 1, 0)[mod] + 0.7 * rng.normal(0, 1, (n, 100))
    pri = {"J": J, "Jf": J, "S": S, "ase": ase, "ase_ix": {t: i for i, t in enumerate(targets) if not unpred[i]},
           "n_go": cnt, "n_universe": n, "n_broad": 0}
    D = {"U": U, "Y": Y, "y_ctrl": y_ctrl, "pri": pri, "E_idx": np.where(U.E.values)[0],
         "AT_idx": np.where(U.always_train.values)[0], "genes": genes,
         "truth": {"a": a, "module": mod, "orphan": orphan, "unpredictable_at": unpred}}
    if halves:
        D["YA"], D["YB"] = YA, YB
    return D
