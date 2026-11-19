# p1_run.py -- run one audit split for mean + linear baselines (protocol v0.2).
# ARCHIVED: v0.2 is no longer reproduced (decision 2026-10-02); kept as the historical
# implementation. The v0.4 P1 runner will be a separate file.
# Requires in namespace: pb (conditions x genes pseudobulk), U (unit table for E),
# JE (GO Jaccard, E x E), Smat (STRING exp+db score, E x E), audit_lib functions.
import os, stat, json
import numpy as np, pandas as pd

ALPHAS = (0.05, 0.10, 0.20)
KS = (20, 50, 200)


def predictions_and_reference(PB, folds):
    """Everything that may depend on expression is computed here from the TRAIN fold only."""
    tr = folds["train"]
    Ytr = PB[U.row.values[tr]]
    ref = train_reference(Ytr)
    models = {"mean": fit_mean(Ytr),
              "linear": fit_linear(Ytr, U.gidx.values[tr], K=10, lam=0.1)}
    ev = np.concatenate([folds["cal"], folds["test"]])
    preds = {m: f(U.gidx.values[ev]) for m, f in models.items()}
    return ref, preds, ev


def prior_distances(folds, ev):
    tr = folds["train"]
    return d_go(JE[np.ix_(ev, tr)], top=20), d_string(Smat[np.ix_(ev, tr)])


def leakage_asserts(folds):
    allidx = np.concatenate(list(folds.values()))
    assert len(allidx) == len(U) and len(np.unique(allidx)) == len(U), "L1/L10: folds not a partition"
    fold_of = np.empty(len(U), object)
    for k, v in folds.items():
        fold_of[v] = k
    # L1: one row per target gene, so each target in exactly one fold
    assert U.target.is_unique
    # L2: every TSS group lies in a single fold
    g = pd.DataFrame({"grp": U.tss_group.values, "fold": fold_of}).groupby("grp").fold.nunique()
    assert (g == 1).all(), "L2: TSS group split across folds"


def canary_check(folds, ref, preds, ev, dG, dS, seed):
    """L3/L4: corrupt cal+test expression arbitrarily; train-derived objects must be unchanged."""
    rng = np.random.default_rng(seed + 999)
    PBc = pb.copy()
    rows = U.row.values[np.concatenate([folds["cal"], folds["test"], folds["val"]])]
    PBc[rows] = rng.normal(5.0, 3.0, size=PBc[rows].shape).astype(PBc.dtype)
    ref2, preds2, ev2 = predictions_and_reference(PBc, folds)
    dG2, dS2 = prior_distances(folds, ev2)
    ok = (np.array_equal(ev, ev2)
          and all(np.allclose(ref[k], ref2[k]) for k in ("mu", "sigma"))
          and np.array_equal(ref["panel"], ref2["panel"])
          and all(np.allclose(preds[m], preds2[m]) for m in preds)
          and np.allclose(dG, dG2) and np.allclose(dS, dS2))
    assert ok, "L3/L4 canary failed: a train-derived quantity depends on held-out expression"
    return ok


def write_readonly_folds(folds, split_name, outdir=None):
    # Paths are resolved from config.toml via paths.py (no hard-coded drive letters).
    if outdir is None:
        from paths import SPLITS
        outdir = str(SPLITS / "v0.2")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{split_name}.json")
    if os.path.exists(path):
        os.chmod(path, stat.S_IWRITE | stat.S_IREAD)
    with open(path, "w") as fh:
        json.dump({k: U.target.values[v].tolist() for k, v in folds.items()}, fh)
    os.chmod(path, stat.S_IREAD)  # L10: read-only once written
    return path


def run_split(i, audit="A", plant=None, permute_d=False, folds=None, check=True):
    seed = MASTER_SEED + i
    if folds is None:
        folds = random_group_split(U.tss_group.values, seed)
    if check:
        leakage_asserts(folds)
    ref, preds, ev = predictions_and_reference(pb, folds)
    dG, dS = prior_distances(folds, ev)
    if check:
        canary_check(folds, ref, preds, ev, dG, dS, seed)
        write_readonly_folds(folds, f"{audit}_{i:03d}")
    ncal = len(folds["cal"])
    is_cal = np.r_[np.ones(ncal, bool), np.zeros(len(folds["test"]), bool)]
    Ytrue = pb[U.row.values[ev]]
    tg = U.gidx.values[ev]
    binG, binS = make_bins(dG), make_bins(dS)
    rng = np.random.default_rng(seed + 12345)
    rows = []
    for m, Yp in preds.items():
        Yp = Yp.copy()
        if plant is not None and m == "linear":
            # positive control: noise growing with D_GO rank, applied to TEST units only
            u = pd.Series(dG[~is_cal]).rank(pct=True).values
            Yp[~is_cal] += plant * u[:, None] * ref["sigma"][None, :] * rng.standard_normal(Yp[~is_cal].shape)
        sc = {f"S1_k{k}": score_S1(Ytrue, Yp, ref, tg, k) for k in KS}
        sc["S2"] = score_S2(Ytrue, Yp, ref, tg)
        dG_use = dG.copy()
        if permute_d:
            t = np.where(~is_cal)[0]
            dG_use[t] = dG[rng.permutation(t)]
        bG = make_bins(dG_use) if permute_d else binG
        for sname, s in sc.items():
            for a in ALPHAS:
                q = conformal_qhat(s[is_cal], a)
                for j in np.where(~is_cal)[0]:
                    rows.append((audit, i, m, sname, a, U.target.values[ev[j]], dG_use[j], dS[j], bG[j], binS[j],
                                 s[j], q, bool(s[j] <= q), ncal))
    return pd.DataFrame(rows, columns=["audit", "split", "model", "score", "alpha", "target", "D_GO", "D_STRING",
                                       "bin_GO", "bin_STRING", "score_value", "qhat", "covered", "n_cal"])
