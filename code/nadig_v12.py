# nadig_v12.py -- protocol combo_crosscell_v1.2 section 3.4 (as amended by v1.2.1) on Nadig et al. 2025 HepG2 and
# Jurkat (CRISPRi, essential library; scPerturb v1.4 files). Prediction: linear-PCA R not elevated.
#   python nadig_v12.py c0 <hepg2|jurkat>       pseudobulk + halves (one backed pass), unit table, read-only splits
#   python nadig_v12.py run <hepg2|jurkat>      per-split mean model + linear-PCA, conformal sets, controls (resumable)
#   python nadig_v12.py summary <hepg2|jurkat>  pooled statistics, TSS-group cluster bootstrap, pre-specified judgement
# Rules reused unchanged: per-cell transform x = log1p(1e4 * counts / L) (pseudobulk_v06; L = obs ncounts per
# v1.2.1 section 5), half split pseudobulk_v06.half_labels, units / E / TSS groups as data_v07.build_units (unit
# metadata aggregated from single-cell obs, protocol 3.4), splits_v06.repeated_group_kfold, linear-PCA as
# pipeline_v05 / eval_v06.refit_predict (K on validation from {5,10,20,30,50}, lambda 0.1), statistics and
# judgement functions of norman_v12 (same constants). No algorithm or parameter is changed.
import json
import os
import pickle
import stat
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED, RESULTS, SPLITS  # noqa: E402
from models_v04 import fit_linear, fit_mean, gene_pca_embedding, score_S1, train_reference  # noqa: E402
from methods_v06 import qhat  # noqa: E402
from eval_v06 import tiers_of  # noqa: E402
from splits_v06 import repeated_group_kfold  # noqa: E402
from pipeline_v05 import _pick  # noqa: E402
from pseudobulk_v06 import half_labels  # noqa: E402
import data_v07 as dv  # noqa: E402
import norman_v12 as nv  # noqa: E402
import resume_v07 as rs  # noqa: E402
from safe_io import create_new  # noqa: E402

DL = Path(r"G:\biology\scperturb_13350497")
DS_ID = {"hepg2": 1, "jurkat": 2}
K_GRID, LAM, K_SCORE, ALPHA, N_PERM = nv.K_GRID, nv.LAM, nv.K_SCORE, nv.ALPHA, nv.N_PERM
BLOCK = 2000
MODELS = ("mean", "linear_pca")


def dirs(ds):
    return PROCESSED / "v1.2" / ds, SPLITS / "v1.2" / ds, RESULTS / "v1.2" / ds


def log(ds, msg):
    od = dirs(ds)[2]
    od.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with open(od / "run.log", "a", encoding="utf8") as fh:
        fh.write(line + "\n")


def col(g, k):
    x = g[k]
    if isinstance(x, h5py.Group):
        c = np.array([v.decode() if isinstance(v, bytes) else v for v in x["categories"][:]], dtype=object)
        return c[x["codes"][:]]
    v = x[:]
    return np.array([a.decode() if isinstance(a, bytes) else a for a in v], dtype=object) if v.dtype.kind in "SO" else v


def save_npy_new(arr, path):
    if path.exists():
        raise FileExistsError(path)
    with open(path, "xb") as fh:
        np.save(fh, arr)
    os.chmod(path, stat.S_IREAD)


# ------------------------------------------------------------------ C0
def c0(ds):
    pdir, sdir, odir = dirs(ds)
    if (pdir / "c0_build.json").exists():
        log(ds, "C0 already complete; nothing to do")
        return
    pdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with h5py.File(DL / f"NadigOConner2024_{ds}.h5ad", "r") as h:
        X = h["X"]
        n, G = X.shape
        o = h["obs"]
        gt = col(o, "gene_transcript")
        L = col(o, "ncounts").astype(np.float64)
        genes = col(h["var"], "gene_name")
        rows, inv = np.unique(gt, return_inverse=True)
        R = len(rows)
        seed = MASTER_SEED + 15000 + DS_ID[ds]
        half = half_labels(inv, R, seed)
        ncell = np.bincount(inv, minlength=R)
        nA = np.bincount(inv[half == 0], minlength=R)
        nB = ncell - nA
        acc = {k: np.zeros((R, G)) for k in ("full", "A", "B")}
        for s in range(0, n, BLOCK):
            e = min(n, s + BLOCK)
            x = np.log1p(X[s:e].astype(np.float64) * (1e4 / L[s:e])[:, None])
            ii = inv[s:e]
            for k, m in (("full", np.ones(e - s, bool)), ("A", half[s:e] == 0), ("B", half[s:e] == 1)):
                ind = sps.csr_matrix((np.ones(int(m.sum())), (ii[m], np.where(m)[0])), shape=(R, e - s))
                acc[k] += ind @ x
    assert (nA > 0).all() or True
    pb = acc["full"] / ncell[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        pA = acc["A"] / nA[:, None]
        pB = acc["B"] / nB[:, None]
    recon = np.nanmax(np.abs((np.nan_to_num(pA) * nA[:, None] + np.nan_to_num(pB) * nB[:, None]) / ncell[:, None] - pb))
    save_npy_new(pb, pdir / "pseudobulk_logmean.npy")
    save_npy_new(pA, pdir / "pseudobulk_logmean_A.npy")
    save_npy_new(pB, pdir / "pseudobulk_logmean_B.npy")
    create_new(pdir / "pseudobulk_rows.csv", pd.DataFrame({"gene_transcript": rows, "n_cells": ncell, "n_A": nA,
                                                          "n_B": nB}).to_csv(index=False))
    create_new(pdir / "pseudobulk_genes.csv", pd.DataFrame({"gene_name": genes}).to_csv(index=False))
    # unit table: data_v07.build_units rules with obs-aggregated metadata (protocol 3.4)
    ob = pd.DataFrame({"gene_transcript": rows, "ncells": ncell})
    U, info = build_units_sc(ob, genes, pdir / "coords_ensembl_extra.csv", ds)
    create_new(pdir / "units_all_targets.csv", U.to_csv(index=False))
    # read-only splits
    Uu = U[U.E | U.always_train].reset_index(drop=True)
    e_idx = np.where(Uu.E.values)[0]
    sbase = MASTER_SEED + 16000 + DS_ID[ds]
    sp = repeated_group_kfold(Uu.tss_group.values[e_idx], sbase)
    sdir.mkdir(parents=True, exist_ok=True)
    for i, s in enumerate(sp):
        rec = {"split": i, "round": s["round"], "fold": s["fold"], "seed": s["seed"],
               **{k: [str(t) for t in Uu.target.values[e_idx[s[k]]]] for k in ("train", "val", "cal", "test")}}
        f = create_new(sdir / f"split_{i:02d}.json", json.dumps(rec))
        os.chmod(f, stat.S_IREAD)
    info.update({"cells": int(n), "genes": int(G), "rows": int(R), "library_size": "obs ncounts (v1.2.1 section 5)",
                 "half_seed": seed, "half_reconstruction_max_abs_diff": float(recon), "split_seed_base": sbase,
                 "n_splits": len(sp), "test_size_range": [min(len(s["test"]) for s in sp), max(len(s["test"]) for s in sp)],
                 "cal_size": len(sp[0]["cal"]), "val_size": len(sp[0]["val"]), "seconds": round(time.time() - t0, 1)})
    create_new(pdir / "c0_build.json", json.dumps(info, indent=1, default=str))
    odir.mkdir(parents=True, exist_ok=True)
    log(ds, "C0 " + json.dumps({k: info[k] for k in ("E", "E_tss_groups", "always_train", "cells", "rows", "coord_source",
                                                     "half_reconstruction_max_abs_diff", "test_size_range")}, default=str))


def build_units_sc(ob, readout, coords_cache, ds):
    """data_v07.build_units (frozen rules) with ob / readout from the single-cell file instead of a bulk file."""
    orig = dv.read_bulk_meta
    dv.read_bulk_meta = lambda _p: (ob.copy(), np.asarray(readout))
    try:
        U, _r, info = dv.build_units(None, coords_cache, (), lambda m: log(ds, m))
    finally:
        dv.read_bulk_meta = orig
    return U, info


# ------------------------------------------------------------------ data
def load(ds):
    pdir, sdir, _ = dirs(ds)
    U = pd.read_csv(pdir / "units_all_targets.csv")
    U = U[U.E | U.always_train].reset_index(drop=True)
    rows = pd.read_csv(pdir / "pseudobulk_rows.csv")
    genes = pd.read_csv(pdir / "pseudobulk_genes.csv").gene_name.values
    rix = {g: i for i, g in enumerate(rows.gene_transcript)}
    U["row"] = U.gene_transcript.map(rix)
    assert U.row.notna().all()
    pb = np.load(pdir / "pseudobulk_logmean.npy", mmap_mode="r")
    ctrl = np.where(rows.gene_transcript.str.contains("non-targeting", case=False))[0]
    w = rows.n_cells.values[ctrl].astype(float)
    y_ctrl = (np.asarray(pb[ctrl]) * w[:, None]).sum(0) / w.sum()
    r = U.row.values.astype(int)
    Y = np.asarray(pb[r])
    YA = np.asarray(np.load(pdir / "pseudobulk_logmean_A.npy", mmap_mode="r")[r])
    YB = np.asarray(np.load(pdir / "pseudobulk_logmean_B.npy", mmap_mode="r")[r])
    gcols = {}
    for j, g in enumerate(genes):
        gcols.setdefault(g, []).append(j)
    U["tcols"] = [np.array(gcols.get(t, []), dtype=int) for t in U.target]
    return U, Y, YA, YB, y_ctrl


def compute_split(ds, i, U, Y, YA, YB, y_ctrl):
    sdir = dirs(ds)[1]
    rec = json.loads((sdir / f"split_{i:02d}.json").read_text())
    pos = {t: j for j, t in enumerate(U.target.values)}
    ix = {k: np.array([pos[t] for t in rec[k]], dtype=int) for k in ("train", "val", "cal", "test")}
    pool = np.concatenate([np.where(U.always_train.values)[0], ix["train"]])
    held = np.concatenate([ix["val"], ix["cal"], ix["test"]])
    assert not set(held) & set(pool), "held-out unit in training pool"
    assert U.E.values[held].all()
    held_groups = set(U.tss_group.values[held])
    assert not (set(U.tss_group.values[ix["train"]]) & set(U.tss_group.values[ix["test"]])), "TSS group shared train/test"
    ref_rows = pool[U.ncells.values[pool] >= 30]
    ref = train_reference(Y[ref_rows])
    mu = fit_mean(Y[ref_rows])(1)[0]
    tc = lambda idx: [U.tcols.values[j] for j in idx]  # noqa: E731
    G = gene_pca_embedding(np.vstack([Y[pool], y_ctrl[None, :]]), K=max(K_GRID), method="gram")
    fit = pool[np.array([len(U.tcols.values[j]) > 0 for j in pool])]
    preds = {"mean": np.tile(mu, (len(held), 1))}
    pcaK = {}
    for K in K_GRID:
        GK = G[:, :K]
        P = np.vstack([GK[U.tcols.values[j][0]] for j in fit] + [np.zeros(K)])
        f = fit_linear(np.vstack([Y[fit], y_ctrl]), P, GK, LAM)
        pcaK[K] = f(np.vstack([GK[U.tcols.values[j][0]] for j in held]))
    nv_, nc = len(ix["val"]), len(ix["cal"])
    errK = {K: float(score_S1(Y[ix["val"]], pcaK[K][:nv_], ref, tc(ix["val"]), K_SCORE).mean()) for K in K_GRID}
    Kc = _pick(errK, K_GRID)
    preds["linear_pca"] = pcaK[Kc]
    ct = np.concatenate([ix["cal"], ix["test"]])
    sl = slice(nv_, None)
    S = {m: score_S1(Y[ct], preds[m][sl], ref, tc(ct), K_SCORE) for m in MODELS}
    q = {m: qhat(S[m][:nc], ALPHA) for m in MODELS}
    eall = np.where(U.E.values)[0]
    E_all = score_S1(Y[eall], np.tile(mu, (len(eall), 1)), ref, tc(eall), K_SCORE)
    tier_all = tiers_of(E_all, 5)
    ep = {j: k for k, j in enumerate(eall)}
    te = ix["test"]
    rows = pd.DataFrame({"split": i, "unit": U.target.values[te], "tss_group": U.tss_group.values[te],
                         "n_cells": U.ncells.values[te], "E": E_all[[ep[j] for j in te]], "tier": tier_all[[ep[j] for j in te]]})
    assert np.allclose(rows.E.values, S["mean"][nc:]), "mean-model score must equal E_true"
    for m in MODELS:
        rows[f"S_{m}"] = S[m][nc:]
        rows[f"q_{m}"] = q[m]
        rows[f"cov_{m}"] = S[m][nc:] <= q[m]
    ht = slice(nv_ + nc, None)
    rows["E_B"] = score_S1(YB[te], np.tile(mu, (len(te), 1)), ref, tc(te), K_SCORE)
    rows["S_ctrl_B"] = score_S1(YB[te], YA[te], ref, tc(te), K_SCORE)
    rows["S_linear_pca_B"] = score_S1(YB[te], preds["linear_pca"][ht], ref, tc(te), K_SCORE)
    permE = {}
    for r in range(N_PERM):
        p = nv.derange(len(eall), np.random.default_rng([nv.SEED_PERM_E, int(i), r]))
        Ep = E_all[p]
        tp = tiers_of(Ep, 5)
        permE[r] = (Ep[[ep[j] for j in te]], tp[[ep[j] for j in te]])
    permP = {}
    Pt = preds["linear_pca"][ht]
    for r in range(N_PERM):
        p = nv.derange(len(te), np.random.default_rng([nv.SEED_PERM_PRED, int(i), r]))
        permP[r] = {"linear_pca": score_S1(Y[te], Pt[p], ref, tc(te), K_SCORE)}
    return {"split": i, "K": Kc, "errK": errK, "rows": rows, "permE": permE, "permP": permP,
            "n": {k: len(v) for k, v in ix.items()}, "n_pool": len(pool), "n_fit_pca": len(fit),
            "n_held_groups": len(held_groups)}


def run(ds):
    od = dirs(ds)[2]
    U, Y, YA, YB, y_ctrl = load(ds)
    cd = od / "cache"
    cd.mkdir(parents=True, exist_ok=True)
    for i in range(21):
        f = cd / f"split_{i:02d}.pkl"
        if rs.is_done(f):
            continue
        t0 = time.time()
        c = compute_split(ds, i, U, Y, YA, YB, y_ctrl)
        rs.pickle_atomic(c, f)
        rs.mark_done(f, time.time() - t0)
        log(ds, f"split {i:02d} done K={c['K']} test={c['n']['test']} cal={c['n']['cal']} {time.time() - t0:.1f}s")


def summary(ds):
    pdir, _, od = dirs(ds)
    files = [od / "cache" / f"split_{i:02d}.pkl" for i in range(21)]
    assert all(rs.is_done(f) for f in files), "not all splits complete"
    cs = [pickle.load(open(f, "rb")) for f in files]
    D = pd.concat([c["rows"] for c in cs], ignore_index=True)
    D.to_parquet(od / "test_rows.parquet")
    cl = D.tss_group.values
    W = nv.boot_weights(len(np.unique(cl)))
    res = {"protocol": "combo_crosscell_v1.2 (+v1.2.1)", "dataset": f"Nadig 2025 {ds} (scPerturb NadigOConner2024_{ds})",
           "n_test_rows": int(len(D)), "n_clusters_tss": int(len(np.unique(cl))), "K_chosen": [c["K"] for c in cs],
           "n_per_split": cs[0]["n"], "thresholds": {"R_lo": nv.R_LO, "R_hi": nv.R_HI}, "models": {}}
    for m in MODELS:
        st = nv.stats_RCV(D[f"S_{m}"].values, D.E.values, cl, W) if m != "mean" else {"R": 0.0, "CV": 0.0, "RER": 0.0}
        if m != "mean":
            t5 = D.tier.values == 5
            st5 = nv.stats_RCV(D[f"S_{m}"].values[t5], D.E.values[t5], cl[t5], nv.boot_weights(len(np.unique(cl[t5]))))
            st["RER_T5"], st["RER_T5_CI"] = st5["RER"], st5["RER_CI"]
            st["verdict"] = nv.judge(st)
        st["coverage_boot"] = nv.cov_stats(D[f"cov_{m}"].values, D.tier.values, cl, W)
        st["coverage_gee"] = nv.gee_cov(D[f"cov_{m}"].values, D.tier.values, cl)
        st["qhat_mean"] = float(D.groupby("split")[f"q_{m}"].first().mean())
        st["rel_width_median_by_tier"] = {str(t): float((D[f"q_{m}"] / D.E)[D.tier == t].median()) for t in range(1, 6)}
        res["models"][m] = st
    df = pd.DataFrame({"cl": cl, "t5": (D.tier == 5).astype(float)})
    df["d"] = (D["cov_linear_pca"].astype(float) - D["cov_mean"].astype(float)) * df.t5
    M, _ = nv.cluster_sums(df, "cl", ["d", "t5"])
    b = (W @ M[:, 0]) / (W @ M[:, 1])
    res["models"]["linear_pca"]["dT5_vs_mean"] = {"est": float(M[:, 0].sum() / M[:, 1].sum()),
                                                  "CI": [float(np.quantile(b, .025)), float(np.quantile(b, .975))]}
    nc1 = {}
    for m in MODELS:
        Rs, gaps = [], []
        for r in range(N_PERM):
            Ep = np.concatenate([c["permE"][r][0] for c in cs])
            tp = np.concatenate([c["permE"][r][1] for c in cs])
            Rs.append(float(1 - np.corrcoef(np.log(Ep), np.log(D[f"S_{m}"].values))[0, 1] ** 2))
            cv = D[f"cov_{m}"].values
            gaps.append(float(cv[tp == 5].mean() - cv[tp == 1].mean()))
        nc1[m] = {"R_mean": float(np.mean(Rs)), "R_range": [min(Rs), max(Rs)], "T5_minus_T1_mean": float(np.mean(gaps)),
                  "T5_minus_T1_range": [min(gaps), max(gaps)], "pass_a": bool(-0.10 <= float(np.mean(gaps)) <= 0.10)}
    nc1["linear_pca"]["pass_b"] = bool(abs(nc1["linear_pca"]["R_mean"] - nc1["mean"]["R_mean"]) < 0.05)
    nc1["pass"] = bool(nc1["mean"]["pass_a"] and nc1["linear_pca"]["pass_a"] and nc1["linear_pca"]["pass_b"])
    res["negative_control_permuted_E"] = nc1
    vs, Rs, RERs = [], [], []
    for r in range(N_PERM):
        Sp = np.concatenate([c["permP"][r]["linear_pca"] for c in cs])
        st = nv.stats_RCV(Sp, D.E.values, cl, W)
        vs.append(nv.judge(st)); Rs.append(st["R"]); RERs.append(st["RER"])
    nc2 = {"linear_pca": {"verdict_counts": {v: vs.count(v) for v in sorted(set(vs))}, "R_mean": float(np.mean(Rs)),
                          "R_range": [min(Rs), max(Rs)], "RER_mean": float(np.mean(RERs)), "pass": vs.count("R clearly elevated") == 0}}
    res["negative_control_permuted_predictions"] = nc2
    pc = {}
    for m, c_ in (("split_half_control", "S_ctrl_B"), ("linear_pca", "S_linear_pca_B")):
        st = nv.stats_RCV(D[c_].values, D.E_B.values, cl, W)
        st["verdict_rule"] = nv.judge(st)
        pc[m] = st
    pc["threshold_reachable"] = bool(pc["split_half_control"]["R"] >= nv.R_HI and pc["split_half_control"]["RER_CI"][0] > 0)
    res["positive_control_half_basis"] = pc
    v = res["models"]["linear_pca"]["verdict"]
    res["judgement"] = {"primary_linear_pca": v, "prediction": "R not elevated", "prediction_holds": v == "R not elevated",
                        "negative_control_1_pass": nc1["pass"], "negative_control_2_pass": nc2["linear_pca"]["pass"],
                        "R_verdicts_trusted": bool(nc1["pass"])}
    if v == "R not elevated":
        res["judgement"]["side_reading_T5_vs_mean"] = res["models"]["linear_pca"]["dT5_vs_mean"]
    create_new(od / "nadig_v12_result.json", json.dumps(res, indent=1, default=float))
    rs.mark_done(od / "nadig_v12_result.json", 0.0)
    print(json.dumps(res["judgement"], default=float))


if __name__ == "__main__":
    {"c0": c0, "run": run, "summary": summary}[sys.argv[1]](sys.argv[2])
