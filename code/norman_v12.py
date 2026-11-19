# norman_v12.py -- protocol combo_crosscell_v1.2 on Norman 2019 (combinatorial CRISPRa, K562).
#   python norman_v12.py c0      unit table, readout genes, E units (section 2.1, 3.1); no score is computed
#   python norman_v12.py splits  21 unseen-combination splits (section 3.2), written read-only
#   python norman_v12.py run     per-split models, scores, conformal sets and controls (resumable: one cache file
#                                per split with a sha256 completion marker; completed splits are skipped)
#   python norman_v12.py summary pooled statistics, cluster bootstrap, controls and the pre-specified judgement
# Frozen functions reused unchanged: models_v04 (train_reference, fit_mean, gene_pca_embedding, fit_linear,
# score_S1), methods_v06.qhat, eval_v06 (tiers_of, _seed), splits_v06.repeated_group_kfold, pipeline_v05._pick,
# resume_v07 (atomic writes and markers).
import json
import os
import stat
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED, RESULTS, SPLITS  # noqa: E402
from models_v04 import fit_linear, fit_mean, gene_pca_embedding, score_S1, train_reference  # noqa: E402
from methods_v06 import qhat  # noqa: E402
from eval_v06 import _seed, tiers_of  # noqa: E402
from splits_v06 import repeated_group_kfold  # noqa: E402
from pipeline_v05 import _pick  # noqa: E402
import resume_v07 as rs  # noqa: E402

ND = PROCESSED / "v1.2" / "norman"
SD = SPLITS / "v1.2" / "norman"
OD = RESULTS / "v1.2" / "norman"
K_GRID = (5, 10, 20, 30, 50)
LAM = 0.1
K_SCORE = 50
ALPHA = 0.10
N_PERM = 20
N_BOOT = 2000
READ_THR = 0.1
R_LO, R_HI = 0.112, 0.224
SEED_SPLIT = MASTER_SEED + 14000
SEED_PERM_E = MASTER_SEED + 13000
SEED_PERM_PRED = MASTER_SEED + 12500
SEED_BOOT = MASTER_SEED + 17000
MODELS = ("mean", "additive", "linear_pca")


def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with open(OD / "run.log", "a", encoding="utf8") as fh:
        fh.write(line + "\n")


# ------------------------------------------------------------------ data
def load():
    U = pd.read_csv(ND / "units_c0.csv")
    read = np.load(ND / "readout_idx.npy")
    F = np.load(ND / "pb_full.npy", mmap_mode="r")
    A = np.load(ND / "pb_A.npy", mmap_mode="r")
    B = np.load(ND / "pb_B.npy", mmap_mode="r")
    Y, YA, YB = (np.ascontiguousarray(M[:, read]) for M in (F, A, B))
    ci = int(np.where(U.unit.values == "control")[0][0])
    U["tcols"] = [np.array(json.loads(t), dtype=int) for t in U.tcols_json]
    return U, Y, YA, YB, ci


def c0():
    OD.mkdir(parents=True, exist_ok=True)
    U = pd.read_csv(ND / "units.csv")
    genes = pd.read_csv(ND / "genes.csv").gene.values
    cm = np.load(ND / "ctrl_mean_counts.npy")
    gi = {g: i for i, g in enumerate(genes)}
    keep = cm >= READ_THR
    targets = sorted({g for u in U.unit if u != "control" for g in u.split("_")})
    missing = [t for t in targets if t not in gi]
    for t in targets:
        if t in gi:
            keep[gi[t]] = True
    read = np.where(keep)[0]
    col = {genes[j]: c for c, j in enumerate(read)}
    tc = []
    for u in U.unit:
        tc.append([] if u == "control" else [col[g] for g in u.split("_") if g in col])
    U["tcols_json"] = [json.dumps(t) for t in tc]
    U["all_targets_readout"] = [u != "control" and len(t) == len(u.split("_")) for u, t in zip(U.unit, tc)]
    singles = set(U.unit[(U.n_genes_in_unit == 1) & (U.n_cells >= 30)])
    U["E"] = [(n == 2) and (c >= 30) and all(g in singles for g in u.split("_")) and ok
              for u, n, c, ok in zip(U.unit, U.n_genes_in_unit, U.n_cells, U.all_targets_readout)]
    U["pool_always"] = (U.unit != "control") & (U.n_cells >= 30) & ~U.E
    rs.write_bytes_atomic(ND / "units_c0.csv", U.to_csv(index=False).encode("utf8"))
    with open(ND / "readout_idx.npy.part", "wb") as fh:
        np.save(fh, read)
    os.replace(ND / "readout_idx.npy.part", ND / "readout_idx.npy")
    info = {"n_units": int(len(U)), "n_singles": int((U.n_genes_in_unit == 1).sum()),
            "n_doubles": int((U.n_genes_in_unit == 2).sum()), "n_E": int(U.E.sum()),
            "n_always_train": int(U.pool_always.sum()), "n_ctrl_cells": int(U.n_cells[U.unit == "control"].iloc[0]),
            "n_readout_threshold": int((cm >= READ_THR).sum()), "n_readout": int(len(read)),
            "targets": len(targets), "targets_missing_from_file": missing,
            "doubles_excluded_missing_target": [u for u, n, ok in zip(U.unit, U.n_genes_in_unit, U.all_targets_readout)
                                                if n == 2 and not ok],
            "E_cells_median": float(U.n_cells[U.E].median()), "E_cells_min": int(U.n_cells[U.E].min()),
            "underpowered_flag": bool(U.E.sum() < 60)}
    rs.write_bytes_atomic(OD / "c0_summary.json", json.dumps(info, indent=1).encode("utf8"))
    print(json.dumps(info))


def splits():
    U = pd.read_csv(ND / "units_c0.csv")
    e_idx = np.where(U.E.values)[0]
    sp = repeated_group_kfold(U.unit.values[e_idx], SEED_SPLIT)
    SD.mkdir(parents=True, exist_ok=True)
    for i, s in enumerate(sp):
        f = SD / f"split_{i:02d}.json"
        rec = {"split": i, "round": s["round"], "fold": s["fold"], "seed": s["seed"],
               **{k: U.unit.values[e_idx[s[k]]].tolist() for k in ("train", "val", "cal", "test")}}
        if f.exists():
            assert json.loads(f.read_text()) == rec, f"existing split {i} differs"
            continue
        f.write_text(json.dumps(rec))
        os.chmod(f, stat.S_IREAD)
    print(len(sp), "splits;", "test sizes", sorted({len(s["test"]) for s in sp}), "cal", sorted({len(s["cal"]) for s in sp}),
          "val", sorted({len(s["val"]) for s in sp}))


def derange(n, rng):
    if n < 2:
        raise ValueError("cannot derange fewer than 2 items")
    while True:
        p = rng.permutation(n)
        if not (p == np.arange(n)).any():
            return p


# ------------------------------------------------------------------ one split
def compute_split(i, U, Y, YA, YB, ci):
    rec = json.loads((SD / f"split_{i:02d}.json").read_text())
    pos = {u: j for j, u in enumerate(U.unit.values)}
    ix = {k: np.array([pos[u] for u in rec[k]]) for k in ("train", "val", "cal", "test")}
    pool = np.concatenate([np.where(U.pool_always.values)[0], ix["train"]])
    held = np.concatenate([ix["val"], ix["cal"], ix["test"]])
    # leakage asserts
    assert not set(held) & set(pool), "held-out combination in training pool"
    assert len(set(ix["val"]) | set(ix["cal"]) | set(ix["test"])) == len(held)
    assert U.E.values[held].all() and (U.n_genes_in_unit.values[pool] >= 1).all()
    y_ctrl = Y[ci]
    ref_rows = pool[U.n_cells.values[pool] >= 30]
    ref = train_reference(Y[ref_rows])
    mu = fit_mean(Y[ref_rows])(1)[0]
    tc = lambda idx: [U.tcols.values[j] for j in idx]
    # additive: singles of each held-out double come from the training pool
    sidx = lambda u: pos[u]
    for j in held:
        for g in U.unit.values[j].split("_"):
            assert sidx(g) in set(pool), "single of a held-out double missing from the pool"
    add = np.vstack([Y[sidx(a)] + Y[sidx(b)] - y_ctrl for a, b in (U.unit.values[j].split("_") for j in held)])
    # linear-PCA with additive embedding of doubles (protocol 3.3 item 3)
    G = gene_pca_embedding(np.vstack([Y[pool], y_ctrl[None, :]]), K=max(K_GRID), method="gram")
    fit = pool[U.all_targets_readout.values[pool]]
    preds = {"mean": np.tile(mu, (len(held), 1)), "additive": add}
    pcaK = {}
    for K in K_GRID:
        GK = G[:, :K]
        P = np.vstack([GK[U.tcols.values[j]].sum(0) for j in fit] + [np.zeros(K)])
        f = fit_linear(np.vstack([Y[fit], y_ctrl]), P, GK, LAM)
        pcaK[K] = f(np.vstack([GK[U.tcols.values[j]].sum(0) for j in held]))
    nv, nc = len(ix["val"]), len(ix["cal"])
    errK = {K: float(score_S1(Y[ix["val"]], pcaK[K][:nv], ref, tc(ix["val"]), K_SCORE).mean()) for K in K_GRID}
    Kc = _pick(errK, K_GRID)
    preds["linear_pca"] = pcaK[Kc]
    ct = np.concatenate([ix["cal"], ix["test"]])
    sl = slice(nv, None)
    S = {m: score_S1(Y[ct], preds[m][sl], ref, tc(ct), K_SCORE) for m in MODELS}
    q = {m: qhat(S[m][:nc], ALPHA) for m in MODELS}
    # effect strength and tiers over all E units of the split
    eall = np.where(U.E.values)[0]
    E_all = score_S1(Y[eall], np.tile(mu, (len(eall), 1)), ref, tc(eall), K_SCORE)
    tier_all = tiers_of(E_all, 5)
    ep = {j: k for k, j in enumerate(eall)}
    te = ix["test"]
    rows = pd.DataFrame({"split": i, "unit": U.unit.values[te], "n_cells": U.n_cells.values[te],
                         "E": E_all[[ep[j] for j in te]], "tier": tier_all[[ep[j] for j in te]]})
    assert np.allclose(rows.E.values, S["mean"][nc:]), "mean-model score must equal E_true"
    for m in MODELS:
        rows[f"S_{m}"] = S[m][nc:]
        rows[f"q_{m}"] = q[m]
        rows[f"cov_{m}"] = S[m][nc:] <= q[m]
    # half-data version (truth = half B; control predictor = half A)
    ht = slice(nv + nc, None)
    rows["E_B"] = score_S1(YB[te], np.tile(mu, (len(te), 1)), ref, tc(te), K_SCORE)
    rows["S_ctrl_B"] = score_S1(YB[te], YA[te], ref, tc(te), K_SCORE)
    for m in ("additive", "linear_pca"):
        rows[f"S_{m}_B"] = score_S1(YB[te], preds[m][ht], ref, tc(te), K_SCORE)
    # negative control 1: permuted effect strength (E over the split's E units), re-tiered
    permE = {}
    for r in range(N_PERM):
        p = derange(len(eall), np.random.default_rng([SEED_PERM_E, int(i), r]))
        Ep = E_all[p]
        tp = tiers_of(Ep, 5)
        permE[r] = (Ep[[ep[j] for j in te]], tp[[ep[j] for j in te]])
    # negative control 2: predictions permuted among test units
    permP = {}
    Pt = {m: preds[m][ht] for m in ("additive", "linear_pca")}
    for r in range(N_PERM):
        p = derange(len(te), np.random.default_rng([SEED_PERM_PRED, int(i), r]))
        permP[r] = {m: score_S1(Y[te], Pt[m][p], ref, tc(te), K_SCORE) for m in Pt}
    return {"split": i, "K": Kc, "errK": errK, "rows": rows, "permE": permE, "permP": permP,
            "n": {k: len(v) for k, v in ix.items()}, "n_pool": len(pool), "n_fit_pca": len(fit)}


def run():
    OD.mkdir(parents=True, exist_ok=True)
    U, Y, YA, YB, ci = load()
    cd = OD / "cache"
    cd.mkdir(exist_ok=True)
    for i in range(21):
        f = cd / f"split_{i:02d}.pkl"
        if rs.is_done(f):
            continue
        t0 = time.time()
        c = compute_split(i, U, Y, YA, YB, ci)
        rs.pickle_atomic(c, f)
        rs.mark_done(f, time.time() - t0)
        log(f"split {i:02d} done K={c['K']} test={c['n']['test']} cal={c['n']['cal']} {time.time() - t0:.1f}s")


# ------------------------------------------------------------------ statistics
def cluster_sums(df, cl, cols):
    g = df.groupby(cl, sort=True)
    return g[cols].sum().values, g.ngroups


def boot_weights(n_cl):
    rng = np.random.default_rng(SEED_BOOT)
    return rng.multinomial(n_cl, np.full(n_cl, 1.0 / n_cl), size=N_BOOT).astype(float)


def stats_RCV(S, E, cl, W):
    """R, CV and RER (point and bootstrap) for score S against effect strength E; clusters cl."""
    x, y, r = np.log(E), np.log(S), S / E
    df = pd.DataFrame({"cl": cl, "n": 1.0, "x": x, "y": y, "xx": x * x, "yy": y * y, "xy": x * y, "r": r, "rr": r * r,
                       "S": S, "E": E})
    M, _ = cluster_sums(df, "cl", ["n", "x", "y", "xx", "yy", "xy", "r", "rr", "S", "E"])
    def f(w):
        t = w @ M if w.ndim == 2 else M.sum(0)
        n, sx, sy, sxx, syy, sxy, sr, srr, sS, sE = (t[..., k] for k in range(10))
        cxy = sxy / n - sx / n * sy / n
        vx = sxx / n - (sx / n) ** 2
        vy = syy / n - (sy / n) ** 2
        R = 1 - cxy ** 2 / (vx * vy)
        mr = sr / n
        CV = np.sqrt(np.maximum(srr / n - mr ** 2, 0) * n / (n - 1)) / mr
        RER = 1 - sS / sE
        return R, CV, RER
    pt = f(np.ones(0))
    bt = f(W)
    out = {}
    for k, nm in enumerate(("R", "CV", "RER")):
        out[nm] = float(pt[k])
        out[nm + "_CI"] = [float(np.quantile(bt[k], 0.025)), float(np.quantile(bt[k], 0.975))]
    return out


def judge(st):
    acc = st["RER_CI"][0] > 0
    if acc and st["R"] >= R_HI and st["R_CI"][0] > R_LO:
        return "R clearly elevated"
    if st["R"] >= R_HI and not acc:
        return "R elevated without accuracy gain"
    if st["R_CI"][1] < R_HI:
        return "R not elevated"
    return "uncertain"


def cov_stats(cov, tier, cl, W, ref_cov=None):
    df = pd.DataFrame({"cl": cl, "c": cov.astype(float), "one": 1.0})
    out = {}
    for t in ("all", 1, 2, 3, 4, 5):
        m = np.ones(len(df), bool) if t == "all" else tier == t
        d = df.assign(c=df.c * m, one=df.one * m)
        M, _ = cluster_sums(d, "cl", ["c", "one"])
        pt = M[:, 0].sum() / M[:, 1].sum()
        b = (W @ M[:, 0]) / (W @ M[:, 1])
        out[str(t)] = {"cov": float(pt), "CI": [float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))], "n": int(m.sum())}
    return out


def gee_cov(cov, tier, cl):
    import statsmodels.api as sm
    from scipy.special import expit
    out = {}
    for t in ("all", 1, 2, 3, 4, 5):
        m = np.ones(len(cov), bool) if t == "all" else tier == t
        y = cov[m].astype(float)
        if y.min() == y.max():
            out[str(t)] = {"cov": float(y.mean()), "CI": None, "note": "all rows identical; GEE undefined"}
            continue
        g = sm.GEE(y, np.ones((m.sum(), 1)), groups=pd.factorize(cl[m])[0], family=sm.families.Binomial(),
                   cov_struct=sm.cov_struct.Independence()).fit()
        b, se = float(g.params[0]), float(g.bse[0])
        out[str(t)] = {"cov": float(expit(b)), "CI": [float(expit(b - 1.96 * se)), float(expit(b + 1.96 * se))]}
    return out


def summary():
    files = [OD / "cache" / f"split_{i:02d}.pkl" for i in range(21)]
    assert all(rs.is_done(f) for f in files), "not all splits complete"
    import pickle
    cs = [pickle.load(open(f, "rb")) for f in files]
    D = pd.concat([c["rows"] for c in cs], ignore_index=True)
    D.to_parquet(OD / "test_rows.parquet")
    cl = D.unit.values
    ucl = np.unique(cl)
    W = boot_weights(len(ucl))
    res = {"protocol": "combo_crosscell_v1.2", "dataset": "Norman 2019 (scPerturb NormanWeissman2019_filtered)",
           "n_test_rows": int(len(D)), "n_clusters": int(len(ucl)), "K_chosen": [c["K"] for c in cs],
           "n_per_split": cs[0]["n"], "thresholds": {"R_lo": R_LO, "R_hi": R_HI}}
    res["models"] = {}
    for m in MODELS:
        st = stats_RCV(D[f"S_{m}"].values, D.E.values, cl, W) if m != "mean" else {"R": 0.0, "CV": 0.0, "RER": 0.0}
        t5 = D.tier.values == 5
        if m != "mean":
            st5 = stats_RCV(D[f"S_{m}"].values[t5], D.E.values[t5], cl[t5],
                            boot_weights(len(np.unique(cl[t5]))))
            st["RER_T5"], st["RER_T5_CI"] = st5["RER"], st5["RER_CI"]
            st["verdict"] = judge(st)
        st["coverage_boot"] = cov_stats(D[f"cov_{m}"].values, D.tier.values, cl, W)
        st["coverage_gee"] = gee_cov(D[f"cov_{m}"].values, D.tier.values, cl)
        st["qhat_mean"] = float(D.groupby("split")[f"q_{m}"].first().mean())
        st["rel_width_median_by_tier"] = {str(t): float((D[f"q_{m}"] / D.E)[D.tier == t].median()) for t in range(1, 6)}
        res["models"][m] = st
    # Delta T5 vs mean model (paired, same rows, same bootstrap)
    df = pd.DataFrame({"cl": cl, "t5": (D.tier == 5).astype(float)})
    for m in ("additive", "linear_pca"):
        df["d"] = (D[f"cov_{m}"].astype(float) - D["cov_mean"].astype(float)) * df.t5
        M, _ = cluster_sums(df, "cl", ["d", "t5"])
        pt = M[:, 0].sum() / M[:, 1].sum()
        b = (W @ M[:, 0]) / (W @ M[:, 1])
        res["models"][m]["dT5_vs_mean"] = {"est": float(pt), "CI": [float(np.quantile(b, .025)), float(np.quantile(b, .975))]}
    # negative control 1: permuted effect strength
    nc1 = {}
    for m in MODELS:
        Rs, gaps = [], []
        for r in range(N_PERM):
            Ep = np.concatenate([c["permE"][r][0] for c in cs])
            tp = np.concatenate([c["permE"][r][1] for c in cs])
            x, y = np.log(Ep), np.log(D[f"S_{m}"].values)
            Rs.append(float(1 - np.corrcoef(x, y)[0, 1] ** 2))
            cv = D[f"cov_{m}"].values
            gaps.append(float(cv[tp == 5].mean() - cv[tp == 1].mean()))
        nc1[m] = {"R_mean": float(np.mean(Rs)), "R_range": [min(Rs), max(Rs)], "T5_minus_T1_mean": float(np.mean(gaps)),
                  "T5_minus_T1_range": [min(gaps), max(gaps)]}
    for m in ("additive", "linear_pca"):
        nc1[m]["pass_a"] = bool(-0.10 <= nc1[m]["T5_minus_T1_mean"] <= 0.10)
        nc1[m]["pass_b"] = bool(abs(nc1[m]["R_mean"] - nc1["mean"]["R_mean"]) < 0.05)
    nc1["mean"]["pass_a"] = bool(-0.10 <= nc1["mean"]["T5_minus_T1_mean"] <= 0.10)
    nc1["pass"] = bool(all(nc1[m]["pass_a"] for m in MODELS) and all(nc1[m]["pass_b"] for m in ("additive", "linear_pca")))
    res["negative_control_permuted_E"] = nc1
    # negative control 2: permuted predictions
    nc2 = {}
    for m in ("additive", "linear_pca"):
        vs, Rs, RERs = [], [], []
        for r in range(N_PERM):
            Sp = np.concatenate([c["permP"][r][m] for c in cs])
            st = stats_RCV(Sp, D.E.values, cl, W)
            vs.append(judge(st)); Rs.append(st["R"]); RERs.append(st["RER"])
        nc2[m] = {"verdict_counts": {v: vs.count(v) for v in sorted(set(vs))}, "R_mean": float(np.mean(Rs)),
                  "R_range": [min(Rs), max(Rs)], "RER_mean": float(np.mean(RERs)),
                  "pass": vs.count("R clearly elevated") == 0}
    res["negative_control_permuted_predictions"] = nc2
    # positive control: split-half predictor, half-data basis
    pc = {}
    for m, col in (("split_half_control", "S_ctrl_B"), ("additive", "S_additive_B"), ("linear_pca", "S_linear_pca_B")):
        st = stats_RCV(D[col].values, D.E_B.values, cl, W)
        st["verdict_rule"] = judge(st)
        pc[m] = st
    pc["threshold_reachable"] = bool(pc["split_half_control"]["R"] >= R_HI and pc["split_half_control"]["RER_CI"][0] > 0)
    res["positive_control_half_basis"] = pc
    # pre-specified judgement
    trust_elev = nc1["pass"] and all(nc2[m]["pass"] for m in nc2)
    res["judgement"] = {"primary_additive": res["models"]["additive"]["verdict"],
                        "secondary_linear_pca": res["models"]["linear_pca"]["verdict"],
                        "negative_control_1_pass": nc1["pass"],
                        "negative_control_2_pass": {m: nc2[m]["pass"] for m in nc2},
                        "R_verdicts_trusted": bool(nc1["pass"]),
                        "elevated_verdict_trusted": bool(trust_elev),
                        "underpowered": bool(json.loads((OD / "c0_summary.json").read_text())["underpowered_flag"])}
    rs.write_bytes_atomic(OD / "norman_v12_result.json", json.dumps(res, indent=1, default=float).encode("utf8"))
    rs.mark_done(OD / "norman_v12_result.json", 0.0)
    print(json.dumps(res["judgement"]))


if __name__ == "__main__":
    {"c0": c0, "splits": splits, "run": run, "summary": summary}[sys.argv[1]]()
