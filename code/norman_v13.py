# norman_v13.py -- protocol norman_ablation_bootstrap_v1.3 (post hoc; frozen before computation).
#   python norman_v13.py ablation   section A: unseen-component ablation (gene-level folds), resumable per split
#   python norman_v13.py summary    section A statistics (dyadic gene bootstrap)
#   python norman_v13.py bootstrap  section B: gene-level resampling of the v1.2 test rows
# Reuses the v1.2 data and the frozen functions unchanged.
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, RESULTS  # noqa: E402
from models_v04 import fit_linear, fit_mean, gene_pca_embedding, score_S1, train_reference  # noqa: E402
import norman_v12 as v12  # noqa: E402
import resume_v07 as rs  # noqa: E402

OD = RESULTS / "v1.3" / "norman"
N_FOLD, N_ROUND = 5, 3
K_MAIN, K_SENS = 50, 10
LAM, K_SCORE = v12.LAM, v12.K_SCORE
N_BOOT = 2000
SEED_FOLD = MASTER_SEED + 19000
SEED_BOOT_A = MASTER_SEED + 19500
SEED_BOOT_TWO = MASTER_SEED + 18000
SEED_BOOT_ONE = MASTER_SEED + 18001
R_LO, R_HI = v12.R_LO, v12.R_HI


def log(msg):
    OD.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with open(OD / "run.log", "a", encoding="utf8") as fh:
        fh.write(line + "\n")


def gene_folds(genes):
    out = []
    for r in range(N_ROUND):
        rng = np.random.default_rng(SEED_FOLD + r)
        perm = rng.permutation(np.asarray(sorted(genes), dtype=object))
        for f, chunk in enumerate(np.array_split(perm, N_FOLD)):
            out.append({"round": r, "fold": f, "held": sorted(chunk.tolist())})
    return out


# ------------------------------------------------------------------ section A
def ablation_split(i, spec, U, Y, ci):
    held = set(spec["held"])
    units = U.unit.values
    genes_of = [[] if u == "control" else u.split("_") for u in units]
    touch = np.array([any(g in held for g in gs) for gs in genes_of])
    pool = np.where((U.E.values | U.pool_always.values) & ~touch)[0]
    test = np.where(U.E.values & touch)[0]
    n_unseen = np.array([sum(g in held for g in genes_of[j]) for j in test])
    # leakage asserts
    assert not any(g in held for j in pool for g in genes_of[j]), "held-out gene in training pool"
    assert not set(pool) & set(test)
    y_ctrl = Y[ci]
    ref_rows = pool[U.n_cells.values[pool] >= 30]
    ref = train_reference(Y[ref_rows])
    mu = fit_mean(Y[ref_rows])(1)[0]
    tc = [U.tcols.values[j] for j in test]
    E = score_S1(Y[test], np.tile(mu, (len(test), 1)), ref, tc, K_SCORE)
    G = gene_pca_embedding(np.vstack([Y[pool], y_ctrl[None, :]]), K=K_MAIN, method="gram")
    fit = pool[U.all_targets_readout.values[pool]]
    S = {}
    for K in (K_MAIN, K_SENS):
        GK = G[:, :K]
        P = np.vstack([GK[U.tcols.values[j]].sum(0) for j in fit] + [np.zeros(K)])
        f = fit_linear(np.vstack([Y[fit], y_ctrl]), P, GK, LAM)
        pred = f(np.vstack([GK[U.tcols.values[j]].sum(0) for j in test]))
        S[f"linear_pca_K{K}"] = score_S1(Y[test], pred, ref, tc, K_SCORE)
    pos = {u: j for j, u in enumerate(units)}
    part = np.full(len(test), np.nan)
    one = np.where(n_unseen == 1)[0]
    if len(one):
        seen = [[g for g in genes_of[test[k]] if g not in held][0] for k in one]
        assert all(pos[g] in set(pool) for g in seen), "seen single missing from pool"
        pp = np.vstack([Y[pos[g]] for g in seen])
        part[one] = score_S1(Y[test[one]], pp, ref, [tc[k] for k in one], K_SCORE)
    rows = pd.DataFrame({"split": i, "round": spec["round"], "fold": spec["fold"], "unit": units[test],
                         "n_unseen": n_unseen, "E": E, "S_partial_additive": part, **{f"S_{k}": v for k, v in S.items()}})
    assert np.isfinite(rows.E).all() and (rows.E > 0).all()
    return {"split": i, "spec": spec, "rows": rows, "n_pool": int(len(pool)), "n_test": int(len(test))}


def ablation():
    U, Y, _, _, ci = v12.load()
    genes = sorted({g for u in U.unit.values[U.E.values] for g in u.split("_")})
    specs = gene_folds(genes)
    OD.mkdir(parents=True, exist_ok=True)
    rs.write_bytes_atomic(OD / "gene_folds.json", json.dumps(specs, indent=1).encode("utf8"))
    cd = OD / "cache"
    cd.mkdir(exist_ok=True)
    for i, spec in enumerate(specs):
        f = cd / f"abl_{i:02d}.pkl"
        if rs.is_done(f):
            continue
        t0 = time.time()
        c = ablation_split(i, spec, U, Y, ci)
        rs.pickle_atomic(c, f)
        rs.mark_done(f, time.time() - t0)
        log(f"ablation split {i:02d} round {spec['round']} fold {spec['fold']} pool={c['n_pool']} test={c['n_test']} "
            f"({(c['rows'].n_unseen == 2).sum()} two-unseen) {time.time() - t0:.1f}s")


# ------------------------------------------------------------------ resampling weights and weighted statistics
def gene_weights(units, seed, mode):
    """Row weights (B x n) from resampling genes. mode 'two': m_A * m_B; 'one': (m_A + m_B) / 2."""
    pairs = [u.split("_") for u in units]
    genes = sorted({g for p in pairs for g in p})
    gi = {g: k for k, g in enumerate(genes)}
    a = np.array([gi[p[0]] for p in pairs])
    b = np.array([gi[p[1]] for p in pairs])
    rng = np.random.default_rng(seed)
    M = rng.multinomial(len(genes), np.full(len(genes), 1.0 / len(genes)), size=N_BOOT).astype(float)
    return M[:, a] * M[:, b] if mode == "two" else (M[:, a] + M[:, b]) / 2.0


def w_stats(S, E, W):
    """Weighted R, CV, RER for each row of W (B x n); W = ones(1, n) gives the point estimate."""
    x, y, r = np.log(E), np.log(S), S / E
    n = W.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mx, my = W @ x / n, W @ y / n
        cxy = W @ (x * y) / n - mx * my
        vx = W @ (x * x) / n - mx ** 2
        vy = W @ (y * y) / n - my ** 2
        R = 1 - cxy ** 2 / (vx * vy)
        mr = W @ r / n
        CV = np.sqrt(np.maximum(W @ (r * r) / n - mr ** 2, 0) * n / (n - 1)) / mr
        RER = 1 - (W @ S) / (W @ E)
    return {"R": R, "CV": CV, "RER": RER}


def ci(v):
    v = v[np.isfinite(v)]
    return [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))], int(len(v))


def judge(st):
    return v12.judge(st)


# ------------------------------------------------------------------ section A summary
def summary():
    files = [OD / "cache" / f"abl_{i:02d}.pkl" for i in range(N_ROUND * N_FOLD)]
    assert all(rs.is_done(f) for f in files), "not all ablation splits complete"
    D = pd.concat([pickle.load(open(f, "rb"))["rows"] for f in files], ignore_index=True)
    D.to_parquet(OD / "ablation_rows.parquet")
    res = {"protocol": "norman_ablation_bootstrap_v1.3", "section": "A", "evidence": "post hoc",
           "n_splits": len(files), "thresholds": {"R_lo": R_LO, "R_hi": R_HI}, "categories": {}}
    for cat in (1, 2):
        d = D[D.n_unseen == cat].reset_index(drop=True)
        W = gene_weights(d.unit.values, SEED_BOOT_A + cat, "two")
        one = np.ones((1, len(d)))
        out = {"n_rows": int(len(d)), "n_pairs": int(d.unit.nunique()),
               "underpowered": bool(cat == 2 and d.unit.nunique() < 40), "models": {}}
        for m in ("linear_pca_K50", "linear_pca_K10", "partial_additive"):
            col = f"S_{m}"
            if d[col].isna().all():
                continue
            pt = w_stats(d[col].values, d.E.values, one)
            bt = w_stats(d[col].values, d.E.values, W)
            st = {}
            for k in ("R", "CV", "RER"):
                st[k] = float(pt[k][0])
                st[k + "_CI"], st[k + "_nvalid"] = ci(bt[k])
            st["verdict_rule"] = judge(st)
            out["models"][m] = st
        res["categories"][f"{cat}_unseen"] = out
    m = res["categories"]["2_unseen"]["models"]["linear_pca_K50"]
    if res["categories"]["2_unseen"]["underpowered"]:
        interp = "underpowered (fewer than 40 distinct two-unseen pairs)"
    elif m["R_CI"][1] < R_HI:
        interp = "returns to screen range"
    elif m["RER_CI"][0] > 0 and m["R"] >= R_HI and m["R_CI"][0] > R_LO:
        interp = "still elevated"
    else:
        interp = "uncertain"
    res["interpretation_A5"] = interp
    rs.write_bytes_atomic(OD / "ablation_result.json", json.dumps(res, indent=1).encode("utf8"))
    rs.mark_done(OD / "ablation_result.json", 0.0)
    print(json.dumps({"interpretation": interp,
                      **{c: {"n_pairs": v["n_pairs"], **{mm: {k: round(s[k], 3) if isinstance(s[k], float) else s[k]
                                                              for k in ("R", "R_CI", "RER", "RER_CI", "verdict_rule")}
                                                         for mm, s in v["models"].items()}}
                         for c, v in res["categories"].items()}}, default=str))


# ------------------------------------------------------------------ section B
def bootstrap():
    D = pd.read_parquet(RESULTS / "v1.2" / "norman" / "test_rows.parquet")
    one = np.ones((1, len(D)))
    Ws = {"two_way": gene_weights(D.unit.values, SEED_BOOT_TWO, "two"),
          "one_way": gene_weights(D.unit.values, SEED_BOOT_ONE, "one")}
    tier = D.tier.values
    res = {"protocol": "norman_ablation_bootstrap_v1.3", "section": "B", "evidence": "post hoc",
           "n_rows": int(len(D)), "n_pairs": int(D.unit.nunique()), "n_genes": int(len({g for u in D.unit for g in u.split("_")})),
           "models": {}}
    for m in v12.MODELS:
        st = {"point": {}, "two_way": {}, "one_way": {}}
        cov = D[f"cov_{m}"].values.astype(float)
        for t in ("all", 1, 2, 3, 4, 5):
            msk = np.ones(len(D), bool) if t == "all" else tier == t
            st["point"][f"cov_{t}"] = float(cov[msk].mean())
            for nm, W in Ws.items():
                v = (W[:, msk] @ cov[msk]) / W[:, msk].sum(1)
                st[nm][f"cov_{t}_CI"], _ = ci(v)
        if m != "mean":
            pt = w_stats(D[f"S_{m}"].values, D.E.values, one)
            for k in ("R", "CV", "RER"):
                st["point"][k] = float(pt[k][0])
            d = (cov - D.cov_mean.values.astype(float)) * (tier == 5)
            st["point"]["dT5"] = float(d.sum() / (tier == 5).sum())
            for nm, W in Ws.items():
                bt = w_stats(D[f"S_{m}"].values, D.E.values, W)
                for k in ("R", "CV", "RER"):
                    st[nm][k + "_CI"], st[nm][k + "_nvalid"] = ci(bt[k])
                v = (W @ d) / (W @ (tier == 5).astype(float))
                st[nm]["dT5_CI"], _ = ci(v)
            j = {"R": st["point"]["R"], "R_CI": st["two_way"]["R_CI"], "RER": st["point"]["RER"], "RER_CI": st["two_way"]["RER_CI"]}
            st["post_hoc_verdict_two_way"] = judge(j)
            j1 = {"R": st["point"]["R"], "R_CI": st["one_way"]["R_CI"], "RER": st["point"]["RER"], "RER_CI": st["one_way"]["RER_CI"]}
            st["post_hoc_verdict_one_way"] = judge(j1)
        res["models"][m] = st
    res["robust_two_way"] = bool(all(res["models"][m]["post_hoc_verdict_two_way"] == "R clearly elevated" for m in ("additive", "linear_pca")))
    rs.write_bytes_atomic(OD / "bootstrap_result.json", json.dumps(res, indent=1).encode("utf8"))
    rs.mark_done(OD / "bootstrap_result.json", 0.0)
    print(json.dumps({"robust_two_way": res["robust_two_way"],
                      **{m: {"two_way_R_CI": [round(x, 3) for x in res["models"][m]["two_way"]["R_CI"]],
                             "one_way_R_CI": [round(x, 3) for x in res["models"][m]["one_way"]["R_CI"]],
                             "two_way_dT5_CI": [round(x, 3) for x in res["models"][m]["two_way"]["dT5_CI"]],
                             "verdict2": res["models"][m]["post_hoc_verdict_two_way"]} for m in ("additive", "linear_pca")}}))


if __name__ == "__main__":
    OD.mkdir(parents=True, exist_ok=True)
    {"ablation": ablation, "summary": summary, "bootstrap": bootstrap}[sys.argv[1]]()
