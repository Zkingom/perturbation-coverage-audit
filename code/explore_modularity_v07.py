# explore_modularity_v07.py -- POST-CONFIRMATION EXPLORATORY analysis (protocol reports/protocol/explore_modularity_v1.0.md).
# Does network modularity of effect strength, and the split-half reliability of E_true, explain the predictability of
# the prior-neighbour proxy h on essential (dev, RPE1) vs non-essential (GWPS) genes? Reads saved caches / rows only;
# no model is refitted, no confirmatory result is touched.
#   python explore_modularity_v07.py [--boot 2000] [--datasets dev rpe1 gwps]
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED, RESULTS  # noqa: E402

OUT = RESULTS / "explore_modularity_v07"
HCFG, MODEL, K = "sim_prior|ridge", "linear_pca", 50
STRING_THR = 0.7
SEED_OFF = {"dev": 0, "rpe1": 10, "gwps": 20}


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


# ------------------------------------------------------------------ inputs

def dataset_cfg(ds):
    if ds == "dev":
        res = RESULTS / "d1_v0.7dev"
        return {"res": res, "cache": RESULTS / "d1_v0.6.1" / "cache",
                "ref_json": res / "confirm_v07_on_dev" / "confirm_v07_dev.json",
                "rand": res / "sens_rand_match_rows_linear_pca.parquet",
                "string": PROCESSED / "priors" / "string_expdb_v0.4_universe.parquet"}
    res = RESULTS / "conf_v0.7" / ds
    return {"res": res, "cache": res / "cache", "ref_json": res / "confirm_v07_result.json",
            "rand": res / "confirm_rand_rows.parquet",
            "string": PROCESSED / "v0.7" / ds / "string_expdb_universe.parquet"}


def load_D(ds):
    import lowmem_v07 as lm
    if ds == "dev":
        U = pd.read_csv(PROCESSED / "v0.3_figshare_k562ess" / "p0v3_units_all_targets.csv")
        return lm.load_dataset_lowmem(U, PROCESSED / "v0.6" / "k562ess",
                                      PROCESSED / "priors" / "string_expdb_v0.4_universe.parquet")
    import confirm_run_v07 as cr
    return cr.load_lowmem(ds)


def gene_table(ds, cfg):
    """Gene-level h, E_true, tier (analysis_v06.gene_proxy) + split-half E_A, E_B; AT-unit E for rho_net."""
    import analysis_v06 as an
    import proxy_v06 as px
    pt = pd.read_parquet(cfg["res"] / "proxy_rows.parquet",
                         filters=[("kind", "==", "h"), ("model", "==", MODEL), ("proxy", "==", HCFG)])
    pt = pt.rename(columns={"value": "h"})
    ref_rho = json.loads(Path(cfg["ref_json"]).read_text(encoding="utf8"))["P"]["rho"]
    chk = an.predictability(pt, boot=False)["rho"]
    assert chk == ref_rho, f"{ds}: rho {chk} != confirmatory P.rho {ref_rho}"
    D = load_D(ds)
    Ud = D["U"]
    pos = {t: i for i, t in enumerate(Ud.target.values)}
    rows, at_rows = [], []
    files = sorted(Path(cfg["cache"]).glob("split_*.pkl"))
    assert len(files) == 21, len(files)
    for f in files:
        c = pickle.load(open(f, "rb"))
        idx = np.array([pos[t] for t in c["target_calt"]])
        ref = {"mu": c["ref_mu"], "sigma": c["ref_sigma"]}
        tc = [Ud.tcols.values[i] for i in idx]
        EA = px.e_true(D["YA"][idx], ref, tc)
        assert np.array_equal(EA, c["EA_calt"]), f"{f.name}: E_A recomputation differs"
        if f.name == "split_00.pkl":
            assert np.array_equal(px.e_true(D["Y"][idx], ref, tc), c["E_calt"]), "E_true recomputation differs"
        EB = px.e_true(D["YB"][idx], ref, tc)
        nt = len(c["test"])
        sl = slice(len(idx) - nt, len(idx))
        rows.append(pd.DataFrame({"split": c["split_id"], "target": c["target_calt"][sl], "E_calt": c["E_calt"][sl],
                                  "EA": EA[sl], "EB": EB[sl]}))
        at = np.asarray(c["at"])
        assert Ud.always_train.values[at].all(), "AT index mismatch"
        at_rows.append(pd.DataFrame({"target": Ud.target.values[at], "E": c["E_at"]}))
    del D
    hr = pd.concat(rows, ignore_index=True)
    m = pt[["split", "target", "E_true"]].merge(hr, on=["split", "target"], how="left", validate="one_to_one")
    assert m.EA.notna().all() and np.array_equal(m.E_true.values, m.E_calt.values), f"{ds}: test-row join failed"
    g = an.gene_proxy(pt)
    ab = hr.groupby("target").agg(EA=("EA", "mean"), EB=("EB", "mean"))
    g = g.join(ab, on="target")
    assert g.EA.notna().all()
    at = pd.concat(at_rows, ignore_index=True).groupby("target").E.mean()
    return g, at, pt


# ------------------------------------------------------------------ annotations

def annotations(genes):
    import enrich_v07dev as en
    gene_terms, syn, _, _ = en.load_go()
    reactome, corum = en.load_reactome(), en.load_corum()
    known = set(gene_terms) | set().union(*reactome.values()) | set().union(*corum.values())
    gmap, how = en.map_symbols(sorted(genes), known, syn)
    members = set().union(*corum.values())
    corum_state = {x: ("unmapped" if x not in gmap else ("member" if gmap[x] in members else "non")) for x in genes}
    info = pd.read_csv(en.DL / "string_v12.0" / "9606.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 1])
    in_string = set(info.iloc[:, 1])
    deg = en.string_degree(sorted(genes))[STRING_THR]
    deg = {x: (int(deg[x]) if x in in_string else None) for x in genes}
    return corum_state, deg, {"corum_symbol_mapping": how, "n_not_in_string": int(sum(v is None for v in deg.values()))}


def tertile(deg_series, cuts):
    q1, q2 = cuts
    return np.where(deg_series <= q1, "low", np.where(deg_series <= q2, "mid", "high"))


# ------------------------------------------------------------------ statistics

def stats(x):
    n = len(x)
    out = {"n": int(n)}
    if n < 10:
        return out
    logE = np.log(x.E.values)
    out.update({"frac_T5": float(x.T5.mean()), "E_sd": float(np.std(x.E, ddof=1)),
                "E_iqr": float(np.subtract(*np.percentile(x.E, [75, 25]))), "logE_sd": float(np.std(logE, ddof=1)),
                "logE_iqr": float(np.subtract(*np.percentile(logE, [75, 25])))})
    rho = spearmanr(x.h, x.E).correlation
    r = spearmanr(x.EA, x.EB).correlation
    rl = pearsonr(np.log(x.EA), np.log(x.EB))[0]
    rel = 2 * r / (1 + r)
    out.update({"rho": float(rho), "auroc_T5": float(roc_auc_score(x.T5, x.h)) if x.T5.nunique() == 2 else np.nan,
                "r_AB": float(r), "rel_SB": float(rel), "r_AB_logpearson": float(rl),
                "rel_SB_logpearson": float(2 * rl / (1 + rl)),
                "rho_corr": float(rho / np.sqrt(rel)) if rel > 0 else np.nan})
    return out


BOOT_KEYS = ("rho", "auroc_T5", "r_AB", "rel_SB", "rho_corr", "r_AB_logpearson", "rel_SB_logpearson")


def boot_groups(g, groups, n_boot, seed):
    """groups: name -> boolean mask over g rows. One TSS-group resample per draw (cluster_boot procedure), statistics
    recomputed inside every group on the same resample. Returns {group: {key: array}}."""
    rng = np.random.default_rng(seed)
    gs = g.tss_group.unique()
    posd = g.groupby("tss_group", sort=False).indices
    res = {k: {s: np.full(n_boot, np.nan) for s in BOOT_KEYS} for k in groups}
    for b in range(n_boot):
        pick = rng.choice(gs, len(gs), replace=True)
        ix = np.concatenate([posd[x] for x in pick])
        gb = g.take(ix).reset_index(drop=True)
        for k, msk in groups.items():
            st = stats(gb[msk[ix]])
            for s in BOOT_KEYS:
                if s in st:
                    res[k][s][b] = st[s]
    return res


def ci(a):
    a = a[np.isfinite(a)]
    return [float(np.quantile(a, 0.025)), float(np.quantile(a, 0.975))] if len(a) > 10 else [np.nan, np.nan]


def boot_p(d):
    d = d[np.isfinite(d)]
    return float(min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))) if len(d) else np.nan


def rho_net(g, at_E, string_path, seed, n_boot):
    st = pd.read_parquet(string_path)
    st = st[st.score >= STRING_THR]
    logE = pd.concat([np.log(g.set_index("target").E), np.log(at_E)])
    logE = logE[~logE.index.duplicated()]
    nb = {}
    for a, b, s in st[["g1", "g2", "score"]].itertuples(index=False):
        if a in logE.index and b in logE.index and a != b:
            nb.setdefault(a, {})[b] = max(nb.get(a, {}).get(b, 0), s)
            nb.setdefault(b, {})[a] = max(nb.get(b, {}).get(a, 0), s)
    own, nbm, tss = [], [], []
    for t, tg in zip(g.target, g.tss_group):
        d = nb.get(t)
        if d:
            w = np.array(list(d.values()))
            v = logE.loc[list(d.keys())].values
            own.append(logE[t])
            nbm.append(float((w * v).sum() / w.sum()))
            tss.append(tg)
    x = pd.DataFrame({"own": own, "nb": nbm, "tss_group": tss})
    point = float(spearmanr(x.own, x.nb).correlation)
    rng = np.random.default_rng(seed)
    gs = x.tss_group.unique()
    posd = x.groupby("tss_group", sort=False).indices
    vals = np.empty(n_boot)
    for b in range(n_boot):
        ix = np.concatenate([posd[k] for k in rng.choice(gs, len(gs), replace=True)])
        xb = x.take(ix)
        vals[b] = spearmanr(xb.own, xb.nb).correlation
    return {"rho_net": point, "ci": ci(vals), "n_genes_with_neighbours": int(len(x)), "n_genes": int(len(g))}, vals


def split_half_oracle(res):
    import analysis_v06 as an
    flt = [("model", "==", MODEL), ("k", "==", K), ("mode", "==", "nominal"), ("perm", "==", -1),
           ("method", "in", ["M0", "oracle"]), ("variant", "in", ["full", "half"])]
    d = pd.read_parquet(res / "coverage_rows.parquet", filters=flt)
    out = {}
    for var in ("full", "half"):
        for meth, hy in (("M0", None), ("oracle", "0.0")):
            x = d[(d.variant == var) & (d.method == meth)]
            if hy is not None:
                x = x[x.hyper == hy]
            x = x.sort_values(["split", "target"])
            gt = an.gee_tiers(x)
            out[f"{var}|{meth}"] = {"n_rows": int(len(x)), "T5": gt[5], "marginal": float(x.covered.mean())}
    return out


def group_E1(rand_path, targets_by_group):
    import analysis_v06 as an
    d = pd.read_parquet(rand_path, filters=[("model", "==", MODEL), ("k", "==", K), ("method", "in", ["M0", "M1"])])
    if "variant" in d:
        d = d[d.variant == "full"]
    if "perm" in d:
        d = d[d.perm == -1]
    d0 = d[d.method == "M0"]
    d1 = d[(d.method == "M1") & (d.hyper.astype(str) == "0.0")]
    out = {}
    for k, tg in targets_by_group.items():
        a = d1[d1.target.isin(tg)].sort_values(["split", "target"])
        b = d0[d0.target.isin(tg)].sort_values(["split", "target"])
        out[k] = an.paired_T5(a, b) if (b.tier == 5).sum() >= 20 else {"n_T5_rows": int((b.tier == 5).sum())}
    return out


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--datasets", nargs="*", default=["dev", "rpe1", "gwps"])
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    G, AT, PT = {}, {}, {}
    for ds in a.datasets:
        log(f"{ds}: gene table")
        G[ds], AT[ds], PT[ds] = gene_table(ds, dataset_cfg(ds))
    allg = set().union(*[set(G[ds].target) for ds in G])
    log(f"annotations for {len(allg)} genes")
    corum, deg, ann_info = annotations(allg)
    res = {"protocol": "reports/protocol/explore_modularity_v1.0.md", "label": "POST-CONFIRMATION EXPLORATORY",
           "annotation_info": ann_info, "n_boot": a.boot, "datasets": {}}
    cuts = {}
    for ds in G:
        g = G[ds]
        g["corum"] = g.target.map(corum)
        g["deg"] = g.target.map(deg)
        dm = g.deg.notna()
        cuts[ds] = tuple(float(v) for v in np.quantile(g.deg[dm].astype(float), [1 / 3, 2 / 3]))
    fixed = cuts.get("rpe1")
    BOOTS = {}
    rows = []
    for ds in G:
        g = G[ds].reset_index(drop=True)
        dm = g.deg.notna().values
        dd = g.deg.fillna(-1).astype(float).values
        t_own = np.where(dm, tertile(dd, cuts[ds]), "unmapped")
        groups = {"all": np.ones(len(g), bool), "corum_member": (g.corum == "member").values,
                  "corum_non": (g.corum == "non").values}
        for lv in ("low", "mid", "high"):
            groups[f"deg_{lv}"] = t_own == lv
        if fixed is not None:
            t_fix = np.where(dm, tertile(dd, fixed), "unmapped")
            for lv in ("low", "mid", "high"):
                groups[f"degfixRPE1_{lv}"] = t_fix == lv
        log(f"{ds}: bootstrap ({a.boot} draws, {len(groups)} groups)")
        bt = boot_groups(g, groups, a.boot, MASTER_SEED + 9100 + SEED_OFF[ds])
        BOOTS[ds] = bt
        gres = {}
        for k, msk in groups.items():
            st = stats(g[msk])
            for s in BOOT_KEYS:
                if s in st:
                    st[s + "_ci"] = ci(bt[k][s])
            gres[k] = st
            rows.append({"dataset": ds, "group": k, **{kk: v for kk, v in st.items() if not kk.endswith("_ci")},
                         **{kk: str(v) for kk, v in st.items() if kk.endswith("_ci")}})
        contr = {}
        for nm, (p, q) in {"corum_member_minus_non": ("corum_member", "corum_non"),
                           "deg_high_minus_low": ("deg_high", "deg_low")}.items():
            for s in ("rho", "auroc_T5"):
                dlt = bt[p][s] - bt[q][s]
                contr[f"{nm}|{s}"] = {"diff": gres[p][s] - gres[q][s], "ci": ci(dlt), "p_boot": boot_p(dlt)}
        log(f"{ds}: rho_net, split-half oracle, group E1")
        rn, rn_b = rho_net(g, AT[ds], dataset_cfg(ds)["string"], MASTER_SEED + 9100 + SEED_OFF[ds] + 1, a.boot)
        BOOTS[ds]["_rho_net"] = rn_b
        tg = {k: set(g.target[msk]) for k, msk in groups.items()}
        res["datasets"][ds] = {"n_genes": int(len(g)), "deg_tertile_cuts": cuts[ds],
                               "n_unmapped": {"corum": int((g.corum == "unmapped").sum()), "string": int((~dm).sum())},
                               "groups": gres, "contrasts": contr, "rho_net": rn,
                               "split_half_oracle": split_half_oracle(dataset_cfg(ds)["res"]),
                               "group_E1_rand_matched": group_E1(dataset_cfg(ds)["rand"], tg)}
    # Holm over the 8 primary contrasts (RPE1, GWPS)
    fam = [(ds, k) for ds in ("rpe1", "gwps") if ds in res["datasets"] for k in res["datasets"][ds]["contrasts"]]
    ps = np.array([res["datasets"][ds]["contrasts"][k]["p_boot"] for ds, k in fam])
    order = np.argsort(ps)
    adj = np.empty(len(ps))
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - r) * ps[i]))
        adj[i] = run
    for (ds, k), v in zip(fam, adj):
        res["datasets"][ds]["contrasts"][k]["p_holm"] = float(v)
    # cross-dataset (independent bootstraps): RPE1 - GWPS
    if "rpe1" in BOOTS and "gwps" in BOOTS:
        cx = {}
        dr = BOOTS["rpe1"]["_rho_net"] - BOOTS["gwps"]["_rho_net"]
        cx["rho_net"] = {"diff": res["datasets"]["rpe1"]["rho_net"]["rho_net"] - res["datasets"]["gwps"]["rho_net"]["rho_net"],
                         "ci": ci(dr), "p_boot": boot_p(dr)}
        for k in ["all", "degfixRPE1_low", "degfixRPE1_mid", "degfixRPE1_high"]:
            for s in ("rho", "rho_corr", "rel_SB"):
                dlt = BOOTS["rpe1"][k][s] - BOOTS["gwps"][k][s]
                cx[f"{k}|{s}"] = {"diff": res["datasets"]["rpe1"]["groups"][k].get(s, np.nan)
                                  - res["datasets"]["gwps"]["groups"][k].get(s, np.nan), "ci": ci(dlt)}
        res["cross_rpe1_minus_gwps"] = cx
        # pre-specified rule 8d (GWPS, all genes)
        ga = res["datasets"]["gwps"]["groups"]["all"]
        ruleN = ga["rel_SB_ci"][1] < 0.6 and ga["rho_corr"] >= 0.30
        ruleI = ga["rel_SB_ci"][0] >= 0.6 or ga["rho_corr_ci"][1] < 0.20
        res["rule_8d"] = {"rule_N_noise": bool(ruleN), "rule_I_information": bool(ruleI) and not ruleN,
                          "undetermined": bool(not ruleN and not ruleI)}
    (out_dir / "explore_results.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf8")
    pd.DataFrame(rows).to_csv(out_dir / "explore_group_table.csv", index=False)
    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
