# enrich_v07dev.py -- protocol v0.6 section 9 functional enrichment (development set only, descriptive) and the
# user-requested check (2026-10-10) "are strong-effect perturbations concentrated in protein complexes, and do
# they have more STRING neighbours?". Nothing here feeds any criterion or proxy.
#   python enrich_v07dev.py download   -> annotation files into <DL>/annotations_2026-10-10 (read-only) + manifest
#   python enrich_v07dev.py run        -> results/enrich_v07dev/
#
# Gene sets (read from saved results, nothing recomputed):
#   G1_core      genes miscovered in every audit-A test appearance under all three primary G1 models (93)
#   G1_extended  ... under at least one of the three primary G1 models (the protocol said "about 120"; the
#                exact count is reported)
#   D1_T5        D1 (v0.6.2) gene-level effect-strength tier 5 (median tier over the 3 rounds; M0, linear-PCA, k=50)
#   D1_always    D1 genes miscovered by M0 (linear-PCA, k=50, nominal) in all 3 rounds
#   controls     for each set, the same number of E genes with the highest mean E_true (D1)
# Background: the E genes (1,786). Statistics: one-sided Fisher exact test, terms with 5-500 background genes,
# BH FDR < 0.05 within each source; redundancy reduction: greedy by p, drop a term whose overlap genes have
# Jaccard > 0.5 with an already kept term.
import gzip
import hashlib
import io
import json
import os
import pickle
import shutil
import stat
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import RESULTS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DL = Path(os.environ.get("PERTAUDIT_DL", r"G:\biology"))
# v0.7 (user decision 2026-10-10): GO OBO and GAF from ONE dated GO release. 2026-05-19 is the most recent
# release that ships goa_human.gaf.gz in the legacy layout (later releases moved to annotations/gaf/).
# Reactome and CORUM are copied (hash-verified) from the 2026-10-10 download.
GO_RELEASE = "2026-05-19"
ANN = DL / f"annotations_go-{GO_RELEASE}"
ANN_PREV = DL / "annotations_2026-10-10"
OUT = RESULTS / "enrich_v07"
URLS = {
    "go-basic.obo": f"https://release.geneontology.org/{GO_RELEASE}/ontology/go-basic.obo",
    "goa_human.gaf.gz": f"https://release.geneontology.org/{GO_RELEASE}/annotations/goa_human.gaf.gz",
    "go_release-date.json": f"https://release.geneontology.org/{GO_RELEASE}/metadata/release-date.json",
    "ReactomePathways.gmt.zip": "https://reactome.org/download/current/ReactomePathways.gmt.zip",
    "corum_human_current.txt": ("https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/"
                                "download_current_file?file_id=human&file_format=txt"),
}
# CORUM's server omits its intermediate certificate; the issuer's published intermediate is added to the CA
# bundle so that TLS verification stays on.
HARICA_INTERMEDIATE = "http://crt.harica.gr/HARICA-GEANT-TLS-R1.cer"
MIN_SIZE, MAX_SIZE, FDR, JACC = 5, 500, 0.05, 0.5
STRING_THR = (0.4, 0.7)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def download():
    import ssl
    import certifi
    ANN.mkdir(parents=True, exist_ok=True)
    der = requests.get(HARICA_INTERMEDIATE, timeout=60).content
    pem = der.decode() if der.startswith(b"-----") else ssl.DER_cert_to_PEM_cert(der)
    bundle = ANN / "_ca_bundle_certifi_plus_harica_geant.pem"
    if not bundle.exists():
        bundle.write_text(Path(certifi.where()).read_text() + "\n" + pem)
    man = {"downloaded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "files": {}}
    prev = json.loads((ANN_PREV / "annotation_manifest.json").read_text())["files"]
    for name, url in URLS.items():
        f = ANN / name
        if not f.exists() and name in ("ReactomePathways.gmt.zip", "corum_human_current.txt"):
            shutil.copy2(ANN_PREV / name, f)                       # copy, then verify
            assert sha256(f) == prev[name]["sha256"], f"copied {name} differs from the 2026-10-10 download"
            os.chmod(f, stat.S_IREAD)
            man["files"][name] = {**prev[name], "copied_from": str(ANN_PREV / name)}
            continue
        if not f.exists():
            r = requests.get(url, timeout=600, verify=str(bundle) if "helmholtz" in url else True)
            r.raise_for_status()
            tmp = f.with_suffix(f.suffix + ".part")
            tmp.write_bytes(r.content)
            tmp.replace(f)
            os.chmod(f, stat.S_IREAD)
            man["files"][name] = {"url": url, "bytes": f.stat().st_size, "sha256": sha256(f),
                                  "last_modified": r.headers.get("Last-Modified")}
        else:
            man["files"][name] = {"url": url, "bytes": f.stat().st_size, "sha256": sha256(f), "last_modified": "existing"}
    man["go_release"] = (ANN / "go_release-date.json").read_text().strip()
    with gzip.open(ANN / "goa_human.gaf.gz", "rt", encoding="utf8") as fh:
        man["gaf_header"] = [ln.strip() for ln, _ in zip(fh, range(40)) if ln.startswith("!")][:12]
    with open(ANN / "go-basic.obo", encoding="utf8") as fh:
        man["obo_header"] = [ln.strip() for ln, _ in zip(fh, range(10))]
    (ANN / "annotation_manifest.json").write_text(json.dumps(man, indent=1), encoding="utf8")
    print(json.dumps({k: (v["bytes"], v["sha256"][:12]) for k, v in man["files"].items()}))


# ------------------------------------------------------------------ annotations

def load_go():
    parents, names, ns = {}, {}, {}
    cur = None
    with open(ANN / "go-basic.obo", encoding="utf8") as fh:
        for ln in fh:
            ln = ln.strip()
            if ln == "[Term]":
                cur = {"is_obsolete": False, "par": []}
            elif ln.startswith("[") and ln.endswith("]"):
                cur = None
            elif cur is not None and ln == "" and "id" in cur:
                if not cur["is_obsolete"]:
                    parents[cur["id"]] = cur["par"]
                    names[cur["id"]] = cur.get("name", "")
                    ns[cur["id"]] = cur.get("namespace", "")
                cur = None
            elif cur is not None:
                if ln.startswith("id: "):
                    cur["id"] = ln[4:]
                elif ln.startswith("name: "):
                    cur["name"] = ln[6:]
                elif ln.startswith("namespace: "):
                    cur["namespace"] = ln[11:]
                elif ln.startswith("is_a: "):
                    cur["par"].append(ln[6:].split(" ")[0])
                elif ln.startswith("relationship: part_of "):
                    cur["par"].append(ln.split(" ")[2])
                elif ln == "is_obsolete: true":
                    cur["is_obsolete"] = True
    anc = {}

    def ancestors(t):
        if t in anc:
            return anc[t]
        s = {t}
        for p in parents.get(t, []):
            s |= ancestors(p)
        anc[t] = s
        return s
    sys.setrecursionlimit(10000)
    gene_terms, syn = {}, {}
    with gzip.open(ANN / "goa_human.gaf.gz", "rt", encoding="utf8") as fh:
        for ln in fh:
            if ln.startswith("!"):
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 15 or "NOT" in c[3] or c[6] == "ND" or c[4] not in parents:
                continue
            g = c[2]
            gene_terms.setdefault(g, set()).update(ancestors(c[4]))
            for s in c[10].split("|"):
                if s:
                    syn.setdefault(s, set()).add(g)
    return gene_terms, syn, names, ns


def load_reactome():
    z = zipfile.ZipFile(ANN / "ReactomePathways.gmt.zip")
    txt = z.read(z.namelist()[0]).decode("utf8")
    sets = {}
    for ln in txt.splitlines():
        c = ln.split("\t")
        if len(c) > 2:
            sets[f"{c[1]} {c[0]}"] = set(c[2:])
    return sets


def load_corum():
    d = pd.read_csv(ANN / "corum_human_current.txt", sep="\t", dtype=str)
    d = d[d.organism.str.lower() == "human"]
    sets = {}
    for _, r in d.iterrows():
        genes = {g.strip() for g in str(r.subunits_gene_name).split(";") if g.strip() and g.strip() != "None"}
        if genes:
            sets[f"CORUM:{r.complex_id} {r.complex_name}"] = genes
    return sets


def map_symbols(E_genes, known, syn):
    """Screen target symbol -> annotation symbol: exact, else a unique GAF synonym."""
    m, how = {}, {"exact": 0, "synonym": 0, "unmapped": 0}
    for g in E_genes:
        if g in known:
            m[g] = g
            how["exact"] += 1
        elif g in syn and len(syn[g]) == 1:
            m[g] = next(iter(syn[g]))
            how["synonym"] += 1
        else:
            how["unmapped"] += 1
    return m, how


# ------------------------------------------------------------------ statistics

def enrich(query, background, term_genes, source, names=None):
    q, bg = set(query) & set(background), set(background)
    rows = []
    for t, gs in term_genes.items():
        gb = gs & bg
        if not (MIN_SIZE <= len(gb) <= MAX_SIZE):
            continue
        a = len(gb & q)
        if a == 0:
            continue
        b, c_ = len(q) - a, len(gb) - a
        d = len(bg) - a - b - c_
        orr, p = fisher_exact([[a, b], [c_, d]], alternative="greater")
        rows.append({"source": source, "term": t, "name": (names or {}).get(t, t), "n_bg": len(gb), "n_query": len(q),
                     "overlap": a, "expected": len(q) * len(gb) / len(bg), "odds_ratio": orr, "p": p,
                     "genes": ";".join(sorted(gb & q))})
    df = pd.DataFrame(rows)
    if len(df):
        df["fdr"] = multipletests(df.p, method="fdr_bh")[1]
        df = df.sort_values("p").reset_index(drop=True)
        keep, kept = [], []
        for i, r in df[df.fdr < FDR].iterrows():
            gs = set(r.genes.split(";"))
            if all(len(gs & k) / len(gs | k) <= JACC for k in kept):
                keep.append(i)
                kept.append(gs)
        df["representative"] = df.index.isin(keep)
    return df


def gene_sets(U):
    import analysis_v06 as an  # noqa: F401
    gi = pd.read_parquet(RESULTS / "g1_v0.5" / "analysis" / "glmm_input.parquet")
    a = gi[gi.audit == "A"]
    cnt = a.groupby(["model", "target"]).miscover.agg(["mean", "size"]).reset_index()
    per = {m: set(cnt[(cnt.model == m) & (cnt["mean"] == 1)].target) for m in ("mean", "linear_pca", "linear_go")}
    core = set.intersection(*per.values())
    ext = set.union(*per.values())
    cov = pd.read_parquet(RESULTS / "d1_v0.6.1" / "coverage_rows.parquet",
                          columns=["model", "k", "variant", "mode", "method", "perm", "target", "tier", "covered", "E_true"])
    d = cov[(cov.model == "linear_pca") & (cov.k == 50) & (cov.variant == "full") & (cov["mode"] == "nominal")
            & (cov.method == "M0") & (cov.perm == -1)]
    g = d.groupby("target").agg(tier=("tier", "median"), cov=("covered", "mean"), E=("E_true", "mean"), n=("covered", "size"))
    t5 = set(g.index[np.floor(g.tier) == 5])
    always = set(g.index[g["cov"] == 0])
    sets = {"G1_core": core, "G1_extended": ext, "D1_T5": t5, "D1_always": always}
    ranked = g.sort_values("E", ascending=False).index.tolist()
    for k in list(sets):
        sets[f"{k}__control_topE"] = set(ranked[:len(sets[k])])
    info = {k: len(v) for k, v in sets.items()}
    info["G1_per_model_always"] = {m: len(v) for m, v in per.items()}
    info["G1_min_appearances_of_core_genes"] = int(cnt[cnt.target.isin(core)]["size"].min())
    info["overlap_D1_T5_with_its_control"] = len(t5 & sets["D1_T5__control_topE"])
    return sets, info, g


def string_degree(E_genes):
    """Genome-wide STRING v12 degree from the detailed links: partners with experimental or database channel
    score >= threshold (the channels used for D_STRING). Returns {threshold: Series gene -> degree}."""
    base = DL / "string_v12.0"
    info = pd.read_csv(base / "9606.protein.info.v12.0.txt.gz", sep="\t", usecols=[0, 1])
    info.columns = ["protein", "gene"]
    p2g = dict(zip(info.protein, info.gene))
    out = {}
    cnt = {thr: {} for thr in STRING_THR}
    for ch in pd.read_csv(base / "9606.protein.links.detailed.v12.0.txt.gz", sep=" ", chunksize=2_000_000,
                          usecols=["protein1", "protein2", "experimental", "database"]):
        sc = np.maximum(ch.experimental.values, ch.database.values) / 1000.0
        for thr in STRING_THR:                       # links are listed in both directions
            vc = ch.protein1[sc >= thr].value_counts()
            for k, v in vc.items():
                cnt[thr][k] = cnt[thr].get(k, 0) + int(v)
    for thr in STRING_THR:
        deg = pd.Series(cnt[thr]).groupby(pd.Series(cnt[thr]).index.map(p2g)).sum()
        out[thr] = deg.reindex(list(E_genes)).fillna(0).astype(int)
    return out


def hypothesis(g, corum, gmap, pri_S, U, deg):
    """D1 T5 genes vs the other E genes."""
    E = g.index.tolist()
    t5 = (np.floor(g.tier) == 5).values
    in_c = np.array([any(gmap.get(x, x) in s for s in corum.values()) for x in E])
    n_c = np.array([sum(gmap.get(x, x) in s for s in corum.values()) for x in E])
    res = {"n_E": len(E), "n_T5": int(t5.sum()), "n_unmapped_for_corum": int(sum(x not in gmap for x in E))}
    tab = [[int((in_c & t5).sum()), int((~in_c & t5).sum())], [int((in_c & ~t5).sum()), int((~in_c & ~t5).sum())]]
    orr, p = fisher_exact(tab, alternative="two-sided")
    res["corum_member"] = {"frac_T5": float(in_c[t5].mean()), "frac_other": float(in_c[~t5].mean()), "table": tab,
                           "odds_ratio": float(orr), "p_two_sided": float(p)}
    res["corum_n_complexes"] = {"median_T5": float(np.median(n_c[t5])), "median_other": float(np.median(n_c[~t5])),
                                "mean_T5": float(n_c[t5].mean()), "mean_other": float(n_c[~t5].mean()),
                                "p_mannwhitney": float(mannwhitneyu(n_c[t5], n_c[~t5]).pvalue)}
    # STRING neighbours within the screen's target set (the matrix used for D_STRING) and genome-wide
    idx = {t: i for i, t in enumerate(U.target.values)}
    rows = np.array([idx[x] for x in E])
    S = pri_S[np.ix_(rows, np.arange(pri_S.shape[1]))].astype(float)
    S[np.arange(len(rows)), rows] = 0
    for thr in STRING_THR:
        dg = (S >= thr).sum(1)
        res[f"string_degree_screen_ge{thr}"] = {"median_T5": float(np.median(dg[t5])), "median_other": float(np.median(dg[~t5])),
                                                "mean_T5": float(dg[t5].mean()), "mean_other": float(dg[~t5].mean()),
                                                "p_mannwhitney": float(mannwhitneyu(dg[t5], dg[~t5]).pvalue)}
        dgw = deg[thr].reindex(E).values
        res[f"string_degree_genome_ge{thr}"] = {"median_T5": float(np.median(dgw[t5])), "median_other": float(np.median(dgw[~t5])),
                                                "mean_T5": float(dgw[t5].mean()), "mean_other": float(dgw[~t5].mean()),
                                                "p_mannwhitney": float(mannwhitneyu(dgw[t5], dgw[~t5]).pvalue)}
    # adjusted view (descriptive): T5 ~ CORUM member + log1p(screen degree >= 0.4) + log1p(number of GO annotations)
    import statsmodels.api as sm
    ngo = np.array([U.n_go.values[idx[x]] if "n_go" in U else np.nan for x in E])
    X = pd.DataFrame({"corum": in_c.astype(float), "log_deg": np.log1p((S >= 0.4).sum(1)), "log_ngo": np.log1p(ngo)})
    X = (X - X.mean()) / X.std()
    X.insert(0, "const", 1.0)
    try:
        r = sm.Logit(t5.astype(float), X.values).fit(disp=0)
        res["logit_T5"] = {c: {"coef_per_SD": float(b), "lo": float(b - 1.96 * s), "hi": float(b + 1.96 * s), "p": float(pv)}
                           for c, b, s, pv in zip(X.columns, r.params, r.bse, r.pvalues)}
    except Exception as e:                                                  # noqa: BLE001
        res["logit_T5"] = {"error": str(e)}
    return res


def run():
    import pipeline_v05 as pl
    OUT.mkdir(parents=True, exist_ok=True)
    U, Y, y_ctrl, genes, pri = pl.load()
    with open(pl.raw("gene2go_all.pkl"), "rb") as fh:
        g2g = pickle.load(fh)
    U = U.assign(n_go=[len(g2g.get(t, ())) for t in U.target.values])
    sets, info, g = gene_sets(U)
    E = g.index.tolist()
    assert len(E) == 1786, len(E)
    gene_terms, syn, names, ns = load_go()
    reactome, corum = load_reactome(), load_corum()
    known = set(gene_terms) | set().union(*reactome.values()) | set().union(*corum.values())
    gmap, how = map_symbols(E, known, syn)
    bg = [gmap[x] for x in E if x in gmap]
    go_sets = {}
    for gg, ts in gene_terms.items():
        for t in ts:
            go_sets.setdefault(t, set()).add(gg)
    go_by_ns = {n_: {t: s for t, s in go_sets.items() if ns.get(t) == n_} for n_ in
                ("biological_process", "cellular_component", "molecular_function")}
    tables = []
    for sname, s in sets.items():
        q = [gmap[x] for x in s if x in gmap]
        for src, tg, nm in [("GO:BP", go_by_ns["biological_process"], names), ("GO:CC", go_by_ns["cellular_component"], names),
                            ("GO:MF", go_by_ns["molecular_function"], names), ("Reactome", reactome, None), ("CORUM", corum, None)]:
            df = enrich(q, bg, tg, src, nm)
            if len(df):
                tables.append(df.assign(gene_set=sname))
    T = pd.concat(tables, ignore_index=True)
    T.to_csv(OUT / "enrichment_all_terms.csv", index=False)
    rep = T[T.representative].copy()
    rep.drop(columns=["genes"]).to_csv(OUT / "enrichment_representative.csv", index=False)
    summary = {sname: {src: int(((T.gene_set == sname) & (T.source == src) & (T.fdr < FDR)).sum())
                       for src in ("GO:BP", "GO:CC", "GO:MF", "Reactome", "CORUM")} for sname in sets}
    deg = string_degree(E)
    hyp = hypothesis(g, corum, gmap, pri["S"], U, deg)
    out = {"gene_set_sizes": info, "symbol_mapping": how, "n_background_mapped": len(bg),
           "n_significant_terms_fdr05": summary, "hypothesis_T5_vs_other": hyp,
           "annotation_manifest": json.loads((ANN / "annotation_manifest.json").read_text())["files"],
           "parameters": {"min_size": MIN_SIZE, "max_size": MAX_SIZE, "fdr": FDR, "jaccard_redundancy": JACC,
                          "go_propagation": "is_a + part_of; NOT qualifiers and ND evidence excluded",
                          "string_degree": "partners with max(experimental, database) >= threshold"}}
    (OUT / "enrichment_summary.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf8")
    (OUT / "gene_sets.json").write_text(json.dumps({k: sorted(v) for k, v in sets.items()}, indent=1), encoding="utf8")
    print(json.dumps({"sizes": info, "mapping": how, "n_sig": summary}, default=float))


if __name__ == "__main__":
    {"download": download, "run": run}[sys.argv[1]]()
