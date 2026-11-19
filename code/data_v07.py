# data_v07.py -- dataset-generic unit table, priors and loader (protocol v0.7; confirmation sets RPE1 and GWPS).
# Same rules as p0_build.py (unit table), build_priors_v04.string_universe (STRING universe) and
# pipeline_v05.load + pipeline_v06.load_dev_v06 (loader); only the input paths are parameters.
# Validated on K562 essential: rebuilt units and loaded objects are identical to the frozen dev versions
# (tests_v06.t_v07_data).
#
# Response isolation (v0.7 section 3, C0): build_units reads ONLY obs gene_transcript and num_cells_filtered and
# the readout gene names of the publisher's bulk file. The response-derived obs columns (energy_test_p_value,
# fold_expr, ...) are never read.
import json
import pickle
import time
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests

from audit_lib import tss_groups
from models_v04 import go_incidence, jaccard_all
from paths import PROCESSED, raw

MIN_CELLS = 30
TSS_MERGE_BP = 1000
GO_TERM_MAX = 500
COORDS_ARCHIVE = PROCESSED / "v0.3_figshare_k562ess" / "p0v3_units_all_targets.csv"   # Ensembl 116, GRCh38.p14
ENSEMBL = "https://rest.ensembl.org"


def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def read_bulk_meta(path):
    with h5py.File(path, "r") as fb:
        ob = pd.DataFrame({"gene_transcript": _dec(fb["obs/gene_transcript"][:]),
                           "ncells": fb["obs/num_cells_filtered"][:]})
        vcat = _dec(fb["var/__categories/gene_name"][:])
        readout = vcat[fb["var/gene_name"][:]]
    return ob, readout


def gears_universe():
    with open(raw("gene2go_all.pkl"), "rb") as fh:
        gene2go = pickle.load(fh)
    with open(raw("essential_all_data_pert_genes.pkl"), "rb") as fh:
        ess = pickle.load(fh)
    return gene2go, [g for g in ess if g in gene2go]


def ensembl_coords(ensgs, symbols, cache_csv, log=print):
    """GRCh38 coordinates of ENSG ids missing from the archive. Cached in cache_csv (with the Ensembl release),
    so a rerun never re-queries. Fallback: lookup by symbol; still missing -> chr 'NA_<target>' (P0 v0.3 rule)."""
    cache_csv = Path(cache_csv)
    if cache_csv.exists():
        return pd.read_csv(cache_csv)
    hdr = {"Content-Type": "application/json", "Accept": "application/json"}
    got = {}
    ids = sorted(set(ensgs))

    def post(batch):
        """POST lookup; on repeated server errors split the batch (a single failing id -> treated as not found,
        then looked up by symbol below). Post-freeze code fix 2026-10-10 (HTTP 500 on one batch)."""
        for _ in range(4):
            try:
                r = requests.post(f"{ENSEMBL}/lookup/id", headers=hdr, json={"ids": batch}, timeout=180)
                if r.ok:
                    got.update({k: v for k, v in r.json().items() if v})
                    return
            except requests.RequestException:
                pass
            time.sleep(5)
        if len(batch) > 1:
            h = len(batch) // 2
            post(batch[:h])
            post(batch[h:])
    for s in range(0, len(ids), 100):
        post(ids[s:s + 100])
    sym = dict(zip(ensgs, symbols))
    how = {}
    for e in ids:
        if e in got:
            how[e] = "id"
            continue
        r = requests.get(f"{ENSEMBL}/lookup/symbol/homo_sapiens/{sym[e]}", headers={"Accept": "application/json"},
                         timeout=60)
        if r.ok:
            got[e], how[e] = r.json(), "symbol"
        else:
            how[e] = "missing"
    rel = requests.get(f"{ENSEMBL}/info/data/", headers={"Accept": "application/json"}, timeout=60).json()
    asm = requests.get(f"{ENSEMBL}/info/assembly/homo_sapiens", headers={"Accept": "application/json"},
                       timeout=60).json().get("assembly_name")
    rows = []
    for e in ids:
        v = got.get(e)
        if v:
            rows.append({"ensg": e, "chr": str(v["seq_region_name"]), "strand": {1: "+", -1: "-"}[v["strand"]],
                         "start": int(v["start"]), "end": int(v["end"]), "how": how[e]})
        else:
            rows.append({"ensg": e, "chr": f"NA_{sym[e]}", "strand": "+", "start": 0, "end": 0, "how": "missing"})
    df = pd.DataFrame(rows)
    df["ensembl_release"] = json.dumps(rel)
    df["assembly"] = asm
    cache_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache_csv, index=False)
    log(f"Ensembl: {len(ids)} ids, by id {sum(v == 'id' for v in how.values())}, by symbol "
        f"{sum(v == 'symbol' for v in how.values())}, missing {sum(v == 'missing' for v in how.values())}, "
        f"release {rel}, assembly {asm}")
    return df


def build_units(bulk_path, coords_cache_csv, exclude_targets=(), log=print):
    """Unit table of one screen (rules of p0_build.py). exclude_targets: targets removed before anything else
    (GWPS: the K562 essential targets never enter any fold, including training)."""
    ob, readout = read_bulk_meta(bulk_path)
    gene2go, universe = gears_universe()
    pert_names = set(universe)
    parts = ob.gene_transcript.str.split("_", n=3, expand=True)
    ob["gid"], ob["target"], ob["tss"], ob["ensg"] = parts[0], parts[1], parts[2], parts[3]
    ob["is_ctrl"] = ob.gene_transcript.str.contains("non-targeting", case=False)
    P = ob[~ob.is_ctrl].copy()
    n_excluded_rows = int(P.target.isin(set(exclude_targets)).sum())
    P = P[~P.target.isin(set(exclude_targets))]
    U = P.sort_values("ncells", ascending=False, kind="quicksort").drop_duplicates("target").copy()
    Rset = set(readout)
    U["in_R"] = U.target.isin(Rset)
    U["in_GO"] = U.target.isin(pert_names)
    U["cells_ok"] = U.ncells >= MIN_CELLS
    U["E"] = U.in_R & U.in_GO & U.cells_ok
    U["E_noGO"] = U.in_R & ~U.in_GO & U.cells_ok
    arch = pd.read_csv(COORDS_ARCHIVE, usecols=["ensg", "chr", "strand", "start", "end"]).drop_duplicates("ensg")
    arch["how"] = "archive_ensembl116"
    miss = sorted(set(U.ensg) - set(arch.ensg))
    if miss:
        sym = U.set_index("ensg").target.to_dict()
        extra = ensembl_coords(miss, [sym[e] for e in miss], coords_cache_csv, log)
        c = pd.concat([arch, extra[["ensg", "chr", "strand", "start", "end", "how"]]], ignore_index=True)
    else:
        c = arch
    c = c.drop_duplicates("ensg").set_index("ensg")
    for col in ("chr", "strand", "start", "end"):
        U[col] = U.ensg.map(c[col])
    U["coord_source"] = U.ensg.map(c["how"])
    assert U.chr.notna().all(), "coordinates missing for some targets"
    U = U.reset_index(drop=True)
    U["tss_group"] = tss_groups(U.target.values, U.chr.values.astype(str), U.strand.values,
                                U.start.values.astype(np.int64), U.end.values.astype(np.int64), TSS_MERGE_BP)
    grp = U.groupby("tss_group").agg(nE=("E", "sum"), n=("E", "size"))
    mixed = grp[(grp.nE > 0) & (grp.n > grp.nE)].index
    U["always_train"] = ~U.E & ~U.tss_group.isin(mixed)
    U["dropped_tss_conflict"] = ~U.E & U.tss_group.isin(mixed)
    assert U.target.is_unique
    info = {"rows_bulk": int(len(ob)), "control_rows": int(ob.is_ctrl.sum()), "excluded_rows": n_excluded_rows,
            "secondary_rows": int(len(P) - len(U)), "targets": int(len(U)), "readout_genes": int(len(readout)),
            "in_R": int(U.in_R.sum()), "in_GO": int(U.in_GO.sum()), "cells_ge30": int(U.cells_ok.sum()),
            "E": int(U.E.sum()), "E_noGO": int(U.E_noGO.sum()), "E_tss_groups": int(U[U.E].tss_group.nunique()),
            "always_train": int(U.always_train.sum()), "dropped_tss_conflict": int(U.dropped_tss_conflict.sum()),
            "coord_source": U.coord_source.value_counts().to_dict(),
            "ncells_E": {"min": int(U[U.E].ncells.min()), "median": float(U[U.E].ncells.median()),
                         "max": int(U[U.E].ncells.max())} if U.E.any() else None}
    return U, readout, info


def string_universe(targets, out_parquet):
    """STRING v12 experimental + database channels within a target universe (build_priors_v04 rule)."""
    info = pd.read_csv(raw("9606.protein.info.v12.0.txt.gz"), sep="\t")
    info.columns = [c.lstrip("#") for c in info.columns]
    p2g = dict(zip(info.string_protein_id, info.preferred_name))
    univ = set(targets)
    keep = []
    for ch in pd.read_csv(raw("9606.protein.links.detailed.v12.0.txt.gz"), sep=" ", chunksize=2_000_000,
                          usecols=["protein1", "protein2", "experimental", "database"]):
        g1 = ch.protein1.map(p2g)
        g2 = ch.protein2.map(p2g)
        m = g1.isin(univ) & g2.isin(univ) & ((ch.experimental > 0) | (ch.database > 0))
        keep.append(pd.DataFrame({"g1": g1[m], "g2": g2[m], "exp": ch.experimental[m] / 1000, "db": ch.database[m] / 1000}))
    st = pd.concat(keep, ignore_index=True)
    prior = 0.041

    def nopr(s):
        return np.clip((s - prior) / (1 - prior), 0, None)
    comb = 1 - (1 - nopr(st.exp)) * (1 - nopr(st.db))
    st["score"] = comb + prior * (1 - comb)
    st = st.groupby(["g1", "g2"], as_index=False).score.max()
    Path(out_parquet).parent.mkdir(parents=True, exist_ok=True)
    st.to_parquet(out_parquet)
    return len(st)


def load_dataset(U_all, pb_dir, string_parquet):
    """pipeline_v05.load + pipeline_v06.load_dev_v06 with parameterised inputs. Returns the D dict of run_dataset."""
    pb_dir = Path(pb_dir)
    pb = np.load(pb_dir / "pseudobulk_logmean.npy").astype(np.float64)
    rows = pd.read_csv(pb_dir / "pseudobulk_rows.csv")
    genes = pd.read_csv(pb_dir / "pseudobulk_genes.csv").gene_name.values
    rix = {g: i for i, g in enumerate(rows.gene_transcript)}
    U = U_all.copy()
    assert U.gene_transcript.isin(rix).all(), "unit rows missing from pseudobulk"
    ctrl_rows = np.where(rows.gene_transcript.str.contains("non-targeting", case=False))[0]
    w = rows.n_cells.values[ctrl_rows].astype(np.float64)
    assert len(ctrl_rows) > 0 and np.isfinite(w).all() and (w > 0).all(), "control cell counts invalid"
    y_ctrl = (pb[ctrl_rows] * w[:, None]).sum(0) / w.sum()
    U = U[U.E | U.always_train].reset_index(drop=True)
    U["row"] = U.gene_transcript.map(rix)
    U["n_cells_sc"] = rows.n_cells.values[U.row.values]
    gcols = {}
    for j, g in enumerate(genes):
        gcols.setdefault(g, []).append(j)
    U["tcols"] = [np.array(gcols.get(t, []), dtype=int) for t in U.target]
    Y = pb[U.row.values]
    del pb
    gene2go, universe = gears_universe()
    cnt = Counter(t for g in set(universe) for t in gene2go[g])
    broad = frozenset(t for t, c in cnt.items() if c > GO_TERM_MAX)
    ase = np.load(PROCESSED / "priors" / "go_ase100.npz", allow_pickle=True)
    ase_ix = {g: i for i, g in enumerate(ase["genes"])}
    st = pd.read_parquet(string_parquet)
    pos = {g: i for i, g in enumerate(U.target)}
    S = np.zeros((len(U), len(U)))
    st = st[st.g1.isin(pos) & st.g2.isin(pos)]
    for a, b, s in st.itertuples(index=False):
        i, j = pos[a], pos[b]
        S[i, j] = S[j, i] = max(S[i, j], s)
    np.fill_diagonal(S, 0.0)
    pri = {"J": jaccard_all(go_incidence(U.target.values, gene2go)),
           "Jf": jaccard_all(go_incidence(U.target.values, gene2go, broad)),
           "S": S, "ase": ase["ase"], "ase_ix": ase_ix,
           "n_universe": len(set(universe)), "n_broad": len(broad)}
    r = U.row.values
    YA = np.load(pb_dir / "pseudobulk_logmean_A.npy").astype(np.float64)[r]
    YB = np.load(pb_dir / "pseudobulk_logmean_B.npy").astype(np.float64)[r]
    pri["n_go"] = np.array([len(gene2go.get(t, ())) for t in U.target.values], float)
    return {"U": U, "Y": Y, "YA": YA, "YB": YB, "y_ctrl": y_ctrl, "pri": pri, "genes": genes,
            "E_idx": np.where(U.E.values)[0], "AT_idx": np.where(U.always_train.values)[0]}
