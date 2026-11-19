# gears_bundle_v07.py -- build the upload bundle for the GEARS robustness check (protocol v0.6 section 10, inherited by
# v0.7): RPE1 essential, the 21 frozen confirmation splits, GEARS 0.1.2 defaults.
#   python gears_bundle_v07.py build      -> <ROOT>/gears_bundle_v07/ (data, priors, splits, code, SHA256SUMS)
# Data rule (protocol): at most 16 cells per perturbation (uniform without replacement, fixed seed); all non-targeting
# control cells. Expression transform identical to pseudobulk_v06: log1p(1e4 * counts / UMI_count) in float64, stored
# as float32 CSR. Units: E units (all in the GEARS GO graph) + always-train units in the GEARS GO graph (always-train
# units outside the graph cannot be represented by GEARS and are excluded; listed in bundle_info.json).
import hashlib
import json
import os
import pickle
import shutil
import sys
import tarfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MASTER_SEED, PROCESSED, ROOT, SPLITS, raw  # noqa: E402
from pseudobulk_v04 import read_obs_column, var_gene_names  # noqa: E402

CAP = 16
SEED = MASTER_SEED + 9700
OUT = ROOT / "gears_bundle_v07"
PB = PROCESSED / "v0.7" / "rpe1"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def var_column(f, name):
    """A var column of an old-style (anndata <0.8) h5ad: codes into var/__categories/<name>, or plain strings."""
    dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])   # noqa: E731
    v = f["var"][name][:]
    if "__categories" in f["var"] and name in f["var"]["__categories"]:
        return dec(f["var"]["__categories"][name][:])[v]
    return dec(v)


def build():
    import anndata as ad
    import confirm_run_v07 as cr
    t0 = time.time()
    for d in ("data", "priors", "splits", "code"):
        (OUT / d).mkdir(parents=True, exist_ok=True)
    U = pd.read_csv(PB / "units_all_targets.csv")
    ess = pickle.load(open(raw("essential_all_data_pert_genes.pkl"), "rb"))
    g2g = pickle.load(open(raw("gene2go_all.pkl"), "rb"))
    graph = set(g for g in ess if g in g2g)
    E = U[U.E]
    assert E.target.isin(graph).all(), "an E unit is outside the GEARS GO graph"
    AT = U[U.always_train]
    AT_in, AT_out = AT[AT.target.isin(graph)], AT[~AT.target.isin(graph)]
    units = pd.concat([E, AT_in], ignore_index=True)
    assert not units.target.duplicated().any()
    sc_path = cr.resolve("rpe1_raw_singlecell_01.h5ad")
    rng = np.random.default_rng(SEED)
    with h5py.File(sc_path, "r") as f:
        X = f["X"]
        n, G = X.shape
        gt = np.asarray(read_obs_column(f, "gene_transcript")).astype(str)
        umi = np.asarray(read_obs_column(f, "UMI_count")).astype(np.float64)
        bc = np.asarray(read_obs_column(f, "cell_barcode")).astype(str)
        gname = np.asarray(var_gene_names(f)).astype(str)
        gid = var_column(f, "gene_id")
        pbg = pd.read_csv(PB / "pseudobulk_genes.csv").gene_name.values.astype(str)
        assert (gname == pbg).all(), "gene order differs from the pseudobulk"
        ctrl = np.flatnonzero(pd.Series(gt).str.contains("non-targeting", case=False).values)
        sel, cond, unit_n = [ctrl], ["ctrl"] * len(ctrl), [0] * len(ctrl)
        idx_by_gt = pd.Series(np.arange(n)).groupby(gt).apply(np.asarray)
        for _, u in units.iterrows():
            cells = idx_by_gt[u.gene_transcript]
            assert len(cells) == int(u.ncells), f"{u.target}: {len(cells)} cells vs {u.ncells}"
            k = min(CAP, len(cells))
            pick = np.sort(rng.choice(cells, k, replace=False))
            sel.append(pick)
            cond += [f"{u.target}+ctrl"] * k
            unit_n += [int(u.ncells)] * k
        sel = np.concatenate(sel)
        order = np.argsort(sel, kind="stable")
        rows = sel[order]
        blocks = []
        CH = 4000
        for s in range(0, n, CH):
            e = min(n, s + CH)
            m = (rows >= s) & (rows < e)
            if not m.any():
                continue
            r = rows[m]
            c = X[s:e][r - s].astype(np.float64)
            blocks.append(sp.csr_matrix(np.log1p(c * (1e4 / umi[r])[:, None]).astype(np.float32)))
            if (s // CH) % 10 == 0:
                print(time.strftime("%H:%M:%S"), f"rows {e}/{n}", flush=True)
        Xs = sp.vstack(blocks, format="csr")
        inv = np.empty_like(order)
        inv[order] = np.arange(len(order))
        Xs = Xs[inv]                                    # back to the (ctrl, unit by unit) order of `sel`
    obs = pd.DataFrame({"condition": cond, "cell_type": "rpe1", "gene_transcript": gt[sel], "cell_barcode": bc[sel],
                        "n_cells_unit_total": unit_n}, index=[f"c{i}" for i in range(len(sel))])
    var = pd.DataFrame({"gene_name": gname}, index=pd.Index(gid, name="gene_id"))
    assert var.index.is_unique
    A = ad.AnnData(X=Xs, obs=obs, var=var)
    A.obs["condition"] = A.obs["condition"].astype("category")
    tmp = OUT / "data" / "rpe1_cap16.h5ad.part"
    A.write_h5ad(tmp, compression="gzip")
    os.replace(tmp, OUT / "data" / "rpe1_cap16.h5ad")
    # sanity: for 5 units, the mean over ALL their cells must reproduce the stored pseudobulk row
    pbr = pd.read_csv(PB / "pseudobulk_rows.csv")
    pbx = np.load(PB / "pseudobulk_logmean.npy", mmap_mode="r")
    chk = []
    with h5py.File(sc_path, "r") as f:
        for t in units.target.values[[0, 100, 500, 1000, len(units) - 1]]:
            u = units[units.target == t].iloc[0]
            cells = idx_by_gt[u.gene_transcript]
            x = np.vstack([f["X"][int(i)] for i in cells]).astype(np.float64)
            m_ = np.log1p(x * (1e4 / umi[cells])[:, None]).mean(0)
            j = int(np.flatnonzero(pbr.gene_transcript.values == u.gene_transcript)[0])
            chk.append({"target": t, "max_abs_diff_vs_pseudobulk": float(np.abs(m_.astype(np.float32) - pbx[j]).max())})
    # priors (GEARS reads these from its data_path; pre-placing them avoids every dataverse download)
    shutil.copy2(raw("gene2go_all.pkl"), OUT / "priors" / "gene2go_all.pkl")
    shutil.copy2(raw("essential_all_data_pert_genes.pkl"), OUT / "priors" / "essential_all_data_pert_genes.pkl")
    with tarfile.open(raw("go_essential_all.tar.gz")) as tf:
        tf.extractall(OUT / "priors")
    assert (OUT / "priors" / "go_essential_all" / "go_essential_all.csv").exists()
    # splits: the 21 frozen files (copied unchanged) + GEARS set2conditions derived from them
    sdir = SPLITS / "v0.7_conf" / "rpe1"
    at_conds = [f"{t}+ctrl" for t in AT_in.target]
    for p in sorted(sdir.glob("conf_split_*.json")):
        shutil.copy2(p, OUT / "splits" / p.name)
        s = json.loads(p.read_text())
        s2c = {"train": [f"{t}+ctrl" for t in s["train"]] + at_conds + ["ctrl"],
               "val": [f"{t}+ctrl" for t in s["val"]],
               "test": [f"{t}+ctrl" for t in s["cal"] + s["test"]]}
        allc = set(s2c["train"]) | set(s2c["val"]) | set(s2c["test"])
        assert len(allc) == sum(len(v) for v in s2c.values()), "a condition is in two folds"
        assert set(A.obs.condition.cat.categories) == allc, "split conditions differ from the data"
        (OUT / "splits" / f"set2conditions_{s['split']:02d}.json").write_text(json.dumps(s2c, indent=0))
    for c in ("gears_server_run.py", "setup_gears_env.sh", "requirements_gears.txt", "README_server.md"):
        if (Path(__file__).parent / "gears_server" / c).exists():
            shutil.copy2(Path(__file__).parent / "gears_server" / c, OUT / "code" / c)
    info = {"protocol": "v0.6 section 10 (inherited by v0.7)", "cap_cells_per_perturbation": CAP, "seed": SEED,
            "source": str(sc_path.name), "source_sha256_from_manifest": "see manifests/raw_manifest_additions_v0.6.json",
            "n_cells": int(A.n_obs), "n_ctrl": int(len(ctrl)), "n_genes": int(A.n_vars), "nnz": int(Xs.nnz),
            "n_E": int(len(E)), "n_AT_in_graph": int(len(AT_in)), "AT_excluded_not_in_GEARS_graph": AT_out.target.tolist(),
            "n_AT_below_cap": int((AT_in.ncells < CAP).sum()), "n_E_below_cap": int((E.ncells < CAP).sum()),
            "duplicate_gene_names": var.gene_name[var.gene_name.duplicated(keep=False)].tolist(),
            "pseudobulk_reproduction_check": chk, "seconds": round(time.time() - t0)}
    (OUT / "bundle_info.json").write_text(json.dumps(info, indent=1))
    write_sums()
    print(json.dumps({k: v for k, v in info.items() if k != "AT_excluded_not_in_GEARS_graph"}), flush=True)


def write_sums():
    lines = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.name != "SHA256SUMS":
            lines.append(f"{sha256(p)}  {p.relative_to(OUT).as_posix()}")
    (OUT / "SHA256SUMS").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    {"build": build, "sums": write_sums}[sys.argv[1]]()
