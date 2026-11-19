# confirm_run_v07.py -- protocol v0.7 confirmation-phase entry point. NO selection is made on a confirmation set:
# every choice comes from FROZEN (identical to the v0.7 freeze manifest; asserted at start).
#
#   python confirm_run_v07.py c0      --dataset rpe1 [--sc PATH]   C0: file check, pseudobulk + bulk QC, P0 unit
#                                                                    table (no response-derived fields), priors,
#                                                                    read-only splits, response-free power
#   python confirm_run_v07.py stage1  --dataset rpe1               C1 (later, needs user OK): stage-1 caches
#   python confirm_run_v07.py ext     --dataset rpe1               C1: ext arrays (MD-*, Munc-*; B_boot = 50)
#   python confirm_run_v07.py analyze --dataset rpe1               C1: run_dataset(frozen) + confirm_v07.analyze
#   python confirm_run_v07.py run --dataset gwps [--only 0]        C1 resumable (resume_v07; 2026-10-12): stage1, ext,
#                                                                    stage2 per split with completion markers, then
#                                                                    finish; pauses before the 23:30 power cut
#
# Outputs: data/processed/v0.7/<dataset>/ (inputs), splits/v0.7_conf/<dataset>/ (read-only),
#          results/conf_v0.7/<dataset>/ (results). GWPS uses --sc (file on the user's exFAT drive) and --passes 2.
import argparse
import hashlib
import json
import os
import stat
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import MANIFESTS, PROCESSED, RESULTS, SPLITS, raw  # noqa: E402

FROZEN = {"h": "sim_prior|ridge", "d": "ridge", "s": "ridge",
          "hypers": {"M1": "0.0", "M2": "3", "M3": "quad", "Mdist-N": "1.0", "Mdist-M": "3", "Mscore-N": "0.0",
                     "Mcomb-N": "0.0", "Mcomb-3": "lin2", "oracle": "0.0", "MD-GO": "5", "MD-STR": "5",
                     "Munc-N": "1.0", "Munc-M": "3"},
          "primary": "M1", "best_dist": "Mdist-M"}
DATASETS = {
    "rpe1": {"sc": "rpe1_raw_singlecell_01.h5ad", "bulk": "rpe1_raw_bulk_01.h5ad", "passes": 1, "exclude": None},
    "gwps": {"sc": None, "bulk": "K562_gwps_raw_bulk_01.h5ad", "passes": 2, "exclude": "k562_essential"},
}
DEV_N_E = 1786
DEV_SE = {"E1_step1": 0.020434, "E2_step1": 0.023378, "E1_step2": 0.015927, "E2_step2": 0.013659}   # dev SEs under randomised matching (same file as confirm_v07.DEV_ESTIMATES)
LOW_POWER_STOP = 0.5           # C0 stop rule (v0.7 section 3): step-2 power at the dev point estimate
MIN_FREE_GB = 2.5


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def avail_gb():
    import psutil
    return psutil.virtual_memory().available / 1024 ** 3


def dirs(ds):
    return PROCESSED / "v0.7" / ds, SPLITS / "v0.7_conf" / ds, RESULTS / "conf_v0.7" / ds


def logger(path):
    path.parent.mkdir(parents=True, exist_ok=True)

    def log(msg):
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
        print(line, flush=True)
        with open(path, "a", encoding="utf8") as fh:
            fh.write(line + "\n")
    return log


def check_frozen():
    fz = json.loads((MANIFESTS / "freeze_v0.7.json").read_text(encoding="utf8"))
    assert fz["frozen_selection"] == FROZEN, "FROZEN differs from the v0.7 freeze manifest"
    return fz


def resolve(name):
    """raw_manifest.json first; files registered later (v0.6 downloads) are in raw_manifest_additions_v0.6.json,
    whose 'location' is stored relative to the download root DL (post-freeze code fix 2026-10-10)."""
    from paths import DL, RAW_MANIFEST
    try:
        return raw(name)
    except KeyError:
        add = json.loads((RAW_MANIFEST.parent / "raw_manifest_additions_v0.6.json").read_text(encoding="utf8"))
        loc = Path(add[name]["location"])
        try:
            return DL / loc.relative_to(Path(r"G:\\biology"))
        except ValueError:
            return loc


def check_file(name, fz, path=None, log=print):
    p = Path(path) if path else resolve(name)
    ref = fz["data_files"][name]
    size = p.stat().st_size
    assert size == ref["size_bytes"], f"{name}: size {size} != {ref['size_bytes']}"
    h = sha256(p)
    assert h == ref["sha256"], f"{name}: sha256 mismatch"
    log(f"{name}: size and sha256 match the freeze manifest")
    return p


def response_free_power(n_e, dev_est, margin1=0.05, margin2=0.03):
    """Power of each endpoint at 1.0 / 0.75 / 0.5 x the dev estimate; SE scaled by sqrt(DEV_N_E / n_E)
    (gene-level clustering; uses only the number of E units)."""
    from confirm_v07 import power_at_se
    f = np.sqrt(DEV_N_E / n_e)
    out = {}
    for k, se in DEV_SE.items():
        mg = margin1 if k.endswith("step1") else margin2
        out[k] = {"se_scaled": se * f,
                  **{f"{x:.2f}x_dev": power_at_se(x * dev_est[k], se * f, mg) for x in (1.0, 0.75, 0.5)}}
    return out


C0_MIN_MINUTES = 180          # pseudobulk (two passes over the 65.8 GB GWPS file) only with >= 3 h to the power cut
EST_S = {"stage1": 600, "ext": 300, "stage2": 1200, "finish": 3600}   # first-split estimates (s); later: median


def _pseudobulk_valid(pdir, log):
    """pseudobulk_build.json present AND every output file listed there has the recorded sha256."""
    bj = pdir / "pseudobulk_build.json"
    if not bj.exists():
        return False
    try:
        rec = json.loads(bj.read_text(encoding="utf8"))["outputs_sha256"]
    except (ValueError, KeyError):
        return False
    for name, h in rec.items():
        q = pdir / name
        if not q.exists() or sha256(q) != h:
            log(f"pseudobulk output {name} missing or differs from pseudobulk_build.json -> rebuild")
            return False
    log("pseudobulk outputs verified against pseudobulk_build.json; not rebuilt")
    return True


def load_lowmem(ds):
    import lowmem_v07 as lm
    pdir, _, _ = dirs(ds)
    U = pd.read_csv(pdir / "units_all_targets.csv")
    D = lm.load_dataset_lowmem(U, pdir, pdir / "string_expdb_universe.parquet")
    lm.install(dirs(ds)[2] / "_pred_memmap")      # round 2: chunked row-wise / per-gene computations (same values)
    return D


def do_run(ds, only=None, finish_f1=False, finish_f2=False, finish_f2c=False):
    """Resumable confirmation run (resume_v07): stage1 -> ext -> stage2 per split -> finish. Exit 10 = paused
    before the power cut (rerun the same command later). With --only, runs stage1/ext/stage2 for those splits
    and reports timing and peak memory (no finish)."""
    import resume_v07 as rs
    from eval_v06 import CONF_SEED_BASE
    check_frozen()
    _, _, odir = dirs(ds)
    log = logger(odir / "c1_resume.log")
    log(f"available memory at start {avail_gb():.2f} GB; process memory cap {rs.MEM_CAP_GB} GB; "
        f"{rs.minutes_to_cut():.0f} min to the power cut")
    if avail_gb() < MIN_FREE_GB:
        log(f"STOP: available memory {avail_gb():.2f} GB < {MIN_FREE_GB} GB")
        return 3
    peak = rs.start_guards(log)
    t = {}
    D = load_lowmem(ds)                       # lowmem_v07: identical values, smaller memory (2026-10-12)
    check_splits(ds, D)
    for name, fn in (("stage1", lambda: rs.stage1(D, CONF_SEED_BASE, odir / "cache", log, EST_S["stage1"], only)),
                     ("ext", lambda: rs.ext(D, odir / "cache", odir / "ext_cache", log, EST_S["ext"], only))):
        t0 = time.time()
        rc = fn()
        t[name] = time.time() - t0
        if rc:
            log(f"{name}: paused (exit {rc}); peak RSS {peak['gb']:.2f} GB")
            return rc
    del D
    t0 = time.time()
    rc = rs.stage2(FROZEN, odir / "cache", odir / "ext_cache", odir / "stage2_parts", log, EST_S["stage2"], only)
    t["stage2"] = time.time() - t0
    if rc:
        return rc
    if only is None:
        t0 = time.time()
        if finish_f2c:
            import finish_f2_v07 as f2
            rc = f2.finish_f2c(FROZEN, odir, odir / "cache", odir / "ext_cache", odir / "stage2_parts", log, EST_S["finish"])
        elif finish_f2:
            import finish_f2_v07 as f2
            rc = f2.finish_f2(FROZEN, odir, odir / "cache", odir / "ext_cache", odir / "stage2_parts", log, EST_S["finish"])
        else:
            rc = rs.finish(FROZEN, odir, odir / "cache", odir / "ext_cache", odir / "stage2_parts", log, EST_S["finish"],
                           f1_intern=finish_f1)
        t["finish"] = time.time() - t0
    rs.write_bytes_atomic(odir / f"run_status{'_only' if only else ''}.json",
                          json.dumps({"seconds": t, "peak_rss_gb": peak["gb"], "only": only, "rc": rc,
                                      "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=1).encode())
    log(f"run done rc={rc} seconds={json.dumps({k: round(v) for k, v in t.items()})} peak RSS {peak['gb']:.2f} GB")
    return rc


def do_c0(ds, sc_path=None, chunk_rows=4000):
    import data_v07 as dv
    import pseudobulk_v06 as pb6
    from confirm_v07 import DEV_ESTIMATES
    from eval_v06 import CONF_SEED_BASE
    from splits_v06 import repeated_group_kfold
    cfg = DATASETS[ds]
    pdir, sdir, odir = dirs(ds)
    log = logger(odir / "c0.log")
    fz = check_frozen()
    import resume_v07 as rs
    if avail_gb() < MIN_FREE_GB:
        log(f"STOP: available memory {avail_gb():.2f} GB < {MIN_FREE_GB} GB; close other programs")
        return 3
    log(f"available memory at start {avail_gb():.2f} GB; process memory cap {rs.MEM_CAP_GB} GB")
    peak = rs.start_guards(log)
    pdir.mkdir(parents=True, exist_ok=True)
    pb_ok = _pseudobulk_valid(pdir, log)
    if not pb_ok and rs.minutes_to_cut() < C0_MIN_MINUTES:
        log(f"PAUSE: pseudobulk not built and only {rs.minutes_to_cut():.0f} min to the power cut "
            f"(< {C0_MIN_MINUTES}); nothing started")
        return rs.PAUSED
    sc_name = cfg["sc"] or "K562_gwps_raw_singlecell_01.h5ad"
    sc = check_file(sc_name, fz, sc_path, log)
    bulk = check_file(cfg["bulk"], fz, None, log)
    # 1. pseudobulk (+ split halves) and single-cell vs bulk QC
    if not pb_ok:
        for q in pdir.glob("pseudobulk_*"):          # an interrupted earlier attempt: everything is rebuilt
            rs._writable(q)
        r = pb6.build(sc, pdir, bulk, cfg["passes"], chunk_rows, odir / "c0_pseudobulk_status.json", MIN_FREE_GB, 1.0,
                      log=log)
        if r["status"] != "done":
            log(f"STOP: pseudobulk {r['status']}")
            return 3
    pbinfo = json.loads((pdir / "pseudobulk_build.json").read_text(encoding="utf8"))
    for p in pdir.glob("pseudobulk_*"):
        os.chmod(p, stat.S_IREAD)
    # 2. P0 unit table (gene names and cell counts only)
    excl = ()
    if cfg["exclude"] == "k562_essential":
        excl = tuple(pd.read_csv(dv.COORDS_ARCHIVE).target.values)
    U, readout, info = dv.build_units(bulk, pdir / "coords_ensembl_extra.csv", excl, log)
    rs._writable(pdir / "units_all_targets.csv")
    rs._writable(pdir / "string_expdb_universe.parquet")
    rs.write_bytes_atomic(pdir / "units_all_targets.csv", U.to_csv(index=False).encode("utf8"))
    rows = pd.read_csv(pdir / "pseudobulk_rows.csv")
    nsc = U.gene_transcript.map(dict(zip(rows.gene_transcript, rows.n_cells)))
    info["ncells_bulk_equals_singlecell_frac"] = float((nsc.values == U.ncells.values).mean())
    info["readout_genes_identical_bulk_vs_singlecell"] = bool(
        (pd.read_csv(pdir / "pseudobulk_genes.csv").gene_name.values == readout).all())
    # 3. priors coverage (prior-only)
    n_edges = dv.string_universe(U.target.values, pdir / "string_expdb_universe.parquet")
    E = U[U.E]
    ase = np.load(PROCESSED / "priors" / "go_ase100.npz", allow_pickle=True)
    st = pd.read_parquet(pdir / "string_expdb_universe.parquet")
    pool_targets = set(U[U.E | U.always_train].target)
    st = st[st.g1.isin(pool_targets) & st.g2.isin(pool_targets)]
    deg = pd.concat([st.g1, st.g2]).value_counts()
    info["priors"] = {"string_edges_universe": int(n_edges),
                      "E_in_go_embedding": float(E.target.isin(set(ase["genes"])).mean()),
                      "E_with_string_neighbour_in_pool": float(E.target.map(deg).fillna(0).gt(0).mean())}
    # 4. read-only splits (before any response is used)
    Uu = U[U.E | U.always_train].reset_index(drop=True)
    e_idx = np.where(Uu.E.values)[0]
    sp = repeated_group_kfold(Uu.tss_group.values[e_idx], CONF_SEED_BASE)
    sdir.mkdir(parents=True, exist_ok=True)
    sh = {}
    for i, s in enumerate(sp):
        f = sdir / f"conf_split_{i:02d}.json"
        rec = {"split": i, "round": s["round"], "fold": s["fold"], "seed": s["seed"],
               **{k: [str(t) for t in Uu.target.values[e_idx[s[k]]]] for k in ("train", "val", "cal", "test")}}
        try:
            old = json.loads(f.read_text(encoding="utf8")) if f.exists() else None
        except ValueError:                        # truncated by an interruption: never completed, rewrite
            old = None
        if old is None:
            rs.write_bytes_atomic(f, json.dumps(rec).encode("utf8"))
            os.chmod(f, stat.S_IREAD)
        else:
            assert old == rec, f"existing split file {f.name} differs"
        sh[f.name] = sha256(f)
    info["splits"] = {"seed_base": int(CONF_SEED_BASE), "n": len(sp), "sha256": sh,
                      "test_size_range": [min(len(s["test"]) for s in sp), max(len(s["test"]) for s in sp)],
                      "cal_size": len(sp[0]["cal"])}
    # 5. response-free power
    pw = response_free_power(int(U.E.sum()), DEV_ESTIMATES)
    low = min(pw["E1_step2"]["1.00x_dev"], pw["E2_step2"]["1.00x_dev"])
    info["power"] = pw
    info["expected_T5_test_rows"] = int(round(3 * U.E.sum() * 0.2))
    info["c0_stop_low_power"] = bool(low < LOW_POWER_STOP)
    info["pseudobulk"] = {k: pbinfo[k] for k in ("cells", "genes", "rows", "passes", "seconds", "qc_vs_bulk")}
    info["inputs_sha256"] = {"units_all_targets.csv": sha256(pdir / "units_all_targets.csv"),
                             "string_expdb_universe.parquet": sha256(pdir / "string_expdb_universe.parquet"),
                             **pbinfo["outputs_sha256"]}
    for p in (pdir / "units_all_targets.csv", pdir / "string_expdb_universe.parquet"):
        os.chmod(p, stat.S_IREAD)
    info["memory_peak_rss_gb"] = round(peak["gb"], 2)
    rs.write_bytes_atomic(odir / "c0_report.json", json.dumps(info, indent=1, default=str).encode("utf8"))
    log(f"C0 done: E={info['E']} TSS groups={info['E_tss_groups']} QC pass={pbinfo['qc_vs_bulk']['pass']} "
        f"step-2 power at dev estimate={low:.2f} stop_low_power={info['c0_stop_low_power']}")
    return 0 if pbinfo["qc_vs_bulk"]["pass"] else 4


def load(ds):
    import data_v07 as dv
    pdir, _, _ = dirs(ds)
    U = pd.read_csv(pdir / "units_all_targets.csv")
    return dv.load_dataset(U, pdir, pdir / "string_expdb_universe.parquet")


def check_splits(ds, D):
    from eval_v06 import CONF_SEED_BASE
    from splits_v06 import repeated_group_kfold
    _, sdir, _ = dirs(ds)
    U, e_idx = D["U"], D["E_idx"]
    sp = repeated_group_kfold(U.tss_group.values[e_idx], CONF_SEED_BASE)
    for i, s in enumerate(sp):
        rec = json.loads((sdir / f"conf_split_{i:02d}.json").read_text(encoding="utf8"))
        for k in ("train", "val", "cal", "test"):
            assert rec[k] == [str(t) for t in U.target.values[e_idx[s[k]]]], f"split {i} {k} differs from C0 file"


def do_stage1(ds):
    import pipeline_v06 as p6
    from eval_v06 import CONF_SEED_BASE
    check_frozen()
    _, _, odir = dirs(ds)
    log = logger(odir / "c1_stage1.log")
    D = load(ds)
    check_splits(ds, D)
    p6.stage1(D, CONF_SEED_BASE, odir / "cache", 3, log)
    return 0


def do_ext(ds):
    import ext_v07dev as ex
    import pipeline_v06 as p6
    check_frozen()
    _, _, odir = dirs(ds)
    log = logger(odir / "c1_ext.log")
    D = load(ds)
    ed = odir / "ext_cache"
    ed.mkdir(parents=True, exist_ok=True)
    for c in p6.load_caches(odir / "cache"):
        f = ed / f"ext_split_{c['split_id']:02d}.npz"
        if f.exists():
            continue
        e = ex.compute_ext(c, D, log=log)
        np.savez(f, **e)
    return 0


def do_analyze(ds):
    import confirm_v07 as cf
    import ext_v07dev as ex
    import pipeline_v06 as p6
    from eval_v06 import CONF_SEED_BASE
    check_frozen()
    _, _, odir = dirs(ds)
    log = logger(odir / "c1_analyze.log")
    D = load(ds)
    check_splits(ds, D)
    summ, cov, _ = p6.run_dataset(D, CONF_SEED_BASE, odir, n_rounds=3, frozen=FROZEN, log=log,
                                   ext_dir=odir / "ext_cache")
    caches = p6.load_caches(odir / "cache")
    for c in caches:
        z = np.load(odir / "ext_cache" / f"ext_split_{c['split_id']:02d}.npz")
        ex.attach(c, {k: z[k] for k in z.files})
    out, rm = cf.analyze(caches, FROZEN, log=log, cov_grid=cov)
    rm.to_parquet(odir / "confirm_rand_rows.parquet")
    cf.write(out, odir / "confirm_v07_result.json")
    log(f"classification: {out['classification']}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["c0", "stage1", "ext", "analyze", "run"])
    ap.add_argument("--only", type=int, nargs="*", default=None)
    ap.add_argument("--mem-cap", type=float, default=None, help="process RSS cap in GB for 'run' (default 4.0)")
    ap.add_argument("--finish-f1", action="store_true", help="finish with F1 string interning (validated on dev)")
    ap.add_argument("--finish-f2", action="store_true", help="finish with the F2 column store (finish_f2_v07)")
    ap.add_argument("--finish-f2c", action="store_true", help="finish with the F2c code store (finish_f2_v07)")
    ap.add_argument("--dataset", required=True, choices=list(DATASETS))
    ap.add_argument("--sc", default=None)
    ap.add_argument("--chunk-rows", type=int, default=4000)
    a = ap.parse_args()
    if a.mode == "c0":
        return do_c0(a.dataset, a.sc, a.chunk_rows)
    if a.mode == "run":
        if a.mem_cap is not None:
            import resume_v07 as rs
            assert 0 < a.mem_cap <= 4.8, "user rule: the cap may be raised to at most 4.8 GB"
            rs.MEM_CAP_GB = a.mem_cap
        return do_run(a.dataset, a.only, a.finish_f1, a.finish_f2, a.finish_f2c)
    return {"stage1": do_stage1, "ext": do_ext, "analyze": do_analyze}[a.mode](a.dataset)


if __name__ == "__main__":
    raise SystemExit(main())
