# d1_run_v061.py -- protocol v0.6.2 (file name kept from v0.6.1), development phase D1 on K562 essential (never touches confirmation data).
#   python d1_run_v061.py pseudobulk          D1a: split-half pseudobulk of the K562 essential single-cell file
#                                              (pseudobulk_v06, 1 pass) + identity check of the full matrix vs v0.4
#   python d1_run_v061.py stages234 [--time-one]  stages 2-4 on cached splits + sens_v07 analyses
#   python d1_run_v061.py stage1 [--only I]   D1 stage 1: compute_split for the 21 dev splits (or split I only),
#                                              cached in results/d1_v0.6.1/cache (same file names as
#                                              pipeline_v06.stage1, so later stages reuse them)
# Every mode samples memory: process RSS every 0.5 s, Windows peak working set, minimum system available.
# Stops (exit code 3) before starting if available memory < MIN_START_GB, and between splits if < MIN_RUN_GB.
import argparse
import hashlib
import json
import shutil
import os
import pickle
import stat
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PROCESSED, RESULTS, SPLITS, raw  # noqa: E402

MIN_START_GB = 3.5
MIN_RUN_GB = 1.0
PB_DIR = PROCESSED / "v0.6" / "k562ess"
OUT = RESULTS / "d1_v0.6.1"
MANIFESTS = RESULTS.parent / "manifests"
SPLIT_DIR = SPLITS / "v0.6.1_dev"


class MemMonitor:
    def __init__(self):
        self.proc = psutil.Process()
        self.peak_rss = 0
        self.min_avail = float("inf")
        self._stop = threading.Event()
        self.t = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while not self._stop.is_set():
            self.peak_rss = max(self.peak_rss, self.proc.memory_info().rss)
            self.min_avail = min(self.min_avail, psutil.virtual_memory().available)
            time.sleep(0.5)

    def start(self):
        self.t.start()
        return self

    def report(self):
        mi = self.proc.memory_info()
        return {"peak_rss_sampled_gb": self.peak_rss / 1024 ** 3,
                "peak_working_set_gb": getattr(mi, "peak_wset", float("nan")) / 1024 ** 3,
                "min_system_available_gb": self.min_avail / 1024 ** 3,
                "available_now_gb": psutil.virtual_memory().available / 1024 ** 3}

    def stop(self):
        self._stop.set()


def log_to(path):
    path.parent.mkdir(parents=True, exist_ok=True)

    def log(msg):
        line = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + str(msg)
        print(line, flush=True)
        with open(path, "a", encoding="utf8") as fh:
            fh.write(line + "\n")
    return log


def avail_gb():
    return psutil.virtual_memory().available / 1024 ** 3


def threads_env():
    return {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")}


def do_pseudobulk():
    import pseudobulk_v06 as pb
    log = log_to(OUT / "d1a_pseudobulk.log")
    if avail_gb() < MIN_START_GB:
        log(f"STOP: available memory {avail_gb():.2f} GB < {MIN_START_GB} GB")
        return 3
    mon = MemMonitor().start()
    r = pb.build(raw("K562_essential_raw_singlecell_01.h5ad"), PB_DIR, raw("K562_essential_raw_bulk_01.h5ad"),
                 passes=1, chunk_rows=4000, status=OUT / "d1a_status.json", min_free=2.5, pause_free=1.0, log=log)
    if r["status"] != "done":
        log(f"STOP: pseudobulk status {r['status']}")
        return 2
    # identity check of the full matrix against the frozen v0.4 representation
    v4 = PROCESSED / "v0.4"
    rows4 = pd.read_csv(v4 / "pseudobulk_rows.csv")
    rows6 = pd.read_csv(PB_DIR / "pseudobulk_rows.csv")
    g4 = pd.read_csv(v4 / "pseudobulk_genes.csv").gene_name.values
    g6 = pd.read_csv(PB_DIR / "pseudobulk_genes.csv").gene_name.values
    X4 = np.load(v4 / "pseudobulk_logmean.npy", mmap_mode="r")
    X6 = np.load(PB_DIR / "pseudobulk_logmean.npy", mmap_mode="r")
    A = np.load(PB_DIR / "pseudobulk_logmean_A.npy", mmap_mode="r")
    B = np.load(PB_DIR / "pseudobulk_logmean_B.npy", mmap_mode="r")
    ix6 = {g: i for i, g in enumerate(rows6.gene_transcript)}
    same_rows = set(rows4.gene_transcript) == set(rows6.gene_transcript)
    j = np.array([ix6[g] for g in rows4.gene_transcript])
    maxdiff = 0.0
    recon = 0.0
    nA, nB = rows6.n_A.values.astype(float), rows6.n_B.values.astype(float)
    for s in range(0, len(j), 256):
        jj = j[s:s + 256]
        maxdiff = max(maxdiff, float(np.abs(np.asarray(X4[s:s + 256], np.float64) - np.asarray(X6[jj], np.float64)).max()))
        w = (nA[jj][:, None] * np.asarray(A[jj], np.float64) + nB[jj][:, None] * np.asarray(B[jj], np.float64)) / \
            (nA[jj] + nB[jj])[:, None]
        recon = max(recon, float(np.abs(w - np.asarray(X6[jj], np.float64)).max()))
    chk = {"rows_identical_sets": bool(same_rows), "genes_identical_order": bool((g4 == g6).all()),
           "n_cells_identical": bool((rows4.n_cells.values == rows6.n_cells.values[j]).all()),
           "max_abs_diff_full_v06_vs_v04": maxdiff, "max_abs_diff_weighted_halves_vs_full": recon,
           "halves_min_cells": int(min(nA.min(), nB.min())), "halves_balanced": bool((np.abs(nA - nB) <= 1).all()),
           "qc_vs_bulk": r.get("qc_vs_bulk"), "memory": mon.report(), "threads": threads_env()}
    chk["pass"] = bool(same_rows and chk["genes_identical_order"] and chk["n_cells_identical"] and maxdiff <= 1e-5
                       and recon <= 1e-5 and (r.get("qc_vs_bulk") or {}).get("pass", False))
    (OUT / "d1a_check.json").write_text(json.dumps(chk, indent=1, default=float), encoding="utf8")
    for p in PB_DIR.glob("pseudobulk_*"):
        os.chmod(p, stat.S_IREAD)
    log(json.dumps(chk, default=float))
    mon.stop()
    return 0 if chk["pass"] else 2


def write_splits(splits, U, E_idx):
    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    for i, sp in enumerate(splits):
        f = SPLIT_DIR / f"dev_split_{i:02d}.json"
        d = {"split": i, "round": sp["round"], "fold": sp["fold"], "seed": sp["seed"],
             **{k: U.target.values[E_idx[sp[k]]].tolist() for k in ("train", "val", "cal", "test")}}
        if f.exists():
            assert json.loads(f.read_text()) == d, f"existing split file {f.name} differs"
            continue
        f.write_text(json.dumps(d), encoding="utf8")
        os.chmod(f, stat.S_IREAD)


def at_definition_check(i, c, D, log):
    """v0.6.2: E_true of always-train (AT) and learner units recomputed from the written definition using
    readout gene NAMES (independent of U.tcols): mean-model S1(k=50), excluding every readout column whose
    gene name equals the target. Writes counts and the maximum deviation; no per-unit values are saved."""
    genes = pd.read_csv(PB_DIR / "pseudobulk_genes.csv").gene_name.values.astype(object)
    U, Y = D["U"], D["Y"]
    out = {"split": i}
    for nm, idx, ref_vals in (("at", c["at"], c["E_at"]), ("learner", c["lrn"], c["E_pool"])):
        s = (Y[idx] - c["ref_mu"]) / c["ref_sigma"]
        e = np.empty(len(idx))
        n_excl = np.empty(len(idx), int)
        for r, tg in enumerate(U.target.values[idx]):
            a = np.abs(s[r]).copy()
            m = genes == tg
            n_excl[r] = int(m.sum())
            a[m] = -np.inf
            top = np.argsort(-a)[:50]
            e[r] = np.sqrt((s[r, top] ** 2).mean())
        out[nm] = {"n": int(len(idx)), "n_target_in_readout": int((n_excl > 0).sum()),
                   "n_target_not_in_readout": int((n_excl == 0).sum()),
                   "n_excluded_columns_counts": {int(k): int(v) for k, v in zip(*np.unique(n_excl, return_counts=True))},
                   "n_nonfinite": int((~np.isfinite(ref_vals)).sum()),
                   "max_abs_diff_vs_cache": float(np.abs(e - ref_vals).max()),
                   "n_cells_lt30": int((U.ncells.values[idx] < 30).sum())}
    out["pass"] = bool(all(out[k]["max_abs_diff_vs_cache"] <= 1e-10 and out[k]["n_nonfinite"] == 0 for k in ("at", "learner")))
    (OUT / f"at_definition_split{i:02d}.json").write_text(json.dumps(out, indent=1), encoding="utf8")
    log(f"E_true definition check split {i:02d}: {json.dumps(out)}")
    if not out["pass"]:
        raise RuntimeError(f"E_true definition check failed for split {i:02d}")


def do_stage1(only):
    import eval_v06 as ev
    import pipeline_v06 as p6
    from splits_v06 import repeated_group_kfold
    log = log_to(OUT / "d1_stage1.log")
    if avail_gb() < MIN_START_GB:
        log(f"STOP: available memory {avail_gb():.2f} GB < {MIN_START_GB} GB")
        return 3
    mon = MemMonitor().start()
    t0 = time.time()
    D = p6.load_dev_v06(PB_DIR)
    log(f"loaded dev data: units={len(D['U'])} E={len(D['E_idx'])} AT={len(D['AT_idx'])} genes={D['Y'].shape[1]} "
        f"in {time.time() - t0:.0f}s; memory {json.dumps(mon.report())}; threads {threads_env()}")
    splits = repeated_group_kfold(D["U"].tss_group.values[D["E_idx"]], ev.DEV_SEED_BASE, n_rounds=3)
    assert len(splits) == 21
    write_splits(splits, D["U"], D["E_idx"])
    cache = OUT / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    todo = [only] if only is not None else range(len(splits))
    per = []
    for i in todo:
        f = cache / f"split_{i:02d}.pkl"
        if f.exists():
            log(f"split {i:02d} cached, skipped")
            continue
        if avail_gb() < MIN_RUN_GB:
            log(f"STOP before split {i:02d}: available memory {avail_gb():.2f} GB < {MIN_RUN_GB} GB")
            return 3
        ts = time.time()
        c = ev.compute_split(i, splits[i], D)
        tmp = f.with_suffix(".tmp")
        with open(tmp, "wb") as fh:
            pickle.dump(c, fh)
        tmp.replace(f)
        at_definition_check(i, c, D, log)
        m = mon.report()
        per.append({"split": i, "round": splits[i]["round"], "fold": splits[i]["fold"], "seconds": time.time() - ts,
                    "K": c["chosen"]["K"], "d": c["chosen"]["d"], "n_test": len(c["test"]), "n_cal": len(c["cal"]),
                    "n_val": len(c["val"]), "n_pool": len(c["pool"]), "peak_rss_gb_so_far": m["peak_rss_sampled_gb"],
                    "min_avail_gb_so_far": m["min_system_available_gb"], "cache_mb": f.stat().st_size / 1e6})
        log(json.dumps(per[-1], default=float))
    meta = {"mode": "only" if only is not None else "all", "splits_run": per, "memory": mon.report(),
            "seconds_total": time.time() - t0, "threads": threads_env()}
    name = f"stage1_meta_split{only:02d}.json" if only is not None else "stage1_meta_all.json"
    (OUT / name).write_text(json.dumps(meta, indent=1, default=float), encoding="utf8")
    log(f"done: {json.dumps(meta['memory'])}")
    mon.stop()
    return 0


def do_stages234(time_one):
    """Stages 2-4 of pipeline_v06.run_dataset on the cached dev splits (stage 1 is not recomputed), then the
    v0.7-style sensitivity analyses and M_dist strata (sens_v07). --time-one: same on split 00 only, in a
    separate folder, to time the run and measure memory before the full run."""
    import shutil
    import eval_v06 as ev
    import pipeline_v06 as p6
    import sens_v07 as sv
    out = OUT / "_timing_split00" if time_one else OUT
    log = log_to(out / "d1_stages234.log")
    if avail_gb() < MIN_START_GB:
        log(f"STOP: available memory {avail_gb():.2f} GB < {MIN_START_GB} GB")
        return 3
    if time_one:
        (out / "cache").mkdir(parents=True, exist_ok=True)
        if not (out / "cache" / "split_00.pkl").exists():
            shutil.copy2(OUT / "cache" / "split_00.pkl", out / "cache" / "split_00.pkl")
    n_cache = len(list((out / "cache").glob("split_*.pkl")))
    assert n_cache == (1 if time_one else 21), f"expected {(1 if time_one else 21)} caches, found {n_cache}"
    p6.stage1 = lambda *a, **k: None        # stage 1 is complete (verified caches); never recomputed here
    mon = MemMonitor().start()
    t0 = time.time()
    D = p6.load_dev_v06(PB_DIR)
    summ, cov, pt = p6.run_dataset(D, ev.DEV_SEED_BASE, out, n_rounds=3, log=log)
    t_main = time.time() - t0
    cov_mem = float(cov.memory_usage(deep=True).sum() / 1024 ** 3)
    n_rows = int(len(cov))
    del cov, pt, D
    log(f"stages 2-4 done in {t_main:.0f}s; coverage rows {n_rows} ({cov_mem:.2f} GB in memory); {json.dumps(mon.report())}")
    t1 = time.time()
    sv.run_all(out, p6.load_caches(out / "cache"), log=log)
    t_sens = time.time() - t1
    meta = {"mode": "time_one" if time_one else "all", "seconds_stages234": t_main, "seconds_sens": t_sens,
            "coverage_rows": n_rows, "coverage_rows_mem_gb": cov_mem, "memory": mon.report(), "threads": threads_env(),
            "selection": summ["selection"]}
    (out / "stages234_meta.json").write_text(json.dumps(meta, indent=1, default=float), encoding="utf8")
    log(f"done: {json.dumps(meta['memory'])} sens {t_sens:.0f}s")
    mon.stop()
    return 0


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def do_ext():
    """v0.7 development: ext arrays (raw-distance Mondrian inputs + bootstrap model-uncertainty proxy) for the 21
    verified stage-1 caches -> results/d1_v0.7dev/ext_cache/ext_split_NN.npz (resumable)."""
    import ext_v07dev as ex
    import pipeline_v06 as p6
    out = RESULTS / "d1_v0.7dev"
    ed = out / "ext_cache"
    ed.mkdir(parents=True, exist_ok=True)
    log = log_to(out / "d1_ext.log")
    if avail_gb() < MIN_START_GB:
        log(f"STOP: available memory {avail_gb():.2f} GB < {MIN_START_GB} GB")
        return 3
    man = json.loads((MANIFESTS / "d1_stage1_cache_v0.6.2.json").read_text(encoding="utf8"))
    want = {r["split"]: r["sha256"] for r in man["files"]}
    cd = out / "cache"
    cd.mkdir(exist_ok=True)
    for i, h in want.items():               # copy the verified stage-1 caches (never modified)
        src, dst = OUT / "cache" / f"split_{i:02d}.pkl", cd / f"split_{i:02d}.pkl"
        if not dst.exists():
            shutil.copy2(src, dst)
        assert sha256(dst) == h, f"cache {i:02d} hash differs from manifest"
    log("21 stage-1 caches copied and verified against manifests/d1_stage1_cache_v0.6.2.json")
    mon = MemMonitor().start()
    D = p6.load_dev_v06(PB_DIR)
    t0 = time.time()
    for c in p6.load_caches(cd):
        f = ed / f"ext_split_{c['split_id']:02d}.npz"
        if f.exists():
            continue
        t1 = time.time()
        e = ex.compute_ext(c, D, log=log)
        tmp = ed / f"ext_split_{c['split_id']:02d}.tmp.npz"
        np.savez(tmp, **e)
        tmp.replace(f)
        log(f"split {c['split_id']:02d} {time.time() - t1:.0f}s {json.dumps(mon.report())}")
    log(f"ext done in {time.time() - t0:.0f}s")
    mon.stop()
    return 0


def do_v07dev():
    """v0.7 development rerun of stages 2-4 (+ sens_v07) with the additional baselines."""
    import eval_v06 as ev
    import pipeline_v06 as p6
    import sens_v07 as sv
    out = RESULTS / "d1_v0.7dev"
    log = log_to(out / "d1_v07dev_stages234.log")
    if avail_gb() < MIN_START_GB:
        log(f"STOP: available memory {avail_gb():.2f} GB < {MIN_START_GB} GB")
        return 3
    n_ext = len(list((out / "ext_cache").glob("ext_split_*.npz")))
    assert len(list((out / "cache").glob("split_*.pkl"))) == 21 and n_ext == 21, f"caches/ext incomplete ({n_ext})"
    p6.stage1 = lambda *a, **k: None
    mon = MemMonitor().start()
    t0 = time.time()
    D = p6.load_dev_v06(PB_DIR)
    summ, cov, pt = p6.run_dataset(D, ev.DEV_SEED_BASE, out, n_rounds=3, log=log, ext_dir=out / "ext_cache")
    t_main = time.time() - t0
    n_rows = int(len(cov))
    del cov, pt, D
    log(f"stages 2-4 done in {t_main:.0f}s; coverage rows {n_rows}; {json.dumps(mon.report())}")
    t1 = time.time()
    caches = p6.load_caches(out / "cache")
    import ext_v07dev as ex
    for c in caches:
        z = np.load(out / "ext_cache" / f"ext_split_{c['split_id']:02d}.npz")
        ex.attach(c, {k: z[k] for k in z.files})
    sv.run_all(out, caches, log=log)
    meta = {"mode": "v07dev", "seconds_stages234": t_main, "seconds_sens": time.time() - t1, "coverage_rows": n_rows,
            "memory": mon.report(), "threads": threads_env(), "selection": summ["selection"]}
    (out / "stages234_meta.json").write_text(json.dumps(meta, indent=1, default=float), encoding="utf8")
    log(f"done: {json.dumps(meta['memory'])}")
    mon.stop()
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["pseudobulk", "stage1", "stages234", "ext", "v07dev"])
    ap.add_argument("--only", type=int, default=None)
    ap.add_argument("--time-one", action="store_true")
    a = ap.parse_args()
    if a.mode == "pseudobulk":
        raise SystemExit(do_pseudobulk())
    if a.mode == "stage1":
        raise SystemExit(do_stage1(a.only))
    if a.mode == "ext":
        raise SystemExit(do_ext())
    if a.mode == "v07dev":
        raise SystemExit(do_v07dev())
    raise SystemExit(do_stages234(a.time_one))
