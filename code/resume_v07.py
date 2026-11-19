# resume_v07.py -- crash- and power-cut-safe execution of the v0.7 confirmation stages (engineering only,
# 2026-10-12, user request for GWPS). NO algorithm, parameter or random-number change: the same frozen functions
# (eval_v06.compute_split, ext_v07dev.compute_ext, eval_v06.apply_methods / proxy_table, pipeline_v06.summarise,
# confirm_v07.analyze) are called with the same arguments in the same order as pipeline_v06.run_dataset /
# confirm_run_v07; only how outputs are written and when work is started changes.
#
# Per-split outputs are written to <file>.tmp, flushed + fsync'ed, then atomically renamed; afterwards a completion
# marker <file>.done.json (sha256, bytes, seconds) is written the same way. On restart a split counts as done only
# if file AND marker exist and the file's sha256 equals the marker; anything else (missing marker, truncated or
# stale file, leftover .tmp) is recomputed and overwritten (no file is deleted).
# Guards: (1) a split / block is started only if it can finish before the daily power cut (23:30 local, margin 15 min,
# duration = 1.5 x median of completed splits or the given estimate); (2) a watchdog ends the process if its RSS
# exceeds MEM_CAP_GB (4 GB); (3) the machine is kept awake while running (Windows SetThreadExecutionState).
import hashlib
import io
import json
import os
import pickle
import stat
import threading
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

POWER_CUT = (23, 30)          # local time, power back 06:30
POWER_BACK = (6, 30)
MARGIN_MIN = 15
SAFETY = 1.5             # user 2026-10-12: >= 1.5 x single-split duration before the cut
MEM_CAP_GB = 4.0
PAUSED = 10                   # exit code: paused before the power cut (rerun later; nothing lost)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _writable(path):
    if path.exists() and not os.access(path, os.W_OK):
        os.chmod(path, stat.S_IWRITE | stat.S_IREAD)


def _commit(tmp, path):
    _writable(path)
    os.replace(tmp, path)


def write_bytes_atomic(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    _commit(tmp, path)


def pickle_atomic(obj, path):
    """Same bytes as pickle.dump(obj, fh) (the call used by pipeline_v06.stage1), plus fsync + atomic rename."""
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh)
        fh.flush()
        os.fsync(fh.fileno())
    _commit(tmp, path)


def npz_atomic(arrays, path):
    buf = io.BytesIO()
    np.savez(buf, **arrays)
    write_bytes_atomic(path, buf.getvalue())


def mark_done(path, seconds, **extra):
    path = Path(path)
    rec = {"file": path.name, "sha256": sha256(path), "bytes": path.stat().st_size, "seconds": float(seconds),
           "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **extra}
    write_bytes_atomic(path.with_name(path.name + ".done.json"), json.dumps(rec, indent=1).encode("utf8"))
    return rec


def is_done(path):
    """Marker dict if <path> and its marker exist and the sha256 matches; otherwise None."""
    path = Path(path)
    mk = path.with_name(path.name + ".done.json")
    if not (path.exists() and mk.exists()):
        return None
    try:
        rec = json.loads(mk.read_text(encoding="utf8"))
    except (ValueError, OSError):
        return None
    if rec.get("bytes") != path.stat().st_size or rec.get("sha256") != sha256(path):
        return None
    return rec


def minutes_to_cut(now=None):
    """Minutes until the next power cut; 0 inside the cut window."""
    now = now or datetime.now()
    cut = now.replace(hour=POWER_CUT[0], minute=POWER_CUT[1], second=0, microsecond=0)
    back = now.replace(hour=POWER_BACK[0], minute=POWER_BACK[1], second=0, microsecond=0)
    if now >= cut or now < back:
        return 0.0
    return (cut - now).total_seconds() / 60


def can_start(est_seconds, what, log):
    left = minutes_to_cut() - MARGIN_MIN
    need = SAFETY * float(est_seconds) / 60
    ok = left >= need
    log(f"time check for {what}: {left:.0f} min usable before the power cut, need {need:.0f} min -> "
        f"{'start' if ok else 'PAUSE (not started)'}")
    return ok


def start_guards(log, cap_gb=None, every_s=2.0):
    """Memory watchdog (RSS of this process > cap -> log and hard exit 7; in-progress output is only a .tmp file
    and is recomputed on restart) and keep-awake for the lifetime of the process."""
    try:
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)
    except Exception:
        pass
    cap_gb = MEM_CAP_GB if cap_gb is None else cap_gb
    proc = psutil.Process()
    peak = {"gb": 0.0}

    def watch():
        while True:
            rss = proc.memory_info().rss / 1024 ** 3
            peak["gb"] = max(peak["gb"], rss)
            if rss > cap_gb:
                log(f"STOP: process RSS {rss:.2f} GB > cap {cap_gb} GB (memory watchdog); partial output discarded")
                os._exit(7)
            time.sleep(every_s)
    threading.Thread(target=watch, daemon=True).start()
    return peak


def intern_strings(df, canon):
    """F1 (2026-10-13, engineering only): make every all-str object column of df reference one shared Python
    str object per distinct value (via pd.factorize + take). Values, dtypes, column order and index are
    unchanged; only object identity (and hence memory) changes. Columns holding any non-str non-NA value are
    left untouched (factorize would merge e.g. 1 and True)."""
    for c in df.columns:
        a = df[c].values
        if a.dtype != object:
            continue
        nn = pd.notna(a)
        if not all(type(v) is str for v in a[nn]):
            continue
        codes, uniq = pd.factorize(a, use_na_sentinel=True)
        u = np.empty(len(uniq), dtype=object)
        for j, v in enumerate(uniq):
            u[j] = canon.setdefault(v, v)
        out = u[codes] if len(u) else a.copy()
        na = codes < 0
        if na.any():
            out[na] = a[na]
        df[c] = out
    return df


# ------------------------------------------------------------------ stages

def stage1(D, seed_base, cache_dir, log, est_s=None, only=None):
    """= pipeline_v06.stage1 (same splits, same compute_split calls, same pickle bytes) with markers + guards."""
    import eval_v06 as ev
    from splits_v06 import repeated_group_kfold
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    splits = repeated_group_kfold(D["U"].tss_group.values[D["E_idx"]], seed_base, n_rounds=3)
    durs = []
    for i, sp in enumerate(splits):
        if only is not None and i not in only:
            continue
        f = cache_dir / f"split_{i:02d}.pkl"
        rec = is_done(f)
        if rec:
            durs.append(rec["seconds"])
            continue
        if f.exists():
            log(f"stage1 split {i:02d}: file without a valid completion marker -> recomputed")
        est = float(np.median(durs)) if durs else est_s
        if est is not None and not can_start(est, f"stage1 split {i:02d}", log):
            return PAUSED
        t0 = time.time()
        c = ev.compute_split(i, sp, D)
        pickle_atomic(c, f)
        mark_done(f, time.time() - t0)
        durs.append(time.time() - t0)
        log(f"stage1 split {i:02d} round {sp['round']} fold {sp['fold']} K={c['chosen']['K']} d={c['chosen']['d']} "
            f"n_test={len(c['test'])} {time.time() - t0:.1f}s")
        del c
    return 0


def _cache_files(cache_dir, n=21):
    fs = [Path(cache_dir) / f"split_{i:02d}.pkl" for i in range(n)]
    bad = [f.name for f in fs if not is_done(f)]
    assert not bad, f"stage-1 caches missing or incomplete: {bad}"
    return fs


def ext(D, cache_dir, ext_dir, log, est_s=None, only=None):
    """= confirm_run_v07.do_ext (compute_ext per cache, same arrays) with markers + guards; one cache in memory."""
    import ext_v07dev as ex
    ext_dir = Path(ext_dir)
    ext_dir.mkdir(parents=True, exist_ok=True)
    durs = []
    for f in _cache_files(cache_dir) if only is None else [Path(cache_dir) / f"split_{i:02d}.pkl" for i in only]:
        i = int(f.stem.split("_")[1])
        out = ext_dir / f"ext_split_{i:02d}.npz"
        rec = is_done(out)
        if rec:
            durs.append(rec["seconds"])
            continue
        est = float(np.median(durs)) if durs else est_s
        if est is not None and not can_start(est, f"ext split {i:02d}", log):
            return PAUSED
        t0 = time.time()
        with open(f, "rb") as fh:
            c = pickle.load(fh)
        e = ex.compute_ext(c, D, log=log)
        npz_atomic(e, out)
        mark_done(out, time.time() - t0)
        durs.append(time.time() - t0)
        del c, e
    return 0


def _load_cache_with_ext(f, ext_dir):
    import ext_v07dev as ex
    with open(f, "rb") as fh:
        c = pickle.load(fh)
    z = np.load(Path(ext_dir) / f"ext_split_{c['split_id']:02d}.npz")
    ex.attach(c, {k: z[k] for k in z.files})
    return c


def stage2(frozen, cache_dir, ext_dir, parts_dir, log, est_s=None, only=None):
    """Stage 2 of pipeline_v06.run_dataset in frozen (confirmation) mode, one split at a time. Per split the same
    calls in the same order: proxy_table, then the main apply_methods block, then the M_perm block (which
    run_dataset runs after all splits; apply_methods is a pure function of (cache, arguments) apart from the
    process-global M3 fallback counter, whose per-split increments are stored and summed)."""
    import eval_v06 as ev
    import methods_v06 as mv
    import pipeline_v06 as p6
    parts_dir = Path(parts_dir)
    parts_dir.mkdir(parents=True, exist_ok=True)
    sel = {"h": tuple(frozen["h"].split("|")), "d": frozen["d"], "s": frozen["s"]}
    P, prim = ev.PRIMARY_MODEL, frozen["primary"]
    prim_h = frozen["hypers"][prim]
    durs = []
    files = _cache_files(cache_dir) if only is None else [Path(cache_dir) / f"split_{i:02d}.pkl" for i in only]
    for f in files:
        i = int(f.stem.split("_")[1])
        out = parts_dir / f"stage2_split_{i:02d}.pkl"
        rec = is_done(out)
        if rec:
            durs.append(rec["seconds"])
            continue
        est = float(np.median(durs)) if durs else est_s
        if est is not None and not can_start(est, f"stage2 split {i:02d}", log):
            return PAUSED
        t0 = time.time()
        c = _load_cache_with_ext(f, ext_dir)
        before = dict(mv.M3_FALLBACKS)
        pt = ev.proxy_table(c)
        main = []
        for m in ev.MODELS:
            for k in ev.KS:
                main.append(ev.apply_methods(c, m, k, sel))
            main.append(ev.apply_methods(c, m, ev.PRIMARY_K, sel, variant="half"))
        for m in (P, ev.SECONDARY_MODEL):
            main.append(ev.apply_methods(c, m, ev.PRIMARY_K, sel, mode="matched"))
        main.append(ev.apply_methods(c, P, ev.PRIMARY_K, sel, mode="matched_cal"))
        perm = []
        for r in range(ev.N_PERM):
            for mode in ("nominal", "matched"):
                perm.append(ev.apply_methods(c, P, ev.PRIMARY_K, sel, grid=[(prim, p6._hyper_value(prim, prim_h))],
                                             mode=mode, perm_r=r))
        fb = {k: int(mv.M3_FALLBACKS.get(k, 0) - before.get(k, 0)) for k in mv.M3_FALLBACKS}
        part = {"split_id": c["split_id"], "pt": pt, "main": main, "perm": perm, "m3_fallbacks": fb,
                "covariates": c["covariates"].assign(split=c["split_id"]),
                "K": c["chosen"]["K"], "d": c["chosen"]["d"]}
        pickle_atomic(part, out)
        mark_done(out, time.time() - t0, rows=int(sum(len(x) for x in main + perm)))
        durs.append(time.time() - t0)
        log(f"stage2 split {i:02d}: {sum(len(x) for x in main + perm)} rows {time.time() - t0:.0f}s")
        del c, part, main, perm, pt
    return 0


def finish(frozen, odir, cache_dir, ext_dir, parts_dir, log, est_s=None, f1_intern=False):
    """Concatenate the per-split parts exactly as run_dataset does, then summarise and confirm_v07.analyze.
    Not resumable inside (a single block); started only if it can finish before the power cut."""
    import confirm_v07 as cf
    import methods_v06 as mv
    import pipeline_v06 as p6
    odir = Path(odir)
    if est_s is not None and not can_start(est_s, "finish (summarise + confirmatory analysis)", log):
        return PAUSED
    t0 = time.time()
    files = _cache_files(cache_dir)
    pf = [Path(parts_dir) / f"stage2_split_{i:02d}.pkl" for i in range(len(files))]
    bad = [p.name for p in pf if not is_done(p)]
    assert not bad, f"stage-2 parts missing or incomplete: {bad}"
    parts = []
    canon = {}
    for p in pf:
        with open(p, "rb") as fh:
            x = pickle.load(fh)
        if f1_intern:                               # F1: shared string objects, identical values
            for key in ("main", "perm"):
                x[key] = [intern_strings(m, canon) for m in x[key]]
            x["pt"] = intern_strings(x["pt"], canon)
            x["covariates"] = intern_strings(x["covariates"], canon)
        parts.append(x)
        del x
    hcfg, dln, sln = frozen["h"], frozen["d"], frozen["s"]
    log(f"selected proxies: h={hcfg} d={dln} s={sln}")
    pt = pd.concat([x["pt"] for x in parts], ignore_index=True)
    pt.to_parquet(odir / "proxy_rows.parquet")
    cov = pd.concat([m for x in parts for m in x["main"]], ignore_index=True)
    log(f"stage 2 rows {len(cov)} (from {len(parts)} per-split parts)")
    hyp = {m: {"chosen": h} for m, h in frozen["hypers"].items()}
    prim, best_dist = frozen["primary"], frozen["best_dist"]
    prim_h = hyp[prim]["chosen"]
    log(f"primary {prim} {prim_h}; best M_dist {best_dist} {hyp[best_dist]['chosen']}")
    perm = [m for x in parts for m in x["perm"]]
    cov = pd.concat([cov] + perm, ignore_index=True)
    cov.to_parquet(odir / "coverage_rows.parquet")
    covs = pd.concat([x["covariates"] for x in parts], ignore_index=True)
    covs.to_parquet(odir / "covariates_analysis_only.parquet")
    mv.M3_FALLBACKS.clear()
    for x in parts:
        for k, v in x["m3_fallbacks"].items():
            mv.M3_FALLBACKS[k] = mv.M3_FALLBACKS.get(k, 0) + v
    p6.summarise(odir, cov, pt, covs, hcfg, None, dln, None, sln, None, hyp, prim, prim_h, best_dist,
                 [x["K"] for x in parts], [x["d"] for x in parts])
    del parts, perm, covs, pt
    caches = [_load_cache_with_ext(f, ext_dir) for f in files]
    out, rm = cf.analyze(caches, frozen, log=log, cov_grid=cov)
    rm.to_parquet(odir / "confirm_rand_rows.parquet")
    cf.write(out, odir / "confirm_v07_result.json")
    mark_done(odir / "confirm_v07_result.json", time.time() - t0)
    log(f"classification: {out['classification']}")
    return 0
