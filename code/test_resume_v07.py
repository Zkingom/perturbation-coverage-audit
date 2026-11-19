# test_resume_v07.py -- validation of resume_v07 on the K562 dev set (2026-10-12). Interrupted-then-resumed runs
# must reproduce the one-shot results exactly. Uses dev data only.
#   python test_resume_v07.py <workdir>
# Reference one-shot outputs: results/d1_v0.6.1/cache/split_00.pkl (stage 1, pipeline_v06.stage1),
# results/d1_v0.7dev/ext_cache/ext_split_00.npz (ext), results/d1_v0.7dev/{coverage_rows,proxy_rows,
# covariates_analysis_only}.parquet (run_dataset, 21 splits; rows of splits 0 and 1 compared).
import json
import pickle
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
import resume_v07 as rs  # noqa: E402
from paths import PROCESSED, RESULTS  # noqa: E402

FROZEN = None
RES = {}


def check(name, cond, detail=None):
    RES[name] = {"pass": bool(cond), "detail": detail}
    print(("PASS " if cond else "FAIL ") + name + ("" if detail is None else f"  {detail}"), flush=True)


def child(kind, wd, only):
    import eval_v06 as ev
    import pipeline_v06 as p6
    from confirm_run_v07 import FROZEN as FZ
    wd = Path(wd)
    log = lambda m: print(m, flush=True)  # noqa: E731
    if kind == "stage2":
        return rs.stage2(FZ, wd / "cache", wd / "ext", wd / "parts", log, None, only)
    D = p6.load_dev_v06(PROCESSED / "v0.6" / "k562ess")
    if kind == "stage1":
        return rs.stage1(D, ev.DEV_SEED_BASE, wd / "cache", log, None, only)
    return rs.ext(D, wd / "cache", wd / "ext", log, None, only)


def run_child(kind, wd, only, kill_when=None, timeout=7200):
    cmd = [sys.executable, __file__, "child", kind, str(wd)] + [str(i) for i in only]
    pr = subprocess.Popen(cmd, stdout=open(Path(wd) / f"child_{kind}_{int(time.time())}.txt", "w"),
                          stderr=subprocess.STDOUT)
    t0 = time.time()
    while pr.poll() is None and time.time() - t0 < timeout:
        if kill_when is not None and kill_when(time.time() - t0):
            psutil.Process(pr.pid).kill()          # TerminateProcess: abrupt stop, as in a power cut
            pr.wait()
            return "killed"
        time.sleep(2)
    return pr.returncode


def main(wd):
    wd = Path(wd)
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "cache").mkdir(exist_ok=True)
    (wd / "ext").mkdir(exist_ok=True)
    ref1 = RESULTS / "d1_v0.6.1" / "cache" / "split_00.pkl"
    # T1 stage 1: kill mid-computation, resume; then truncate the finished file (marker kept), resume again
    f0 = wd / "cache" / "split_00.pkl"
    r = run_child("stage1", wd, [0], kill_when=lambda t: t > 60)
    check("T1a stage1 killed mid-split leaves no completed split", r == "killed" and rs.is_done(f0) is None, r)
    r = run_child("stage1", wd, [0])
    check("T1b resumed stage1 split 00 byte-identical to one-shot pipeline_v06.stage1",
          r == 0 and f0.read_bytes() == ref1.read_bytes(), r)
    b = f0.read_bytes()
    f0.write_bytes(b[:len(b) // 3])
    check("T1c truncated file with stale marker is detected as incomplete", rs.is_done(f0) is None)
    r = run_child("stage1", wd, [0])
    check("T1d recomputed split 00 byte-identical again", r == 0 and f0.read_bytes() == ref1.read_bytes(), r)
    # T2 ext: kill, resume, compare arrays bitwise with the one-shot npz (npz zip timestamps differ by design)
    e0 = wd / "ext" / "ext_split_00.npz"
    r = run_child("ext", wd, [0], kill_when=lambda t: t > 45)
    check("T2a ext killed mid-split leaves no completed file", r == "killed" and rs.is_done(e0) is None, r)
    r = run_child("ext", wd, [0])
    za, zb = np.load(e0), np.load(RESULTS / "d1_v0.7dev" / "ext_cache" / "ext_split_00.npz")
    same = sorted(za.files) == sorted(zb.files) and all(
        za[k].dtype == zb[k].dtype and za[k].shape == zb[k].shape and za[k].tobytes() == zb[k].tobytes() for k in za.files)
    check("T2b resumed ext arrays bitwise identical to one-shot ext", r == 0 and same, r)
    # T3 stage 2 for splits 0 and 1 (split 1 inputs copied from the verified one-shot caches)
    man = {x["split"]: x["sha256"] for x in json.loads((RESULTS.parent / "manifests" / "d1_stage1_cache_v0.6.2.json").read_text())["files"]}
    f1 = wd / "cache" / "split_01.pkl"
    shutil.copy2(RESULTS / "d1_v0.6.1" / "cache" / "split_01.pkl", f1)
    assert rs.sha256(f1) == man[1]
    rs.mark_done(f1, 0)
    shutil.copy2(RESULTS / "d1_v0.7dev" / "ext_cache" / "ext_split_01.npz", wd / "ext" / "ext_split_01.npz")
    rs.mark_done(wd / "ext" / "ext_split_01.npz", 0)
    p0 = wd / "parts" / "stage2_split_00.pkl"
    state = {"t": None}

    def kill_in_split1(t):
        if rs.is_done(p0) and state["t"] is None:
            state["t"] = t
        return state["t"] is not None and t > state["t"] + 90
    r = run_child("stage2", wd, [0, 1], kill_when=kill_in_split1)
    check("T3a stage2 killed during split 01 after split 00 completed",
          r == "killed" and rs.is_done(p0) is not None and rs.is_done(wd / "parts" / "stage2_split_01.pkl") is None, r)
    r = run_child("stage2", wd, [0, 1])
    parts = [pickle.load(open(wd / "parts" / f"stage2_split_{i:02d}.pkl", "rb")) for i in (0, 1)]
    cov = pd.concat([m for x in parts for m in x["main"]], ignore_index=True)
    cov = pd.concat([cov] + [m for x in parts for m in x["perm"]], ignore_index=True)
    ref = pd.read_parquet(RESULTS / "d1_v0.7dev" / "coverage_rows.parquet", filters=[("split", "in", [0, 1])])
    ref = ref.reset_index(drop=True)
    cov_rt = pd.read_parquet(_rt(cov, wd / "cov_rt.parquet"))
    check("T3b resumed stage-2 rows (splits 00, 01) identical to the one-shot run_dataset rows",
          r == 0 and len(cov) == len(ref) and list(cov_rt.dtypes) == list(ref.dtypes) and cov_rt.equals(ref),
          {"rows": len(cov), "ref_rows": len(ref)})
    pt = pd.concat([x["pt"] for x in parts], ignore_index=True)
    ptr = pd.read_parquet(RESULTS / "d1_v0.7dev" / "proxy_rows.parquet")
    ptr = ptr[ptr.split.isin([0, 1])].reset_index(drop=True) if "split" in ptr else None
    pt_rt = pd.read_parquet(_rt(pt, wd / "pt_rt.parquet"))
    check("T3c proxy rows identical", ptr is not None and pt_rt.equals(ptr), None if ptr is None else len(ptr))
    cv = pd.concat([x["covariates"] for x in parts], ignore_index=True)
    cvr = pd.read_parquet(RESULTS / "d1_v0.7dev" / "covariates_analysis_only.parquet")
    cvr = cvr[cvr.split.isin([0, 1])].reset_index(drop=True)
    check("T3d covariate rows identical", pd.read_parquet(_rt(cv, wd / "cv_rt.parquet")).equals(cvr))
    (wd / "test_resume_results.json").write_text(json.dumps(RES, indent=1, default=str), encoding="utf8")
    print(f"{sum(v['pass'] for v in RES.values())}/{len(RES)} passed", flush=True)


def _rt(df, path):
    df.to_parquet(path)
    return path


if __name__ == "__main__":
    if sys.argv[1] == "child":
        raise SystemExit(child(sys.argv[2], sys.argv[3], [int(x) for x in sys.argv[4:]]))
    main(sys.argv[1])
