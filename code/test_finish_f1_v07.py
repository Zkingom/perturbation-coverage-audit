# test_finish_f1_v07.py -- validation of the F1 finish (string interning) on the K562 dev set (2026-10-13).
# Builds dev stage-2 parts with the validated resumable stage2 (frozen v0.7 selection), then runs the finish twice
# in separate processes: A = frozen finish logic, B = F1 (intern_strings). Acceptance (user): summary.json and
# confirm_v07_result.json byte-identical; also coverage/proxy/covariate/rand-row tables equal on read-back.
#   python test_finish_f1_v07.py prep <workdir>          (dev stage-2 parts; resumable)
#   python test_finish_f1_v07.py child <A|B> <workdir>
#   python test_finish_f1_v07.py compare <workdir>
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
import resume_v07 as rs  # noqa: E402
from paths import MANIFESTS, RESULTS  # noqa: E402


def prep(wd):
    from confirm_run_v07 import FROZEN
    wd = Path(wd)
    (wd / "cache").mkdir(parents=True, exist_ok=True)
    (wd / "ext").mkdir(exist_ok=True)
    man = {x["split"]: x["sha256"] for x in json.loads((MANIFESTS / "d1_stage1_cache_v0.6.2.json").read_text())["files"]}
    for i in range(21):
        f = wd / "cache" / f"split_{i:02d}.pkl"
        if not rs.is_done(f):
            shutil.copy2(RESULTS / "d1_v0.6.1" / "cache" / f.name, f)
            assert rs.sha256(f) == man[i], f"dev cache {i} hash differs from manifest"
            rs.mark_done(f, 0, source="results/d1_v0.6.1/cache (manifest-verified)")
        e = wd / "ext" / f"ext_split_{i:02d}.npz"
        if not rs.is_done(e):
            shutil.copy2(RESULTS / "d1_v0.7dev" / "ext_cache" / e.name, e)
            rs.mark_done(e, 0, source="results/d1_v0.7dev/ext_cache")
    log = lambda m: print(m, flush=True)  # noqa: E731
    return rs.stage2(FROZEN, wd / "cache", wd / "ext", wd / "parts", log, None, None)


def child(which, wd):
    from confirm_run_v07 import FROZEN
    wd = Path(wd)
    od = wd / f"finish_{which}"
    od.mkdir(exist_ok=True)
    log = lambda m: print(m, flush=True)  # noqa: E731
    t0 = time.time()
    rc = rs.finish(FROZEN, od, wd / "cache", wd / "ext", wd / "parts", log, None, f1_intern=(which == "B"))
    mi = psutil.Process().memory_info()
    (od / "mem.json").write_text(json.dumps({"rc": rc, "seconds": time.time() - t0, "peak_wset_gb": mi.peak_wset / 1024 ** 3,
                                             "peak_pagefile_gb": mi.peak_pagefile / 1024 ** 3}))
    return rc


def compare(wd):
    wd = Path(wd)
    A, B = wd / "finish_A", wd / "finish_B"
    res = {}
    for f in ("summary.json", "confirm_v07_result.json"):
        res[f + " byte-identical"] = (A / f).read_bytes() == (B / f).read_bytes()
    for f in ("coverage_rows.parquet", "proxy_rows.parquet", "covariates_analysis_only.parquet", "confirm_rand_rows.parquet"):
        a, b = pd.read_parquet(A / f), pd.read_parquet(B / f)
        res[f + " equal"] = bool(a.equals(b) and list(a.dtypes) == list(b.dtypes))
        del a, b
    res["memory"] = {w: json.loads((wd / f"finish_{w}" / "mem.json").read_text()) for w in "AB"}
    res["pass"] = all(v for k, v in res.items() if k != "memory")
    (wd / "test_finish_f1_results.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res), flush=True)


if __name__ == "__main__":
    m, wd = sys.argv[1], sys.argv[-1]
    if m == "prep":
        raise SystemExit(prep(wd))
    if m == "child":
        raise SystemExit(child(sys.argv[2], wd))
    if m == "main":
        for w in "AB":
            r = subprocess.run([sys.executable, __file__, "child", w, wd])
            assert r.returncode == 0, f"finish {w} failed"
        compare(wd)
    else:
        compare(wd)
