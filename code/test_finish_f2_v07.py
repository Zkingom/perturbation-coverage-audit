# test_finish_f2_v07.py -- validation of F2 (finish_f2_v07) on the K562 dev set (2026-10-14). Uses the dev stage-2
# parts and the frozen-finish outputs (finish_A) already produced by test_finish_f1_v07.py in the same workdir.
#   python test_finish_f2_v07.py child <workdir>   ;   python test_finish_f2_v07.py main <workdir>
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))


def child(wd, variant="f2"):
    import finish_f2_v07 as f2
    from confirm_run_v07 import FROZEN
    wd = Path(wd)
    od = wd / ("finish_C" if variant == "f2" else "finish_E")
    od.mkdir(exist_ok=True)
    t0 = time.time()
    fn = f2.finish_f2 if variant == "f2" else f2.finish_f2c
    rc = fn(FROZEN, od, wd / "cache", wd / "ext", wd / "parts", lambda m: print(m, flush=True), None,
                      colstore_dir=od / "_f2_colstore")
    mi = psutil.Process().memory_info()
    (od / "mem.json").write_text(json.dumps({"rc": rc, "seconds": time.time() - t0, "peak_wset_gb": mi.peak_wset / 1024 ** 3}))
    return rc


def compare(wd, variant="f2"):
    wd = Path(wd)
    A, B = wd / "finish_A", wd / ("finish_C" if variant == "f2" else "finish_E")
    res = {}
    for f in ("summary.json", "confirm_v07_result.json"):
        res[f + " byte-identical"] = (A / f).read_bytes() == (B / f).read_bytes()
    for f in ("coverage_rows.parquet", "proxy_rows.parquet", "covariates_analysis_only.parquet", "confirm_rand_rows.parquet"):
        a, b = pd.read_parquet(A / f), pd.read_parquet(B / f)
        res[f + " equal"] = bool(a.equals(b) and list(a.dtypes) == list(b.dtypes) and list(a.columns) == list(b.columns))
        del a, b
    res["memory"] = {"A_frozen": json.loads((A / "mem.json").read_text()), variant: json.loads((B / "mem.json").read_text())}
    res["pass"] = all(v for k, v in res.items() if k != "memory")
    (wd / ("test_finish_f2_results.json" if variant == "f2" else "test_finish_f2c_results.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res), flush=True)


if __name__ == "__main__":
    # python test_finish_f2_v07.py main <workdir> [f2|f2c] ; child <workdir> [f2|f2c]
    var = sys.argv[3] if len(sys.argv) > 3 else "f2"
    if sys.argv[1] == "child":
        raise SystemExit(child(sys.argv[2], var))
    r = subprocess.run([sys.executable, __file__, "child", sys.argv[2], var])
    assert r.returncode == 0, "finish failed"
    compare(sys.argv[2], var)
