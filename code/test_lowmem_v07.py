# test_lowmem_v07.py -- validation of lowmem_v07 against the frozen loader on K562 dev split 00 (2026-10-12).
#   python test_lowmem_v07.py child <frozen|lowmem> <outdir>    (one process per implementation, peak memory)
#   python test_lowmem_v07.py main <outdir>
import json
import pickle
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PROCESSED, RESULTS  # noqa: E402

MODES = ("frozen", "lowmem2")     # round 2 (2026-10-12): frozen loader vs lowmem loader + install()


def child(mode, out):
    import data_v07 as dv
    import eval_v06 as ev
    import ext_v07dev as ex
    import lowmem_v07 as lm
    from splits_v06 import repeated_group_kfold
    out = Path(out)
    pr = psutil.Process()
    U = pd.read_csv(PROCESSED / "v0.3_figshare_k562ess" / "p0v3_units_all_targets.csv")
    args = (U, PROCESSED / "v0.6" / "k562ess", PROCESSED / "priors" / "string_expdb_v0.4_universe.parquet")
    t0 = time.time()
    D = dv.load_dataset(*args) if mode == "frozen" else lm.load_dataset_lowmem(*args)
    if mode == "lowmem2":
        lm.install(out / "_pred_memmap")
    after_load = pr.memory_info().rss / 1024 ** 3
    sp = repeated_group_kfold(D["U"].tss_group.values[D["E_idx"]], ev.DEV_SEED_BASE)[0]
    t1 = time.time()
    c = ev.compute_split(0, sp, D)
    with open(out / f"split_00_{mode}.pkl", "wb") as fh:
        pickle.dump(c, fh)
    t2 = time.time()
    e = ex.compute_ext(c, D, log=lambda *_: None)
    np.savez(out / f"ext_split_00_{mode}.npz", **e)
    t3 = time.time()
    mi = pr.memory_info()
    rec = {"mode": mode, "rss_after_load_gb": after_load, "peak_wset_gb": mi.peak_wset / 1024 ** 3,
           "peak_pagefile_gb": mi.peak_pagefile / 1024 ** 3, "load_s": t1 - t0, "stage1_s": t2 - t1, "ext_s": t3 - t2}
    (out / f"mem_{mode}.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec), flush=True)


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    for mode in MODES:
        r = subprocess.run([sys.executable, __file__, "child", mode, str(out)], capture_output=True, text=True)
        (out / f"child_{mode}.txt").write_text(r.stdout + r.stderr)
        assert r.returncode == 0, f"{mode} child failed: {r.stderr[-2000:]}"
    ref = RESULTS / "d1_v0.6.1" / "cache" / "split_00.pkl"
    for mode in MODES:
        res[f"stage1_{mode}_bytes_equal_reference"] = (out / f"split_00_{mode}.pkl").read_bytes() == ref.read_bytes()
    zr = np.load(RESULTS / "d1_v0.7dev" / "ext_cache" / "ext_split_00.npz")
    for mode in MODES:
        z = np.load(out / f"ext_split_00_{mode}.npz")
        res[f"ext_{mode}_arrays_bitwise_equal_reference"] = sorted(z.files) == sorted(zr.files) and all(
            z[k].dtype == zr[k].dtype and z[k].shape == zr[k].shape and z[k].tobytes() == zr[k].tobytes() for k in z.files)
    res["memory"] = {m: json.loads((out / f"mem_{m}.json").read_text()) for m in MODES}
    res["pass"] = all(v for k, v in res.items() if k not in ("memory", "pass"))
    (out / "test_lowmem_results.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "child":
        child(sys.argv[2], sys.argv[3])
    else:
        main(sys.argv[2])
