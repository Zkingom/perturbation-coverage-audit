# confirm_primary_v07.py -- execution-order helper (2026-10-10, user-requested before a power cut; no frozen file
# changed). Runs ONLY the confirmatory analysis of protocol v0.7 section 4 (frozen confirm_v07.analyze, frozen
# selection, randomised width matching, M_perm) on the existing C1 stage-1 caches + ext arrays, before the full
# descriptive run_dataset (confirm_run_v07.py analyze). cov_grid=None: the grid-matching sensitivity block is left
# to the full run, which recomputes the primary blocks; both must agree exactly.
#   python confirm_primary_v07.py rpe1
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import confirm_run_v07 as cr  # noqa: E402
import confirm_v07 as cf  # noqa: E402
import ext_v07dev as ex  # noqa: E402
import pipeline_v06 as p6  # noqa: E402


def main(ds):
    cr.check_frozen()
    _, _, odir = cr.dirs(ds)
    log = cr.logger(odir / "c1_primary.log")
    import resume_v07 as rs                       # 2026-10-15 (GWPS): watchdog + keep-awake, as in the 'run' mode
    rs.MEM_CAP_GB = 4.8
    log("primary helper start")
    rs.start_guards(log)
    D = cr.load_lowmem(ds)                        # identical values (lowmem_v07, validated); only for check_splits
    cr.check_splits(ds, D)
    del D
    log("split files verified")
    caches = p6.load_caches(odir / "cache")
    log(f"loaded {len(caches)} caches")
    assert len(caches) == 21, f"expected 21 caches, found {len(caches)}"
    for c in caches:
        z = np.load(odir / "ext_cache" / f"ext_split_{c['split_id']:02d}.npz")
        ex.attach(c, {k: z[k] for k in z.files})
    out, rm = cf.analyze(caches, cr.FROZEN, log=log, cov_grid=None)
    rm.to_parquet(odir / "confirm_rand_rows_primary.parquet")
    cf.write(out, odir / "confirm_v07_result_primary.json")
    log(f"classification: {out['classification']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
