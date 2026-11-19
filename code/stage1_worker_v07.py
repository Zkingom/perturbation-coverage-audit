# stage1_worker_v07.py -- execution helper (2026-10-10, user-approved speed-up; no frozen file changed).
# Computes stage-1 caches for a given list of split ids with the FROZEN eval_v06.compute_split, writing
# <cache>/split_NN.pkl through a worker-specific temp file + atomic rename. Splits that already exist are skipped,
# so it can run alongside confirm_run_v07.py stage1 (which also skips existing files).
#   python stage1_worker_v07.py <dev|rpe1> <cache_dir> <id> [<id> ...]
import os
import pickle
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import eval_v06 as ev  # noqa: E402
from splits_v06 import repeated_group_kfold  # noqa: E402


def main():
    ds, cache, ids = sys.argv[1], Path(sys.argv[2]), [int(x) for x in sys.argv[3:]]
    if ds == "dev":
        import pipeline_v06 as p6
        from paths import PROCESSED
        D = p6.load_dev_v06(PROCESSED / "v0.6" / "k562ess")
        base = ev.DEV_SEED_BASE
    else:
        import confirm_run_v07 as cr
        cr.check_frozen()
        D = cr.load(ds)
        cr.check_splits(ds, D)
        base = ev.CONF_SEED_BASE
    sp = repeated_group_kfold(D["U"].tss_group.values[D["E_idx"]], base, n_rounds=3)
    cache.mkdir(parents=True, exist_ok=True)
    for i in ids:
        f = cache / f"split_{i:02d}.pkl"
        if f.exists():
            continue
        c = ev.compute_split(i, sp[i], D)
        tmp = cache / f"split_{i:02d}.w{os.getpid()}.tmp"
        with open(tmp, "wb") as fh:
            pickle.dump(c, fh)
        if not f.exists():
            os.replace(tmp, f)
        else:
            tmp.unlink()
        print(f"worker {os.getpid()} split {i:02d} done", flush=True)


if __name__ == "__main__":
    main()
