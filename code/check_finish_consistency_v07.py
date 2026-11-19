# check_finish_consistency_v07.py -- after the GWPS finish: confirmatory outputs vs per-split outputs (2026-10-13).
#   python check_finish_consistency_v07.py <dataset>
# 1) coverage_rows.parquet restricted to split i equals split i's stage-2 part (main rows, then M_perm rows), for all i;
# 2) proxy / covariate rows per split equal the parts;
# 3) confirm_rand_rows.parquet per split equals confirm_v07.rand_rows recomputed from that split's cache (splits 0, 10, 20);
# 4) counts: n_test_rows of the result = sum of test sizes over caches = primary matched M0 rows in coverage_rows.
import json
import pickle
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import resume_v07 as rs  # noqa: E402
from confirm_run_v07 import FROZEN, dirs  # noqa: E402


def main(ds):
    import confirm_v07 as cf
    import eval_v06 as ev
    _, _, od = dirs(ds)
    res = {"per_split": {}}
    sel = {"h": tuple(FROZEN["h"].split("|")), "d": FROZEN["d"], "s": FROZEN["s"]}
    grid = [(m, cf._hv(m, h)) for m, h in FROZEN["hypers"].items()]
    n_main_total = sum(json.loads((od / "stage2_parts" / f"stage2_split_{i:02d}.pkl.done.json").read_text())["rows"]
                       for i in range(21))
    # row offsets of main / perm blocks in the frozen concatenation order
    main_lens, perm_lens = [], []
    for i in range(21):
        with open(od / "stage2_parts" / f"stage2_split_{i:02d}.pkl", "rb") as fh:
            x = pickle.load(fh)
        main_lens.append(sum(len(m) for m in x["main"]))
        perm_lens.append(sum(len(m) for m in x["perm"]))
        del x
    assert sum(main_lens) + sum(perm_lens) == n_main_total
    pf = pq.ParquetFile(od / "coverage_rows.parquet")
    assert pf.metadata.num_rows == n_main_total, (pf.metadata.num_rows, n_main_total)
    offs_m = [sum(main_lens[:i]) for i in range(21)]
    offs_p = [sum(main_lens) + sum(perm_lens[:i]) for i in range(21)]
    pt_all = pd.read_parquet(od / "proxy_rows.parquet")
    cv_all = pd.read_parquet(od / "covariates_analysis_only.parquet")
    rr_all = pd.read_parquet(od / "confirm_rand_rows.parquet")
    for i in range(21):
        with open(od / "stage2_parts" / f"stage2_split_{i:02d}.pkl", "rb") as fh:
            x = pickle.load(fh)
        main = pd.concat(x["main"], ignore_index=True)
        perm = pd.concat(x["perm"], ignore_index=True)
        got_m = _rows(od, offs_m[i], main_lens[i])
        got_p = _rows(od, offs_p[i], perm_lens[i])
        r = {"main_equal": _eq(got_m, main), "perm_equal": _eq(got_p, perm),
             "proxy_equal": _eq(pt_all[pt_all.split == x["split_id"]].reset_index(drop=True), x["pt"].reset_index(drop=True)),
             "covariates_equal": _eq(cv_all[cv_all.split == x["split_id"]].reset_index(drop=True),
                                     x["covariates"].reset_index(drop=True))}
        if i in (0, 10, 20):
            c = rs._load_cache_with_ext(od / "cache" / f"split_{i:02d}.pkl", od / "ext_cache")
            rr = cf.rand_rows(c, ev.PRIMARY_MODEL, ev.PRIMARY_K, sel, grid)
            r["rand_rows_equal"] = _eq(rr_all[rr_all.split == i].reset_index(drop=True), rr.reset_index(drop=True))
            del c, rr
        res["per_split"][i] = r
        del x, main, perm, got_m, got_p
    out = json.loads((od / "confirm_v07_result.json").read_text())
    n_test = 0
    for i in range(21):
        with open(od / "cache" / f"split_{i:02d}.pkl", "rb") as fh:
            n_test += len(pickle.load(fh)["test"])
    m0 = pd.read_parquet(od / "coverage_rows.parquet", columns=["model", "k", "variant", "mode", "perm", "method"],
                         filters=[("model", "==", ev.PRIMARY_MODEL), ("k", "==", ev.PRIMARY_K), ("variant", "==", "full"),
                                  ("mode", "==", "matched"), ("perm", "==", -1), ("method", "==", "M0")])
    res["counts"] = {"n_test_rows_result": out["n_test_rows"], "sum_test_sizes": n_test, "grid_matched_M0_rows": len(m0)}
    res["pass"] = all(all(v.values()) for v in res["per_split"].values()) and \
        out["n_test_rows"] == n_test == len(m0)
    (od / "finish_consistency_check.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({"pass": res["pass"], "counts": res["counts"],
                      "failed": {k: v for k, v in res["per_split"].items() if not all(v.values())}}), flush=True)


def _rows(od, start, n):
    pf = pq.ParquetFile(od / "coverage_rows.parquet")
    out, pos = [], 0
    for g in range(pf.num_row_groups):
        m = pf.metadata.row_group(g).num_rows
        if pos + m > start and pos < start + n:
            t = pf.read_row_group(g).to_pandas()
            a, b = max(0, start - pos), min(m, start + n - pos)
            out.append(t.iloc[a:b])
        pos += m
    return pd.concat(out, ignore_index=True)


def _eq(a, b):
    """Exact value equality; dtype differences caused by concatenation across splits (e.g. int -> float upcast)
    are tolerated, values must still be identical (check_exact)."""
    a, b = a.reset_index(drop=True), b.reset_index(drop=True)
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        return False
    if a.equals(b):
        return True
    try:
        pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=True)
        return True
    except AssertionError:
        return False


if __name__ == "__main__":
    main(sys.argv[1])
