# confounding_v07_run2.py -- low-memory variant of confounding_v07_run.py (execution helper, 2026-10-15; no frozen
# file changed). sens_v07.confounding_v07 uses only the coverage rows with k = PRIMARY_K, mode = nominal, method = M0,
# model in LINEAR, perm = -1 (variant full / half). Instead of loading the whole coverage table, those rows are read
# with a pyarrow row filter (row order preserved). The confounding code itself is called unchanged.
#   python confounding_v07_run2.py <results_dir> [--out <json>]
# Validation: on RPE1 the output must be byte-identical to confounding_v07.json written by confounding_v07_run.py.
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import eval_v06 as ev  # noqa: E402
import sens_v07 as sv  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    d = Path(a.results_dir)
    flt = [("k", "==", ev.PRIMARY_K), ("mode", "==", "nominal"), ("method", "==", "M0"), ("perm", "==", -1),
           ("model", "in", list(ev.LINEAR)), ("variant", "in", ["full", "half"])]
    cov = pd.read_parquet(d / "coverage_rows.parquet", filters=flt)
    covs = pd.read_parquet(d / "covariates_analysis_only.parquet")
    print(f"rows read: {len(cov)}", flush=True)
    out = {"full": sv.confounding_v07(cov, covs), "half": sv.confounding_v07(cov, covs, half=True)}
    p = Path(a.out) if a.out else d / "confounding_v07.json"
    p.write_text(json.dumps(out, indent=1, default=sv._js), encoding="utf8")
    print("confounding done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
