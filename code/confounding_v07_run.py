# confounding_v07_run.py -- runs the frozen v0.7 section-5 confounding analysis (sens_v07.confounding_v07:
# logistic coefficients per SD, TSS-group cluster bootstrap B=500 seed master_seed+9000, Firth sensitivity and a
# Firth bootstrap when > 5% of plain fits fail/separate) on a dataset's saved coverage rows. Execution helper,
# 2026-10-11; no frozen file changed.   python confounding_v07_run.py <results_dir>
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sens_v07 as sv  # noqa: E402


def main(d):
    d = Path(d)
    cov = pd.read_parquet(d / "coverage_rows.parquet")
    covs = pd.read_parquet(d / "covariates_analysis_only.parquet")
    out = {"full": sv.confounding_v07(cov, covs), "half": sv.confounding_v07(cov, covs, half=True)}
    (d / "confounding_v07.json").write_text(json.dumps(out, indent=1, default=sv._js), encoding="utf8")
    print("confounding done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
