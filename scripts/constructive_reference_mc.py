# constructive_reference_mc.py -- Monte Carlo check of the constructive T5 reference (Methods, "Constructive reference").
# Setting: the conformity score equals the stratifying variable E (a strictly increasing function suffices), E ~ log-normal
# (meanlog 1.2, sdlog 0.35) for the n_units of a dataset; per replicate, a calibration set of n_cal and a test set of n_test
# units are drawn without replacement from those units; split-conformal threshold q = ceil((n_cal+1)(1-alpha))-th smallest
# calibration score; T5 = test units above the 80th percentile of E over all units; statistic = fraction of T5 test units
# with score <= q. 4,000 replicates per dataset, one generator with seed 20261001, datasets in the order dev, rpe1, gwps.
# Sizes are those of one confirmation split (reports/review/constructive_reference_simulation.json).
#   python scripts/constructive_reference_mc.py   -> results/constructive_reference_mc/mc_result.json
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SIZES = {"dev": (267, 255, 1786), "rpe1": (262, 251, 1751), "gwps": (852, 812, 5683)}
SEED, REPS, ALPHA = 20261001, 4000, 0.1


def sim(rng, n_cal, n_test, n_all, reps=REPS, alpha=ALPHA):
    out = np.empty(reps)
    for i in range(reps):
        z = rng.lognormal(1.2, 0.35, n_all)
        cal = rng.choice(z, n_cal, replace=False)
        te = rng.choice(z, n_test, replace=False)
        q = np.sort(cal)[int(np.ceil((n_cal + 1) * (1 - alpha))) - 1]
        t5 = te[te > np.quantile(z, 0.8)]
        out[i] = (t5 <= q).mean() if len(t5) else np.nan
    return out


def main():
    rng = np.random.default_rng(SEED)
    res = {}
    for ds, (nc, nt, na) in SIZES.items():
        x = sim(rng, nc, nt, na)
        res[ds] = {"n_cal": nc, "n_test_per_split": nt, "n_units": na, "sim_mean": float(np.nanmean(x)),
                   "mc_se": float(np.nanstd(x, ddof=1) / np.sqrt(np.isfinite(x).sum())),
                   "sim_ci_single_split": [float(v) for v in np.nanpercentile(x, [2.5, 97.5])]}
    out = ROOT / "results" / "constructive_reference_mc"
    out.mkdir(parents=True, exist_ok=True)
    (out / "mc_result.json").write_text(json.dumps({"seed": SEED, "reps": REPS, "alpha": ALPHA, "results": res}, indent=1))
    print(json.dumps(res))


if __name__ == "__main__":
    main()
