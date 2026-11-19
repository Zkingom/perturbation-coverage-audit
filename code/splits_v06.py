# splits_v06.py -- protocol v0.6 section 5: repeated grouped 7-fold.
# Round r (seed = base_seed + r): TSS groups are shuffled and dealt to 7 test folds, each new group going
# to the fold with the fewest units so far (balanced sizes; a group never spans folds). For test fold f
# the remaining units are split at random by TSS group into val (floor(0.10 n)), cal (floor(0.15 n))
# and train (rest); n = number of E units. 3 rounds x 7 folds = 21 splits; every unit is a test unit
# exactly once per round.
import numpy as np

N_FOLDS = 7
N_ROUNDS = 3
VAL_FRAC, CAL_FRAC = 0.10, 0.15


def _fill(group_id, idx, quotas, rng):
    """Assign whole groups of the units `idx` to folds in `quotas` order (dict name -> size);
    the remainder is 'train'. Same greedy rule as audit_lib.random_group_split."""
    names = list(quotas)
    fold = np.full(len(idx), "train", dtype=object)
    counts = dict.fromkeys(names, 0)
    k = 0
    for g in rng.permutation(np.unique(group_id[idx])):
        m = np.where(group_id[idx] == g)[0]
        while k < len(names) and counts[names[k]] + len(m) > quotas[names[k]]:
            k += 1
        if k == len(names):
            break
        fold[m] = names[k]
        counts[names[k]] += len(m)
    return {nm: idx[fold == nm] for nm in names + ["train"]}


def repeated_group_kfold(group_id, base_seed, n_rounds=N_ROUNDS, n_folds=N_FOLDS):
    """Returns a list of dicts {round, fold, seed, train, val, cal, test} with indices into group_id."""
    group_id = np.asarray(group_id)
    n = len(group_id)
    quotas = {"cal": int(np.floor(CAL_FRAC * n)), "val": int(np.floor(VAL_FRAC * n))}
    out = []
    for r in range(n_rounds):
        seed = base_seed + r
        rng = np.random.default_rng(seed)
        groups = rng.permutation(np.unique(group_id))
        size = np.zeros(n_folds, int)
        gfold = {}
        for g in groups:
            f = int(np.argmin(size))
            gfold[g] = f
            size[f] += int((group_id == g).sum())
        unit_fold = np.array([gfold[g] for g in group_id])
        for f in range(n_folds):
            test = np.where(unit_fold == f)[0]
            rest = np.where(unit_fold != f)[0]
            sub = _fill(group_id, rest, quotas, np.random.default_rng([seed, f]))
            out.append({"round": r, "fold": f, "seed": seed, "test": test, **sub})
    return out
