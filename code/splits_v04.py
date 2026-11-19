# splits_v04.py -- audit B (structured) split for protocol v0.4.2.
# Replaces audit_lib.cluster_holdout_split (v0.2), which was never executed and fails: it passes
# fractions without a "test" key to random_group_split, and its fractions did not refer to the
# total number of units.
import numpy as np


def cluster_holdout_split(group_id, cluster_id, seed, test_frac=0.15, val_frac=0.10, cal_frac=0.15):
    """Whole prior-defined clusters go to test, added in random order until >= test_frac * n units;
    any TSS group touching test is moved entirely to test. The remaining units are split at random
    by TSS group into val / cal (quotas = fraction x n, n = all units) and train (the rest)."""
    rng = np.random.default_rng(seed)
    group_id = np.asarray(group_id)
    cluster_id = np.asarray(cluster_id)
    n = len(group_id)
    n_test = int(np.floor(test_frac * n))
    test = np.zeros(n, bool)
    for c in rng.permutation(np.unique(cluster_id)):
        if test.sum() >= n_test:
            break
        test[cluster_id == c] = True
    for g in np.unique(group_id[test]):
        test[group_id == g] = True
    rest = np.where(~test)[0]
    quota = {"cal": int(np.floor(cal_frac * n)), "val": int(np.floor(val_frac * n))}
    fold = np.full(len(rest), "train", dtype=object)
    counts = {"cal": 0, "val": 0}
    order = ["cal", "val"]
    k = 0
    for g in rng.permutation(np.unique(group_id[rest])):
        members = np.where(group_id[rest] == g)[0]
        while k < len(order) and counts[order[k]] + len(members) > quota[order[k]]:
            k += 1
        if k == len(order):
            break
        fold[members] = order[k]
        counts[order[k]] += len(members)
    out = {nm: rest[fold == nm] for nm in ("train", "val", "cal")}
    out["test"] = np.where(test)[0]
    return out


def cluster_holdout_split_C(group_id, cluster_id, seed, test_frac=0.15, val_frac=0.10, cal_frac=0.15):
    """Audit B, stopping rule C (protocol v0.5). Clusters are visited in random order; a cluster
    (together with every TSS group it touches) is added to test only if the resulting test size is at
    least as close to test_frac * n as the current size; stop once test >= test_frac * n or all
    clusters were visited. The remaining units are split at random by TSS group into val / cal
    (quotas = fraction x n, n = all units) and train (the rest)."""
    rng = np.random.default_rng(seed)
    group_id = np.asarray(group_id)
    cluster_id = np.asarray(cluster_id)
    n = len(group_id)
    target = test_frac * n
    test = np.zeros(n, bool)
    for c in rng.permutation(np.unique(cluster_id)):
        add = test | (cluster_id == c)
        for g in np.unique(group_id[add & ~test]):
            add |= group_id == g
        if abs(add.sum() - target) <= abs(test.sum() - target):
            test = add
        if test.sum() >= target:
            break
    rest = np.where(~test)[0]
    quota = {"cal": int(np.floor(cal_frac * n)), "val": int(np.floor(val_frac * n))}
    fold = np.full(len(rest), "train", dtype=object)
    counts = {"cal": 0, "val": 0}
    order = ["cal", "val"]
    k = 0
    for g in rng.permutation(np.unique(group_id[rest])):
        members = np.where(group_id[rest] == g)[0]
        while k < len(order) and counts[order[k]] + len(members) > quota[order[k]]:
            k += 1
        if k == len(order):
            break
        fold[members] = order[k]
        counts[order[k]] += len(members)
    out = {nm: rest[fold == nm] for nm in ("train", "val", "cal")}
    out["test"] = np.where(test)[0]
    return out
