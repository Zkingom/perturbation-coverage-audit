# p0_build.py -- P0 unit table for the Replogle 2022 K562 essential screen (protocol v0.3/v0.4).
# Consolidates the P0 notebook cells of 2026-10-01 into one script; the algorithm and all
# parameters are unchanged. Gene coordinates are NOT re-queried from Ensembl (release 116 was
# used originally); they are read from the archived unit table so that later Ensembl releases
# cannot change the TSS grouping.
import pickle
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import PROCESSED, raw  # noqa: E402
from audit_lib import tss_groups  # noqa: E402

UNITS_ARCHIVE = PROCESSED / "v0.3_figshare_k562ess" / "p0v3_units_all_targets.csv"
MIN_CELLS = 30
TSS_MERGE_BP = 1000


def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def read_bulk_obs_var(path):
    """Read perturbation metadata and readout gene names from the Figshare raw-bulk h5ad."""
    with h5py.File(path, "r") as fb:
        ob = pd.DataFrame({
            "gene_transcript": _dec(fb["obs/gene_transcript"][:]),
            "ncells": fb["obs/num_cells_filtered"][:],
            "core_control": fb["obs/core_control"][:],
            "energy_p": fb["obs/energy_test_p_value"][:],
            "fold_expr": fb["obs/fold_expr"][:],
        })
        vcat = _dec(fb["var/__categories/gene_name"][:])
        readout = vcat[fb["var/gene_name"][:]]
    return ob, readout


def build_units(coords_source=UNITS_ARCHIVE):
    ob, readout = read_bulk_obs_var(raw("K562_essential_raw_bulk_01.h5ad"))
    with open(raw("gene2go_all.pkl"), "rb") as fh:
        gene2go = pickle.load(fh)
    with open(raw("essential_all_data_pert_genes.pkl"), "rb") as fh:
        ess = pickle.load(fh)
    pert_names = set(g for g in ess if g in gene2go)          # GEARS 0.1.2 set_pert_genes()

    parts = ob.gene_transcript.str.split("_", n=3, expand=True)
    ob["gid"], ob["target"], ob["tss"], ob["ensg"] = parts[0], parts[1], parts[2], parts[3]
    ob["is_ctrl"] = ob.gene_transcript.str.contains("non-targeting", case=False)
    P = ob[~ob.is_ctrl].copy()
    # one row per target gene: the TSS row with most cells; other rows never enter any fold
    U = P.sort_values("ncells", ascending=False).drop_duplicates("target").copy()
    R = set(readout)
    U["in_R"] = U.target.isin(R)
    U["in_GO"] = U.target.isin(pert_names)
    U["cells_ok"] = U.ncells >= MIN_CELLS
    U["E"] = U.in_R & U.in_GO & U.cells_ok
    U["E_noGO"] = U.in_R & ~U.in_GO & U.cells_ok
    U["sig"] = U.energy_p < 0.05

    # coordinates (Ensembl release 116, GRCh38.p14) from the archived table, keyed by ENSG
    c = pd.read_csv(coords_source, usecols=["ensg", "chr", "strand", "start", "end"]).drop_duplicates("ensg")
    c = c.set_index("ensg")
    U["chr"] = U.ensg.map(c["chr"])
    U["strand"] = U.ensg.map(c["strand"])
    U["start"] = U.ensg.map(c["start"])
    U["end"] = U.ensg.map(c["end"])
    assert U.chr.notna().all(), "coordinates missing for some targets"
    U = U.reset_index(drop=True)
    U["tss_group"] = tss_groups(U.target.values, U.chr.values, U.strand.values,
                                U.start.values.astype(np.int64), U.end.values.astype(np.int64), TSS_MERGE_BP)
    grp = U.groupby("tss_group").agg(nE=("E", "sum"), n=("E", "size"))
    mixed = grp[(grp.nE > 0) & (grp.n > grp.nE)].index
    U["always_train"] = ~U.E & ~U.tss_group.isin(mixed)
    U["dropped_tss_conflict"] = ~U.E & U.tss_group.isin(mixed)
    assert U.target.is_unique
    return U, readout, n_secondary_rows(P, U)


def n_secondary_rows(P, U):
    return len(P) - len(U)


if __name__ == "__main__":
    U, readout, n_sec = build_units()
    out = PROCESSED / "v0.3_figshare_k562ess" / "p0_units_rebuilt.csv"
    U.to_csv(out, index=False)
    print({"E": int(U.E.sum()), "tss_groups_E": int(U[U.E].tss_group.nunique()),
           "always_train": int(U.always_train.sum()), "dropped": int(U.dropped_tss_conflict.sum()),
           "secondary_rows": n_sec, "readout_genes": len(readout)})
