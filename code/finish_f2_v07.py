# finish_f2_v07.py -- F2 low-memory finish for the v0.7 confirmation runs (engineering only, 2026-10-14; user rule:
# F1 exceeded 4.8 GB on GWPS -> F2, acceptance = byte-identical summary.json and confirm_v07_result.json on dev).
#
# What changes: only how the stage-2 coverage table is held in memory.
#   * The table is never materialised as one DataFrame. Each column is concatenated on its own through the SAME
#     pandas DataFrame-concat path as the frozen finish (first the main frames, then [that] + the M_perm frames),
#     so every column gets the identical values and dtype. Numeric / bool columns are stored as .npy files on F:
#     (memory-mapped, float64 / int64 / bool unchanged); string columns stay in RAM as object arrays with shared
#     string objects (resume_v07.intern_strings, validated as F1).
#   * Every subset the frozen pipeline_v06.summarise / confirm_v07.analyze takes from the full table
#     (cov[mask], cov[cov.perm == -1].groupby(...)) is rebuilt from the column store with the same boolean mask,
#     evaluated on the same key columns, and with the same index labels, column order and dtypes. All analysis
#     functions (analysis_v06.*, confirm_v07.*) are the frozen ones and are called in the same order with the same
#     arguments; summary.json is written with the same json.dumps call.
#   * coverage_rows.parquet is written in row batches with an explicit schema (reads back to the same table).
# No algorithm, parameter or random-number change.
import json
import os
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

import resume_v07 as rs

KEYS = ["model", "k", "variant", "mode", "method", "hyper", "perm"]
BATCH = 500_000


class ColumnStore:
    """Column-wise copy of the frozen `cov` table (identical values, dtypes, column order, RangeIndex labels)."""

    def __init__(self, cols, arrays, n):
        self.cols, self.a, self.n = cols, arrays, n

    def series(self, c):
        return pd.Series(np.asarray(self.a[c]), name=c)

    def take(self, pos):
        pos = np.asarray(pos, dtype=np.int64)
        data = {c: np.asarray(self.a[c][pos]) for c in self.cols}
        return pd.DataFrame(data, columns=self.cols, index=pd.Index(pos, dtype="int64"))

    def keyframe(self):
        return pd.DataFrame({c: np.asarray(self.a[c]) for c in KEYS}, columns=KEYS)


def build_store(parts_files, colstore_dir, log):
    """Load the per-split parts once (strings interned), concatenate column by column exactly as the frozen
    finish does, return (store, pt, covs, fallbacks, K, d)."""
    colstore_dir = Path(colstore_dir)
    colstore_dir.mkdir(parents=True, exist_ok=True)
    mains, perms, pts, covss, fb, Ks, ds = [], [], [], [], [], [], []
    canon = {}
    for p in parts_files:
        with open(p, "rb") as fh:
            x = pickle.load(fh)
        mains += [rs.intern_strings(m, canon) for m in x["main"]]
        perms += [rs.intern_strings(m, canon) for m in x["perm"]]
        pts.append(x["pt"])
        covss.append(x["covariates"])
        fb.append(x["m3_fallbacks"])
        Ks.append(x["K"])
        ds.append(x["d"])
        del x
    cols = list(mains[0].columns)
    for m in mains + perms:
        assert list(m.columns) == cols, "stage-2 frames differ in columns; F2 column-wise concat not applicable"
    n_main = sum(len(m) for m in mains)
    log(f"F2: {len(mains)} main + {len(perms)} perm frames, {n_main} main rows, {len(cols)} columns")
    arrays = {}
    for c in cols:
        s1 = pd.concat([m[[c]] for m in mains], ignore_index=True)          # = frozen first concat, column c
        s = pd.concat([s1] + [m[[c]] for m in perms], ignore_index=True)[c]  # = frozen second concat, column c
        del s1
        v = s.to_numpy()
        if v.dtype == object:
            nn = pd.notna(v)
            assert all(type(t) is str for t in v[nn]), f"column {c}: non-str objects"
            arrays[c] = rs.intern_strings(pd.DataFrame({c: v}), canon)[c].to_numpy()
        else:
            f = colstore_dir / f"{c}.npy"
            np.save(f, v)
            arrays[c] = np.load(f, mmap_mode="r")
        del s, v
    n = len(arrays[cols[0]])
    del mains, perms
    pt = pd.concat(pts, ignore_index=True)
    covs = pd.concat(covss, ignore_index=True)
    return ColumnStore(cols, arrays, n), pt, covs, fb, Ks, ds


def write_parquet(store, path):
    fields = []
    for c in store.cols:
        a = store.a[c]
        if a.dtype == object:
            fields.append(pa.field(c, pa.string()))
        else:
            fields.append(pa.field(c, pa.from_numpy_dtype(a.dtype)))
    schema = pa.schema(fields)
    tmp = Path(str(path) + ".tmp")
    with pq.ParquetWriter(tmp, schema) as w:
        for a in range(0, store.n, BATCH):
            b = min(a + BATCH, store.n)
            df = pd.DataFrame({c: np.asarray(store.a[c][a:b]) for c in store.cols}, columns=store.cols)
            w.write_table(pa.Table.from_pandas(df, schema=schema, preserve_index=False))
    with open(tmp, "rb+") as fh:
        os.fsync(fh.fileno())
    rs._commit(tmp, Path(path))


def summarise_f2(out, store, pt, covs, hcfg, rho_h, dln, rho_d, sln, rho_s, hyp, prim, prim_h, best_dist,
                 chosen_K, chosen_d):
    """pipeline_v06.summarise with every full-table subset taken from the column store (same masks, same order)."""
    import pipeline_v06 as p6
    an, ev = p6.an, p6.ev
    out = Path(out)
    P = ev.PRIMARY_MODEL
    Kf = store.keyframe()
    pos = lambda m: np.flatnonzero(m.to_numpy())                                      # noqa: E731
    pm = lambda mode, model=P: store.take(pos((Kf.model == model) & (Kf.k == ev.PRIMARY_K)  # noqa: E731
                                              & (Kf.variant == "full") & (Kf["mode"] == mode) & (Kf.perm == -1)))
    chosen = lambda d, meth: d[(d.method == meth) & (d.hyper == hyp[meth]["chosen"])]        # noqa: E731
    summ = {"selection": {"h": hcfg, "rho_h": rho_h, "d": dln, "rho_d": rho_d, "s": sln, "rho_s": rho_s,
                          "hypers": {m: v["chosen"] for m, v in hyp.items()}, "hyper_detail": hyp,
                          "primary": [prim, prim_h], "best_dist": [best_dist, hyp[best_dist]["chosen"]]},
            "chosen_K": chosen_K, "chosen_d": chosen_d}
    tabs = {}
    for (m, k, var, mode, meth, hy), sub in Kf[Kf.perm == -1].groupby(["model", "k", "variant", "mode", "method", "hyper"]):
        tabs[f"{m}|{k}|{var}|{mode}|{meth}|{hy}"] = an.tier_table(store.take(sub.index.to_numpy()))
    summ["tier_tables"] = tabs
    summ["monotone_all"] = bool(store.series("monotone").all())
    import methods_v06 as mv
    summ["m3_solver_fallbacks"] = dict(mv.M3_FALLBACKS)
    m_rows = pm("matched")
    m0n = pm("nominal")
    pmc = {"nominal": m0n, "matched": m_rows}
    summ["dgo_strata"] = {f"{mode}|{meth}": an.strata_table(chosen(pmc[mode], meth) if meth != "M0"
                                                              else pmc[mode][pmc[mode].method == "M0"], "D_GO")
                          for mode in ("nominal", "matched")
                          for meth in ["M0", "M1", prim, "Mscore-N"] + [b for b in an.BASELINES if b in hyp]}
    d0m = m_rows[m_rows.method == "M0"].sort_values(["split", "target"])
    comp = {}
    for meth in [prim, best_dist, "Mcomb-N", "Mcomb-3", "Mscore-N", "oracle"]:
        dm = chosen(m_rows, meth).sort_values(["split", "target"])
        comp[f"{meth}_vs_M0"] = an.paired_T5(dm, d0m)
    dbd = chosen(m_rows, best_dist).sort_values(["split", "target"])
    for meth in [prim, "Mcomb-N", "Mcomb-3", "Mscore-N"]:
        comp[f"{meth}_vs_{best_dist}"] = an.paired_T5(chosen(m_rows, meth).sort_values(["split", "target"]), dbd)
    bases = [b for b in an.BASELINES if b in hyp]
    for b in bases:
        db = chosen(m_rows, b).sort_values(["split", "target"])
        comp[f"{b}_vs_M0"] = an.paired_T5(db, d0m)
        for meth in sorted({prim, "Mscore-N"}):
            if meth != b:
                comp[f"{meth}_vs_{b}"] = an.paired_T5(chosen(m_rows, meth).sort_values(["split", "target"]), db)
    summ["matched_paired_T5"] = comp
    dec = {}
    for meth in [prim, "Mscore-N"] + bases:
        dec[f"{meth}_vs_M0"] = an.decision_paired(chosen(m_rows, meth), d0m)
    for meth in sorted({prim, "Mscore-N"}):
        dec[f"{meth}_vs_{best_dist}"] = an.decision_paired(chosen(m_rows, meth), chosen(m_rows, best_dist))
    summ["decision_paired_top10_matched"] = dec
    rows_lin = m0n.iloc[:0]
    for m in ev.LINEAR:
        x = pm("nominal", m)
        rows_lin = pd.concat([rows_lin, x[x.method == "M0"]])
        del x
    summ["confounding"] = an.confounding(an.confounding_frame(rows_lin, covs))
    hr = store.take(pos((Kf.k == ev.PRIMARY_K) & (Kf.variant == "half") & (Kf["mode"] == "nominal")
                        & (Kf.method == "M0") & Kf.model.isin(ev.LINEAR) & (Kf.perm == -1)))
    summ["confounding_half"] = an.confounding(an.confounding_frame(hr, covs, half=True))
    del hr
    summ["incremental"] = an.incremental(rows_lin, covs)
    pth = pt[(pt.kind == "h") & (pt.proxy == hcfg)].rename(columns={"value": "h"})
    for mode in ("nominal", "matched"):
        cm = store.take(pos((Kf.variant == "full") & (Kf["mode"] == mode)))
        summ[f"criteria_{mode}"] = an.criteria(cm, pth, prim, prim_h, P, ev.PRIMARY_K)
        del cm
    summ["decision_metrics_matched"] = {"M0": an.decision_metrics(d0m), prim: an.decision_metrics(chosen(m_rows, prim))}
    summ["gene_level_matched"] = {"M0": an.gene_level(d0m), prim: an.gene_level(chosen(m_rows, prim))}
    summ["tier_stability_frac_genes_changing"] = an.tier_stability(m0n[m0n.method == "M0"])
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=p6._js), encoding="utf8")
    return summ


def finish_f2(frozen, odir, cache_dir, ext_dir, parts_dir, log, est_s=None, colstore_dir=None):
    """resume_v07.finish with the F2 column store; same outputs (summary.json, confirm_v07_result.json, parquets)."""
    import confirm_v07 as cf
    import methods_v06 as mv
    import pipeline_v06 as p6
    ev = p6.ev
    odir = Path(odir)
    if est_s is not None and not rs.can_start(est_s, "finish F2 (summarise + confirmatory analysis)", log):
        return rs.PAUSED
    t0 = time.time()
    files = rs._cache_files(cache_dir)
    pf = [Path(parts_dir) / f"stage2_split_{i:02d}.pkl" for i in range(len(files))]
    bad = [p.name for p in pf if not rs.is_done(p)]
    assert not bad, f"stage-2 parts missing or incomplete: {bad}"
    hcfg, dln, sln = frozen["h"], frozen["d"], frozen["s"]
    log(f"selected proxies: h={hcfg} d={dln} s={sln}")
    store, pt, covs, fb, Ks, ds = build_store(pf, colstore_dir or (odir / "_f2_colstore"), log)
    pt.to_parquet(odir / "proxy_rows.parquet")
    log(f"F2: column store built, {store.n} rows")
    hyp = {m: {"chosen": h} for m, h in frozen["hypers"].items()}
    prim, best_dist = frozen["primary"], frozen["best_dist"]
    prim_h = hyp[prim]["chosen"]
    log(f"primary {prim} {prim_h}; best M_dist {best_dist} {hyp[best_dist]['chosen']}")
    write_parquet(store, odir / "coverage_rows.parquet")
    covs.to_parquet(odir / "covariates_analysis_only.parquet")
    mv.M3_FALLBACKS.clear()
    for x in fb:
        for k, v in x.items():
            mv.M3_FALLBACKS[k] = mv.M3_FALLBACKS.get(k, 0) + v
    summarise_f2(odir, store, pt, covs, hcfg, None, dln, None, sln, None, hyp, prim, prim_h, best_dist, Ks, ds)
    log("F2: summary.json written")
    del covs, pt
    Kf = store.keyframe()
    g = store.take(np.flatnonzero(((Kf.model == ev.PRIMARY_MODEL) & (Kf.k == ev.PRIMARY_K) & (Kf.variant == "full")
                                   & (Kf["mode"] == "matched") & (Kf.perm == -1)).to_numpy()))
    del Kf, store
    caches = [rs._load_cache_with_ext(f, ext_dir) for f in files]
    out, rm = cf.analyze(caches, frozen, log=log, cov_grid=g)
    rm.to_parquet(odir / "confirm_rand_rows.parquet")
    cf.write(out, odir / "confirm_v07_result.json")
    rs.mark_done(odir / "confirm_v07_result.json", time.time() - t0)
    log(f"classification: {out['classification']}")
    return 0



# ------------------------------------------------------------------ F2 revision 2 ("F2c", 2026-10-15)
# Memory profile on dev (F2): the plateau comes from the object-dtype key frame (7 columns x all rows) kept during
# summarise, and the peak from criteria() receiving all full/<mode> rows (58% of the table). F2c changes only the
# representation, never a value:
#   * string columns are held as int32 codes + one array of the distinct (interned) str objects; take() returns
#     uniques[codes[pos]], i.e. the same str values;
#   * the key masks of summarise are evaluated on codes (== value  <=>  code == code(value); NaN -> code -1 never
#     matches), and the tier-table groupby is reproduced from codes: groups ordered like pandas groupby(sort=True)
#     (sorted unique values per key, lexicographic over keys, NaN keys dropped), rows inside a group in original
#     order;
#   * criteria() receives the rows of its own first filter (model == P, k == K, method in {M0, primary}) out of the
#     full/<mode> rows; every selection inside criteria() is a sub-filter of that, so it sees the same rows, in the
#     same order, with the same index labels;
#   * on Windows the process working set is trimmed (EmptyWorkingSet) after large blocks, so memory-mapped file
#     pages of the numeric columns are not counted by the RSS watchdog.


def _trim():
    try:
        import ctypes
        ctypes.windll.psapi.EmptyWorkingSet(ctypes.windll.kernel32.GetCurrentProcess())
    except Exception:
        pass


class CodeStore(ColumnStore):
    def __init__(self, base):
        self.cols, self.n = base.cols, base.n
        self.a, self.codes, self.uniq, self.na = {}, {}, {}, {}
        for c in base.cols:
            v = base.a[c]
            if isinstance(v, np.ndarray) and v.dtype == object:
                codes, uniq = pd.factorize(v, use_na_sentinel=True)
                u = np.empty(len(uniq), dtype=object)
                u[:] = list(uniq)
                self.codes[c] = codes.astype(np.int32)
                self.uniq[c] = u
                na = np.flatnonzero(codes < 0)                 # missing values: keep the original objects
                self.na[c] = (na, v[na].copy())
            else:
                self.a[c] = v
        del base

    def column(self, c, pos=None):
        if c in self.codes:
            pos_all = np.arange(self.n) if pos is None else np.asarray(pos, dtype=np.int64)
            cd = self.codes[c][pos_all]
            out = np.empty(len(cd), dtype=object)
            ok = cd >= 0
            out[ok] = self.uniq[c][cd[ok]]
            if (~ok).any():
                na_pos, na_obj = self.na[c]
                out[~ok] = na_obj[np.searchsorted(na_pos, pos_all[~ok])]
            return out
        return np.asarray(self.a[c] if pos is None else self.a[c][pos])

    def take(self, pos):
        pos = np.asarray(pos, dtype=np.int64)
        return pd.DataFrame({c: self.column(c, pos) for c in self.cols}, columns=self.cols,
                            index=pd.Index(pos, dtype="int64"))

    def series(self, c):
        return pd.Series(self.column(c), name=c)

    def eq(self, c, value):
        if c in self.codes:
            hit = np.flatnonzero(self.uniq[c] == value) if len(self.uniq[c]) else []
            return self.codes[c] == hit[0] if len(hit) else np.zeros(self.n, bool)
        return np.asarray(self.a[c]) == value

    def isin(self, c, values):
        idx = [i for i, u in enumerate(self.uniq[c]) if u in set(values)]
        return np.isin(self.codes[c], np.array(idx, dtype=np.int32))

    def groups(self, keys, row_mask):
        """(key tuple, positions) in the order of DataFrame.groupby(keys, sort=True, dropna=True) on the masked rows."""
        rows = np.flatnonzero(row_mask)
        ranks, vals = [], []
        keep = np.ones(len(rows), bool)
        for c in keys:
            if c in self.codes:
                cd = self.codes[c][rows]
                keep &= cd >= 0
                order = sorted(range(len(self.uniq[c])), key=lambda i: self.uniq[c][i])
                rank_of = np.empty(len(order), np.int64)
                rank_of[order] = np.arange(len(order))
                r = np.where(cd >= 0, rank_of[np.maximum(cd, 0)], -1)
                vals.append(np.array([self.uniq[c][i] for i in order], dtype=object))
            else:
                x = np.asarray(self.a[c])[rows]
                u, r = np.unique(x, return_inverse=True)
                vals.append(u)
            ranks.append(r)
        rows = rows[keep]
        ranks = [r[keep] for r in ranks]
        o = np.lexsort(ranks[::-1])                                   # stable: original order inside each group
        R = np.stack([r[o] for r in ranks], axis=1)
        brk = np.flatnonzero(np.any(np.diff(R, axis=0) != 0, axis=1)) + 1
        for s, e in zip(np.r_[0, brk], np.r_[brk, len(o)]):
            key = tuple(vals[j][R[s, j]] for j in range(len(keys)))
            yield key, rows[o[s:e]]


def summarise_f2c(out, store, pt, covs, hcfg, rho_h, dln, rho_d, sln, rho_s, hyp, prim, prim_h, best_dist,
                  chosen_K, chosen_d):
    """pipeline_v06.summarise on a CodeStore (same subsets, same calls, same order)."""
    import pipeline_v06 as p6
    an, ev = p6.an, p6.ev
    out = Path(out)
    P = ev.PRIMARY_MODEL
    S = store
    kk, pp = np.asarray(S.a["k"]), np.asarray(S.a["perm"])
    pm = lambda mode, model=P: S.take(np.flatnonzero(S.eq("model", model) & (kk == ev.PRIMARY_K)   # noqa: E731
                                                     & S.eq("variant", "full") & S.eq("mode", mode) & (pp == -1)))
    chosen = lambda d, meth: d[(d.method == meth) & (d.hyper == hyp[meth]["chosen"])]        # noqa: E731
    summ = {"selection": {"h": hcfg, "rho_h": rho_h, "d": dln, "rho_d": rho_d, "s": sln, "rho_s": rho_s,
                          "hypers": {m: v["chosen"] for m, v in hyp.items()}, "hyper_detail": hyp,
                          "primary": [prim, prim_h], "best_dist": [best_dist, hyp[best_dist]["chosen"]]},
            "chosen_K": chosen_K, "chosen_d": chosen_d}
    tabs = {}
    for i, ((m, k, var, mode, meth, hy), g) in enumerate(S.groups(["model", "k", "variant", "mode", "method", "hyper"],
                                                                 pp == -1)):
        tabs[f"{m}|{k}|{var}|{mode}|{meth}|{hy}"] = an.tier_table(S.take(g))
        if i % 50 == 49:
            _trim()
    summ["tier_tables"] = tabs
    _trim()
    summ["monotone_all"] = bool(S.series("monotone").all())
    import methods_v06 as mv
    summ["m3_solver_fallbacks"] = dict(mv.M3_FALLBACKS)
    m_rows = pm("matched")
    m0n = pm("nominal")
    pmc = {"nominal": m0n, "matched": m_rows}
    summ["dgo_strata"] = {f"{mode}|{meth}": an.strata_table(chosen(pmc[mode], meth) if meth != "M0"
                                                              else pmc[mode][pmc[mode].method == "M0"], "D_GO")
                          for mode in ("nominal", "matched")
                          for meth in ["M0", "M1", prim, "Mscore-N"] + [b for b in an.BASELINES if b in hyp]}
    d0m = m_rows[m_rows.method == "M0"].sort_values(["split", "target"])
    comp = {}
    for meth in [prim, best_dist, "Mcomb-N", "Mcomb-3", "Mscore-N", "oracle"]:
        dm = chosen(m_rows, meth).sort_values(["split", "target"])
        comp[f"{meth}_vs_M0"] = an.paired_T5(dm, d0m)
    dbd = chosen(m_rows, best_dist).sort_values(["split", "target"])
    for meth in [prim, "Mcomb-N", "Mcomb-3", "Mscore-N"]:
        comp[f"{meth}_vs_{best_dist}"] = an.paired_T5(chosen(m_rows, meth).sort_values(["split", "target"]), dbd)
    bases = [b for b in an.BASELINES if b in hyp]
    for b in bases:
        db = chosen(m_rows, b).sort_values(["split", "target"])
        comp[f"{b}_vs_M0"] = an.paired_T5(db, d0m)
        for meth in sorted({prim, "Mscore-N"}):
            if meth != b:
                comp[f"{meth}_vs_{b}"] = an.paired_T5(chosen(m_rows, meth).sort_values(["split", "target"]), db)
    summ["matched_paired_T5"] = comp
    dec = {}
    for meth in [prim, "Mscore-N"] + bases:
        dec[f"{meth}_vs_M0"] = an.decision_paired(chosen(m_rows, meth), d0m)
    for meth in sorted({prim, "Mscore-N"}):
        dec[f"{meth}_vs_{best_dist}"] = an.decision_paired(chosen(m_rows, meth), chosen(m_rows, best_dist))
    summ["decision_paired_top10_matched"] = dec
    rows_lin = m0n.iloc[:0]
    for m in ev.LINEAR:
        x = pm("nominal", m)
        rows_lin = pd.concat([rows_lin, x[x.method == "M0"]])
        del x
    summ["confounding"] = an.confounding(an.confounding_frame(rows_lin, covs))
    hr = S.take(np.flatnonzero((kk == ev.PRIMARY_K) & S.eq("variant", "half") & S.eq("mode", "nominal")
                               & S.eq("method", "M0") & S.isin("model", ev.LINEAR) & (pp == -1)))
    summ["confounding_half"] = an.confounding(an.confounding_frame(hr, covs, half=True))
    del hr
    _trim()
    summ["incremental"] = an.incremental(rows_lin, covs)
    pth = pt[(pt.kind == "h") & (pt.proxy == hcfg)].rename(columns={"value": "h"})
    for mode in ("nominal", "matched"):
        # rows of criteria()'s own first filter, out of cov[(variant == "full") & (mode == mode)]
        cm = S.take(np.flatnonzero(S.eq("variant", "full") & S.eq("mode", mode) & S.eq("model", P)
                                   & (kk == ev.PRIMARY_K) & S.isin("method", ["M0", prim])))
        summ[f"criteria_{mode}"] = an.criteria(cm, pth, prim, prim_h, P, ev.PRIMARY_K)
        del cm
        _trim()
    summ["decision_metrics_matched"] = {"M0": an.decision_metrics(d0m), prim: an.decision_metrics(chosen(m_rows, prim))}
    summ["gene_level_matched"] = {"M0": an.gene_level(d0m), prim: an.gene_level(chosen(m_rows, prim))}
    summ["tier_stability_frac_genes_changing"] = an.tier_stability(m0n[m0n.method == "M0"])
    (out / "summary.json").write_text(json.dumps(summ, indent=1, default=p6._js), encoding="utf8")
    return summ


def finish_f2c(frozen, odir, cache_dir, ext_dir, parts_dir, log, est_s=None, colstore_dir=None):
    """finish_f2 with the CodeStore representation (F2c)."""
    import confirm_v07 as cf
    import methods_v06 as mv
    import pipeline_v06 as p6
    ev = p6.ev
    odir = Path(odir)
    if est_s is not None and not rs.can_start(est_s, "finish F2c (summarise + confirmatory analysis)", log):
        return rs.PAUSED
    t0 = time.time()
    files = rs._cache_files(cache_dir)
    pf = [Path(parts_dir) / f"stage2_split_{i:02d}.pkl" for i in range(len(files))]
    bad = [p.name for p in pf if not rs.is_done(p)]
    assert not bad, f"stage-2 parts missing or incomplete: {bad}"
    hcfg, dln, sln = frozen["h"], frozen["d"], frozen["s"]
    log(f"selected proxies: h={hcfg} d={dln} s={sln}")
    base, pt, covs, fb, Ks, ds = build_store(pf, colstore_dir or (odir / "_f2_colstore"), log)
    pt.to_parquet(odir / "proxy_rows.parquet")
    write_parquet(base, odir / "coverage_rows.parquet")
    store = CodeStore(base)
    del base
    _trim()
    log(f"F2c: code store built, {store.n} rows")
    hyp = {m: {"chosen": h} for m, h in frozen["hypers"].items()}
    prim, best_dist = frozen["primary"], frozen["best_dist"]
    prim_h = hyp[prim]["chosen"]
    log(f"primary {prim} {prim_h}; best M_dist {best_dist} {hyp[best_dist]['chosen']}")
    covs.to_parquet(odir / "covariates_analysis_only.parquet")
    mv.M3_FALLBACKS.clear()
    for x in fb:
        for k, v in x.items():
            mv.M3_FALLBACKS[k] = mv.M3_FALLBACKS.get(k, 0) + v
    summarise_f2c(odir, store, pt, covs, hcfg, None, dln, None, sln, None, hyp, prim, prim_h, best_dist, Ks, ds)
    log("F2c: summary.json written")
    del covs, pt
    kk, pp = np.asarray(store.a["k"]), np.asarray(store.a["perm"])
    g = store.take(np.flatnonzero(store.eq("model", ev.PRIMARY_MODEL) & (kk == ev.PRIMARY_K) & store.eq("variant", "full")
                                  & store.eq("mode", "matched") & (pp == -1)))
    del store, kk, pp
    _trim()
    caches = [rs._load_cache_with_ext(f, ext_dir) for f in files]
    out, rm = cf.analyze(caches, frozen, log=log, cov_grid=g)
    rm.to_parquet(odir / "confirm_rand_rows.parquet")
    cf.write(out, odir / "confirm_v07_result.json")
    rs.mark_done(odir / "confirm_v07_result.json", time.time() - t0)
    log(f"classification: {out['classification']}")
    return 0
