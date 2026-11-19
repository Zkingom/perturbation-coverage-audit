# m3_crosscheck_v061.py -- protocol v0.6.1 section 11c: our exact one-LP M3 vs conditionalconformal 0.0.5.
#   python m3_crosscheck_v061.py make    <dir>   (env pert-audit)  synthetic cal/test sets, Phi matrices, U
#   python m3_crosscheck_v061.py ref     <dir>   (env condconf-ref) reference thresholds from the package
#   python m3_crosscheck_v061.py compare <dir>   (env pert-audit)  our thresholds + acceptance criteria
# Both implementations read the same Phi matrices, so the basis functions are identical.
import json
import sys
from pathlib import Path

import numpy as np

ALPHA = 0.10
N_CAL, N_TEST = 267, 60
TOL_99, TOL_ALL, TOL_BISECT = 1e-6, 1e-4, 2e-3
MAX_ERR_FRAC = 0.05


def make(d):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import methods_v06 as mv
    d.mkdir(parents=True, exist_ok=True)
    meta = []
    for seed in (0, 1, 2):
        rng = np.random.default_rng(seed)
        hc = rng.uniform(0.5, 3, N_CAL)
        dc = hc * np.exp(rng.normal(0, 0.5, N_CAL))
        s = hc * rng.exponential(1, N_CAL)
        ht = np.linspace(0.3, 3.5, N_TEST)
        dt = ht * np.exp(rng.normal(0, 0.5, N_TEST))
        ho = rng.uniform(0.5, 3, 800)
        do = ho * np.exp(rng.normal(0, 0.5, 800))
        U = rng.uniform(-ALPHA, 1 - ALPHA, N_TEST)
        for kind in ("quad", "ns3", "ns5", "lin2"):
            if kind == "lin2":
                f = mv.phi_lin2(ho, do)
                Pc, Pt = f(np.column_stack([hc, dc])), f(np.column_stack([ht, dt]))
            else:
                f = mv.phi_basis(kind, ho)
                Pc, Pt = f(hc), f(ht)
            np.savez(d / f"case_s{seed}_{kind}.npz", s=s, Pc=Pc, Pt=Pt, U=U)
            meta.append(f"case_s{seed}_{kind}.npz")
    (d / "cases.json").write_text(json.dumps(meta))


def ref(d):
    import importlib.metadata as md
    from conditionalconformal import CondConf
    out = {"versions": {p: md.version(p) for p in ("conditionalconformal", "numpy", "scipy", "cvxpy")}, "cases": {}}
    for name in json.loads((d / "cases.json").read_text()):
        z = np.load(d / name)
        s, Pc, Pt, U = z["s"], z["Pc"], z["Pt"], z["U"]
        PHI = np.vstack([Pc, Pt])
        phi_fn = lambda x: PHI[np.asarray(x).reshape(-1).astype(int)]          # noqa: E731
        cc = CondConf(score_fn=lambda x, y: y, Phi_fn=phi_fn)
        cc.setup_problem(np.arange(len(s)).reshape(-1, 1).astype(float), s)
        ex, bs, err = [], [], []
        for i in range(len(Pt)):
            x = np.array([[float(len(s) + i)]])
            try:
                ex.append(float(np.ravel(cc.predict(1 - ALPHA, x, lambda c, x: c, threshold=U[i], exact=True))[0]))
                err.append("")
            except Exception as e:                                          # noqa: BLE001
                ex.append(np.nan)
                err.append(repr(e)[:200])
            try:
                bs.append(float(np.ravel(cc.predict(1 - ALPHA, x, lambda c, x: c, threshold=U[i], exact=False))[0]))
            except Exception as e:                                          # noqa: BLE001
                bs.append(np.nan)
        out["cases"][name] = {"exact": ex, "bisect": bs, "errors": err}
    (d / "ref_thresholds.json").write_text(json.dumps(out, default=float))


def compare(d):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import methods_v06 as mv
    R = json.loads((d / "ref_thresholds.json").read_text())
    rows, n_err, n_tot = [], 0, 0
    inf_ok = True
    rel_ex, rel_bs = [], []
    for name, v in R["cases"].items():
        z = np.load(d / name)
        ours = np.array([mv.m3_threshold(z["Pc"], z["s"], z["Pt"][i], z["U"][i], ALPHA) for i in range(len(z["U"]))])
        ex = np.array(v["exact"], float)
        bs = np.array(v["bisect"], float)
        for i in range(len(ours)):
            n_tot += 1
            if np.isnan(ex[i]):
                n_err += 1
                if np.isfinite(bs[i]) and np.isfinite(ours[i]):
                    rel_bs.append(abs(ours[i] - bs[i]) / max(1, abs(bs[i])))
                elif not (np.isinf(bs[i]) and np.isinf(ours[i]) and np.sign(bs[i]) == np.sign(ours[i])):
                    rel_bs.append(np.inf)
                continue
            if np.isinf(ex[i]):
                inf_ok &= bool(np.isinf(ours[i]) and np.sign(ours[i]) == np.sign(ex[i]))
                continue
            rel_ex.append(abs(ours[i] - ex[i]) / max(1, abs(ex[i])) if np.isfinite(ours[i]) else np.inf)
        rows.append({"case": name, "max_rel_exact": float(np.nanmax([abs(o - e) / max(1, abs(e)) for o, e in zip(ours, ex)
                                                                        if np.isfinite(e) and np.isfinite(o)] or [0])),
                     "n_pkg_errors": int(np.isnan(ex).sum()), "n_inf_ref": int(np.isinf(ex).sum())})
    rel_ex = np.array(rel_ex)
    crit = {"c1_frac_within_1e-6": float((rel_ex <= TOL_99).mean()) if len(rel_ex) else None,
            "c1_max_rel": float(rel_ex.max()) if len(rel_ex) else None,
            "c1_pass": bool(len(rel_ex) and (rel_ex <= TOL_99).mean() >= 0.99 and rel_ex.max() <= TOL_ALL),
            "c2_infinite_agree": inf_ok,
            "c3_pkg_error_frac": n_err / n_tot,
            "c3_bisect_max_rel": float(max(rel_bs)) if rel_bs else None,
            "c3_pass": bool(n_err / n_tot <= MAX_ERR_FRAC and (not rel_bs or max(rel_bs) <= TOL_BISECT)),
            "n_points": n_tot, "n_compared_exact": int(len(rel_ex))}
    crit["pass"] = bool(crit["c1_pass"] and crit["c2_infinite_agree"] and crit["c3_pass"])
    out = {"versions_ref": R["versions"], "criteria": crit, "per_case": rows}
    (d / "m3_crosscheck_result.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(crit))


if __name__ == "__main__":
    {"make": make, "ref": ref, "compare": compare}[sys.argv[1]](Path(sys.argv[2]))
