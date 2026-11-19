# gears_server_run.py -- GEARS 0.1.2 robustness run on RPE1 essential (protocol v0.6 section 10, inherited by v0.7).
# Runs on the rented GPU server inside the bundle directory (venv from setup_gears_env.sh).
#   python code/gears_server_run.py prepare            once: GEARS new_data_process (DE genes, per-cell graphs)
#   python code/gears_server_run.py run --splits 0     smoke test (one split)
#   python code/gears_server_run.py run                all 21 splits; finished splits are verified and skipped
#   python code/gears_server_run.py verify             re-hash every finished split, write results/SHA256SUMS
# GEARS settings are the protocol's: hidden_size = 64, epochs = 20, uncertainty = False, all other model_initialize
# defaults (num_similar_genes_go_graph = 20, num_similar_genes_co_express_graph = 20, coexpress_threshold = 0.4); best
# epoch = lowest validation top-20 DE MSE (GEARS default); batch sizes 32 / 128 and custom-split fields as in the
# Ahlmann-Eltze et al. 2025 benchmark script (run_gears.py). Each split has its own GEARS dataset directory so the
# co-expression graph (cached by GEARS under a file name that does not contain the split) is never shared.
# Resumability: a split's outputs are written to *.tmp, fsync'ed and renamed; done.json (sha256 of every output) is
# written last. A split without a valid done.json is recomputed from scratch.
import argparse
import hashlib
import json
import os
import random
import shutil
import sys
import threading
import time
from pathlib import Path

import numpy as np

B = Path(__file__).resolve().parents[1]
WORK = B / "work" / "gears_data"
RES = B / "results"
DS = "rpe1_cap16"
PREP_SEED = 20261001 + 9701
HP = {"hidden_size": 64, "epochs": 20, "uncertainty": False, "batch_size": 32, "test_batch_size": 128}


def log(msg, f=None):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    if f is not None:
        with open(f, "a") as fh:
            fh.write(line + "\n")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def atomic_bytes(path, data):
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "wb") as fh:
        fh.write(data)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def atomic_np(path, arr):
    tmp = Path(str(path) + ".tmp.npy")
    with open(tmp, "wb") as fh:
        np.save(fh, arr)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


class Peak:
    """Samples host RSS (process + children) every second; GPU peaks come from torch."""
    def __init__(self):
        import psutil
        self.p, self.peak, self.stop = psutil.Process(), 0, False
        threading.Thread(target=self._run, daemon=True).start()

    def _run(self):
        while not self.stop:
            try:
                r = self.p.memory_info().rss + sum(c.memory_info().rss for c in self.p.children(recursive=True))
                self.peak = max(self.peak, r)
            except Exception:
                pass
            time.sleep(1)


def seed_all(s):
    import torch
    random.seed(s)
    np.random.seed(s)
    torch.manual_seed(s)
    torch.cuda.manual_seed_all(s)


def place_priors():
    WORK.mkdir(parents=True, exist_ok=True)
    for f in ("gene2go_all.pkl", "essential_all_data_pert_genes.pkl"):
        if not (WORK / f).exists():
            shutil.copy2(B / "priors" / f, WORK / f)
    if not (WORK / "go_essential_all").exists():
        shutil.copytree(B / "priors" / "go_essential_all", WORK / "go_essential_all")


def prepare():
    import scanpy as sc
    from gears import PertData
    import gears.version
    assert gears.version.__version__ == "0.1.2"
    place_priors()
    done = WORK / DS / "prepare_done.json"
    if done.exists():
        d = json.loads(done.read_text())
        if all(sha256(WORK / DS / k) == v for k, v in d["sha256"].items()):
            log("prepare: already done and verified")
            return 0
        log("prepare: marker does not match the files -> rebuilding")
    shutil.rmtree(WORK / DS, ignore_errors=True)
    t0 = time.time()
    pk = Peak()
    adata = sc.read_h5ad(B / "data" / "rpe1_cap16.h5ad")
    seed_all(PREP_SEED)
    pd_ = PertData(str(WORK))                       # priors are pre-placed: no dataverse download
    pd_.new_data_process(dataset_name=DS, adata=adata)
    files = ["perturb_processed.h5ad", "data_pyg/cell_graphs.pkl"]
    rec = {"sha256": {f: sha256(WORK / DS / f) for f in files}, "seconds": round(time.time() - t0),
           "peak_rss_gb": round(pk.peak / 2**30, 2), "seed": PREP_SEED}
    atomic_bytes(done, json.dumps(rec, indent=1).encode())
    log(f"prepare done {json.dumps(rec)}")
    return 0


def split_done(i):
    d = RES / f"split_{i:02d}"
    m = d / "done.json"
    if not m.exists():
        return None
    rec = json.loads(m.read_text())
    try:
        ok = all(sha256(d / f) == h for f, h in rec["sha256"].items())
    except FileNotFoundError:
        ok = False
    return rec if ok else None


def run_split(i, seed_mode, logf):
    import torch
    from gears import GEARS, PertData
    t = {}
    t0 = time.time()
    out = RES / f"split_{i:02d}"
    shutil.rmtree(out, ignore_errors=True)          # anything without a valid done.json is discarded
    out.mkdir(parents=True)
    s2c = json.loads((B / "splits" / f"set2conditions_{i:02d}.json").read_text())
    assert not (set(s2c["train"]) & set(s2c["val"])) and not (set(s2c["train"]) & set(s2c["test"])) \
        and not (set(s2c["val"]) & set(s2c["test"])), "folds overlap"
    sdir = WORK / f"rpe1_s{i:02d}"                  # own dataset dir per split -> own co-expression cache
    shutil.rmtree(sdir, ignore_errors=True)
    (sdir / "data_pyg").mkdir(parents=True)
    os.symlink(WORK / DS / "perturb_processed.h5ad", sdir / "perturb_processed.h5ad")
    os.symlink(WORK / DS / "data_pyg" / "cell_graphs.pkl", sdir / "data_pyg" / "cell_graphs.pkl")
    if seed_mode == "benchmark":                    # run_gears.py: np.random.seed(1); torch seed 0 is set by `import gears`
        seed = 1
        np.random.seed(1)
        torch.manual_seed(0)
        torch.cuda.manual_seed_all(0)
    else:
        seed = json.loads((B / "splits" / f"conf_split_{i:02d}.json").read_text())["seed"]
        seed_all(seed)
    torch.cuda.reset_peak_memory_stats()
    pk = Peak()
    pdat = PertData(str(WORK))
    pdat.load(data_path=str(sdir))
    present = set(pdat.adata.obs.condition.unique())
    assert set(s2c["train"]) | set(s2c["val"]) | set(s2c["test"]) == present, "conditions differ from the data"
    pdat.set2conditions = s2c
    pdat.split = "custom"
    pdat.subgroup = None
    pdat.seed = 1
    pdat.train_gene_set_size = 0.75
    pdat.get_dataloader(batch_size=HP["batch_size"], test_batch_size=HP["test_batch_size"])
    t["load_s"] = round(time.time() - t0)
    model = GEARS(pdat, device="cuda")
    model.model_initialize(hidden_size=HP["hidden_size"], uncertainty=HP["uncertainty"])
    t1 = time.time()
    model.train(epochs=HP["epochs"])
    t["train_s"] = round(time.time() - t1)
    t2 = time.time()
    conds = s2c["val"] + s2c["test"]
    genes = [c.split("+")[0] for c in conds]
    pred = model.predict([[g] for g in genes])
    P = np.vstack([np.asarray(pred[g], dtype=np.float32) for g in genes])
    t["predict_s"] = round(time.time() - t2)
    gnames = pdat.adata.var.gene_name.astype(str).tolist()
    gids = pdat.adata.var.index.astype(str).tolist()
    atomic_np(out / "pred.npy", P)
    atomic_bytes(out / "pred_targets.json", json.dumps({"val": [c.split("+")[0] for c in s2c["val"]],
                                                         "test_cal_and_test": [c.split("+")[0] for c in s2c["test"]],
                                                         "row_order": genes}).encode())
    atomic_bytes(out / "genes.json", json.dumps({"gene_name": gnames, "gene_id": gids}).encode())
    tmpm = out / "model"
    model.save_model(str(tmpm))
    pk.stop = True
    rec = {"split": i, "seed_mode": seed_mode, "seed": seed, "hp": HP, "timing_s": t, "total_s": round(time.time() - t0),
           "peak_rss_gb": round(pk.peak / 2**30, 2),
           "peak_gpu_alloc_gb": round(torch.cuda.max_memory_allocated() / 2**30, 2),
           "peak_gpu_reserved_gb": round(torch.cuda.max_memory_reserved() / 2**30, 2),
           "gpu": torch.cuda.get_device_name(0), "n_pred": int(P.shape[0]), "n_genes": int(P.shape[1]),
           "finite": bool(np.isfinite(P).all())}
    assert rec["finite"], "non-finite predictions"
    files = ["pred.npy", "pred_targets.json", "genes.json"] + [str(p.relative_to(out)) for p in sorted(tmpm.rglob("*")) if p.is_file()]
    rec["sha256"] = {f: sha256(out / f) for f in files}
    atomic_bytes(out / "done.json", json.dumps(rec, indent=1).encode())
    log(f"split {i:02d} done {json.dumps({k: v for k, v in rec.items() if k != 'sha256'})}", logf)
    return rec


def run(splits, seed_mode):
    if not (WORK / DS / "prepare_done.json").exists():
        log("run: prepare first")
        return 2
    RES.mkdir(parents=True, exist_ok=True)
    logf = RES / "run.log"
    for i in splits:
        if split_done(i):
            log(f"split {i:02d}: verified done, skipped", logf)
            continue
        log(f"split {i:02d}: start", logf)
        run_split(i, seed_mode, logf)
    verify()
    return 0


def verify():
    lines, n = [], 0
    for i in range(21):
        rec = split_done(i)
        if rec:
            n += 1
            lines += [f"{h}  split_{i:02d}/{f}" for f, h in rec["sha256"].items()]
            lines.append(f"{sha256(RES / f'split_{i:02d}' / 'done.json')}  split_{i:02d}/done.json")
    if lines:
        atomic_bytes(RES / "SHA256SUMS", ("\n".join(lines) + "\n").encode())
    log(f"verify: {n}/21 splits done and hash-verified")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["prepare", "run", "verify"])
    ap.add_argument("--splits", type=int, nargs="*", default=list(range(21)))
    ap.add_argument("--seed-mode", choices=["benchmark", "split"], default="benchmark",
                    help="benchmark: numpy seed 1, torch seed 0 for every split (as run_gears.py + gears import); split: the split's own seed")
    a = ap.parse_args()
    if a.mode == "prepare":
        sys.exit(prepare())
    if a.mode == "verify":
        sys.exit(verify())
    sys.exit(run(a.splits, a.seed_mode))
