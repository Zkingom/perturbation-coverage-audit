"""audit_pack.py — 给本机审核用的材料包（只读，不改任何结果）。
产出 /root/autodl-tmp/gears/gears_audit_v07.tar.gz (+ .sha256)，内容：
  runtime_patches/        实际注入的补丁源码（含 sitecustomize.py）
  run_all.sh              实际启动脚本（补丁如何注入：PYTHONPATH=runtime_patches -> sitecustomize）
  test_mem_patch.py/.out  evaluate 等价性测试及输出
  coexp_manifest.tsv      每个 split -> 它的共表达网络缓存文件、SHA256、训练扰动数、训练集指纹
  netcache/               全部共表达网络缓存文件 + GO 网络 pickle
  pip_freeze.txt          venv 的 pip freeze
  bundle_check.txt        sha256sum -c SHA256SUMS 结果（证明 code/ 等未改动）
  logs_runtime/           run_all.log、mem_run.txt、各 split 的 pre_XX.txt
"""
import hashlib, json, os, shutil, subprocess, tarfile
from pathlib import Path

B = Path("/root/autodl-tmp/gears/gears_bundle_v07")
G = Path("/root/autodl-tmp/gears")
A = Path("/root/autodl-tmp")
OUT = G / "AUDIT"
shutil.rmtree(OUT, ignore_errors=True); OUT.mkdir()


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


shutil.copytree(G / "runtime_patches", OUT / "runtime_patches", ignore=shutil.ignore_patterns("__pycache__"))
for f in ["run_all.sh", "test_mem_patch.py", "test_mem_patch.out", "run_all.log", "mem_run.txt"]:
    if (A / f).exists():
        shutil.copy2(A / f, OUT / f)
(OUT / "logs_runtime").mkdir()
for p in sorted((B / "logs").glob("pre_*.txt")):
    shutil.copy2(p, OUT / "logs_runtime" / p.name)
if (G / "netcache").exists():
    shutil.copytree(G / "netcache", OUT / "netcache")
if (G / "POSTCHECK").exists():
    shutil.copytree(G / "POSTCHECK", OUT / "POSTCHECK")
for f in ["postchecks.py", "audit_pack.py", "test_analysis_patch.py", "test_analysis_patch.out", "debug_run.txt", "debug_mem.txt", "debug_run.sh", "final_summary.py"]:
    if (G / f).exists():
        shutil.copy2(G / f, OUT / f)

rows = ["split\ttrain_n\tval_n\ttest_n\tdisjoint\ttrain_sha1_16\tcoexp_file\texists\tsha256"]
names = {}
for i in range(21):
    s2c = json.loads((B / "splits" / f"set2conditions_{i:02d}.json").read_text())
    tr, va, te = set(s2c["train"]), set(s2c["val"]), set(s2c["test"])
    h = hashlib.sha1("\n".join(sorted(s2c["train"])).encode()).hexdigest()[:16]
    fn = f"coexp_{h}_custom_1_0.75_0.4_20.csv"   # 与 mem_patch.coexp_cache_name 相同的命名规则
    p = G / "netcache" / fn
    names.setdefault(fn, []).append(i)
    rows.append(f"{i:02d}\t{len(tr)}\t{len(va)}\t{len(te)}\t{not (tr & va) and not (tr & te) and not (va & te)}\t{h}\t{fn}\t{p.exists()}\t{sha(p) if p.exists() else ''}")
shared = {k: v for k, v in names.items() if len(v) > 1}
rows.append(f"# distinct coexp files: {len(names)} / 21 ; files shared by >1 split: {shared if shared else 'none'}")
(OUT / "coexp_manifest.tsv").write_text("\n".join(rows) + "\n")

(OUT / "pip_freeze.txt").write_text(subprocess.run([str(B / "venv_gears/bin/python"), "-m", "pip", "freeze"],
                                                   capture_output=True, text=True).stdout)
r = subprocess.run(["sha256sum", "-c", "SHA256SUMS", "--quiet"], cwd=B, capture_output=True, text=True)
(OUT / "bundle_check.txt").write_text(f"rc={r.returncode}\n{r.stdout}{r.stderr}" + ("BUNDLE OK (code/, data/, priors/, splits/ unchanged)\n" if r.returncode == 0 else ""))

tgz = G / "gears_audit_v07.tar.gz"
with tarfile.open(tgz, "w:gz") as t:
    t.add(OUT, arcname="AUDIT")
(G / "gears_audit_v07.tar.gz.sha256").write_text(f"{sha(tgz)}  gears_audit_v07.tar.gz\n")
print("AUDIT PACKED", tgz, os.path.getsize(tgz))
