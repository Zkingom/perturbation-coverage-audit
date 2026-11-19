"""final_summary.py — 全部结束后生成一页总结 /root/autodl-tmp/gears/FINAL_SUMMARY.txt（只读，不改任何结果）。"""
import hashlib, json, os, time
from pathlib import Path

G = Path("/root/autodl-tmp/gears"); B = G / "gears_bundle_v07"; A = Path("/root/autodl-tmp")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def rd(p):
    try:
        return Path(p).read_text()
    except Exception:
        return "(missing)"


L = [f"GEARS v0.7 run summary — generated {time.strftime('%F %T')}", ""]
done, missing = [], []
for i in range(21):
    d = B / "results" / f"split_{i:02d}" / "done.json"
    if d.exists():
        j = json.loads(d.read_text())
        done.append(i)
        L.append(f"split {i:02d}: DONE  epochs={j['hp']['epochs']} hidden={j['hp']['hidden_size']} unc={j['hp']['uncertainty']} "
                 f"bs={j['hp']['batch_size']} seed={j['seed']} finite={j['finite']} n_pred={j['n_pred']} "
                 f"peak_rss={j['peak_rss_gb']}GB total={round(j['total_s']/60)}min")
    else:
        missing.append(i)
        L.append(f"split {i:02d}: NOT DONE")
L += ["", f"DONE {len(done)}/21; missing: {missing if missing else 'none'}", ""]
L += ["--- run_all_timing.txt ---", rd(B / "logs" / "run_all_timing.txt").strip(), ""]
L += ["--- verify (code/gears_server_run.py verify) ---", rd(B / "logs" / "verify.txt").strip()[-1500:], ""]
L += ["--- bundle check ---", rd(G / "AUDIT" / "bundle_check.txt").strip(), ""]
L += ["--- postchecks ---", rd(G / "POSTCHECK" / "postchecks.txt").strip(), ""]
L += ["--- coexp manifest (last line) ---", rd(G / "AUDIT" / "coexp_manifest.tsv").strip().splitlines()[-1] if (G / "AUDIT" / "coexp_manifest.tsv").exists() else "(missing)", ""]
L.append("--- files to download ---")
for p in [B / "gears_results_v07.tar.gz", B / "gears_results_v07.tar.gz.sha256", G / "gears_audit_v07.tar.gz", G / "gears_audit_v07.tar.gz.sha256"]:
    L.append(f"{p}  {'%.1f MB' % (p.stat().st_size / 1e6) if p.exists() else 'MISSING'}" + (f"  sha256={sha(p)}" if p.exists() and p.suffix == '.gz' else ""))
L += ["", "READY TO DOWNLOAD" if (B / "gears_results_v07.tar.gz").exists() and (G / "gears_audit_v07.tar.gz").exists() else "NOT READY"]
(G / "FINAL_SUMMARY.txt").write_text("\n".join(L) + "\n")
print("\n".join(L))
