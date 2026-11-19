# download_resumable.py -- resumable HTTP download with MD5 verification (engineering utility).
# Usage: python download_resumable.py <out_dir> <name>=<md5> [<name>=<md5> ...]
# Files come from Zenodo record 13350497 (scPerturb v1.4). Partial data is kept in <name>.part and resumed with
# HTTP Range requests after an interruption (network or power cut). After completion the MD5 is checked; on a
# mismatch the file is renamed to <name>.bad and the script exits with code 2. Progress goes to download_log.txt.
import hashlib
import os
import sys
import time
from pathlib import Path

import requests

BASE = "https://zenodo.org/records/13350497/files/{}?download=1"
CHUNK = 1 << 20


def log(out, msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with open(out / "download_log.txt", "a", encoding="utf8") as fh:
        fh.write(line + "\n")


def md5_of(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def fetch(out, name, md5):
    final = out / name
    if final.exists():
        if md5_of(final) == md5:
            log(out, f"{name}: already complete, MD5 ok")
            return True
        log(out, f"{name}: existing file has wrong MD5; renaming to .bad")
        os.replace(final, out / (name + ".bad"))
    part = out / (name + ".part")
    url = BASE.format(name)
    total = int(requests.head(url, allow_redirects=True, timeout=60).headers.get("Content-Length", 0))
    fails = 0
    while True:
        have = part.stat().st_size if part.exists() else 0
        if total and have >= total:
            break
        try:
            with requests.get(url, stream=True, timeout=120, headers={"Range": f"bytes={have}-"}) as r:
                if r.status_code not in (200, 206):
                    raise IOError(f"HTTP {r.status_code}")
                if r.status_code == 200 and have > 0:
                    have = 0
                    part.unlink()
                mode = "ab" if have else "wb"
                t0, n, last = time.time(), 0, time.time()
                with open(part, mode) as fh:
                    for ch in r.iter_content(CHUNK):
                        fh.write(ch)
                        n += len(ch)
                        if time.time() - last > 60:
                            fh.flush()
                            os.fsync(fh.fileno())
                            done = have + n
                            log(out, f"{name}: {done / 1e6:.0f}/{total / 1e6:.0f} MB ({100 * done / max(total, 1):.1f}%), "
                                     f"{n / 1e6 / (time.time() - t0):.2f} MB/s")
                            last = time.time()
                    fh.flush()
                    os.fsync(fh.fileno())
            fails = 0
        except Exception as e:  # network error: wait and resume from the bytes already on disk
            fails += 1
            log(out, f"{name}: interrupted ({type(e).__name__}: {str(e)[:120]}); retry {fails} in 30 s")
            if fails > 40:
                log(out, f"{name}: giving up after {fails} failures")
                return False
            time.sleep(30)
    got = md5_of(part)
    if got != md5:
        os.replace(part, out / (name + ".bad"))
        log(out, f"{name}: MD5 MISMATCH (got {got}, expected {md5}); file renamed to .bad")
        return False
    os.replace(part, final)
    log(out, f"{name}: complete, {final.stat().st_size:,} bytes, MD5 ok")
    return True


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    ok = True
    for spec in sys.argv[2:]:
        name, md5 = spec.split("=")
        ok = fetch(out, name, md5) and ok
    log(out, "ALL DONE" if ok else "FINISHED WITH ERRORS")
    sys.exit(0 if ok else 2)


if __name__ == "__main__":
    main()
