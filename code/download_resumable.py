# download_resumable.py -- parallel, resumable HTTP Range download of one Figshare file.
# (The previous single-stream version is kept in code/archive/download_resumable_v1_single_stream.py.)
#
# Usage:
#   python download_resumable.py --file-id 35775507 --dest <absolute path of final file>
#          [--threads 16] [--chunk-mib 32] [--limit-bytes N]
#
# Design
# - Expected size and MD5 come from the Figshare API (article 20029387 v1); the size is also checked
#   against the Content-Range total reported by S3.
# - Every request goes to the ndownloader URL again, so each one gets a fresh pre-signed S3 redirect
#   (the signature is only valid for ~10 s, and is checked when the request starts).
# - Missing byte ranges are cut into chunks and handed out in ascending order to N worker threads.
#   Each chunk is written at its own offset in <dest>.part. Ascending order keeps all writes close to
#   the end of the file, so the file system never has to zero-fill a large gap (relevant on NTFS/exFAT).
# - Finished byte intervals are recorded (after fsync) in <dest>.part.state.json, so a restart only
#   fetches what is missing; the chunk size can differ between runs.
# - A failed chunk is retried with exponential back-off (HTTP 429/5xx, timeouts, short reads).
# - Every 60 s one line with progress, speed over the last minute, mean speed and ETA goes to
#   <dest>.download.log.
# - Creating the file <dest>.stop makes all workers finish their current chunk and exit cleanly.
# - When the whole file is present: size, MD5 and SHA256 are computed; <dest>.part is renamed to
#   <dest> only if the size matches and the MD5 matches the Figshare MD5 (when one is published).
import argparse
import hashlib
import json
import os
import shutil
import threading
import time
from pathlib import Path

import requests

ARTICLE = 20029387
API = f"https://api.figshare.com/v2/articles/{ARTICLE}"
NDL = "https://ndownloader.figshare.com/files/{fid}"
READ_TIMEOUT = 60
MAX_TRIES_PER_CHUNK = 30
LOG_EVERY = 60


def utc():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def keep_awake(on):
    """Windows only: ask the OS not to enter idle sleep while this process runs (the display may
    still turn off). Has no effect on lid-close or manual sleep."""
    if os.name != "nt":
        return False
    import ctypes
    ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001
    flags = ES_CONTINUOUS | (ES_SYSTEM_REQUIRED if on else 0)
    return bool(ctypes.windll.kernel32.SetThreadExecutionState(flags))


def write_status(path, d):
    """Atomically write a small JSON status file (used to monitor from another machine/sandbox)."""
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(d), encoding="utf8")
    os.replace(tmp, path)


class Logger:
    def __init__(self, path):
        self.path, self.lock = path, threading.Lock()

    def __call__(self, msg):
        line = f"{utc()} {msg}"
        with self.lock, open(self.path, "a", encoding="utf8") as fh:
            fh.write(line + "\n")
        print(line, flush=True)


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def missing(done, total):
    gaps, pos = [], 0
    for a, b in merge(done):
        if a > pos:
            gaps.append((pos, a))
        pos = max(pos, b)
    if pos < total:
        gaps.append((pos, total))
    return gaps


class State:
    """Completed byte intervals [start, end), persisted atomically."""

    def __init__(self, path, total, fid):
        self.path, self.lock = path, threading.Lock()
        if path.exists():
            s = json.loads(path.read_text(encoding="utf8"))
            assert s["total"] == total and s["file_id"] == fid, "state file belongs to another file"
            self.done = [list(x) for x in s["done"]]
        else:
            self.done = []
        self.total, self.fid = total, fid

    def add(self, a, b):
        with self.lock:
            self.done = merge(self.done + [[a, b]])
            tmp = self.path.with_suffix(".tmp")
            tmp.write_text(json.dumps({"file_id": self.fid, "total": self.total, "done": self.done}), encoding="utf8")
            os.replace(tmp, self.path)

    def nbytes(self):
        with self.lock:
            return sum(b - a for a, b in self.done)


def fetch_chunk(sess, url, part, a, b, counter):
    """Download bytes [a, b) to offset a of part. Returns on success, raises on failure."""
    r = sess.get(url, headers={"Range": f"bytes={a}-{b - 1}"}, stream=True, timeout=(30, READ_TIMEOUT))
    try:
        if r.status_code != 206:
            raise RuntimeError(f"HTTP {r.status_code}")
        cr = r.headers.get("Content-Range", "")
        if not cr.startswith(f"bytes {a}-{b - 1}/"):
            raise RuntimeError(f"unexpected Content-Range {cr!r}")
        got = 0
        try:
            with open(part, "r+b") as fo:
                fo.seek(a)
                for blk in r.iter_content(1024 * 1024):
                    fo.write(blk)
                    got += len(blk)
                    counter.add(len(blk))
                fo.flush()
                os.fsync(fo.fileno())
            if got != b - a:
                raise RuntimeError(f"short read {got} of {b - a}")
        except Exception:
            counter.add(-got)       # bytes of a failed attempt do not count as progress
            raise
    finally:
        r.close()


class Counter:
    def __init__(self):
        self.n, self.lock = 0, threading.Lock()

    def add(self, k):
        with self.lock:
            self.n += k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file-id", type=int, required=True)
    ap.add_argument("--dest", required=True, help="absolute path of the final file")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--chunk-mib", type=float, default=32)
    ap.add_argument("--limit-bytes", type=int, default=None, help="only fetch bytes [0, limit) (speed test)")
    ap.add_argument("--status-file", default=None, help="JSON status mirrored every minute (optional)")
    args = ap.parse_args()
    status_path = Path(args.status_file) if args.status_file else None

    dest = Path(args.dest)
    part = dest.with_name(dest.name + ".part")
    statef = dest.with_name(dest.name + ".part.state.json")
    stopf = dest.with_name(dest.name + ".stop")
    log = Logger(dest.with_name(dest.name + ".download.log"))
    dest.parent.mkdir(parents=True, exist_ok=True)

    meta = [f for f in requests.get(API, timeout=60).json()["files"] if f["id"] == args.file_id][0]
    total = int(meta["size"])
    md5_pub = meta.get("computed_md5") or meta.get("supplied_md5") or ""
    url = NDL.format(fid=args.file_id)
    probe = requests.get(url, headers={"Range": "bytes=0-0"}, timeout=60)
    s3_total = int(probe.headers["Content-Range"].split("/")[-1])
    assert probe.status_code == 206 and s3_total == total, f"size mismatch: API {total}, S3 {s3_total}"

    state = State(statef, total, args.file_id)
    if not part.exists():
        free = shutil.disk_usage(dest.parent).free
        need = total - state.nbytes()
        assert args.limit_bytes or free > need + 2 * 1024**3, f"not enough free space: {free} < {need} + 2 GiB"
        part.touch()
    target = min(total, args.limit_bytes) if args.limit_bytes else total
    chunk = int(args.chunk_mib * 1024 * 1024)
    jobs = []
    for g0, g1 in missing(state.done, target):
        for a in range(g0, g1, chunk):
            jobs.append((a, min(a + chunk, g1)))
    log(f"start pid={os.getpid()} file={meta['name']} total={total} md5_published={md5_pub or 'none'} "
        f"target={target} already={state.nbytes()} chunks={len(jobs)} chunk={chunk} threads={args.threads}")

    lock = threading.Lock()
    queue = list(reversed(jobs))           # pop() hands out the lowest offset first
    failed = []
    counter = Counter()
    stop_evt = threading.Event()

    def worker():
        sess = requests.Session()
        while not stop_evt.is_set():
            with lock:
                if not queue:
                    return
                a, b = queue.pop()
            for k in range(1, MAX_TRIES_PER_CHUNK + 1):
                if stopf.exists():
                    stop_evt.set()
                    with lock:
                        queue.append((a, b))
                    return
                try:
                    fetch_chunk(sess, url, part, a, b, counter)
                    state.add(a, b)
                    break
                except Exception as e:
                    wait = min(120, 2 ** min(k, 7))
                    log(f"chunk {a}-{b} try {k} failed: {type(e).__name__}: {str(e)[:150]}; retry in {wait}s")
                    sess.close()
                    sess = requests.Session()
                    time.sleep(wait)
            else:
                with lock:
                    failed.append((a, b))

    awake = keep_awake(True)
    log(f"keep_awake={awake}")
    t0 = time.time()
    base = state.nbytes()
    st = {"file": meta["name"], "file_id": args.file_id, "dest": str(dest), "pid": os.getpid(), "total": total,
          "target": target, "started_utc": utc(), "phase": "downloading", "failed_chunks": 0}
    write_status(status_path, {**st, "updated_utc": utc(), "have": base, "last_min_MBps": None, "mean_MBps": None})
    ths = [threading.Thread(target=worker, daemon=True) for _ in range(args.threads)]
    for t in ths:
        t.start()
    last_t, last_n = t0, 0
    while any(t.is_alive() for t in ths):
        time.sleep(1)
        now = time.time()
        if now - last_t >= LOG_EVERY:
            n = counter.n
            rate, mean = (n - last_n) / (now - last_t), n / (now - t0)
            have = base + n
            eta = (target - have) / mean / 3600 if mean > 0 else float("inf")
            log(f"progress {have / 1e9:.3f}/{target / 1e9:.3f} GB ({100 * have / target:.2f}%) "
                f"last-min {rate / 1e6:.2f} MB/s mean {mean / 1e6:.2f} MB/s ETA {eta:.2f} h")
            write_status(status_path, {**st, "updated_utc": utc(), "have": have, "last_min_MBps": rate / 1e6,
                                       "mean_MBps": mean / 1e6, "eta_h": eta, "failed_chunks": len(failed)})
            last_t, last_n = now, n
    dt = time.time() - t0
    n = counter.n
    log(f"end-of-run fetched={n} bytes in {dt:.1f}s mean={n / dt / 1e6:.2f} MB/s failed_chunks={len(failed)} "
        f"stopped={stop_evt.is_set()} complete_bytes={state.nbytes()}/{total}")
    summary = {"fetched_bytes": n, "seconds": dt, "mean_MBps": n / dt / 1e6, "failed_chunks": failed,
               "stopped": stop_evt.is_set(), "complete_bytes": state.nbytes(), "total": total}
    write_status(status_path, {**st, "updated_utc": utc(), "have": state.nbytes(), "phase": "verifying"
                               if state.nbytes() == total and not failed else "stopped", **summary})
    if state.nbytes() == total and not failed:
        m, s = hashlib.md5(), hashlib.sha256()
        with open(part, "rb") as fh:
            for blk in iter(lambda: fh.read(16 * 1024 * 1024), b""):
                m.update(blk)
                s.update(blk)
        size = part.stat().st_size
        ok = size == total and (not md5_pub or m.hexdigest() == md5_pub)
        log(f"verify size={size} md5={m.hexdigest()} md5_published={md5_pub or 'none'} sha256={s.hexdigest()} ok={ok}")
        summary.update({"size": size, "md5": m.hexdigest(), "sha256": s.hexdigest(), "verified": ok})
        if ok:
            os.replace(part, dest)
            statef.unlink()
    write_status(status_path, {**st, "updated_utc": utc(), "have": state.nbytes(),
                               "phase": "done" if summary.get("verified") else ("verify_failed" if "verified" in summary
                                                                                else "stopped"), **summary})
    keep_awake(False)
    print(json.dumps(summary), flush=True)
    return 0 if not failed and summary.get("verified", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
