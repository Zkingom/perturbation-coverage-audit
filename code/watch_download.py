# watch_download.py -- periodic check of a download status file written by download_resumable.py.
# Returns (stops watching) as soon as something needs attention:
#   finished / verify failed / stopped, failed chunks, no status update for STALL_MIN minutes while this
#   watcher itself stayed awake, status never appears within WAIT_START_H hours, or MAX_H hours elapse.
# If the watcher's own sleep overran by a lot, the machine was asleep; that is logged, not treated as a stall.
import json
import time
from datetime import datetime, timezone
from pathlib import Path

INTERVAL_MIN = 15
STALL_MIN = 30
WAIT_START_H = 3
MAX_H = 14


def parse(ts):
    return datetime.strptime(ts, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc).timestamp()


def watch(status_path, log_path):
    status_path, log_path = Path(status_path), Path(log_path)
    t_start = time.time()

    def log(msg):
        line = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + msg
        with open(log_path, "a", encoding="utf8") as fh:
            fh.write(line + "\n")
        return line

    log(f"watch start; status={status_path}")
    last = None
    grace = False                     # True right after a detected machine sleep
    while True:
        t_before = time.time()
        if status_path.exists():
            s = json.loads(status_path.read_text(encoding="utf8"))
            age_min = (time.time() - parse(s["updated_utc"])) / 60
            have, tot = s.get("have", 0), s.get("target", 1)
            msg = (f"check phase={s['phase']} {have / 1e9:.2f}/{tot / 1e9:.2f} GB ({100 * have / tot:.1f}%) "
                   f"last-min={s.get('last_min_MBps')} mean={s.get('mean_MBps')} eta_h={s.get('eta_h')} "
                   f"failed={s.get('failed_chunks')} status_age_min={age_min:.1f} pid={s.get('pid')}")
            last = log(msg)
            if s["phase"] in ("done", "verify_failed", "stopped"):
                return {"reason": s["phase"], "status": s, "last": last}
            fc = s.get("failed_chunks")
            if (isinstance(fc, list) and fc) or (isinstance(fc, int) and fc > 0):
                return {"reason": "failed_chunks", "status": s, "last": last}
            if age_min > STALL_MIN and not grace:
                return {"reason": "stalled", "status": s, "last": last}
        else:
            last = log("check: no status file yet")
            if time.time() - t_start > WAIT_START_H * 3600:
                return {"reason": "never_started", "last": last}
        if time.time() - t_start > MAX_H * 3600:
            return {"reason": "max_time", "last": last}
        time.sleep(INTERVAL_MIN * 60)
        overrun = (time.time() - t_before) / 60 - INTERVAL_MIN
        grace = overrun > 5
        if grace:
            log(f"note: watcher slept {overrun:.0f} min longer than planned (machine sleep?)")
