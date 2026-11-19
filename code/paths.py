"""Central path resolution. No drive letters appear anywhere else in the code.

ROOT  = directory containing config.toml (override: env PERTAUDIT_ROOT)
DL    = download_root from config.toml   (override: env PERTAUDIT_DL)
"""
import os
import tomllib
from pathlib import Path

ROOT = Path(os.environ.get("PERTAUDIT_ROOT", Path(__file__).resolve().parents[1]))
with open(ROOT / "config.toml", "rb") as _fh:
    CONFIG = tomllib.load(_fh)
if CONFIG.get("project_root", ".") != ".":
    ROOT = (ROOT / CONFIG["project_root"]).resolve()
DL = Path(os.environ.get("PERTAUDIT_DL", CONFIG["download_root"]))
MASTER_SEED = int(CONFIG["master_seed"])
PROTOCOL_VERSION = str(CONFIG["protocol_version"])

RAW_MANIFEST = ROOT / "data" / "raw" / "raw_manifest.json"
PROCESSED = ROOT / "data" / "processed"
SPLITS = ROOT / "splits"
RESULTS = ROOT / "results"
REPORTS = ROOT / "reports"
MANIFESTS = ROOT / "manifests"


def raw(name):
    """Absolute path of a raw file registered in data/raw/raw_manifest.json (paths stored relative to DL)."""
    import json
    entry = json.loads(RAW_MANIFEST.read_text(encoding="utf8"))["files"][name]
    return DL / entry["relpath"]
