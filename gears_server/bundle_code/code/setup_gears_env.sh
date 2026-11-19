#!/usr/bin/env bash
# setup_gears_env.sh -- one-shot install + self-check for the GEARS robustness run (Linux, NVIDIA GPU).
#   bash setup_gears_env.sh [PIP_INDEX_URL]
# Creates ./venv_gears (Python 3.10) next to this bundle. torch 2.3.0 + CUDA 12.1 (as in gears_env2.yml). If the image
# already has exactly torch 2.3.0 with CUDA 12.1, it is reused (--system-site-packages); otherwise the cu121 wheel is
# installed. Writes logs/setup_selfcheck.json. Re-running is safe.
set -euo pipefail
B="$(cd "$(dirname "$0")/.." && pwd)"; cd "$B"; mkdir -p logs
IDX="${1:-https://pypi.tuna.tsinghua.edu.cn/simple}"
PY="$(command -v python3.10 || true)"
if [ -z "$PY" ]; then
  if command -v conda >/dev/null; then conda create -y -p "$B/py310" python=3.10.14 >/dev/null; PY="$B/py310/bin/python"; else echo "need python3.10 or conda"; exit 2; fi
fi
HAS_TORCH=$($PY -c "import torch;print(torch.__version__.split('+')[0]=='2.3.0' and torch.version.cuda=='12.1')" 2>/dev/null || echo False)
if [ "$HAS_TORCH" = "True" ]; then $PY -m venv --system-site-packages venv_gears; else $PY -m venv venv_gears; fi
source venv_gears/bin/activate
python -m pip install -q --upgrade "pip<25" -i "$IDX"
if ! python -c "import torch,sys;sys.exit(0 if torch.__version__.split('+')[0]=='2.3.0' else 1)" 2>/dev/null; then
  pip install -q torch==2.3.0 --index-url https://download.pytorch.org/whl/cu121 || pip install -q torch==2.3.0 -i "$IDX"
fi
pip install -q -r code/requirements_gears.txt -i "$IDX"
python - <<'EOF'
import json, platform, subprocess, importlib.metadata as md, torch, numpy as np
out = {"python": platform.python_version(), "torch": torch.__version__, "cuda": torch.version.cuda,
       "cuda_available": torch.cuda.is_available(),
       "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
       "gpu_mem_gb": round(torch.cuda.get_device_properties(0).total_memory / 2**30, 1) if torch.cuda.is_available() else None,
       "packages": {p: md.version(p) for p in ("cell-gears", "torch-geometric", "numpy", "pandas", "scanpy", "anndata", "scipy", "scikit-learn", "networkx", "dcor")}}
import gears.version; assert gears.version.__version__ == "0.1.2"
from gears import PertData, GEARS  # noqa: F401
x = torch.randn(256, 256, device="cuda"); out["gpu_matmul_ok"] = bool(torch.isfinite(x @ x).all())
try: out["host_mem_gb"] = round(int(open("/proc/meminfo").read().split()[1]) / 2**20, 1)
except Exception: pass
out["cpu"] = subprocess.run(["bash", "-c", "lscpu | grep 'Model name' | head -1"], capture_output=True, text=True).stdout.strip()
out["nproc"] = int(subprocess.run(["nproc"], capture_output=True, text=True).stdout)
expected = {"cell-gears": "0.1.2", "torch-geometric": "2.5.3", "numpy": "1.26.4", "pandas": "2.2.2", "scanpy": "1.10.1", "anndata": "0.10.7"}
out["pins_ok"] = all(out["packages"][k] == v for k, v in expected.items()) and out["torch"].startswith("2.3.0")
out["pass"] = bool(out["cuda_available"] and out["gpu_matmul_ok"] and out["pins_ok"])
json.dump(out, open("logs/setup_selfcheck.json", "w"), indent=1); print(json.dumps(out, indent=1))
assert out["pass"], "self-check failed"
EOF
cd "$B" && sha256sum -c SHA256SUMS --quiet && echo "bundle SHA256 OK"
echo "SETUP PASS"
