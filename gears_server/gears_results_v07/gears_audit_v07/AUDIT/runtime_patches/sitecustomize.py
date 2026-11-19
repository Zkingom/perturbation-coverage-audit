# Runtime-only patches for the frozen GEARS bundle (no file under code/ is modified).
# Loaded automatically by Python because this directory is on PYTHONPATH (run_all.sh sets it).
#   fast_predict_patch : predict() input construction without per-perturbation AnnData slicing (same RNG calls, same inputs)
#   mem_patch          : similarity networks read from precomputed cache (GEARS's own cache-hit path) + lean evaluate()
import sys
try:
    import fast_predict_patch, mem_patch  # noqa: F401
except Exception as e:  # never silently run half-patched
    print(f"[sitecustomize] FAILED to load runtime patches: {e!r}", file=sys.stderr, flush=True)
    raise
