# safe_io.py -- create-only file writes (engineering rule from 2026-11-04, after an accidental overwrite of a
# frozen file). create_new() refuses to write if the path already exists; use it for every new file.
from pathlib import Path


class PathExistsError(FileExistsError):
    pass


def create_new(path, data, encoding="utf8"):
    """Write a NEW file. Raises PathExistsError if anything already exists at `path` (never overwrites)."""
    p = Path(path)
    if p.exists():
        raise PathExistsError(f"refusing to overwrite existing file: {p}")
    p.parent.mkdir(parents=True, exist_ok=True)
    mode = "xb" if isinstance(data, (bytes, bytearray)) else "x"
    with open(p, mode, **({} if mode == "xb" else {"encoding": encoding})) as fh:
        fh.write(data)
    return p
