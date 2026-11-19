# vcs_freeze_tags.py -- git history and freeze tags for the project (dulwich; no git binary on this machine).
# Repository: <ROOT>/_vcs_repo.git, a BARE repository (no working tree), so that no git operation can ever write to,
# delete or check out files in the project directory. Clone it with `git clone F:\codex_in_this\claude_science\_vcs_repo.git`
# or inspect it with `git --git-dir=_vcs_repo.git --work-tree=. status`.
#   python vcs_freeze_tags.py init       create the repository and the retroactive freeze commits + annotated tags
#   python vcs_freeze_tags.py snapshot "message"   commit the current tracked tree on top of main
#   python vcs_freeze_tags.py tag <name> <manifest>  commit the files of a new freeze manifest and tag it
# Retroactive tags: each tag points to a commit whose tree contains the freeze manifest and every tracked file whose
# CURRENT sha256 is recorded in that manifest. Files whose recorded hash differs from the current file (superseded by a
# later frozen version) cannot be reconstructed and are listed in the tag message. Commit times are the actual time
# of tagging; the freeze time is stated in the tag message.
import hashlib
import json
import os
import stat
import sys
import time
from pathlib import Path

from dulwich.objects import Blob, Commit, Tag, Tree
from dulwich.repo import Repo

ROOT = Path(__file__).resolve().parent.parent
VCS = ROOT / "_vcs_repo.git"
IDENT = b"pert-audit project <project@localhost>"
TRACK_DIRS = ("code", "scripts", "reports", "manifests", "manuscript_v1.0", "manuscript_v1.1", "review_packet",
              "gears_bundle_v07/code")
TRACK_TOP = ("handover.md",)
TRACK_EXT = {".py", ".R", ".r", ".sh", ".ps1", ".md", ".json", ".csv", ".txt", ".tex", ".png", ".pdf", ".toml"}
EXCLUDE_PARTS = {"fonts", "_tmp", "_tectonic_cache", "_home", "__pycache__", "cache", "ext_cache", "archive_unused"}
MAX_BYTES = 20 * 1024 * 1024
FREEZES = [  # (tag, manifest)
    ("v0.5", "freeze_v0.5.json"), ("v0.5.1", "freeze_v0.5.1.json"), ("v0.6", "freeze_v0.6.json"),
    ("v0.6.1", "freeze_v0.6.1.json"), ("v0.6.2", "freeze_v0.6.2.json"), ("v0.7", "freeze_v0.7.json"),
    ("gears_coupling-v1.0", "freeze_gears_coupling_v1.0.json"), ("gears_coupling-v1.1", "freeze_gears_coupling_v1.1.json"),
    ("v1.2", "freeze_v1.2.json"), ("v1.2.1", "freeze_v1.2.1.json"), ("v1.3", "freeze_v1.3.json"), ("v1.4", "freeze_v1.4.json")]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def tracked_files():
    out = [ROOT / f for f in TRACK_TOP if (ROOT / f).is_file()]
    for d in TRACK_DIRS:
        for p in (ROOT / d).rglob("*"):
            if p.is_file() and p.suffix in TRACK_EXT and not (set(p.relative_to(ROOT).parts) & EXCLUDE_PARTS) \
                    and p.stat().st_size <= MAX_BYTES:
                out.append(p)
    return sorted(set(out))


def _tree_dict(repo, tid):
    d = {}
    for e in repo[tid].items():
        d[e.path.decode()] = _tree_dict(repo, e.sha) if e.mode == stat.S_IFDIR else e.sha
    return d


def _build(repo, d):
    t = Tree()
    for name, v in d.items():
        if isinstance(v, dict):
            t.add(name.encode(), stat.S_IFDIR, _build(repo, v).id)
        else:
            t.add(name.encode(), 0o100644, v)
    repo.object_store.add_object(t)
    return t


def commit(repo, files, message, replace=False):
    try:
        parent = repo.refs[b"refs/heads/main"]
        d = {} if replace else _tree_dict(repo, repo[parent].tree)
        parents = [parent]
    except KeyError:
        d, parents = {}, []
    for p in files:
        parts = p.relative_to(ROOT).as_posix().split("/")
        b = Blob.from_string(p.read_bytes())
        repo.object_store.add_object(b)
        cur = d
        for q in parts[:-1]:
            cur = cur.setdefault(q, {})
        cur[parts[-1]] = b.id
    c = Commit()
    c.tree = _build(repo, d).id
    c.author = c.committer = IDENT
    c.commit_time = c.author_time = int(time.time())
    c.commit_timezone = c.author_timezone = 0
    c.encoding = b"UTF-8"
    c.message = message.encode("utf8")
    c.parents = parents
    repo.object_store.add_object(c)
    repo.refs[b"refs/heads/main"] = c.id
    return c.id


def annotated_tag(repo, name, target, message):
    ref = b"refs/tags/" + name.encode()
    if ref in repo.refs:
        raise SystemExit(f"tag {name} already exists; refusing to move it")
    t = Tag()
    t.tagger = IDENT
    t.message = message.encode("utf8")
    t.name = name.encode()
    t.object = (Commit, target)
    t.tag_time = int(time.time())
    t.tag_timezone = 0
    repo.object_store.add_object(t)
    repo.refs[ref] = t.id
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t.tag_time))


def freeze_commit_and_tag(repo, tag, manifest, cur_hash):
    mp = ROOT / "manifests" / manifest
    text = mp.read_text(encoding="utf8")
    j = json.loads(text)
    rec = set(__import__("re").findall(r"[0-9a-f]{64}", text))
    files = [mp] + [p for p, h in cur_hash.items() if h in rec]
    frozen = j.get("frozen_utc", "unknown")
    names_rec = set()

    def walk(o, key=""):
        if isinstance(o, dict):
            for k, v in o.items():
                if isinstance(v, str) and len(v) == 64:
                    names_rec.add((k, v))
                elif isinstance(v, dict) and isinstance(v.get("sha256"), str):
                    names_rec.add((v.get("path") or k, v["sha256"]))
                walk(v, k)
    walk(j)
    have = {h for h in cur_hash.values()}
    missing = sorted({Path(n).name for n, h in names_rec if h not in have})
    cid = commit(repo, files, f"Freeze {tag} ({manifest}); frozen {frozen}. Retroactive commit of hash-verified files.")
    msg = (f"Freeze {tag}\nmanifest: manifests/{manifest}\nfrozen_utc: {frozen}\n"
           f"retroactive: tag created {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}; the commit contains the files "
           f"whose current sha256 equals a hash recorded in the manifest ({len(files) - 1} files + manifest).\n"
           f"recorded but not included (current file differs or file is outside the tracked tree, e.g. raw data): "
           f"{', '.join(missing) if missing else 'none'}\n")
    when = annotated_tag(repo, tag, cid, msg)
    return {"tag": tag, "manifest": manifest, "frozen_utc": frozen, "tag_utc": when, "commit": cid.decode(),
            "n_files": len(files) - 1, "not_included": missing}


def main():
    mode = sys.argv[1]
    if mode == "init":
        if VCS.exists():
            raise SystemExit(f"{VCS} already exists; refusing to re-initialise")
        VCS.mkdir()
        repo = Repo.init_bare(str(VCS))
        files = tracked_files()
        cur = {p: sha256(p) for p in files}
        log = [freeze_commit_and_tag(repo, t, m, cur) for t, m in FREEZES]
        snap = commit(repo, files, "Snapshot of the tracked project tree at repository creation", replace=False)
        out = {"repo": str(VCS), "tags": log, "snapshot_commit": snap.decode(), "n_tracked": len(files)}
        (VCS / "freeze_tags_init.json").write_text(json.dumps(out, indent=1), encoding="utf8")
        print(json.dumps(out, indent=1))
    elif mode == "snapshot":
        repo = Repo(str(VCS))
        print(commit(repo, tracked_files(), sys.argv[2]).decode())
    elif mode == "tag":
        repo = Repo(str(VCS))
        files = tracked_files()
        cur = {p: sha256(p) for p in files}
        print(json.dumps(freeze_commit_and_tag(repo, sys.argv[2], sys.argv[3], cur)))


if __name__ == "__main__":
    main()
