#!/usr/bin/env python3
"""Hash every file under ROOT (optionally excluding sub-directories) or compare with a snapshot.

    python scripts/snapshot_tree.py ROOT --exclude DIR [--exclude DIR] --out snapshot.json
    python scripts/snapshot_tree.py ROOT --exclude DIR --compare snapshot.json
Exit status 0 only if the comparison finds no changed, missing or added file.
"""
import argparse, hashlib, json, os, sys
from pathlib import Path


def snapshot(root: Path, exclude: set[str]) -> dict:
    out = {}
    for dp, dns, fns in os.walk(root):
        dns[:] = [d for d in dns if d not in exclude and d != "__pycache__"]
        for f in fns:
            if f == ".DS_Store":
                continue
            p = Path(dp) / f
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for b in iter(lambda: fh.read(1 << 22), b""):
                    h.update(b)
            out[str(p.relative_to(root))] = h.hexdigest()
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("root"); ap.add_argument("--exclude", action="append", default=[])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--out"); g.add_argument("--compare")
    a = ap.parse_args()
    cur = snapshot(Path(a.root).resolve(), set(a.exclude))
    if a.out:
        Path(a.out).write_text(json.dumps(cur, indent=0, sort_keys=True))
        print(f"{len(cur)} files hashed")
        return 0
    old = json.loads(Path(a.compare).read_text())
    changed = sorted(k for k in old if k in cur and cur[k] != old[k])
    missing = sorted(set(old) - set(cur)); added = sorted(set(cur) - set(old))
    print(json.dumps({"files_in_snapshot": len(old), "unchanged": len(old) - len(changed) - len(missing),
                      "changed": changed, "missing": missing, "added": added}, indent=1))
    return 0 if not (changed or missing or added) else 1


if __name__ == "__main__":
    sys.exit(main())
