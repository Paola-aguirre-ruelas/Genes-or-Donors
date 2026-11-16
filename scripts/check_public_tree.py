#!/usr/bin/env python3
"""Static checks for files that would be published (git-tracked or untracked-but-not-ignored).

Fails on: expression data or archives, files larger than 5 MB, credential-like strings,
absolute home-directory paths, e-mail addresses, and terms from an optional private deny-list.

The deny-list is never stored in the repository. It is read, one term per line ('#' starts a
comment), from ``--deny-file PATH``, from the file named by ``$CGP_PRIVATE_DENY_FILE``, or from
``.private-deny-terms`` in the repository root (git-ignored). With ``--history`` the same checks are
applied to every blob, path and commit message reachable in the Git history.

Run from the repository root:  python scripts/check_public_tree.py . [--history] [--deny-file PATH]
"""
import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

DATA_SUFFIXES = {".h5ad", ".h5", ".loom", ".mtx", ".rds", ".zarr", ".npz", ".zip", ".tar", ".gz"}
BINARY_SUFFIXES = {".png", ".pdf", ".jpg", ".jpeg"}
MAX_BYTES = 5 * 1024 * 1024
PATTERNS = {
    "home path": re.compile(r"/(Users|home)/[A-Za-z0-9._-]+/"),
    "e-mail": re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}"),
    "credential": re.compile(r"(api[_-]?key|secret|token|password)\s*[:=]\s*['\"][^'\"]{8,}", re.I),
    "private key": re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
}
# GitHub no-reply addresses, the AI co-author trailer and the approved public research contact.
ALLOWED_EMAIL = re.compile(r"^(noreply@anthropic\.com|[0-9]+\+[A-Za-z0-9-]+@users\.noreply\.github\.com"
                           r"|paolajanetaguirre@gmail\.com)$|\.invalid$")
DEFAULT_DENY_FILE = ".private-deny-terms"
LOCAL_ONLY = {DEFAULT_DENY_FILE, "historical_dirs.local.toml", "paths.local.toml"}


def git(root: Path, *args: str) -> bytes:
    return subprocess.run(["git", *args], cwd=root, capture_output=True, check=True).stdout


def tracked_files(root: Path) -> list[Path]:
    try:
        out = git(root, "ls-files", "--cached", "--others", "--exclude-standard").decode().split("\n")
        return [root / f for f in out if f]
    except (subprocess.CalledProcessError, FileNotFoundError):
        ignore = {".git", ".venv", "runs", "data", "__pycache__", "build", "dist"}
        return [p for p in root.rglob("*") if p.is_file() and not set(p.relative_to(root).parts) & ignore]


def load_deny(root: Path, explicit: str | None) -> tuple[list[bytes], str]:
    candidate = explicit or os.environ.get("CGP_PRIVATE_DENY_FILE") or str(root / DEFAULT_DENY_FILE)
    path = Path(candidate)
    if not path.is_file():
        return [], "no private deny-list supplied (generic checks only)"
    terms = [ln.split("#", 1)[0].strip() for ln in path.read_text(encoding="utf-8").splitlines()]
    return [t.encode() for t in terms if t], f"private deny-list: {sum(1 for t in terms if t)} terms"


def scan_text(raw: bytes, where: str, deny: list[bytes], problems: list[str]) -> None:
    for i, token in enumerate(deny):
        escaped = "".join(f"\\x{b:02x}" for b in token).encode()   # also catch \xNN-escaped spellings
        if (re.search(rb"(?<![A-Za-z])" + re.escape(token) + rb"(?![a-z])", raw)
                or re.search(re.escape(escaped), raw, re.I)):
            problems.append(f"{where}: contains private deny-list term #{i + 1}")
    text = raw.decode("utf-8", errors="ignore")
    for name, rx in PATTERNS.items():
        for m in rx.finditer(text):
            if name == "e-mail" and ALLOWED_EMAIL.search(m.group(0)):
                continue
            problems.append(f"{where}: {name}: {m.group(0)[:60]}")


def check_tree(root: Path, deny: list[bytes], problems: list[str]) -> int:
    files = tracked_files(root)
    for p in files:
        rel = p.relative_to(root)
        if rel.name in LOCAL_ONLY:
            problems.append(f"{rel}: local-only file would be published (must stay git-ignored)")
            continue
        if p.suffix.lower() in DATA_SUFFIXES:
            problems.append(f"{rel}: data/archive file type")
        if p.stat().st_size > MAX_BYTES:
            problems.append(f"{rel}: larger than 5 MB")
        scan_text(str(rel).encode(), f"path {rel}", deny, problems)
        if p.suffix.lower() not in BINARY_SUFFIXES:
            scan_text(p.read_bytes(), str(rel), deny, problems)
    return len(files)


def check_history(root: Path, deny: list[bytes], problems: list[str]) -> int:
    objects = git(root, "rev-list", "--all", "--objects").decode().splitlines()
    for line in objects:
        sha, _, path = line.partition(" ")
        kind = git(root, "cat-file", "-t", sha).decode().strip()
        label = f"history {sha[:10]} {path}".rstrip()
        if path:
            scan_text(path.encode(), f"{label} (path)", deny, problems)
            if Path(path).name in LOCAL_ONLY:
                problems.append(f"{label}: local-only file committed")
            if Path(path).suffix.lower() in DATA_SUFFIXES:
                problems.append(f"{label}: data/archive file type")
        if kind == "commit" or (kind == "blob" and Path(path).suffix.lower() not in BINARY_SUFFIXES):
            scan_text(git(root, "cat-file", "-p", sha), label, deny, problems)
    return len(objects)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root", nargs="?", default=".")
    ap.add_argument("--history", action="store_true", help="also scan every object reachable in the Git history")
    ap.add_argument("--deny-file", help="private deny-list (one term per line; never commit it)")
    a = ap.parse_args()
    root = Path(a.root).resolve()
    deny, deny_note = load_deny(root, a.deny_file)
    problems: list[str] = []
    n_files = check_tree(root, deny, problems)
    note = f"{n_files} files"
    if a.history:
        note += f", {check_history(root, deny, problems)} history objects"
    for line in problems:
        print("FAIL", line)
    print(f"public-tree check: {'PASS' if not problems else 'FAIL'} ({len(problems)} problems; {note}; {deny_note})")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
