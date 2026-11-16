"""Hashing and atomic, write-once file helpers.

Every artifact written by the pipeline goes through these helpers so that an
interrupted process never leaves a half-written file under its final name, and
completed artifacts are never silently overwritten.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
from pathlib import Path
from typing import Any, Iterable

import numpy as np


class IntegrityError(RuntimeError):
    """Raised when a verification step fails. The pipeline stops; nothing is repaired."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise IntegrityError("STOP: " + message)


def sha256_file(path: str | os.PathLike, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_lines(lines: Iterable[str]) -> str:
    """Hash of a newline-terminated list file (the historical cell-list convention)."""
    return hashlib.sha256(("\n".join(lines) + "\n").encode("utf-8")).hexdigest()


def json_sha(obj: Any) -> str:
    """Historical provenance hash: json.dumps(obj, sort_keys=True) with default separators."""
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def _atomic_write(path: Path, data: bytes, overwrite: bool) -> None:
    path = Path(path)
    if path.exists() and not overwrite:
        raise IntegrityError(f"STOP: refusing to overwrite existing artifact {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def write_bytes(path: str | os.PathLike, data: bytes, *, overwrite: bool = False) -> None:
    _atomic_write(Path(path), data, overwrite)


def json_bytes(obj: Any) -> bytes:
    return (json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False, default=_default) + "\n").encode("utf-8")


def _default(o: Any):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"not JSON serializable: {type(o)}")


def write_json(path: str | os.PathLike, obj: Any, *, overwrite: bool = False) -> None:
    _atomic_write(Path(path), json_bytes(obj), overwrite)


def read_json(path: str | os.PathLike) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_npz(path: str | os.PathLike, *, overwrite: bool = False, **arrays: np.ndarray) -> None:
    buf = io.BytesIO()
    np.savez_compressed(buf, **arrays)
    _atomic_write(Path(path), buf.getvalue(), overwrite)


def write_lines(path: str | os.PathLike, lines: Iterable[str], *, overwrite: bool = False) -> None:
    _atomic_write(Path(path), ("\n".join(lines) + "\n").encode("utf-8"), overwrite)
