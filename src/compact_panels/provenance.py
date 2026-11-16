"""Software, platform and resource provenance."""
from __future__ import annotations

import datetime as dt
import importlib.metadata as md
import platform
import resource
import sys

PACKAGES = ("numpy", "scipy", "pandas", "scikit-learn", "anndata", "h5py", "scanpy", "numba", "threadpoolctl", "matplotlib")


def maxrss() -> int:
    """Peak resident set size of this process in bytes (macOS reports bytes, Linux KiB)."""
    r = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return r if sys.platform == "darwin" else r * 1024


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def versions() -> dict:
    out = {"python": sys.version.split()[0]}
    for p in PACKAGES:
        try:
            out[p] = md.version(p)
        except md.PackageNotFoundError:
            out[p] = None
    try:
        out["compact_panels"] = md.version("compact-gene-panels")
    except md.PackageNotFoundError:
        from . import __version__
        out["compact_panels"] = __version__
    return out


def thread_pools() -> list[dict]:
    try:
        from threadpoolctl import threadpool_info
        return [{k: v for k, v in d.items() if k in ("user_api", "internal_api", "num_threads", "version", "threading_layer")}
                for d in threadpool_info()]
    except Exception:  # noqa: BLE001
        return []


def environment() -> dict:
    return {"versions": versions(), "platform": platform.platform(), "machine": platform.machine(),
            "thread_pools": thread_pools()}


def check_pinned(expected: dict) -> list[str]:
    """Return mismatches between installed and pinned versions (empty list = all pinned)."""
    got = versions()
    return [f"{k}: expected {v}, found {got.get(k)}" for k, v in expected.items() if got.get(k) != v]
