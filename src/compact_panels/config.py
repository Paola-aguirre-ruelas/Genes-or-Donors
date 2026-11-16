"""Versioned configuration.

Two kinds of TOML files are used:

* a *study configuration* (``configs/study.toml``) holding every scientific parameter:
  classes, eligibility rules, seeds, caps, selector and classifier settings, the
  historical fold-specific conventions and the expected real-data invariants;
* a *paths file* (``configs/paths.local.toml``, never committed) telling the pipeline
  where the verified inputs, optional historical reference tree and run directory are.

A study configuration may declare ``extends = "other.toml"``; tables are deep-merged.
Relative paths in a paths file are resolved against the paths file's own directory.
"""
from __future__ import annotations

import copy
import os
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .io_utils import IntegrityError, json_sha


def _deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_toml(path: str | os.PathLike) -> dict:
    path = Path(path).resolve()
    with open(path, "rb") as f:
        data = tomllib.load(f)
    parent = data.pop("extends", None)
    resets = data.pop("reset", [])
    if parent:
        base = load_toml(path.parent / parent)
        for key in resets:          # drop inherited tables (e.g. real-data invariants) before merging
            base.pop(key, None)
        data = _deep_merge(base, data)
    return data


REQUIRED_SECTIONS = ("meta", "data", "cohort", "partition", "withholding", "subsample", "hvg",
                     "selector", "conventions", "classifier", "evaluation", "random_reference", "stability")


@dataclass(frozen=True)
class StudyConfig:
    raw: dict
    path: Path

    def __getitem__(self, key: str) -> Any:
        return self.raw[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.raw.get(key, default)

    @property
    def sha256(self) -> str:
        return json_sha(self.raw)

    @property
    def classes(self) -> list[str]:
        return list(self.raw["cohort"]["classes"])

    @property
    def k(self) -> int:
        return int(self.raw["selector"]["k"])

    @property
    def n_folds(self) -> int:
        return int(self.raw["partition"]["n_splits"])

    @property
    def seed(self) -> int:
        return int(self.raw["partition"]["seed"])

    def convention(self, fold: int) -> dict:
        """Historical fold-specific conventions (fold 0 used a separate recovery script)."""
        conv = self.raw["conventions"]
        special = {int(f) for f in conv.get("fold0_style_folds", [])}
        base = dict(conv["default"])
        if fold in special:
            base.update(conv["fold0_style"])
        return base

    @property
    def expected(self) -> dict:
        return self.raw.get("expected", {})

    @property
    def inputs(self) -> dict:
        return self.raw.get("inputs", {})


def load_study(path: str | os.PathLike) -> StudyConfig:
    path = Path(path).resolve()
    raw = load_toml(path)
    missing = [s for s in REQUIRED_SECTIONS if s not in raw]
    if missing:
        raise IntegrityError(f"STOP: study configuration {path.name} lacks sections {missing}")
    if raw["conventions"].get("default", {}).get("search_variant") not in ("cv_halving", "fold0_recovery"):
        raise IntegrityError("STOP: unknown default search variant")
    return StudyConfig(raw=raw, path=path)


@dataclass(frozen=True)
class Paths:
    h5ad: Path
    identity_mapping: Path
    qc_exclusion: Path
    run_dir: Path
    legacy_root: Path | None
    source: Path | None

    def as_public_dict(self) -> dict:
        """Names only: absolute local paths are never written into public artifacts."""
        return {"h5ad": self.h5ad.name, "identity_mapping": self.identity_mapping.name,
                "qc_exclusion": self.qc_exclusion.name, "run_dir": self.run_dir.name,
                "legacy_root_configured": self.legacy_root is not None}


def load_paths(path: str | os.PathLike | None = None, **overrides: str | None) -> Paths:
    """Load the paths file (``--paths``, else ``$CGP_PATHS``) and apply CLI overrides."""
    data: dict = {}
    src = None
    path = path or os.environ.get("CGP_PATHS")
    if path:
        src = Path(path).resolve()
        with open(src, "rb") as f:
            data = tomllib.load(f).get("paths", {})
    for k, v in overrides.items():
        if v is not None:
            data[k] = v
    base = src.parent if src else Path.cwd()

    def resolve(key: str, required: bool = True) -> Path | None:
        v = data.get(key)
        if v in (None, ""):
            if required:
                raise IntegrityError(f"STOP: path '{key}' not configured (paths file or CLI flag)")
            return None
        p = Path(os.path.expanduser(str(v)))
        return (p if p.is_absolute() else base / p).resolve()

    return Paths(h5ad=resolve("h5ad"), identity_mapping=resolve("identity_mapping"),
                 qc_exclusion=resolve("qc_exclusion"), run_dir=resolve("run_dir"),
                 legacy_root=resolve("legacy_root", required=False), source=src)
