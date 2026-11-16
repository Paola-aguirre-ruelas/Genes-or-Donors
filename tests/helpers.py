"""Shared fixtures for the synthetic test-suite (no real data is used anywhere)."""
from __future__ import annotations

import functools
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from compact_panels.config import load_paths, load_study

REPO = Path(__file__).resolve().parents[1]
SYNTH_CONFIG = REPO / "configs" / "synthetic.toml"
STUDY_CONFIG = REPO / "configs" / "study.toml"


def synth_cfg():
    return load_study(SYNTH_CONFIG)


def study_cfg():
    return load_study(STUDY_CONFIG)


@functools.lru_cache(maxsize=1)
def synthetic_dataset() -> dict:
    """Generate the synthetic dataset once per test process."""
    from compact_panels.synthetic import make_synthetic
    tmp = Path(tempfile.mkdtemp(prefix="cgp_synth_"))
    info = make_synthetic(tmp / "data")
    info["root"] = tmp
    return info


def synthetic_paths(run_name: str):
    info = synthetic_dataset()
    run = info["root"] / run_name
    return load_paths(None, h5ad=str(info["h5ad"]), identity_mapping=str(info["identity_mapping"]),
                      qc_exclusion=str(info["qc_exclusion"]), run_dir=str(run))


def toy_obs(n_ind=12, classes=("A", "B", "C"), per=25, sites=("S1", "S2", "S3"), seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_ind):
        for c in classes:
            for j in range(per + int(rng.integers(0, 5))):
                rows.append({"name": f"c{i:02d}_{c}_{j:03d}", "individual_id": f"I{i:02d}", "author_cell_type": c,
                             "Site": sites[i % len(sites)], "patient_id": f"P{i:02d}", "Status": "Healthy"})
    return pd.DataFrame(rows).set_index("name")
