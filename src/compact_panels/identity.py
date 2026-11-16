"""Biological identity mapping and cohort construction.

The H5AD's ``patient_id`` denotes a *sample* label. Some biological individuals contributed
several samples under different labels, so all isolation is enforced on the published
biological individual identifier obtained through a frozen mapping table
(columns ``h5ad_patient_id`` and ``individual_id``).

Cohort steps (historical order):
1. map patient_id -> individual_id (every patient_id must be mapped);
2. remove cells with the excluded stimulation status (LPS);
3. remove exactly the recorded QC-excluded cells (every listed cell must be present);
4. restrict to the five frozen classes and count cells per individual x class over all
   remaining samples and collection days;
5. eligible individuals: >= ``min_evaluable_classes`` classes with >= ``min_cells_per_class`` cells.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

MAPPING_COLUMNS = ("h5ad_patient_id", "individual_id")


def load_identity_mapping(path) -> pd.DataFrame:
    m = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in MAPPING_COLUMNS if c not in m.columns]
    if missing:
        raise ValueError(f"identity mapping lacks columns {missing}")
    if m.h5ad_patient_id.duplicated().any():
        raise ValueError("identity mapping has duplicated h5ad_patient_id rows")
    return m


def load_cell_list(path) -> list[str]:
    with open(path, encoding="utf-8") as f:
        return f.read().split()


def mapping_pairs_sha(mapping: pd.DataFrame) -> str:
    """Content hash of the (h5ad_patient_id, individual_id) pairs, independent of extra columns."""
    import hashlib
    lines = sorted(f"{p}\t{i}" for p, i in zip(mapping.h5ad_patient_id.astype(str), mapping.individual_id.astype(str)))
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def mapping_from_sdrf(sdrf_path, patient_ids) -> pd.DataFrame:
    """Derive patient_id -> individual_id from the public ArrayExpress SDRF (E-MTAB-10026).

    ``Source Name`` identifies the sample; Cambridge sources are named ``BGCVnn_<patient_id>``,
    all others equal the H5AD patient_id. ``Characteristics[individual]`` is the published
    biological individual."""
    s = pd.read_csv(sdrf_path, sep="\t", dtype=str)
    s.columns = [c.strip() for c in s.columns]
    s = s[["Source Name", "Characteristics[individual]"]].drop_duplicates()
    key: dict[str, str] = {}
    for sn, ind in zip(s["Source Name"], s["Characteristics[individual]"]):
        pid = sn.split("_", 1)[1] if sn.startswith("BGCV") else sn
        if key.get(pid, ind) != ind:
            raise ValueError(f"SDRF assigns sample {pid} to several individuals")
        key[pid] = ind
    pids = sorted(set(patient_ids))
    missing = [p for p in pids if p not in key]
    if missing:
        raise KeyError(f"patient_id values absent from the SDRF: {missing}")
    return pd.DataFrame({"h5ad_patient_id": pids, "individual_id": [key[p] for p in pids]})


def qc_exclusion_from_rule(obs: pd.DataFrame, rule: dict) -> list[str]:
    """Cells matching the recorded QC rule (all columns must match; ``status_is_not`` negated)."""
    mask = pd.Series(True, index=obs.index)
    for col, val in rule.items():
        if col == "status_is_not":
            mask &= obs["Status"] != val
        else:
            mask &= obs[col] == val
    return sorted(obs.index[mask])


def map_individuals(obs: pd.DataFrame, mapping: pd.DataFrame, patient_column: str = "patient_id") -> pd.Series:
    """patient_id -> individual_id via the frozen mapping. Raises if any patient_id is unmapped."""
    m = dict(zip(mapping.h5ad_patient_id.astype(str), mapping.individual_id.astype(str)))
    missing = sorted(set(obs[patient_column]) - set(m))
    if missing:
        raise KeyError(f"patient_id values not in the frozen mapping: {missing}")
    return obs[patient_column].map(m)


def apply_exclusion(pool: pd.DataFrame, names) -> pd.DataFrame:
    """Remove exactly the listed cells; raises if a listed cell is not in the pool."""
    s = set(names)
    if not s <= set(pool.index):
        raise ValueError(f"{len(s - set(pool.index))} listed cells are not in the non-excluded-status pool")
    return pool[~pool.index.isin(s)]


def eligibility(pool: pd.DataFrame, classes, min_cells: int, min_classes: int, cell_type_column="author_cell_type"):
    """Counts n_{d,c} over all remaining collections (frozen classes only), k_d and the eligible list."""
    five = pool[pool[cell_type_column].isin(classes)]
    counts = pd.crosstab(five.individual_id, five[cell_type_column]).reindex(columns=list(classes), fill_value=0)
    counts = counts.reindex(sorted(pool.individual_id.unique()), fill_value=0)
    k = (counts >= min_cells).sum(axis=1)
    return counts, k, sorted(k[k >= min_classes].index)


@dataclass
class Cohort:
    obs: pd.DataFrame            # all cells, with individual_id
    nonexcluded: pd.DataFrame    # after removing the excluded status
    pool: pd.DataFrame           # after QC exclusion
    counts: pd.DataFrame
    k: pd.Series
    eligible: list[str]
    analysis: pd.DataFrame       # eligible individuals x frozen classes
    site_of: dict[str, str]
    excluded_cells: set[str]


def build_cohort(obs: pd.DataFrame, mapping: pd.DataFrame, excluded_cells: list[str], cfg) -> Cohort:
    d = cfg["data"]
    classes = cfg.classes
    obs = obs.copy()
    obs["individual_id"] = map_individuals(obs, mapping, d["patient_column"])
    nonexcluded = obs[obs[d["status_column"]] != d["excluded_status"]].copy()
    pool = apply_exclusion(nonexcluded, excluded_cells)
    c = cfg["cohort"]
    counts, k, E = eligibility(pool, classes, int(c["min_cells_per_class"]), int(c["min_evaluable_classes"]),
                               d["cell_type_column"])
    analysis = pool[pool.individual_id.isin(E) & pool[d["cell_type_column"]].isin(classes)]
    sites = pool.groupby("individual_id")[d["site_column"]].agg(lambda x: sorted(set(x)))
    multi_site = {i: s for i, s in sites.items() if len(s) != 1}
    if multi_site:
        raise ValueError(f"individuals mapped to more than one site: {multi_site}")
    site_of = {i: s[0] for i, s in sites.items()}
    return Cohort(obs=obs, nonexcluded=nonexcluded, pool=pool, counts=counts, k=k, eligible=E,
                  analysis=analysis, site_of=site_of, excluded_cells=set(excluded_cells))
