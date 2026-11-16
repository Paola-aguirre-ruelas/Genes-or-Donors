"""Donor-disjoint folds, within-training withholding and training subsamples.

Historical conventions reproduced exactly:

* folds: ``StratifiedKFold(n_splits, shuffle=True, random_state=seed)`` over the
  lexicographically sorted eligible individuals, stratified by site only; if a fold's
  training individuals lack a frozen class, one redraw with the next unused seed;
* Stage B withholding for test fold k: ``default_rng([seed, k])``; for each training
  individual (lexicographic) and class (frozen order) with n >= 1 cells, draw
  m = floor(0.2 n + 1/2) cells without replacement from that stratum's cells sorted by
  name, calling ``rng.choice`` for every stratum (also when m = 0);
* subsamples from the fitting cells with independent generators ``[seed, k, 1]``
  (selector, cap 50) and ``[seed, k, 2]`` (classifier, cap 200), same iteration rules.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from fractions import Fraction

import numpy as np
import pandas as pd

from .io_utils import json_sha, sha256_lines


# ------------------------------------------------------------------ folds
def make_folds(individuals: list[str], sites: list[str], n_splits: int, shuffle: bool, seed: int):
    from sklearn.model_selection import StratifiedKFold
    if list(individuals) != sorted(individuals) or len(set(individuals)) != len(individuals):
        raise ValueError("individuals must be unique and sorted (lexicographic Python str order)")
    skf = StratifiedKFold(n_splits=n_splits, shuffle=shuffle, random_state=seed)
    return list(skf.split(np.zeros(len(individuals)), np.asarray(sites)))


def canonical_assignment(fold_individuals: list[list[str]], seed: int, n_splits: int) -> tuple[str, str]:
    payload = {"seed": seed, "n_splits": n_splits, "unit": "individual_id",
               "test_folds": [sorted(f) for f in fold_individuals]}
    s = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return s, hashlib.sha256(s.encode()).hexdigest()


@dataclass
class FoldAssignment:
    folds: list[list[str]]       # test individuals of each fold, sorted
    seed_used: int
    redraws: list[dict]
    canonical: str
    canonical_sha256: str


def assign_folds(cohort, cfg) -> FoldAssignment:
    p = cfg["partition"]
    E = cohort.eligible
    sites = [cohort.site_of[d] for d in E]
    classes = set(cfg.classes)
    ctype = cfg["data"]["cell_type_column"]
    seed = int(p["seed"])
    used, redraws = seed, []
    while True:
        folds = make_folds(E, sites, int(p["n_splits"]), bool(p["shuffle"]), used)
        cover = [set(cohort.analysis[cohort.analysis.individual_id.isin([E[j] for j in tr])][ctype]) >= classes
                 for tr, _ in folds]
        if all(cover) or used != seed:
            break
        redraws.append({"seed": used, "class_coverage_per_fold": cover})
        used = int(p["redraw_seed"])
    if not all(cover):
        raise ValueError("a fold's training individuals lack a frozen class even after the redraw")
    fold_inds = [sorted(E[j] for j in te) for _, te in folds]
    canon, sha = canonical_assignment(fold_inds, used, int(p["n_splits"]))
    return FoldAssignment(fold_inds, used, redraws, canon, sha)


# ------------------------------------------------------------------ stratified draws
def n_withheld(n: int, fraction: str | float = "0.2") -> int:
    """floor(fraction * n + 1/2) in exact rational arithmetic; for 0.2 equals (2n + 5) // 10."""
    if n < 0:
        raise ValueError("negative stratum size")
    f = Fraction(str(fraction))
    return int((f * n + Fraction(1, 2)).__floor__())


def stratified_draw(cells: pd.DataFrame, classes, seed, m_of_n, cell_type_column="author_cell_type") -> tuple[list[str], list[dict]]:
    """One generator consumed sequentially over individuals (lexicographic) x classes (frozen order)."""
    rng = np.random.default_rng(seed)
    out, rows = [], []
    for d in sorted(cells.individual_id.unique()):
        sub = cells[cells.individual_id == d]
        for c in classes:
            names = sorted(sub.index[sub[cell_type_column] == c])
            n = len(names)
            if n == 0:
                continue
            m = m_of_n(n)
            idx = rng.choice(n, size=m, replace=False)
            out += [names[i] for i in idx]
            rows.append({"individual_id": d, "class": c, "n": n, "m": m})
    return sorted(out), rows


def withhold(train: pd.DataFrame, classes, seed, fraction="0.2", cell_type_column="author_cell_type"):
    return stratified_draw(train, classes, seed, lambda n: n_withheld(n, fraction), cell_type_column)


def draw_subsample(cells: pd.DataFrame, classes, cap: int, seed, cell_type_column="author_cell_type"):
    return stratified_draw(cells, classes, seed, lambda n: min(cap, n), cell_type_column)


# ------------------------------------------------------------------ fold plan
@dataclass
class FoldPlan:
    fold: int
    train: list[str]
    test: list[str]
    fitting: list[str]
    withheld: list[str]
    selector: list[str]
    classifier: list[str]
    train_people: list[str]
    test_people: list[str]
    seeds: dict = field(default_factory=dict)

    LISTS = ("train", "test", "fitting", "withheld", "selector", "classifier", "train_people", "test_people")

    def provenance(self) -> dict:
        return {k: {"size": len(getattr(self, k)), "sha256": sha256_lines(getattr(self, k))} for k in self.LISTS}

    def provenance_sha(self) -> str:
        return json_sha(self.provenance())


def build_fold_plan(cohort, assignment: FoldAssignment, fold: int, cfg) -> FoldPlan:
    d = cfg["data"]
    ctype = d["cell_type_column"]
    classes = cfg.classes
    seed = int(cfg["partition"]["seed"]) if not assignment.redraws else assignment.seed_used
    analysis = cohort.analysis
    test_people = sorted(assignment.folds[fold])
    train_people = sorted(set().union(*(set(f) for i, f in enumerate(assignment.folds) if i != fold)))
    if set(train_people) & set(test_people):
        raise ValueError("biological identity leakage between training and test individuals")
    test = sorted(analysis.index[analysis.individual_id.isin(test_people)])
    train = sorted(set(analysis.index) - set(test))
    obs = cohort.obs
    tm, sm = obs.loc[train], obs.loc[test]
    if set(tm.individual_id) != set(train_people) or set(sm.individual_id) != set(test_people):
        raise ValueError("cell-level individual_id does not match the individual-level partition")
    if set(tm[d["patient_column"]]) & set(sm[d["patient_column"]]):
        raise ValueError("a sample label (patient_id) occurs on both sides")
    for name, df in (("train", tm), ("test", sm)):
        if (df[d["status_column"]] == d["excluded_status"]).any():
            raise ValueError(f"excluded-status cells in {name}")
        if not set(df[ctype]) <= set(classes):
            raise ValueError(f"non-frozen classes in {name}")
        if set(df.index) & cohort.excluded_cells:
            raise ValueError(f"QC-excluded cells in {name}")
    stage_b_seed = [seed, fold]
    withheld, _ = withhold(tm, classes, stage_b_seed, cfg["withholding"]["fraction"], ctype)
    fitting = sorted(set(train) - set(withheld))
    if len(fitting) + len(withheld) != len(train) or set(fitting) & set(withheld):
        raise ValueError("Stage B partition failed")
    fm = obs.loc[fitting]
    if set(fm.individual_id) != set(train_people):
        raise ValueError("a training individual lacks fitting cells")
    s = cfg["subsample"]
    sel_seed = [seed, fold, int(s["selector_stream"])]
    clf_seed = [seed, fold, int(s["classifier_stream"])]
    selector, _ = draw_subsample(fm, classes, int(s["selector_cap"]), sel_seed, ctype)
    classifier, _ = draw_subsample(fm, classes, int(s["classifier_cap"]), clf_seed, ctype)
    if not (set(selector) <= set(fitting) and set(classifier) <= set(fitting)):
        raise ValueError("subsample leakage outside the fitting cells")
    for name, lst in (("selector", selector), ("classifier", classifier)):
        if set(fm.loc[lst][ctype]) != set(classes):
            raise ValueError(f"{name} subsample lacks a frozen class")
    return FoldPlan(fold, train, test, fitting, withheld, selector, classifier, train_people, test_people,
                    seeds={"stage_b": stage_b_seed, "selector": sel_seed, "classifier": clf_seed})
