"""Read-only input validation and dry run. Performs no model fitting and no held-out read.

Checks, in order:
  C1  configuration loads; pinned software versions present
  C2  input files exist; SHA-256 equal to the configuration pins
  C3  H5AD structure: CSR X of the expected shape/dtype; required obs columns; guard self-test
  C4  cohort reconstruction equals the expected invariants (cells, individuals, sites, classes)
  C5  fold assignment: canonical SHA-256, sizes, site balance, no redraw
  C6  per-fold isolation: individuals, sample labels and cells disjoint; withheld/fitting
      partition of training; subsamples inside fitting cells
  C7  optional historical comparison: every referenced artifact SHA-verified; fold plans equal
  C8  run directory isolated from the historical tree and writable
  C9  (--deep) training-only HVG / standardization replay against the historical checkpoints
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from . import data, preprocessing
from .io_utils import require, write_json
from .pipeline import Logger, Study
from .provenance import check_pinned, environment, maxrss, utc_now
from .reference import EXACT, Legacy, compare_plan, compare_preprocessing


class Checks:
    def __init__(self, log):
        self.items: list[dict] = []
        self.log = log

    def add(self, cid: str, desc: str, passed: bool, **detail) -> bool:
        self.items.append({"id": cid, "check": desc, "passed": bool(passed), **detail})
        self.log(f"[{'PASS' if passed else 'FAIL'}] {cid} {desc}" + ("" if passed else f" {detail}"))
        return bool(passed)

    @property
    def ok(self) -> bool:
        return all(i["passed"] for i in self.items)


def _inside(child: Path, parent: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def run_validation(cfg, paths, *, refs_path: Path | None, acceptance: dict | None, deep: bool = False,
                   write_report: bool = True) -> dict:
    log = Logger(None)
    ck = Checks(log)
    exp = cfg.expected
    # C1
    pins = cfg.get("software", {})
    mism = check_pinned(pins) if pins else []
    ck.add("C1", "configuration loaded and pinned software versions installed", not mism, config_sha256=cfg.sha256,
           mismatches=mism)
    # C2
    try:
        st = Study(cfg, paths, log, verify_inputs=True)
        ck.add("C2", "input files present and SHA-256 equal to the configuration pins", True, sha256=st.input_hashes)
    except Exception as e:  # noqa: BLE001
        ck.add("C2", "input files present and SHA-256 equal to the configuration pins", False, error=str(e))
        return _finish(ck, cfg, paths, write_report)
    # C3
    with data.GuardedH5AD(paths.h5ad) as g:
        info = g.x_info()
    obs, var_names, accessed = st.metadata()
    d = cfg["data"]
    cols = [d[k] for k in ("patient_column", "site_column", "status_column", "cell_type_column")]
    shape_ok = True
    if "n_cells" in exp:
        shape_ok = info["shape"] == (int(exp["n_cells"]), int(exp["n_genes"]))
    ck.add("C3", "X is CSR float32 of the expected shape; obs columns present; only obs/var read; forbidden groups refused",
           info["encoding_type"] == "csr_matrix" and info["dtype"] == "float32" and shape_ok
           and all(c in obs.columns for c in cols) and set(a.split("/")[0] for a in accessed) <= {"obs", "var"}
           and data.guard_selftest(paths.h5ad) and obs.index.is_unique and var_names.is_unique, x=info)
    # C4
    co = st.cohort()
    n_sites = {s: sum(co.site_of[i] == s for i in co.eligible) for s in sorted(set(co.site_of[i] for i in co.eligible))}
    by_class = {c: int((co.analysis[d["cell_type_column"]] == c).sum()) for c in cfg.classes}
    observed = {"n_individuals_mapped": int(co.obs.individual_id.nunique()),
                "n_excluded_status_cells": int(len(co.obs) - len(co.nonexcluded)),
                "n_pool_cells": int(len(co.pool)), "n_eligible": len(co.eligible),
                "eligible_by_site": n_sites, "n_analysis_cells": int(len(co.analysis)),
                "analysis_cells_by_class": by_class}
    want = {k: exp[k] for k in observed if k in exp}
    ck.add("C4", "cohort reconstruction equals the expected invariants", all(observed[k] == want[k] for k in want),
           observed=observed, expected=want)
    # C5
    asg = st.assignment()
    sizes = [len(f) for f in asg.folds]
    c5 = not asg.redraws and len(set().union(*map(set, asg.folds))) == len(co.eligible) and sum(sizes) == len(co.eligible)
    if "fold_assignment_sha256" in exp:
        c5 = c5 and asg.canonical_sha256 == exp["fold_assignment_sha256"]
    ck.add("C5", "fold assignment reproduces the canonical SHA-256; every individual in exactly one fold; no redraw",
           c5, canonical_sha256=asg.canonical_sha256, fold_sizes=sizes)
    # C6
    plans = {}
    iso_ok = True
    detail = {}
    for k in range(cfg.n_folds):
        p = st.plan(k)
        plans[k] = p
        tm, sm = co.obs.loc[p.train], co.obs.loc[p.test]
        ok = (not set(p.train_people) & set(p.test_people) and not set(p.train) & set(p.test)
              and not set(tm[d["patient_column"]]) & set(sm[d["patient_column"]])
              and set(p.fitting) | set(p.withheld) == set(p.train) and not set(p.fitting) & set(p.withheld)
              and set(p.selector) <= set(p.fitting) and set(p.classifier) <= set(p.fitting)
              and set(co.obs.loc[p.fitting].individual_id) == set(p.train_people))
        if "test_cells_per_fold" in exp:
            ok = ok and len(p.test) == int(exp["test_cells_per_fold"][k])
        iso_ok &= ok
        detail[k] = {x: len(getattr(p, x)) for x in p.LISTS}
    ck.add("C6", "per-fold isolation: individuals, sample labels and cells disjoint; Stage B partition; subsamples within fitting",
           iso_ok, sizes=detail)
    # C7
    legacy = None
    if paths.legacy_root and refs_path:
        legacy = Legacy(paths.legacy_root, refs_path)
        ver = legacy.verify_all()
        ck.add("C7a", "historical reference artifacts present and SHA-256 verified", all(ver.values()),
               n_artifacts=len(ver), failed=[k for k, v in ver.items() if not v])
        ck.add("C7b", "historical canonical fold assignment equals the reconstruction",
               legacy.canonical_sha() == asg.canonical_sha256)
        cmp = {k: compare_plan(plans[k].provenance(), legacy.plan_lists(k)) for k in plans}
        ck.add("C7c", "every fold plan (train, test, fitting, withheld, selector, classifier, individuals) equals the historical lists",
               all(c["status"] == EXACT for c in cmp.values()),
               per_fold={k: c["status"] for k, c in cmp.items()})
    else:
        ck.add("C7", "historical comparison skipped (no legacy_root configured)", True, skipped=True)
    # C8
    rd = paths.run_dir
    isolated = True
    if paths.legacy_root:
        from .reference import load_dir_aliases
        aliases = load_dir_aliases(refs_path) if refs_path else {}
        protected = [paths.legacy_root / x for x in aliases.values()]
        isolated = not any(_inside(rd, p) for p in protected)
    rd.mkdir(parents=True, exist_ok=True)
    ck.add("C8", "run directory is isolated from historical result directories and writable",
           isolated and os.access(rd, os.W_OK), run_dir=rd.name)
    # C9
    if deep:
        for k in range(cfg.n_folds):
            p = plans[k]
            blk = st.read_rows(p.fitting, allowed=set(p.fitting))
            try:
                st.read_rows(p.test[:1], allowed=set(p.fitting))
                refused = False
            except data.LeakageGuardError:
                refused = True
            conv = cfg.convention(k)
            hv = cfg["hvg"]
            pre, F = preprocessing.fit_preprocessing(blk.matrix, blk.var_names, int(hv["n_top_genes"]), hv["flavor"],
                                                     conv["hvg_convention"])
            del blk, F
            new = {"candidate_var_idx": pre.candidate_var_idx, "genes": pre.genes, "means": pre.means, "sds": pre.sds}
            if legacy:
                c = compare_preprocessing(new, legacy.preprocessing(k), acceptance or {"preprocessing_atol": 1e-12})
                ck.add(f"C9.{k}", f"fold {k}: training-only HVG and standardization replay ({conv['hvg_convention']} convention) "
                       "matches the historical checkpoint; held-out row read refused", c["status"] != "DISCREPANCY" and refused,
                       comparison=c, peak_rss_bytes=maxrss())
            else:
                ck.add(f"C9.{k}", f"fold {k}: training-only preprocessing computed; held-out row read refused", refused)
    return _finish(ck, cfg, paths, write_report)


def _finish(ck: Checks, cfg, paths, write_report: bool) -> dict:
    rep = {"validated_utc": utc_now(), "status": "PASS" if ck.ok else "FAIL", "config_sha256": cfg.sha256,
           "config_file": cfg.path.name, "checks": ck.items, "environment": environment(),
           "models_fitted": 0, "held_out_expression_read": False}
    if write_report:
        out = paths.run_dir / "validation" / f"validation_{utc_now().replace(':', '')}.json"
        write_json(out, rep)
        rep["report_file"] = str(out.relative_to(paths.run_dir))
    return rep
