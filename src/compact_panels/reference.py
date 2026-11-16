"""Historical reference artifacts and component-level comparisons.

The historical tree (``legacy_root``) is opened read-only. Every referenced file is
SHA-256-verified against ``configs/historical_references.toml`` before it is used.

Agreement classes (tolerances in ``configs/acceptance.toml``):

* ``EXACT`` -- identical gene sets / cell lists and bit-identical numbers;
* ``NUMERICAL_AGREEMENT`` -- identical discrete outcomes (lists, panels, predictions-derived
  scores) with continuous values inside the declared tolerance;
* ``DISCREPANCY`` -- anything else.
"""
from __future__ import annotations

import json
import tomllib
from pathlib import Path

import numpy as np

from .io_utils import IntegrityError, read_json, require, sha256_file, sha256_lines

EXACT, NUMERIC, DISCREPANCY, MISSING = "EXACT", "NUMERICAL_AGREEMENT", "DISCREPANCY", "NOT_AVAILABLE"


def worst(statuses) -> str:
    order = {EXACT: 0, NUMERIC: 1, MISSING: 2, DISCREPANCY: 3}
    statuses = list(statuses)
    return max(statuses, key=lambda s: order[s]) if statuses else MISSING


ALIASES_LOCAL = "historical_dirs.local.toml"


def load_dir_aliases(refs_path: Path) -> dict[str, str]:
    """Map the ``@alias`` directory tokens of the public reference file to local directory names.

    The mapping lives in ``historical_dirs.local.toml`` next to the reference file (git-ignored; see
    ``historical_dirs.example.toml``) or in the file named by ``$CGP_HISTORICAL_DIRS``. The public
    repository carries only the aliases, never the private directory names.
    """
    import os
    cand = os.environ.get("CGP_HISTORICAL_DIRS") or str(Path(refs_path).parent / ALIASES_LOCAL)
    if not Path(cand).is_file():
        return {}
    with open(cand, "rb") as f:
        return {str(k): str(v) for k, v in tomllib.load(f).get("dirs", {}).items()}


class Legacy:
    def __init__(self, root: Path, refs_path: Path):
        self.root = Path(root)
        with open(refs_path, "rb") as f:
            self.refs = tomllib.load(f)
        self.pins: dict[str, str] = self.refs.get("sha256", {})
        self.aliases = load_dir_aliases(refs_path)
        self._verified: set[str] = set()

    def resolve(self, rel: str) -> Path:
        """``@alias/sub/path`` -> ``root/<local dir>/sub/path``; pins stay keyed by the alias form."""
        head, _, tail = rel.partition("/")
        if head.startswith("@"):
            name = head[1:]
            require(name in self.aliases, f"historical directory alias @{name} is not configured "
                                          f"(see configs/historical_dirs.example.toml)")
            return self.root / self.aliases[name] / tail
        return self.root / rel

    def path(self, rel: str) -> Path:
        p = self.resolve(rel)
        if rel not in self._verified:
            require(p.is_file(), f"historical artifact missing: {rel}")
            want = self.pins.get(rel)
            require(want is not None, f"historical artifact {rel} has no pinned SHA-256")
            got = sha256_file(p)
            require(got == want, f"historical artifact changed: {rel}")
            self._verified.add(rel)
        return p

    def fold(self, fold: int) -> dict:
        return self.refs["fold0"] if fold == 0 else {k: v.format(k=fold) for k, v in self.refs["folds"].items()}

    def verify_all(self) -> dict:
        out = {}
        for rel in sorted(self.pins):
            try:
                self.path(rel)
                out[rel] = True
            except IntegrityError:
                out[rel] = False
        return out

    # ------------------------------------------------------------------ loaders
    def plan_lists(self, fold: int) -> dict:
        """{list_name: {"size", "sha256"}} for the historical fold plan."""
        f = self.fold(fold)
        if fold == 0:
            out = {}
            for key in ("train", "test", "fitting", "withheld", "selector", "classifier", "train_people", "test_people"):
                names = self.path(f[key]).read_text().split()
                out[key] = {"size": len(names), "sha256": sha256_lines(names)}
            return out
        return read_json(self.path(f["provenance"]))["lists"]

    def preprocessing(self, fold: int) -> dict:
        f = self.fold(fold)
        p = self.path(f["preprocessing"])  # SHA verified first; the fold-0 file holds legacy object arrays
        with np.load(p, allow_pickle=(fold == 0)) as z:
            genes = z["candidate_genes"] if fold == 0 else z["genes"]
            return {"candidate_var_idx": z["candidate_var_idx"].astype(np.int64), "genes": genes.astype(str),
                    "means": z["means"].astype(np.float64), "sds": z["sds"].astype(np.float64)}

    def checkpoint(self, fold: int, arm: str) -> dict:
        with np.load(self.path(self.fold(fold)[f"checkpoint_{arm}"]), allow_pickle=False) as z:
            return {"meta": json.loads(str(z["meta"].item())), "coef": z["coef"].copy(), "best_norms": z["best_norms"].copy()}

    def panel(self, fold: int, arm: str) -> dict:
        return read_json(self.path(self.fold(fold)[f"panel_{arm}"]))

    def model(self, fold: int, arm: str) -> dict:
        with np.load(self.path(self.fold(fold)[f"model_{arm}"]), allow_pickle=False) as z:
            idx = z["candidate_indices"] if "candidate_indices" in z.files else z["panel_idx"]
            return {"coef": z["coef"].copy(), "intercept": z["intercept"].copy(), "classes": z["classes"].astype(str),
                    "panel_idx": idx.astype(np.int64)}

    def evaluation(self, fold: int) -> dict:
        return read_json(self.path(self.fold(fold)["evaluation"]))

    def random_dir(self, fold: int) -> str:
        return self.refs["random"]["fold_dir"].format(k=fold)

    def random_manifest(self, fold: int) -> dict:
        return read_json(self.path(self.random_dir(fold) + "/manifest.json"))

    def random_model(self, fold: int, draw: int) -> dict:
        with np.load(self.path(self.random_dir(fold) + f"/model_{draw:02d}.npz"), allow_pickle=False) as z:
            return {"coef": z["coef"].copy(), "intercept": z["intercept"].copy(), "classes": z["classes"].astype(str),
                    "cols": z["cols"].copy()}

    def random_evaluation(self, fold: int) -> dict:
        return read_json(self.path(self.random_dir(fold) + "/evaluation.json"))

    def stability(self, key: str) -> dict:
        return read_json(self.path(self.refs["stability"][key]))

    def canonical_sha(self) -> str:
        import hashlib
        return hashlib.sha256(self.path(self.refs["splits"]["canonical"]).read_text().strip().encode()).hexdigest()


# ====================================================================== comparisons
def _maxabs(a, b) -> float:
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        return float("inf")
    return float(np.max(np.abs(a - b))) if a.size else 0.0


def _num_status(a, b, atol: float, discrete_equal: bool = True) -> tuple[str, float]:
    a_, b_ = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    d = _maxabs(a_, b_)
    if not discrete_equal:
        return DISCREPANCY, d
    if a_.shape == b_.shape and np.array_equal(a_, b_):
        return EXACT, 0.0
    return (NUMERIC if d <= atol else DISCREPANCY), d


def compare_plan(new: dict, old: dict) -> dict:
    lists = {k: {"new": new[k], "old": old[k], "equal": new[k] == old[k]} for k in old}
    return {"component": "fold_plan", "status": EXACT if all(v["equal"] for v in lists.values()) else DISCREPANCY,
            "lists": lists}


def compare_preprocessing(new: dict, old: dict, tol: dict) -> dict:
    cand_eq = np.array_equal(new["candidate_var_idx"], old["candidate_var_idx"]) and np.array_equal(new["genes"], old["genes"])
    s_mu, d_mu = _num_status(new["means"], old["means"], float(tol["preprocessing_atol"]), cand_eq)
    s_sd, d_sd = _num_status(new["sds"], old["sds"], float(tol["preprocessing_atol"]), cand_eq)
    return {"component": "preprocessing", "status": worst([s_mu, s_sd]) if cand_eq else DISCREPANCY,
            "hvg_candidates_identical": bool(cand_eq), "means_max_abs_diff": d_mu, "sds_max_abs_diff": d_sd}


def compare_search(new_ckpt: dict, new_panel: dict, old_ckpt: dict, old_panel: dict, tol: dict) -> dict:
    nf, of = new_ckpt["meta"]["fits"], old_ckpt["meta"]["fits"]
    n = min(len(nf), len(of))
    alphas_new = [f["alpha"] for f in nf]
    alphas_old = [f["alpha"] for f in of]
    nnz_eq = [f["n_nonzero"] for f in nf] == [f["n_nonzero"] for f in of]
    iter_eq = [f["n_iter"] for f in nf] == [f["n_iter"] for f in of]
    same_len = len(nf) == len(of)
    alpha_rel = max((abs(a - b) / abs(b) for a, b in zip(alphas_new[:n], alphas_old[:n])), default=0.0)
    gap_rel = max((abs(a["dual_gap"] - b["dual_gap"]) / max(abs(b["dual_gap"]), 1e-300) for a, b in zip(nf[:n], of[:n])), default=0.0)
    genes_eq = sorted(new_panel["panel_genes"]) == sorted(old_panel["panel_genes"])
    idx_eq = sorted(new_panel["panel_idx"]) == sorted(old_panel["panel_idx"])
    alpha_final_rel = abs(new_panel["alpha_final"] - old_panel["alpha_final"]) / abs(old_panel["alpha_final"])
    coef_d = _maxabs(new_ckpt["coef"], old_ckpt["coef"])
    discrete = same_len and nnz_eq and genes_eq and idx_eq and new_panel["exact64"] == old_panel["exact64"]
    if discrete and alphas_new == alphas_old and coef_d == 0.0 and iter_eq:
        status = EXACT
    elif discrete and alpha_rel <= float(tol["alpha_rtol"]) and coef_d <= float(tol["coef_atol"]):
        status = NUMERIC
    else:
        status = DISCREPANCY
    return {"component": "panel_search", "status": status, "n_fits_new": len(nf), "n_fits_old": len(of),
            "alpha_trace_bit_identical": alphas_new == alphas_old, "alpha_max_rel_diff": alpha_rel,
            "nnz_trace_identical": nnz_eq, "n_iter_trace_identical": iter_eq, "dual_gap_max_rel_diff": gap_rel,
            "final_coef_max_abs_diff": coef_d, "panel_genes_identical": genes_eq, "panel_idx_identical": idx_eq,
            "exact64_new": new_panel["exact64"], "exact64_old": old_panel["exact64"],
            "discarded_new": new_panel["discarded"], "discarded_old": old_panel.get("discarded"),
            "alpha_final_rel_diff": alpha_final_rel, "jaccard": len(set(new_panel["panel_genes"]) & set(old_panel["panel_genes"]))
            / len(set(new_panel["panel_genes"]) | set(old_panel["panel_genes"]))}


def compare_model(new: dict, old: dict, tol: dict) -> dict:
    discrete = np.array_equal(new["classes"], old["classes"]) and np.array_equal(np.sort(new["panel_idx"]), np.sort(old["panel_idx"]))
    s1, d1 = _num_status(new["coef"], old["coef"], float(tol["classifier_atol"]), discrete)
    s2, d2 = _num_status(new["intercept"], old["intercept"], float(tol["classifier_atol"]), discrete)
    return {"component": "classifier", "status": worst([s1, s2]), "classes_and_genes_identical": bool(discrete),
            "coef_max_abs_diff": d1, "intercept_max_abs_diff": d2}


def compare_scores(new_rows: list[dict], old_rows: list[dict], keys: list[str], tol: dict) -> dict:
    o = {str(r["individual_id"]): r for r in old_rows}
    n = {str(r["individual_id"]): r for r in new_rows}
    same_people = set(o) == set(n)
    worst_diff, statuses = 0.0, []
    for d in sorted(set(o) & set(n)):
        for k in keys:
            s, diff = _num_status(n[d][k], o[d][k], float(tol["score_atol"]))
            statuses.append(s)
            worst_diff = max(worst_diff, diff)
    status = worst(statuses) if same_people else DISCREPANCY
    return {"status": status, "individuals_identical": same_people, "n_individuals": len(n), "max_abs_diff": worst_diff}
