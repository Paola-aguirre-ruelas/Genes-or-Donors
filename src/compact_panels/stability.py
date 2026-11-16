"""Training-only numerical stability checks (historical exploratory controls).

* G3 fixed-alpha check: one cold (non-warm-started) MultiTaskLasso fit at the panel's
  selected alpha with tol = 1e-5 and max_iter = 10000. Exact-64 panels are *stable* if the
  new support has 64 genes and Jaccard >= 0.95 with the panel; fallback panels if the
  top-64-by-norm of the new fit has Jaccard >= 0.95.
* G4 identity-permutation sanity check: on fold 0, donor labels of the selector rows are
  permuted with ``default_rng(20261115)``; a cold-fit alpha search (tol 1e-5) for the
  w = 1 target with permuted donors selects a panel compared with the w = 0 panel
  (Jaccard >= 0.90 expected when donor labels carry no structure).
* H1 diagnosis (fold 3, w = 1): cold fits at four fixed alphas around the fallback alpha.

These reproduce the historical sensitivity checks; they are not decoy or confirmatory
controls.
"""
from __future__ import annotations

import time
import warnings

import numpy as np

K_DEFAULT = 64


def jaccard(a, b) -> float:
    aa, bb = set(a), set(b)
    if not (aa | bb):
        raise ValueError("Jaccard undefined for two empty sets")
    return len(aa & bb) / len(aa | bb)


def top_k_positive(norms: np.ndarray, k: int) -> list[int]:
    """Top-k among strictly positive rows; ties by candidate index (historical review helper)."""
    positive = np.flatnonzero(norms > 0)
    if len(positive) < k:
        raise ValueError(f"cannot select {k} genes from {len(positive)} active rows")
    by_norm = np.lexsort((positive, -norms[positive]))
    return sorted(positive[by_norm[:k]].tolist())


def cold_fit(X, Y, alpha, *, tol: float, max_iter: int):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import MultiTaskLasso
    if not (alpha > 0 and tol > 0):
        raise ValueError("invalid solver settings")
    est = MultiTaskLasso(alpha=float(alpha), fit_intercept=True, max_iter=max_iter, tol=tol, warm_start=False,
                         selection="cyclic", copy_X=True)
    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as found:
        warnings.simplefilter("always")
        est.fit(X, Y)
    runtime = time.perf_counter() - t0
    cw = [str(w.message) for w in found if issubclass(w.category, ConvergenceWarning)]
    norms = np.linalg.norm(est.coef_, axis=0)
    selected = np.flatnonzero(norms > 0)
    yy = Y - Y.mean(axis=0)
    threshold = tol * float(np.sum(yy * yy)) / len(Y)
    gap, iters = float(est.dual_gap_), int(est.n_iter_)
    converged = bool(not cw and iters < max_iter and gap <= threshold * (1 + 1e-8) + 1e-12)
    if iters == 0 and len(selected) > 0:
        converged = False
    return norms, {"alpha": float(alpha), "nnz": int(len(selected)), "n_iter": iters, "dual_gap": gap,
                   "gap_threshold": threshold, "converged": converged, "convergence_warnings": cw,
                   "seconds": round(runtime, 3)}


def judge_fixed_alpha(reference_idx, reference_exact64: bool, support, top, threshold: float):
    base = set(reference_idx)
    raw_j = jaccard(base, set(support))
    top_j = jaccard(base, top) if top is not None else None
    if reference_exact64:
        status = "STABLE_AT_FIXED_ALPHA" if len(support) == len(base) and raw_j >= threshold else "REVIEW_REQUIRED"
    else:
        status = "FALLBACK_TOP64_STABLE_AT_FIXED_ALPHA" if top_j is not None and top_j >= threshold else "REVIEW_REQUIRED"
    return status, raw_j, top_j


def advance_cold_search(state: dict, *, alpha: float, rec: dict, norms: np.ndarray, k: int,
                        max_steps: int, extra_steps: int, trigger: int, lower_frac: float):
    """Historical cold-search decision rule (G4). Returns (state, selected_idx or None, exact)."""
    step = rec["step"]
    if rec["nnz"] > k and (state["best_alpha_over_k"] is None or alpha > state["best_alpha_over_k"]):
        state["best_alpha_over_k"] = alpha
        state["best_norms_over_k"] = norms.tolist()
    if rec["nnz"] == k:
        return state, sorted(np.flatnonzero(norms > 0).tolist()), True
    if rec["nnz"] > k:
        state["low"] = alpha
    else:
        state["high"] = alpha
    if step >= max_steps and state["best_norms_over_k"] is not None:
        best = np.asarray(state["best_norms_over_k"])
        if int(np.sum(best > 0)) - k <= trigger or step == max_steps + extra_steps:
            return state, top_k_positive(best, k), False
    if step == max_steps + extra_steps:
        state["status"] = "STOP_REVIEW"
        return state, None, False
    if state["best_alpha_over_k"] is None:
        cand = max(alpha / 2, state["alpha_max"] * lower_frac)
        if not cand < alpha:
            raise RuntimeError("cannot descend above k active genes")
        state["next_alpha"] = cand
    else:
        if not state["high"] > state["low"] > 0:
            raise RuntimeError("invalid alpha bracket")
        state["next_alpha"] = float(np.sqrt(state["low"] * state["high"]))
    return state, None, False
