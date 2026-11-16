"""MultiTaskLasso panel selection with the historical alpha-search state machines.

One warm-started ``MultiTaskLasso`` per condition. A fit is *converged* when
``n_iter_ < max_iter``, no ConvergenceWarning was raised and
``dual_gap_ <= tol * ||Yc||_F^2 / n``; a failed fit is retried once from the same
starting coefficients with ``max_iter = 10000``.

Two historical search variants are implemented without modification:

``fold0_recovery`` (the separate script that produced the fold-0 panels)
    bracket ``[lo, alpha_max]`` with ``lo = 0.01 alpha_max`` for w = 0 (a prior bound from the
    earlier frozen pilot) and ``1e-4 alpha_max`` for w = 1; first alpha ``0.1 alpha_max``;
    w = 0 bisects geometrically from the second fit; w = 1 halves alpha until a fit has
    more than k active genes, then bisects. At fit 25 the fallback (k largest coefficient
    row norms of the largest alpha with > k genes) is taken if it discards <= 8 genes,
    otherwise bisection continues to fit 50. With no fit above k at fit 25 the run stops.

``cv_halving`` (folds 1-3)
    bracket ``[1e-4 alpha_max, alpha_max]``, first alpha ``0.1 alpha_max``, halving until a
    fit exceeds k, then geometric bisection; same fallback rule from fit 25, but without a
    fit above k the search simply continues and stops with an error after fit 50.

The pure function :func:`advance` encodes each decision so that recorded historical
alpha traces can be replayed and checked without fitting anything.
"""
from __future__ import annotations

import json
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .io_utils import IntegrityError, require, sha256_file, write_json, write_npz
from .provenance import maxrss

VARIANTS = ("fold0_recovery", "cv_halving")


class SearchStop(RuntimeError):
    """The search reached a state that requires review (not converged, memory, no bracket)."""


def alpha_max(X: np.ndarray, Y: np.ndarray) -> float:
    Xc = X - X.mean(axis=0)
    Yc = Y - Y.mean(axis=0)
    return float(np.max(np.linalg.norm(Xc.T @ Yc, axis=1)) / X.shape[0])


def top_k_by_norm(norms: np.ndarray, k: int) -> np.ndarray:
    """Indices of the k largest row norms; ties broken by lower gene index; returned sorted."""
    order = np.lexsort((np.arange(len(norms)), -norms))
    return np.sort(order[:k])


# ====================================================================== estimator wrapper
class Selector:
    """Warm-started MultiTaskLasso with per-fit convergence records and one retry."""

    def __init__(self, X: np.ndarray, Y: np.ndarray, *, tol: float, max_iter: int, retry_max_iter: int,
                 label: str, log=print):
        from sklearn.linear_model import MultiTaskLasso
        self.X, self.Y, self.tol, self.label = X, Y, tol, label
        self.max_iter, self.retry_max_iter = int(max_iter), int(retry_max_iter)
        self.n = X.shape[0]
        Yc = Y - Y.mean(axis=0)
        self.gap_threshold = tol * float(np.sum(Yc * Yc)) / self.n
        self.est = MultiTaskLasso(alpha=1.0, fit_intercept=True, max_iter=self.max_iter, tol=tol, warm_start=True,
                                  selection="cyclic", copy_X=True)
        self.fits: list[dict] = []
        self.norms: dict[float, np.ndarray] = {}
        self.log = log

    def _attempt(self, alpha: float, max_iter: int) -> dict:
        from sklearn.exceptions import ConvergenceWarning
        self.est.set_params(alpha=alpha, max_iter=max_iter)
        t0 = time.perf_counter()
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter("always")
            self.est.fit(self.X, self.Y)
        el = time.perf_counter() - t0
        conv_warn = [str(w.message)[:200] for w in W if issubclass(w.category, ConvergenceWarning)]
        other = [f"{w.category.__name__}: {str(w.message)[:200]}" for w in W if not issubclass(w.category, ConvergenceWarning)]
        n_iter, gap = int(self.est.n_iter_), float(self.est.dual_gap_)
        ok = (n_iter < max_iter) and not conv_warn and (gap <= self.gap_threshold)
        return {"alpha": alpha, "max_iter": max_iter, "n_iter": n_iter, "dual_gap": gap,
                "gap_threshold": self.gap_threshold, "convergence_warnings": conv_warn, "other_warnings": other,
                "converged": bool(ok), "seconds": round(el, 3)}

    def fit_at(self, alpha: float, step: int, phase: str) -> int:
        start = self.est.coef_.copy() if hasattr(self.est, "coef_") else None
        rec = self._attempt(alpha, self.max_iter)
        rec.update(step=step, phase=phase, retried=False)
        if not rec["converged"]:
            first = dict(rec)
            if start is None:
                del self.est.coef_
            else:
                self.est.coef_ = start.copy()
            rec = self._attempt(alpha, self.retry_max_iter)
            rec.update(step=step, phase=phase, retried=True, first_attempt=first)
            self.est.set_params(max_iter=self.max_iter)
        norms = np.linalg.norm(self.est.coef_, axis=0)
        nz = int(np.sum(np.any(self.est.coef_ != 0.0, axis=0)))
        rec["n_nonzero"] = nz
        self.norms[alpha] = norms
        self.fits.append(rec)
        self.log(f"  [{self.label}] {phase} step {step:2d} alpha={alpha:.6e} nnz={nz:4d} n_iter={rec['n_iter']} "
                 f"gap={rec['dual_gap']:.3e} thr={self.gap_threshold:.3e} conv={rec['converged']} "
                 f"retried={rec['retried']} {rec['seconds']:.1f}s")
        return nz


# ====================================================================== state machine
@dataclass
class Outcome:
    panel_idx: list[int]
    exact: bool
    discarded: int
    alpha_final: float


def initial_state(variant: str, arm: str, amax: float, sel_cfg: dict, conv: dict) -> dict:
    if variant not in VARIANTS:
        raise ValueError(f"unknown search variant {variant!r}")
    lower = float(sel_cfg["lower_frac"])
    skip_halving = False
    if variant == "fold0_recovery" and arm == "w0":
        lower = float(conv["w0_lower_frac"])
        skip_halving = True
    return {"status": "SEARCH", "variant": variant, "arm": arm, "alpha_max": float(amax),
            "lo": float(lower * amax), "hi": float(amax),
            "next_alpha": float(float(sel_cfg["first_alpha_frac"]) * amax),
            "skip_halving": skip_halving, "fits": [], "best_alpha": None, "best_nnz": None}


def advance(state: dict, *, alpha: float, nz: int, norms: np.ndarray | None, best_norms: np.ndarray | None,
            step: int, sel_cfg: dict) -> tuple[dict, np.ndarray | None, Outcome | None]:
    """Apply the historical decision rule after the fit at ``alpha`` (fit number ``step``)."""
    K = int(sel_cfg["k"])
    max_steps, extra, trigger = int(sel_cfg["max_steps"]), int(sel_cfg["extra_steps"]), int(sel_cfg["fallback_trigger"])
    variant = state["variant"]
    if nz > K and (state["best_alpha"] is None or alpha > state["best_alpha"]):
        state["best_alpha"], state["best_nnz"] = alpha, nz
        best_norms = None if norms is None else norms.copy()
    if nz == K:
        panel = [] if norms is None else np.flatnonzero(norms > 0).tolist()
        return state, best_norms, Outcome(panel, True, 0, float(alpha))
    if nz > K:
        state["lo"] = alpha
    else:
        state["hi"] = alpha
    if step >= max_steps:
        if variant == "fold0_recovery":
            if state["best_alpha"] is None:
                raise SearchStop(f"{state['arm']}: no fit above {K} genes; cannot apply fallback")
            discarded = int(state["best_nnz"] - K)
            if step >= max_steps + extra or discarded <= trigger:
                panel = [] if best_norms is None else top_k_by_norm(best_norms, K).tolist()
                return state, best_norms, Outcome(panel, False, discarded, float(state["best_alpha"]))
        elif state["best_alpha"] is not None:
            discarded = int(state["best_nnz"] - K)
            if discarded <= trigger or step == max_steps + extra:
                panel = [] if best_norms is None else top_k_by_norm(best_norms, K).tolist()
                return state, best_norms, Outcome(panel, False, discarded, float(state["best_alpha"]))
    if not state["hi"] > state["lo"] > 0:
        raise SearchStop("degenerate alpha bracket")
    if state["best_alpha"] is None and not state["skip_halving"]:
        nxt = float(max(alpha / 2, state["lo"]))
        if not nxt < alpha:
            raise SearchStop("unable to bracket the target gene count")
        state["next_alpha"] = nxt
    else:
        state["next_alpha"] = float(np.sqrt(state["lo"] * state["hi"]))
    return state, best_norms, None


def replay_trace(variant: str, arm: str, amax: float, trace: list[tuple[float, int]], sel_cfg: dict, conv: dict):
    """Re-derive the alpha sequence from recorded (alpha, nnz) pairs; returns the derived alphas and outcome."""
    state = initial_state(variant, arm, amax, sel_cfg, conv)
    derived, outcome = [], None
    for step, (alpha_rec, nz) in enumerate(trace, start=1):
        alpha = state["next_alpha"]
        derived.append(alpha)
        state, _, outcome = advance(state, alpha=alpha, nz=int(nz), norms=None, best_norms=None, step=step,
                                    sel_cfg=sel_cfg)
        if outcome is not None:
            break
    return derived, outcome, state


# ====================================================================== checkpointed search
def _save_checkpoint(path: Path, state: dict, coef: np.ndarray, best_norms: np.ndarray) -> None:
    write_npz(path, overwrite=True, meta=np.array(json.dumps(state, sort_keys=True)), coef=coef,
              best_norms=best_norms)


def run_search(X: np.ndarray, Y: np.ndarray, *, arm: str, fold: int, genes: np.ndarray, cfg, dest: Path,
               max_new_fits: int, provenance: dict, log=print) -> int:
    """Resumable search for one condition. Writes checkpoint_<arm>.npz after every fit and
    panel_<arm>.json at completion. Returns the number of new fits."""
    sel = cfg["selector"]
    conv = cfg.convention(fold)
    variant = conv["search_variant"]
    ckpt = dest / f"checkpoint_{arm}.npz"
    panel_path = dest / f"panel_{arm}.json"
    amax = alpha_max(X, Y)
    require(np.isfinite(amax) and amax > 0, "invalid alpha_max")
    selector = Selector(X, Y, tol=float(sel["tol"]), max_iter=int(sel["max_iter"]),
                        retry_max_iter=int(sel["retry_max_iter"]), label=f"fold{fold}/{arm}", log=log)
    n_genes = X.shape[1]
    if ckpt.exists():
        with np.load(ckpt, allow_pickle=False) as z:
            state = json.loads(str(z["meta"].item()))
            coef = z["coef"].copy()
            best_norms = z["best_norms"].copy()
        require(state.get("provenance") == provenance and np.isclose(state["alpha_max"], amax, rtol=1e-12, atol=0)
                and coef.shape == (Y.shape[1], n_genes) and state["variant"] == variant,
                f"fold{fold}/{arm}: incompatible checkpoint")
        if state["status"] == "COMPLETE":
            require(panel_path.exists(), f"fold{fold}/{arm}: checkpoint COMPLETE but panel missing")
            log(f"fold{fold}/{arm}: panel already complete")
            return 0
        require(state["status"] == "SEARCH", f"fold{fold}/{arm}: checkpoint status {state['status']} requires review")
        if state["fits"]:
            selector.est.coef_ = coef
        selector.fits = list(state["fits"])
        log(f"fold{fold}/{arm}: resumed after {len(state['fits'])} fits")
    else:
        require(not panel_path.exists(), f"fold{fold}/{arm}: panel exists without checkpoint")
        state = initial_state(variant, arm, amax, sel, conv)
        state["provenance"] = provenance
        best_norms = np.zeros(n_genes, dtype=np.float64)
        log(f"fold{fold}/{arm}: variant={variant} alpha_max={amax:.8g} first alpha={state['next_alpha']:.8g} "
            f"bracket=[{state['lo']:.8g},{state['hi']:.8g}]")
    K = int(sel["k"])
    limit = int(sel["max_steps"]) + int(sel["extra_steps"])
    used = 0
    while used < max_new_fits:
        step = len(state["fits"]) + 1
        if variant == "cv_halving" and step > limit:
            state["status"] = "FIT_LIMIT_REVIEW_REQUIRED"
            _save_checkpoint(ckpt, state, selector.est.coef_, best_norms)
            raise SearchStop(f"fold{fold}/{arm}: exceeded the maximum of {limit} fits")
        alpha = float(state["next_alpha"])
        nz = selector.fit_at(alpha, step, "search")
        used += 1
        rec = selector.fits[-1]
        rec["peak_rss_bytes_after_fit"] = int(maxrss())
        state["fits"] = list(selector.fits)
        if not rec["converged"]:
            state["status"] = "NONCONVERGED_REVIEW_REQUIRED"
            _save_checkpoint(ckpt, state, selector.est.coef_, best_norms)
            raise SearchStop(f"fold{fold}/{arm}: fit did not converge; state saved")
        if rec["peak_rss_bytes_after_fit"] > int(sel["memory_gate_bytes"]):
            state["status"] = "MEMORY_GATE_REVIEW_REQUIRED"
            _save_checkpoint(ckpt, state, selector.est.coef_, best_norms)
            raise SearchStop(f"fold{fold}/{arm}: resident memory above the gate")
        state, best_norms_new, outcome = advance(state, alpha=alpha, nz=nz, norms=selector.norms[alpha],
                                                 best_norms=best_norms, step=step, sel_cfg=sel)
        if best_norms_new is not None:
            best_norms = best_norms_new
        if outcome is not None:
            require(len(outcome.panel_idx) == K and len(set(outcome.panel_idx)) == K, "selected panel has the wrong size")
            write_json(panel_path, {
                "fold": fold, "arm": arm, "variant": variant, "k": K,
                "panel_idx": [int(i) for i in outcome.panel_idx],
                "panel_genes": [str(genes[i]) for i in outcome.panel_idx],
                "exact64": bool(outcome.exact), "discarded": int(outcome.discarded),
                "alpha_final": float(outcome.alpha_final), "alpha_max": float(amax), "n_fits": step,
                "all_converged": all(f["converged"] for f in state["fits"]),
                "n_retries": int(sum(bool(f.get("retried")) for f in state["fits"])),
                "fit_seconds_sum": float(sum(f["seconds"] for f in state["fits"])),
                "provenance": provenance})
            state["status"] = "COMPLETE"
            log(f"fold{fold}/{arm}: PANEL {K} genes exact={outcome.exact} discarded={outcome.discarded} fits={step}")
        _save_checkpoint(ckpt, state, selector.est.coef_, best_norms)
        if outcome is not None:
            return used
    log(f"fold{fold}/{arm}: {used} new fits; checkpoint saved, rerun to resume")
    return used


def load_panel(dest: Path, arm: str) -> dict:
    p = dest / f"panel_{arm}.json"
    require(p.exists(), f"missing panel {p}")
    return json.loads(p.read_text())


def panel_sha(dest: Path, arm: str) -> str:
    return sha256_file(dest / f"panel_{arm}.json")
