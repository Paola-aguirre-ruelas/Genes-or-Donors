"""Downstream cell-class classifier on a 64-gene panel.

Multinomial logistic regression, L2 penalty, C = 1, lbfgs, max_iter = 1000, intercept,
no class weights, trained on the standardized classifier subsample of the fitting cells.
Models are stored as plain arrays (no pickled estimators). Prediction is the arg-max of
``X @ coef.T + intercept`` over the stored class order.
"""
from __future__ import annotations

import time
import warnings

import numpy as np


def hyperparameters(cfg) -> dict:
    c = cfg["classifier"]
    return {"l1_ratio": 0, "C": float(c["C"]), "solver": str(c["solver"]), "max_iter": int(c["max_iter"]),
            "fit_intercept": True, "class_weight": None}


def train_logistic(X: np.ndarray, y: np.ndarray, cfg) -> dict:
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    hp = hyperparameters(cfg)
    model = LogisticRegression(**hp)
    t0 = time.perf_counter()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter("always")
        model.fit(np.ascontiguousarray(X), y)
    seconds = time.perf_counter() - t0
    cw = [str(w.message) for w in W if issubclass(w.category, ConvergenceWarning)]
    other = [f"{w.category.__name__}: {str(w.message)[:200]}" for w in W if not issubclass(w.category, ConvergenceWarning)]
    n_iter = int(np.max(model.n_iter_))
    return {"coef": model.coef_.astype(np.float64), "intercept": model.intercept_.astype(np.float64),
            "classes": model.classes_.astype(str), "n_iter": n_iter,
            "converged": bool(not cw and n_iter < hp["max_iter"]), "convergence_warnings": cw,
            "other_warnings": other, "seconds": round(seconds, 3), "hyperparameters": hp}


def predict(coef: np.ndarray, intercept: np.ndarray, classes: np.ndarray, X: np.ndarray) -> np.ndarray:
    logits = X @ coef.T + intercept
    if not np.all(np.isfinite(logits)):
        raise ValueError("non-finite class scores")
    return np.asarray(classes)[np.argmax(logits, axis=1)]
