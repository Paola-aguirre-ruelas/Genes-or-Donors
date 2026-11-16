"""Donor-level evaluation on held-out biological individuals.

For each held-out individual d, the evaluable classes C_d are those with at least 20 true
test cells; d must have at least 3. Macro-F1 averages the per-class F1 over C_d, but
*every* test cell of d is counted when forming true positives, false positives and false
negatives, so a non-evaluable cell predicted as an evaluable class is a false positive.
Fold and study summaries are donor-equal means; the primary contrast is the mean over
individuals of macro-F1(w=0) - macro-F1(w=1).
"""
from __future__ import annotations

import numpy as np


def evaluable_classes(truth: np.ndarray, classes, min_cells: int) -> list[str]:
    return [c for c in classes if int(np.sum(truth == c)) >= min_cells]


def score_individual(y_true, y_pred, evaluable, classes, min_cells: int = 20, min_classes: int = 3):
    truth = np.asarray(y_true, dtype=str)
    prediction = np.asarray(y_pred, dtype=str)
    if truth.shape != prediction.shape or truth.ndim != 1:
        raise ValueError("prediction length mismatch")
    if len(set(evaluable)) < min_classes or not set(evaluable) <= set(classes):
        raise ValueError("incorrect evaluable classes")
    if not (set(truth) <= set(classes) and set(prediction) <= set(classes)):
        raise ValueError("unexpected labels")
    values = {}
    for c in evaluable:
        tp = int(np.sum((truth == c) & (prediction == c)))
        fp = int(np.sum((truth != c) & (prediction == c)))
        fn = int(np.sum((truth == c) & (prediction != c)))
        if tp + fn < min_cells:
            raise ValueError("an F1 class has fewer true cells than the minimum")
        values[c] = (2.0 * tp) / (2 * tp + fp + fn)
    return float(np.mean(list(values.values()))), values


def score_donors(truth: np.ndarray, donors: np.ndarray, predictions: dict[str, np.ndarray], people: list[str],
                 classes, min_cells: int, min_classes: int) -> list[dict]:
    """Per-individual macro-F1 for each prediction set (keys are arm labels)."""
    rows = []
    for donor in people:
        mask = donors == donor
        t = truth[mask]
        ev = evaluable_classes(t, classes, min_cells)
        if len(ev) < min_classes:
            raise ValueError(f"held-out individual not evaluable ({len(ev)} classes)")
        row = {"individual_id": donor, "n_test_cells": int(mask.sum()), "n_evaluable_classes": len(ev),
               "evaluable_classes": "|".join(ev)}
        for label, pred in predictions.items():
            row[f"macro_f1_{label}"], _ = score_individual(t, pred[mask], ev, classes, min_cells, min_classes)
        rows.append(row)
    return rows
