"""Supervision targets for the multi-task selector.

Y_A is the cell-class block (5 columns), Y_B the donor-identity block (one column per
training individual, lexicographic). Each block is one-hot, column-centred and scaled
so that ||block||_F^2 = n. Condition w uses Y = [Y_A, sqrt(w) * Y_B].

Historical note: the w = 0 panels were fitted with the five-column Y_A alone (not a
71-column Y with a zero block); w = 1 uses all 5 + 66 = 71 columns.
"""
from __future__ import annotations

import numpy as np


def block(labels, columns) -> np.ndarray:
    B = (np.asarray(labels)[:, None] == np.asarray(columns)[None, :]).astype(np.float64)
    if not np.all(B.sum(axis=1) == 1):
        raise ValueError("label outside the column set")
    B -= B.mean(axis=0)
    tv = B.var(axis=0, ddof=0).sum()
    return B / np.sqrt(tv)


def make_y(cell_class, individual, individuals, classes, w: float) -> np.ndarray:
    return np.hstack([block(cell_class, list(classes)), np.sqrt(w) * block(individual, list(individuals))])


def targets_for_arm(arm: str, cell_class, individual, individuals, classes) -> np.ndarray:
    if arm == "w0":
        return block(cell_class, list(classes))
    if arm == "w1":
        return make_y(cell_class, individual, individuals, classes, 1.0)
    raise ValueError(f"unknown arm {arm!r}")
