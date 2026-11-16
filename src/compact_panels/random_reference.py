"""Expression-matched random 64-gene reference panels.

Rule (fixed before the historical random-reference scores were computed):
rank the candidate genes by increasing training-only mean expression (ties by candidate
index), cut them into equal bins, count how many genes of the fold's w = 0 panel fall in
each bin, and draw that many genes per bin uniformly without replacement. Each draw uses
its own generator ``default_rng([seed, fold, draw, salt])``.
"""
from __future__ import annotations

import numpy as np


def random_panels(means, ref_idx, fold: int, *, n_genes: int, k: int, bins: int, draws: int, seed: int, salt: int):
    means = np.asarray(means, dtype=np.float64)
    ref = np.asarray(ref_idx, dtype=np.int64)
    if n_genes % bins:
        raise ValueError("number of candidate genes must be divisible by the number of bins")
    width = n_genes // bins
    if means.shape != (n_genes,) or not np.all(np.isfinite(means)):
        raise ValueError("invalid training means")
    if ref.shape != (k,) or len(set(ref.tolist())) != k or not np.all((ref >= 0) & (ref < n_genes)):
        raise ValueError("invalid reference panel")
    order = np.lexsort((np.arange(n_genes), means))
    rank_bin = np.empty(n_genes, dtype=np.int64)
    rank_bin[order] = np.arange(n_genes) // width
    counts = np.bincount(rank_bin[ref], minlength=bins)
    if counts.sum() != k or counts.max() > width:
        raise ValueError("invalid bin occupancy")
    panels = []
    for draw in range(draws):
        rng = np.random.default_rng([seed, fold, draw, salt])
        selected = []
        for b, n in enumerate(counts):
            if not n:
                continue
            pool = order[b * width:(b + 1) * width]
            selected.extend(rng.choice(pool, size=int(n), replace=False).tolist())
        idx = sorted(selected)
        if len(idx) != k or len(set(idx)) != k:
            raise ValueError("random panel duplicate or size error")
        if not np.array_equal(np.bincount(rank_bin[idx], minlength=bins), counts):
            raise ValueError("mean-expression histogram differs from the reference")
        panels.append([int(i) for i in idx])
    return panels, counts.tolist()
