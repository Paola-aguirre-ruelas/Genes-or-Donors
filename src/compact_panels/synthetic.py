"""Synthetic single-cell dataset with the structure of the real input (no real data).

It reproduces the features the pipeline must handle: several sites, individuals with
several sample labels (multi-ID), an excluded-status (LPS) individual, a recorded QC
exclusion list, individuals that fail the class-count eligibility rule, five cell
classes with class-specific marker genes, and individual-specific expression offsets.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sps

CLASSES = ("NK_16hi", "CD4.Naive", "CD4.CM", "CD8.Naive", "CD8.TE")


def make_synthetic(out_dir: Path, *, seed: int = 20260101, n_per_site: int = 8, n_genes: int = 400,
                   classes=CLASSES) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    sites = ["SiteA", "SiteB", "SiteC"]
    individuals, rows, mapping = [], [], []
    n_cls = len(classes)
    markers = {c: np.arange(i * 8, i * 8 + 8) for i, c in enumerate(classes)}
    base = rng.gamma(2.0, 0.5, size=n_genes)
    donor_genes = np.arange(8 * n_cls, 8 * n_cls + 40)
    idx = 0
    for s_i, site in enumerate(sites):
        for j in range(n_per_site):
            ind = f"IND{s_i}{j:02d}"
            individuals.append((ind, site))
            labels = [f"S{s_i}{j:02d}"] if j != 1 else [f"S{s_i}{j:02d}a", f"S{s_i}{j:02d}b"]  # multi-ID individuals
            for lab in labels:
                mapping.append({"h5ad_patient_id": lab, "individual_id": ind})
            offset = rng.normal(0, 0.35, size=n_genes)
            offset[donor_genes] = rng.normal(0, 2.5, size=len(donor_genes))   # strong individual-specific genes
            for c_i, c in enumerate(classes):
                n = int(rng.integers(30, 46))
                if j == 7 and c_i >= 2:      # one individual per site fails eligibility (2 classes only)
                    n = 6
                for _ in range(n):
                    lab = labels[idx % len(labels)]
                    day = "D0" if idx % 3 else "D7"
                    rows.append({"name": f"cell{idx:06d}-{lab}", "patient_id": lab, "Site": site, "Status": "Covid" if j % 2 else "Healthy",
                                 "author_cell_type": c, "Collection_Day": day, "ind": ind, "offset_key": len(individuals) - 1,
                                 "class_idx": c_i})
                    idx += 1
            individuals[-1] = (ind, site, offset)
    # an individual whose cells are all LPS-stimulated (removed before eligibility)
    lps_label = "S9LPS"
    mapping.append({"h5ad_patient_id": lps_label, "individual_id": "IND9LPS"})
    individuals.append(("IND9LPS", "SiteA", rng.normal(0, 0.35, size=n_genes)))
    for c_i, c in enumerate(classes):
        for _ in range(25):
            rows.append({"name": f"cell{idx:06d}-{lps_label}", "patient_id": lps_label, "Site": "SiteA", "Status": "LPS",
                         "author_cell_type": c, "Collection_Day": "D0", "ind": "IND9LPS", "offset_key": len(individuals) - 1,
                         "class_idx": c_i})
            idx += 1
    obs = pd.DataFrame(rows).set_index("name")
    n = len(obs)
    lam = np.empty((n, n_genes))
    for i, r in enumerate(obs.itertuples()):
        mu = base.copy()
        mu[markers[classes[r.class_idx]]] += 1.2
        lam[i] = np.exp(np.log(mu) + individuals[r.offset_key][2] * 0.5)
    counts = rng.poisson(lam).astype(np.float64)
    lib = counts.sum(axis=1, keepdims=True)
    X = np.log1p(counts / np.maximum(lib, 1) * 1e4).astype(np.float32)
    X[X < 0.35] = 0  # sparsify
    import anndata as ad
    o = obs[["patient_id", "Site", "Status", "author_cell_type", "Collection_Day"]].copy()
    adata = ad.AnnData(X=sps.csr_matrix(X), obs=o, var=pd.DataFrame(index=[f"GENE{g:04d}" for g in range(n_genes)]))
    h5ad = out_dir / "synthetic.h5ad"
    adata.write_h5ad(h5ad)
    pd.DataFrame(mapping).to_csv(out_dir / "identity_mapping.tsv", sep="\t", index=False, lineterminator="\n")
    # QC exclusion: the D7 cells of one sample label of the first individual of SiteB
    q = obs[(obs.patient_id == "S100") & (obs.Collection_Day == "D7")].index.tolist()
    (out_dir / "qc_exclusion_cells.txt").write_text("\n".join(sorted(q)) + "\n")
    return {"h5ad": h5ad, "identity_mapping": out_dir / "identity_mapping.tsv",
            "qc_exclusion": out_dir / "qc_exclusion_cells.txt", "n_cells": n, "n_genes": n_genes}
