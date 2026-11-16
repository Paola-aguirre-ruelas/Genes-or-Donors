"""Training-only highly variable gene (HVG) selection and standardization.

Both are estimated from the fitting cells of the fold only. Two historical call
conventions of ``scanpy.pp.highly_variable_genes(flavor='seurat')`` are preserved:

* ``"pilot"`` (fold 0): ``AnnData(X=csr, var=index, uns={'log1p': {'base': None}})``,
  ``inplace=False``, and a STOP if two genes tie at the cutoff;
* ``"cv"`` (folds 1-3): ``AnnData(X=csr.copy())``, ``inplace=True``, STOP unless exactly
  ``n_top_genes`` genes are flagged.

Candidate genes keep ascending H5AD column order. Standardization uses float64 per-gene
means and standard deviations (two-pass, ddof = 0); a zero SD stops the run.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sps


@dataclass
class Preprocessing:
    candidate_var_idx: np.ndarray
    genes: np.ndarray
    means: np.ndarray
    sds: np.ndarray
    hvg_convention: str
    hvg_cutoff_margin: float | None
    warnings: list[str]


def select_hvg(csr: sps.csr_matrix, var_names: pd.Index, n_top: int, flavor: str, convention: str):
    import scanpy as sc
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter("always")
        if convention == "pilot":
            A = sc.AnnData(X=csr, var=pd.DataFrame(index=var_names), uns={"log1p": {"base": None}})
            hv = sc.pp.highly_variable_genes(A, flavor=flavor, n_top_genes=n_top, inplace=False)
            del A
            mask = hv["highly_variable"].to_numpy()
            dn = hv["dispersions_norm"].to_numpy(dtype=float)
            fin = np.sort(dn[~np.isnan(dn)])[::-1]
            if int(mask.sum()) != n_top or not fin[n_top - 1] > fin[n_top]:
                raise ValueError(f"HVG tie or count mismatch at the cutoff ({int(mask.sum())} flagged)")
            margin = float(fin[n_top - 1] - fin[n_top])
        elif convention == "cv":
            a = sc.AnnData(X=csr.copy())
            sc.pp.highly_variable_genes(a, flavor=flavor, n_top_genes=n_top, subset=False, inplace=True)
            mask = a.var["highly_variable"].to_numpy(dtype=bool)
            del a
            margin = None
        else:
            raise ValueError(f"unknown HVG convention {convention!r}")
    cand = np.flatnonzero(mask)
    if len(cand) != n_top or len(np.unique(cand)) != n_top:
        raise ValueError(f"HVG ambiguity: {len(cand)} genes flagged, expected {n_top}")
    msgs = [f"{w.category.__name__}: {str(w.message)[:200]}" for w in W]
    return cand, margin, msgs


def fit_standardization(F: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mu = F.mean(axis=0)
    sd = np.sqrt(((F - mu) ** 2).mean(axis=0))
    if not (np.all(sd > 0) and np.all(np.isfinite(sd))):
        raise ValueError("zero or invalid gene variance among candidate genes")
    return mu, sd


def fit_preprocessing(csr_fit: sps.csr_matrix, var_names: pd.Index, n_top: int, flavor: str, convention: str):
    """Returns the Preprocessing record and the dense float64 candidate matrix of the fitting cells."""
    if str(csr_fit.dtype) != "float32":
        raise ValueError(f"unexpected expression dtype {csr_fit.dtype}")
    cand, margin, msgs = select_hvg(csr_fit, var_names, n_top, flavor, convention)
    F = csr_fit[:, cand].toarray().astype(np.float64)
    mu, sd = fit_standardization(F)
    genes = np.asarray(var_names[cand], dtype=str)
    if len(set(genes.tolist())) != n_top:
        raise ValueError("non-unique candidate gene names")
    return Preprocessing(cand, genes, mu, sd, convention, margin, msgs), F


def standardize_rows(F: np.ndarray, rows: np.ndarray, mu: np.ndarray, sd: np.ndarray, fortran: bool) -> np.ndarray:
    X = (F[rows] - mu) / sd
    X = np.asfortranarray(X) if fortran else np.ascontiguousarray(X)
    if not np.all(np.isfinite(X)):
        raise ValueError("non-finite standardized values")
    return X


def positions(sorted_names: list[str], subset: list[str]) -> np.ndarray:
    arr = np.asarray(sorted_names)
    pos = np.searchsorted(arr, np.asarray(subset))
    if not np.array_equal(arr[pos], np.asarray(subset)):
        raise ValueError("subset rows not located in the sorted fitting list")
    return pos
