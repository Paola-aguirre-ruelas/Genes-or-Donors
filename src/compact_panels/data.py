"""Guarded, read-only access to the processed H5AD file.

The reader only ever opens the ``obs`` and ``var`` groups and the three CSR components
of ``X``. Expression rows are returned only for cell names on an explicit allow-list,
so a training stage cannot read withheld or held-out cells even by mistake. Every
HDF5 access is logged and can be audited.

The row loader reproduces the historical implementation exactly (block-wise reads in
file order, then restoration of the requested order).
"""
from __future__ import annotations

from dataclasses import dataclass

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sps

METADATA_GROUPS = {"obs", "var"}
MATRIX_KEYS = {"X", "X/data", "X/indices", "X/indptr"}
FORBIDDEN_GROUPS = ("layers", "raw", "obsm", "obsp", "uns", "varm", "varp")
STRING_COLUMNS = ("patient_id", "Site", "Status", "author_cell_type", "Collection_Day")


class LeakageGuardError(PermissionError):
    pass


class GuardedH5AD:
    """h5py wrapper restricted to obs, var and the CSR parts of X."""

    def __init__(self, path, *, allow_matrix: bool = True):
        self._f = h5py.File(path, "r")
        self.accessed: list[str] = []
        self.allowed = METADATA_GROUPS | (MATRIX_KEYS if allow_matrix else set())

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self) -> None:
        self._f.close()

    def __getitem__(self, key: str):
        k = key.strip("/")
        top = k.split("/")[0]
        self.accessed.append(k)
        if k not in self.allowed and not (top in METADATA_GROUPS and top in self.allowed):
            raise LeakageGuardError(f"leakage guard: access to '{key}' is forbidden")
        return self._f[k]

    # ------------------------------------------------------------------ metadata
    def read_obs(self, string_columns=STRING_COLUMNS) -> pd.DataFrame:
        from anndata.io import read_elem
        obs = read_elem(self["obs"])
        obs.index = obs.index.astype(str)
        for c in string_columns:
            if c in obs.columns:
                obs[c] = obs[c].astype(str)
        return obs

    def read_var_names(self) -> pd.Index:
        from anndata.io import read_elem
        return pd.Index(read_elem(self["var"]).index.astype(str))

    # ------------------------------------------------------------------ matrix
    def x_info(self) -> dict:
        a = self["X"].attrs
        return {"encoding_type": str(a["encoding-type"]), "shape": tuple(int(s) for s in a["shape"]),
                "dtype": str(self["X/data"].dtype)}

    def indptr(self) -> np.ndarray:
        return self["X/indptr"][:]

    def load_rows(self, names: list[str], row_of: dict[str, int], indptr: np.ndarray, n_genes: int,
                  allowed_rows: set[str], block: int = 2000) -> sps.csr_matrix:
        bad = [n for n in names if n not in allowed_rows]
        if bad:
            raise LeakageGuardError(f"leakage guard: {len(bad)} requested rows are outside the allowed set (e.g. {bad[0]})")
        r = np.array([row_of[n] for n in names], dtype=np.int64)
        order = np.argsort(r, kind="stable")
        rs = r[order]
        D, I = self["X/data"], self["X/indices"]
        datas, inds, lens = [], [], []
        for b in range(0, len(rs), block):
            blk = rs[b:b + block]
            lo, hi = int(indptr[blk[0]]), int(indptr[blk[-1] + 1])
            d, ix = D[lo:hi], I[lo:hi]
            for row in blk:
                s, e = int(indptr[row]) - lo, int(indptr[row + 1]) - lo
                datas.append(d[s:e]); inds.append(ix[s:e]); lens.append(e - s)
        ptr = np.zeros(len(rs) + 1, dtype=np.int64)
        ptr[1:] = np.cumsum(lens)
        M = sps.csr_matrix((np.concatenate(datas), np.concatenate(inds), ptr), shape=(len(rs), n_genes))
        inv = np.argsort(order, kind="stable")
        return M[inv]


def guard_selftest(path) -> bool:
    """True if every forbidden group (and X when matrix access is disabled) raises."""
    with GuardedH5AD(path, allow_matrix=False) as g:
        for grp in FORBIDDEN_GROUPS + ("X",):
            try:
                g[grp]
                return False
            except LeakageGuardError:
                pass
            except KeyError:
                return False
    return True


def load_metadata(path) -> tuple[pd.DataFrame, pd.Index, list[str]]:
    """Read obs and var only (no expression values)."""
    with GuardedH5AD(path, allow_matrix=False) as g:
        obs = g.read_obs()
        var_names = g.read_var_names()
        accessed = list(g.accessed)
    return obs, var_names, accessed


@dataclass
class ExpressionBlock:
    matrix: sps.csr_matrix
    var_names: pd.Index
    accessed: list[str]


def read_expression(path, names: list[str], obs_index: pd.Index, *, allowed: set[str]) -> ExpressionBlock:
    """Read the CSR rows of ``names`` (all must be in ``allowed``) in the given order."""
    if len(names) != len(set(names)):
        raise LeakageGuardError("duplicate row names requested")
    with GuardedH5AD(path) as g:
        info = g.x_info()
        if info["encoding_type"] != "csr_matrix":
            raise ValueError(f"X must be CSR, found {info['encoding_type']}")
        indptr = g.indptr()
        var_names = g.read_var_names()
        if len(obs_index) != info["shape"][0] or len(var_names) != info["shape"][1]:
            raise ValueError("obs/var sizes do not match X")
        row_of = {n: i for i, n in enumerate(obs_index)}
        M = g.load_rows(list(names), row_of, indptr, info["shape"][1], allowed)
        return ExpressionBlock(M, var_names, list(g.accessed))
