from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components


PERMUTATION_SCOPES = {"whole_slide", "within_connected_components", "within_user_strata"}


def permutation_groups(connectivities: sparse.spmatrix, scope: str, strata=None) -> tuple[np.ndarray, np.ndarray, int]:
    scope = str(scope or "whole_slide")
    if scope not in PERMUTATION_SCOPES:
        raise ValueError("permutation scope must be whole_slide, within_connected_components, or within_user_strata")
    n = int(connectivities.shape[0])
    valid = np.ones(n, dtype=bool)
    if scope == "whole_slide":
        groups = np.zeros(n, dtype=int)
    elif scope == "within_connected_components":
        _, groups = connected_components(connectivities > 0, directed=False, return_labels=True)
        groups = np.asarray(groups, dtype=int)
    else:
        if strata is None:
            raise ValueError("within_user_strata requires a categorical adata.obs column")
        series = pd.Series(np.asarray(strata, dtype=object))
        valid = series.notna().to_numpy(dtype=bool)
        groups = np.full(n, -1, dtype=int)
        codes, _ = pd.factorize(series[valid].astype(str), sort=True)
        groups[valid] = codes
    number = int(len(np.unique(groups[valid]))) if valid.any() else 0
    return groups, valid, number


def permute_within_groups(values, groups: np.ndarray, rng: np.random.Generator, valid_mask=None) -> np.ndarray:
    result = np.asarray(values).copy()
    groups = np.asarray(groups, dtype=int)
    valid = np.ones(len(result), dtype=bool) if valid_mask is None else np.asarray(valid_mask, dtype=bool).copy()
    for group in np.unique(groups[valid]):
        indices = np.flatnonzero(valid & (groups == group))
        if len(indices) > 1:
            result[indices] = rng.permutation(result[indices])
    return result
