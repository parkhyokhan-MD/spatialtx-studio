from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from ..advanced_analysis import _bh_adjust
from .metadata import FDR_SCOPE, PERMUTATION_LIMITATION, json_safe
from .permutation import permutation_groups, permute_within_groups


def _edge_data(connectivities: sparse.spmatrix) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    upper = sparse.triu(connectivities, k=1).tocoo()
    weights = np.asarray(upper.data, dtype=float)
    weights[~np.isfinite(weights)] = 0.0
    return np.asarray(upper.row, dtype=int), np.asarray(upper.col, dtype=int), weights


def _statistic(x: np.ndarray, y: np.ndarray, i: np.ndarray, j: np.ndarray, w: np.ndarray, symmetric: bool) -> tuple[float, int]:
    valid = np.isfinite(x[i]) & np.isfinite(y[j]) & np.isfinite(w) & (w > 0)
    if symmetric:
        valid = valid & np.isfinite(y[i]) & np.isfinite(x[j])
    if not valid.any():
        return np.nan, 0
    if symmetric:
        values = 0.5 * (x[i[valid]] * y[j[valid]] + y[i[valid]] * x[j[valid]])
    else:
        values = x[i[valid]] * y[j[valid]]
    weights = w[valid]
    return float(np.sum(weights * values) / np.sum(weights)), int(valid.sum())


def _empirical_p(observed: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not len(null) or not np.isfinite(observed):
        return np.nan
    center = float(np.mean(null))
    extreme = int(np.sum(np.abs(null - center) >= abs(observed - center)))
    return float((extreme + 1) / (len(null) + 1))


def continuous_edge_association(
    connectivities: sparse.spmatrix,
    x,
    y,
    *,
    statistic_name: str,
    x_source: str,
    y_source: str,
    permutations: int = 999,
    seed: int = 0,
    symmetric: bool = True,
    permutation_scope: str = "whole_slide",
    strata=None,
    tissue_mask=None,
    active_graph: str = "",
    graph_parameters: dict | None = None,
    stratification_column: str = "",
) -> dict:
    if int(permutations) < 1:
        raise ValueError("permutations must be at least 1")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) != connectivities.shape[0] or len(y) != connectivities.shape[0]:
        raise ValueError("continuous variables must have one value per graph node")
    i, j, w = _edge_data(connectivities)
    groups, scope_valid, number_of_strata = permutation_groups(connectivities, permutation_scope, strata)
    valid_observations = scope_valid & np.isfinite(x) & np.isfinite(y)
    if tissue_mask is not None:
        valid_observations &= np.asarray(tissue_mask, dtype=bool)
    x_active = x.copy()
    y_active = y.copy()
    x_active[~valid_observations] = np.nan
    y_active[~valid_observations] = np.nan
    observed, valid_edges = _statistic(x_active, y_active, i, j, w, symmetric=symmetric)
    null: list[float] = []
    if valid_edges:
        rng = np.random.default_rng(seed)
        for _ in range(int(permutations)):
            permuted = permute_within_groups(y_active, groups, rng, valid_observations)
            value, _ = _statistic(x_active, permuted, i, j, w, symmetric=symmetric)
            if np.isfinite(value):
                null.append(value)
    null_values = np.asarray(null, dtype=float)
    mean = float(np.mean(null_values)) if len(null_values) else np.nan
    sd = float(np.std(null_values, ddof=1)) if len(null_values) > 1 else np.nan
    z = (observed - mean) / sd if np.isfinite(sd) and sd > 0 and np.isfinite(observed) else np.nan
    return {
        "statistic": statistic_name,
        "statistic_name": statistic_name,
        "field_x": x_source,
        "field_y": y_source,
        "X_source": x_source,
        "Y_source": y_source,
        "status": "ok" if valid_edges and np.isfinite(observed) else "insufficient_edges",
        "raw_weighted_statistic": observed if valid_edges else np.nan,
        "permutation_null_mean": mean,
        "permutation_null_sd": sd,
        "z_score": z,
        "empirical_p": _empirical_p(observed, null_values),
        "valid_observation_count": int(valid_observations.sum()),
        "valid_edge_count": int(valid_edges),
        "missing_value_count": int(len(x) - valid_observations.sum()),
        "active_graph": active_graph,
        "graph_parameters": json_safe(graph_parameters or {}),
        "permutation_scope": permutation_scope,
        "stratification_column": stratification_column if permutation_scope == "within_user_strata" else "",
        "number_of_strata": number_of_strata,
        "permutations": int(permutations),
        "seed": int(seed),
        "tissue_only_restriction": tissue_mask is not None,
        "symmetric": bool(symmetric),
        "fdr_scope": FDR_SCOPE,
        "permutation_limitation": PERMUTATION_LIMITATION,
    }


def continuous_edge_statistics(
    connectivities: sparse.spmatrix,
    fields: dict[str, np.ndarray],
    *,
    permutations: int = 999,
    seed: int = 0,
    symmetric: bool = True,
    custom_pairs: list[tuple[str, str, str]] | None = None,
    permutation_scope: str = "whole_slide",
    strata=None,
    tissue_mask=None,
    active_graph: str = "",
    graph_parameters: dict | None = None,
    stratification_column: str = "",
) -> pd.DataFrame:
    """Graph-weighted continuous edge-interaction statistics."""
    if int(permutations) < 1:
        raise ValueError("permutations must be at least 1")
    requested = [
        ("I_CS", "C", "S"),
        ("I_RR", "R", "R"),
        ("I_RH", "R", "H"),
        ("I_RV", "R", "V"),
        ("I_HV", "H", "V"),
    ]
    requested.extend(custom_pairs or [])
    rows: list[dict] = []
    for name, x_key, y_key in requested:
        if x_key not in fields or y_key not in fields:
            rows.append({
                "statistic": name,
                "field_x": x_key,
                "field_y": y_key,
                "status": "missing_field",
                "raw_weighted_statistic": np.nan,
                "valid_edge_count": 0,
                "permutations": int(permutations),
                "seed": int(seed),
            })
            continue
        rows.append(continuous_edge_association(
            connectivities,
            fields[x_key],
            fields[y_key],
            statistic_name=name,
            x_source=x_key,
            y_source=y_key,
            permutations=permutations,
            seed=seed,
            symmetric=symmetric,
            permutation_scope=permutation_scope,
            strata=strata,
            tissue_mask=tissue_mask,
            active_graph=active_graph,
            graph_parameters=graph_parameters,
            stratification_column=stratification_column,
        ))
    table = pd.DataFrame(rows)
    table["fdr_bh"] = _bh_adjust(table["empirical_p"].to_numpy()) if "empirical_p" in table else np.nan
    table["q_value_bh_fdr"] = table["fdr_bh"]
    table["fdr_scope"] = FDR_SCOPE
    table["permutation_limitation"] = PERMUTATION_LIMITATION
    if len(table) and "z_score" in table:
        table["observed_vs_null_direction"] = "neutral_or_inconclusive"
        significant = pd.to_numeric(table["fdr_bh"], errors="coerce") <= 0.05
        table.loc[significant & (table["z_score"] > 0), "observed_vs_null_direction"] = "enriched"
        table.loc[significant & (table["z_score"] < 0), "observed_vs_null_direction"] = "depleted"
    return table
