from __future__ import annotations

from itertools import combinations_with_replacement

import numpy as np
import pandas as pd
from scipy import sparse

from ..advanced_analysis import _bh_adjust
from .metadata import FDR_SCOPE, PERMUTATION_LIMITATION
from .permutation import permutation_groups, permute_within_groups


def _edge_arrays(connectivities: sparse.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    upper = sparse.triu(connectivities, k=1).tocoo()
    return np.asarray(upper.row, dtype=int), np.asarray(upper.col, dtype=int)


def _empirical_p(observed: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not len(null) or not np.isfinite(observed):
        return np.nan
    center = float(np.mean(null))
    extreme = int(np.sum(np.abs(null - center) >= abs(observed - center)))
    return float((extreme + 1) / (len(null) + 1))


def _zscore(observed: float, null: np.ndarray) -> tuple[float, float, float]:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not len(null):
        return np.nan, np.nan, np.nan
    mean = float(np.mean(null))
    sd = float(np.std(null, ddof=1)) if len(null) > 1 else np.nan
    z = (float(observed) - mean) / sd if np.isfinite(sd) and sd > 0 else np.nan
    return mean, sd, z


def _state_counts(labels: np.ndarray, edges_i: np.ndarray, edges_j: np.ndarray, states: list[str]) -> dict[tuple[str, str], int]:
    counts = {pair: 0 for pair in combinations_with_replacement(states, 2)}
    for i, j in zip(edges_i, edges_j):
        a, b = str(labels[i]), str(labels[j])
        if a not in states or b not in states:
            continue
        pair = tuple(sorted((a, b), key=states.index))
        counts[pair] = counts.get(pair, 0) + 1
    return counts


def categorical_neighborhood_enrichment(
    connectivities: sparse.spmatrix,
    labels,
    *,
    permutations: int = 999,
    seed: int = 0,
    states: list[str] | None = None,
    valid_mask=None,
    permutation_scope: str = "whole_slide",
    strata=None,
    tissue_mask=None,
) -> pd.DataFrame:
    """Permutation-based categorical neighborhood enrichment for one sample."""
    if int(permutations) < 1:
        raise ValueError("permutations must be at least 1")
    labels = np.asarray(labels, dtype=object)
    groups, scope_valid, number_of_strata = permutation_groups(connectivities, permutation_scope, strata)
    combined_valid = scope_valid.copy()
    if valid_mask is not None:
        combined_valid &= np.asarray(valid_mask, dtype=bool)
    if tissue_mask is not None:
        combined_valid &= np.asarray(tissue_mask, dtype=bool)
    if not combined_valid.all():
        labels = labels.copy()
        labels[~combined_valid] = None
    valid = pd.notna(labels)
    if states is None:
        states = [str(x) for x in pd.unique(labels[valid]) if x is not None]
    states = [str(s) for s in states]
    edges_i, edges_j = _edge_arrays(connectivities)
    keep_edges = valid[edges_i] & valid[edges_j]
    edges_i, edges_j = edges_i[keep_edges], edges_j[keep_edges]
    observed = _state_counts(labels, edges_i, edges_j, states)
    rng = np.random.default_rng(seed)
    null_counts = {pair: [] for pair in observed}
    for _ in range(int(permutations)):
        permuted = permute_within_groups(labels, groups, rng, valid)
        counts = _state_counts(permuted, edges_i, edges_j, states)
        for pair, count in counts.items():
            null_counts[pair].append(count)
    rows: list[dict] = []
    for pair, count in observed.items():
        null = np.asarray(null_counts[pair], dtype=float)
        mean, sd, z = _zscore(float(count), null)
        rows.append({
            "label_a": pair[0],
            "label_b": pair[1],
            "observed_edge_count": int(count),
            "expected_edge_count": mean,
            "permutation_sd": sd,
            "enrichment_z": z,
            "empirical_p": _empirical_p(float(count), null),
            "observed_expected_ratio": float(count / mean) if np.isfinite(mean) and mean > 0 else np.nan,
            "n_valid_spots": int(valid.sum()),
            "n_valid_graph_edges": int(len(edges_i)),
            "permutations": int(permutations),
            "seed": int(seed),
            "permutation_scope": permutation_scope,
            "stratification_column": "user_selected" if permutation_scope == "within_user_strata" else "",
            "number_of_strata": number_of_strata,
            "tissue_only_restriction": tissue_mask is not None,
            "interpretation": (
                "neighborhood enrichment"
                if np.isfinite(z) and z >= 0
                else "neighborhood depletion"
                if np.isfinite(z)
                else "not_interpretable"
            ),
        })
    table = pd.DataFrame(rows)
    table["fdr_bh"] = _bh_adjust(table["empirical_p"].to_numpy()) if len(table) else []
    table["q_value_bh_fdr"] = table["fdr_bh"]
    table["fdr_scope"] = FDR_SCOPE
    table["permutation_limitation"] = PERMUTATION_LIMITATION
    if len(table):
        table["association_direction"] = "neutral_or_inconclusive"
        significant = pd.to_numeric(table["fdr_bh"], errors="coerce") <= 0.05
        table.loc[significant & (table["enrichment_z"] > 0), "association_direction"] = "enriched"
        table.loc[significant & (table["enrichment_z"] < 0), "association_direction"] = "depleted"
    return table


def _neighbor_count(mask_a: np.ndarray, mask_b: np.ndarray, edges_i: np.ndarray, edges_j: np.ndarray) -> int:
    if mask_a is mask_b:
        return int(np.sum(mask_a[edges_i] & mask_a[edges_j]))
    return int(np.sum((mask_a[edges_i] & mask_b[edges_j]) | (mask_b[edges_i] & mask_a[edges_j])))


def binary_mask_association(
    connectivities: sparse.spmatrix,
    masks: dict[str, np.ndarray],
    *,
    pairs: list[tuple[str, str]] | None = None,
    permutations: int = 999,
    seed: int = 0,
    valid_mask=None,
    permutation_scope: str = "whole_slide",
    strata=None,
    tissue_mask=None,
) -> pd.DataFrame:
    """Test same-spot overlap and neighboring-spot association between binary masks."""
    if int(permutations) < 1:
        raise ValueError("permutations must be at least 1")
    names = list(masks)
    if pairs is None:
        pairs = [(a, b) for i, a in enumerate(names) for b in names[i:]]
    groups, scope_valid, number_of_strata = permutation_groups(connectivities, permutation_scope, strata)
    combined_valid = scope_valid.copy()
    if valid_mask is not None:
        combined_valid &= np.asarray(valid_mask, dtype=bool)
    if tissue_mask is not None:
        combined_valid &= np.asarray(tissue_mask, dtype=bool)
    clean_masks: dict[str, np.ndarray] = {}
    for name, values in masks.items():
        arr = np.asarray(values, dtype=bool)
        arr = arr & combined_valid
        clean_masks[name] = arr
    edges_i, edges_j = _edge_arrays(connectivities)
    if not combined_valid.all():
        valid = combined_valid
        keep = valid[edges_i] & valid[edges_j]
        edges_i, edges_j = edges_i[keep], edges_j[keep]
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for label_a, label_b in pairs:
        a = clean_masks[label_a]
        b = clean_masks[label_b]
        valid = combined_valid
        same_overlap = int(np.sum(a & b & valid))
        neighbor_obs = _neighbor_count(a, b, edges_i, edges_j)
        same_null: list[int] = []
        neighbor_null: list[int] = []
        for _ in range(int(permutations)):
            permuted_b = permute_within_groups(b, groups, rng, valid)
            same_null.append(int(np.sum(a & permuted_b & valid)))
            neighbor_null.append(_neighbor_count(a, permuted_b, edges_i, edges_j))
        for mode, observed, null_values in (
            ("same_spot_overlap", same_overlap, np.asarray(same_null, dtype=float)),
            ("neighboring_spot_association", neighbor_obs, np.asarray(neighbor_null, dtype=float)),
        ):
            mean, sd, z = _zscore(float(observed), null_values)
            rows.append({
                "mask_a": label_a,
                "mask_b": label_b,
                "mode": mode,
                "observed_count": int(observed),
                "expected_count": mean,
                "permutation_sd": sd,
                "association_z": z,
                "empirical_p": _empirical_p(float(observed), null_values),
                "observed_expected_ratio": float(observed / mean) if np.isfinite(mean) and mean > 0 else np.nan,
                "n_valid_spots": int(valid.sum()),
                "n_valid_graph_edges": int(len(edges_i)),
                "permutations": int(permutations),
                "seed": int(seed),
                "permutation_scope": permutation_scope,
                "stratification_column": "user_selected" if permutation_scope == "within_user_strata" else "",
                "number_of_strata": number_of_strata,
                "tissue_only_restriction": tissue_mask is not None,
                "direction": (
                    "association"
                    if np.isfinite(z) and z >= 0
                    else "avoidance_or_depletion"
                    if np.isfinite(z)
                    else "not_interpretable"
                ),
            })
    table = pd.DataFrame(rows)
    table["fdr_bh"] = _bh_adjust(table["empirical_p"].to_numpy()) if len(table) else []
    table["q_value_bh_fdr"] = table["fdr_bh"]
    table["fdr_scope"] = FDR_SCOPE
    table["permutation_limitation"] = PERMUTATION_LIMITATION
    if len(table):
        table["association_direction"] = "neutral_or_inconclusive"
        significant = pd.to_numeric(table["fdr_bh"], errors="coerce") <= 0.05
        table.loc[significant & (table["association_z"] > 0), "association_direction"] = "enriched"
        table.loc[significant & (table["association_z"] < 0), "association_direction"] = "depleted"
    return table
