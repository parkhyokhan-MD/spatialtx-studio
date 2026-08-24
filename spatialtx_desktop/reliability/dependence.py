from __future__ import annotations

from collections.abc import Mapping, Sequence
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from .models import AxisReliabilityResult, ReliabilityConfig


DEPENDENCE_COLUMNS = [
    "sample_id",
    "axis_i",
    "axis_j",
    "dependence_type",
    "metric_i",
    "metric_j",
    "pearson_correlation",
    "spearman_correlation",
    "valid_spot_count",
    "missing_undefined_fraction",
    "permutation_p_value",
    "bh_fdr",
    "metric_inference_eligible",
    "metric_inference_qc_reason",
    "qc_status",
    "permutation_scope",
]


def bh_adjust(pvalues: Sequence[float]) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(values))
    if not len(finite_indices):
        return adjusted
    finite_values = values[finite_indices]
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    corrected = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    restored = np.empty_like(corrected)
    restored[order] = np.clip(corrected, 0.0, 1.0)
    adjusted[finite_indices] = restored
    return adjusted


def _correlations(left: np.ndarray, right: np.ndarray) -> tuple[float, float, np.ndarray]:
    valid = np.isfinite(left) & np.isfinite(right)
    if valid.sum() < 3:
        return np.nan, np.nan, valid
    x, y = left[valid], right[valid]
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan, valid
    pearson = float(np.corrcoef(x, y)[0, 1])
    spearman = float(stats.spearmanr(x, y).statistic)
    return pearson, spearman, valid


def _permute_within(values: np.ndarray, groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    output = values.copy()
    for group in np.unique(groups):
        indices = np.flatnonzero(groups == group)
        output[indices] = values[rng.permutation(indices)]
    return output


def _permutation_p(
    left: np.ndarray,
    right: np.ndarray,
    valid: np.ndarray,
    config: ReliabilityConfig,
    rng: np.random.Generator,
    groups: np.ndarray | None,
) -> float:
    if valid.sum() < config.minimum_valid_spots:
        return np.nan
    x, y = left[valid], right[valid]
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    observed = abs(float(np.corrcoef(x, y)[0, 1]))
    local_groups = np.zeros(len(x), dtype=int) if groups is None else np.asarray(groups)[valid]
    exceed = 0
    for _ in range(int(config.permutation_iterations)):
        permuted = _permute_within(y, local_groups, rng)
        candidate = abs(float(np.corrcoef(x, permuted)[0, 1]))
        exceed += int(np.isfinite(candidate) and candidate >= observed)
    return float((exceed + 1) / (int(config.permutation_iterations) + 1))


def _direction_inference_support(
    result: AxisReliabilityResult,
    config: ReliabilityConfig,
) -> tuple[bool, str]:
    """Apply the same conservative sample-level Direction gate used by pair exports."""

    total_n = int(len(result.valid_input))
    valid_input_n = int(np.asarray(result.valid_input, dtype=bool).sum())
    defined_n = int((
        np.asarray(result.valid_input, dtype=bool)
        & np.asarray(result.direction_defined, dtype=bool)
        & np.isfinite(np.asarray(result.direction_D, dtype=float))
    ).sum())
    valid_input_fraction = float(valid_input_n / total_n) if total_n else 0.0
    defined_fraction = float(defined_n / valid_input_n) if valid_input_n else 0.0
    eligible = bool(
        valid_input_n >= int(config.minimum_valid_spots)
        and valid_input_fraction >= float(config.minimum_valid_fraction)
        and defined_n >= int(config.minimum_valid_spots)
        and defined_fraction >= float(config.minimum_valid_fraction)
    )
    if eligible:
        reason = "ok"
    elif valid_input_n < int(config.minimum_valid_spots):
        reason = "upstream_activity_insufficient_valid_spots"
    elif valid_input_fraction < float(config.minimum_valid_fraction):
        reason = "upstream_activity_low_valid_fraction"
    elif defined_n < int(config.minimum_valid_spots):
        reason = "defined_spot_count_below_minimum"
    else:
        reason = "defined_fraction_below_pass_threshold"
    return eligible, reason


def compute_axis_dependence(
    axes: Mapping[str, AxisReliabilityResult],
    config: ReliabilityConfig | Mapping | None = None,
    *,
    sample_id: str = "",
    permutation_groups: Sequence[int] | None = None,
) -> pd.DataFrame:
    """Report inter-axis dependence without transforming or removing axes."""

    cfg = ReliabilityConfig.from_value(config)
    cfg.validate()
    if not cfg.dependence_qc or len(axes) < 2:
        return pd.DataFrame(columns=DEPENDENCE_COLUMNS)
    lengths = {len(result.C) for result in axes.values()}
    if len(lengths) != 1:
        raise ValueError("All axes must contain the same number of spots for dependence QC.")
    n_spots = next(iter(lengths), 0)
    groups = None if permutation_groups is None else np.asarray(permutation_groups)
    if groups is not None and groups.shape != (n_spots,):
        raise ValueError("permutation_groups must have one value per spot.")

    rows: list[dict] = []
    comparisons: list[
        tuple[
            str,
            str,
            str,
            str,
            str,
            np.ndarray,
            np.ndarray,
            AxisReliabilityResult,
            AxisReliabilityResult,
        ]
    ] = []
    for axis_i, axis_j in combinations(sorted(axes), 2):
        first, second = axes[axis_i], axes[axis_j]
        comparisons.extend((
            (axis_i, axis_j, "direction_dependence", "direction_D", "direction_D", first.direction_D, second.direction_D, first, second),
            (axis_i, axis_j, "activity_dependence", "activity_A", "activity_A", first.activity_A, second.activity_A, first, second),
            (axis_i, axis_j, "balance_dependence", "balance_B", "balance_B", first.balance_B, second.balance_B, first, second),
            (axis_i, axis_j, "cross_dependence", "direction_D", "activity_A", first.direction_D, second.activity_A, first, second),
            (axis_j, axis_i, "cross_dependence", "direction_D", "activity_A", second.direction_D, first.activity_A, second, first),
        ))
    for index, (
        axis_i,
        axis_j,
        kind,
        metric_i,
        metric_j,
        left,
        right,
        result_i,
        result_j,
    ) in enumerate(comparisons):
        pearson, spearman, valid = _correlations(left, right)
        rng = np.random.default_rng(np.random.SeedSequence([int(cfg.seed), index]))
        support_checks: list[tuple[bool, str]] = []
        if metric_i == "direction_D":
            support_checks.append(_direction_inference_support(result_i, cfg))
        if metric_j == "direction_D":
            support_checks.append(_direction_inference_support(result_j, cfg))
        metric_inference_eligible = all(check[0] for check in support_checks)
        support_reasons = [reason for eligible, reason in support_checks if not eligible]
        metric_inference_qc_reason = (
            "ok" if metric_inference_eligible else ";".join(dict.fromkeys(support_reasons))
        )
        p_value = (
            _permutation_p(left, right, valid, cfg, rng, groups)
            if metric_inference_eligible
            else np.nan
        )
        valid_count = int(valid.sum())
        if not metric_inference_eligible:
            qc_status = "insufficient_direction_defined_support"
        elif valid_count < int(cfg.minimum_valid_spots):
            qc_status = "insufficient_valid_spots"
        elif not np.isfinite(pearson) or not np.isfinite(spearman):
            qc_status = "undefined_constant_or_invalid"
        elif max(abs(pearson), abs(spearman)) >= float(cfg.dependence_abs_correlation_warning):
            qc_status = "warning_high_dependence"
        else:
            qc_status = "ok"
        rows.append({
            "sample_id": sample_id,
            "axis_i": axis_i,
            "axis_j": axis_j,
            "dependence_type": kind,
            "metric_i": metric_i,
            "metric_j": metric_j,
            "pearson_correlation": pearson,
            "spearman_correlation": spearman,
            "valid_spot_count": valid_count,
            "missing_undefined_fraction": float(1.0 - valid_count / n_spots) if n_spots else np.nan,
            "permutation_p_value": p_value,
            "bh_fdr": np.nan,
            "metric_inference_eligible": metric_inference_eligible,
            "metric_inference_qc_reason": metric_inference_qc_reason,
            "qc_status": qc_status,
            "permutation_scope": "within_supplied_groups" if groups is not None else "within_sample_global",
        })
    table = pd.DataFrame(rows, columns=DEPENDENCE_COLUMNS)
    table["bh_fdr"] = bh_adjust(table["permutation_p_value"].to_numpy(dtype=float))
    return table


def direction_dependence_matrix(
    table: pd.DataFrame,
    axis_names: Sequence[str],
) -> pd.DataFrame:
    axes = [str(axis) for axis in axis_names]
    matrix = pd.DataFrame(np.nan, index=axes, columns=axes, dtype=float)
    for axis in axes:
        matrix.loc[axis, axis] = 1.0
    if table.empty:
        return matrix
    rows = table.loc[table["dependence_type"].eq("direction_dependence")]
    for _, row in rows.iterrows():
        value = float(row["pearson_correlation"])
        matrix.loc[str(row["axis_i"]), str(row["axis_j"])] = value
        matrix.loc[str(row["axis_j"]), str(row["axis_i"])] = value
    return matrix
