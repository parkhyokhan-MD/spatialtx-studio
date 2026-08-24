from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .models import AxisReliabilityResult, ReliabilityConfig


def _one_dimensional(values, label: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{label} must be a one-dimensional score array.")
    return array


def compute_axis_reliability(
    C,
    S,
    config: ReliabilityConfig | Mapping | None = None,
    *,
    axis: str = "FRAME2.6_CS",
    activity_C=None,
    activity_S=None,
    balance_score_source: str = "legacy_signed_cs",
    balance_score_domain: str = "signed",
    activity_score_source: str = "explicit_nonnegative_program_abundance",
    activity_score_domain: str = "nonnegative",
    activity_source_transformations: str = "explicit input; no clipping or shift",
    activity_source_version: str = "v0.65-nonnegative-program-mean-v1",
) -> AxisReliabilityResult:
    """Compute legacy Balance and nonnegative Activity quantities separately.

    ``C`` and ``S`` are the preserved signed v0.6 fields and are used only for
    legacy Balance.  ``activity_C`` and ``activity_S`` are explicit
    pre-centering program abundances used only for Activity, Direction, and
    Co-activation.  The latter default to C/S solely for the standalone public
    calculation API; the application integration always supplies them
    explicitly.  No clipping, shifting, normalization, or missing-value
    replacement occurs.
    """

    cfg = ReliabilityConfig.from_value(config)
    cfg.validate()
    c_values = _one_dimensional(C, "C")
    s_values = _one_dimensional(S, "S")
    if c_values.shape != s_values.shape:
        raise ValueError("C and S score arrays must have the same shape.")
    activity_c_values = _one_dimensional(
        c_values if activity_C is None else activity_C,
        "activity_C",
    )
    activity_s_values = _one_dimensional(
        s_values if activity_S is None else activity_S,
        "activity_S",
    )
    if activity_c_values.shape != c_values.shape or activity_s_values.shape != c_values.shape:
        raise ValueError("Legacy and Activity score arrays must have the same shape.")

    n = len(c_values)
    balance_finite_input = np.isfinite(c_values) & np.isfinite(s_values)
    finite_input = np.isfinite(activity_c_values) & np.isfinite(activity_s_values)
    nonnegative_input = finite_input & (activity_c_values >= 0) & (activity_s_values >= 0)
    valid_input = nonnegative_input.copy()

    balance = np.full(n, np.nan, dtype=float)
    balance[balance_finite_input] = (
        c_values[balance_finite_input] - s_values[balance_finite_input]
    )

    activity_balance = np.full(n, np.nan, dtype=float)
    activity_balance[valid_input] = (
        activity_c_values[valid_input] - activity_s_values[valid_input]
    )

    activity = np.full(n, np.nan, dtype=float)
    activity[valid_input] = activity_c_values[valid_input] + activity_s_values[valid_input]

    ca_strength = np.full(n, np.nan, dtype=float)
    ca_strength[valid_input] = (
        activity[valid_input] - np.abs(activity_balance[valid_input])
    )

    direction_defined = valid_input & (activity > float(cfg.epsilon))
    direction = np.full(n, np.nan, dtype=float)
    direction[direction_defined] = activity_balance[direction_defined] / (
        activity[direction_defined] + float(cfg.epsilon)
    )
    ca_fraction = np.full(n, np.nan, dtype=float)
    ca_fraction[direction_defined] = ca_strength[direction_defined] / (
        activity[direction_defined] + float(cfg.epsilon)
    )

    status = np.full(n, "invalid_nonfinite_input", dtype=object)
    status[finite_input & ~nonnegative_input] = "invalid_negative_score"
    status[valid_input] = "valid_continuous"
    status[valid_input & ~direction_defined] = "undefined_direction_zero_activity"

    states = np.full(n, None, dtype=object)
    if cfg.classification_enabled:
        activity_threshold = float(cfg.activity_threshold)
        direction_threshold = float(cfg.direction_threshold)
        low = valid_input & np.isfinite(activity) & (activity <= activity_threshold)
        active = direction_defined & ~low
        c_dominant = active & (direction >= direction_threshold)
        s_dominant = active & (direction <= -direction_threshold)
        coactive = active & ~(c_dominant | s_dominant)
        states[low] = "low_activity"
        states[c_dominant] = "c_dominant_active"
        states[s_dominant] = "s_dominant_active"
        states[coactive] = "active_coactivation_candidate"
        classified = np.fromiter((value is not None for value in states), dtype=bool, count=n)
        status[classified] = "valid_classified"

    return AxisReliabilityResult(
        axis=str(axis),
        C=c_values.copy(),
        S=s_values.copy(),
        balance_B=balance,
        balance_finite_input=balance_finite_input,
        activity_C=activity_c_values.copy(),
        activity_S=activity_s_values.copy(),
        activity_balance=activity_balance,
        activity_A=activity,
        direction_D=direction,
        ca_strength=ca_strength,
        ca_fraction=ca_fraction,
        reliability_state=states,
        status=status,
        finite_input=finite_input,
        nonnegative_input=nonnegative_input,
        valid_input=valid_input,
        direction_defined=direction_defined,
        balance_score_source=str(balance_score_source),
        balance_score_domain=str(balance_score_domain),
        activity_score_source=str(activity_score_source),
        activity_score_domain=str(activity_score_domain),
        activity_source_transformations=str(activity_source_transformations),
        activity_source_version=str(activity_source_version),
    )


def compute_reliability_axes(
    axis_scores: Mapping[str, tuple[Sequence[float], Sequence[float]]],
    config: ReliabilityConfig | Mapping | None = None,
    *,
    activity_axis_scores: Mapping[str, tuple[Sequence[float], Sequence[float]]] | None = None,
    source_metadata: Mapping[str, str] | None = None,
) -> dict[str, AxisReliabilityResult]:
    """Apply the same pure reliability calculation to explicit paired-pole axes."""

    cfg = ReliabilityConfig.from_value(config)
    cfg.validate()
    activity_scores = activity_axis_scores or axis_scores
    metadata = dict(source_metadata or {})
    missing = set(axis_scores) - set(activity_scores)
    if missing:
        raise ValueError(f"Activity scores are missing for axes: {', '.join(sorted(missing))}")
    return {
        str(axis): compute_axis_reliability(
            C,
            S,
            cfg,
            axis=str(axis),
            activity_C=activity_scores[axis][0],
            activity_S=activity_scores[axis][1],
            **metadata,
        )
        for axis, (C, S) in axis_scores.items()
    }


def axis_spot_frame(
    result: AxisReliabilityResult,
    *,
    pair_label: str,
    sample_role: str,
    sample_file: str,
    spot_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    n = len(result.C)
    identifiers = list(spot_ids) if spot_ids is not None else [str(index) for index in range(n)]
    if len(identifiers) != n:
        raise ValueError("spot_ids must have the same length as the score arrays.")
    return pd.DataFrame({
        "pair_label": pair_label,
        "sample_role": sample_role,
        "sample_file": sample_file,
        "axis": result.axis,
        "spot_id": identifiers,
        "C_input": result.C,
        "S_input": result.S,
        "balance_score_source": result.balance_score_source,
        "balance_score_domain": result.balance_score_domain,
        "balance_B": result.balance_B,
        "activity_C_input": result.activity_C,
        "activity_S_input": result.activity_S,
        "activity_balance": result.activity_balance,
        "activity_score_source": result.activity_score_source,
        "activity_score_domain": result.activity_score_domain,
        "activity_source_transformations": result.activity_source_transformations,
        "activity_source_version": result.activity_source_version,
        "activity_A": result.activity_A,
        "direction_D": result.direction_D,
        "CA_strength": result.ca_strength,
        "CA_fraction": result.ca_fraction,
        "reliability_state": result.reliability_state,
        "reliability_status": result.status,
        "finite_input": result.finite_input,
        "balance_finite_input": result.balance_finite_input,
        "nonnegative_input": result.nonnegative_input,
        "valid_input": result.valid_input,
        "direction_defined": result.direction_defined,
    })


def summarize_axis(result: AxisReliabilityResult) -> dict:
    def median(values: np.ndarray) -> float:
        finite = np.isfinite(values)
        return float(np.median(values[finite])) if finite.any() else np.nan

    n = len(result.C)
    valid_count = int(result.valid_input.sum())
    direction_count = int(result.direction_defined.sum())
    classified = np.fromiter(
        (value is not None for value in result.reliability_state), dtype=bool, count=n
    )
    classified_count = int(classified.sum())
    summary = {
        "axis": result.axis,
        "n_spots": n,
        "n_finite_input": int(result.finite_input.sum()),
        "n_valid_nonnegative_input": valid_count,
        "n_direction_defined": direction_count,
        "finite_input_fraction": float(result.finite_input.mean()) if n else np.nan,
        "valid_input_fraction": float(result.valid_input.mean()) if n else np.nan,
        "direction_defined_fraction": float(result.direction_defined.mean()) if n else np.nan,
        "undefined_fraction": float((~result.direction_defined).mean()) if n else np.nan,
        "negative_input_fraction": float((result.finite_input & ~result.nonnegative_input).mean()) if n else np.nan,
        "nonfinite_input_fraction": float((~result.finite_input).mean()) if n else np.nan,
        "median_balance_B": median(result.balance_B),
        "median_activity_balance": median(result.activity_balance),
        "median_activity_A": median(result.activity_A),
        "median_direction_D": median(result.direction_D),
        "median_CA_strength": median(result.ca_strength),
        "median_CA_fraction": median(result.ca_fraction),
        "activity_min": float(np.nanmin(result.activity_A)) if np.isfinite(result.activity_A).any() else np.nan,
        "activity_max": float(np.nanmax(result.activity_A)) if np.isfinite(result.activity_A).any() else np.nan,
        "classified_spot_count": classified_count,
    }
    for state in (
        "low_activity",
        "c_dominant_active",
        "s_dominant_active",
        "active_coactivation_candidate",
    ):
        count = int(np.sum(result.reliability_state == state))
        summary[f"{state}_count"] = count
        summary[f"{state}_fraction"] = (
            float(count / classified_count) if classified_count else np.nan
        )
    return summary
