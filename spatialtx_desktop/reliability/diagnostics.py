from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd


def _one_dimensional(values: Sequence[float], label: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{label} must be a one-dimensional score array.")
    return array


def _quantile(values: np.ndarray, q: float) -> float:
    finite = values[np.isfinite(values)]
    return float(np.quantile(finite, q)) if len(finite) else np.nan


def score_domain_diagnostic_row(
    *,
    pair_label: str,
    sample_role: str,
    sample_file: str,
    score_role: str,
    C: Sequence[float],
    S: Sequence[float],
    score_source: str,
    score_domain: str,
    transformation_history: str,
    source_version: str,
) -> dict:
    c_values = _one_dimensional(C, "C")
    s_values = _one_dimensional(S, "S")
    if c_values.shape != s_values.shape:
        raise ValueError("C and S diagnostic arrays must have the same shape.")
    total = int(len(c_values))
    finite = np.isfinite(c_values) & np.isfinite(s_values)
    either_negative = finite & ((c_values < 0) | (s_values < 0))
    both_nonnegative = finite & (c_values >= 0) & (s_values >= 0)
    nonfinite = ~finite
    first_negative_stage = (
        "per_gene_zscore_across_spots"
        if score_role == "legacy_balance" and either_negative.any()
        else "activity_source_input_contract_failed"
        if score_role == "activity" and either_negative.any()
        else "not_observed"
    )
    return {
        "pair_label": str(pair_label),
        "sample_role": str(sample_role),
        "sample_file": str(sample_file),
        "score_role": str(score_role),
        "score_source": str(score_source),
        "score_domain": str(score_domain),
        "source_version": str(source_version),
        "transformation_history": str(transformation_history),
        "first_negative_stage": first_negative_stage,
        "total_spots": total,
        "finite_spots": int(finite.sum()),
        "C_min": _quantile(c_values, 0.0),
        "C_q01": _quantile(c_values, 0.01),
        "C_median": _quantile(c_values, 0.5),
        "C_q99": _quantile(c_values, 0.99),
        "C_max": _quantile(c_values, 1.0),
        "S_min": _quantile(s_values, 0.0),
        "S_q01": _quantile(s_values, 0.01),
        "S_median": _quantile(s_values, 0.5),
        "S_q99": _quantile(s_values, 0.99),
        "S_max": _quantile(s_values, 1.0),
        "either_negative_count": int(either_negative.sum()),
        "either_negative_fraction": float(either_negative.mean()) if total else np.nan,
        "both_nonnegative_count": int(both_nonnegative.sum()),
        "both_nonnegative_fraction": float(both_nonnegative.mean()) if total else np.nan,
        "nonfinite_count": int(nonfinite.sum()),
        "nonfinite_fraction": float(nonfinite.mean()) if total else np.nan,
    }


def build_score_domain_diagnostic(
    *,
    pair_label: str,
    sample_role: str,
    sample_file: str,
    legacy_C: Sequence[float],
    legacy_S: Sequence[float],
    activity_C: Sequence[float],
    activity_S: Sequence[float],
    source_metadata: Mapping[str, object],
) -> pd.DataFrame:
    metadata = dict(source_metadata)
    rows = [
        score_domain_diagnostic_row(
            pair_label=pair_label,
            sample_role=sample_role,
            sample_file=sample_file,
            score_role="legacy_balance",
            C=legacy_C,
            S=legacy_S,
            score_source=str(metadata.get("balance_score_source", "legacy_signed_cs")),
            score_domain=str(metadata.get("balance_score_domain", "signed")),
            transformation_history=str(metadata.get("balance_source_transformations", "")),
            source_version="v0.6-preserved",
        ),
        score_domain_diagnostic_row(
            pair_label=pair_label,
            sample_role=sample_role,
            sample_file=sample_file,
            score_role="activity",
            C=activity_C,
            S=activity_S,
            score_source=str(metadata.get("activity_score_source", "unavailable")),
            score_domain=str(metadata.get("activity_score_domain", "unknown")),
            transformation_history=str(metadata.get("activity_source_transformations", "")),
            source_version=str(metadata.get("activity_source_version", "unversioned")),
        ),
    ]
    return pd.DataFrame(rows)
