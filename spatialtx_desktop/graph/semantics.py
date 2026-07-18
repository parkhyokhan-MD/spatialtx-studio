from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .metadata import json_safe


VARIABLE_MODES = {
    "categorical_state",
    "binary_mask",
    "continuous_score",
    "proportion_composition",
}


def infer_variable_mode(values: pd.Series) -> str:
    valid = values.dropna()
    if pd.api.types.is_bool_dtype(valid.dtype):
        return "binary_mask"
    numeric = pd.to_numeric(valid, errors="coerce")
    numeric_complete = len(valid) == 0 or numeric.notna().all()
    if numeric_complete:
        unique = pd.unique(numeric)
        if len(unique) <= 2 and set(float(value) for value in unique).issubset({0.0, 1.0}):
            return "binary_mask"
        if len(numeric) and float(numeric.min()) >= 0 and float(numeric.max()) <= 1:
            return "proportion_composition"
        return "continuous_score"
    return "categorical_state"


def describe_obs_variable(
    adata,
    column: str,
    *,
    confirmed_mode: str | None = None,
    missing_handling: str = "exclude",
) -> dict[str, Any]:
    if column not in adata.obs:
        raise ValueError(f"adata.obs column not found: {column}")
    values = adata.obs[column]
    inferred = infer_variable_mode(values)
    mode = confirmed_mode or inferred
    if mode not in VARIABLE_MODES:
        raise ValueError("analysis mode must be categorical_state, binary_mask, continuous_score, or proportion_composition")
    valid = values.dropna()
    metadata: dict[str, Any] = {
        "source_column": column,
        "inferred_data_type": inferred,
        "user_confirmed_analysis_mode": mode,
        "missing_value_handling": missing_handling,
        "number_of_valid_observations": int(len(valid)),
        "missing_value_count": int(values.isna().sum()),
        "public_unit": "spot-level composition" if mode == "proportion_composition" else "spot-level state" if mode in {"categorical_state", "binary_mask"} else "continuous spot-level score",
    }
    if mode in {"categorical_state", "binary_mask"}:
        metadata["category_counts"] = {str(key): int(value) for key, value in valid.astype(str).value_counts().items()}
        metadata["value_range"] = None
    else:
        numeric = pd.to_numeric(valid, errors="coerce").dropna()
        metadata["value_range"] = [float(numeric.min()), float(numeric.max())] if len(numeric) else [np.nan, np.nan]
        metadata["category_counts"] = None
    existing = dict(adata.uns.get("spatialtx_variable_semantics", {}))
    existing[column] = json_safe(metadata)
    adata.uns["spatialtx_variable_semantics"] = existing
    return json_safe(metadata)


def semantics_table(records: list[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame(records)
