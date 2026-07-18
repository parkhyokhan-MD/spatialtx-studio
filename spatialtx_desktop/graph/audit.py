from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from .. import __version__
from .metadata import json_safe


AUDIT_SCHEMA_VERSION = "0.4"


def _summary(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {"minimum": np.nan, "q05": np.nan, "median": np.nan, "mean": np.nan, "q95": np.nan, "maximum": np.nan}
    return {
        "minimum": float(np.min(values)),
        "q05": float(np.quantile(values, 0.05)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "q95": float(np.quantile(values, 0.95)),
        "maximum": float(np.max(values)),
    }


def _matrix_diagnostics(X, n_obs: int, n_vars: int) -> tuple[str, str, np.ndarray, np.ndarray, float]:
    if sparse.issparse(X):
        matrix = X.tocsr()
        totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
        detected = np.asarray(matrix.getnnz(axis=1)).ravel().astype(float)
        data = np.asarray(matrix.data, dtype=float)
        zero_fraction = 1.0 - (matrix.nnz / max(1, n_obs * n_vars))
        storage = "sparse"
    else:
        matrix = np.asarray(X)
        totals = np.nansum(matrix, axis=1).astype(float) if matrix.ndim == 2 else np.asarray([], dtype=float)
        detected = np.sum(np.isfinite(matrix) & (matrix != 0), axis=1).astype(float) if matrix.ndim == 2 else np.asarray([], dtype=float)
        data = np.asarray(matrix, dtype=float).ravel()
        zero_fraction = float(np.mean(data == 0)) if len(data) else np.nan
        storage = "dense"
    finite = data[np.isfinite(data)]
    nonzero = finite[finite != 0]
    if not len(finite):
        state = "unknown"
    elif np.nanmin(finite) < 0:
        state = "centered_or_scaled"
    else:
        probe = nonzero if len(nonzero) else finite
        if len(probe) > 100000:
            step = max(1, len(probe) // 100000)
            probe = probe[::step][:100000]
        integer_fraction = float(np.mean(np.isclose(probe, np.round(probe), atol=1e-8))) if len(probe) else 0.0
        if integer_fraction >= 0.995:
            state = "raw_counts"
        elif float(np.nanmax(finite)) <= 30 and np.nanmin(finite) >= 0:
            state = "log1p_normalized"
        else:
            state = "normalized_unknown_method"
    return storage, state, totals, detected, float(zero_fraction)


def _spatial_payload(adata) -> tuple[bool, bool]:
    spatial = adata.uns.get("spatial", {}) if hasattr(adata, "uns") else {}
    scalefactors = False
    images = False
    if isinstance(spatial, dict):
        for value in spatial.values():
            if isinstance(value, dict):
                scalefactors = scalefactors or bool(value.get("scalefactors"))
                images = images or bool(value.get("images"))
    return scalefactors, images


def _platform_guess(adata, has_spatial: bool, scalefactors: bool) -> str:
    obs_columns = {str(column).lower() for column in getattr(adata, "obs", pd.DataFrame()).columns}
    uns_keys = {str(key).lower() for key in getattr(adata, "uns", {})}
    import_meta = getattr(adata, "uns", {}).get("spatialtx_import", {})
    source_text = str(import_meta.get("source_format", "")).lower() if isinstance(import_meta, dict) else ""
    if {"array_row", "array_col"}.issubset(obs_columns) or scalefactors or "visium" in source_text:
        return "visium"
    imaging_terms = {"cell_id", "x_centroid", "y_centroid", "transcript_count"}
    if obs_columns & imaging_terms or uns_keys & {"xenium", "cosmx", "merscope", "merfish"}:
        return "imaging_or_cell_based"
    if has_spatial:
        return "generic_spot_based"
    return "unknown"


def audit_input(
    adata,
    *,
    input_filename: str | Path | None = None,
    coordinate_unit: str = "native",
    coordinate_scale: float | None = None,
    scale_source: str = "",
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Audit one AnnData object without normalizing or mutating its expression matrix."""
    n_obs, n_vars = int(adata.n_obs), int(adata.n_vars)
    storage, expression_state, totals, detected, zero_fraction = _matrix_diagnostics(adata.X, n_obs, n_vars)
    coords = None
    coordinate_source = "unavailable"
    if hasattr(adata, "obsm") and "spatial" in adata.obsm:
        try:
            coords = np.asarray(adata.obsm["spatial"], dtype=float)
            coordinate_source = 'adata.obsm["spatial"]'
        except Exception:
            coords = None
    missing_coordinates = coords is None or coords.ndim != 2 or coords.shape != (n_obs, 2) or coords.size == 0
    nonfinite_coordinates = bool(not missing_coordinates and not np.isfinite(coords).all())
    duplicate_coordinates = 0
    if not missing_coordinates and not nonfinite_coordinates:
        duplicate_coordinates = int(n_obs - len(np.unique(coords, axis=0)))
    scalefactors, images = _spatial_payload(adata)
    obs = getattr(adata, "obs", pd.DataFrame(index=np.arange(n_obs)))
    tissue_mask_available = "in_tissue" in obs
    tissue_count = None
    if tissue_mask_available:
        tissue_values = pd.to_numeric(obs["in_tissue"], errors="coerce")
        tissue_count = int((tissue_values == 1).sum())
    platform = _platform_guess(adata, not missing_coordinates, scalefactors)
    if expression_state in {"raw_counts", "log1p_normalized"}:
        detection_interpretation = "expression_greater_than_zero_is_detection_fraction"
    elif expression_state == "centered_or_scaled":
        detection_interpretation = "detection_fraction_unavailable_positive_value_fraction_only"
    else:
        detection_interpretation = "detection_fraction_requires_verified_nonnegative_scale"
    unit = str(coordinate_unit or "native").strip().lower()
    calibrated = bool(str(scale_source).strip()) and (
        unit in {"micrometer", "micrometre", "um", "µm"}
        or (coordinate_scale is not None and float(coordinate_scale) > 0)
    )
    if calibrated:
        unit = "micrometer"
    warnings: list[str] = []
    if expression_state == "unknown":
        warnings.append("expression preprocessing state could not be determined")
    warnings.append("cross-sample score magnitudes may not be directly comparable")
    if missing_coordinates:
        warnings.append("spatial coordinates are unavailable")
    if unit in {"native", "unknown", ""}:
        warnings.append("coordinate unit is unknown")
    if nonfinite_coordinates or duplicate_coordinates:
        warnings.append("duplicate or non-finite coordinates were detected")
    if not tissue_mask_available:
        warnings.append("tissue mask is unavailable")
        if not missing_coordinates:
            warnings.append("spatial analysis may include off-tissue observations")
    var_names = [str(value) for value in adata.var_names]
    obs_names = [str(value) for value in adata.obs_names]
    audit = {
        "application_version": __version__,
        "audit_schema_version": AUDIT_SCHEMA_VERSION,
        "input_filename": str(Path(input_filename).resolve()) if input_filename else "",
        "platform_guess": platform,
        "expression_state_guess": expression_state,
        "expression_scale_guess": expression_state,
        "detection_metric_interpretation": detection_interpretation,
        "n_obs": n_obs,
        "n_vars": n_vars,
        "sparse_or_dense_X": storage,
        "coordinate_source": coordinate_source,
        "coordinate_unit": unit or "native",
        "coordinate_scale": float(coordinate_scale) if coordinate_scale is not None else None,
        "scale_source": str(scale_source or ""),
        "physical_calibration_available": bool(calibrated),
        "array_row_column_available": bool({"array_row", "array_col"}.issubset(obs.columns)),
        "tissue_mask_available": bool(tissue_mask_available),
        "scalefactor_available": bool(scalefactors),
        "histology_image_available": bool(images),
        "duplicate_gene_names": int(len(var_names) - len(set(var_names))),
        "duplicate_observation_names": int(len(obs_names) - len(set(obs_names))),
        "missing_coordinates": bool(missing_coordinates),
        "non_finite_coordinates": bool(nonfinite_coordinates),
        "duplicate_coordinates": int(duplicate_coordinates),
        "total_count_distribution": _summary(totals),
        "detected_gene_distribution": _summary(detected),
        "zero_expression_fraction": zero_fraction,
        "tissue_spot_count": tissue_count,
        "warnings": warnings,
        "matrix_altered": False,
    }
    adata.uns["spatialtx_input_audit"] = json_safe(audit)
    row = dict(audit)
    row["warnings"] = "; ".join(warnings)
    row["total_count_distribution"] = json_safe(audit["total_count_distribution"])
    row["detected_gene_distribution"] = json_safe(audit["detected_gene_distribution"])
    return json_safe(audit), pd.DataFrame([row])
