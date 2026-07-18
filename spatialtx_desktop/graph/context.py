from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

from .. import __version__
from ..workflow import _dense, _gene_indices, _is_count_like, parse_gene_text
from .metadata import SMOOTHING_LIMITATION, json_safe


DEFAULT_HYPOXIA_GENES = [
    "CA9", "VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "ALDOA", "BNIP3", "NDRG1", "ADM", "EGLN3", "P4HA1"
]

DEFAULT_VASCULAR_PROXY_GENES = [
    "PECAM1", "VWF", "KDR", "FLT1", "ESAM", "ENG", "CDH5", "RAMP2", "EMCN", "CLDN5", "ANGPT2", "PLVAP"
]


@dataclass
class ContextFieldConfig:
    field: str
    genes: Iterable[str] | None = None
    gene_set_name: str | None = None
    score_method: str = "z_score_mean"
    min_coverage: float = 0.25
    allow_low_coverage: bool = False
    smoothing: str = "none"
    high_quantile: float = 0.80
    min_spot_fraction: float = 0.01
    dominant_gene_warning_fraction: float = 0.60
    library_correlation_warning: float = 0.80


def _rank_quantile(values: np.ndarray) -> np.ndarray:
    result = np.full(len(values), np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    ranks = pd.Series(values[finite]).rank(method="average").to_numpy(dtype=float)
    result[finite] = (ranks - 1.0) / max(len(ranks) - 1, 1)
    return result


def _zscore_columns(matrix: np.ndarray) -> np.ndarray:
    mean = np.nanmean(matrix, axis=0)
    sd = np.nanstd(matrix, axis=0)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return (matrix - mean) / sd


def _score_matrix(matrix: np.ndarray, method: str) -> np.ndarray:
    method = (method or "z_score_mean").lower()
    if matrix.ndim == 1:
        matrix = matrix[:, None]
    if method == "raw_mean":
        return np.nanmean(matrix, axis=1)
    if method == "z_score_mean":
        return np.nanmean(_zscore_columns(matrix), axis=1)
    if method == "rank_quantile":
        return _rank_quantile(np.nanmean(_zscore_columns(matrix), axis=1))
    raise ValueError("context score method must be raw_mean, z_score_mean, or rank_quantile")


def _smooth(values: np.ndarray, connectivities) -> np.ndarray:
    from scipy import sparse

    W = connectivities.tocsr().astype(float)
    W = W + sparse.eye(W.shape[0], format="csr")
    row_sum = np.asarray(W.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    return np.asarray(W.dot(values) / row_sum, dtype=float).ravel()


def _default_genes(field: str) -> tuple[str, list[str]]:
    if field == "H":
        return "default_hypoxia_associated_expression", list(DEFAULT_HYPOXIA_GENES)
    if field == "V":
        return "default_endothelial_angiogenic_expression_proxy", list(DEFAULT_VASCULAR_PROXY_GENES)
    raise ValueError("field must be 'H' or 'V'")


def _context_label(field: str) -> str:
    return "H_expr" if field == "H" else "V_expr"


def _safe_corr(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    valid = np.isfinite(left) & np.isfinite(right)
    if valid.sum() < 3 or np.nanstd(left[valid]) == 0 or np.nanstd(right[valid]) == 0:
        return np.nan
    value = np.corrcoef(left[valid], right[valid])[0, 1]
    return float(value) if np.isfinite(value) else np.nan


def _library_metrics(X) -> tuple[np.ndarray, np.ndarray]:
    if hasattr(X, "getnnz"):
        return np.asarray(X.sum(axis=1)).ravel().astype(float), np.asarray(X.getnnz(axis=1)).ravel().astype(float)
    matrix = np.asarray(X, dtype=float)
    return np.nansum(matrix, axis=1), np.sum(np.isfinite(matrix) & (matrix != 0), axis=1).astype(float)


def _expression_scale_guess(X) -> str:
    data = np.asarray(X.data if hasattr(X, "data") and hasattr(X, "getnnz") else X, dtype=float).ravel()
    finite = data[np.isfinite(data)]
    if not len(finite):
        return "unknown"
    if float(np.min(finite)) < 0:
        return "centered_or_scaled"
    probe = finite[finite != 0]
    if len(probe) > 100000:
        probe = probe[:: max(1, len(probe) // 100000)][:100000]
    integer_fraction = float(np.mean(np.isclose(probe, np.round(probe), atol=1e-8))) if len(probe) else 1.0
    if integer_fraction >= 0.995:
        return "raw_counts"
    if float(np.max(finite)) <= 30:
        return "log1p_normalized"
    return "normalized_unknown_method"


def _detection_matrix(adata, indices: list[int], matched: list[str], fallback: np.ndarray) -> tuple[np.ndarray, str, str]:
    layers = getattr(adata, "layers", {})
    for layer_name in ("counts", "raw"):
        if layer_name in layers:
            matrix = _dense(layers[layer_name][:, indices])
            return np.asarray(matrix, dtype=float), f"adata.layers[{layer_name!r}]", "raw_counts"
    raw = getattr(adata, "raw", None)
    if raw is not None:
        lookup = {str(name).upper(): i for i, name in enumerate(raw.var_names)}
        raw_indices = [lookup.get(str(gene).upper()) for gene in matched]
        if all(index is not None for index in raw_indices):
            matrix = _dense(raw.X[:, [int(index) for index in raw_indices]])
            return np.asarray(matrix, dtype=float), "adata.raw.X", _expression_scale_guess(raw.X)
    audit = dict(getattr(adata, "uns", {}).get("spatialtx_input_audit", {}))
    scale = str(audit.get("expression_scale_guess") or audit.get("expression_state_guess") or "")
    return np.asarray(fallback, dtype=float), "adata.X", scale or _expression_scale_guess(adata.X)


def _leave_one_gene_out(matrix: np.ndarray, matched: list[str], method: str, full_field: np.ndarray) -> list[dict]:
    rows: list[dict] = []
    for index, gene in enumerate(matched):
        if matrix.shape[1] <= 1:
            without = np.full(len(full_field), np.nan, dtype=float)
        else:
            without = _score_matrix(np.delete(matrix, index, axis=1), method)
        rows.append({
            "gene": gene,
            "field_correlation_without_gene": _safe_corr(full_field, without),
            "mean_absolute_field_change": (
                float(np.nanmean(np.abs(full_field - without))) if np.isfinite(without).any() else np.nan
            ),
            "variance_change": (
                float(np.nanvar(without) - np.nanvar(full_field)) if np.isfinite(without).any() else np.nan
            ),
        })
    ordered = sorted(
        range(len(rows)),
        key=lambda idx: (
            -rows[idx]["mean_absolute_field_change"]
            if np.isfinite(rows[idx]["mean_absolute_field_change"])
            else float("inf")
        ),
    )
    for rank, index in enumerate(ordered, 1):
        rows[index]["rank"] = rank
    return rows


def add_context_field(
    adata,
    config: ContextFieldConfig,
    graph_connectivities=None,
    *,
    reference_fields: dict[str, np.ndarray] | None = None,
    active_graph: str = "",
) -> tuple[pd.DataFrame, dict]:
    """Add H_expr or V_expr to adata.obs as descriptive context fields only."""
    field = config.field.upper()
    default_name, default_genes = _default_genes(field)
    requested = parse_gene_text(config.genes if config.genes is not None else default_genes)
    gene_set_name = config.gene_set_name or default_name
    label = _context_label(field)
    if not requested:
        raise ValueError(f"{label} gene set is empty.")
    indices, matched, missing = _gene_indices(adata, requested)
    coverage = len(matched) / len(requested) if requested else 0.0
    warnings: list[str] = []
    if coverage < float(config.min_coverage):
        warning = f"low {label} gene coverage: {coverage:.1%}"
        warnings.append(warning)
        if not config.allow_low_coverage:
            raise ValueError(warning)
    if not len(indices):
        raise ValueError(f"No requested genes were found for {label}.")
    matrix = np.asarray(_dense(adata.X[:, indices]), dtype=float)
    detection_matrix, detection_source, expression_scale_guess = _detection_matrix(adata, indices, matched, matrix)
    positive_value_fraction = np.mean(np.isfinite(detection_matrix) & (detection_matrix > 0), axis=0)
    detection_is_interpretable = expression_scale_guess in {"raw_counts", "log1p_normalized"}
    detection_fraction = (
        positive_value_fraction.copy()
        if detection_is_interpretable
        else np.full(len(matched), np.nan, dtype=float)
    )
    detection_interpretation = (
        "expression_greater_than_zero_is_detection_fraction"
        if detection_is_interpretable
        else "detection_fraction_unavailable_positive_value_fraction_only"
    )
    expressed = positive_value_fraction >= float(config.min_spot_fraction)
    expressed_fraction = float(np.mean(expressed)) if len(expressed) else 0.0
    if expressed_fraction < 0.5:
        warnings.append(f"most matched genes are not detectably expressed for {label}")
    if _is_count_like(adata.X):
        matrix = np.log1p(np.maximum(matrix, 0.0))
    values = _score_matrix(np.asarray(matrix, dtype=float), config.score_method)
    leave_one_out = _leave_one_gene_out(matrix, matched, config.score_method, values)
    base_column = f"spatialtx_{label}"
    adata.obs[base_column] = values
    active_values = values
    smoothed_column = ""
    if config.smoothing != "none":
        if graph_connectivities is None:
            raise ValueError(f"{label} smoothing requires an active spatial graph.")
        smoothed = _smooth(values, graph_connectivities)
        smoothed_column = f"spatialtx_{label}_smoothed"
        adata.obs[smoothed_column] = smoothed
        active_values = smoothed
        warnings.append(SMOOTHING_LIMITATION)
    threshold = float(np.nanquantile(active_values[np.isfinite(active_values)], float(config.high_quantile))) if np.isfinite(active_values).any() else np.nan
    high_column = f"spatialtx_{label.replace('_expr', '')}_high"
    adata.obs[high_column] = np.asarray(active_values >= threshold, dtype=bool) if np.isfinite(threshold) else False
    high_fraction = float(np.mean(adata.obs[high_column]))
    field_q05, field_q95 = np.nanquantile(active_values, [0.05, 0.95]) if np.isfinite(active_values).any() else (np.nan, np.nan)
    dynamic_range = float(field_q95 - field_q05) if np.isfinite(field_q05) and np.isfinite(field_q95) else np.nan
    if np.isfinite(dynamic_range) and dynamic_range <= 1e-6:
        warnings.append(f"{label} has very low dynamic range")
    if high_fraction < 0.02 or high_fraction > 0.50:
        warnings.append(f"{label} high-state fraction is unexpectedly small or large: {high_fraction:.1%}")
    contribution = np.nanmean(np.abs(_zscore_columns(matrix)), axis=0)
    contribution_sum = float(np.nansum(contribution))
    dominant_index = int(np.nanargmax(contribution)) if len(contribution) and np.isfinite(contribution).any() else -1
    dominant_fraction = float(contribution[dominant_index] / contribution_sum) if dominant_index >= 0 and contribution_sum > 0 else np.nan
    dominant_gene = matched[dominant_index] if dominant_index >= 0 else ""
    if np.isfinite(dominant_fraction) and dominant_fraction >= float(config.dominant_gene_warning_fraction):
        warnings.append(f"one gene dominates {label}: {dominant_gene} ({dominant_fraction:.1%})")
    total_counts, detected_genes = _library_metrics(adata.X)
    corr_total = _safe_corr(active_values, total_counts)
    corr_detected = _safe_corr(active_values, detected_genes)
    if np.isfinite(corr_total) and abs(corr_total) >= float(config.library_correlation_warning):
        warnings.append(f"{label} is strongly correlated with library size: r={corr_total:.3f}")
    reference_fields = reference_fields or {}
    coverage_table = pd.DataFrame({
        "field": label,
        "gene": requested,
        "status": ["matched" if gene.upper() in {m.upper() for m in matched} else "missing" for gene in requested],
    })
    detection_lookup = {str(gene).upper(): float(value) for gene, value in zip(matched, detection_fraction)}
    positive_lookup = {str(gene).upper(): float(value) for gene, value in zip(matched, positive_value_fraction)}
    leave_one_out_lookup = {str(row["gene"]).upper(): row for row in leave_one_out}
    coverage_table["detection_fraction"] = [detection_lookup.get(gene.upper(), np.nan) for gene in requested]
    coverage_table["positive_value_fraction"] = [positive_lookup.get(gene.upper(), np.nan) for gene in requested]
    coverage_table["detection_metric_interpretation"] = detection_interpretation
    coverage_table["detection_source"] = detection_source
    coverage_table["expressed_above_min_spot_fraction"] = coverage_table["positive_value_fraction"] >= float(config.min_spot_fraction)
    for column in ("field_correlation_without_gene", "mean_absolute_field_change", "variance_change", "rank"):
        coverage_table[column] = [leave_one_out_lookup.get(gene.upper(), {}).get(column, np.nan) for gene in requested]
    metadata = {
        "field": label,
        "public_term": "hypoxia-associated expression field" if field == "H" else "endothelial/angiogenic expression proxy",
        "gene_set_name": gene_set_name,
        "requested_genes": requested,
        "matched_genes": matched,
        "missing_genes": missing,
        "requested_gene_count": len(requested),
        "matched_gene_count": len(matched),
        "coverage_fraction": coverage,
        "coverage_warning": "; ".join(warnings),
        "warnings": warnings,
        "genes_expressed_above_min_spot_fraction": [gene for gene, keep in zip(matched, expressed) if keep],
        "expressed_gene_fraction": expressed_fraction,
        "per_gene_detection_fraction": detection_lookup,
        "per_gene_positive_value_fraction": positive_lookup,
        "expression_scale_guess": expression_scale_guess,
        "detection_source": detection_source,
        "detection_metric_interpretation": detection_interpretation,
        "minimum_spot_fraction": float(config.min_spot_fraction),
        "field_zero_fraction": float(np.mean(np.isclose(active_values, 0.0))) if len(active_values) else np.nan,
        "field_minimum": float(np.nanmin(active_values)) if np.isfinite(active_values).any() else np.nan,
        "field_maximum": float(np.nanmax(active_values)) if np.isfinite(active_values).any() else np.nan,
        "field_median": float(np.nanmedian(active_values)) if np.isfinite(active_values).any() else np.nan,
        "field_dynamic_range": dynamic_range,
        "score_method": config.score_method,
        "smoothing_method": config.smoothing,
        "smoothing_warning": SMOOTHING_LIMITATION if config.smoothing != "none" else "",
        "active_graph_used_for_smoothing": active_graph if config.smoothing != "none" else "",
        "high_state_threshold": threshold,
        "high_state_quantile": float(config.high_quantile),
        "high_state_fraction": high_fraction,
        "correlation_with_total_counts": corr_total,
        "correlation_with_detected_gene_count": corr_detected,
        "correlation_with_C": _safe_corr(active_values, reference_fields["C"]) if "C" in reference_fields else np.nan,
        "correlation_with_S": _safe_corr(active_values, reference_fields["S"]) if "S" in reference_fields else np.nan,
        "correlation_with_R": _safe_corr(active_values, reference_fields["R"]) if "R" in reference_fields else np.nan,
        "dominant_gene": dominant_gene,
        "dominant_gene_contribution_fraction": dominant_fraction,
        "dominant_gene_metric_status": "experimental_legacy_absolute_z_contribution",
        "leave_one_gene_out_qc": leave_one_out,
        "version": __version__,
        "interpretation_limit": (
            "Descriptive expression-derived context field only; it does not modify R(x), Type A/B/C calls, or transition masks."
            if field == "H"
            else "Expression proxy only; not true vascularity, perfusion, vessel density, or functional blood supply."
        ),
        "obs_columns": {
            "base": base_column,
            "smoothed": smoothed_column,
            "high": high_column,
        },
    }
    context_uns = dict(adata.uns.get("spatialtx_context_fields", {}))
    uns_metadata = json_safe(metadata)
    uns_metadata["leave_one_gene_out_qc_json"] = json.dumps(
        uns_metadata.pop("leave_one_gene_out_qc", []), ensure_ascii=False
    )
    context_uns[label] = uns_metadata
    adata.uns["spatialtx_context_fields"] = context_uns
    return coverage_table, json_safe(metadata)
