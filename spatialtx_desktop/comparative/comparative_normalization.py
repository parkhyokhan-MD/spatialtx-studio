from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

from ..workflow import _dense, _gene_indices, _is_count_like
from .metric_registry import METRIC_REGISTRY, metric_definition


PERCENT_REFERENCE_EPSILON = 1e-8
SYMMETRIC_PERCENT_EPSILON = 1e-12
POOLED_HIGH_QUANTILE = 0.90
SCALE_WARNING_BANNER = (
    "Caution: Reference and target differ in valid spot count, tissue extent, or tissue-component count. "
    "Review normalized topology metrics before interpreting raw component counts."
)
CENTERED_HV_WARNING = (
    "H/V centered sample means are expected to be approximately zero after within-sample normalization "
    "and are therefore omitted from comparative interpretation."
)


NORMALIZED_EQUIVALENTS = {
    "n_diffuse_components": "diffuse_components_per_1000_valid_spots",
    "n_small_components": "small_components_per_1000_valid_spots",
    "n_interface_components": "interface_segments_per_1000_valid_spots",
}


def _finite_number(value) -> float:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return float(numeric) if np.isfinite(numeric) else np.nan


def _safe_per_1000(numerator, denominator) -> float:
    numerator = _finite_number(numerator)
    denominator = _finite_number(denominator)
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator <= 0:
        return np.nan
    return float(1000.0 * numerator / denominator)


def _safe_ratio(numerator, denominator) -> float:
    numerator = _finite_number(numerator)
    denominator = _finite_number(denominator)
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator <= 0:
        return np.nan
    return float(numerator / denominator)


def sample_scale_metrics(adata, scored_fields: dict, transition_mask: np.ndarray, graph_result) -> dict:
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    coords = np.asarray(scored_fields["coords"], dtype=float)
    C = np.asarray(scored_fields["C"], dtype=float)
    S = np.asarray(scored_fields["S"], dtype=float)
    R = np.asarray(scored_fields["R"], dtype=float)
    G = np.asarray(scored_fields["G"], dtype=float)
    valid = np.isfinite(C) & np.isfinite(S) & np.isfinite(R) & np.isfinite(G)
    n_total = int(adata.n_obs)
    n_valid = int(valid.sum())
    in_tissue = np.nan
    if "in_tissue" in adata.obs:
        tissue = pd.to_numeric(adata.obs["in_tissue"], errors="coerce")
        if tissue.notna().any():
            in_tissue = int(tissue.eq(1).sum())
    if graph_result is not None and graph_result.connectivities.shape[0] == n_total:
        tissue_components = int(connected_components(graph_result.connectivities > 0, directed=False)[0])
    else:
        tissue_components = np.nan
    extent_x = float(np.ptp(coords[:, 0])) if len(coords) else np.nan
    extent_y = float(np.ptp(coords[:, 1])) if len(coords) else np.nan
    if len(coords) >= 2:
        nearest = np.asarray(cKDTree(coords).query(coords, k=2)[0][:, 1], dtype=float)
        spacing = float(np.nanmean(nearest[np.isfinite(nearest)])) if np.isfinite(nearest).any() else np.nan
    else:
        spacing = np.nan
    return {
        "n_total_spots": n_total,
        "n_valid_spots": n_valid,
        "n_in_tissue_spots": in_tissue,
        "n_transition_spots": int(np.asarray(transition_mask, dtype=bool).sum()),
        "tissue_area_proxy": n_valid,
        "tissue_area_proxy_definition": "n_valid_spots",
        "physical_tissue_area_estimate": np.nan,
        "physical_tissue_area_unit": "unavailable_without_reliable_physical_scaling",
        "tissue_component_count": tissue_components,
        "spatial_extent_x": extent_x,
        "spatial_extent_y": extent_y,
        "spatial_extent_area_proxy": extent_x * extent_y if np.isfinite(extent_x) and np.isfinite(extent_y) else np.nan,
        "spatial_extent_area_proxy_unit": "coordinate_units_squared_not_physical_area",
        "mean_spot_spacing": spacing,
    }


def add_normalized_topology_metrics(metrics: dict) -> dict:
    result = dict(metrics)
    diffuse_components = _finite_number(result.get("n_diffuse_components"))
    interface_components = _finite_number(result.get("n_interface_components"))
    small_fraction = _finite_number(result.get("small_component_fraction"))
    n_small = (
        int(round(diffuse_components * small_fraction))
        if np.isfinite(diffuse_components) and np.isfinite(small_fraction)
        else np.nan
    )
    result["n_small_components"] = n_small
    result["diffuse_components_per_1000_valid_spots"] = _safe_per_1000(
        diffuse_components, result.get("n_valid_spots")
    )
    result["diffuse_components_per_1000_in_tissue_spots"] = _safe_per_1000(
        diffuse_components, result.get("n_in_tissue_spots")
    )
    result["diffuse_components_per_tissue_component"] = _safe_ratio(
        diffuse_components, result.get("tissue_component_count")
    )
    result["small_components_per_1000_valid_spots"] = _safe_per_1000(
        n_small, result.get("n_valid_spots")
    )
    result["transition_components_per_1000_transition_spots"] = _safe_per_1000(
        diffuse_components, result.get("n_transition_spots")
    )
    result["interface_segments_per_1000_valid_spots"] = _safe_per_1000(
        interface_components, result.get("n_valid_spots")
    )
    total_transition_components = (
        diffuse_components + interface_components
        if np.isfinite(diffuse_components) and np.isfinite(interface_components)
        else np.nan
    )
    result["normalized_fragmentation_score"] = _safe_per_1000(
        total_transition_components, result.get("n_valid_spots")
    )
    result["normalized_fragmentation_definition"] = (
        "1000 * (n_diffuse_components + n_interface_components) / max(n_valid_spots, 1)"
    )
    return result


def noncentered_context_values(adata, matched_genes: Iterable[str], expression_scale_guess: str) -> tuple[np.ndarray | None, str]:
    genes = [str(gene) for gene in matched_genes if str(gene).strip()]
    if not genes:
        return None, "unavailable_no_matched_genes"
    if expression_scale_guess not in {"raw_counts", "log1p_normalized"}:
        return None, f"unavailable_expression_scale_{expression_scale_guess or 'unknown'}"
    indices, _matched, _missing = _gene_indices(adata, genes)
    if not indices:
        return None, "unavailable_no_matched_genes"
    matrix = np.asarray(_dense(adata.X[:, indices]), dtype=float)
    if expression_scale_guess == "raw_counts" or _is_count_like(adata.X):
        matrix = np.log1p(np.maximum(matrix, 0.0))
        method = "log1p_counts_then_program_mean"
    else:
        method = "existing_nonnegative_expression_then_program_mean"
    values = np.nanmean(matrix, axis=1)
    return (values if np.isfinite(values).any() else None), method


def raw_context_summary(axis: str, values: np.ndarray | None, transition_mask: np.ndarray, method: str) -> dict:
    if values is None:
        return {
            f"{axis}_raw_available": False,
            f"{axis}_raw_mean": np.nan,
            f"{axis}_raw_median": np.nan,
            f"{axis}_q75": np.nan,
            f"{axis}_q90": np.nan,
            f"{axis}_high_fraction": np.nan,
            f"{axis}_variance": np.nan,
            f"{axis}_MAD": np.nan,
            f"{axis}_transition_median": np.nan,
            f"{axis}_nontransition_median": np.nan,
            f"{axis}_transition_enrichment": np.nan,
            f"{axis}_spatial_variance": np.nan,
            f"{axis}_coefficient_of_variation": np.nan,
            f"{axis}_local_hotspot_fraction": np.nan,
            f"{axis}_pooled_high_threshold": np.nan,
            f"{axis}_high_threshold_method": "unavailable",
            f"{axis}_raw_normalization_method": method,
        }
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    transition = np.asarray(transition_mask, dtype=bool) & finite
    other = ~np.asarray(transition_mask, dtype=bool) & finite
    clean = values[finite]
    transition_median = float(np.nanmedian(values[transition])) if transition.any() else np.nan
    other_median = float(np.nanmedian(values[other])) if other.any() else np.nan
    median = float(np.nanmedian(clean))
    mean = float(np.nanmean(clean))
    nonnegative = bool(np.nanmin(clean) >= 0)
    coefficient = float(np.nanstd(clean) / mean) if nonnegative and abs(mean) > PERCENT_REFERENCE_EPSILON else np.nan
    return {
        f"{axis}_raw_available": True,
        f"{axis}_raw_mean": mean,
        f"{axis}_raw_median": median,
        f"{axis}_q75": float(np.nanquantile(clean, 0.75)),
        f"{axis}_q90": float(np.nanquantile(clean, 0.90)),
        f"{axis}_high_fraction": np.nan,
        f"{axis}_variance": float(np.nanvar(clean)),
        f"{axis}_MAD": float(np.nanmedian(np.abs(clean - median))),
        f"{axis}_transition_median": transition_median,
        f"{axis}_nontransition_median": other_median,
        f"{axis}_transition_enrichment": (
            transition_median - other_median
            if np.isfinite(transition_median) and np.isfinite(other_median)
            else np.nan
        ),
        f"{axis}_spatial_variance": float(np.nanvar(clean)),
        f"{axis}_coefficient_of_variation": coefficient,
        f"{axis}_local_hotspot_fraction": np.nan,
        f"{axis}_pooled_high_threshold": np.nan,
        f"{axis}_high_threshold_method": "pending_pooled_reference_target_q90",
        f"{axis}_raw_normalization_method": method,
    }


def apply_pooled_hv_thresholds(
    sample_metrics: pd.DataFrame,
    fields_by_sample: dict[str, dict[str, np.ndarray]],
    reference: str,
    target: str,
) -> pd.DataFrame:
    result = sample_metrics.copy()
    eligible_ids = set(result.loc[result["group"].isin([reference, target]), "sample_id"].astype(str))
    for axis in ("H", "V"):
        field_name = f"{axis}_expr_raw"
        pooled = [
            np.asarray(fields_by_sample[sample_id][field_name], dtype=float)
            for sample_id in sorted(eligible_ids)
            if sample_id in fields_by_sample and field_name in fields_by_sample[sample_id]
        ]
        finite_pooled = np.concatenate(pooled) if pooled else np.asarray([], dtype=float)
        finite_pooled = finite_pooled[np.isfinite(finite_pooled)]
        threshold = float(np.nanquantile(finite_pooled, POOLED_HIGH_QUANTILE)) if len(finite_pooled) else np.nan
        for index, row in result.iterrows():
            sample_id = str(row["sample_id"])
            fields = fields_by_sample.get(sample_id, {})
            values = np.asarray(fields.get(field_name, []), dtype=float)
            result.at[index, f"{axis}_pooled_high_threshold"] = threshold
            result.at[index, f"{axis}_high_threshold_method"] = (
                "pooled_reference_target_q90" if np.isfinite(threshold) and len(values) else "unavailable_noncentered_scores"
            )
            if not np.isfinite(threshold) or not len(values):
                result.at[index, f"{axis}_high_fraction"] = np.nan
                result.at[index, f"{axis}_local_hotspot_fraction"] = np.nan
                continue
            high = np.isfinite(values) & (values >= threshold)
            result.at[index, f"{axis}_high_fraction"] = float(np.mean(high))
            edge_i = np.asarray(fields.get("context_edge_i", []), dtype=int)
            edge_j = np.asarray(fields.get("context_edge_j", []), dtype=int)
            has_high_neighbor = np.zeros(len(values), dtype=bool)
            if len(edge_i) and len(edge_i) == len(edge_j):
                both = high[edge_i] & high[edge_j]
                has_high_neighbor[edge_i[both]] = True
                has_high_neighbor[edge_j[both]] = True
            result.at[index, f"{axis}_local_hotspot_fraction"] = float(np.mean(high & has_high_neighbor))
    return result


def percent_change(reference: float, target: float) -> tuple[float, str]:
    if not np.isfinite(reference) or not np.isfinite(target):
        return np.nan, "missing"
    if abs(reference) <= PERCENT_REFERENCE_EPSILON:
        return np.nan, "unstable_reference"
    return float(100.0 * (target - reference) / abs(reference)), "ok"


def symmetric_percent_change(reference: float, target: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(target):
        return np.nan
    difference = target - reference
    denominator = abs(target) + abs(reference)
    if denominator <= SYMMETRIC_PERCENT_EPSILON:
        return 0.0 if abs(difference) <= SYMMETRIC_PERCENT_EPSILON else float(
            200.0 * difference / SYMMETRIC_PERCENT_EPSILON
        )
    return float(200.0 * difference / denominator)


def build_metric_change_table(delta_metrics: pd.DataFrame) -> pd.DataFrame:
    table = delta_metrics.copy()
    if table.empty:
        return pd.DataFrame(columns=[
            "metric_name", "display_name", "category", "unit", "reference_value", "target_value",
            "raw_delta", "percent_change", "symmetric_percent_change", "normalization_denominator",
            "normalization_status", "scale_sensitive", "observational_only", "interpretation_flag", "warning",
        ])
    normalized_lookup = {
        (str(row["comparison_id"]), str(row["metric"])): _finite_number(row.get("delta"))
        for _, row in table.iterrows()
    }
    rows: list[dict] = []
    for _, row in table.iterrows():
        name = str(row["metric"])
        definition = metric_definition(name)
        reference = _finite_number(row.get("reference_value"))
        target = _finite_number(row.get("target_value"))
        raw_delta = _finite_number(row.get("delta"))
        ordinary, ordinary_status = percent_change(reference, target)
        symmetric = symmetric_percent_change(reference, target)
        normalized_name = NORMALIZED_EQUIVALENTS.get(name, "")
        normalized_delta = normalized_lookup.get((str(row["comparison_id"]), normalized_name), np.nan)
        warning = ""
        if definition.deprecated:
            interpretation = "non_informative_centered_mean"
            warning = CENTERED_HV_WARNING
        elif not np.isfinite(raw_delta):
            interpretation = "missing"
            warning = "Reference or target value is unavailable."
        elif ordinary_status == "unstable_reference":
            interpretation = "unstable_reference"
            warning = "Ordinary percent change omitted because the reference value is zero or near zero."
        elif definition.scale_sensitive:
            interpretation = "scale_sensitive"
            warning = "Raw value is scale-sensitive; review sample scale and its normalized equivalent when available."
        elif abs(raw_delta) <= 1e-12:
            interpretation = "unchanged"
        else:
            interpretation = "increased" if raw_delta > 0 else "decreased"
        if definition.normalization_denominator:
            normalization_status = f"normalized_by_{definition.normalization_denominator}"
        elif normalized_name and np.isfinite(normalized_delta):
            normalization_status = f"normalized_equivalent:{normalized_name}"
        elif definition.scale_sensitive:
            normalization_status = "normalized_equivalent_unavailable"
        else:
            normalization_status = "not_required"
        rows.append({
            "comparison_id": row.get("comparison_id", ""),
            "pair_id": row.get("pair_id", ""),
            "reference_sample_id": row.get("reference_sample_id", ""),
            "target_sample_id": row.get("target_sample_id", ""),
            "metric_name": name,
            "delta_metric": row.get("delta_metric", definition.resolved_delta_name),
            "display_name": definition.display_name,
            "category": definition.category,
            "unit": definition.unit,
            "reference_value": reference,
            "target_value": target,
            "raw_delta": raw_delta,
            "normalized_delta": normalized_delta,
            "percent_change": ordinary,
            "symmetric_percent_change": symmetric,
            "normalization_denominator": definition.normalization_denominator,
            "normalization_status": normalization_status,
            "scale_sensitive": bool(definition.scale_sensitive),
            "observational_only": bool(definition.observational_only),
            "deprecated": bool(definition.deprecated),
            "interpretation_priority": int(definition.interpretation_priority),
            "plot_group": definition.plot_group,
            "interpretation_flag": interpretation,
            "warning": warning,
            "direction_definition": row.get("direction_definition", "Target - Reference"),
            "status": row.get("status", ""),
        })
    return pd.DataFrame(rows)


def build_sample_scale_table(sample_metrics: pd.DataFrame) -> pd.DataFrame:
    identifiers = [column for column in ("sample_id", "group", "pair_id", "condition", "batch") if column in sample_metrics]
    metrics = [
        definition.internal_name
        for definition in METRIC_REGISTRY
        if definition.plot_group == "sample_scale" and definition.internal_name in sample_metrics
    ]
    extras = [
        column for column in (
            "tissue_area_proxy_definition", "physical_tissue_area_estimate", "physical_tissue_area_unit",
            "spatial_extent_area_proxy_unit",
        ) if column in sample_metrics
    ]
    return sample_metrics[identifiers + metrics + extras].copy()


def build_normalized_metrics_table(sample_metrics: pd.DataFrame) -> pd.DataFrame:
    identifiers = [column for column in ("sample_id", "group", "pair_id", "condition", "batch") if column in sample_metrics]
    metrics = [
        definition.internal_name
        for definition in METRIC_REGISTRY
        if definition.plot_group == "topology_normalized" and definition.internal_name in sample_metrics
    ]
    extras = [column for column in ("normalized_fragmentation_definition",) if column in sample_metrics]
    return sample_metrics[identifiers + metrics + extras].copy()


def build_relative_changes_table(metric_changes: pd.DataFrame) -> pd.DataFrame:
    if metric_changes.empty:
        return metric_changes.copy()
    return metric_changes.loc[
        ~metric_changes["deprecated"].astype(bool),
        [
            "comparison_id", "metric_name", "display_name", "category", "unit", "reference_value",
            "target_value", "raw_delta", "percent_change", "symmetric_percent_change",
            "interpretation_flag", "warning",
        ],
    ].copy()


def _value_for(table: pd.DataFrame, comparison_id: str, metric: str, column: str) -> float:
    selected = table.loc[
        table["comparison_id"].astype(str).eq(comparison_id) & table["metric_name"].eq(metric), column
    ]
    return _finite_number(selected.iloc[0]) if len(selected) else np.nan


def _scale_ratio(left: float, right: float) -> float:
    if not np.isfinite(left) or not np.isfinite(right) or min(abs(left), abs(right)) <= 0:
        return np.nan
    return float(max(abs(left), abs(right)) / min(abs(left), abs(right)))


def build_scale_warnings(metric_changes: pd.DataFrame) -> pd.DataFrame:
    columns = ["comparison_id", "warning_code", "severity", "metrics", "message"]
    if metric_changes.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict] = []
    for comparison_id in metric_changes["comparison_id"].astype(str).drop_duplicates():
        def values(metric: str) -> tuple[float, float]:
            return (
                _value_for(metric_changes, comparison_id, metric, "reference_value"),
                _value_for(metric_changes, comparison_id, metric, "target_value"),
            )

        affected: list[str] = []
        ref_valid, tar_valid = values("n_valid_spots")
        if _scale_ratio(ref_valid, tar_valid) > 1.5:
            affected.append("n_valid_spots")
        ref_area, tar_area = values("tissue_area_proxy")
        if _scale_ratio(ref_area, tar_area) > 1.5:
            affected.append("tissue_area_proxy")
        ref_extent, tar_extent = values("spatial_extent_area_proxy")
        if _scale_ratio(ref_extent, tar_extent) > 1.5:
            affected.append("spatial_extent_area_proxy")
        ref_components, tar_components = values("tissue_component_count")
        if np.isfinite(ref_components) and np.isfinite(tar_components) and abs(tar_components - ref_components) > 2:
            affected.append("tissue_component_count")
        if affected:
            rows.append({
                "comparison_id": comparison_id,
                "warning_code": "substantial_sample_scale_difference",
                "severity": "caution",
                "metrics": ";".join(dict.fromkeys(affected)),
                "message": SCALE_WARNING_BANNER,
            })
        ref_tissue, tar_tissue = values("n_in_tissue_spots")
        if not np.isfinite(ref_tissue) or not np.isfinite(tar_tissue):
            rows.append({
                "comparison_id": comparison_id,
                "warning_code": "in_tissue_count_unavailable",
                "severity": "information",
                "metrics": "n_in_tissue_spots",
                "message": "In-tissue spot counts were unavailable for at least one sample; valid-spot normalization remains available.",
            })
        for raw_name, normalized_name in NORMALIZED_EQUIVALENTS.items():
            raw_change = _value_for(metric_changes, comparison_id, raw_name, "symmetric_percent_change")
            normalized_change = _value_for(metric_changes, comparison_id, normalized_name, "symmetric_percent_change")
            if np.isfinite(raw_change) and np.isfinite(normalized_change) and abs(raw_change) >= 50 and abs(normalized_change) <= 20:
                rows.append({
                    "comparison_id": comparison_id,
                    "warning_code": "raw_normalized_change_discordance",
                    "severity": "caution",
                    "metrics": f"{raw_name};{normalized_name}",
                    "message": (
                        f"Raw {metric_definition(raw_name).display_name.lower()} changed substantially, but "
                        f"{metric_definition(normalized_name).display_name.lower()} changed only modestly. "
                        "The raw difference may partly reflect tissue size or valid-spot differences."
                    ),
                })
    return pd.DataFrame(rows, columns=columns)


def build_hv_summary(metric_changes: pd.DataFrame, sample_metrics: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "comparison_id", "axis", "metric", "reference_value", "target_value", "raw_delta",
        "symmetric_percent_change", "threshold_value", "threshold_method", "normalization_method",
        "observational_only", "warning",
    ]
    rows: list[dict] = []
    selected_names = {
        "H_raw_mean", "H_raw_median", "H_q75", "H_q90", "H_high_fraction", "H_variance", "H_MAD",
        "H_transition_enrichment", "H_spatial_variance", "H_coefficient_of_variation", "H_local_hotspot_fraction",
        "V_raw_mean", "V_raw_median", "V_q75", "V_q90", "V_high_fraction", "V_variance", "V_MAD",
        "V_transition_enrichment", "V_spatial_variance", "V_coefficient_of_variation", "V_local_hotspot_fraction",
        "H_expr_mean", "V_expr_mean",
    }
    for _, row in metric_changes.loc[metric_changes["metric_name"].isin(selected_names)].iterrows():
        name = str(row["metric_name"])
        axis = name[0]
        threshold_values = pd.to_numeric(
            sample_metrics.get(f"{axis}_pooled_high_threshold", pd.Series(index=sample_metrics.index, dtype=float)),
            errors="coerce",
        ).dropna()
        methods = sample_metrics.get(f"{axis}_high_threshold_method", pd.Series(dtype=str)).astype(str)
        normalizations = sample_metrics.get(f"{axis}_raw_normalization_method", pd.Series(dtype=str)).astype(str)
        deprecated = bool(row.get("deprecated", False))
        warning = CENTERED_HV_WARNING if deprecated else str(row.get("warning", ""))
        if not deprecated and not np.isfinite(_finite_number(row.get("raw_delta"))):
            warning = warning or "Non-centered H/V scores were unavailable; this observational metric was omitted from figures."
        rows.append({
            "comparison_id": row.get("comparison_id", ""),
            "axis": axis,
            "metric": name,
            "reference_value": row.get("reference_value", np.nan),
            "target_value": row.get("target_value", np.nan),
            "raw_delta": row.get("raw_delta", np.nan),
            "symmetric_percent_change": row.get("symmetric_percent_change", np.nan),
            "threshold_value": float(threshold_values.iloc[0]) if len(threshold_values) else np.nan,
            "threshold_method": next((value for value in methods if value and value != "nan"), "unavailable"),
            "normalization_method": next((value for value in normalizations if value and value != "nan"), "unavailable"),
            "observational_only": True,
            "warning": warning,
        })
    return pd.DataFrame(rows, columns=columns)
