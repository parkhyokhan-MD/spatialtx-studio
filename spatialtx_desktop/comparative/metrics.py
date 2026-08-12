from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__
from ..graph.builder import GraphBuildConfig, build_spatial_graph
from ..graph.context import (
    ContextFieldConfig,
    add_context_field,
    audit_context_program,
    resolve_context_program,
)
from ..graph.metadata import json_safe
from ..workflow import ScoringOptions, _read_h5ad, inspect_h5ad_memory, score_adata
from .comparative_normalization import (
    add_normalized_topology_metrics,
    noncentered_context_values,
    raw_context_summary,
    sample_scale_metrics,
)
from .metric_registry import DELTA_METRIC_SPECS, GROUP_METRICS
from .models import ComparativeConfig, NO_REGISTRATION_NOTICE, SampleRecord


COMPARATIVE_METRIC_LAYER_SCHEMA = "v0.6-hv-validation-v1"
CONTEXT_STATUS_VALUES = (
    "available",
    "insufficient_gene_coverage",
    "no_matched_genes",
    "unsupported_expression_scale",
    "graph_unavailable",
    "not_requested",
    "calculation_error",
)


def _operational_regime_confidence(
    regime: str,
    interface_fraction: float,
    diffuse_fraction: float,
    transition_burden_score: float,
) -> float:
    """Cautious descriptive confidence using established SpatialTX summary quantities."""
    if regime == "Type_A_candidate":
        value = 0.6 * interface_fraction + 0.4 * transition_burden_score
    else:
        value = 0.6 * diffuse_fraction + 0.4 * transition_burden_score
    return float(value)


def file_sha256(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def _cache_payload(record: SampleRecord, config: ComparativeConfig, input_hash: str) -> dict:
    h_program = resolve_context_program("H", config.h_genes)
    v_program = resolve_context_program("V", config.v_genes)
    return {
        "application_version": __version__,
        "comparative_metric_layer_schema": COMPARATIVE_METRIC_LAYER_SCHEMA,
        "input_hash": input_hash,
        "C_genes": config.c_genes,
        "S_genes": config.s_genes,
        "thresholds": {"C": config.c_q, "S": config.s_q, "G": config.g_q},
        "scoring_options": config.scoring_options,
        "graph_settings": config.graph_settings,
        "H": {
            "enabled": config.enable_h_expr,
            "genes": config.h_genes,
            "effective_genes": h_program["requested_genes"],
            "method": config.context_score_method,
            "smoothing": config.context_smoothing,
            "minimum_coverage": config.context_min_coverage,
        },
        "V": {
            "enabled": config.enable_v_expr,
            "genes": config.v_genes,
            "effective_genes": v_program["requested_genes"],
            "method": config.context_score_method,
            "smoothing": config.context_smoothing,
            "minimum_coverage": config.context_min_coverage,
        },
        "sample_id": record.sample_id,
    }


def cache_key(record: SampleRecord, config: ComparativeConfig, input_hash: str) -> str:
    payload = json.dumps(_cache_payload(record, config, input_hash), sort_keys=True, ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _cache_paths(cache_dir: Path, key: str) -> tuple[Path, Path]:
    return cache_dir / f"{key}.json", cache_dir / f"{key}.npz"


def _load_cache(cache_dir: Path, key: str) -> tuple[dict, dict[str, np.ndarray]] | None:
    metadata_path, fields_path = _cache_paths(cache_dir, key)
    if not metadata_path.is_file() or not fields_path.is_file():
        return None
    try:
        metrics = json.loads(metadata_path.read_text(encoding="utf-8"))
        with np.load(fields_path, allow_pickle=False) as data:
            fields = {name: np.asarray(data[name]) for name in data.files}
        return metrics, fields
    except (OSError, ValueError, json.JSONDecodeError):
        return None


def _write_cache(cache_dir: Path, key: str, metrics: dict, fields: dict[str, np.ndarray]) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    metadata_path, fields_path = _cache_paths(cache_dir, key)
    temporary_json = metadata_path.with_suffix(".tmp.json")
    temporary_npz = fields_path.with_suffix(".tmp.npz")
    temporary_json.write_text(json.dumps(json_safe(metrics), indent=2, ensure_ascii=False), encoding="utf-8")
    np.savez_compressed(temporary_npz, **fields)
    temporary_json.replace(metadata_path)
    temporary_npz.replace(fields_path)


def _active_context_values(adata, meta: dict) -> np.ndarray:
    columns = meta.get("obs_columns", {})
    active = columns.get("smoothed") or columns.get("base")
    return np.asarray(adata.obs[active], dtype=float)


def _context_summary(
    label: str,
    values: np.ndarray,
    meta: dict,
    transition_mask: np.ndarray,
) -> dict:
    transition = np.asarray(transition_mask, dtype=bool)
    other = ~transition
    transition_mean = float(np.nanmean(values[transition])) if transition.any() else np.nan
    other_mean = float(np.nanmean(values[other])) if other.any() else np.nan
    return {
        f"{label}_available": True,
        f"{label}_mean": float(np.nanmean(values)),
        f"{label}_median": float(np.nanmedian(values)),
        f"{label}_std": float(np.nanstd(values)),
        f"{label}_coverage_fraction": meta.get("coverage_fraction", np.nan),
        f"{label}_high_state_fraction": meta.get("high_state_fraction", np.nan),
        f"{label}_correlation_with_C": meta.get("correlation_with_C", np.nan),
        f"{label}_correlation_with_S": meta.get("correlation_with_S", np.nan),
        f"{label}_correlation_with_R": meta.get("correlation_with_R", np.nan),
        f"{label}_transition_mean": transition_mean,
        f"{label}_nontransition_mean": other_mean,
        f"{label}_transition_enrichment": transition_mean - other_mean if np.isfinite(transition_mean) and np.isfinite(other_mean) else np.nan,
        f"{label}_warnings": "; ".join(map(str, meta.get("warnings", []))),
        f"{label}_interpretation": meta.get("interpretation_limit", "Observational expression-derived context only."),
    }


def _missing_context(label: str, warning: str) -> dict:
    return {
        f"{label}_available": False,
        f"{label}_mean": np.nan,
        f"{label}_median": np.nan,
        f"{label}_std": np.nan,
        f"{label}_coverage_fraction": np.nan,
        f"{label}_high_state_fraction": np.nan,
        f"{label}_correlation_with_C": np.nan,
        f"{label}_correlation_with_S": np.nan,
        f"{label}_correlation_with_R": np.nan,
        f"{label}_transition_mean": np.nan,
        f"{label}_nontransition_mean": np.nan,
        f"{label}_transition_enrichment": np.nan,
        f"{label}_warnings": warning,
        f"{label}_interpretation": "Observational expression-derived context unavailable; core C/S results remain unchanged.",
    }


def _context_audit_metrics(
    axis: str,
    audit: dict,
    status: str,
    warning: str,
    raw_normalization_method: str,
) -> dict:
    if status not in CONTEXT_STATUS_VALUES:
        raise ValueError(f"Unsupported context status: {status}")
    payload = dict(audit)
    payload.update({
        "axis": axis,
        "context_available": status == "available",
        "context_status": status,
        "context_warning": warning,
        "raw_normalization_method": raw_normalization_method,
    })
    return {
        f"{axis}_context_available": status == "available",
        f"{axis}_context_status": status,
        f"{axis}_context_warning": warning,
        f"{axis}_context_audit": json_safe(payload),
    }


def analyze_sample(
    record: SampleRecord,
    config: ComparativeConfig,
    cache_dir: Path,
) -> tuple[dict, dict[str, np.ndarray], str]:
    input_hash = file_sha256(record.file_path)
    key = cache_key(record, config, input_hash)
    if config.use_cache:
        cached = _load_cache(cache_dir, key)
        if cached is not None:
            metrics, fields = cached
            metrics = dict(metrics)
            metrics.update(record.to_dict())
            metrics["cache_status"] = "reused"
            return metrics, fields, "reused"

    memory_info = inspect_h5ad_memory(record.file_path)
    adata = _read_h5ad(record.file_path)
    options = ScoringOptions(**config.scoring_options) if config.scoring_options else ScoringOptions()
    base, scored_fields = score_adata(
        adata,
        config.c_genes,
        config.s_genes,
        config.c_q,
        config.s_q,
        config.g_q,
        options=options,
        preflight_info=memory_info,
        source_path=record.file_path,
        sample_name=record.sample_id,
        gene_program_mode="comparative",
    )
    if not scored_fields.get("spatial_available", False):
        raise ValueError("Valid spatial coordinates are required for Comparative Analysis.")
    C = np.asarray(scored_fields["C"], dtype=float)
    S = np.asarray(scored_fields["S"], dtype=float)
    R = np.asarray(scored_fields["R"], dtype=float)
    G = np.asarray(scored_fields["G"], dtype=float)
    transition_mask = np.asarray(scored_fields["interface"], dtype=bool) | np.asarray(scored_fields["diffuse"], dtype=bool)
    regime_confidence = _operational_regime_confidence(
        str(base["regime_label"]),
        float(base["interface_fraction"]),
        float(base["diffuse_fraction"]),
        float(base["transition_burden_score"]),
    )
    metrics = {
        **base,
        **record.to_dict(),
        "input_sha256": input_hash,
        "analysis_status": "ok",
        "C_median": float(np.nanmedian(C)),
        "S_median": float(np.nanmedian(S)),
        "R_median": float(np.nanmedian(R)),
        "R_std": float(np.nanstd(R)),
        "gradient_mean": float(np.nanmean(G)),
        "gradient_q90": float(np.nanquantile(G, 0.90)),
        "localized_interface_fraction": float(base["interface_fraction"]),
        "localized_interface_spots": int(base["interface_spots"]),
        "regime_confidence": float(regime_confidence),
        "registration_status": "not_performed",
        "registration_notice": NO_REGISTRATION_NOTICE,
        "cache_status": "computed",
    }
    field_cache: dict[str, np.ndarray] = {
        "coords": np.asarray(scored_fields["coords"], dtype=float),
        "C": C,
        "S": S,
        "R": R,
        "G": G,
        "interface": np.asarray(scored_fields["interface"], dtype=np.uint8),
        "diffuse": np.asarray(scored_fields["diffuse"], dtype=np.uint8),
    }
    context_warnings: list[str] = []
    context_requested = config.enable_h_expr or config.enable_v_expr
    graph_kwargs = {
        key: value for key, value in config.graph_settings.items()
        if key in GraphBuildConfig.__dataclass_fields__
    }
    graph_result = build_spatial_graph(adata, GraphBuildConfig(**graph_kwargs))
    metrics.update(sample_scale_metrics(adata, scored_fields, transition_mask, graph_result))
    metrics = add_normalized_topology_metrics(metrics)
    if int(graph_result.qc.get("n_edges", 0)) > 0:
        from scipy import sparse

        upper = sparse.triu(graph_result.connectivities, k=1).tocoo()
        field_cache["context_edge_i"] = np.asarray(upper.row, dtype=np.int64)
        field_cache["context_edge_j"] = np.asarray(upper.col, dtype=np.int64)
    elif context_requested:
        graph_result = None
        context_warnings.append("H/V context skipped because graph construction produced no usable edges.")
    if context_requested and graph_result is None:
        if not context_warnings:
            context_warnings.append("H/V context skipped because graph construction produced no usable edges.")
    for field, enabled, genes in (
        ("H", config.enable_h_expr, config.h_genes),
        ("V", config.enable_v_expr, config.v_genes),
    ):
        label = f"{field}_expr"
        field_config = ContextFieldConfig(
            field=field,
            genes=genes,
            score_method=config.context_score_method,
            min_coverage=config.context_min_coverage,
            allow_low_coverage=False,
            smoothing=config.context_smoothing,
        )
        try:
            audit = audit_context_program(adata, field_config)
        except Exception as exc:
            program = resolve_context_program(field, genes)
            audit = {
                **program,
                "matched_genes": [],
                "missing_genes": list(program["requested_genes"]),
                "matched_gene_count": 0,
                "coverage_fraction": 0.0,
                "genes_expressed_above_min_spot_fraction": [],
                "expressed_gene_count": 0,
                "expressed_gene_fraction": 0.0,
                "expression_scale_guess": "unknown",
                "detection_source": "audit_calculation_error",
                "raw_normalization_method": "unavailable_audit_calculation_error",
                "minimum_coverage": float(config.context_min_coverage),
            }
            warning = f"{label} gene audit failed: {exc}"
            metrics.update(_missing_context(label, warning))
            metrics.update(raw_context_summary(field, None, transition_mask, "unavailable_audit_calculation_error"))
            metrics.update(_context_audit_metrics(
                field,
                audit,
                "calculation_error",
                warning,
                "unavailable_audit_calculation_error",
            ))
            context_warnings.append(warning)
            continue
        if not enabled:
            metrics.update(_missing_context(label, "Not requested."))
            metrics.update(raw_context_summary(field, None, transition_mask, "not_requested"))
            metrics.update(_context_audit_metrics(
                field,
                audit,
                "not_requested",
                "Not requested.",
                "not_requested",
            ))
            continue
        if int(audit.get("matched_gene_count", 0)) == 0:
            warning = f"No requested genes were found for {label}."
            metrics.update(_missing_context(label, warning))
            metrics.update(raw_context_summary(field, None, transition_mask, "unavailable_no_matched_genes"))
            metrics.update(_context_audit_metrics(
                field,
                audit,
                "no_matched_genes",
                warning,
                "unavailable_no_matched_genes",
            ))
            context_warnings.append(warning)
            continue
        if float(audit.get("coverage_fraction", 0.0)) < float(config.context_min_coverage):
            warning = (
                f"low {label} gene coverage: "
                f"{float(audit.get('coverage_fraction', 0.0)):.1%}"
            )
            metrics.update(_missing_context(label, warning))
            metrics.update(raw_context_summary(field, None, transition_mask, "unavailable_insufficient_gene_coverage"))
            metrics.update(_context_audit_metrics(
                field,
                audit,
                "insufficient_gene_coverage",
                warning,
                "unavailable_insufficient_gene_coverage",
            ))
            context_warnings.append(warning)
            continue
        if graph_result is None:
            warning = f"{label} unavailable because a usable context graph was not available."
            metrics.update(_missing_context(label, warning))
            metrics.update(raw_context_summary(field, None, transition_mask, "unavailable_graph"))
            metrics.update(_context_audit_metrics(
                field,
                audit,
                "graph_unavailable",
                warning,
                "unavailable_graph",
            ))
            context_warnings.append(warning)
            continue
        try:
            _coverage, meta = add_context_field(
                adata,
                field_config,
                graph_result.connectivities,
                reference_fields={"C": C, "S": S, "R": R},
                active_graph=graph_result.method,
            )
            for key in (
                "gene_set_name", "requested_genes", "matched_genes", "missing_genes",
                "requested_gene_count", "matched_gene_count", "coverage_fraction",
                "genes_expressed_above_min_spot_fraction", "expressed_gene_fraction",
                "expression_scale_guess", "detection_source",
            ):
                if key in meta:
                    audit[key] = meta[key]
            audit["expressed_gene_count"] = len(
                audit.get("genes_expressed_above_min_spot_fraction", [])
            )
            values = _active_context_values(adata, meta)
            field_cache[label] = values
            metrics.update(_context_summary(label, values, meta, transition_mask))
            raw_values, raw_method = noncentered_context_values(
                adata,
                meta.get("matched_genes", []),
                str(meta.get("expression_scale_guess", "")),
            )
            if raw_values is not None:
                field_cache[f"{field}_expr_raw"] = np.asarray(raw_values, dtype=float)
            metrics.update(raw_context_summary(field, raw_values, transition_mask, raw_method))
            if raw_values is None and str(raw_method).startswith("unavailable_expression_scale_"):
                context_status = "unsupported_expression_scale"
                warning = (
                    f"{label} non-centered comparison unavailable for expression scale "
                    f"{audit.get('expression_scale_guess', 'unknown')}."
                )
            elif raw_values is None:
                context_status = "calculation_error"
                warning = f"{label} non-centered context summary could not be calculated."
            else:
                context_status = "available"
                warning = "; ".join(map(str, meta.get("warnings", [])))
            metrics.update(_context_audit_metrics(
                field,
                audit,
                context_status,
                warning,
                raw_method,
            ))
            context_warnings.extend(map(str, meta.get("warnings", [])))
            if context_status != "available":
                context_warnings.append(warning)
        except Exception as exc:
            warning = f"{label} skipped: {exc}"
            metrics.update(_missing_context(label, warning))
            metrics.update(raw_context_summary(field, None, transition_mask, "unavailable_context_error"))
            metrics.update(_context_audit_metrics(
                field,
                audit,
                "calculation_error",
                warning,
                "unavailable_context_error",
            ))
            context_warnings.append(warning)
    combined_warnings = [str(base.get("QC_notes", "")).strip(), *context_warnings]
    metrics["warning_messages"] = "; ".join(item for item in combined_warnings if item)
    metrics["H_V_observational_only"] = True
    metrics["H_V_core_effect"] = "none; H/V do not alter C/S/R, transition masks, or Type A/B/C"
    metrics["H_expr_mean_interpretation_flag"] = "non_informative_centered_mean"
    metrics["V_expr_mean_interpretation_flag"] = "non_informative_centered_mean"
    if config.use_cache:
        _write_cache(cache_dir, key, metrics, field_cache)
    return metrics, field_cache, "computed"


def compute_delta_metrics(sample_metrics: pd.DataFrame, matches: pd.DataFrame) -> pd.DataFrame:
    by_sample = sample_metrics.set_index("sample_id")
    rows: list[dict] = []
    for _, match in matches.iterrows():
        reference_id = str(match["reference_sample_id"])
        target_id = str(match["target_sample_id"])
        reference_row = by_sample.loc[reference_id]
        target_row = by_sample.loc[target_id]
        for delta_name, metric in DELTA_METRIC_SPECS:
            reference_value = pd.to_numeric(pd.Series([reference_row.get(metric)]), errors="coerce").iloc[0]
            target_value = pd.to_numeric(pd.Series([target_row.get(metric)]), errors="coerce").iloc[0]
            available = bool(np.isfinite(reference_value) and np.isfinite(target_value))
            delta = float(target_value - reference_value) if available else np.nan
            rows.append({
                "comparison_id": match["comparison_id"],
                "pair_id": match.get("pair_id", ""),
                "reference_sample_id": reference_id,
                "target_sample_id": target_id,
                "metric": metric,
                "delta_metric": delta_name,
                "reference_value": reference_value,
                "target_value": target_value,
                "delta": delta,
                "direction_definition": "Target - Reference",
                "positive_delta_definition": "Target > Reference",
                "standardized_delta": np.nan,
                "standardized_delta_status": "not_computed_for_pairwise",
                "status": "ok" if available else "unavailable_optional_metric",
            })
    return pd.DataFrame(rows)


def group_mean_deltas(sample_metrics: pd.DataFrame, reference: str, target: str) -> pd.DataFrame:
    rows: list[dict] = []
    for delta_name, metric in DELTA_METRIC_SPECS:
        if metric not in sample_metrics:
            ref = pd.Series(dtype=float)
            tar = pd.Series(dtype=float)
        else:
            ref = pd.to_numeric(sample_metrics.loc[sample_metrics["group"].eq(reference), metric], errors="coerce").dropna()
            tar = pd.to_numeric(sample_metrics.loc[sample_metrics["group"].eq(target), metric], errors="coerce").dropna()
        available = len(ref) > 0 and len(tar) > 0
        reference_value = float(ref.mean()) if len(ref) else np.nan
        target_value = float(tar.mean()) if len(tar) else np.nan
        pooled = pd.concat([ref, tar], ignore_index=True)
        pooled_scale = float(pooled.std(ddof=1)) if len(pooled) > 1 else np.nan
        standardized_delta = (
            float((target_value - reference_value) / pooled_scale)
            if available and np.isfinite(pooled_scale) and pooled_scale > 1e-12
            else np.nan
        )
        rows.append({
            "comparison_id": f"{reference}_vs_{target}_group_mean",
            "pair_id": "",
            "reference_sample_id": f"group:{reference}",
            "target_sample_id": f"group:{target}",
            "metric": metric,
            "delta_metric": delta_name,
            "reference_value": reference_value,
            "target_value": target_value,
            "delta": target_value - reference_value if available else np.nan,
            "direction_definition": "Target group mean - Reference group mean",
            "positive_delta_definition": "Target group mean > Reference group mean",
            "standardized_delta": standardized_delta,
            "standardized_delta_status": "pooled_sample_scale" if np.isfinite(standardized_delta) else "unavailable_pooled_scale",
            "status": "ok" if available else "unavailable_optional_metric",
        })
    return pd.DataFrame(rows)


def regime_transition_table(
    sample_metrics: pd.DataFrame,
    matches: pd.DataFrame | None,
    reference: str,
    target: str,
    confidence_threshold: float,
) -> pd.DataFrame:
    rows: list[dict] = []
    if matches is not None and not matches.empty:
        by_sample = sample_metrics.set_index("sample_id")
        for _, match in matches.iterrows():
            ref = by_sample.loc[str(match["reference_sample_id"])]
            tar = by_sample.loc[str(match["target_sample_id"])]
            ref_conf = float(ref.get("regime_confidence", np.nan))
            tar_conf = float(tar.get("regime_confidence", np.nan))
            uncertain = (
                not np.isfinite(ref_conf)
                or not np.isfinite(tar_conf)
                or min(ref_conf, tar_conf) < confidence_threshold
                or str(ref.get("spatial_qc_status", "")) != "PASS"
                or str(tar.get("spatial_qc_status", "")) != "PASS"
            )
            ref_regime = str(ref.get("regime_label", ""))
            tar_regime = str(tar.get("regime_label", ""))
            rows.append({
                "comparison_id": match["comparison_id"],
                "pair_id": match.get("pair_id", ""),
                "comparison_basis": "matched_sample_operational_change",
                "reference_regime": ref_regime,
                "target_regime": tar_regime,
                "regime_transition": f"{ref_regime} → {tar_regime}",
                "reference_confidence": ref_conf,
                "target_confidence": tar_conf,
                "transition_confidence_flag": "uncertain" if uncertain else "adequate_for_descriptive_comparison",
                "interpretation": "Descriptive operational change; not a validated biological state transition.",
            })
    else:
        for regime in ("Type_A_candidate", "Type_B_candidate", "Type_C_candidate"):
            ref_values = sample_metrics.loc[sample_metrics["group"].eq(reference), "regime_label"].astype(str)
            tar_values = sample_metrics.loc[sample_metrics["group"].eq(target), "regime_label"].astype(str)
            rows.append({
                "comparison_id": f"{reference}_vs_{target}_distribution",
                "pair_id": "",
                "comparison_basis": "unpaired_group_distribution",
                "reference_regime": regime,
                "target_regime": regime,
                "regime_transition": "No direct transition inferred; group distribution comparison",
                "reference_count": int(ref_values.eq(regime).sum()),
                "target_count": int(tar_values.eq(regime).sum()),
                "reference_fraction": float(ref_values.eq(regime).mean()) if len(ref_values) else np.nan,
                "target_fraction": float(tar_values.eq(regime).mean()) if len(tar_values) else np.nan,
                "reference_confidence": np.nan,
                "target_confidence": np.nan,
                "transition_confidence_flag": "not_applicable_unpaired_distribution",
                "interpretation": "Descriptive group distribution; not a validated biological state transition.",
            })
    return pd.DataFrame(rows)
