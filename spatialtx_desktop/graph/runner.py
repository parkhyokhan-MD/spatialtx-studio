from __future__ import annotations

import datetime as dt
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from threading import Event
from typing import Callable

import numpy as np
import pandas as pd

from .. import __version__
from ..gene_program_validation import validate_gene_programs
from ..workflow import _read_h5ad, score_adata
from .audit import audit_input
from .builder import GraphBuildConfig, build_spatial_graph, store_graph
from .context import ContextFieldConfig, add_context_field
from .continuous import continuous_edge_statistics
from .enrichment import binary_mask_association, categorical_neighborhood_enrichment
from .metadata import (
    FDR_SCOPE,
    PERMUTATION_LIMITATION,
    SMOOTHING_LIMITATION,
    SPATIAL_ASSOCIATION_LIMITATION,
    json_safe,
)
from .plotting import plot_context_map, plot_enrichment_heatmap, plot_graph_qc, plot_joint_hv_map
from .robustness import compare_graph_robustness
from .semantics import describe_obs_variable, semantics_table


Progress = Callable[[str], None]
MODULE_VERSION = "0.4"


@dataclass
class SpatialGraphAnalysisConfig:
    graph: GraphBuildConfig
    c_genes: list[str] | None = None
    s_genes: list[str] | None = None
    enable_h: bool = True
    enable_v: bool = True
    h_genes: list[str] | None = None
    v_genes: list[str] | None = None
    h_score_method: str = "z_score_mean"
    v_score_method: str = "z_score_mean"
    h_high_quantile: float = 0.80
    v_high_quantile: float = 0.80
    context_smoothing: str = "none"
    min_gene_coverage: float = 0.25
    allow_low_coverage: bool = False
    context_min_spot_fraction: float = 0.01
    label_source: str = "auto_spatialtx_states"
    label_mode: str = "categorical_state"
    user_mask_a_column: str = ""
    user_mask_b_column: str = ""
    continuous_x_column: str = ""
    continuous_y_column: str = ""
    continuous_x_mode: str = "continuous_score"
    continuous_y_mode: str = "continuous_score"
    permutation_scope: str = "whole_slide"
    stratification_column: str = ""
    tissue_only_restriction: bool = False
    run_graph_robustness: bool = False
    permutations: int = 999
    seed: int = 20260713
    c_q: float = 0.80
    s_q: float = 0.80
    g_q: float = 0.60
    write_annotated_h5ad: bool = False
    max_plot_edges: int = 50000
    plot_edge_seed: int = 42


def _ensure_dirs(run_dir: Path) -> dict[str, Path]:
    paths = {
        "spatial_graph": run_dir / "spatial_graph",
        "input_audit": run_dir / "input_audit",
        "variable_semantics": run_dir / "variable_semantics",
        "neighborhood": run_dir / "neighborhood",
        "context_fields": run_dir / "context_fields",
        "robustness": run_dir / "robustness",
        "figures": run_dir / "figures",
        "annotated": run_dir / "annotated",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def _spatialtx_labels(fields: dict, c_q: float, s_q: float) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    C = np.asarray(fields["C"], dtype=float)
    S = np.asarray(fields["S"], dtype=float)
    interface = np.asarray(fields.get("interface", np.zeros(len(C), dtype=bool)), dtype=bool)
    diffuse = np.asarray(fields.get("diffuse", np.zeros(len(C), dtype=bool)), dtype=bool)
    c_threshold = float(np.nanquantile(C[np.isfinite(C)], c_q)) if np.isfinite(C).any() else np.nan
    s_threshold = float(np.nanquantile(S[np.isfinite(S)], s_q)) if np.isfinite(S).any() else np.nan
    c_high = C >= c_threshold if np.isfinite(c_threshold) else np.zeros(len(C), dtype=bool)
    s_high = S >= s_threshold if np.isfinite(s_threshold) else np.zeros(len(C), dtype=bool)
    labels = np.full(len(C), "balanced_or_transition_like", dtype=object)
    labels[c_high & ~s_high] = "C_dominant"
    labels[s_high & ~c_high] = "S_dominant"
    labels[diffuse] = "diffuse_transition"
    labels[interface] = "localized_interface_like"
    masks = {
        "C_dominant": c_high & ~s_high,
        "S_dominant": s_high & ~c_high,
        "transition_like": interface | diffuse | (c_high & s_high),
        "localized_interface_like": interface,
        "diffuse_transition": diffuse,
    }
    return labels, masks


def _obs_labels(adata, column: str) -> np.ndarray:
    if column not in adata.obs:
        raise ValueError(f"Label source column not found in adata.obs: {column}")
    values = np.asarray(adata.obs[column], dtype=object)
    return values


def _obs_binary_mask(adata, column: str) -> np.ndarray:
    if column not in adata.obs:
        raise ValueError(f"Binary-mask source column not found in adata.obs: {column}")
    values = adata.obs[column]
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False).to_numpy(dtype=bool)
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().any():
        return numeric.fillna(0).to_numpy(dtype=float) > 0
    true_values = {"true", "yes", "y", "1", "positive", "high"}
    return values.fillna("").astype(str).str.strip().str.lower().isin(true_values).to_numpy(dtype=bool)


def _safe_statistic_name(x_source: str, y_source: str) -> str:
    x = re.sub(r"[^A-Za-z0-9]+", "_", x_source).strip("_") or "X"
    y = re.sub(r"[^A-Za-z0-9]+", "_", y_source).strip("_") or "Y"
    return f"I_{x}_{y}"


def _context_masks(adata) -> dict[str, np.ndarray]:
    masks: dict[str, np.ndarray] = {}
    if "spatialtx_H_high" in adata.obs:
        h = np.asarray(adata.obs["spatialtx_H_high"], dtype=bool)
        masks["H_high"] = h
        masks["H_low"] = ~h
    if "spatialtx_V_high" in adata.obs:
        v = np.asarray(adata.obs["spatialtx_V_high"], dtype=bool)
        masks["V_high"] = v
        masks["V_low"] = ~v
    return masks


def _mask_label_source(name: str) -> str:
    if str(name).startswith("obs:"):
        return "user_supplied"
    if str(name).startswith(("H_", "V_")):
        return "expression_derived_context"
    return "spatialtx_derived"


def _interpretation_for_sources(*sources: str) -> str:
    if "spatialtx_derived" in sources:
        return "descriptive_spatial_organization"
    return "exploratory_neighbor_association"


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(json_safe(payload), indent=2, ensure_ascii=False), encoding="utf-8")


def run_spatial_graph_neighborhood_sample(
    source: Path,
    out: dict[str, Path],
    config: SpatialGraphAnalysisConfig,
    progress: Progress | None = None,
) -> dict:
    sample = source.stem
    if progress:
        progress(f"Building graph/context workflow for {sample}")
    adata = _read_h5ad(source)
    audit, audit_table = audit_input(
        adata,
        input_filename=source,
        coordinate_unit=config.graph.coordinate_unit,
        coordinate_scale=config.graph.coordinate_scale,
        scale_source=config.graph.scale_source,
    )
    _write_json(out["input_audit"] / f"input_audit_{sample}.json", audit)
    audit_table.to_csv(out["input_audit"] / f"input_audit_{sample}.csv", index=False)
    c_genes = config.c_genes or []
    s_genes = config.s_genes or []
    if not c_genes or not s_genes:
        raise ValueError("C/S gene programs are required for v0.4 graph labels and C/S edge statistics.")
    metrics, fields = score_adata(
        adata,
        c_genes,
        s_genes,
        config.c_q,
        config.s_q,
        config.g_q,
        source_path=source,
        sample_name=sample,
        gene_program_mode="graph",
    )
    if not fields.get("spatial_available", False):
        raise ValueError("Valid spatial coordinates are required for v0.4 spatial graph analysis.")
    result = build_spatial_graph(adata, config.graph)
    if not result.qc.get("coordinate_valid", False):
        raise ValueError("Valid spatial coordinates are required for graph analysis.")
    if int(result.qc.get("n_edges", 0)) < 1:
        raise ValueError(
            "The selected graph parameters produced no usable spatial edges. "
            "Choose an automatic/larger radius or review the coordinate scale."
        )
    store_graph(adata, result)
    strata = None
    if config.permutation_scope == "within_user_strata":
        if not config.stratification_column:
            raise ValueError("within_user_strata requires a stratification column")
        if config.stratification_column not in adata.obs:
            raise ValueError(f"Stratification column not found in adata.obs: {config.stratification_column}")
        strata = np.asarray(adata.obs[config.stratification_column], dtype=object)
    tissue_mask = None
    if config.tissue_only_restriction:
        if "in_tissue" in adata.obs:
            tissue_mask = pd.to_numeric(adata.obs["in_tissue"], errors="coerce").fillna(0).to_numpy(dtype=float) == 1
        elif progress:
            progress("  Tissue-only restriction requested but in_tissue is unavailable; continuing with audit warning.")
    graph_prefix = f"{sample}_"
    _write_json(out["spatial_graph"] / f"{graph_prefix}graph_metadata.json", {
        "application": "SpatialTX Studio Desktop",
        "application_version": __version__,
        "analysis_module": "Spatial Graph and Neighborhood Analysis",
        "analysis_module_version": MODULE_VERSION,
        "created": dt.datetime.now().isoformat(timespec="seconds"),
        "input_filename": str(source.resolve()),
        "output_path": str(out["spatial_graph"].resolve()),
        "graph": result.metadata,
        "warnings": result.warnings,
        "analysis_status": "experimental",
        "spatial_association_limitation": SPATIAL_ASSOCIATION_LIMITATION,
        "permutation_limitation": PERMUTATION_LIMITATION,
        "fdr_scope": FDR_SCOPE,
        "gene_program_validation": fields.get("gene_program_validation", {}),
    })
    pd.DataFrame([{**result.qc, "sample": sample, "method": result.method}]).to_csv(out["spatial_graph"] / f"{graph_prefix}graph_qc.csv", index=False)
    result.degree_table.to_csv(out["spatial_graph"] / f"{graph_prefix}graph_degree_distribution.csv", index=False)
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    degree = np.asarray((result.connectivities > 0).sum(axis=1)).ravel()
    plot_graph_qc(
        coords,
        result.connectivities,
        degree == 0,
        out["figures"] / f"{graph_prefix}graph_qc.png",
        f"{sample} spatial graph QC",
        max_plot_edges=config.max_plot_edges,
        plot_edge_seed=config.plot_edge_seed,
    )
    context_rows: list[pd.DataFrame] = []
    context_metadata: dict[str, dict] = {}
    context_status_rows: list[dict] = []

    def _optional_context(field: str, genes, method: str, quantile: float, cmap: str, title: str) -> None:
        label = "H_expr" if field == "H" else "V_expr"
        try:
            table, meta = add_context_field(
                adata,
                ContextFieldConfig(
                    field=field,
                    genes=genes,
                    score_method=method,
                    min_coverage=config.min_gene_coverage,
                    allow_low_coverage=config.allow_low_coverage,
                    smoothing=config.context_smoothing,
                    high_quantile=quantile,
                    min_spot_fraction=config.context_min_spot_fraction,
                ),
                result.connectivities,
                reference_fields={"C": fields["C"], "S": fields["S"], "R": fields["R"]},
                active_graph=result.method,
            )
            meta = {"status": "ok", **meta}
            context_rows.append(table)
            context_metadata[label] = meta
            context_status_rows.append(meta)
            base_col = meta["obs_columns"]["base"]
            smoothed_col = meta["obs_columns"]["smoothed"]
            value_col = smoothed_col or base_col
            plot_context_map(coords, np.asarray(adata.obs[value_col], dtype=float), out["figures"] / f"{graph_prefix}{label}_map.png", f"{sample} {title}", cmap)
            if smoothed_col:
                plot_context_map(
                    coords,
                    np.asarray(adata.obs[base_col], dtype=float),
                    out["figures"] / f"{graph_prefix}{label}_unsmoothed_map.png",
                    f"{sample} {title} (unsmoothed)",
                    cmap,
                )
                plot_context_map(
                    coords,
                    np.asarray(adata.obs[smoothed_col], dtype=float),
                    out["figures"] / f"{graph_prefix}{label}_smoothed_map.png",
                    f"{sample} {title} (graph-smoothed; exploratory sensitivity)",
                    cmap,
                )
        except ValueError as exc:
            skipped = {
                "status": "skipped_qc",
                "field": label,
                "public_term": "hypoxia-associated expression field" if field == "H" else "endothelial/angiogenic expression proxy",
                "warnings": [str(exc)],
                "coverage_warning": str(exc),
                "version": __version__,
            }
            context_status_rows.append(skipped)
            _write_json(out["context_fields"] / f"{sample}_{label}_qc.json", skipped)
            if progress:
                progress(f"  {label} skipped by context QC: {exc}")

    if config.enable_h:
        _optional_context("H", config.h_genes, config.h_score_method, config.h_high_quantile, "Purples", "H_expr hypoxia-associated expression field")
    if config.enable_v:
        _optional_context("V", config.v_genes, config.v_score_method, config.v_high_quantile, "Greens", "V_expr endothelial/angiogenic expression proxy")
    if context_rows:
        context_table = pd.concat(context_rows, ignore_index=True)
        context_table.to_csv(out["context_fields"] / f"{graph_prefix}context_field_gene_coverage.csv", index=False)
        loo_columns = [
            "field", "gene", "field_correlation_without_gene",
            "mean_absolute_field_change", "variance_change", "rank",
        ]
        available_loo = [column for column in loo_columns if column in context_table]
        if available_loo:
            context_table[available_loo].dropna(subset=["rank"], how="all").to_csv(
                out["context_fields"] / f"{graph_prefix}context_field_leave_one_gene_out.csv",
                index=False,
            )
    if context_status_rows:
        summary_rows = []
        for meta in context_status_rows:
            key = meta.get("field", "")
            row = {"sample": sample, "field": key, **meta}
            for list_key in ("requested_genes", "matched_genes", "missing_genes", "genes_expressed_above_min_spot_fraction", "warnings"):
                if isinstance(row.get(list_key), list):
                    row[list_key] = ";".join(map(str, row[list_key]))
            if isinstance(row.get("per_gene_detection_fraction"), dict):
                row["per_gene_detection_fraction"] = json.dumps(row["per_gene_detection_fraction"], ensure_ascii=False)
            if isinstance(row.get("per_gene_positive_value_fraction"), dict):
                row["per_gene_positive_value_fraction"] = json.dumps(row["per_gene_positive_value_fraction"], ensure_ascii=False)
            if isinstance(row.get("leave_one_gene_out_qc"), list):
                row["leave_one_gene_out_qc"] = json.dumps(row["leave_one_gene_out_qc"], ensure_ascii=False)
            if isinstance(row.get("obs_columns"), dict):
                row["obs_columns"] = json.dumps(row["obs_columns"], ensure_ascii=False)
            summary_rows.append(row)
            if meta.get("status") == "ok":
                _write_json(out["context_fields"] / f"{sample}_{key}_qc.json", meta)
        pd.DataFrame(summary_rows).to_csv(out["context_fields"] / f"{graph_prefix}context_field_summary.csv", index=False)
    if "spatialtx_H_high" in adata.obs and "spatialtx_V_high" in adata.obs:
        plot_joint_hv_map(
            coords,
            np.asarray(adata.obs["spatialtx_H_high"], dtype=bool),
            np.asarray(adata.obs["spatialtx_V_high"], dtype=bool),
            out["figures"] / f"{graph_prefix}H_V_association_map.png",
            f"{sample} H/V expression-context joint state",
        )
    labels, masks = _spatialtx_labels(fields, config.c_q, config.s_q)
    variable_records: list[dict] = []
    if config.label_source != "auto_spatialtx_states":
        variable_records.append(describe_obs_variable(adata, config.label_source, confirmed_mode=config.label_mode))
        labels = _obs_labels(adata, config.label_source)
    if config.stratification_column:
        variable_records.append(describe_obs_variable(adata, config.stratification_column, confirmed_mode="categorical_state"))
    categorical = categorical_neighborhood_enrichment(
        result.connectivities,
        labels,
        permutations=config.permutations,
        seed=config.seed,
        permutation_scope=config.permutation_scope,
        strata=strata,
        tissue_mask=tissue_mask,
    )
    categorical_source = "spatialtx_derived" if config.label_source == "auto_spatialtx_states" else "user_supplied"
    categorical["label_source"] = categorical_source
    categorical["analysis_interpretation"] = _interpretation_for_sources(categorical_source)
    categorical["stratification_column"] = config.stratification_column if config.permutation_scope == "within_user_strata" else ""
    categorical.insert(0, "sample", sample)
    categorical.to_csv(out["neighborhood"] / f"{graph_prefix}categorical_enrichment.csv", index=False)
    plot_enrichment_heatmap(categorical, out["figures"] / f"{graph_prefix}neighborhood_enrichment_heatmap.png", f"{sample} neighborhood enrichment")
    masks.update(_context_masks(adata))
    if config.user_mask_a_column:
        variable_records.append(describe_obs_variable(adata, config.user_mask_a_column, confirmed_mode="binary_mask"))
        masks[f"obs:{config.user_mask_a_column}"] = _obs_binary_mask(adata, config.user_mask_a_column)
    if config.user_mask_b_column:
        variable_records.append(describe_obs_variable(adata, config.user_mask_b_column, confirmed_mode="binary_mask"))
        masks[f"obs:{config.user_mask_b_column}"] = _obs_binary_mask(adata, config.user_mask_b_column)
    binary_pairs = [
        ("C_dominant", "H_high"), ("S_dominant", "H_high"), ("transition_like", "H_high"),
        ("localized_interface_like", "H_high"), ("diffuse_transition", "H_high"),
        ("C_dominant", "V_high"), ("S_dominant", "V_high"), ("transition_like", "V_high"),
        ("H_high", "V_high"), ("H_high", "V_low"),
    ]
    if config.user_mask_a_column and config.user_mask_b_column:
        binary_pairs.append((f"obs:{config.user_mask_a_column}", f"obs:{config.user_mask_b_column}"))
    binary_pairs = [pair for pair in binary_pairs if pair[0] in masks and pair[1] in masks]
    binary = binary_mask_association(
        result.connectivities,
        masks,
        pairs=binary_pairs,
        permutations=config.permutations,
        seed=config.seed,
        permutation_scope=config.permutation_scope,
        strata=strata,
        tissue_mask=tissue_mask,
    )
    if len(binary):
        binary["label_source_a"] = binary["mask_a"].map(_mask_label_source)
        binary["label_source_b"] = binary["mask_b"].map(_mask_label_source)
        binary["label_source"] = binary["label_source_a"] + "+" + binary["label_source_b"]
        binary["analysis_interpretation"] = [
            _interpretation_for_sources(a, b)
            for a, b in zip(binary["label_source_a"], binary["label_source_b"])
        ]
    binary["stratification_column"] = config.stratification_column if config.permutation_scope == "within_user_strata" else ""
    binary.insert(0, "sample", sample)
    binary.to_csv(out["neighborhood"] / f"{graph_prefix}binary_mask_association.csv", index=False)
    same_spot = binary[binary["mode"].eq("same_spot_overlap")].copy() if "mode" in binary else binary.copy()
    neighboring = binary[binary["mode"].eq("neighboring_spot_association")].copy() if "mode" in binary else binary.copy()
    same_spot.to_csv(out["neighborhood"] / f"{graph_prefix}same_spot_overlap.csv", index=False)
    neighboring.to_csv(out["neighborhood"] / f"{graph_prefix}neighboring_spot_association.csv", index=False)
    continuous_fields = {
        "C": np.asarray(fields["C"], dtype=float),
        "S": np.asarray(fields["S"], dtype=float),
        "R": np.asarray(fields["R"], dtype=float),
    }
    custom_pairs: list[tuple[str, str, str]] = []
    if "H_expr" in context_metadata:
        h_columns = context_metadata["H_expr"]["obs_columns"]
        h_active = h_columns["smoothed"] or h_columns["base"]
        continuous_fields["H"] = np.asarray(adata.obs[h_active], dtype=float)
        continuous_fields["H_unsmoothed"] = np.asarray(adata.obs[h_columns["base"]], dtype=float)
        if h_columns["smoothed"]:
            continuous_fields["H_smoothed"] = np.asarray(adata.obs[h_columns["smoothed"]], dtype=float)
            custom_pairs.extend([
                ("I_RH_unsmoothed", "R", "H_unsmoothed"),
                ("I_RH_smoothed", "R", "H_smoothed"),
            ])
    if "V_expr" in context_metadata:
        v_columns = context_metadata["V_expr"]["obs_columns"]
        v_active = v_columns["smoothed"] or v_columns["base"]
        continuous_fields["V"] = np.asarray(adata.obs[v_active], dtype=float)
        continuous_fields["V_unsmoothed"] = np.asarray(adata.obs[v_columns["base"]], dtype=float)
        if v_columns["smoothed"]:
            continuous_fields["V_smoothed"] = np.asarray(adata.obs[v_columns["smoothed"]], dtype=float)
            custom_pairs.extend([
                ("I_RV_unsmoothed", "R", "V_unsmoothed"),
                ("I_RV_smoothed", "R", "V_smoothed"),
            ])
    if bool(config.continuous_x_column) != bool(config.continuous_y_column):
        raise ValueError("Generic continuous edge association requires both X and Y source columns.")
    if config.continuous_x_column and config.continuous_y_column:
        variable_records.append(describe_obs_variable(adata, config.continuous_x_column, confirmed_mode=config.continuous_x_mode))
        variable_records.append(describe_obs_variable(adata, config.continuous_y_column, confirmed_mode=config.continuous_y_mode))
        continuous_fields[config.continuous_x_column] = pd.to_numeric(adata.obs[config.continuous_x_column], errors="coerce").to_numpy(dtype=float)
        continuous_fields[config.continuous_y_column] = pd.to_numeric(adata.obs[config.continuous_y_column], errors="coerce").to_numpy(dtype=float)
        custom_pairs.append((_safe_statistic_name(config.continuous_x_column, config.continuous_y_column), config.continuous_x_column, config.continuous_y_column))
    continuous = continuous_edge_statistics(
        result.connectivities,
        continuous_fields,
        permutations=config.permutations,
        seed=config.seed,
        custom_pairs=custom_pairs,
        permutation_scope=config.permutation_scope,
        strata=strata,
        tissue_mask=tissue_mask,
        active_graph=result.method,
        graph_parameters=result.metadata,
        stratification_column=config.stratification_column,
    )
    continuous["label_source"] = "expression_derived_fields"
    continuous["analysis_interpretation"] = "exploratory_neighbor_association"
    if len(continuous):
        continuous["uses_graph_smoothed_context"] = [
            bool(config.context_smoothing != "none" and any(
                token in {str(row.get("field_x", "")), str(row.get("field_y", ""))}
                for token in ("H", "V", "H_smoothed", "V_smoothed")
            ))
            for _, row in continuous.iterrows()
        ]
        continuous["context_smoothing_method"] = config.context_smoothing
        continuous["smoothing_warning"] = [
            SMOOTHING_LIMITATION if used else ""
            for used in continuous["uses_graph_smoothed_context"]
        ]
    continuous.insert(0, "sample", sample)
    continuous.to_csv(out["neighborhood"] / f"{graph_prefix}continuous_edge_statistics.csv", index=False)
    if variable_records:
        semantics_table(variable_records).drop_duplicates(subset=["source_column", "user_confirmed_analysis_mode"]).to_csv(
            out["variable_semantics"] / f"variable_semantics_{sample}.csv", index=False
        )
        _write_json(out["variable_semantics"] / f"variable_semantics_{sample}.json", {"sample": sample, "variables": variable_records})
    _write_json(out["neighborhood"] / f"{graph_prefix}permutation_parameters.json", {
        "sample": sample,
        "permutations": config.permutations,
        "seed": config.seed,
        "empirical_p_value": "(extreme_count + 1) / (n_permutations + 1)",
        "slide_wise_permutation": True,
        "label_source": config.label_source,
        "permutation_scope": config.permutation_scope,
        "stratification_column": config.stratification_column if config.permutation_scope == "within_user_strata" else "",
        "tissue_only_restriction_requested": config.tissue_only_restriction,
        "tissue_only_restriction_applied": tissue_mask is not None,
        "fdr_scope": FDR_SCOPE,
        "fdr_not_across": ["all_samples", "all_graph_types", "all_analysis_families", "all_robustness_runs"],
        "permutation_limitation": PERMUTATION_LIMITATION,
        "spatial_association_limitation": SPATIAL_ASSOCIATION_LIMITATION,
        "context_smoothing_method": config.context_smoothing,
        "context_smoothing_warning": SMOOTHING_LIMITATION if config.context_smoothing != "none" else "",
    })
    robustness_graph_csv = ""
    robustness_direction_csv = ""
    if config.run_graph_robustness:
        try:
            graph_robustness, direction_stability = compare_graph_robustness(
                adata,
                config.graph,
                labels,
                permutations=config.permutations,
                seed=config.seed,
                permutation_scope=config.permutation_scope,
                strata=strata,
                tissue_mask=tissue_mask,
            )
            robustness_graph_path = out["robustness"] / f"graph_robustness_summary_{sample}.csv"
            robustness_direction_path = out["robustness"] / f"association_direction_stability_{sample}.csv"
            graph_robustness.to_csv(robustness_graph_path, index=False)
            direction_stability.to_csv(robustness_direction_path, index=False)
            robustness_graph_csv = str(robustness_graph_path)
            robustness_direction_csv = str(robustness_direction_path)
        except Exception as exc:
            if progress:
                progress(f"  Optional graph robustness comparison skipped: {exc}")
    annotated_path = ""
    if config.write_annotated_h5ad:
        annotated_path = str(out["annotated"] / f"{sample}_spatialtx_v0_4_graph_context.h5ad")
        adata.write_h5ad(annotated_path, compression="gzip")
    return {
        "sample": sample,
        "source_h5ad": str(source.resolve()),
        "status": "ok",
        "graph_method": result.method,
        "graph_requested_method": result.metadata.get("requested_method", result.method),
        "graph_effective_method": result.metadata.get("effective_method", result.method),
        "graph_fallback_used": result.metadata.get("fallback_used", False),
        "n_nodes": result.qc.get("n_nodes"),
        "n_edges": result.qc.get("n_edges"),
        "isolated_fraction": result.qc.get("isolated_fraction"),
        "n_connected_components": result.qc.get("n_connected_components"),
        "categorical_csv": str(out["neighborhood"] / f"{graph_prefix}categorical_enrichment.csv"),
        "binary_csv": str(out["neighborhood"] / f"{graph_prefix}binary_mask_association.csv"),
        "same_spot_overlap_csv": str(out["neighborhood"] / f"{graph_prefix}same_spot_overlap.csv"),
        "neighboring_spot_association_csv": str(out["neighborhood"] / f"{graph_prefix}neighboring_spot_association.csv"),
        "continuous_csv": str(out["neighborhood"] / f"{graph_prefix}continuous_edge_statistics.csv"),
        "input_audit_json": str(out["input_audit"] / f"input_audit_{sample}.json"),
        "input_audit_csv": str(out["input_audit"] / f"input_audit_{sample}.csv"),
        "graph_robustness_csv": robustness_graph_csv,
        "association_direction_stability_csv": robustness_direction_csv,
        "annotated_h5ad": annotated_path,
        "warnings": "; ".join(result.warnings),
    }


def run_spatial_graph_neighborhood_batch(
    paths: list[Path],
    output_root: str | Path,
    c_genes: list[str],
    s_genes: list[str],
    config: SpatialGraphAnalysisConfig,
    progress: Progress | None = None,
    cancel_event: Event | None = None,
) -> tuple[Path, pd.DataFrame]:
    if not paths:
        raise ValueError("Select at least one h5ad sample.")
    if int(config.permutations) < 1:
        raise ValueError("Neighborhood permutations must be at least 1.")
    if int(config.max_plot_edges) < 1:
        raise ValueError("Maximum plotted edges must be at least 1.")
    for name, value in (("H high quantile", config.h_high_quantile), ("V high quantile", config.v_high_quantile)):
        if not 0 < float(value) < 1:
            raise ValueError(f"{name} must be between 0 and 1.")
    if not 0 <= float(config.min_gene_coverage) <= 1:
        raise ValueError("Minimum context-gene coverage must be between 0 and 1.")
    if not 0 <= float(config.context_min_spot_fraction) <= 1:
        raise ValueError("Minimum context-gene spot fraction must be between 0 and 1.")
    if config.permutation_scope not in {"whole_slide", "within_connected_components", "within_user_strata"}:
        raise ValueError("Unsupported permutation scope.")
    if config.permutation_scope == "within_user_strata" and not config.stratification_column:
        raise ValueError("within_user_strata requires a stratification column.")
    validation = validate_gene_programs(c_genes, s_genes, mode="graph")
    c_genes = validation.normalized_c_genes
    s_genes = validation.normalized_s_genes
    config.c_genes = list(c_genes)
    config.s_genes = list(s_genes)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_dir = Path(output_root).expanduser().resolve() / f"spatial_graph_neighborhood_{stamp}"
    out = _ensure_dirs(run_dir)
    _write_json(run_dir / "run_parameters.json", {
        "application": "SpatialTX Studio Desktop",
        "application_version": __version__,
        "analysis_module": "Spatial Graph and Neighborhood Analysis",
        "analysis_module_version": MODULE_VERSION,
        "created": dt.datetime.now().isoformat(timespec="seconds"),
        "graph_parameters": asdict(config.graph),
        "context_parameters": asdict(config),
        "Cx_genes": list(c_genes),
        "Sx_genes": list(s_genes),
        "gene_program_validation": validation.to_provenance(),
        "output_path": str(run_dir.resolve()),
        "analysis_status": "experimental",
        "spatial_association_limitation": SPATIAL_ASSOCIATION_LIMITATION,
        "permutation_limitation": PERMUTATION_LIMITATION,
        "fdr_scope": FDR_SCOPE,
        "context_smoothing_warning": SMOOTHING_LIMITATION if config.context_smoothing != "none" else "",
    })
    rows: list[dict] = []
    for number, source in enumerate(paths, 1):
        if cancel_event is not None and cancel_event.is_set():
            rows.append({"sample": source.stem, "source_h5ad": str(source.resolve()), "status": "cancelled"})
            break
        if progress:
            progress(f"[{number}/{len(paths)}] Spatial graph and neighborhood: {source.name}")
        try:
            rows.append(run_spatial_graph_neighborhood_sample(Path(source), out, config, progress))
        except Exception as exc:
            rows.append({"sample": Path(source).stem, "source_h5ad": str(Path(source).resolve()), "status": f"error: {exc}"})
            if progress:
                progress(f"  Error: {exc}")
    manifest = pd.DataFrame(rows)
    manifest.to_csv(run_dir / "spatial_graph_neighborhood_manifest.csv", index=False)
    ok = manifest[manifest["status"].eq("ok")] if "status" in manifest else pd.DataFrame()
    if len(ok):
        ok.to_csv(run_dir / "combined_cohort_summary.csv", index=False)
    if progress:
        progress(f"Completed: {run_dir}")
    return run_dir, manifest
