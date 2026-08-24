from __future__ import annotations

import datetime as dt
import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__
from ..graph.metadata import json_safe
from .matching import infer_manifest_batch_mode, paired_group_matches, pairwise_match
from .comparative_normalization import (
    CENTERED_HV_WARNING,
    apply_pooled_hv_thresholds,
    build_hv_summary,
    build_metric_change_table,
    build_normalized_metrics_table,
    build_relative_changes_table,
    build_sample_scale_table,
    build_scale_warnings,
)
from .metrics import (
    analyze_sample,
    compute_delta_metrics,
    group_mean_deltas,
    regime_transition_table,
)
from .models import (
    EXPLORATORY_NOTICE,
    NO_REGISTRATION_NOTICE,
    ComparativeCancelled,
    ComparativeConfig,
    ComparativeRunResult,
    Progress,
    SampleRecord,
)
from .plotting import generate_comparative_figures
from .reporting import rules_based_summary, write_html_report, write_pdf_report
from .statistics import comparative_group_statistics
from .validation import preflight_records, valid_records_from_preflight, validate_record_structure


OUTPUT_TABLE_SCHEMAS = {
    "comparative_sample_metrics.csv": ["sample_id", "group", "analysis_status"],
    "comparative_delta_metrics.csv": ["comparison_id", "metric", "delta"],
    "comparative_group_statistics.csv": ["metric", "status"],
    "comparative_regime_transitions.csv": ["comparison_id", "regime_transition"],
    "comparative_warnings.csv": ["warning_scope", "warning_code", "message"],
    "comparative_run_manifest.csv": ["sample_id", "status", "error"],
    "comparative_sample_scale.csv": ["sample_id", "group", "n_valid_spots"],
    "comparative_metric_change_table.csv": ["metric_name", "reference_value", "target_value", "raw_delta"],
    "comparative_normalized_metrics.csv": ["sample_id", "group"],
    "comparative_relative_changes.csv": ["metric_name", "symmetric_percent_change"],
    "comparative_scale_warnings.csv": ["comparison_id", "warning_code", "message"],
    "comparative_HV_summary.csv": ["axis", "metric", "observational_only"],
}


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(payload), indent=2, ensure_ascii=False), encoding="utf-8")


def _git_commit(start: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(start), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return "unavailable"
    return result.stdout.strip() if result.returncode == 0 and result.stdout.strip() else "unavailable"


def _environment_summary() -> dict:
    packages = {}
    for name in ("anndata", "numpy", "pandas", "scipy", "matplotlib", "h5py"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = "not installed"
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "processor": platform.processor(),
        "packages": packages,
        "cpu_only_compatible": True,
        "gpu_required": False,
    }


def _empty_table(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=columns)


def _write_csv(path: Path, table: pd.DataFrame, required_columns: list[str]) -> None:
    output = table.copy()
    for column in required_columns:
        if column not in output:
            output[column] = pd.Series(dtype=object)
    output.to_csv(path, index=False)


def _common_metadata(config: ComparativeConfig, effective_mode: str, timestamp: str, commit: str) -> dict:
    return {
        "spatialtx_version": f"v{__version__}",
        "analysis_timestamp": timestamp,
        "comparison_mode": config.mode,
        "effective_statistical_design": effective_mode,
        "reference_definition": config.reference,
        "target_definition": config.target,
        "seed": int(config.seed),
        "C_gene_program_json": json.dumps(config.c_genes, ensure_ascii=False),
        "S_gene_program_json": json.dumps(config.s_genes, ensure_ascii=False),
        "graph_settings_json": json.dumps(json_safe(config.graph_settings), ensure_ascii=False, sort_keys=True),
        "thresholds_json": json.dumps({"C": config.c_q, "S": config.s_q, "G": config.g_q}, sort_keys=True),
        "scoring_options_json": json.dumps(json_safe(config.scoring_options), ensure_ascii=False, sort_keys=True),
        "fdr_method": config.fdr_method,
        "H_V_observational_only": True,
        "git_commit": commit,
        "registration_status": "not_performed",
        "registration_notice": NO_REGISTRATION_NOTICE,
        "analysis_scope": "exploratory_non_diagnostic",
    }


def _annotate(table: pd.DataFrame, metadata: dict) -> pd.DataFrame:
    result = table.copy()
    for name, value in reversed(list(metadata.items())):
        if name not in result:
            result.insert(0, name, value)
    return result


def _sample_warning_rows(sample_metrics: pd.DataFrame) -> list[dict]:
    rows: list[dict] = []
    for _, sample in sample_metrics.iterrows():
        messages = [item.strip() for item in str(sample.get("warning_messages", "")).split(";") if item.strip()]
        for message in dict.fromkeys(messages):
            rows.append({
                "warning_scope": "sample",
                "sample_id": sample.get("sample_id", ""),
                "warning_code": "sample_or_context_qc",
                "message": message,
            })
    return rows


def _write_partial(run_dir: Path, manifest: pd.DataFrame, sample_rows: list[dict], log_lines: list[str]) -> None:
    _write_csv(
        run_dir / "comparative_run_manifest.csv",
        manifest,
        OUTPUT_TABLE_SCHEMAS["comparative_run_manifest.csv"],
    )
    _write_csv(
        run_dir / "comparative_sample_metrics.csv",
        pd.DataFrame(sample_rows),
        OUTPUT_TABLE_SCHEMAS["comparative_sample_metrics.csv"],
    )
    logs = run_dir / "comparative_logs"
    logs.mkdir(parents=True, exist_ok=True)
    (logs / "run.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")


def _prepare_manifest(preflight: pd.DataFrame) -> pd.DataFrame:
    manifest = preflight.copy()
    manifest["status"] = np.where(manifest["validation_status"].eq("ok"), "pending", "failed_validation")
    manifest["error"] = manifest["validation_error"]
    manifest["input_sha256"] = ""
    manifest["cache_status"] = ""
    return manifest


def _resolve_effective_mode(config: ComparativeConfig, records: list[SampleRecord]) -> str:
    if config.mode == "manifest_batch":
        return infer_manifest_batch_mode(records, config.reference, config.target)
    return config.mode


def _matched_rows_for_mode(
    effective_mode: str,
    records: list[SampleRecord],
    config: ComparativeConfig,
) -> pd.DataFrame | None:
    if effective_mode == "pairwise":
        return pairwise_match(records, config.reference, config.target)
    if effective_mode == "paired":
        return paired_group_matches(records, config.reference, config.target)
    return None


def _group_summary(table: pd.DataFrame, group: str) -> pd.DataFrame:
    return table.loc[table["group"].eq(group)].copy()


def run_comparative_analysis(
    records: list[SampleRecord],
    config: ComparativeConfig,
    output_root: str | Path,
    progress: Progress | None = None,
    cancel_event=None,
) -> ComparativeRunResult:
    """Run sample-summary comparison without spot-wise subtraction or registration assumptions."""
    validate_record_structure(records, config)
    if progress:
        progress("Validating comparative H5AD inputs...")
    preflight = preflight_records(records, config)
    valid_records = valid_records_from_preflight(records, preflight)
    if len(valid_records) < 2:
        errors = "; ".join(
            f"{row['sample_id']}: {row['validation_error']}"
            for _, row in preflight.loc[~preflight["validation_status"].eq("ok")].iterrows()
        )
        raise ValueError(f"Fewer than two valid spatial H5AD samples remain. {errors}")
    if config.mode == "pairwise" and len(valid_records) != len(records):
        errors = "; ".join(
            f"{row['sample_id']}: {row['validation_error']}"
            for _, row in preflight.loc[~preflight["validation_status"].eq("ok")].iterrows()
        )
        raise ValueError(f"Pairwise comparison cannot continue with an invalid sample. {errors}")
    effective_mode = _resolve_effective_mode(config, valid_records)
    initial_matches = _matched_rows_for_mode(effective_mode, valid_records, config)

    output_root = Path(output_root).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    timestamp = dt.datetime.now().isoformat(timespec="seconds")
    run_dir = output_root / f"comparative_analysis_{stamp}"
    figures_dir = run_dir / "comparative_figures"
    logs_dir = run_dir / "comparative_logs"
    cache_dir = output_root / ".spatialtx_comparative_cache"
    figures_dir.mkdir(parents=True)
    logs_dir.mkdir(parents=True)
    commit = _git_commit(Path(__file__).resolve().parents[3])
    log_lines = [
        f"SpatialTX Studio v{__version__} Comparative Analysis",
        f"Timestamp: {timestamp}",
        f"Mode: {config.mode}; effective design: {effective_mode}",
        NO_REGISTRATION_NOTICE,
    ]
    manifest = _prepare_manifest(preflight)
    warning_rows: list[dict] = [{
        "warning_scope": "run",
        "sample_id": "",
        "warning_code": "no_coordinate_registration",
        "message": NO_REGISTRATION_NOTICE,
    }]
    for _, row in preflight.loc[~preflight["validation_status"].eq("ok")].iterrows():
        warning_rows.append({
            "warning_scope": "sample_validation",
            "sample_id": row["sample_id"],
            "warning_code": "invalid_sample_retained_in_manifest",
            "message": row["validation_error"],
        })
    if len(records) > int(config.large_batch_warning_count):
        warning_rows.append({
            "warning_scope": "run",
            "sample_id": "",
            "warning_code": "large_batch",
            "message": f"Large batch ({len(records)} samples); review memory and runtime before repeating.",
        })
    sample_rows: list[dict] = []
    fields_by_sample: dict[str, dict[str, np.ndarray]] = {}
    for number, record in enumerate(valid_records, 1):
        if cancel_event is not None and cancel_event.is_set():
            pending = manifest["status"].eq("pending")
            manifest.loc[pending, "status"] = "cancelled"
            manifest.loc[pending, "error"] = "Cancelled before sample analysis."
            log_lines.append("Cancellation requested; partial results written safely.")
            _write_partial(run_dir, manifest, sample_rows, log_lines)
            raise ComparativeCancelled(f"Comparative analysis cancelled. Partial results: {run_dir}")
        message = f"[{number}/{len(valid_records)}] Scoring {record.sample_id} with the canonical C/S engine..."
        log_lines.append(message)
        if progress:
            progress(message)
        try:
            metrics, fields, cache_status = analyze_sample(record, config, cache_dir)
            sample_rows.append(metrics)
            fields_by_sample[record.sample_id] = fields
            mask = manifest["sample_id"].astype(str).eq(record.sample_id)
            manifest.loc[mask, "status"] = "ok"
            manifest.loc[mask, "error"] = ""
            manifest.loc[mask, "input_sha256"] = metrics.get("input_sha256", "")
            manifest.loc[mask, "cache_status"] = cache_status
        except Exception as exc:
            mask = manifest["sample_id"].astype(str).eq(record.sample_id)
            manifest.loc[mask, "status"] = "failed_analysis"
            manifest.loc[mask, "error"] = str(exc)
            warning_rows.append({
                "warning_scope": "sample_analysis",
                "sample_id": record.sample_id,
                "warning_code": "analysis_failure",
                "message": str(exc),
            })
            log_lines.append(f"ERROR {record.sample_id}: {exc}")
        _write_partial(run_dir, manifest, sample_rows, log_lines)

    sample_metrics = pd.DataFrame(sample_rows)
    if len(sample_metrics) < 2:
        raise ValueError(f"Fewer than two samples completed comparative scoring. Partial results: {run_dir}")
    successful_ids = set(sample_metrics["sample_id"].astype(str))
    successful_records = [record for record in valid_records if record.sample_id in successful_ids]
    reference_records = [record for record in successful_records if record.group == config.reference]
    target_records = [record for record in successful_records if record.group == config.target]
    if not reference_records or not target_records:
        raise ValueError(f"At least one successful sample is required in both reference and target groups. Partial results: {run_dir}")
    sample_metrics = apply_pooled_hv_thresholds(
        sample_metrics,
        fields_by_sample,
        config.reference,
        config.target,
    )
    matches = initial_matches
    if matches is not None:
        complete = matches["reference_sample_id"].isin(successful_ids) & matches["target_sample_id"].isin(successful_ids)
        for _, dropped in matches.loc[~complete].iterrows():
            warning_rows.append({
                "warning_scope": "matching",
                "sample_id": "",
                "warning_code": "incomplete_pair_after_analysis_failure",
                "message": f"Pair/comparison {dropped['comparison_id']} was excluded from paired statistics because one member failed.",
            })
        matches = matches.loc[complete].reset_index(drop=True)
        if matches.empty:
            raise ValueError(f"No complete successful pair remains. Partial results: {run_dir}")

    if effective_mode in {"pairwise", "paired"}:
        assert matches is not None
        delta_metrics = compute_delta_metrics(sample_metrics, matches)
    else:
        delta_metrics = group_mean_deltas(sample_metrics, config.reference, config.target)
    metric_change_table = build_metric_change_table(delta_metrics)
    sample_scale = build_sample_scale_table(sample_metrics)
    normalized_metrics = build_normalized_metrics_table(sample_metrics)
    relative_changes = build_relative_changes_table(metric_change_table)
    scale_warnings = build_scale_warnings(metric_change_table)
    hv_summary = build_hv_summary(metric_change_table, sample_metrics)
    stats_design = "paired" if effective_mode == "paired" else "unpaired"
    group_statistics = comparative_group_statistics(
        sample_metrics,
        config,
        matches if stats_design == "paired" else None,
        effective_mode=stats_design,
    )
    if effective_mode == "pairwise" and len(group_statistics):
        group_statistics["comparison_design"] = "pairwise_descriptive"
    transitions = regime_transition_table(
        sample_metrics,
        matches if effective_mode in {"pairwise", "paired"} else None,
        config.reference,
        config.target,
        config.low_regime_confidence_threshold,
    )
    for group in (config.reference, config.target):
        count = int(sample_metrics["group"].eq(group).sum())
        if count < 3:
            warning_rows.append({
                "warning_scope": "group",
                "sample_id": "",
                "warning_code": "small_sample_size",
                "message": f"{group} has n={count}; group inference is unstable for n < 3.",
            })
    warning_rows.extend(_sample_warning_rows(sample_metrics))
    for _, warning in scale_warnings.iterrows():
        warning_rows.append({
            "warning_scope": "comparison_scale",
            "sample_id": "",
            "warning_code": warning.get("warning_code", "sample_scale"),
            "message": warning.get("message", ""),
        })
    if config.enable_h_expr or config.enable_v_expr:
        warning_rows.append({
            "warning_scope": "observational_context",
            "sample_id": "",
            "warning_code": "centered_hv_means_non_informative",
            "message": CENTERED_HV_WARNING,
        })

    reference_sample = reference_records[0].sample_id
    target_sample = target_records[0].sample_id
    if progress:
        progress("Generating fresh comparative figures...")
    figures = generate_comparative_figures(
        sample_metrics,
        delta_metrics,
        transitions,
        fields_by_sample,
        reference_sample,
        target_sample,
        figures_dir,
        config,
        metric_change_table=metric_change_table,
        sample_scale=sample_scale,
        scale_warnings=scale_warnings,
        hv_summary=hv_summary,
        effective_mode=effective_mode,
    )
    summary_text = rules_based_summary(
        sample_metrics,
        delta_metrics,
        transitions,
        config,
        effective_mode,
        metric_change_table=metric_change_table,
        scale_warnings=scale_warnings,
        hv_summary=hv_summary,
    )
    metadata = _common_metadata(config, effective_mode, timestamp, commit)
    environment = _environment_summary()
    failed_samples = manifest.loc[
        ~manifest["status"].eq("ok"), ["sample_id", "status", "error"]
    ].to_dict("records")
    metadata.update({
        "input_paths_json": json.dumps(
            [str(record.file_path) for record in records], ensure_ascii=False
        ),
        "input_hashes_json": json.dumps(
            {
                str(row["sample_id"]): str(row["input_sha256"])
                for _, row in manifest.iterrows()
                if str(row.get("input_sha256", ""))
            },
            ensure_ascii=False,
            sort_keys=True,
        ),
        "warnings_json": json.dumps(json_safe(warning_rows), ensure_ascii=False),
        "failed_samples_json": json.dumps(json_safe(failed_samples), ensure_ascii=False),
        "software_environment_json": json.dumps(json_safe(environment), ensure_ascii=False, sort_keys=True),
    })
    sample_metrics = _annotate(sample_metrics, metadata)
    delta_metrics = _annotate(delta_metrics, metadata)
    group_statistics = _annotate(group_statistics, metadata)
    transitions = _annotate(transitions, metadata)
    sample_scale = _annotate(sample_scale, metadata)
    metric_change_table = _annotate(metric_change_table, metadata)
    normalized_metrics = _annotate(normalized_metrics, metadata)
    relative_changes = _annotate(relative_changes, metadata)
    scale_warnings = _annotate(scale_warnings, metadata)
    hv_summary = _annotate(hv_summary, metadata)
    warnings_table = _annotate(pd.DataFrame(warning_rows), metadata)
    manifest = _annotate(manifest, metadata)

    tables = {
        "comparative_sample_metrics.csv": sample_metrics,
        "comparative_delta_metrics.csv": delta_metrics,
        "comparative_group_statistics.csv": group_statistics,
        "comparative_regime_transitions.csv": transitions,
        "comparative_warnings.csv": warnings_table,
        "comparative_run_manifest.csv": manifest,
        "comparative_sample_scale.csv": sample_scale,
        "comparative_metric_change_table.csv": metric_change_table,
        "comparative_normalized_metrics.csv": normalized_metrics,
        "comparative_relative_changes.csv": relative_changes,
        "comparative_scale_warnings.csv": scale_warnings,
        "comparative_HV_summary.csv": hv_summary,
    }
    for filename, table in tables.items():
        _write_csv(run_dir / filename, table, OUTPUT_TABLE_SCHEMAS[filename])
    parameters = {
        **metadata,
        "application": "SpatialTX Studio",
        "analysis_module": "Comparative Spatial Transition Analysis",
        "analysis_status": "development_exploratory",
        "config": config.to_dict(),
        "C_gene_program": config.c_genes,
        "S_gene_program": config.s_genes,
        "graph_settings": config.graph_settings,
        "thresholds": {"C": config.c_q, "S": config.s_q, "G": config.g_q},
        "inputs": [
            {
                **record.to_dict(),
                "input_sha256": str(
                    manifest.loc[manifest["sample_id"].astype(str).eq(record.sample_id), "input_sha256"].iloc[0]
                ),
                "status": str(manifest.loc[manifest["sample_id"].astype(str).eq(record.sample_id), "status"].iloc[0]),
            }
            for record in records
        ],
        "failed_samples": failed_samples,
        "warnings": warning_rows,
        "software_environment": environment,
        "spotwise_subtraction_performed": False,
        "H_V_observational_only": True,
        "candidate_discovery_added": False,
        "ligand_receptor_analysis_added": False,
        "QUBO_added_to_comparative_module": False,
        "AI_interpretation_added": False,
        "multi_axis_modeling_added": False,
        "exploratory_notice": EXPLORATORY_NOTICE,
    }
    _write_json(run_dir / "comparative_parameters.json", parameters)

    if effective_mode == "pairwise":
        _write_csv(run_dir / "sample_A_summary.csv", _group_summary(sample_metrics, config.reference), ["sample_id"])
        _write_csv(run_dir / "sample_B_summary.csv", _group_summary(sample_metrics, config.target), ["sample_id"])
        wide = delta_metrics.pivot_table(index="comparison_id", columns="delta_metric", values="delta", aggfunc="first").reset_index()
        _write_csv(run_dir / "pairwise_delta_summary.csv", _annotate(wide, metadata), ["comparison_id"])
    else:
        _write_csv(run_dir / "group_A_summary.csv", _group_summary(sample_metrics, config.reference), ["sample_id"])
        _write_csv(run_dir / "group_B_summary.csv", _group_summary(sample_metrics, config.target), ["sample_id"])
        effects = group_statistics[[column for column in group_statistics if column in {
            *metadata.keys(), "metric", "effect_size", "effect_size_method", "ci_95_low", "ci_95_high", "status"
        }]].copy()
        fdr = group_statistics[[column for column in group_statistics if column in {
            *metadata.keys(), "metric", "p_value", "adjusted_p_value_bh", "significant_fdr_0_05", "fdr_scope", "status"
        }]].copy()
        _write_csv(run_dir / "group_effect_sizes.csv", effects, ["metric", "effect_size"])
        _write_csv(run_dir / "group_fdr_results.csv", fdr, ["metric", "adjusted_p_value_bh"])

    report_tables = {
        "Sample metrics": sample_metrics,
        "Metric changes": delta_metrics,
        "Group statistics": group_statistics,
        "Operational regime transitions": transitions,
        "Warnings": warnings_table,
        "Sample scale": sample_scale,
        "Metric change table": metric_change_table,
        "Normalized topology metrics": normalized_metrics,
        "Relative changes": relative_changes,
        "H/V observational summary": hv_summary,
        "Run manifest": manifest,
    }
    write_html_report(run_dir / "comparative_summary_report.html", summary_text, report_tables, figures, config)
    write_pdf_report(run_dir / "comparative_summary_report.pdf", summary_text, figures, config, len(sample_metrics))
    log_lines.append(f"Completed: {run_dir}")
    log_lines.append(f"Inputs: {metadata['input_paths_json']}")
    log_lines.append(f"Input hashes: {metadata['input_hashes_json']}")
    log_lines.append(f"C gene program: {metadata['C_gene_program_json']}")
    log_lines.append(f"S gene program: {metadata['S_gene_program_json']}")
    log_lines.append(f"Graph settings: {metadata['graph_settings_json']}")
    log_lines.append(f"Thresholds: {metadata['thresholds_json']}; seed: {config.seed}")
    log_lines.append(f"Warnings: {metadata['warnings_json']}")
    log_lines.append(f"Failed samples: {metadata['failed_samples_json']}")
    log_lines.append(f"Environment: {metadata['software_environment_json']}")
    log_lines.append(EXPLORATORY_NOTICE)
    (logs_dir / "run.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    if progress:
        progress(f"Comparative analysis completed: {run_dir}")
    return ComparativeRunResult(
        run_dir=run_dir,
        sample_metrics=sample_metrics,
        delta_metrics=delta_metrics,
        group_statistics=group_statistics,
        regime_transitions=transitions,
        warnings=warnings_table,
        run_manifest=manifest,
        figures=figures,
        summary_text=summary_text,
        effective_mode=effective_mode,
        sample_scale=sample_scale,
        metric_change_table=metric_change_table,
        normalized_metrics=normalized_metrics,
        relative_changes=relative_changes,
        scale_warnings=scale_warnings,
        hv_summary=hv_summary,
    )
