"""Canonical FRAME2.6/C-S CLI wrapper.

The former independent CLI implementation is retained in ``legacy_frame26.py``
for audit history only.  Public CLI execution delegates to the same in-memory
engine used by SpatialTX Studio Desktop's Main Mapper.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from spatialtx_desktop.gene_program_validation import validate_gene_programs
from spatialtx_desktop.workflow import (
    ScoringOptions,
    _read_h5ad,
    parse_gene_text,
    save_spatial_map,
    score_adata,
)

from .gene_program import select_gene_programs


C_GENES_FIXED = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG"]
B_GENES_FIXED = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN"]


@dataclass(frozen=True)
class Frame26Config:
    """Compatibility view of public FRAME2.6 thresholds.

    Legacy grid-only parameters remain readable in run metadata, but the
    canonical engine controls scoring, smoothing, gradients, masks, and regime
    calls.
    """

    sigma: float = 0.8
    qscale_low: float = 0.05
    qscale_high: float = 0.95
    c_q: float = 0.80
    b_q: float = 0.80
    g_q: float = 0.60
    transition_radius_spots: int = 2
    diffuse_exclude_interface_buffer: int = 1
    eps: float = 1e-9


def frame26_config_from_dict(config: dict | None) -> Frame26Config:
    if not config:
        return Frame26Config()
    frame_cfg = config.get("frame26", {})
    thresholds = config.get("thresholds", {})
    smoothing = config.get("smoothing", {})
    return Frame26Config(
        sigma=float(frame_cfg.get("sigma", smoothing.get("sigma", Frame26Config.sigma))),
        qscale_low=float(frame_cfg.get("qscale_low", thresholds.get("qscale_low", Frame26Config.qscale_low))),
        qscale_high=float(frame_cfg.get("qscale_high", thresholds.get("qscale_high", Frame26Config.qscale_high))),
        c_q=float(frame_cfg.get("c_q", thresholds.get("c_quantile", Frame26Config.c_q))),
        b_q=float(frame_cfg.get("b_q", thresholds.get("b_quantile", Frame26Config.b_q))),
        g_q=float(frame_cfg.get("g_q", thresholds.get("g_quantile", Frame26Config.g_q))),
        transition_radius_spots=int(
            frame_cfg.get(
                "transition_radius_spots",
                config.get("transition_zone", {}).get("radius_spots", Frame26Config.transition_radius_spots),
            )
        ),
        diffuse_exclude_interface_buffer=int(
            frame_cfg.get("diffuse_exclude_interface_buffer", Frame26Config.diffuse_exclude_interface_buffer)
        ),
        eps=float(frame_cfg.get("eps", Frame26Config.eps)),
    )


def _scoring_options_from_dict(config: dict | None) -> ScoringOptions:
    config = config or {}
    scoring = dict(config.get("scoring_options", {}))
    smoothing = config.get("smoothing", {})
    smoothing_mode = str(scoring.get("smoothing_mode", smoothing.get("mode", "none"))).lower()
    if smoothing_mode not in {"none", "knn_mean", "gaussian"}:
        # Historical values such as "standard" described the removed legacy
        # grid engine and must not silently select a different canonical path.
        smoothing_mode = "none"
    normalization_mode = str(scoring.get("normalization_mode", "raw_mean")).lower()
    if normalization_mode not in {"raw_mean", "z_score", "rank_quantile"}:
        normalization_mode = "raw_mean"
    return ScoringOptions(
        smoothing_mode=smoothing_mode,
        smoothing_k=int(scoring.get("smoothing_k", 6)),
        gaussian_sigma=float(scoring.get("gaussian_sigma", smoothing.get("sigma", 0.0)) or 0.0),
        normalization_mode=normalization_mode,
        perturbation_check=bool(scoring.get("perturbation_check", False)),
        parameter_log_export=bool(scoring.get("parameter_log_export", True)),
    )


def _selected_programs(adata, config: dict | None, mode: str) -> tuple[list[str], list[str], str]:
    if config:
        selected_c, selected_s, note = select_gene_programs(adata, config, mode)
    else:
        names = {str(name).upper(): str(name) for name in adata.var_names}
        fixed = validate_gene_programs(C_GENES_FIXED, B_GENES_FIXED, mode="fixed")
        selected_c = [gene for gene in fixed.normalized_c_genes if gene in names]
        selected_s = [gene for gene in fixed.normalized_s_genes if gene in names]
        note = "fixed_gene_program"
    return parse_gene_text(selected_c), parse_gene_text(selected_s), str(note)


def _legacy_compatible_summary(
    metrics: dict,
    fields: dict,
    cfg: Frame26Config,
    mode: str,
    selected_c: list[str],
    selected_s: list[str],
    selection_note: str,
) -> pd.DataFrame:
    C = np.asarray(fields["C"], dtype=float)
    S = np.asarray(fields["S"], dtype=float)
    G = np.asarray(fields["G"], dtype=float)
    high_g = np.asarray(fields.get("high_g", np.zeros(len(C), dtype=bool)), dtype=bool)
    legacy = {
        "mode": mode,
        "sigma": cfg.sigma,
        "qscale_low": cfg.qscale_low,
        "qscale_high": cfg.qscale_high,
        "c_q": cfg.c_q,
        "b_q": cfg.b_q,
        "g_q": cfg.g_q,
        "transition_radius_spots": cfg.transition_radius_spots,
        "c_threshold": float(np.nanquantile(C, cfg.c_q)),
        "b_threshold": float(np.nanquantile(S, cfg.b_q)),
        "g_threshold": float(np.nanquantile(G, cfg.g_q)) if np.isfinite(G).any() else np.nan,
        "interface_spot_count": metrics.get("interface_spots"),
        "transition_spot_count": metrics.get("diffuse_spots"),
        "diffuse_transition_spot_count": metrics.get("diffuse_spots"),
        "transition_fraction": metrics.get("diffuse_fraction"),
        "diffuse_transition_fraction": metrics.get("diffuse_fraction"),
        "g_high_fraction": float(np.mean(high_g)) if len(high_g) else np.nan,
        "r_variance": float(metrics.get("R_sd", np.nan)) ** 2,
        "diffuse_fraction": metrics.get("diffuse_fraction"),
        "regime_confidence": np.nan,
        "warnings": metrics.get("QC_notes", ""),
        "selection_note": selection_note,
        "selected_C": ",".join(selected_c),
        "selected_B": ",".join(selected_s),
        "legacy_transition_alias_note": (
            "transition_spot_count aliases the canonical diffuse-transition count; "
            "the retired independent grid-shell engine is not recomputed"
        ),
        "common_engine": "spatialtx_desktop.workflow.score_adata",
    }
    return pd.DataFrame([{**metrics, **legacy}])


def _annotate_compatibility_columns(adata, fields: dict) -> None:
    C = np.asarray(fields["C"], dtype=float)
    S = np.asarray(fields["S"], dtype=float)
    R = np.asarray(fields["R"], dtype=float)
    G = np.asarray(fields["G"], dtype=float)
    interface = np.asarray(fields["interface"], dtype=bool)
    diffuse = np.asarray(fields["diffuse"], dtype=bool)
    high_g = np.asarray(fields.get("high_g", np.zeros(len(C), dtype=bool)), dtype=bool)
    adata.obs["frame26_C"] = C
    adata.obs["frame26_B"] = S
    adata.obs["frame26_R"] = R
    adata.obs["frame26_C_smooth"] = C
    adata.obs["frame26_B_smooth"] = S
    adata.obs["frame26_G"] = G
    adata.obs["frame26_GR"] = G
    adata.obs["frame26_interface_like"] = interface.astype(int)
    adata.obs["frame26_diffuse_transition"] = diffuse.astype(int)
    adata.obs["frame26_high_grad"] = high_g.astype(int)


def write_frame26_qc(path: Path, summary: pd.DataFrame) -> None:
    row = summary.iloc[0]
    spatial_available = row.get("spatial_results_status") == "available"
    qc = [
        {"check": "analysis", "status": "ok", "detail": "Canonical C/S balance-field workflow"},
        {"check": "common_engine", "status": "ok", "detail": str(row["common_engine"])},
        {"check": "gene_selection", "status": "ok", "detail": str(row["selection_note"])},
        {"check": "regime_label", "status": "ok" if spatial_available else "warning", "detail": str(row["regime_label"])},
        {"check": "spatial_results", "status": "ok" if spatial_available else "warning", "detail": str(row.get("spatial_qc_message", ""))},
        {"check": "research_use_only", "status": "notice", "detail": "Not intended for clinical diagnosis or treatment decision-making."},
    ]
    pd.DataFrame(qc).to_csv(path, index=False)


def run_frame26(
    input_path: str | Path,
    sample_name: str,
    output_dir: str | Path,
    mode: str = "fixed",
    cfg: Frame26Config | None = None,
    config: dict | None = None,
):
    """Run the CLI through the canonical Main Mapper scoring engine."""
    cfg = cfg or frame26_config_from_dict(config)
    output_dir = Path(output_dir)
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    adata = _read_h5ad(input_path)
    selected_c, selected_s, selection_note = _selected_programs(adata, config, mode)
    metrics, fields = score_adata(
        adata,
        selected_c,
        selected_s,
        cfg.c_q,
        cfg.b_q,
        cfg.g_q,
        options=_scoring_options_from_dict(config),
        source_path=input_path,
        sample_name=sample_name,
        gene_program_mode=mode,
    )
    _annotate_compatibility_columns(adata, fields)
    summary = _legacy_compatible_summary(metrics, fields, cfg, mode, selected_c, selected_s, selection_note)

    summary.to_csv(output_dir / "metrics.csv", index=False)
    summary.to_csv(output_dir / "frame26_summary.csv", index=False)
    write_frame26_qc(output_dir / "qc_report.csv", summary)
    pd.DataFrame(
        [{"program": "C", "canonical_program": "C", "gene": gene} for gene in selected_c]
        + [{"program": "B", "canonical_program": "S", "gene": gene} for gene in selected_s]
    ).to_csv(output_dir / "selected_genes.csv", index=False)
    if fields.get("spatial_available", False):
        interface_path = save_spatial_map(figures_dir / "interface_map.png", metrics, fields)
        shutil.copyfile(interface_path, figures_dir / "transition_zone_map.png")
    else:
        (figures_dir / "SPATIAL_RESULTS_UNAVAILABLE.txt").write_text(
            str(metrics.get("spatial_qc_message", "Spatial results unavailable.")), encoding="utf-8"
        )
    with (output_dir / "run_log.txt").open("w", encoding="utf-8") as handle:
        for key, value in summary.iloc[0].items():
            handle.write(f"{key}={value}\n")
    (output_dir / "engine_metadata.json").write_text(
        json.dumps(
            {
                "common_engine": "spatialtx_desktop.workflow.score_adata",
                "legacy_engine_status": "audit_only_not_called",
                "legacy_module": "spatialtx_studio.legacy_frame26",
                "C_definition": "canonical Main Mapper C field",
                "S_definition": "canonical Main Mapper S field; legacy B aliases are retained in output names",
                "R_definition": "C-S",
                "gene_program_validation": fields.get("gene_program_validation", {}),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return summary, adata
