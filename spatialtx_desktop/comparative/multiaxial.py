from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__


SITE_COMPARABILITY_VALUES = ("same_site", "different_site", "unknown_site")
SITE_SHIFT_WARNING = (
    "SITE-SHIFT WARNING: The paired specimens originate from different anatomical sites. "
    "Observed differences may reflect site-specific tissue composition or microenvironmental variation "
    "in addition to treatment-associated differences."
)


def _finite(value) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return np.nan
    return numeric if np.isfinite(numeric) else np.nan


def normalize_site_comparability(value: str | None) -> str:
    normalized = str(value or "unknown_site").strip().casefold().replace(" ", "_")
    if normalized not in SITE_COMPARABILITY_VALUES:
        raise ValueError(
            "Site comparability must be same_site, different_site, or unknown_site."
        )
    return normalized


def interpretation_confidence(comparability: str, site_comparability: str, status: str) -> str:
    if str(status) != "PASS" or str(comparability) == "Low":
        return "LOW"
    if str(comparability) == "Caution" or site_comparability == "different_site":
        return "CAUTION"
    return "GOOD"


def _change_word(delta, tolerance: float) -> str:
    value = _finite(delta)
    if not np.isfinite(value):
        return "not available"
    if value > tolerance:
        return "increased"
    if value < -tolerance:
        return "decreased"
    return "remained broadly stable"


def _context_axis_sentence(axis: str, row: dict) -> str:
    label = (
        "Hypoxia-associated expression context"
        if axis == "H"
        else "Vascular-associated expression proxy"
    )
    pre_status = str(row.get(f"pre_{axis}_context_status", "available"))
    post_status = str(row.get(f"post_{axis}_context_status", "available"))
    median_change = _change_word(row.get(f"delta_{axis}"), 0.05)
    high_change = _change_word(row.get(f"delta_{axis}_high_fraction"), 0.01)
    q90_change = _change_word(row.get(f"delta_{axis}_q90"), 0.05)
    if median_change == "not available":
        return f"{axis}: not available (Pre={pre_status}, Post={post_status})."
    if median_change == "remained broadly stable" and high_change not in {
        "not available",
        "remained broadly stable",
    }:
        return (
            f"{label} median remained broadly stable; the pair-pooled upper-tail high-context fraction "
            f"{high_change}."
        )
    details = [f"median {median_change}"]
    if q90_change != "not available":
        details.append(f"q90 {q90_change}")
    if high_change != "not available":
        details.append(f"pair-pooled high-context fraction {high_change}")
    return f"{label}: " + "; ".join(details) + "."


def context_change_lines(row: dict) -> tuple[str, str]:
    return _context_axis_sentence("H", row), _context_axis_sentence("V", row)


def spatial_change_sentence(row: dict) -> str:
    regime_pre = str(row.get("regime_pre", "")).strip()
    regime_post = str(row.get("regime_post", "")).strip()
    interface_delta = _finite(row.get("delta_interface_fraction"))
    diffuse_delta = _finite(row.get("delta_diffuse_fraction"))
    if regime_pre == "Type_A_candidate" and regime_post == "Type_B_candidate":
        return (
            "Localized interface-like organization was no longer assigned as the broad operational regime, "
            "while diffuse transition structure remained present."
        )
    if regime_pre and regime_pre == regime_post:
        if interface_delta < -0.01 and diffuse_delta > 0.01:
            return (
                "The localized-interface component decreased while diffuse transition increased, consistent "
                "with redistribution of transition organization without a regime-level change."
            )
        return "The localized transition regime remained broadly stable; metric changes are descriptive."
    if regime_pre and regime_post:
        return (
            f"The operational regime changed from {regime_pre} to {regime_post}; this remains a descriptive "
            "candidate transition."
        )
    return "The operational regime comparison was not available."


def site_warning_text(site_comparability: str) -> str:
    if site_comparability == "different_site":
        return SITE_SHIFT_WARNING
    if site_comparability == "unknown_site":
        return "Site comparability: unknown."
    return "Site comparability: same anatomical site reported."


def qc_aware_interpretation(row: dict) -> tuple[str, str]:
    comparability = str(row.get("comparability", "Low"))
    site = normalize_site_comparability(row.get("site_comparability", "unknown_site"))
    confidence = interpretation_confidence(comparability, site, str(row.get("status", "")))
    h_sentence, v_sentence = context_change_lines(row)
    sentences = [spatial_change_sentence(row), h_sentence, v_sentence]
    if comparability == "Low":
        sentences.append(
            "The pair has low technical/spatial comparability. Observed changes should not be attributed "
            "directly to treatment without independent validation."
        )
    elif comparability == "Caution":
        sentences.append("Comparability cautions should be reviewed before biological attribution.")
    else:
        sentences.append(
            "No configured comparability warning was detected; this does not establish biological causality."
        )
    sentences.append(site_warning_text(site))
    return confidence, " ".join(sentences)


_QC_LABELS = {
    "n_spots": "Spot count",
    "n_features": "Feature count",
    "median_detected_genes_per_spot": "Detected genes/spot",
    "median_total_counts": "Median counts",
    "q25_total_counts": "Q25 counts",
    "q75_total_counts": "Q75 counts",
    "in_tissue_fraction": "In-tissue fraction",
    "low_quality_fraction": "Low-quality fraction",
    "spatial_extent_area_proxy": "Spatial extent",
    "tissue_component_count": "Tissue components",
    "valid_spot_fraction": "Valid-spot fraction",
    "C_gene_coverage": "C gene coverage",
    "S_gene_coverage": "S gene coverage",
}


def _mismatch_display(row: pd.Series) -> str:
    metric = str(row.get("qc_metric", ""))
    label = _QC_LABELS.get(metric, metric.replace("_", " ").title())
    value = _finite(row.get("comparison_value"))
    if not np.isfinite(value):
        return str(row.get("reason", "")).strip() or f"{label}: unavailable"
    if metric in {
        "n_spots", "n_features", "median_detected_genes_per_spot", "median_total_counts",
        "q25_total_counts", "q75_total_counts", "spatial_extent_area_proxy", "tissue_component_count",
    }:
        return f"{label} {value:.2f}x"
    return f"{label} difference {value:.3f}"


def _primary_mismatches(qc: pd.DataFrame, limit: int = 3) -> list[str]:
    if qc.empty:
        return []
    selected = qc.loc[
        qc["severity"].astype(str).isin(("low", "caution"))
        & ~qc["category"].astype(str).isin(("pair_validation", "site_metadata", "composition_proxy_secondary"))
    ].copy()
    if selected.empty:
        return []
    selected["_severity_order"] = selected["severity"].map({"low": 0, "caution": 1}).fillna(2)
    selected["_primary_order"] = (~selected["primary_for_classification"].fillna(False).astype(bool)).astype(int)
    selected["_source_order"] = np.arange(len(selected))
    selected = selected.sort_values(["_severity_order", "_primary_order", "_source_order"])
    return [_mismatch_display(row) for _, row in selected.head(limit).iterrows()]


def build_comparative_qc_summary(
    pair_results: pd.DataFrame,
    comparability_qc: pd.DataFrame,
    interpretation: pd.DataFrame,
) -> pd.DataFrame:
    interpretation_by_pair = {
        str(row["pair_label"]): row for _, row in interpretation.iterrows()
    }
    rows: list[dict] = []
    for _, pair in pair_results.iterrows():
        label = str(pair.get("pair_label", ""))
        qc = (
            comparability_qc.loc[comparability_qc["pair_label"].astype(str).eq(label)].copy()
            if not comparability_qc.empty and "pair_label" in comparability_qc
            else pd.DataFrame()
        )
        mismatches = _primary_mismatches(qc)
        detail = interpretation_by_pair.get(label, {})
        site = normalize_site_comparability(pair.get("site_comparability", "unknown_site"))
        rows.append({
            "pair_id": label,
            "sample_pre": Path(str(pair.get("pre_file", ""))).stem,
            "sample_post": Path(str(pair.get("post_file", ""))).stem,
            "comparability_status": pair.get("comparability", "Low"),
            "comparability_reason": pair.get("comparability_reasons", ""),
            "primary_mismatch_1": mismatches[0] if len(mismatches) > 0 else "",
            "primary_mismatch_2": mismatches[1] if len(mismatches) > 1 else "",
            "primary_mismatch_3": mismatches[2] if len(mismatches) > 2 else "",
            "primary_mismatch_summary": "; ".join(mismatches),
            "site_comparability": site,
            "site_shift_warning": SITE_SHIFT_WARNING if site == "different_site" else "",
            "interpretation_confidence": detail.get(
                "interpretation_confidence",
                interpretation_confidence(str(pair.get("comparability", "Low")), site, str(pair.get("status", ""))),
            ),
        })
    return pd.DataFrame(rows)


def build_multiaxial_pair_summary(
    pair_results: pd.DataFrame,
    interpretation: pd.DataFrame,
) -> pd.DataFrame:
    interpretation_by_pair = {
        str(row["pair_label"]): row for _, row in interpretation.iterrows()
    }
    rows: list[dict] = []
    for _, pair in pair_results.iterrows():
        label = str(pair.get("pair_label", ""))
        detail = interpretation_by_pair.get(label, {})
        rows.append({
            "pair_id": label,
            "sample_pre": Path(str(pair.get("pre_file", ""))).stem,
            "sample_post": Path(str(pair.get("post_file", ""))).stem,
            "regime_pre": pair.get("regime_pre", ""),
            "regime_post": pair.get("regime_post", ""),
            "regime_change": pair.get("regime_transition", ""),
            "interface_pre": pair.get("pre_interface_fraction", np.nan),
            "interface_post": pair.get("post_interface_fraction", np.nan),
            "delta_interface": pair.get("delta_interface_fraction", np.nan),
            "diffuse_pre": pair.get("pre_diffuse_fraction", np.nan),
            "diffuse_post": pair.get("post_diffuse_fraction", np.nan),
            "delta_diffuse": pair.get("delta_diffuse_fraction", np.nan),
            "burden_pre": pair.get("pre_transition_burden", np.nan),
            "burden_post": pair.get("post_transition_burden", np.nan),
            "delta_burden": pair.get("delta_transition_burden", np.nan),
            "cs_balance_pre": pair.get("pre_R", np.nan),
            "cs_balance_post": pair.get("post_R", np.nan),
            "delta_cs_balance": pair.get("delta_R", np.nan),
            "H_pre": pair.get("pre_H", np.nan),
            "H_post": pair.get("post_H", np.nan),
            "delta_H": pair.get("delta_H", np.nan),
            "H_raw_mean_pre": pair.get("pre_H_raw_mean", np.nan),
            "H_raw_mean_post": pair.get("post_H_raw_mean", np.nan),
            "delta_H_raw_mean": pair.get("delta_H_raw_mean", np.nan),
            "H_q75_pre": pair.get("pre_H_q75", np.nan),
            "H_q75_post": pair.get("post_H_q75", np.nan),
            "delta_H_q75": pair.get("delta_H_q75", np.nan),
            "H_q90_pre": pair.get("pre_H_q90", np.nan),
            "H_q90_post": pair.get("post_H_q90", np.nan),
            "delta_H_q90": pair.get("delta_H_q90", np.nan),
            "H_high_fraction_pre": pair.get("pre_H_high_fraction", np.nan),
            "H_high_fraction_post": pair.get("post_H_high_fraction", np.nan),
            "delta_H_high_fraction": pair.get("delta_H_high_fraction", np.nan),
            "H_local_fraction_pre": pair.get("pre_H_local_hotspot_fraction", np.nan),
            "H_local_fraction_post": pair.get("post_H_local_hotspot_fraction", np.nan),
            "delta_H_local_fraction": pair.get("delta_H_local_hotspot_fraction", np.nan),
            "H_transition_enrichment_pre": pair.get("pre_H_transition_enrichment", np.nan),
            "H_transition_enrichment_post": pair.get("post_H_transition_enrichment", np.nan),
            "delta_H_transition_enrichment": pair.get("delta_H_transition_enrichment", np.nan),
            "H_coefficient_of_variation_pre": pair.get("pre_H_coefficient_of_variation", np.nan),
            "H_coefficient_of_variation_post": pair.get("post_H_coefficient_of_variation", np.nan),
            "H_pair_pooled_q90": pair.get("H_pair_pooled_q90", np.nan),
            "H_gene_coverage_pre": pair.get("pre_H_gene_coverage", np.nan),
            "H_gene_coverage_post": pair.get("post_H_gene_coverage", np.nan),
            "H_context_status_pre": pair.get("pre_H_context_status", "calculation_error"),
            "H_context_status_post": pair.get("post_H_context_status", "calculation_error"),
            "V_pre": pair.get("pre_V", np.nan),
            "V_post": pair.get("post_V", np.nan),
            "delta_V": pair.get("delta_V", np.nan),
            "V_raw_mean_pre": pair.get("pre_V_raw_mean", np.nan),
            "V_raw_mean_post": pair.get("post_V_raw_mean", np.nan),
            "delta_V_raw_mean": pair.get("delta_V_raw_mean", np.nan),
            "V_q75_pre": pair.get("pre_V_q75", np.nan),
            "V_q75_post": pair.get("post_V_q75", np.nan),
            "delta_V_q75": pair.get("delta_V_q75", np.nan),
            "V_q90_pre": pair.get("pre_V_q90", np.nan),
            "V_q90_post": pair.get("post_V_q90", np.nan),
            "delta_V_q90": pair.get("delta_V_q90", np.nan),
            "V_high_fraction_pre": pair.get("pre_V_high_fraction", np.nan),
            "V_high_fraction_post": pair.get("post_V_high_fraction", np.nan),
            "delta_V_high_fraction": pair.get("delta_V_high_fraction", np.nan),
            "V_local_fraction_pre": pair.get("pre_V_local_hotspot_fraction", np.nan),
            "V_local_fraction_post": pair.get("post_V_local_hotspot_fraction", np.nan),
            "delta_V_local_fraction": pair.get("delta_V_local_hotspot_fraction", np.nan),
            "V_transition_enrichment_pre": pair.get("pre_V_transition_enrichment", np.nan),
            "V_transition_enrichment_post": pair.get("post_V_transition_enrichment", np.nan),
            "delta_V_transition_enrichment": pair.get("delta_V_transition_enrichment", np.nan),
            "V_coefficient_of_variation_pre": pair.get("pre_V_coefficient_of_variation", np.nan),
            "V_coefficient_of_variation_post": pair.get("post_V_coefficient_of_variation", np.nan),
            "V_pair_pooled_q90": pair.get("V_pair_pooled_q90", np.nan),
            "V_gene_coverage_pre": pair.get("pre_V_gene_coverage", np.nan),
            "V_gene_coverage_post": pair.get("post_V_gene_coverage", np.nan),
            "V_context_status_pre": pair.get("pre_V_context_status", "calculation_error"),
            "V_context_status_post": pair.get("post_V_context_status", "calculation_error"),
            "comparability_status": pair.get("comparability", "Low"),
            "comparability_reason": pair.get("comparability_reasons", ""),
            "site_comparability": pair.get("site_comparability", "unknown_site"),
            "interpretation_confidence": detail.get("interpretation_confidence", "LOW"),
            "interpretation_note": detail.get("qc_aware_interpretation_note", ""),
            "status": pair.get("status", ""),
            "error": pair.get("error", ""),
        })
    return pd.DataFrame(rows)


def _format_value(value: float) -> str:
    return f"{value:+.4g}" if np.isfinite(value) else "NA"


def _label_bars(ax, bars, values: list[float]) -> None:
    finite = [abs(value) for value in values if np.isfinite(value)]
    pad = max(finite, default=1.0) * 0.035
    for bar, value in zip(bars, values):
        if not np.isfinite(value):
            continue
        ax.text(
            value + (pad if value >= 0 else -pad),
            bar.get_y() + bar.get_height() / 2,
            _format_value(value),
            ha="left" if value >= 0 else "right",
            va="center",
            fontsize=8,
        )


def _center_zero_axis(ax, values: list[float]) -> None:
    finite = [abs(value) for value in values if np.isfinite(value)]
    limit = max(finite, default=1.0)
    if limit <= 0:
        limit = 1.0
    ax.set_xlim(-1.50 * limit, 1.50 * limit)


def _context_plot_layers(summary: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return the two independently scaled H/V layers used by the overview figure."""
    return {
        "median": (
            pd.to_numeric(summary.get("delta_H"), errors="coerce").to_numpy(dtype=float),
            pd.to_numeric(summary.get("delta_V"), errors="coerce").to_numpy(dtype=float),
        ),
        "pair_pooled_high_fraction": (
            pd.to_numeric(summary.get("delta_H_high_fraction"), errors="coerce").to_numpy(dtype=float),
            pd.to_numeric(summary.get("delta_V_high_fraction"), errors="coerce").to_numpy(dtype=float),
        ),
    }


def plot_multiaxial_overview(summary: pd.DataFrame, output: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    passed = summary.loc[summary.get("status", pd.Series(dtype=str)).eq("PASS")].copy()
    fig_height = max(7.0, 1.20 * max(1, len(passed)))
    fig = plt.figure(figsize=(20.5, fig_height), layout="constrained")
    outer = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.0, 1.0, 1.08),
        height_ratios=(1.0, 0.055),
        wspace=0.18,
    )
    balance_ax = fig.add_subplot(outer[0, 0])
    spatial_ax = fig.add_subplot(outer[0, 1])
    context_grid = outer[0, 2].subgridspec(2, 1, hspace=0.24)
    context_median_ax = fig.add_subplot(context_grid[0, 0])
    context_high_ax = fig.add_subplot(context_grid[1, 0])
    footer_ax = fig.add_subplot(outer[1, :])
    footer_ax.axis("off")
    all_axes = (balance_ax, spatial_ax, context_median_ax, context_high_ax)
    labels = [
        f"{row['pair_id']} [{row['comparability_status']}]"
        + (" [SITE SHIFT]" if row.get("site_comparability") == "different_site" else "")
        for _, row in passed.iterrows()
    ]
    positions = np.arange(len(passed), dtype=float)

    def setup_axis(ax, title: str, xlabel: str) -> None:
        ax.set_title(title, weight="bold")
        ax.axvline(0, color="#111827", linewidth=0.9)
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", alpha=0.2)

    if passed.empty:
        for ax in all_axes:
            ax.text(0.5, 0.5, "No completed pairs", ha="center", va="center")
            ax.axis("off")
    else:
        balance = pd.to_numeric(passed["delta_cs_balance"], errors="coerce").to_numpy(dtype=float)
        balance_bars = balance_ax.barh(positions, np.nan_to_num(balance, nan=0.0), color="#1d4ed8", alpha=0.86)
        balance_ax.set_yticks(positions, labels)
        balance_ax.invert_yaxis()
        setup_axis(balance_ax, "C/S balance", "Raw Delta R = Post - Pre")
        _center_zero_axis(balance_ax, balance.tolist())
        _label_bars(balance_ax, balance_bars, balance.tolist())

        spatial_specs = (
            ("delta_interface", "Interface", "#047857"),
            ("delta_diffuse", "Diffuse", "#0f766e"),
            ("delta_burden", "Burden", "#b45309"),
        )
        offsets = (-0.24, 0.0, 0.24)
        spatial_values: list[float] = []
        for (column, label, color), offset in zip(spatial_specs, offsets):
            values = pd.to_numeric(passed[column], errors="coerce").to_numpy(dtype=float)
            spatial_values.extend(values.tolist())
            bars = spatial_ax.barh(positions + offset, np.nan_to_num(values, nan=0.0), height=0.21, label=label, color=color, alpha=0.82)
            _label_bars(spatial_ax, bars, values.tolist())
        spatial_ax.set_yticks(positions, labels)
        spatial_ax.invert_yaxis()
        setup_axis(spatial_ax, "Spatial organization", "Raw Delta = Post - Pre")
        _center_zero_axis(spatial_ax, spatial_values)
        spatial_ax.legend(loc="best", fontsize=8)

        context_layers = _context_plot_layers(passed)

        def plot_context_layer(ax, values_by_axis, title: str, xlabel: str) -> None:
            context_specs = (
                (values_by_axis[0], "H context", "#7e22ce"),
                (values_by_axis[1], "V context", "#15803d"),
            )
            context_available = False
            context_values: list[float] = []
            for (values, label, color), offset in zip(context_specs, (-0.16, 0.16)):
                context_values.extend(values.tolist())
                context_available = context_available or np.isfinite(values).any()
                bars = ax.barh(
                    positions + offset,
                    np.nan_to_num(values, nan=0.0),
                    height=0.28,
                    label=label,
                    color=color,
                    alpha=0.82,
                )
                _label_bars(ax, bars, values.tolist())
            ax.set_yticks(positions, labels)
            ax.invert_yaxis()
            setup_axis(ax, title, xlabel)
            _center_zero_axis(ax, context_values)
            if context_available:
                ax.legend(loc="best", fontsize=8)
            else:
                ax.text(0.5, 0.5, "H/V not available", transform=ax.transAxes, ha="center", va="center")

        plot_context_layer(
            context_median_ax,
            context_layers["median"],
            "Parallel context axes — Median Δ",
            "Raw median Delta = Post - Pre",
        )
        plot_context_layer(
            context_high_ax,
            context_layers["pair_pooled_high_fraction"],
            "Pair-pooled high-context fraction Δ",
            "Raw fraction Delta = Post - Pre (shared within-pair q90)",
        )

    fig.suptitle("Multiaxial Change Profile — raw values, no composite response score", fontsize=14, weight="bold")
    footer_ax.text(
        0.0,
        0.5,
        f"SpatialTX Studio v{__version__} | H/V panels separately show raw-median Delta and pair-pooled high-context fraction Delta (shared within-pair q90). No composite response score; H/V do not modify C/S, R, masks, or Type A/B/C.",
        fontsize=7,
        color="#4b5563",
        ha="left",
        va="center",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output
