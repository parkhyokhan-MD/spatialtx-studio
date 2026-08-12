from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__
from ..graph.metadata import json_safe
from .metric_registry import (
    METRIC_REGISTRY,
    heatmap_metrics,
    metric_definition,
    metrics_for_plot_group,
)
from .models import EXPLORATORY_NOTICE, NO_REGISTRATION_NOTICE, ComparativeConfig


PRIMARY_SPATIAL_STATE_METRICS = (
    "localized_interface_fraction",
    "diffuse_fraction",
    "transition_burden_score",
    "R_crossing_fraction",
    "largest_diffuse_component_ratio",
    "small_component_fraction",
)

TOPOLOGY_COMPONENT_COMPLEXITY_METRICS = (
    "diffuse_components_per_1000_valid_spots",
    "diffuse_components_per_tissue_component",
    "small_components_per_1000_valid_spots",
    "transition_components_per_1000_transition_spots",
    "normalized_fragmentation_score",
)


def _plot_modules():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _save(
    fig,
    path: Path,
    config: ComparativeConfig,
    figure_type: str,
    extra: dict | None = None,
    *,
    effective_mode: str | None = None,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = effective_mode or config.mode
    fig.text(
        0.01,
        0.005,
        (
            f"SpatialTX Studio v{__version__} | mode={mode} | {config.reference} -> {config.target} | "
            f"{NO_REGISTRATION_NOTICE} | Exploratory, non-diagnostic"
        ),
        fontsize=7,
        color="#4b5563",
    )
    fig.savefig(path, dpi=160, bbox_inches="tight")
    metadata = {
        "spatialtx_version": __version__,
        "figure_type": figure_type,
        "comparison_mode": config.mode,
        "effective_comparison_mode": mode,
        "reference": config.reference,
        "target": config.target,
        "delta_direction": "Target - Reference",
        "seed": config.seed,
        "thresholds": {"C": config.c_q, "S": config.s_q, "G": config.g_q},
        "C_gene_program": config.c_genes,
        "S_gene_program": config.s_genes,
        "scoring_options": config.scoring_options,
        "graph_settings": config.graph_settings,
        "H_V_observational_only": True,
        "spotwise_subtraction_performed": False,
        "registration_notice": NO_REGISTRATION_NOTICE,
        "interpretation_notice": EXPLORATORY_NOTICE,
        **(extra or {}),
    }
    path.with_suffix(path.suffix + ".metadata.json").write_text(
        json.dumps(json_safe(metadata), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    _plot_modules().close(fig)
    return path


def _format_value(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if abs(value) <= 1e-12:
        return "0"
    if abs(value) >= 1000 or abs(value) < 0.001:
        return f"{value:.3g}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _label_horizontal_bars(ax, bars, values: list[float]) -> None:
    finite = [abs(value) for value in values if np.isfinite(value)]
    pad = max(finite, default=1.0) * 0.025
    for bar, value in zip(bars, values):
        if not np.isfinite(value):
            continue
        inside = abs(value) > 4 * pad
        position = value - pad if value > 0 and inside else value + pad if value < 0 and inside else value + (pad if value >= 0 else -pad)
        alignment = "right" if value > 0 and inside else "left" if value < 0 and inside else "left" if value >= 0 else "right"
        ax.text(
            position,
            bar.get_y() + bar.get_height() / 2,
            _format_value(value),
            va="center",
            ha=alignment,
            fontsize=8,
        )


def _mean_change_rows(
    metric_changes: pd.DataFrame,
    definitions,
    *,
    value_column: str = "raw_delta",
) -> pd.DataFrame:
    rows: list[dict] = []
    for definition in definitions:
        selected = metric_changes.loc[
            metric_changes["metric_name"].eq(definition.internal_name)
            & metric_changes["status"].eq("ok")
        ]
        if value_column not in selected:
            continue
        values = pd.to_numeric(selected.get(value_column), errors="coerce").dropna()
        if not len(values):
            continue
        reference = pd.to_numeric(selected.get("reference_value"), errors="coerce").dropna()
        target = pd.to_numeric(selected.get("target_value"), errors="coerce").dropna()
        row = {
            "metric_name": definition.internal_name,
            "display_name": definition.display_name,
            "unit": definition.unit,
            value_column: float(values.mean()),
            "reference_value": float(reference.mean()) if len(reference) else np.nan,
            "target_value": float(target.mean()) if len(target) else np.nan,
        }
        if value_column != "raw_delta":
            raw_values = pd.to_numeric(selected.get("raw_delta"), errors="coerce").dropna()
            row["raw_delta"] = float(raw_values.mean()) if len(raw_values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _metric_definitions(names: tuple[str, ...]):
    definitions = {definition.internal_name: definition for definition in METRIC_REGISTRY}
    return tuple(definitions[name] for name in names)


def _plot_overview_panel(ax, data: pd.DataFrame, value_column: str, title: str, x_label: str) -> list[dict]:
    ax.set_title(title, fontsize=11, weight="bold")
    if data.empty:
        ax.text(0.5, 0.5, "No available metrics", ha="center", va="center")
        ax.axis("off")
        return []
    raw_values = data[value_column].astype(float).tolist()
    values = [0.0 if abs(value) <= 1e-12 else value for value in raw_values]
    positions = np.arange(len(data))
    bars = ax.barh(
        positions,
        values,
        color=["#b91c1c" if value < 0 else "#0369a1" for value in values],
        alpha=0.86,
    )
    ax.set_yticks(positions, data["display_name"])
    ax.invert_yaxis()
    ax.axvline(0, color="#111827", linewidth=0.9)
    ax.set_xlabel(x_label)
    ax.grid(axis="x", alpha=0.2)
    ax.margins(x=0.26)
    _label_horizontal_bars(ax, bars, values)
    return data.to_dict("records")


def _plot_two_panel_metric_overview(
    metric_changes: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    value_column: str,
    figure_type: str,
    overall_title: str,
    x_label: str,
    effective_mode: str | None = None,
) -> Path | None:
    panel_a = _mean_change_rows(
        metric_changes,
        _metric_definitions(PRIMARY_SPATIAL_STATE_METRICS),
        value_column=value_column,
    )
    panel_b = _mean_change_rows(
        metric_changes,
        _metric_definitions(TOPOLOGY_COMPONENT_COMPLEXITY_METRICS),
        value_column=value_column,
    )
    if value_column != "raw_delta" and panel_a.empty and panel_b.empty:
        return None

    plt = _plot_modules()
    fig, axes = plt.subplots(1, 2, figsize=(17, 7.2), sharex=False)
    displayed = {
        "panel_a": _plot_overview_panel(
            axes[0],
            panel_a,
            value_column,
            "Panel A — Primary spatial-state summary metrics",
            x_label,
        ),
        "panel_b": _plot_overview_panel(
            axes[1],
            panel_b,
            value_column,
            "Panel B — Topology / component complexity metrics",
            x_label,
        ),
    }
    fig.suptitle(overall_title, fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0.055, 1, 0.94), w_pad=4.0)
    return _save(
        fig,
        output,
        config,
        figure_type,
        {
            "display_mode": value_column,
            "panel_titles": {
                "panel_a": "Primary spatial-state summary metrics",
                "panel_b": "Topology / component complexity metrics",
            },
            "panel_metric_names": {
                "panel_a": list(PRIMARY_SPATIAL_STATE_METRICS),
                "panel_b": list(TOPOLOGY_COMPONENT_COMPLEXITY_METRICS),
            },
            "displayed_values": displayed,
            "x_axes_shared": False,
            "numeric_value_labels": True,
            "raw_values_preserved_in_csv_and_metadata": True,
            "tissue_denominator_note": (
                "The existing diffuse_components_per_tissue_component metric is retained; "
                "it is not labeled as physical tissue area."
            ),
        },
        effective_mode=effective_mode,
    )


def _plot_change_group(
    metric_changes: pd.DataFrame,
    plot_group: str,
    output: Path,
    config: ComparativeConfig,
    title: str,
    figure_type: str,
    *,
    note: str = "",
    effective_mode: str | None = None,
) -> Path:
    plt = _plot_modules()
    definitions = metrics_for_plot_group(plot_group)
    data = _mean_change_rows(metric_changes, definitions)
    units = list(dict.fromkeys(data["unit"].astype(str))) if not data.empty else ["value"]
    fig, axes = plt.subplots(
        1,
        len(units),
        figsize=(max(7.5, 6.0 * len(units)), max(4.2, 0.52 * max(len(data), 1))),
        squeeze=False,
    )
    metadata_rows: list[dict] = []
    for ax, unit in zip(axes[0], units):
        subset = data.loc[data["unit"].eq(unit)] if not data.empty else data
        if subset.empty:
            ax.text(0.5, 0.5, "No available metrics", ha="center", va="center")
            ax.axis("off")
            continue
        raw_values = subset["raw_delta"].astype(float).tolist()
        values = [0.0 if abs(value) <= 1e-12 else value for value in raw_values]
        positions = np.arange(len(subset))
        colors = ["#b91c1c" if value < 0 else "#0369a1" for value in values]
        bars = ax.barh(positions, values, color=colors, alpha=0.86)
        ax.axvline(0, color="#111827", linewidth=0.9)
        ax.set_yticks(positions, subset["display_name"])
        ax.invert_yaxis()
        ax.set_xlabel(f"Delta = Target - Reference ({unit})")
        ax.set_title(unit)
        ax.grid(axis="x", alpha=0.2)
        ax.margins(x=0.22)
        _label_horizontal_bars(ax, bars, values)
        metadata_rows.extend(subset.to_dict("records"))
    subtitle = title + (f"\n{note}" if note else "")
    fig.suptitle(subtitle, fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.91 if note else 0.94))
    return _save(
        fig,
        output,
        config,
        figure_type,
        {
            "displayed_values": metadata_rows,
            "undefined_metrics_omitted": True,
            "numeric_value_labels": True,
            "display_zero_tolerance": 1e-12,
            "raw_values_preserved_in_csv_and_metadata": True,
        },
        effective_mode=effective_mode,
    )


def plot_metric_change(
    metric_changes: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str | None = None,
) -> Path:
    """Raw-delta overview with independent primary and topology axes."""
    result = _plot_two_panel_metric_overview(
        metric_changes,
        output,
        config,
        value_column="raw_delta",
        figure_type="compatibility_primary_metric_overview",
        overall_title="Comparative spatial change overview — raw delta (Target - Reference)",
        x_label="Raw delta = Target - Reference",
        effective_mode=effective_mode,
    )
    assert result is not None
    return result


def plot_standardized_metric_change(
    metric_changes: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str | None = None,
) -> Path | None:
    """Optional group-only overview using pooled-sample-scale standardized deltas."""
    return _plot_two_panel_metric_overview(
        metric_changes,
        output,
        config,
        value_column="standardized_delta",
        figure_type="standardized_two_panel_metric_overview",
        overall_title="Comparative spatial change overview — standardized change",
        x_label="Standardized delta (pooled sample SD)",
        effective_mode=effective_mode,
    )


def plot_group_distributions(sample_metrics: pd.DataFrame, output: Path, config: ComparativeConfig, *, effective_mode: str | None = None) -> Path:
    plt = _plot_modules()
    metrics = [
        metric for metric in ("R_mean", "R_std", "localized_interface_fraction", "diffuse_fraction", "transition_burden_score")
        if metric in sample_metrics and pd.to_numeric(sample_metrics[metric], errors="coerce").notna().any()
    ]
    if not metrics:
        metrics = ["n_spots"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(max(9, 3.2 * len(metrics)), 5), squeeze=False)
    groups = [config.reference, config.target]
    for ax, metric in zip(axes[0], metrics):
        data = [
            pd.to_numeric(sample_metrics.loc[sample_metrics["group"].eq(group), metric], errors="coerce").dropna().to_numpy()
            for group in groups
        ]
        if all(len(values) for values in data):
            parts = ax.violinplot(data, positions=[1, 2], showmeans=True, showmedians=True)
            for body, color in zip(parts["bodies"], ("#64748b", "#0284c7")):
                body.set_facecolor(color)
                body.set_alpha(0.55)
            for index, values in enumerate(data, 1):
                ax.scatter(np.full(len(values), index), values, s=18, color="#111827", alpha=0.7)
        ax.set_xticks([1, 2], groups, rotation=20)
        ax.set_title(metric_definition(metric).display_name)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Sample-level group distributions")
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    return _save(fig, output, config, "group_distribution", effective_mode=effective_mode)


def plot_regime_transitions(
    transitions: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str,
) -> Path:
    plt = _plot_modules()
    fig, ax = plt.subplots(figsize=(8, 5.8))
    visualization_mode = "unavailable"
    if transitions.empty:
        ax.text(0.5, 0.5, "No regime transition data", ha="center", va="center")
        ax.axis("off")
    elif effective_mode == "pairwise" and len(transitions) == 1:
        visualization_mode = "pairwise_transition_card"
        row = transitions.iloc[0]
        reference_regime = str(row.get("reference_regime", "Unavailable"))
        target_regime = str(row.get("target_regime", "Unavailable"))
        confidence_flag = str(row.get("transition_confidence_flag", "uncertain"))
        if confidence_flag == "uncertain":
            status = "Uncertain operational comparison"
            color = "#b45309"
        elif reference_regime == target_regime:
            status = "Stable operational regime"
            color = "#047857"
        else:
            status = "Changed operational regime label"
            color = "#1d4ed8"
        ax.text(0.18, 0.62, reference_regime, ha="center", va="center", fontsize=13,
                bbox={"boxstyle": "round,pad=0.6", "facecolor": "#e2e8f0", "edgecolor": "#64748b"})
        ax.annotate("", xy=(0.72, 0.62), xytext=(0.32, 0.62), arrowprops={"arrowstyle": "->", "lw": 2, "color": "#334155"})
        ax.text(0.82, 0.62, target_regime, ha="center", va="center", fontsize=13,
                bbox={"boxstyle": "round,pad=0.6", "facecolor": "#dbeafe", "edgecolor": "#0284c7"})
        ax.text(0.18, 0.45, f"Confidence: {_format_value(float(row.get('reference_confidence', np.nan)))}", ha="center")
        ax.text(0.82, 0.45, f"Confidence: {_format_value(float(row.get('target_confidence', np.nan)))}", ha="center")
        ax.text(0.5, 0.27, status, ha="center", fontsize=12, weight="bold", color=color)
        ax.text(0.5, 0.15, "Descriptive operational classification only; not a biological state transition.", ha="center", fontsize=9)
        ax.set_title("Pairwise operational regime comparison")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
    elif transitions["comparison_basis"].eq("unpaired_group_distribution").all():
        visualization_mode = "unpaired_regime_distribution"
        labels = list(transitions["reference_regime"])
        x = np.arange(len(labels))
        ax.bar(x - 0.18, transitions["reference_fraction"], width=0.36, label=config.reference, color="#64748b")
        ax.bar(x + 0.18, transitions["target_fraction"], width=0.36, label=config.target, color="#0284c7")
        ax.set_xticks(x, labels, rotation=20, ha="right")
        ax.set_ylabel("Sample fraction")
        ax.legend()
        ax.set_title("Operational regime distribution comparison")
    else:
        visualization_mode = "matched_group_3x3_transition_matrix"
        order = ["Type_A_candidate", "Type_B_candidate", "Type_C_candidate"]
        matrix = pd.crosstab(transitions["reference_regime"], transitions["target_regime"]).reindex(
            index=order, columns=order, fill_value=0
        )
        row_totals = matrix.sum(axis=1).replace(0, np.nan)
        row_percent = matrix.div(row_totals, axis=0) * 100.0
        image = ax.imshow(matrix.to_numpy(dtype=float), cmap="Blues")
        for row_index in range(len(order)):
            for column_index in range(len(order)):
                count = int(matrix.iloc[row_index, column_index])
                percentage = row_percent.iloc[row_index, column_index]
                label = f"{count}\n{percentage:.1f}%" if np.isfinite(percentage) else f"{count}\nNA"
                ax.text(column_index, row_index, label, ha="center", va="center")
        uncertain = int(transitions.get("transition_confidence_flag", pd.Series(dtype=str)).eq("uncertain").sum())
        missing = int(
            (~transitions["reference_regime"].isin(order) | ~transitions["target_regime"].isin(order)).sum()
        )
        ax.set_xticks(range(len(order)), order, rotation=20, ha="right")
        ax.set_yticks(range(len(order)), order)
        ax.set_xlabel("Target operational regime")
        ax.set_ylabel("Reference operational regime")
        ax.set_title(f"Matched operational regime matrix: n={len(transitions)}, uncertain={uncertain}, missing={missing}\nCells show count and reference-row percentage")
        fig.colorbar(image, ax=ax, fraction=0.045, label="Matched-pair count")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return _save(
        fig,
        output,
        config,
        "operational_regime_transition",
        {"visualization_mode": visualization_mode},
        effective_mode=effective_mode,
    )


def plot_standardized_heatmap(
    sample_metrics: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str | None = None,
) -> Path:
    plt = _plot_modules()
    definitions = [definition for definition in heatmap_metrics() if definition.internal_name in sample_metrics]
    metrics = [definition.internal_name for definition in definitions]
    numeric = sample_metrics[metrics].apply(pd.to_numeric, errors="coerce") if metrics else pd.DataFrame(index=sample_metrics.index)
    available = [metric for metric in metrics if numeric[metric].notna().any()]
    definitions = [definition for definition in definitions if definition.internal_name in available]
    numeric = numeric[available]
    if numeric.empty:
        standard = pd.DataFrame(np.zeros((len(sample_metrics), 1)), columns=["No available metrics"])
        labels = ["No available metrics"]
        categories = ["Unavailable"]
    else:
        filled = numeric.fillna(numeric.mean())
        standard = ((filled - filled.mean()) / filled.std(ddof=0).replace(0, 1)).fillna(0)
        labels = [definition.display_name for definition in definitions]
        categories = [definition.category for definition in definitions]
    fig, ax = plt.subplots(figsize=(max(10, len(labels) * 0.56), max(4, len(sample_metrics) * 0.42 + 1.0)))
    image = ax.imshow(standard.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(labels)), labels, rotation=55, ha="right")
    ax.set_yticks(range(len(sample_metrics)), sample_metrics["sample_id"].astype(str))
    boundaries = [index for index in range(1, len(categories)) if categories[index] != categories[index - 1]]
    for boundary in boundaries:
        ax.axvline(boundary - 0.5, color="#111827", linewidth=1.3)
    starts = [0, *boundaries]
    ends = [boundary - 1 for boundary in boundaries] + [len(categories) - 1]
    for start, end in zip(starts, ends):
        ax.text((start + end) / 2, -0.85, categories[start], ha="center", va="bottom", fontsize=8, weight="bold")
    ax.set_title("Category-grouped standardized display heatmap")
    fig.colorbar(image, ax=ax, label="Within-run metric-wise z-score (visualization only)", fraction=0.025)
    fig.subplots_adjust(left=0.10, right=0.92, bottom=0.38, top=0.86)
    return _save(
        fig,
        output,
        config,
        "standardized_sample_metric_heatmap",
        {
            "visualization_transform": "metric-wise within-run z-score; missing values filled with metric mean",
            "raw_export_values_replaced": False,
            "metric_order": available,
            "metric_categories": categories,
        },
        effective_mode=effective_mode,
    )


def plot_side_by_side_maps(
    fields_by_sample: dict[str, dict[str, np.ndarray]],
    reference_sample: str,
    target_sample: str,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str | None = None,
) -> Path:
    from matplotlib.lines import Line2D

    plt = _plot_modules()
    ref = fields_by_sample[reference_sample]
    tar = fields_by_sample[target_sample]
    finite = np.r_[np.asarray(ref["R"], dtype=float), np.asarray(tar["R"], dtype=float)]
    limit = float(np.nanmax(np.abs(finite))) if np.isfinite(finite).any() else 1.0
    limit = max(limit, np.finfo(float).eps)
    # Reserve an independent footer band for the legend and the common run notice.
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.0))
    for ax, sample, fields in ((axes[0], reference_sample, ref), (axes[1], target_sample, tar)):
        coords = np.asarray(fields["coords"], dtype=float)
        image = ax.scatter(
            coords[:, 0], coords[:, 1], c=np.asarray(fields["R"], dtype=float), s=12,
            cmap="coolwarm", vmin=-limit, vmax=limit, linewidths=0,
        )
        mask = np.asarray(fields["interface"], dtype=bool) | np.asarray(fields["diffuse"], dtype=bool)
        if mask.any():
            ax.scatter(coords[mask, 0], coords[mask, 1], facecolors="none", edgecolors="#111827", s=22, linewidths=0.45)
        ax.set_title(f"{sample}\nR(x)=C(x)-S(x)")
        ax.set_aspect("equal", adjustable="datalim")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03, label="Fill color: R = C - S")
    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#111827",
               markersize=7, label="Outline: transition candidate (localized interface-like or diffuse mask)"),
    ]
    fig.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.065),
        frameon=False,
        fontsize=9,
    )
    fig.suptitle("Side-by-side selected sample R maps (no registration or spot-wise subtraction)")
    fig.tight_layout(rect=(0, 0.16, 1, 0.93))
    return _save(
        fig,
        output,
        config,
        "side_by_side_R_maps",
        {
            "reference_sample": reference_sample,
            "target_sample": target_sample,
            "fill_definition": "R = C - S",
            "outline_definition": "union of existing localized interface-like and diffuse transition masks",
            "spotwise_subtraction": False,
        },
        effective_mode=effective_mode,
    )


def plot_sample_scale(
    metric_changes: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str | None = None,
) -> Path:
    definitions = tuple(
        definition for definition in metrics_for_plot_group("sample_scale")
        if definition.internal_name in {
            "n_valid_spots", "n_in_tissue_spots", "tissue_component_count", "tissue_area_proxy",
            "spatial_extent_area_proxy",
        }
    )
    data = _mean_change_rows(metric_changes, definitions)
    plt = _plot_modules()
    units = list(dict.fromkeys(data["unit"].astype(str))) if not data.empty else ["value"]
    fig, axes = plt.subplots(1, len(units), figsize=(max(8, 5.2 * len(units)), 5), squeeze=False)
    for ax, unit in zip(axes[0], units):
        subset = data.loc[data["unit"].eq(unit)] if not data.empty else data
        if subset.empty:
            ax.text(0.5, 0.5, "No available sample-scale metrics", ha="center", va="center")
            ax.axis("off")
            continue
        x = np.arange(len(subset))
        ref = subset["reference_value"].astype(float).to_numpy()
        tar = subset["target_value"].astype(float).to_numpy()
        left = ax.bar(x - 0.19, ref, width=0.38, label=config.reference, color="#64748b")
        right = ax.bar(x + 0.19, tar, width=0.38, label=config.target, color="#0284c7")
        ax.set_xticks(x, subset["display_name"], rotation=25, ha="right")
        ax.set_ylabel(unit)
        ax.grid(axis="y", alpha=0.2)
        ax.bar_label(left, labels=[_format_value(value) for value in ref], fontsize=7, padding=2)
        ax.bar_label(right, labels=[_format_value(value) for value in tar], fontsize=7, padding=2)
        ax.legend()
    fig.suptitle("Reference and target sample-scale context")
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    return _save(fig, output, config, "sample_scale_context", effective_mode=effective_mode)


def plot_relative_changes(
    metric_changes: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str | None = None,
) -> Path:
    ordered = []
    for definition in METRIC_REGISTRY:
        if definition.deprecated or definition.observational_only or definition.scale_sensitive or definition.interpretation_priority > 3:
            continue
        selected = metric_changes.loc[metric_changes["metric_name"].eq(definition.internal_name)]
        values = pd.to_numeric(selected.get("symmetric_percent_change"), errors="coerce").dropna()
        if len(values):
            ordered.append((definition, float(values.mean())))
    plt = _plot_modules()
    fig, ax = plt.subplots(figsize=(10, max(4.5, len(ordered) * 0.34)))
    if not ordered:
        ax.text(0.5, 0.5, "No stable eligible relative changes", ha="center", va="center")
        ax.axis("off")
    else:
        values = [value for _definition, value in ordered]
        labels = [definition.display_name for definition, _value in ordered]
        bars = ax.barh(np.arange(len(values)), values, color=["#b91c1c" if value < 0 else "#0369a1" for value in values])
        ax.set_yticks(np.arange(len(values)), labels)
        ax.invert_yaxis()
        ax.axvline(0, color="#111827", linewidth=0.9)
        ax.set_xlabel("Symmetric percent change (%)")
        ax.set_xlim(-210, 210)
        ax.grid(axis="x", alpha=0.2)
        _label_horizontal_bars(ax, bars, values)
    ax.set_title("Relative changes (undefined and unstable ordinary percent changes omitted)")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return _save(
        fig,
        output,
        config,
        "symmetric_percent_changes",
        {"transform": "200 * (target-reference) / (abs(target)+abs(reference))", "undefined_metrics_omitted": True},
        effective_mode=effective_mode,
    )


def plot_hv_summary(
    hv_summary: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str | None = None,
) -> Path:
    plt = _plot_modules()
    panels = [
        ("raw_median", "Raw-scale median"),
        ("q90", "90th percentile"),
        ("high_fraction", "Pooled-threshold high fraction"),
        ("transition_enrichment", "Transition enrichment (median difference)"),
        ("spatial_variance", "Spatial variance"),
    ]
    fig, axes = plt.subplots(1, len(panels), figsize=(17, 4.8), squeeze=False)
    valid_panel_count = 0
    for ax, (suffix, title) in zip(axes[0], panels):
        values: list[float] = []
        labels: list[str] = []
        colors: list[str] = []
        for axis in ("H", "V"):
            name = f"{axis}_{suffix}"
            selected = hv_summary.loc[hv_summary["metric"].eq(name)]
            if selected.empty:
                continue
            ref = pd.to_numeric(selected["reference_value"], errors="coerce").dropna()
            tar = pd.to_numeric(selected["target_value"], errors="coerce").dropna()
            if len(ref):
                values.append(float(ref.mean())); labels.append(f"{axis}\n{config.reference}"); colors.append("#64748b")
            if len(tar):
                values.append(float(tar.mean())); labels.append(f"{axis}\n{config.target}"); colors.append("#0284c7")
        if not values:
            ax.text(0.5, 0.5, "Unavailable", ha="center", va="center")
            ax.axis("off")
            continue
        valid_panel_count += 1
        bars = ax.bar(np.arange(len(values)), values, color=colors)
        ax.set_xticks(np.arange(len(values)), labels)
        ax.axhline(0, color="#111827", linewidth=0.8)
        ax.grid(axis="y", alpha=0.2)
        ax.bar_label(bars, labels=[_format_value(value) for value in values], fontsize=7, padding=2)
        ax.set_title(title, fontsize=9)
    fig.suptitle(
        "H/V observational context summary\nPooled reference/target thresholds; H/V do not influence transition detection",
        fontsize=12,
    )
    if valid_panel_count == 0:
        fig.text(0.5, 0.16, "Non-centered H/V scores were unavailable; no values were fabricated.", ha="center", color="#b45309")
    fig.tight_layout(rect=(0, 0.06, 1, 0.89))
    return _save(
        fig,
        output,
        config,
        "observational_H_V_summary",
        {
            "observational_only": True,
            "centered_means_excluded": True,
            "high_fraction_threshold_method": "pooled_reference_target_q90",
            "transition_enrichment_definition": "median inside transition mask - median outside transition mask",
        },
        effective_mode=effective_mode,
    )


def plot_summary_card(
    metric_changes: pd.DataFrame,
    transitions: pd.DataFrame,
    scale_warnings: pd.DataFrame,
    hv_summary: pd.DataFrame,
    output: Path,
    config: ComparativeConfig,
    *,
    effective_mode: str,
) -> Path:
    import textwrap

    plt = _plot_modules()
    fig, ax = plt.subplots(figsize=(10, 6.4))
    ax.axis("off")

    def delta(metric: str) -> str:
        values = pd.to_numeric(metric_changes.loc[metric_changes["metric_name"].eq(metric), "raw_delta"], errors="coerce").dropna()
        return _format_value(float(values.mean())) if len(values) else "NA"

    transition = "Unavailable"
    if len(transitions):
        if effective_mode == "pairwise" and len(transitions) == 1:
            transition = str(transitions.iloc[0].get("regime_transition", "Unavailable"))
        else:
            transition = f"{len(transitions)} matched/distribution rows; see regime figure"
    scale_message = (
        str(scale_warnings.loc[scale_warnings["severity"].eq("caution"), "message"].iloc[0])
        if len(scale_warnings) and scale_warnings["severity"].eq("caution").any()
        else "No configured substantial sample-scale difference was detected."
    )
    hv_available = bool(
        len(hv_summary)
        and pd.to_numeric(hv_summary.loc[~hv_summary["metric"].isin(["H_expr_mean", "V_expr_mean"]), "raw_delta"], errors="coerce").notna().any()
    )
    lines = [
        ("Comparison", f"{config.reference} -> {config.target} ({effective_mode})"),
        ("Transition burden change", delta("transition_burden_score")),
        ("Localized interface-like fraction change", delta("localized_interface_fraction")),
        ("Diffuse fraction change", delta("diffuse_fraction")),
        ("Normalized diffuse component density change", delta("diffuse_components_per_1000_valid_spots")),
        ("Operational regime comparison", transition),
        ("H/V observational context", "Available; see H/V summary" if hv_available else "Non-centered summaries unavailable; no values fabricated"),
        ("Scale context", scale_message),
    ]
    ax.text(0.04, 0.94, "SpatialTX comparative summary card", fontsize=17, weight="bold", va="top")
    y = 0.84
    for label, value in lines:
        ax.text(0.05, y, label, fontsize=10, weight="bold", va="top", color="#334155")
        wrapped_value = "\n".join(textwrap.wrap(value, width=68)) or value
        ax.text(0.50, y, wrapped_value, fontsize=10, va="top")
        y -= 0.085 + 0.04 * max(0, wrapped_value.count("\n"))
    ax.text(0.05, 0.07, "Exploratory, non-diagnostic. No direct spatial registration or spot-wise subtraction.", fontsize=10, color="#9a3412", weight="bold")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return _save(fig, output, config, "pairwise_or_group_summary_card", effective_mode=effective_mode)


def _copy_compatibility_figure(source: Path, target: Path, figure_type: str) -> Path:
    shutil.copy2(source, target)
    source_metadata = source.with_suffix(source.suffix + ".metadata.json")
    payload = json.loads(source_metadata.read_text(encoding="utf-8"))
    payload["figure_type"] = figure_type
    payload["compatibility_alias_of"] = source.name
    target.with_suffix(target.suffix + ".metadata.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return target


def generate_comparative_figures(
    sample_metrics: pd.DataFrame,
    delta_metrics: pd.DataFrame,
    transitions: pd.DataFrame,
    fields_by_sample: dict[str, dict[str, np.ndarray]],
    reference_sample: str,
    target_sample: str,
    output_dir: Path,
    config: ComparativeConfig,
    *,
    metric_change_table: pd.DataFrame,
    sample_scale: pd.DataFrame,
    scale_warnings: pd.DataFrame,
    hv_summary: pd.DataFrame,
    effective_mode: str,
) -> list[Path]:
    del delta_metrics, sample_scale  # Raw tables remain exported; redesigned plots use the enriched change table.
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = [
        plot_summary_card(
            metric_change_table, transitions, scale_warnings, hv_summary,
            output_dir / "comparative_summary_card.png", config, effective_mode=effective_mode,
        ),
        _plot_change_group(
            metric_change_table, "program", output_dir / "comparative_program_score_changes.png", config,
            "Program score changes", "program_score_changes", effective_mode=effective_mode,
        ),
        _plot_change_group(
            metric_change_table, "transition", output_dir / "comparative_transition_changes.png", config,
            "Gradient and transition changes", "transition_metric_changes", effective_mode=effective_mode,
        ),
        _plot_change_group(
            metric_change_table, "graph", output_dir / "comparative_graph_changes.png", config,
            "Graph and adjacency changes", "graph_metric_changes", effective_mode=effective_mode,
        ),
        _plot_change_group(
            metric_change_table, "topology_raw", output_dir / "comparative_topology_raw_counts.png", config,
            "Topology raw-count changes", "topology_raw_count_changes",
            note="Raw counts; sensitive to tissue size, valid spot count, and tissue fragmentation.",
            effective_mode=effective_mode,
        ),
        _plot_change_group(
            metric_change_table, "topology_normalized", output_dir / "comparative_topology_normalized.png", config,
            "Normalized topology changes", "topology_normalized_changes", effective_mode=effective_mode,
        ),
        plot_sample_scale(metric_change_table, output_dir / "comparative_sample_scale.png", config, effective_mode=effective_mode),
        plot_relative_changes(metric_change_table, output_dir / "comparative_relative_changes.png", config, effective_mode=effective_mode),
        plot_standardized_heatmap(sample_metrics, output_dir / "comparative_standardized_heatmap.png", config, effective_mode=effective_mode),
        plot_regime_transitions(transitions, output_dir / "comparative_regime_transitions.png", config, effective_mode=effective_mode),
        plot_side_by_side_maps(
            fields_by_sample, reference_sample, target_sample,
            output_dir / "comparative_side_by_side_maps.png", config, effective_mode=effective_mode,
        ),
        plot_group_distributions(sample_metrics, output_dir / "comparative_group_distributions.png", config, effective_mode=effective_mode),
        plot_hv_summary(hv_summary, output_dir / "comparative_HV_summary.png", config, effective_mode=effective_mode),
    ]
    figures.append(plot_metric_change(
        metric_change_table, output_dir / "comparative_metric_changes.png", config, effective_mode=effective_mode
    ))
    standardized_overview = plot_standardized_metric_change(
        metric_change_table,
        output_dir / "comparative_metric_changes_standardized.png",
        config,
        effective_mode=effective_mode,
    )
    if standardized_overview is not None:
        figures.append(standardized_overview)
    figures.append(_copy_compatibility_figure(
        output_dir / "comparative_HV_summary.png",
        output_dir / "comparative_H_V_context.png",
        "compatibility_observational_H_V_summary",
    ))
    return figures
