from __future__ import annotations

import html
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__
from .models import EXPLORATORY_NOTICE, NO_REGISTRATION_NOTICE, ComparativeConfig


def _mean_delta(delta_metrics: pd.DataFrame, name: str) -> float:
    values = pd.to_numeric(
        delta_metrics.loc[delta_metrics["delta_metric"].eq(name) & delta_metrics["status"].eq("ok"), "delta"],
        errors="coerce",
    ).dropna()
    return float(values.mean()) if len(values) else np.nan


def rules_based_summary(
    sample_metrics: pd.DataFrame,
    delta_metrics: pd.DataFrame,
    transitions: pd.DataFrame,
    config: ComparativeConfig,
    effective_mode: str,
    *,
    metric_change_table: pd.DataFrame | None = None,
    scale_warnings: pd.DataFrame | None = None,
    hv_summary: pd.DataFrame | None = None,
) -> str:
    clauses: list[str] = []
    metric_change_table = metric_change_table if metric_change_table is not None else pd.DataFrame()
    scale_warnings = scale_warnings if scale_warnings is not None else pd.DataFrame()
    hv_summary = hv_summary if hv_summary is not None else pd.DataFrame()

    def enriched_delta(name: str) -> float:
        if metric_change_table.empty:
            return np.nan
        values = pd.to_numeric(
            metric_change_table.loc[metric_change_table["metric_name"].eq(name), "raw_delta"], errors="coerce"
        ).dropna()
        return float(values.mean()) if len(values) else np.nan

    normalized_diffuse = enriched_delta("diffuse_components_per_1000_valid_spots")
    diffuse = _mean_delta(delta_metrics, "delta_diffuse_fraction")
    localized = _mean_delta(delta_metrics, "delta_localized_interface_fraction")
    burden = _mean_delta(delta_metrics, "delta_transition_burden_score")
    if np.isfinite(normalized_diffuse):
        clauses.append(
            f"a {'higher' if normalized_diffuse > 0 else 'lower' if normalized_diffuse < 0 else 'similar'} "
            "normalized diffuse-component density"
        )
    if np.isfinite(diffuse):
        clauses.append(f"a {'higher' if diffuse > 0 else 'lower' if diffuse < 0 else 'similar'} diffuse transition fraction")
    if np.isfinite(localized):
        clauses.append(
            f"a {'higher' if localized > 0 else 'lower' if localized < 0 else 'similar'} localized interface-like fraction"
        )
    if np.isfinite(burden):
        clauses.append(f"a {'higher' if burden > 0 else 'lower' if burden < 0 else 'similar'} transition burden score")
    comparison = ", ".join(clauses) if clauses else "no uniformly available primary transition metric change"
    text = (
        f"Compared with the reference definition ({config.reference}), the target definition ({config.target}) showed {comparison}."
    )
    direct = transitions[
        transitions.get("comparison_basis", pd.Series(index=transitions.index, dtype=str)).eq(
            "matched_sample_operational_change"
        )
    ]
    if len(direct) == 1:
        text += f" The operational regime comparison was {direct.iloc[0]['regime_transition']}."
    elif len(direct) > 1:
        changed = int((direct["reference_regime"] != direct["target_regime"]).sum())
        text += f" Operational regime labels differed in {changed} of {len(direct)} matched comparisons."
    else:
        text += " Unpaired operational regime distributions were compared without inferring sample-level state transitions."
    uncertain = int(transitions.get("transition_confidence_flag", pd.Series(dtype=str)).eq("uncertain").sum())
    if uncertain:
        text += f" {uncertain} regime comparison(s) were flagged as uncertain because confidence or QC was limited."
    group_sizes = sample_metrics.groupby("group")["sample_id"].nunique().to_dict()
    small = [f"{group} n={count}" for group, count in group_sizes.items() if count < 3]
    if small:
        text += " Small sample-size warning: " + ", ".join(small) + "."
    cautions = (
        scale_warnings.loc[scale_warnings.get("severity", pd.Series(index=scale_warnings.index, dtype=str)).eq("caution"), "message"]
        if not scale_warnings.empty and "message" in scale_warnings
        else pd.Series(dtype=str)
    )
    if len(cautions):
        text += " " + " ".join(dict.fromkeys(str(message) for message in cautions if str(message).strip()))
        raw_diffuse = enriched_delta("n_diffuse_components")
        if np.isfinite(raw_diffuse):
            text += (
                f" The raw diffuse-component count changed by {raw_diffuse:.3g}; this raw count is reported only "
                "with the preceding scale warning."
            )
    if sample_metrics.get("H_expr_available", pd.Series(False, index=sample_metrics.index)).astype(bool).any() or sample_metrics.get(
        "V_expr_available", pd.Series(False, index=sample_metrics.index)
    ).astype(bool).any():
        valid_hv = (
            not hv_summary.empty
            and pd.to_numeric(
                hv_summary.loc[~hv_summary.get("metric", pd.Series(dtype=str)).isin(["H_expr_mean", "V_expr_mean"]), "raw_delta"],
                errors="coerce",
            ).notna().any()
        )
        text += (
            " H/V non-centered distribution and transition-enrichment summaries were available as observational context only;"
            if valid_hv
            else " H/V centered means were omitted because within-sample centering constrains them near zero;"
        )
        text += " H/V does not alter C/S/R, transition masks, or operational regimes."
    text += f" {NO_REGISTRATION_NOTICE}"
    text += " Statistical significance alone must not be interpreted as biological importance."
    text += f" {EXPLORATORY_NOTICE}"
    return text


def write_html_report(
    path: Path,
    summary_text: str,
    tables: dict[str, pd.DataFrame],
    figures: list[Path],
    config: ComparativeConfig,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    sections: list[str] = []
    for title, table in tables.items():
        sections.append(f"<h2>{html.escape(title)}</h2>")
        if table.empty:
            sections.append("<p>No rows available.</p>")
        else:
            sections.append(table.to_html(index=False, border=0, classes="dataframe", max_rows=500))
    figure_html = "".join(
        f'<figure><img src="comparative_figures/{html.escape(figure.name)}" alt="{html.escape(figure.stem)}">'
        f"<figcaption>{html.escape(figure.stem)}</figcaption></figure>"
        for figure in figures
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SpatialTX comparative summary</title>
<style>
body {{ font-family: Segoe UI, Arial, sans-serif; margin: 2rem; color: #172033; }}
.notice {{ border-left: 5px solid #b45309; background: #fff7ed; padding: 1rem; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.82rem; margin-bottom: 1.5rem; }}
th, td {{ border: 1px solid #dbe3ed; padding: 0.35rem; text-align: left; }}
th {{ background: #eef4f8; position: sticky; top: 0; }}
figure {{ margin: 1.5rem 0; }} img {{ max-width: 100%; height: auto; }}
</style>
</head>
<body>
<h1>SpatialTX Studio v{html.escape(__version__)} — Comparative Spatial Transition Analysis</h1>
<p><strong>Mode:</strong> {html.escape(effective_mode_label(config.mode))} &nbsp; <strong>Reference:</strong> {html.escape(config.reference)} &nbsp; <strong>Target:</strong> {html.escape(config.target)}</p>
<div class="notice"><p>{html.escape(summary_text)}</p><p>{html.escape(NO_REGISTRATION_NOTICE)}</p><p>{html.escape(EXPLORATORY_NOTICE)}</p></div>
<h2>Figures</h2>{figure_html}
{''.join(sections)}
</body></html>
"""
    path.write_text(document, encoding="utf-8")
    return path


def effective_mode_label(mode: str) -> str:
    return {
        "pairwise": "Pairwise",
        "paired": "Paired groups",
        "unpaired": "Unpaired groups",
        "manifest_batch": "Manifest batch",
    }.get(mode, mode)


def write_pdf_report(
    path: Path,
    summary_text: str,
    figures: list[Path],
    config: ComparativeConfig,
    sample_count: int,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.08, 0.94, f"SpatialTX Studio v{__version__}", fontsize=16, weight="bold")
        fig.text(0.08, 0.90, "Comparative Spatial Transition Analysis", fontsize=14)
        fig.text(0.08, 0.86, f"Mode: {effective_mode_label(config.mode)}", fontsize=10)
        fig.text(0.08, 0.83, f"Reference: {config.reference}    Target: {config.target}    Samples: {sample_count}", fontsize=10)
        wrapped = "\n".join(textwrap.wrap(summary_text, width=96))
        fig.text(0.08, 0.76, wrapped, fontsize=9, va="top", linespacing=1.4)
        fig.text(0.08, 0.12, NO_REGISTRATION_NOTICE, fontsize=9, color="#7c2d12")
        fig.text(0.08, 0.08, EXPLORATORY_NOTICE, fontsize=9, color="#7c2d12")
        plt.axis("off")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        for figure_path in figures:
            try:
                image = mpimg.imread(figure_path)
            except (OSError, ValueError):
                continue
            fig, ax = plt.subplots(figsize=(11.69, 8.27))
            ax.imshow(image)
            ax.axis("off")
            ax.set_title(figure_path.stem, fontsize=11)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    return path
