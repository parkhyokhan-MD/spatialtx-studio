from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd

from ..graph.metadata import json_safe
from .core import compute_axis_reliability, summarize_axis
from .dependence import bh_adjust
from .models import (
    RELIABILITY_SCHEMA_VERSION,
    AxisReliabilityResult,
    ReliabilityConfig,
    ReliabilityResult,
)


RELIABILITY_PAIR_METRICS = {
    "B": "balance_B",
    "A": "activity_A",
    "D": "direction_D",
    "CA_strength": "ca_strength",
    "CA_fraction": "ca_fraction",
}

DESCRIPTIVE_INFERENCE_WARNING = (
    "Descriptive spot-distribution comparison of unregistered slides. "
    "Not specimen-level inference and not evidence of treatment effect."
)


def _stable_seed_offset(*parts: str) -> int:
    digest = hashlib.sha256("\x1f".join(parts).encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little", signed=False)


def _finite(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return array[np.isfinite(array)]


def _median(values: np.ndarray) -> float:
    finite = _finite(values)
    return float(np.median(finite)) if len(finite) else np.nan


def _bootstrap_delta_ci(
    pre: np.ndarray,
    post: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    pre_values, post_values = _finite(pre), _finite(post)
    if len(pre_values) < 2 or len(post_values) < 2:
        return np.nan, np.nan
    deltas = np.empty(int(iterations), dtype=float)
    for index in range(int(iterations)):
        pre_sample = pre_values[rng.integers(0, len(pre_values), len(pre_values))]
        post_sample = post_values[rng.integers(0, len(post_values), len(post_values))]
        deltas[index] = np.median(post_sample) - np.median(pre_sample)
    return float(np.quantile(deltas, 0.025)), float(np.quantile(deltas, 0.975))


def _label_permutation_p(
    pre: np.ndarray,
    post: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> float:
    pre_values, post_values = _finite(pre), _finite(post)
    if len(pre_values) < 2 or len(post_values) < 2:
        return np.nan
    observed = abs(float(np.median(post_values) - np.median(pre_values)))
    pooled = np.concatenate((pre_values, post_values))
    n_pre = len(pre_values)
    exceed = 0
    for _ in range(int(iterations)):
        shuffled = pooled[rng.permutation(len(pooled))]
        candidate = abs(float(np.median(shuffled[n_pre:]) - np.median(shuffled[:n_pre])))
        exceed += int(candidate >= observed)
    return float((exceed + 1) / (int(iterations) + 1))


def _sample_validity(
    result: AxisReliabilityResult,
    config: ReliabilityConfig,
) -> tuple[str, str, int, float]:
    total = int(len(result.valid_input))
    count = int(result.valid_input.sum())
    fraction = float(count / total) if total else 0.0
    count_failed = count < int(config.minimum_valid_spots)
    fraction_failed = fraction < float(config.warning_valid_fraction)
    if count_failed and fraction_failed:
        status = "qc_fail_insufficient_valid_spots_and_fraction"
    elif count_failed:
        status = "qc_fail_insufficient_valid_spots"
    elif fraction_failed:
        status = "qc_fail_insufficient_valid_fraction"
    elif fraction < float(config.minimum_valid_fraction):
        status = "warning_low_valid_fraction"
    else:
        status = "valid"
    reason = (
        f"valid={count}/{total} ({fraction:.1%}); minimum_valid_spots="
        f"{int(config.minimum_valid_spots)}; warning_valid_fraction="
        f"{float(config.warning_valid_fraction):.1%}; minimum_valid_fraction="
        f"{float(config.minimum_valid_fraction):.1%}; status={status}"
    )
    return status, reason, count, fraction


def _shared_source(pre_value: str, post_value: str) -> str:
    return pre_value if pre_value == post_value else f"Pre={pre_value}; Post={post_value}"


def _legacy_balance_bh_fdr_compatibility(
    pair_label: str,
    pre: AxisReliabilityResult,
    post: AxisReliabilityResult,
    config: ReliabilityConfig,
) -> float:
    """Reproduce the v0.65-v1 mixed-metric B correction without exporting it as Activity.

    v1 corrected B together with A/D/CA p-values computed from the signed C/S
    arrays whenever the old absolute-count-only gate passed.  The corrected
    Activity contract changes those companion p-values and would therefore
    change the established ``delta_B_bh_fdr`` field despite identical B data.
    This isolated compatibility path preserves that one legacy field.  New
    source-correct FDR is exported separately as
    ``delta_B_reliability_v2_bh_fdr``.
    """

    legacy_pre = compute_axis_reliability(pre.C, pre.S, config, axis=pre.axis)
    legacy_post = compute_axis_reliability(post.C, post.S, config, axis=post.axis)
    old_pair_domain_valid = (
        int(legacy_pre.valid_input.sum()) >= int(config.minimum_valid_spots)
        and int(legacy_post.valid_input.sum()) >= int(config.minimum_valid_spots)
    )
    p_values: list[float] = []
    for metric_name, attribute in RELIABILITY_PAIR_METRICS.items():
        if metric_name != "B" and not old_pair_domain_valid:
            p_values.append(np.nan)
            continue
        offset = _stable_seed_offset(pair_label, pre.axis, metric_name)
        seed_sequence = np.random.SeedSequence([int(config.seed), offset])
        _bootstrap_seed, permutation_seed = seed_sequence.spawn(2)
        p_values.append(_label_permutation_p(
            getattr(legacy_pre, attribute),
            getattr(legacy_post, attribute),
            int(config.permutation_iterations),
            np.random.default_rng(permutation_seed),
        ))
    return float(bh_adjust(p_values)[0])


def build_pair_summary(
    pair_label: str,
    pre: AxisReliabilityResult,
    post: AxisReliabilityResult,
    config: ReliabilityConfig | Mapping | None = None,
    *,
    pre_coverage: Mapping[str, float] | None = None,
    post_coverage: Mapping[str, float] | None = None,
) -> pd.DataFrame:
    cfg = ReliabilityConfig.from_value(config)
    cfg.validate()
    if pre.axis != post.axis:
        raise ValueError("Pre and Post reliability results must describe the same axis.")
    pre_validity, pre_reason, pre_valid_count, pre_valid_fraction = _sample_validity(pre, cfg)
    post_validity, post_reason, post_valid_count, post_valid_fraction = _sample_validity(post, cfg)
    sample_validities = (pre_validity, post_validity)
    qc_fail = any(status.startswith("qc_fail_") for status in sample_validities)
    warning = any(status == "warning_low_valid_fraction" for status in sample_validities)
    if qc_fail:
        failed = [status for status in sample_validities if status.startswith("qc_fail_")]
        pair_validity = failed[0] if len(set(failed)) == 1 else "qc_fail_pre_or_post"
    elif warning:
        pair_validity = "warning_low_valid_fraction"
    else:
        pair_validity = "valid"
    pair_reason = f"Pre: {pre_reason} | Post: {post_reason}"
    row: dict = {
        "pair_label": pair_label,
        "axis": pre.axis,
        "score_domain": cfg.score_domain,
        "balance_score_source": _shared_source(
            pre.balance_score_source, post.balance_score_source
        ),
        "balance_score_domain": _shared_source(
            pre.balance_score_domain, post.balance_score_domain
        ),
        "activity_score_source": _shared_source(
            pre.activity_score_source, post.activity_score_source
        ),
        "activity_score_domain": _shared_source(
            pre.activity_score_domain, post.activity_score_domain
        ),
        "activity_source_transformations": _shared_source(
            pre.activity_source_transformations,
            post.activity_source_transformations,
        ),
        "activity_source_version": _shared_source(
            pre.activity_source_version, post.activity_source_version
        ),
        "epsilon": float(cfg.epsilon),
        "classification_enabled": bool(cfg.classification_enabled),
        "activity_threshold": cfg.activity_threshold,
        "direction_threshold": cfg.direction_threshold,
        "resampling_scope": "spot_distribution_descriptive_unregistered_slides",
        "specimen_level_inference": False,
        "inference_level": "spot_distribution_descriptive",
        "registered_spots": False,
        "biological_replicate_inference": False,
        "treatment_effect_claim_allowed": False,
        "inference_warning": DESCRIPTIVE_INFERENCE_WARNING,
        "minimum_valid_spots": int(cfg.minimum_valid_spots),
        "minimum_valid_fraction": float(cfg.minimum_valid_fraction),
        "warning_valid_fraction": float(cfg.warning_valid_fraction),
        "pre_total_spot_count": int(len(pre.valid_input)),
        "post_total_spot_count": int(len(post.valid_input)),
        "pre_valid_spot_count": pre_valid_count,
        "post_valid_spot_count": post_valid_count,
        "pre_valid_input_fraction": pre_valid_fraction,
        "post_valid_input_fraction": post_valid_fraction,
        "pre_validity_reason": pre_reason,
        "post_validity_reason": post_reason,
        "pair_validity_reason": pair_reason,
        "pre_direction_defined_spot_count": int(pre.direction_defined.sum()),
        "post_direction_defined_spot_count": int(post.direction_defined.sum()),
    }
    row["pre_score_validity"] = pre_validity
    row["post_score_validity"] = post_validity
    row["pair_score_validity"] = pair_validity
    row["balance_validity"] = "valid_for_finite_signed_inputs"
    row["activity_summary_included_in_conclusion"] = not qc_fail
    pre_activity_balance = _median(pre.activity_balance)
    post_activity_balance = _median(post.activity_balance)
    row["pre_activity_balance"] = pre_activity_balance
    row["post_activity_balance"] = post_activity_balance
    row["delta_activity_balance"] = (
        post_activity_balance - pre_activity_balance
        if np.isfinite(pre_activity_balance) and np.isfinite(post_activity_balance)
        else np.nan
    )
    p_columns: list[str] = []
    for metric_name, attribute in RELIABILITY_PAIR_METRICS.items():
        pre_values = getattr(pre, attribute)
        post_values = getattr(post, attribute)
        pre_median = _median(pre_values)
        post_median = _median(post_values)
        delta = post_median - pre_median if np.isfinite(pre_median) and np.isfinite(post_median) else np.nan
        offset = _stable_seed_offset(pair_label, pre.axis, metric_name)
        seed_sequence = np.random.SeedSequence([int(cfg.seed), offset])
        bootstrap_seed, permutation_seed = seed_sequence.spawn(2)
        ci_low, ci_high = _bootstrap_delta_ci(
            pre_values,
            post_values,
            int(cfg.bootstrap_iterations),
            np.random.default_rng(bootstrap_seed),
        )
        p_value = _label_permutation_p(
            pre_values,
            post_values,
            int(cfg.permutation_iterations),
            np.random.default_rng(permutation_seed),
        )
        row.update({
            f"pre_{metric_name}": pre_median,
            f"post_{metric_name}": post_median,
            f"delta_{metric_name}": delta,
            f"delta_{metric_name}_bootstrap_ci_low": ci_low,
            f"delta_{metric_name}_bootstrap_ci_high": ci_high,
            f"delta_{metric_name}_permutation_p_value": p_value,
            f"delta_{metric_name}_bh_fdr": np.nan,
        })
        p_columns.append(f"delta_{metric_name}_permutation_p_value")
    q_values = bh_adjust([row[column] for column in p_columns])
    for metric_name, q_value in zip(RELIABILITY_PAIR_METRICS, q_values):
        row[f"delta_{metric_name}_bh_fdr"] = q_value
    row["delta_B_reliability_v2_bh_fdr"] = row["delta_B_bh_fdr"]
    row["delta_B_bh_fdr"] = _legacy_balance_bh_fdr_compatibility(
        pair_label,
        pre,
        post,
        cfg,
    )
    row["delta_B_bh_fdr_scope"] = (
        "legacy_v0.65-v1_mixed_metric_compatibility; signed C/S used only to preserve this "
        "existing multiplicity field, never as the v2 Activity source"
    )
    row["delta_B_reliability_v2_bh_fdr_scope"] = (
        "v0.65-v2 source-correct B/A/D/CA descriptive metric family"
    )

    pre_summary, post_summary = summarize_axis(pre), summarize_axis(post)
    for name in (
        "valid_input_fraction",
        "direction_defined_fraction",
        "undefined_fraction",
        "negative_input_fraction",
        "nonfinite_input_fraction",
        "low_activity_fraction",
        "c_dominant_active_fraction",
        "s_dominant_active_fraction",
        "active_coactivation_candidate_fraction",
    ):
        pre_value, post_value = pre_summary[name], post_summary[name]
        row[f"pre_{name}"] = pre_value
        row[f"post_{name}"] = post_value
        row[f"delta_{name}"] = (
            post_value - pre_value
            if np.isfinite(pre_value) and np.isfinite(post_value)
            else np.nan
        )
    for pole in ("C", "S"):
        pre_value = float((pre_coverage or {}).get(pole, np.nan))
        post_value = float((post_coverage or {}).get(pole, np.nan))
        row[f"pre_{pole}_gene_coverage"] = pre_value
        row[f"post_{pole}_gene_coverage"] = post_value
        row[f"delta_{pole}_gene_coverage"] = (
            post_value - pre_value
            if np.isfinite(pre_value) and np.isfinite(post_value)
            else np.nan
        )
    return pd.DataFrame([row])


def _write_csv(path: Path, table: pd.DataFrame, columns: list[str] | None = None) -> Path:
    output = table.copy()
    if output.empty and columns is not None:
        output = pd.DataFrame(columns=columns)
    output.to_csv(path, index=False)
    return path


def plot_dependence_heatmap(matrix: pd.DataFrame, output: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    if matrix.empty or len(matrix) < 2:
        ax.axis("off")
        axis_text = ", ".join(map(str, matrix.index)) if len(matrix) else "none"
        ax.text(
            0.5,
            0.55,
            "Axis-dependence heatmap not applicable",
            ha="center",
            va="center",
            fontsize=15,
            weight="bold",
        )
        ax.text(
            0.5,
            0.43,
            f"Paired-pole reliability axes available: {axis_text}\nAt least two explicit paired-pole axes are required.",
            ha="center",
            va="center",
            fontsize=10,
            color="#475569",
        )
    else:
        image = ax.imshow(matrix.to_numpy(dtype=float), vmin=-1.0, vmax=1.0, cmap="coolwarm")
        ax.set_xticks(range(len(matrix.columns)), labels=matrix.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(matrix.index)), labels=matrix.index)
        ax.set_title("Direction dependence — Pearson correlation")
        fig.colorbar(image, ax=ax, label="correlation", shrink=0.85)
        for row in range(len(matrix.index)):
            for column in range(len(matrix.columns)):
                value = matrix.iat[row, column]
                label = "NA" if not np.isfinite(value) else f"{value:.2f}"
                ax.text(column, row, label, ha="center", va="center", fontsize=8)
    fig.text(
        0.5,
        0.015,
        "QC display only; axes are not orthogonalized, transformed, or removed.",
        ha="center",
        fontsize=8,
        color="#475569",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_valid_fraction_qc(
    pair_summary: pd.DataFrame,
    config: ReliabilityConfig,
    output: Path,
) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    if pair_summary.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "Reliability valid-fraction QC unavailable", ha="center", va="center")
    else:
        x = np.arange(len(pair_summary), dtype=float)
        width = 0.36
        pre = pd.to_numeric(pair_summary["pre_valid_input_fraction"], errors="coerce")
        post = pd.to_numeric(pair_summary["post_valid_input_fraction"], errors="coerce")
        ax.bar(x - width / 2, pre, width, label="Pre", color="#2563eb")
        ax.bar(x + width / 2, post, width, label="Post", color="#0d9488")
        ax.axhline(
            float(config.minimum_valid_fraction),
            color="#15803d",
            linestyle="--",
            label=f"valid threshold ({float(config.minimum_valid_fraction):.0%})",
        )
        ax.axhline(
            float(config.warning_valid_fraction),
            color="#b45309",
            linestyle=":",
            label=f"QC-fail threshold ({float(config.warning_valid_fraction):.0%})",
        )
        ax.set_xticks(x, pair_summary["pair_label"].astype(str), rotation=25, ha="right")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Valid nonnegative Activity input fraction")
        ax.set_title("v0.65 Reliability source-domain QC")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(axis="y", alpha=0.2)
    fig.text(
        0.5,
        0.01,
        DESCRIPTIVE_INFERENCE_WARNING,
        ha="center",
        fontsize=7,
        color="#475569",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def write_reliability_sidecars(
    run_dir: str | Path,
    result: ReliabilityResult,
    config: ReliabilityConfig | Mapping,
) -> list[Path]:
    cfg = ReliabilityConfig.from_value(config)
    cfg.validate()
    root = Path(run_dir)
    files = [
        _write_csv(root / "reliability_spot_results.csv", result.spot_results),
        _write_csv(root / "reliability_pair_summary.csv", result.pair_summary),
        _write_csv(root / "reliability_gene_coverage.csv", result.gene_coverage),
        _write_csv(root / "cross_exclusivity_audit.csv", result.cross_exclusivity_audit),
        _write_csv(root / "axis_dependence_long.csv", result.axis_dependence_long),
        _write_csv(
            root / "reliability_score_domain_diagnostic.csv",
            result.score_domain_diagnostic,
        ),
    ]
    diagnostic_json_path = root / "reliability_score_domain_diagnostic.json"
    diagnostic_json_path.write_text(
        json.dumps(
            json_safe({
                **dict(result.score_domain_diagnostic_metadata),
                "records": result.score_domain_diagnostic.to_dict(orient="records"),
            }),
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    files.append(diagnostic_json_path)
    matrix_path = root / "axis_dependence_matrix.csv"
    result.axis_dependence_matrix.to_csv(matrix_path, index=True, index_label="axis")
    files.append(matrix_path)
    heatmap_path = root / "axis_dependence_heatmap.png"
    files.append(plot_dependence_heatmap(result.axis_dependence_matrix, heatmap_path))
    valid_fraction_path = root / "reliability_valid_fraction_qc.png"
    files.append(plot_valid_fraction_qc(result.pair_summary, cfg, valid_fraction_path))

    qc_payload = {
        "reliability_schema_version": RELIABILITY_SCHEMA_VERSION,
        "configuration": cfg.to_dict(),
        "score_domain_interpretation": (
            "Legacy Balance uses preserved signed v0.6 C/S. Activity, Direction, and Co-activation "
            "use a separate pre-centering nonnegative program abundance. Neither source is clipped, "
            "shifted, min-max scaled, or overwritten."
        ),
        "classification_mode": "classified" if cfg.classification_enabled else "continuous",
        "state_fraction_denominator": "classified_valid_spots",
        "dependence_policy": "reported_only_no_axis_transformation",
        "resampling_scope": "spot_distribution_descriptive_unregistered_slides",
        "specimen_level_inference": False,
        "inference_level": "spot_distribution_descriptive",
        "registered_spots": False,
        "biological_replicate_inference": False,
        "treatment_effect_claim_allowed": False,
        "inference_warning": DESCRIPTIVE_INFERENCE_WARNING,
        "H_V_reinterpreted_as_paired_pole_axes": False,
        "axis_weights": {
            str(axis): 1.0 for axis in result.axis_dependence_matrix.index
        },
        "pan_cancer_weight": 1.0,
        "axis_count": int(len(result.axis_dependence_matrix)),
        "spot_row_count": int(len(result.spot_results)),
        "pair_summary_row_count": int(len(result.pair_summary)),
        "gene_coverage_row_count": int(len(result.gene_coverage)),
        "cross_exclusivity_hard_error_count": int(
            result.cross_exclusivity_audit.get("severity", pd.Series(dtype=str)).eq("hard_error").sum()
        ),
        "dependence_warning_count": int(
            result.axis_dependence_long.get("qc_status", pd.Series(dtype=str)).eq("warning_high_dependence").sum()
        ),
        **dict(result.qc),
        "output_files": [path.name for path in files] + ["reliability_qc.json"],
    }
    qc_path = root / "reliability_qc.json"
    qc_path.write_text(
        json.dumps(json_safe(qc_payload), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    files.append(qc_path)
    result.files = files
    return files
