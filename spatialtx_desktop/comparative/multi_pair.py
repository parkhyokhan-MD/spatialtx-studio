from __future__ import annotations

import datetime as dt
import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from threading import Event
from typing import Callable

import numpy as np
import pandas as pd

from .. import __version__
from ..graph.context import resolve_context_program
from ..graph.metadata import json_safe
from ..workflow import _read_h5ad
from .comparative_normalization import POOLED_HIGH_QUANTILE, apply_pooled_hv_thresholds, percent_change
from .metrics import COMPARATIVE_METRIC_LAYER_SCHEMA, analyze_sample, file_sha256
from .models import ComparativeConfig, SampleRecord
from .multiaxial import (
    SITE_COMPARABILITY_VALUES,
    build_comparative_qc_summary,
    build_multiaxial_pair_summary,
    normalize_site_comparability,
    plot_multiaxial_overview,
    qc_aware_interpretation,
    site_warning_text,
)


MULTI_PAIR_NOTICE = (
    "Observed Pre/Post differences are descriptive spatial changes. They do not establish treatment response, "
    "therapeutic efficacy, drug sensitivity, or clinical benefit."
)
MAX_MULTI_PAIRS = 6


@dataclass(frozen=True)
class PairSpec:
    label: str
    pre_path: Path
    post_path: Path
    site_comparability: str = "unknown_site"

    def normalized(self, fallback_index: int = 1) -> "PairSpec":
        return PairSpec(
            self.label.strip() or f"Pair_{fallback_index}",
            Path(self.pre_path).expanduser().resolve(),
            Path(self.post_path).expanduser().resolve(),
            normalize_site_comparability(self.site_comparability),
        )


@dataclass(frozen=True)
class MultiPairMetric:
    internal_name: str
    export_name: str
    display_name: str
    unit: str
    direction_tolerance: float
    percent_meaningful: bool = True


# One registry drives calculation, export ordering, GUI columns, and figure selection.
# Future axes can be appended here without changing the pair runner contract.
MULTI_PAIR_METRICS: tuple[MultiPairMetric, ...] = (
    MultiPairMetric("C_median", "C", "C program median", "score", 0.05, False),
    MultiPairMetric("S_median", "S", "S program median", "score", 0.05, False),
    MultiPairMetric("R_median", "R", "R = C - S median", "score", 0.05, False),
    MultiPairMetric("C_mean", "C_mean", "C program mean (compatibility)", "score", 0.05, False),
    MultiPairMetric("S_mean", "S_mean", "S program mean (compatibility)", "score", 0.05, False),
    MultiPairMetric("R_mean", "R_mean", "R mean (compatibility)", "score", 0.05, False),
    MultiPairMetric("localized_interface_fraction", "interface_fraction", "Localized interface-like fraction", "fraction", 0.01),
    MultiPairMetric("diffuse_fraction", "diffuse_fraction", "Diffuse transition fraction", "fraction", 0.01),
    MultiPairMetric("transition_burden_score", "transition_burden", "Transition burden score", "score", 0.02),
    MultiPairMetric("adj_same_fraction", "adj_same_fraction", "Same-side adjacency fraction", "fraction", 0.01),
    MultiPairMetric("adj_zero_fraction", "adj_zero_fraction", "Near-zero adjacency fraction", "fraction", 0.01),
    MultiPairMetric("adj_opposite_fraction", "adj_opposite_fraction", "Opposite-side adjacency fraction", "fraction", 0.01),
    MultiPairMetric("interface_fragmentation_index", "interface_fragmentation", "Interface fragmentation index", "ratio", 0.02),
    MultiPairMetric("largest_diffuse_component_ratio", "largest_diffuse_component_ratio", "Largest diffuse-component ratio", "ratio", 0.02),
    MultiPairMetric("small_component_fraction", "small_component_fraction", "Small-component fraction", "fraction", 0.01),
    MultiPairMetric(
        "diffuse_components_per_1000_valid_spots",
        "diffuse_components_per_1000_valid_spots",
        "Diffuse components per 1,000 valid spots",
        "per 1,000 spots",
        0.5,
    ),
    MultiPairMetric("H_raw_median", "H", "Hypoxia-associated expression context median", "expression score", 0.05),
    MultiPairMetric("H_raw_mean", "H_raw_mean", "H raw expression-context mean", "expression score", 0.05),
    MultiPairMetric("H_q75", "H_q75", "H expression-context q75", "expression score", 0.05),
    MultiPairMetric("H_q90", "H_q90", "H expression-context q90", "expression score", 0.05),
    MultiPairMetric("H_high_fraction", "H_high_fraction", "H pair-pooled high-context fraction", "fraction", 0.01),
    MultiPairMetric("H_local_hotspot_fraction", "H_local_hotspot_fraction", "H high-context local fraction", "fraction", 0.01),
    MultiPairMetric("H_transition_enrichment", "H_transition_enrichment", "H transition-context enrichment", "expression score", 0.05),
    MultiPairMetric("H_coefficient_of_variation", "H_coefficient_of_variation", "H coefficient of variation", "ratio", 0.05),
    MultiPairMetric("V_raw_median", "V", "Vascular-associated expression proxy median", "expression score", 0.05),
    MultiPairMetric("V_raw_mean", "V_raw_mean", "V raw expression-context mean", "expression score", 0.05),
    MultiPairMetric("V_q75", "V_q75", "V expression-context q75", "expression score", 0.05),
    MultiPairMetric("V_q90", "V_q90", "V expression-context q90", "expression score", 0.05),
    MultiPairMetric("V_high_fraction", "V_high_fraction", "V pair-pooled high-context fraction", "fraction", 0.01),
    MultiPairMetric("V_local_hotspot_fraction", "V_local_hotspot_fraction", "V high-context local fraction", "fraction", 0.01),
    MultiPairMetric("V_transition_enrichment", "V_transition_enrichment", "V transition-context enrichment", "expression score", 0.05),
    MultiPairMetric("V_coefficient_of_variation", "V_coefficient_of_variation", "V coefficient of variation", "ratio", 0.05),
)

BALANCE_METRIC_NAMES: tuple[str, ...] = ("C", "S", "R", "C_mean", "S_mean", "R_mean")
CONTEXT_METRIC_NAMES: tuple[str, ...] = (
    "H", "H_raw_mean", "H_q75", "H_q90", "H_high_fraction",
    "H_local_hotspot_fraction", "H_transition_enrichment", "H_coefficient_of_variation",
    "V", "V_raw_mean", "V_q75", "V_q90", "V_high_fraction",
    "V_local_hotspot_fraction", "V_transition_enrichment", "V_coefficient_of_variation",
)
SPATIAL_ORGANIZATION_METRIC_NAMES: tuple[str, ...] = tuple(
    metric.export_name
    for metric in MULTI_PAIR_METRICS
    if metric.export_name not in BALANCE_METRIC_NAMES + CONTEXT_METRIC_NAMES
)


@dataclass
class ComparabilityConfig:
    """Central, auditable thresholds for the specimen comparability gate."""

    spot_count_fold_caution: float = 1.50
    spot_count_fold_low: float = 2.50
    feature_count_fold_caution: float = 1.25
    feature_count_fold_low: float = 1.75
    detected_genes_fold_caution: float = 1.50
    detected_genes_fold_low: float = 2.50
    library_size_fold_caution: float = 2.00
    library_size_fold_low: float = 4.00
    occupancy_difference_caution: float = 0.15
    occupancy_difference_low: float = 0.30
    low_quality_difference_caution: float = 0.15
    low_quality_difference_low: float = 0.30
    valid_spot_difference_caution: float = 0.10
    valid_spot_difference_low: float = 0.25
    extent_area_fold_caution: float = 2.00
    extent_area_fold_low: float = 4.00
    tissue_component_fold_caution: float = 2.00
    tissue_component_fold_low: float = 4.00
    gene_coverage_caution: float = 0.80
    gene_coverage_low: float = 0.50
    composition_proxy_difference_caution: float = 0.25
    composition_proxy_difference_low: float = 0.45
    caution_count_for_low: int = 2

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class PairInterpretationConfig:
    """Auditable thresholds for descriptive pair-level change labels."""

    balance_moderate_abs_delta: float = 0.25
    balance_large_abs_delta: float = 0.75
    interface_moderate_abs_delta: float = 0.03
    interface_large_abs_delta: float = 0.10
    diffuse_moderate_abs_delta: float = 0.05
    diffuse_large_abs_delta: float = 0.15
    transition_burden_moderate_abs_delta: float = 0.10
    transition_burden_large_abs_delta: float = 0.30
    adjacency_moderate_abs_delta: float = 0.05
    adjacency_large_abs_delta: float = 0.15
    fragmentation_moderate_abs_delta: float = 0.10
    fragmentation_large_abs_delta: float = 0.30

    def validate(self) -> None:
        threshold_pairs = (
            ("balance", self.balance_moderate_abs_delta, self.balance_large_abs_delta),
            ("interface", self.interface_moderate_abs_delta, self.interface_large_abs_delta),
            ("diffuse", self.diffuse_moderate_abs_delta, self.diffuse_large_abs_delta),
            ("transition_burden", self.transition_burden_moderate_abs_delta, self.transition_burden_large_abs_delta),
            ("adjacency", self.adjacency_moderate_abs_delta, self.adjacency_large_abs_delta),
            ("fragmentation", self.fragmentation_moderate_abs_delta, self.fragmentation_large_abs_delta),
        )
        for label, moderate, large in threshold_pairs:
            if not 0 <= float(moderate) < float(large):
                raise ValueError(f"{label} interpretation thresholds must satisfy 0 <= moderate < large")

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class MultiPairRunResult:
    run_dir: Path
    pair_results: pd.DataFrame
    comparative_overview: pd.DataFrame
    comparability_qc: pd.DataFrame
    cohort_summary: pd.DataFrame
    figures: list[Path] = field(default_factory=list)
    summary_text: str = ""
    balance_changes: pd.DataFrame = field(default_factory=pd.DataFrame)
    spatial_organization_changes: pd.DataFrame = field(default_factory=pd.DataFrame)
    specimen_reliability: pd.DataFrame = field(default_factory=pd.DataFrame)
    pair_interpretation_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    comparability_details: pd.DataFrame = field(default_factory=pd.DataFrame)
    overview_interpretation: pd.DataFrame = field(default_factory=pd.DataFrame)
    context_changes: pd.DataFrame = field(default_factory=pd.DataFrame)
    context_gene_audit: pd.DataFrame = field(default_factory=pd.DataFrame)
    multiaxial_pair_summary: pd.DataFrame = field(default_factory=pd.DataFrame)
    comparative_qc_summary: pd.DataFrame = field(default_factory=pd.DataFrame)


Progress = Callable[[str], None]
CancelEvent = Event | None


def validate_pair_specs(pairs: list[PairSpec]) -> list[PairSpec]:
    if not 1 <= len(pairs) <= MAX_MULTI_PAIRS:
        raise ValueError(
            f"Multi-Pair Pre/Post requires between 1 and {MAX_MULTI_PAIRS} complete pairs."
        )
    normalized = [pair.normalized(index) for index, pair in enumerate(pairs, start=1)]
    labels = [pair.label.casefold() for pair in normalized]
    if len(labels) != len(set(labels)):
        raise ValueError("Pair labels must be unique.")
    # Per-pair file failures are intentionally handled inside the runner so one
    # unreadable or missing pair never terminates the other valid pairs.
    return normalized


def direction_symbol(delta: float, tolerance: float) -> str:
    if not np.isfinite(delta):
        return "NA"
    if abs(float(delta)) < float(tolerance):
        return "→"
    return "↑" if delta > 0 else "↓"


def _finite(value) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return np.nan
    return numeric if np.isfinite(numeric) else np.nan


def _fold_ratio(first, second) -> float:
    first_value, second_value = _finite(first), _finite(second)
    if not np.isfinite(first_value) or not np.isfinite(second_value):
        return np.nan
    if first_value < 0 or second_value < 0:
        return np.nan
    smaller, larger = min(first_value, second_value), max(first_value, second_value)
    if smaller <= 0:
        return 1.0 if larger <= 0 else np.inf
    return float(larger / smaller)


def _matrix_qc(adata) -> tuple[np.ndarray, np.ndarray]:
    from scipy import sparse

    matrix = adata.X
    totals = np.asarray(matrix.sum(axis=1)).reshape(-1)
    if sparse.issparse(matrix):
        detected = np.asarray((matrix > 0).sum(axis=1)).reshape(-1)
    else:
        values = np.asarray(matrix)
        detected = np.count_nonzero(np.isfinite(values) & (values > 0), axis=1)
    return np.asarray(totals, dtype=float), np.asarray(detected, dtype=float)


def _optional_low_quality_fraction(adata) -> tuple[float, str]:
    for column in ("low_quality", "qc_fail", "is_low_quality", "spatialtx_qc_fail"):
        if column not in adata.obs:
            continue
        values = adata.obs[column]
        if pd.api.types.is_bool_dtype(values):
            return float(values.fillna(False).mean()), column
        numeric = pd.to_numeric(values, errors="coerce")
        if numeric.notna().any():
            return float(numeric.fillna(0).ne(0).mean()), column
    return np.nan, "not_available"


def inspect_h5ad_qc(path: str | Path, c_genes: list[str], s_genes: list[str]) -> dict:
    """Read transparent technical QC inputs without inventing pathology labels."""
    resolved = Path(path).resolve()
    adata = _read_h5ad(resolved)
    totals, detected = _matrix_qc(adata)
    genes = {str(name).strip().upper() for name in adata.var_names}
    c_requested = {str(name).strip().upper() for name in c_genes if str(name).strip()}
    s_requested = {str(name).strip().upper() for name in s_genes if str(name).strip()}
    coords = np.asarray(adata.obsm.get("spatial", []), dtype=float)
    spatial_valid = bool(
        coords.ndim == 2
        and coords.shape[0] == int(adata.n_obs)
        and coords.shape[1] >= 2
        and coords.size
        and np.isfinite(coords[:, :2]).all()
    )
    in_tissue_fraction = np.nan
    if "in_tissue" in adata.obs:
        in_tissue = pd.to_numeric(adata.obs["in_tissue"], errors="coerce")
        if in_tissue.notna().any():
            in_tissue_fraction = float(in_tissue.eq(1).mean())
    low_quality_fraction, low_quality_source = _optional_low_quality_fraction(adata)
    result = {
        "path": str(resolved),
        "input_sha256": file_sha256(resolved),
        "n_spots": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "median_total_counts": float(np.nanmedian(totals)) if len(totals) else np.nan,
        "q25_total_counts": float(np.nanquantile(totals, 0.25)) if len(totals) else np.nan,
        "q75_total_counts": float(np.nanquantile(totals, 0.75)) if len(totals) else np.nan,
        "median_detected_genes_per_spot": float(np.nanmedian(detected)) if len(detected) else np.nan,
        "in_tissue_fraction": in_tissue_fraction,
        "low_quality_fraction": low_quality_fraction,
        "low_quality_source": low_quality_source,
        "spatial_coordinates_valid": spatial_valid,
        "C_gene_coverage": len(c_requested & genes) / max(1, len(c_requested)),
        "S_gene_coverage": len(s_requested & genes) / max(1, len(s_requested)),
        "matrix_values_note": "Observed adata.X values; library-size comparisons depend on preprocessing compatibility.",
    }
    return result


def _comparison_row(
    pair_label: str,
    metric: str,
    category: str,
    pre_value,
    post_value,
    comparison_value,
    availability: str,
    severity: str,
    reason: str,
    primary: bool,
) -> dict:
    return {
        "pair_label": pair_label,
        "qc_metric": metric,
        "category": category,
        "pre_value": pre_value,
        "post_value": post_value,
        "comparison_value": comparison_value,
        "availability_status": availability,
        "severity": severity,
        "primary_for_classification": bool(primary),
        "reason": reason,
    }


def evaluate_comparability(
    pair_label: str,
    pre_qc: dict | None,
    post_qc: dict | None,
    config: ComparabilityConfig,
    pre_metrics: dict | None = None,
    post_metrics: dict | None = None,
    pre_fields: dict | None = None,
    post_fields: dict | None = None,
) -> tuple[str, list[str], pd.DataFrame]:
    rows: list[dict] = []
    if pre_qc is None or post_qc is None:
        missing = "Pre" if pre_qc is None else "Post"
        reason = f"{missing} specimen could not be read for technical comparability QC"
        rows.append(_comparison_row(pair_label, "sample_readability", "technical", np.nan, np.nan, np.nan, "warning", "low", reason, True))
        return "Low", [reason], pd.DataFrame(rows)
    pre_qc = dict(pre_qc)
    post_qc = dict(post_qc)

    def ratio_check(name: str, category: str, caution: float, low: float, primary: bool = True) -> None:
        pre, post = pre_qc.get(name, np.nan), post_qc.get(name, np.nan)
        ratio = _fold_ratio(pre, post)
        if not np.isfinite(ratio):
            rows.append(_comparison_row(pair_label, name, category, pre, post, ratio, "not_available", "not_available", "Optional QC metric unavailable", primary))
            return
        severity = "low" if ratio >= low else "caution" if ratio >= caution else "ok"
        reason = f"{name} mismatch ({ratio:.2f}-fold)" if severity != "ok" else ""
        rows.append(_comparison_row(pair_label, name, category, pre, post, ratio, "warning" if severity != "ok" else "available", severity, reason, primary))

    def difference_check(name: str, category: str, caution: float, low: float, primary: bool = True) -> None:
        pre, post = _finite(pre_qc.get(name)), _finite(post_qc.get(name))
        difference = abs(post - pre) if np.isfinite(pre) and np.isfinite(post) else np.nan
        if not np.isfinite(difference):
            rows.append(_comparison_row(pair_label, name, category, pre, post, difference, "not_available", "not_available", "Optional QC metric unavailable", primary))
            return
        severity = "low" if difference >= low else "caution" if difference >= caution else "ok"
        reason = f"{name} mismatch (absolute difference {difference:.3f})" if severity != "ok" else ""
        rows.append(_comparison_row(pair_label, name, category, pre, post, difference, "warning" if severity != "ok" else "available", severity, reason, primary))

    pre_spatial = bool(pre_qc.get("spatial_coordinates_valid", False))
    post_spatial = bool(post_qc.get("spatial_coordinates_valid", False))
    spatial_ok = pre_spatial and post_spatial
    rows.append(_comparison_row(
        pair_label,
        "spatial_coordinates_valid",
        "technical",
        pre_spatial,
        post_spatial,
        int(spatial_ok),
        "available" if spatial_ok else "warning",
        "ok" if spatial_ok else "low",
        "missing or invalid spatial coordinates" if not spatial_ok else "",
        True,
    ))
    for side in ("C", "S"):
        name = f"{side}_gene_coverage"
        pre, post = _finite(pre_qc.get(name)), _finite(post_qc.get(name))
        minimum = min(pre, post) if np.isfinite(pre) and np.isfinite(post) else np.nan
        if not np.isfinite(minimum):
            severity, availability, reason = "not_available", "not_available", "Required gene coverage unavailable"
        elif minimum < config.gene_coverage_low:
            severity, availability, reason = "low", "warning", f"{side}-side required-gene coverage below {config.gene_coverage_low:.0%}"
        elif minimum < config.gene_coverage_caution:
            severity, availability, reason = "caution", "warning", f"{side}-side required-gene coverage below {config.gene_coverage_caution:.0%}"
        else:
            severity, availability, reason = "ok", "available", ""
        rows.append(_comparison_row(pair_label, name, "technical", pre, post, minimum, availability, severity, reason, True))

    ratio_check("n_spots", "sampling", config.spot_count_fold_caution, config.spot_count_fold_low)
    ratio_check("n_features", "technical", config.feature_count_fold_caution, config.feature_count_fold_low)
    ratio_check("median_detected_genes_per_spot", "technical", config.detected_genes_fold_caution, config.detected_genes_fold_low)
    ratio_check("median_total_counts", "technical", config.library_size_fold_caution, config.library_size_fold_low)
    ratio_check("q25_total_counts", "technical_context", config.library_size_fold_caution, config.library_size_fold_low, primary=False)
    ratio_check("q75_total_counts", "technical_context", config.library_size_fold_caution, config.library_size_fold_low, primary=False)
    difference_check("in_tissue_fraction", "sampling", config.occupancy_difference_caution, config.occupancy_difference_low)
    difference_check("low_quality_fraction", "technical", config.low_quality_difference_caution, config.low_quality_difference_low)

    if pre_metrics is not None and post_metrics is not None:
        for name in ("spatial_extent_area_proxy", "tissue_component_count"):
            pre_qc[name] = pre_metrics.get(name, np.nan)
            post_qc[name] = post_metrics.get(name, np.nan)
        for target, metrics in ((pre_qc, pre_metrics), (post_qc, post_metrics)):
            total = _finite(metrics.get("n_total_spots"))
            valid = _finite(metrics.get("n_valid_spots"))
            target["valid_spot_fraction"] = valid / total if np.isfinite(total) and total > 0 and np.isfinite(valid) else np.nan
        ratio_check("spatial_extent_area_proxy", "geometry_sampling", config.extent_area_fold_caution, config.extent_area_fold_low)
        ratio_check("tissue_component_count", "geometry_sampling", config.tissue_component_fold_caution, config.tissue_component_fold_low)
        difference_check(
            "valid_spot_fraction",
            "technical",
            config.valid_spot_difference_caution,
            config.valid_spot_difference_low,
        )

    # C/S distribution is explicitly secondary context because it may represent true
    # treatment-associated biology. It cannot by itself produce a Low classification.
    if pre_fields is not None and post_fields is not None:
        for label, fields, target in (("pre", pre_fields, pre_qc), ("post", post_fields, post_qc)):
            c_values = np.asarray(fields.get("C", []), dtype=float)
            s_values = np.asarray(fields.get("S", []), dtype=float)
            valid = np.isfinite(c_values) & np.isfinite(s_values)
            target["C_dominant_fraction_composition_proxy"] = float(np.mean(c_values[valid] > s_values[valid])) if valid.any() else np.nan
            target["S_dominant_fraction_composition_proxy"] = float(np.mean(s_values[valid] > c_values[valid])) if valid.any() else np.nan
        difference_check(
            "C_dominant_fraction_composition_proxy",
            "composition_proxy_secondary",
            config.composition_proxy_difference_caution,
            config.composition_proxy_difference_low,
            primary=False,
        )
        difference_check(
            "S_dominant_fraction_composition_proxy",
            "composition_proxy_secondary",
            config.composition_proxy_difference_caution,
            config.composition_proxy_difference_low,
            primary=False,
        )

    table = pd.DataFrame(rows)
    primary_low = table.loc[table["primary_for_classification"] & table["severity"].eq("low")]
    primary_caution = table.loc[table["primary_for_classification"] & table["severity"].eq("caution")]
    secondary_warning = table.loc[~table["primary_for_classification"] & table["severity"].isin(["caution", "low"])]
    if not primary_low.empty or len(primary_caution) >= int(config.caution_count_for_low):
        classification = "Low"
    elif not primary_caution.empty or not secondary_warning.empty:
        classification = "Caution"
    else:
        classification = "Good"
    reasons = [str(value) for value in table.loc[table["severity"].isin(["caution", "low"]), "reason"] if str(value).strip()]
    return classification, list(dict.fromkeys(reasons)), table


def _safe_slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._")
    return slug or "pair"


def _unique_run_dir(output_root: str | Path, tag: str | None = None) -> tuple[Path, str]:
    timestamp = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    stem = _safe_slug(tag or timestamp.replace(":", "").replace("-", ""))
    base = Path(output_root).expanduser().resolve() / "comparative_multi_pair"
    candidate = base / stem
    suffix = 2
    while candidate.exists():
        candidate = base / f"{stem}_{suffix}"
        suffix += 1
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate, timestamp


_NAMED_PAIR_ID = re.compile(
    r"(?:^|[_.-])(patient|subject|sample|case|pair|pt)[_.-]*([a-z0-9]+)(?=$|[_.-])",
    re.IGNORECASE,
)
_TIMEPOINT_SUFFIX_ID = re.compile(
    r"(?:^|[_.-])([a-z]{1,12}[_.-]?\d+)(?=[_.-](?:pre|post|baseline|pretreatment|posttreatment|before|after)(?:$|[_.-]))",
    re.IGNORECASE,
)


def infer_pair_id_from_filename(path: str | Path) -> str:
    """Infer only explicit patient/sample-like identifiers; accessions alone are not pair IDs."""
    stem = Path(path).stem.strip().lower()
    named = _NAMED_PAIR_ID.search(stem)
    if named:
        return f"{named.group(1).lower()}_{named.group(2).lower()}"
    trailing = _TIMEPOINT_SUFFIX_ID.search(stem)
    if trailing:
        candidate = re.sub(r"[_.-]+", "_", trailing.group(1).lower()).strip("_")
        compact = candidate.replace("_", "")
        if re.fullmatch(r"(?:gsm|gse|srr|sra|srx|samn|ena)\d+", compact, re.IGNORECASE):
            return ""
        return candidate
    return ""


def validate_pair_identity(pre_path: str | Path, post_path: str | Path) -> dict[str, str]:
    pre_id = infer_pair_id_from_filename(pre_path)
    post_id = infer_pair_id_from_filename(post_path)
    if pre_id and post_id and pre_id != post_id:
        warning = (
            f"Possible pair-ID mismatch: Pre suggests '{pre_id}' while Post suggests '{post_id}'. "
            "Please confirm this is an intended comparison."
        )
        status = "warning"
    elif pre_id and post_id:
        warning = ""
        status = "matched"
    else:
        warning = ""
        status = "not_available"
    return {
        "pre_inferred_pair_id": pre_id,
        "post_inferred_pair_id": post_id,
        "pair_id_validation": status,
        "pair_id_warning": warning,
    }


def _pair_result_row(pair: PairSpec) -> dict:
    row = {
        "pair_label": pair.label,
        "pre_file": str(pair.pre_path),
        "post_file": str(pair.post_path),
        "site_comparability": pair.site_comparability,
        "status": "ERROR",
        "error": "",
        "comparability": "Low",
        "comparability_reasons": "",
        "regime_pre": "",
        "regime_post": "",
        "regime_transition": "",
        "pattern_pre": "",
        "pattern_post": "",
        "pattern_transition": "",
    }
    row.update(validate_pair_identity(pair.pre_path, pair.post_path))
    return row


def _add_metric_results(row: dict, pre_metrics: dict, post_metrics: dict) -> None:
    for definition in MULTI_PAIR_METRICS:
        pre = _finite(pre_metrics.get(definition.internal_name))
        post = _finite(post_metrics.get(definition.internal_name))
        delta = post - pre if np.isfinite(pre) and np.isfinite(post) else np.nan
        percent, percent_status = percent_change(pre, post) if definition.percent_meaningful else (np.nan, "not_meaningful")
        prefix = definition.export_name
        row[f"pre_{prefix}"] = pre
        row[f"post_{prefix}"] = post
        row[f"delta_{prefix}"] = delta
        row[f"percent_change_{prefix}"] = percent
        row[f"percent_change_status_{prefix}"] = percent_status
        row[f"direction_{prefix}"] = direction_symbol(delta, definition.direction_tolerance)


def _apply_within_pair_context_thresholds(
    pre_metrics: dict,
    post_metrics: dict,
    pre_fields: dict[str, np.ndarray],
    post_fields: dict[str, np.ndarray],
) -> tuple[dict, dict]:
    """Apply one pooled q90 per axis within this Pre/Post pair only."""
    sample_metrics = pd.DataFrame([pre_metrics, post_metrics])
    fields_by_sample = {
        str(pre_metrics["sample_id"]): pre_fields,
        str(post_metrics["sample_id"]): post_fields,
    }
    pooled = apply_pooled_hv_thresholds(
        sample_metrics,
        fields_by_sample,
        "Pre",
        "Post",
    ).set_index("sample_id", drop=False)
    pre_updated = pooled.loc[str(pre_metrics["sample_id"])].to_dict()
    post_updated = pooled.loc[str(post_metrics["sample_id"])].to_dict()
    return pre_updated, post_updated


def _add_context_audit_results(row: dict, pre_metrics: dict, post_metrics: dict) -> None:
    for axis in ("H", "V"):
        pre_audit = dict(pre_metrics.get(f"{axis}_context_audit", {}))
        post_audit = dict(post_metrics.get(f"{axis}_context_audit", {}))
        row[f"pre_{axis}_gene_coverage"] = _finite(pre_audit.get("coverage_fraction"))
        row[f"post_{axis}_gene_coverage"] = _finite(post_audit.get("coverage_fraction"))
        row[f"pre_{axis}_context_status"] = str(
            pre_metrics.get(f"{axis}_context_status", "calculation_error")
        )
        row[f"post_{axis}_context_status"] = str(
            post_metrics.get(f"{axis}_context_status", "calculation_error")
        )
        row[f"pre_{axis}_raw_normalization_method"] = str(
            pre_metrics.get(f"{axis}_raw_normalization_method", "")
        )
        row[f"post_{axis}_raw_normalization_method"] = str(
            post_metrics.get(f"{axis}_raw_normalization_method", "")
        )
        pre_threshold = _finite(pre_metrics.get(f"{axis}_pooled_high_threshold"))
        post_threshold = _finite(post_metrics.get(f"{axis}_pooled_high_threshold"))
        if np.isfinite(pre_threshold) and np.isfinite(post_threshold) and np.isclose(pre_threshold, post_threshold):
            pair_threshold = pre_threshold
        else:
            pair_threshold = np.nan
        row[f"{axis}_pair_pooled_q90"] = pair_threshold
        row[f"{axis}_pair_pooled_threshold_scope"] = "within_pair_pre_plus_post"


def _join_gene_list(value) -> str:
    if isinstance(value, (list, tuple)):
        return ";".join(str(item) for item in value)
    return str(value or "")


def _context_gene_audit_row(
    pair: PairSpec,
    sample_role: str,
    sample_file: Path,
    axis: str,
    metrics: dict | None,
    config: ComparativeConfig,
) -> dict:
    configured_genes = config.h_genes if axis == "H" else config.v_genes
    enabled = config.enable_h_expr if axis == "H" else config.enable_v_expr
    effective = resolve_context_program(axis, configured_genes)
    audit = dict((metrics or {}).get(f"{axis}_context_audit", {}))
    if not audit:
        audit = {
            **effective,
            "matched_genes": [],
            "missing_genes": [],
            "matched_gene_count": 0,
            "coverage_fraction": np.nan,
            "genes_expressed_above_min_spot_fraction": [],
            "expressed_gene_count": 0,
            "expressed_gene_fraction": np.nan,
            "expression_scale_guess": "unknown",
            "detection_source": "sample_analysis_not_completed",
            "raw_normalization_method": "unavailable_sample_analysis_not_completed",
        }
    status = str((metrics or {}).get(
        f"{axis}_context_status",
        "not_requested" if not enabled else "calculation_error",
    ))
    sample_metrics = metrics or {}
    single_sample_warning = str(sample_metrics.get(
        f"{axis}_context_warning",
        "Not requested." if not enabled else "Sample analysis did not complete.",
    ))
    single_sample_high_fraction = _finite(
        sample_metrics.get(f"{axis}_expr_high_state_fraction")
    )
    pair_pooled_high_fraction = _finite(sample_metrics.get(f"{axis}_high_fraction"))
    pair_pooled_q90 = _finite(sample_metrics.get(f"{axis}_pooled_high_threshold"))
    pair_pooled_warning = ""
    if status == "available" and not (
        np.isfinite(pair_pooled_high_fraction) and np.isfinite(pair_pooled_q90)
    ):
        pair_pooled_warning = (
            "Pair-pooled raw-context q90 summary was unavailable despite available context status."
        )

    warning_parts = []
    if single_sample_warning:
        warning_parts.append(
            "Single-sample context QC (legacy within-sample centered-context q80): "
            f"{single_sample_warning}"
        )
    if pair_pooled_warning:
        warning_parts.append(
            "Pair-pooled context QC (within-pair raw-context q90): "
            f"{pair_pooled_warning}"
        )
    expressed = audit.get("genes_expressed_above_min_spot_fraction", [])
    return {
        "pair_label": pair.label,
        "sample_role": sample_role,
        "sample_file": str(sample_file),
        "axis": axis,
        "gene_set_name": audit.get("gene_set_name", effective["gene_set_name"]),
        "gene_set_source": audit.get("source", effective["source"]),
        "requested_gene_count": audit.get("requested_gene_count", effective["requested_gene_count"]),
        "matched_gene_count": audit.get("matched_gene_count", 0),
        "coverage_fraction": audit.get("coverage_fraction", np.nan),
        "requested_genes": _join_gene_list(audit.get("requested_genes", effective["requested_genes"])),
        "matched_genes": _join_gene_list(audit.get("matched_genes", [])),
        "missing_genes": _join_gene_list(audit.get("missing_genes", [])),
        "expressed_gene_count": audit.get("expressed_gene_count", len(expressed)),
        "expressed_gene_fraction": audit.get("expressed_gene_fraction", np.nan),
        "expressed_genes": _join_gene_list(expressed),
        "expression_scale_guess": audit.get("expression_scale_guess", "unknown"),
        "detection_source": audit.get("detection_source", ""),
        "raw_normalization_method": audit.get("raw_normalization_method", ""),
        "context_available": status == "available",
        "context_status": status,
        "single_sample_context_high_fraction": single_sample_high_fraction,
        "single_sample_context_threshold_scope": "within_sample_centered_context_q80",
        "single_sample_context_warning": single_sample_warning,
        "single_sample_context_warning_provenance": "legacy_within_sample_centered_context_q80",
        "pair_pooled_context_high_fraction": pair_pooled_high_fraction,
        "pair_pooled_context_q90_threshold": pair_pooled_q90,
        "pair_pooled_context_threshold_scope": "within_pair_pre_plus_post_raw_context_q90",
        "pair_pooled_context_warning": pair_pooled_warning,
        "pair_pooled_context_warning_provenance": "within_pair_pre_plus_post_raw_context_q90",
        "context_warning": " | ".join(warning_parts),
    }


def _format_delta_summary(row: dict, metric: str) -> str:
    delta = _finite(row.get(f"delta_{metric}"))
    direction = str(row.get(f"direction_{metric}", "NA"))
    return f"{direction} {delta:+.4g}" if np.isfinite(delta) else "NA"


def _overview_row(row: dict) -> dict:
    columns = (
        "pair_label", "status", "error", "comparability", "comparability_reasons",
        "delta_C", "direction_C", "delta_S", "direction_S", "delta_R", "direction_R",
        "delta_interface_fraction", "direction_interface_fraction",
        "delta_diffuse_fraction", "direction_diffuse_fraction",
        "delta_transition_burden", "direction_transition_burden",
        "delta_H", "direction_H", "delta_V", "direction_V", "site_comparability",
        "regime_transition", "pattern_transition",
    )
    overview = {column: row.get(column, np.nan if column.startswith("delta_") else "") for column in columns}
    overview.update({
        "balance_change_summary": "; ".join(
            f"{metric} {_format_delta_summary(row, metric)}" for metric in ("C", "S", "R")
        ),
        "spatial_organization_summary": "; ".join((
            f"Interface {_format_delta_summary(row, 'interface_fraction')}",
            f"Diffuse {_format_delta_summary(row, 'diffuse_fraction')}",
            f"Burden {_format_delta_summary(row, 'transition_burden')}",
        )),
        "derived_pattern_transition": str(row.get("regime_transition", "")),
        "specimen_reliability": str(row.get("comparability", "")),
        "specimen_reliability_reasons": str(row.get("comparability_reasons", "")),
        "context_axes_summary": "; ".join((
            f"H {_format_delta_summary(row, 'H')}",
            f"V {_format_delta_summary(row, 'V')}",
        )),
        "site_comparability_summary": str(row.get("site_comparability", "unknown_site")),
    })
    return overview


def _change_layer_table(pair_results: pd.DataFrame, layer: str, metric_names: tuple[str, ...]) -> pd.DataFrame:
    base_columns = ["pair_label", "status", "error", "pre_file", "post_file", "site_comparability"]
    metric_columns: list[str] = []
    for metric in metric_names:
        metric_columns.extend((
            f"pre_{metric}",
            f"post_{metric}",
            f"delta_{metric}",
            f"percent_change_{metric}",
            f"percent_change_status_{metric}",
            f"direction_{metric}",
        ))
    extra_columns = ["regime_pre", "regime_post", "regime_transition", "pattern_pre", "pattern_post", "pattern_transition"]
    columns = base_columns + metric_columns + (extra_columns if layer == "Spatial organization change" else [])
    table = pair_results.reindex(columns=columns).copy()
    table.insert(1, "result_layer", layer)
    return table


def _specimen_reliability_table(pair_results: pd.DataFrame, qc_table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for _, pair in pair_results.iterrows():
        if not qc_table.empty and "pair_label" in qc_table:
            qc = qc_table.loc[qc_table["pair_label"].eq(pair.get("pair_label"))].copy()
        else:
            qc = pd.DataFrame()
        if qc.empty:
            primary_low_count = primary_caution_count = secondary_warning_count = unavailable_count = 0
        else:
            primary = qc["primary_for_classification"].fillna(False).astype(bool)
            severity = qc["severity"].fillna("").astype(str)
            primary_low_count = int((primary & severity.eq("low")).sum())
            primary_caution_count = int((primary & severity.eq("caution")).sum())
            secondary_warning_count = int((
                ~primary
                & severity.isin(("caution", "low"))
                & ~qc["category"].fillna("").astype(str).eq("pair_validation")
            ).sum())
            unavailable_count = int(qc["availability_status"].fillna("").astype(str).eq("not_available").sum())
        classification = str(pair.get("comparability", "Low"))
        if str(pair.get("status", "")) != "PASS":
            note = "Pair analysis did not complete; review the error and QC details."
        elif classification == "Low":
            note = "Technical or sampling differences may substantially influence the observed changes."
        elif classification == "Caution":
            note = "Review the listed QC differences before interpreting the observed changes."
        else:
            note = "No configured comparability warning was detected; this does not validate biological causality."
        rows.append({
            "pair_label": pair.get("pair_label", ""),
            "result_layer": "Specimen reliability",
            "status": pair.get("status", ""),
            "specimen_reliability": classification,
            "reliability_reasons": pair.get("comparability_reasons", ""),
            "primary_low_count": primary_low_count,
            "primary_caution_count": primary_caution_count,
            "secondary_warning_count": secondary_warning_count,
            "unavailable_qc_count": unavailable_count,
            "interpretation_note": note,
            "pair_id_validation": pair.get("pair_id_validation", "not_available"),
            "pair_id_warning": pair.get("pair_id_warning", ""),
            "site_comparability": pair.get("site_comparability", "unknown_site"),
            "site_warning": site_warning_text(str(pair.get("site_comparability", "unknown_site"))),
            "error": pair.get("error", ""),
            "pre_file": pair.get("pre_file", ""),
            "post_file": pair.get("post_file", ""),
        })
    return pd.DataFrame(rows)


def _qc_reason_summary(qc: pd.DataFrame, categories: tuple[str, ...]) -> str:
    if qc.empty:
        return ""
    selected = qc.loc[
        qc["category"].astype(str).isin(categories)
        & qc["severity"].astype(str).isin(("caution", "low")),
        "reason",
    ]
    return "; ".join(dict.fromkeys(str(value) for value in selected if str(value).strip()))


def _qc_values(qc: pd.DataFrame, metric: str) -> tuple[float, float, float]:
    if qc.empty:
        return np.nan, np.nan, np.nan
    selected = qc.loc[qc["qc_metric"].astype(str).eq(metric)]
    if selected.empty:
        return np.nan, np.nan, np.nan
    row = selected.iloc[0]
    return _finite(row.get("pre_value")), _finite(row.get("post_value")), _finite(row.get("comparison_value"))


def _comparability_details_table(pair_results: pd.DataFrame, qc_table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for _, pair in pair_results.iterrows():
        if not qc_table.empty and "pair_label" in qc_table:
            qc = qc_table.loc[qc_table["pair_label"].eq(pair.get("pair_label"))].copy()
        else:
            qc = pd.DataFrame()
        spot_pre, spot_post, spot_ratio = _qc_values(qc, "n_spots")
        genes_pre, genes_post, genes_ratio = _qc_values(qc, "median_detected_genes_per_spot")
        counts_pre, counts_post, counts_ratio = _qc_values(qc, "median_total_counts")
        q25_pre, q25_post, q25_ratio = _qc_values(qc, "q25_total_counts")
        extent_pre, extent_post, extent_ratio = _qc_values(qc, "spatial_extent_area_proxy")
        components_pre, components_post, components_ratio = _qc_values(qc, "tissue_component_count")
        occupancy_pre, occupancy_post, occupancy_difference = _qc_values(qc, "in_tissue_fraction")
        technical_reason = _qc_reason_summary(qc, ("technical", "technical_context"))
        sampling_reason = _qc_reason_summary(qc, ("sampling", "geometry_sampling"))
        composition_reason = _qc_reason_summary(qc, ("composition_proxy_secondary",))
        sampling_composition = "; ".join(
            value for value in (sampling_reason, composition_reason) if value
        )
        rows.append({
            "pair_label": pair.get("pair_label", ""),
            "comparability": pair.get("comparability", "Low"),
            "spot_count_pre": spot_pre,
            "spot_count_post": spot_post,
            "spot_count_ratio": spot_ratio,
            "detected_gene_median_pre": genes_pre,
            "detected_gene_median_post": genes_post,
            "detected_gene_ratio": genes_ratio,
            "median_counts_pre": counts_pre,
            "median_counts_post": counts_post,
            "median_counts_ratio": counts_ratio,
            "q25_counts_pre": q25_pre,
            "q25_counts_post": q25_post,
            "q25_counts_ratio": q25_ratio,
            "spatial_extent_pre": extent_pre,
            "spatial_extent_post": extent_post,
            "spatial_extent_ratio": extent_ratio,
            "tissue_component_count_pre": components_pre,
            "tissue_component_count_post": components_post,
            "tissue_component_count_difference": (
                abs(components_post - components_pre)
                if np.isfinite(components_pre) and np.isfinite(components_post)
                else np.nan
            ),
            "tissue_component_count_ratio": components_ratio,
            "tissue_occupancy_pre": occupancy_pre,
            "tissue_occupancy_post": occupancy_post,
            "tissue_occupancy_difference": occupancy_difference,
            "technical_mismatch_reason": technical_reason,
            "sampling_mismatch_reason": sampling_reason,
            "composition_proxy_reason": composition_reason,
            "sampling_composition_mismatch_summary": sampling_composition,
            "pair_id_validation": pair.get("pair_id_validation", "not_available"),
            "pair_id_warning": pair.get("pair_id_warning", ""),
            "all_comparability_reasons": pair.get("comparability_reasons", ""),
            "site_comparability": pair.get("site_comparability", "unknown_site"),
            "site_warning": site_warning_text(str(pair.get("site_comparability", "unknown_site"))),
        })
    return pd.DataFrame(rows)


def _balance_change_class(row: dict, config: PairInterpretationConfig) -> tuple[str, str]:
    values = [(name, abs(_finite(row.get(f"delta_{name}")))) for name in ("C", "S", "R")]
    values = [(name, value) for name, value in values if np.isfinite(value)]
    if not values:
        return "Unavailable", "Balance metrics unavailable"
    metric, maximum = max(values, key=lambda item: item[1])
    if maximum >= config.balance_large_abs_delta:
        label = "Large"
    elif maximum >= config.balance_moderate_abs_delta:
        label = "Moderate"
    else:
        label = "Minimal"
    return label, f"max absolute Delta={maximum:.4g} ({metric})"


def _spatial_change_class(row: dict, config: PairInterpretationConfig) -> tuple[str, str]:
    specifications = (
        ("interface_fraction", "interface", config.interface_moderate_abs_delta, config.interface_large_abs_delta),
        ("diffuse_fraction", "diffuse", config.diffuse_moderate_abs_delta, config.diffuse_large_abs_delta),
        ("transition_burden", "transition burden", config.transition_burden_moderate_abs_delta, config.transition_burden_large_abs_delta),
        ("adj_same_fraction", "same-side adjacency", config.adjacency_moderate_abs_delta, config.adjacency_large_abs_delta),
        ("adj_zero_fraction", "near-zero adjacency", config.adjacency_moderate_abs_delta, config.adjacency_large_abs_delta),
        ("adj_opposite_fraction", "opposite-side adjacency", config.adjacency_moderate_abs_delta, config.adjacency_large_abs_delta),
        ("interface_fragmentation", "fragmentation", config.fragmentation_moderate_abs_delta, config.fragmentation_large_abs_delta),
    )
    observed: list[tuple[str, float, float, float]] = []
    for metric, label, moderate, large in specifications:
        value = abs(_finite(row.get(f"delta_{metric}")))
        if np.isfinite(value):
            observed.append((label, value, moderate, large))
    if not observed:
        return "Unavailable", "Spatial organization metrics unavailable"
    large_hits = [(label, value) for label, value, _moderate, large in observed if value >= large]
    moderate_hits = [(label, value) for label, value, moderate, _large in observed if value >= moderate]
    if large_hits:
        label = "Large"
        basis_hits = large_hits
    elif moderate_hits:
        label = "Moderate"
        basis_hits = moderate_hits
    else:
        label = "Minimal"
        basis_hits = [(name, value) for name, value, _moderate, _large in observed]
    basis_hits = sorted(basis_hits, key=lambda item: item[1], reverse=True)[:3]
    basis = "; ".join(f"abs Delta {name}={value:.4g}" for name, value in basis_hits)
    return label, basis


def _regime_preservation(row: dict) -> str:
    pre = str(row.get("regime_pre", "")).strip()
    post = str(row.get("regime_post", "")).strip()
    unavailable = {"", "nan", "Spatial_QC_incomplete"}
    if pre in unavailable or post in unavailable:
        return "not_available"
    return "yes" if pre == post else "no"


def _pair_interpretation_table(
    pair_results: pd.DataFrame,
    comparability_details: pd.DataFrame,
    config: PairInterpretationConfig,
) -> pd.DataFrame:
    detail_by_pair = {
        str(row["pair_label"]): row for _, row in comparability_details.iterrows()
    }
    rows: list[dict] = []
    for _, pair_series in pair_results.iterrows():
        pair = pair_series.to_dict()
        pair_label = str(pair.get("pair_label", ""))
        detail = detail_by_pair.get(pair_label, {})
        balance_class, balance_basis = _balance_change_class(pair, config)
        spatial_class, spatial_basis = _spatial_change_class(pair, config)
        regime_preserved = _regime_preservation(pair)
        comparability = str(pair.get("comparability", "Low"))
        if str(pair.get("status", "")) != "PASS" or spatial_class == "Unavailable":
            structure_preserved = "not_available"
            structure_note = "Spatial structure could not be assessed."
        elif regime_preserved == "no":
            structure_preserved = "no"
            structure_note = "Broad spatial regime changed; this remains a descriptive categorical transition."
        elif spatial_class == "Large":
            structure_preserved = "no"
            structure_note = "Same broad regime may coexist with a substantial internal spatial-organization shift."
        elif spatial_class == "Moderate":
            structure_preserved = "probably"
            structure_note = "Same regime with an internal structural shift."
        elif comparability == "Low":
            structure_preserved = "probably"
            structure_note = "Broad regime and measured spatial structure appear largely preserved, with low comparability."
        else:
            structure_preserved = "yes"
            structure_note = "Broad regime and measured spatial structure are largely preserved."

        confidence = {"Good": "higher-confidence", "Caution": "moderate-confidence", "Low": "low-confidence"}.get(
            comparability,
            "low-confidence",
        )
        if spatial_class == "Large" or regime_preserved == "no":
            status_prefix = "Strong shift"
        elif spatial_class == "Moderate":
            status_prefix = "Possible structural shift"
        elif spatial_class == "Minimal":
            status_prefix = "Stable structure"
        else:
            status_prefix = "Structure unavailable"
        pair_status_flag = f"{status_prefix} / {confidence} pair; balance change: {balance_class}"

        if balance_class == "Minimal" and spatial_class == "Minimal":
            interpretive_flag = "minimal overall observed change"
        elif regime_preserved == "yes" and spatial_class == "Minimal" and balance_class in ("Moderate", "Large"):
            interpretive_flag = "balance shift with preserved structure"
        elif spatial_class == "Large":
            interpretive_flag = "possible strong structural shift"
        elif spatial_class == "Moderate" or regime_preserved == "no":
            interpretive_flag = "possible structural shift"
        else:
            interpretive_flag = "observed quantitative change"
        if comparability == "Low":
            interpretive_flag += ", low comparability"

        caution_message = (
            "CAUTION: Observed pre/post differences may be substantially influenced by specimen-selection, "
            "sampling, or technical differences."
            if comparability == "Low"
            else ""
        )
        note_message = (
            "NOTE: Broad spatial regime is preserved. Observed differences mainly reflect quantitative balance "
            "shifts rather than a categorical regime transition."
            if regime_preserved == "yes" and spatial_class == "Minimal" and balance_class in ("Moderate", "Large")
            else ""
        )
        interpretation_confidence, qc_aware_note = qc_aware_interpretation(pair)
        site_comparability = normalize_site_comparability(
            str(pair.get("site_comparability", "unknown_site"))
        )
        interpretation_sentences = [
            interpretive_flag.capitalize() + ".",
            structure_note,
            qc_aware_note,
        ]
        if caution_message:
            interpretation_sentences.append(caution_message)
        if note_message:
            interpretation_sentences.append(note_message)
        interpretation_text = (
            f"PAIR: {pair_label}\n"
            f"Balance change: {balance_class} ({balance_basis})\n"
            f"Spatial organization change: {spatial_class} ({spatial_basis})\n"
            f"Regime preservation: {pair.get('regime_transition', '')} | preserved={regime_preserved}\n"
            f"Comparability: {comparability}\n"
            f"Interpretation confidence: {interpretation_confidence}\n"
            f"Site comparability: {site_comparability}\n\n"
            "Interpretation:\n" + "\n".join(interpretation_sentences)
        )
        rows.append({
            "pair_label": pair_label,
            "status": pair.get("status", ""),
            "balance_change_class": balance_class,
            "balance_change_basis": balance_basis,
            "spatial_change_class": spatial_class,
            "spatial_change_basis": spatial_basis,
            "regime_transition": pair.get("regime_transition", ""),
            "regime_preserved": regime_preserved,
            "structure_preserved": structure_preserved,
            "structure_preservation_note": structure_note,
            "comparability": comparability,
            "interpretation_confidence": interpretation_confidence,
            "site_comparability": site_comparability,
            "site_shift_warning": site_warning_text(site_comparability),
            "qc_aware_interpretation_note": qc_aware_note,
            "pair_status_flag": pair_status_flag,
            "interpretive_flag": interpretive_flag,
            "technical_mismatch_reason": detail.get("technical_mismatch_reason", ""),
            "sampling_mismatch_reason": detail.get("sampling_mismatch_reason", ""),
            "composition_proxy_reason": detail.get("composition_proxy_reason", ""),
            "pair_id_validation": pair.get("pair_id_validation", "not_available"),
            "pair_id_warning": pair.get("pair_id_warning", ""),
            "caution_message": caution_message,
            "note_message": note_message,
            "interpretation_text": interpretation_text,
            "error": pair.get("error", ""),
        })
    return pd.DataFrame(rows)


def _cohort_summary(overview: pd.DataFrame) -> pd.DataFrame:
    passed = overview.loc[overview["status"].eq("PASS")] if not overview.empty else overview
    total = int(len(passed))
    rows = [
        {"summary_metric": "R increased", "count": int(passed.get("direction_R", pd.Series(dtype=str)).eq("↑").sum()), "eligible_pairs": total},
        {"summary_metric": "R decreased", "count": int(passed.get("direction_R", pd.Series(dtype=str)).eq("↓").sum()), "eligible_pairs": total},
        {"summary_metric": "Diffuse fraction decreased", "count": int(passed.get("direction_diffuse_fraction", pd.Series(dtype=str)).eq("↓").sum()), "eligible_pairs": total},
        {"summary_metric": "Transition burden decreased", "count": int(passed.get("direction_transition_burden", pd.Series(dtype=str)).eq("↓").sum()), "eligible_pairs": total},
        {"summary_metric": "Localized interface increased", "count": int(passed.get("direction_interface_fraction", pd.Series(dtype=str)).eq("↑").sum()), "eligible_pairs": total},
        {
            "summary_metric": "B/C to A descriptive transitions",
            "count": int(passed.get("regime_transition", pd.Series(dtype=str)).astype(str).str.match(r"Type_[BC]_candidate → Type_A_candidate").sum()),
            "eligible_pairs": total,
        },
    ]
    for label in ("Good", "Caution", "Low"):
        rows.append({
            "summary_metric": f"Comparability: {label}",
            "count": int(overview.get("comparability", pd.Series(dtype=str)).eq(label).sum()),
            "eligible_pairs": int(len(overview)),
        })
    return pd.DataFrame(rows)


def _write_csv(path: Path, table: pd.DataFrame) -> None:
    table.to_csv(path, index=False, encoding="utf-8-sig", na_rep="NA")


def _plot_overview(overview: pd.DataFrame, output: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    columns = [
        ("pair_label", "Pair"),
        ("comparability", "L3 Reliability"),
        ("balance_change_class", "L1 Class"),
        ("spatial_change_class", "L2 Class"),
        ("regime_transition", "Regime"),
        ("regime_preserved", "Regime preserved"),
        ("delta_C", "L1 ΔC"),
        ("delta_S", "L1 ΔS"),
        ("delta_R", "L1 ΔR"),
        ("delta_interface_fraction", "L2 ΔInterface"),
        ("delta_diffuse_fraction", "L2 ΔDiffuse"),
        ("delta_transition_burden", "L2 ΔBurden"),
        ("interpretive_flag", "Interpretive flag"),
    ]
    display = overview[[name for name, _label in columns]].copy() if not overview.empty else pd.DataFrame(columns=[name for name, _label in columns])
    for name, _label in columns:
        if name.startswith("delta_"):
            display[name] = display[name].map(
                lambda value: f"{_finite(value):+.3g}" if np.isfinite(_finite(value)) else "NA"
            )
    display = display.fillna("NA").astype(str)
    fig_height = max(3.2, 1.0 + 0.62 * max(1, len(display)))
    fig, ax = plt.subplots(figsize=(20.5, fig_height))
    ax.axis("off")
    table = ax.table(
        cellText=display.to_numpy(),
        colLabels=[label for _name, label in columns],
        cellLoc="center",
        colLoc="center",
        colWidths=[0.055, 0.075, 0.07, 0.07, 0.18, 0.075, 0.055, 0.055, 0.055, 0.075, 0.075, 0.075, 0.23],
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.6)
    header_colors = [
        "#334155", "#b45309", "#1d4ed8", "#047857", "#475569", "#047857",
        "#1d4ed8", "#1d4ed8", "#1d4ed8", "#047857", "#047857", "#047857", "#475569",
    ]
    for column in range(len(columns)):
        table[(0, column)].set_facecolor(header_colors[column])
        table[(0, column)].set_text_props(color="white", weight="bold")
    comparability_column = 1
    colors = {"Good": "#dcfce7", "Caution": "#fef3c7", "Low": "#fee2e2"}
    for row_index, status in enumerate(display.get("comparability", []), start=1):
        table[(row_index, comparability_column)].set_facecolor(colors.get(str(status), "#f3f4f6"))
    ax.set_title("Three-layer Multi-Pair Pre/Post overview", fontsize=15, weight="bold", pad=28)
    fig.text(
        0.5,
        0.91,
        "Layer 1: Balance change  |  Layer 2: Spatial organization change  |  Layer 3: Specimen reliability",
        ha="center",
        fontsize=9,
        color="#334155",
    )
    fig.text(0.01, 0.015, f"SpatialTX Studio v{__version__} | Separate descriptive layers; no composite response score. | {MULTI_PAIR_NOTICE}", fontsize=7, color="#4b5563")
    fig.tight_layout(rect=(0, 0.05, 1, 0.88))
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output


def _plot_pair_panel(row: dict, output: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = [
        ("Layer 1: Balance", "C", "C"),
        ("Layer 1: Balance", "S", "S"),
        ("Layer 1: Balance", "R", "R"),
        ("Layer 2: Spatial", "interface_fraction", "Interface"),
        ("Layer 2: Spatial", "diffuse_fraction", "Diffuse"),
        ("Layer 2: Spatial", "transition_burden", "Burden"),
        ("Context axis", "H", "H context"),
        ("Context axis", "V", "V context"),
    ]
    cells = []
    for layer, key, label in metrics:
        cells.append([
            layer,
            label,
            f"{_finite(row.get(f'pre_{key}')):.4g}" if np.isfinite(_finite(row.get(f"pre_{key}"))) else "NA",
            f"{_finite(row.get(f'post_{key}')):.4g}" if np.isfinite(_finite(row.get(f"post_{key}"))) else "NA",
            f"{_finite(row.get(f'delta_{key}')):+.4g}" if np.isfinite(_finite(row.get(f"delta_{key}"))) else "NA",
            str(row.get(f"direction_{key}", "NA")),
        ])
    fig, ax = plt.subplots(figsize=(10.4, 6.4))
    ax.axis("off")
    table = ax.table(
        cellText=cells,
        colLabels=["Result layer", "Metric", "Pre", "Post", "Delta", "Direction"],
        cellLoc="center",
        colWidths=[0.22, 0.16, 0.14, 0.14, 0.16, 0.14],
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.55)
    for column in range(6):
        table[(0, column)].set_facecolor("#0f3d56")
        table[(0, column)].set_text_props(color="white", weight="bold")
    for row_index in range(1, 4):
        table[(row_index, 0)].set_facecolor("#dbeafe")
    for row_index in range(4, 7):
        table[(row_index, 0)].set_facecolor("#d1fae5")
    for row_index in range(7, 9):
        table[(row_index, 0)].set_facecolor("#f3e8ff")
    ax.set_title(
        f"{row['pair_label']} | Three separate result layers\n"
        f"Layer 1 balance: {row.get('balance_change_class', 'Unavailable')} | "
        f"Layer 2 spatial: {row.get('spatial_change_class', 'Unavailable')} | "
        f"Layer 3 reliability: {row['comparability']}\n"
        f"{row.get('pair_status_flag', '')} | Regime: {row.get('regime_transition', '')}\n"
        f"Interpretation confidence: {row.get('interpretation_confidence', 'LOW')} | "
        f"Site: {row.get('site_comparability', 'unknown_site')}",
        fontsize=14,
        weight="bold",
        pad=15,
    )
    warning_parts: list[str] = []
    if row["comparability"] == "Low":
        warning_parts.append(
            str(row.get("caution_message", "")).strip()
            or "CAUTION: Observed pre/post differences may be substantially influenced by specimen-selection, sampling, or technical differences."
        )
    if row.get("site_comparability") == "different_site":
        warning_parts.append(site_warning_text("different_site"))
    warning = " ".join(warning_parts) or MULTI_PAIR_NOTICE
    fig.text(0.5, 0.055, "Balance, spatial organization, and specimen reliability are not combined into one score.", ha="center", fontsize=8, color="#334155")
    fig.text(0.5, 0.025, warning, ha="center", fontsize=8, color="#9a3412" if warning_parts else "#4b5563")
    fig.tight_layout(rect=(0, 0.10, 1, 0.88))
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output


def run_multi_pair_analysis(
    pairs: list[PairSpec],
    analysis_config: ComparativeConfig,
    output_root: str | Path,
    comparability_config: ComparabilityConfig | None = None,
    progress: Progress | None = None,
    cancel_event: CancelEvent = None,
    run_tag: str | None = None,
    interpretation_config: PairInterpretationConfig | None = None,
) -> MultiPairRunResult:
    pairs = validate_pair_specs(pairs)
    analysis_config.validate()
    if analysis_config.mode != "pairwise":
        raise ValueError("Multi-Pair Pre/Post uses the existing pairwise analysis engine for each independent pair.")
    comparability_config = comparability_config or ComparabilityConfig()
    interpretation_config = interpretation_config or PairInterpretationConfig()
    interpretation_config.validate()
    run_dir, timestamp = _unique_run_dir(output_root, run_tag)
    figures_dir = run_dir / "figures"
    # Match the established Single Pair cache location. Keeping the 64-character
    # cache keys outside the deeper timestamped run directory avoids Windows
    # path-length failures without changing result-folder contents.
    cache_dir = Path(output_root).expanduser().resolve() / ".spatialtx_comparative_cache"
    figures_dir.mkdir(parents=True, exist_ok=True)
    pair_rows: list[dict] = []
    qc_tables: list[pd.DataFrame] = []
    context_gene_audit_rows: list[dict] = []
    figures: list[Path] = []
    input_hashes: dict[str, dict[str, str]] = {}

    for index, pair in enumerate(pairs, start=1):
        if cancel_event is not None and cancel_event.is_set():
            break
        if progress:
            progress(f"Pair {index}/{len(pairs)}: {pair.label} — technical QC")
        row = _pair_result_row(pair)
        pre_qc = post_qc = None
        structural_errors: list[str] = []
        if pair.pre_path == pair.post_path:
            structural_errors.append("Pre and Post refer to the same file")
        if pair.pre_path.suffix.casefold() != ".h5ad" or pair.post_path.suffix.casefold() != ".h5ad":
            structural_errors.append("both inputs must use the .h5ad extension")
        if not pair.pre_path.is_file():
            structural_errors.append(f"Pre H5AD does not exist: {pair.pre_path}")
        if not pair.post_path.is_file():
            structural_errors.append(f"Post H5AD does not exist: {pair.post_path}")
        if structural_errors:
            row["error"] = "; ".join(structural_errors)
        else:
            try:
                pre_qc = inspect_h5ad_qc(pair.pre_path, analysis_config.c_genes, analysis_config.s_genes)
            except Exception as exc:
                row["error"] = f"Pre unreadable or corrupted: {exc}"
            try:
                post_qc = inspect_h5ad_qc(pair.post_path, analysis_config.c_genes, analysis_config.s_genes)
            except Exception as exc:
                row["error"] = "; ".join(filter(None, [row["error"], f"Post unreadable or corrupted: {exc}"]))
        input_hashes[pair.label] = {
            "pre_sha256": str((pre_qc or {}).get("input_sha256", "")),
            "post_sha256": str((post_qc or {}).get("input_sha256", "")),
        }

        pre_metrics = post_metrics = pre_fields = post_fields = None
        if pre_qc is not None and post_qc is not None:
            try:
                if progress:
                    progress(f"Pair {index}/{len(pairs)}: {pair.label} — Pre analysis")
                pre_record = SampleRecord(f"{_safe_slug(pair.label)}__pre", pair.pre_path, "Pre", pair_id=pair.label, condition="Pre")
                post_record = SampleRecord(f"{_safe_slug(pair.label)}__post", pair.post_path, "Post", pair_id=pair.label, condition="Post")
                pre_metrics, pre_fields, _ = analyze_sample(pre_record, analysis_config, cache_dir)
                if progress:
                    progress(f"Pair {index}/{len(pairs)}: {pair.label} — Post analysis")
                post_metrics, post_fields, _ = analyze_sample(post_record, analysis_config, cache_dir)
                pre_metrics, post_metrics = _apply_within_pair_context_thresholds(
                    pre_metrics,
                    post_metrics,
                    pre_fields,
                    post_fields,
                )
                _add_metric_results(row, pre_metrics, post_metrics)
                _add_context_audit_results(row, pre_metrics, post_metrics)
                row.update({
                    "status": "PASS",
                    "error": "",
                    "regime_pre": str(pre_metrics.get("regime_label", "")),
                    "regime_post": str(post_metrics.get("regime_label", "")),
                    "pattern_pre": str(pre_metrics.get("public_transition_pattern", "")),
                    "pattern_post": str(post_metrics.get("public_transition_pattern", "")),
                })
                row["regime_transition"] = f"{row['regime_pre']} → {row['regime_post']}"
                row["pattern_transition"] = f"{row['pattern_pre'] or '-'} → {row['pattern_post'] or '-'}"
            except Exception as exc:
                row["error"] = f"Analysis exception: {exc}"

        classification, reasons, qc_table = evaluate_comparability(
            pair.label,
            pre_qc,
            post_qc,
            comparability_config,
            pre_metrics,
            post_metrics,
            pre_fields,
            post_fields,
        )
        pair_id_status = str(row.get("pair_id_validation", "not_available"))
        pair_id_reason = str(row.get("pair_id_warning", ""))
        pair_id_qc = _comparison_row(
            pair.label,
            "filename_pair_id_check",
            "pair_validation",
            row.get("pre_inferred_pair_id", ""),
            row.get("post_inferred_pair_id", ""),
            1 if pair_id_status == "matched" else 0 if pair_id_status == "warning" else np.nan,
            "warning" if pair_id_status == "warning" else "available" if pair_id_status == "matched" else "not_available",
            "caution" if pair_id_status == "warning" else "ok" if pair_id_status == "matched" else "not_available",
            pair_id_reason or ("Filename-based pair ID was not confidently inferable" if pair_id_status == "not_available" else ""),
            False,
        )
        site_status = normalize_site_comparability(pair.site_comparability)
        site_qc = _comparison_row(
            pair.label,
            "anatomical_site_comparability",
            "site_metadata",
            site_status,
            site_status,
            0 if site_status == "different_site" else 1 if site_status == "same_site" else np.nan,
            "warning" if site_status == "different_site" else "available" if site_status == "same_site" else "not_available",
            "caution" if site_status == "different_site" else "ok" if site_status == "same_site" else "not_available",
            site_warning_text(site_status),
            False,
        )
        qc_table = pd.concat((qc_table, pd.DataFrame([pair_id_qc, site_qc])), ignore_index=True)
        row["comparability"] = classification
        row["comparability_reasons"] = "; ".join(reasons)
        qc_table.insert(1, "pre_file", str(pair.pre_path))
        qc_table.insert(2, "post_file", str(pair.post_path))
        qc_table["comparability"] = classification
        qc_table["classification_reasons"] = "; ".join(reasons)
        qc_tables.append(qc_table)
        for sample_role, sample_file, sample_metrics in (
            ("Pre", pair.pre_path, pre_metrics),
            ("Post", pair.post_path, post_metrics),
        ):
            for axis in ("H", "V"):
                context_gene_audit_rows.append(_context_gene_audit_row(
                    pair,
                    sample_role,
                    sample_file,
                    axis,
                    sample_metrics,
                    analysis_config,
                ))
        pair_rows.append(row)

    pair_results = pd.DataFrame(pair_rows)
    context_gene_audit = pd.DataFrame(context_gene_audit_rows)
    comparability_qc = pd.concat(qc_tables, ignore_index=True) if qc_tables else pd.DataFrame()
    comparability_details = _comparability_details_table(pair_results, comparability_qc)
    pair_interpretation_summary = _pair_interpretation_table(
        pair_results,
        comparability_details,
        interpretation_config,
    )
    overview = pd.DataFrame([_overview_row(row) for row in pair_rows])
    interpretation_columns = [
        "pair_label", "balance_change_class", "balance_change_basis", "spatial_change_class",
        "spatial_change_basis", "regime_preserved", "structure_preserved",
        "structure_preservation_note", "pair_status_flag", "interpretive_flag",
        "technical_mismatch_reason", "sampling_mismatch_reason", "composition_proxy_reason",
        "pair_id_validation", "pair_id_warning", "caution_message", "note_message",
        "interpretation_confidence", "site_shift_warning",
        "qc_aware_interpretation_note",
    ]
    overview = overview.merge(
        pair_interpretation_summary[interpretation_columns],
        on="pair_label",
        how="left",
        validate="one_to_one",
    )
    overview_interpretation = overview[[
        "pair_label", "comparability", "balance_change_class", "spatial_change_class",
        "regime_transition", "regime_preserved", "structure_preserved",
        "delta_C", "delta_S", "delta_R", "delta_interface_fraction",
        "delta_diffuse_fraction", "delta_transition_burden", "delta_H", "delta_V",
        "interpretive_flag", "interpretation_confidence", "site_comparability",
        "site_shift_warning", "qc_aware_interpretation_note",
        "pair_status_flag", "technical_mismatch_reason", "sampling_mismatch_reason",
        "composition_proxy_reason", "pair_id_validation", "pair_id_warning",
    ]].copy()
    overview_interpretation = overview_interpretation.rename(columns={
        "delta_interface_fraction": "delta_interface",
        "delta_diffuse_fraction": "delta_diffuse",
        "delta_transition_burden": "delta_burden",
    })
    balance_changes = _change_layer_table(pair_results, "Balance change", BALANCE_METRIC_NAMES)
    spatial_organization_changes = _change_layer_table(
        pair_results,
        "Spatial organization change",
        SPATIAL_ORGANIZATION_METRIC_NAMES,
    )
    context_changes = _change_layer_table(pair_results, "Context axes", CONTEXT_METRIC_NAMES)
    context_audit_columns = (
        "pre_H_gene_coverage", "post_H_gene_coverage", "pre_H_context_status", "post_H_context_status",
        "pre_H_raw_normalization_method", "post_H_raw_normalization_method", "H_pair_pooled_q90",
        "pre_V_gene_coverage", "post_V_gene_coverage", "pre_V_context_status", "post_V_context_status",
        "pre_V_raw_normalization_method", "post_V_raw_normalization_method", "V_pair_pooled_q90",
    )
    for column in context_audit_columns:
        context_changes[column] = pair_results.get(column, pd.Series(index=pair_results.index, dtype=object))
    specimen_reliability = _specimen_reliability_table(pair_results, comparability_qc)
    comparative_qc_summary = build_comparative_qc_summary(
        pair_results,
        comparability_qc,
        pair_interpretation_summary,
    )
    multiaxial_pair_summary = build_multiaxial_pair_summary(
        pair_results,
        pair_interpretation_summary,
    )
    cohort_summary = _cohort_summary(overview)
    _write_csv(run_dir / "pair_results.csv", pair_results)
    _write_csv(run_dir / "balance_changes.csv", balance_changes)
    _write_csv(run_dir / "spatial_organization_changes.csv", spatial_organization_changes)
    _write_csv(run_dir / "context_changes.csv", context_changes)
    _write_csv(run_dir / "context_gene_audit.csv", context_gene_audit)
    _write_csv(run_dir / "specimen_reliability.csv", specimen_reliability)
    _write_csv(run_dir / "pair_interpretation_summary.csv", pair_interpretation_summary)
    _write_csv(run_dir / "comparability_details.csv", comparability_details)
    _write_csv(run_dir / "overview_interpretation.csv", overview_interpretation)
    _write_csv(run_dir / "comparative_overview.csv", overview)
    _write_csv(run_dir / "comparability_qc.csv", comparability_qc)
    _write_csv(run_dir / "comparative_qc_summary.csv", comparative_qc_summary)
    _write_csv(run_dir / "multiaxial_pair_summary.csv", multiaxial_pair_summary)
    _write_csv(run_dir / "cohort_summary.csv", cohort_summary)
    interpretation_by_pair = {
        str(row["pair_label"]): row.to_dict() for _, row in pair_interpretation_summary.iterrows()
    }
    for _, pair_row in pair_results.loc[pair_results["status"].eq("PASS")].iterrows():
        figure_row = pair_row.to_dict()
        figure_row.update(interpretation_by_pair.get(str(pair_row["pair_label"]), {}))
        figures.append(
            _plot_pair_panel(
                figure_row,
                figures_dir / f"pair_{_safe_slug(str(pair_row['pair_label']))}_metrics.png",
            )
        )
    figures.insert(0, _plot_overview(overview, figures_dir / "multi_pair_comparative_overview.png"))
    figures.insert(1, plot_multiaxial_overview(
        multiaxial_pair_summary,
        figures_dir / "multiaxial_pair_overview.png",
    ))

    metadata = {
        "application": "SpatialTX Studio Desktop",
        "spatialtx_version": __version__,
        "analysis_module": "Multi-Pair Pre/Post Comparative Analysis",
        "created_utc": timestamp,
        "pair_count_requested": len(pairs),
        "pair_count_completed": int(pair_results.get("status", pd.Series(dtype=str)).eq("PASS").sum()),
        "pair_count_failed": int(pair_results.get("status", pd.Series(dtype=str)).eq("ERROR").sum()),
        "analysis_parameters": analysis_config.to_dict(),
        "comparative_metric_layer_schema": COMPARATIVE_METRIC_LAYER_SCHEMA,
        "effective_context_programs": {
            axis: {
                **resolve_context_program(
                    axis,
                    analysis_config.h_genes if axis == "H" else analysis_config.v_genes,
                ),
                "enabled": analysis_config.enable_h_expr if axis == "H" else analysis_config.enable_v_expr,
                "minimum_coverage": float(analysis_config.context_min_coverage),
            }
            for axis in ("H", "V")
        },
        "raw_context_summary_method": (
            "Non-centered matched-gene program mean per spot; raw counts use log1p before the program mean, "
            "while supported nonnegative log1p expression is used as stored."
        ),
        "pair_pooled_high_quantile": POOLED_HIGH_QUANTILE,
        "pair_pooled_threshold_scope": "within_pair_pre_plus_post",
        "context_warning_provenance": {
            "single_sample_context_warning": "legacy_within_sample_centered_context_q80",
            "pair_pooled_context_warning": "within_pair_pre_plus_post_raw_context_q90",
            "interpretation": (
                "A single-sample centered-context warning does not describe the pair-pooled "
                "raw-context high fraction."
            ),
        },
        "direction_tolerances": {metric.export_name: metric.direction_tolerance for metric in MULTI_PAIR_METRICS},
        "comparability_configuration": comparability_config.to_dict(),
        "pair_interpretation_configuration": interpretation_config.to_dict(),
        "selected_pairs": [
            {
                "label": pair.label,
                "pre_file": str(pair.pre_path),
                "post_file": str(pair.post_path),
                "pre_sha256": input_hashes.get(pair.label, {}).get("pre_sha256", ""),
                "post_sha256": input_hashes.get(pair.label, {}).get("post_sha256", ""),
                "site_comparability": pair.site_comparability,
            }
            for pair in pairs
        ],
        "metric_registry": [asdict(metric) for metric in MULTI_PAIR_METRICS],
        "result_layers": {
            "layer_1": "Balance change: C, S, and R = C - S",
            "layer_2": "Spatial organization change: interface, diffuse, adjacency, fragmentation, and topology metrics",
            "layer_3": "Specimen reliability: comparability and technical/sampling QC",
            "parallel_context_axes": (
                "H hypoxia-associated expression context and V endothelial/angiogenic expression proxy; "
                "observational only and not combined with C/S or FRAME2.6."
            ),
            "combination_rule": "The three layers remain separate and are not combined into a response or quality score.",
        },
        "pair_id_validation": (
            "Filename-based patient/sample IDs are checked conservatively. A mismatch is a visible warning, "
            "not an automatic run blocker or a comparability-class override."
        ),
        "delta_definition": "Post - Pre",
        "comparability_role": "Separate specimen technical/sampling context; not confidence in biological truth.",
        "spotwise_subtraction_performed": False,
        "spatial_registration_performed": False,
        "clinical_interpretation_performed": False,
        "composite_response_score_computed": False,
        "H_V_modify_core_fields_or_regimes": False,
        "site_comparability_values": list(SITE_COMPARABILITY_VALUES),
        "notice": MULTI_PAIR_NOTICE,
    }
    (run_dir / "run_metadata.json").write_text(
        json.dumps(json_safe(metadata), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    passed = int(pair_results.get("status", pd.Series(dtype=str)).eq("PASS").sum())
    failed = int(len(pair_results) - passed)
    summary_text = (
        f"Multi-Pair Pre/Post complete: {passed}/{len(pair_results)} pairs passed; {failed} failed.\n"
        "Layer 1 — Balance change: C, S, and R = C - S are reported separately as Pre, Post, and Delta.\n"
        "Layer 2 — Spatial organization change: interface, diffuse, adjacency, fragmentation, and topology metrics are reported separately.\n"
        f"Layer 3 — Specimen reliability: Good {int(overview['comparability'].eq('Good').sum())}, "
        f"Caution {int(overview['comparability'].eq('Caution').sum())}, "
        f"Low {int(overview['comparability'].eq('Low').sum())}. "
        f"Pair-ID warnings: {int(overview['pair_id_validation'].eq('warning').sum())}.\n"
        f"Balance classes — Minimal {int(overview['balance_change_class'].eq('Minimal').sum())}, "
        f"Moderate {int(overview['balance_change_class'].eq('Moderate').sum())}, "
        f"Large {int(overview['balance_change_class'].eq('Large').sum())}. "
        f"Spatial classes — Minimal {int(overview['spatial_change_class'].eq('Minimal').sum())}, "
        f"Moderate {int(overview['spatial_change_class'].eq('Moderate').sum())}, "
        f"Large {int(overview['spatial_change_class'].eq('Large').sum())}.\n"
        f"H available {int(multiaxial_pair_summary['H_pre'].notna().sum())}/{len(multiaxial_pair_summary)}; "
        f"V available {int(multiaxial_pair_summary['V_pre'].notna().sum())}/{len(multiaxial_pair_summary)}. "
        f"Site-shift warnings {int(multiaxial_pair_summary['site_comparability'].eq('different_site').sum())}.\n"
        f"The three layers and H/V context axes are not combined into a single response score. {MULTI_PAIR_NOTICE}"
    )
    return MultiPairRunResult(
        run_dir=run_dir,
        pair_results=pair_results,
        comparative_overview=overview,
        comparability_qc=comparability_qc,
        cohort_summary=cohort_summary,
        figures=figures,
        summary_text=summary_text,
        balance_changes=balance_changes,
        spatial_organization_changes=spatial_organization_changes,
        specimen_reliability=specimen_reliability,
        pair_interpretation_summary=pair_interpretation_summary,
        comparability_details=comparability_details,
        overview_interpretation=overview_interpretation,
        context_changes=context_changes,
        context_gene_audit=context_gene_audit,
        multiaxial_pair_summary=multiaxial_pair_summary,
        comparative_qc_summary=comparative_qc_summary,
    )
