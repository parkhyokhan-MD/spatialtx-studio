from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from ..advanced_analysis import _bh_adjust
from .metrics import GROUP_METRICS
from .models import ComparativeConfig


def _descriptive(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {"n": 0, "mean": np.nan, "median": np.nan, "std": np.nan, "iqr": np.nan}
    q25, q75 = np.quantile(values, [0.25, 0.75])
    return {
        "n": int(len(values)),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
        "iqr": float(q75 - q25),
    }


def _bootstrap_ci(
    reference: np.ndarray,
    target: np.ndarray,
    *,
    paired: bool,
    iterations: int,
    seed: int,
) -> tuple[float, float]:
    reference = np.asarray(reference, dtype=float)
    target = np.asarray(target, dtype=float)
    rng = np.random.default_rng(seed)
    estimates = np.empty(iterations, dtype=float)
    if paired:
        differences = target - reference
        for index in range(iterations):
            sample = rng.choice(differences, size=len(differences), replace=True)
            estimates[index] = np.mean(sample)
    else:
        for index in range(iterations):
            sampled_reference = rng.choice(reference, size=len(reference), replace=True)
            sampled_target = rng.choice(target, size=len(target), replace=True)
            estimates[index] = np.mean(sampled_target) - np.mean(sampled_reference)
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def _paired_effect(reference: np.ndarray, target: np.ndarray) -> tuple[float, str]:
    differences = np.asarray(target, dtype=float) - np.asarray(reference, dtype=float)
    sd = float(np.std(differences, ddof=1)) if len(differences) > 1 else np.nan
    effect = float(np.mean(differences) / sd) if np.isfinite(sd) and sd > 0 else np.nan
    return effect, "paired_standardized_mean_difference_dz"


def _hedges_g(reference: np.ndarray, target: np.ndarray) -> tuple[float, str]:
    n_ref, n_tar = len(reference), len(target)
    if n_ref < 2 or n_tar < 2:
        return np.nan, "hedges_g_target_minus_reference"
    pooled_n = n_ref + n_tar - 2
    pooled_variance = ((n_ref - 1) * np.var(reference, ddof=1) + (n_tar - 1) * np.var(target, ddof=1)) / pooled_n
    if pooled_variance <= 0:
        return np.nan, "hedges_g_target_minus_reference"
    d = (np.mean(target) - np.mean(reference)) / np.sqrt(pooled_variance)
    correction = 1 - 3 / max(4 * (n_ref + n_tar) - 9, 1)
    return float(d * correction), "hedges_g_target_minus_reference"


def _unpaired_rank_effect(reference: np.ndarray, target: np.ndarray, u_target: float) -> tuple[float, str]:
    effect = 2 * float(u_target) / (len(reference) * len(target)) - 1
    return float(effect), "rank_biserial_target_minus_reference"


def _test_metric(
    reference: np.ndarray,
    target: np.ndarray,
    *,
    paired: bool,
    selected_test: str,
) -> tuple[str, float, float, float, str]:
    if paired:
        test_name = "paired_t" if selected_test == "paired_t" else "wilcoxon_signed_rank"
        effect, effect_name = _paired_effect(reference, target)
        if test_name == "paired_t":
            result = stats.ttest_rel(target, reference, nan_policy="omit")
            return test_name, float(result.statistic), float(result.pvalue), effect, effect_name
        differences = target - reference
        if np.allclose(differences, 0):
            return test_name, 0.0, 1.0, effect, effect_name
        result = stats.wilcoxon(target, reference, alternative="two-sided", zero_method="wilcox")
        return test_name, float(result.statistic), float(result.pvalue), effect, effect_name
    test_name = "welch_t" if selected_test == "welch_t" else "mann_whitney_u"
    if test_name == "welch_t":
        result = stats.ttest_ind(target, reference, equal_var=False, nan_policy="omit")
        effect, effect_name = _hedges_g(reference, target)
        return test_name, float(result.statistic), float(result.pvalue), effect, effect_name
    result = stats.mannwhitneyu(target, reference, alternative="two-sided", method="auto")
    effect, effect_name = _unpaired_rank_effect(reference, target, float(result.statistic))
    return test_name, float(result.statistic), float(result.pvalue), effect, effect_name


def comparative_group_statistics(
    sample_metrics: pd.DataFrame,
    config: ComparativeConfig,
    matches: pd.DataFrame | None,
    *,
    effective_mode: str,
) -> pd.DataFrame:
    paired = effective_mode == "paired"
    selected_test = config.statistical_test
    if selected_test == "auto":
        selected_test = "wilcoxon" if paired else "mannwhitney"
    allowed = {"wilcoxon", "paired_t"} if paired else {"mannwhitney", "welch_t"}
    if selected_test not in allowed:
        raise ValueError(f"Statistical test {selected_test!r} is not valid for {effective_mode} mode.")
    by_sample = sample_metrics.set_index("sample_id")
    rows: list[dict] = []
    for metric_index, metric in enumerate(GROUP_METRICS):
        if metric not in sample_metrics:
            continue
        if paired:
            assert matches is not None
            reference = pd.to_numeric(
                pd.Series([by_sample.loc[str(row["reference_sample_id"])].get(metric) for _, row in matches.iterrows()]),
                errors="coerce",
            ).to_numpy(dtype=float)
            target = pd.to_numeric(
                pd.Series([by_sample.loc[str(row["target_sample_id"])].get(metric) for _, row in matches.iterrows()]),
                errors="coerce",
            ).to_numpy(dtype=float)
            finite = np.isfinite(reference) & np.isfinite(target)
            reference, target = reference[finite], target[finite]
        else:
            reference = pd.to_numeric(
                sample_metrics.loc[sample_metrics["group"].eq(config.reference), metric], errors="coerce"
            ).dropna().to_numpy(dtype=float)
            target = pd.to_numeric(
                sample_metrics.loc[sample_metrics["group"].eq(config.target), metric], errors="coerce"
            ).dropna().to_numpy(dtype=float)
        ref_desc, tar_desc = _descriptive(reference), _descriptive(target)
        minimum = 2 if paired else 1
        status = "ok" if len(reference) >= minimum and len(target) >= minimum else "insufficient_values"
        test_name = statistic = p_value = effect_size = np.nan
        effect_name = ""
        ci_low = ci_high = np.nan
        if status == "ok":
            try:
                test_name, statistic, p_value, effect_size, effect_name = _test_metric(
                    reference, target, paired=paired, selected_test=selected_test
                )
                ci_low, ci_high = _bootstrap_ci(
                    reference,
                    target,
                    paired=paired,
                    iterations=int(config.bootstrap_iterations),
                    seed=int(config.seed) + metric_index,
                )
            except (ValueError, FloatingPointError) as exc:
                status = f"test_unavailable: {exc}"
        rows.append({
            "metric": metric,
            "comparison_design": effective_mode,
            "reference_group": config.reference,
            "target_group": config.target,
            "n_reference": ref_desc["n"],
            "reference_mean": ref_desc["mean"],
            "reference_median": ref_desc["median"],
            "reference_std": ref_desc["std"],
            "reference_iqr": ref_desc["iqr"],
            "n_target": tar_desc["n"],
            "target_mean": tar_desc["mean"],
            "target_median": tar_desc["median"],
            "target_std": tar_desc["std"],
            "target_iqr": tar_desc["iqr"],
            "mean_difference_target_minus_reference": (
                tar_desc["mean"] - ref_desc["mean"]
                if np.isfinite(tar_desc["mean"]) and np.isfinite(ref_desc["mean"])
                else np.nan
            ),
            "effect_size": effect_size,
            "effect_size_method": effect_name,
            "confidence_interval_method": "95% bootstrap mean-difference CI",
            "ci_95_low": ci_low,
            "ci_95_high": ci_high,
            "test": test_name,
            "test_statistic": statistic,
            "p_value": p_value,
            "status": status,
            "interpretation": "Exploratory statistical comparison; significance alone does not establish biological importance.",
        })
    table = pd.DataFrame(rows)
    if len(table):
        table["adjusted_p_value_bh"] = _bh_adjust(pd.to_numeric(table["p_value"], errors="coerce").to_numpy())
        table["significant_fdr_0_05"] = table["adjusted_p_value_bh"] < 0.05
        table["fdr_scope"] = "across available comparative metrics within this run"
    return table
