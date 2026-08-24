from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from spatialtx_desktop import __version__
from spatialtx_desktop.comparative.models import ComparativeConfig
from spatialtx_desktop.comparative.multi_pair import PairSpec, run_multi_pair_analysis
from spatialtx_desktop.graph.metadata import json_safe
from spatialtx_desktop.reliability.models import ReliabilityConfig


C_GENES = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG"]
S_GENES = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN"]
PAIR_IDS = (30, 38, 41, 42, 43, 44)
GOLDEN_DELTA_B = {
    "sample_30": 0.0663950512784694,
    "sample_38": 0.037868473546208314,
    "sample_41": -0.030295522043840173,
    "sample_42": 0.11513113040116349,
    "sample_43": -0.22667093441959016,
    "sample_44": -0.06837231767997905,
}
BALANCE_REGRESSION_FIELDS = (
    "pre_B",
    "post_B",
    "delta_B",
    "delta_B_bootstrap_ci_low",
    "delta_B_bootstrap_ci_high",
    "delta_B_permutation_p_value",
    "delta_B_bh_fdr",
    "pre_C_gene_coverage",
    "post_C_gene_coverage",
    "pre_S_gene_coverage",
    "post_S_gene_coverage",
)


def _one_match(data_root: Path, pattern: str) -> Path:
    matches = sorted(data_root.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one file for {pattern!r}; found {len(matches)}")
    return matches[0]


def _write_reliability_regression_artifacts(
    run_dir: Path,
    summary: pd.DataFrame,
    baseline_path: Path | None,
) -> dict:
    fixed_summary = run_dir / "reliability_pair_summary_fixed.csv"
    summary.to_csv(fixed_summary, index=False, encoding="utf-8-sig", na_rep="NA")

    expected_order = list(GOLDEN_DELTA_B)
    actual_order = summary["pair_label"].astype(str).tolist()
    order_pass = actual_order == expected_order
    by_pair = summary.set_index("pair_label")
    diff_rows: list[dict] = []
    for pair_label, reference in GOLDEN_DELTA_B.items():
        actual = float(by_pair.loc[pair_label, "delta_B"]) if pair_label in by_pair.index else np.nan
        passed = bool(np.isclose(actual, reference, rtol=1.0e-12, atol=1.0e-12, equal_nan=True))
        diff_rows.append({
            "comparison": "golden_delta_B",
            "pair_label": pair_label,
            "field": "delta_B",
            "baseline_value": reference,
            "new_value": actual,
            "absolute_difference": abs(actual - reference) if np.isfinite(actual) else np.nan,
            "rtol": 1.0e-12,
            "atol": 1.0e-12,
            "passed": passed,
        })

    baseline_used = baseline_path is not None
    baseline_order_pass = True
    if baseline_path is not None:
        baseline = pd.read_csv(baseline_path)
        baseline_order = baseline["pair_label"].astype(str).tolist()
        baseline_order_pass = baseline_order == actual_order
        baseline_by_pair = baseline.set_index("pair_label")
        for pair_label in actual_order:
            for field in BALANCE_REGRESSION_FIELDS:
                old = pd.to_numeric(
                    pd.Series([baseline_by_pair.loc[pair_label, field]]), errors="coerce"
                ).iloc[0]
                new = pd.to_numeric(pd.Series([by_pair.loc[pair_label, field]]), errors="coerce").iloc[0]
                passed = bool(np.isclose(old, new, rtol=1.0e-12, atol=1.0e-12, equal_nan=True))
                diff_rows.append({
                    "comparison": "preserved_balance_baseline",
                    "pair_label": pair_label,
                    "field": field,
                    "baseline_value": old,
                    "new_value": new,
                    "absolute_difference": abs(new - old) if np.isfinite(old) and np.isfinite(new) else np.nan,
                    "rtol": 1.0e-12,
                    "atol": 1.0e-12,
                    "passed": passed,
                })

    diff = pd.DataFrame(diff_rows)
    diff_path = run_dir / "reliability_regression_diff.csv"
    diff.to_csv(diff_path, index=False, encoding="utf-8-sig", na_rep="NA")
    numeric_pass = bool(diff["passed"].all()) if not diff.empty else False
    required_files = {
        "cross_exclusivity_audit.csv": (run_dir / "cross_exclusivity_audit.csv").is_file(),
        "reliability_gene_coverage.csv": (run_dir / "reliability_gene_coverage.csv").is_file(),
        "run_metadata.json": (run_dir / "run_metadata.json").is_file(),
    }
    report = {
        "status": "PASS" if numeric_pass and order_pass and baseline_order_pass and all(required_files.values()) else "FAIL",
        "rtol": 1.0e-12,
        "atol": 1.0e-12,
        "expected_pair_order": expected_order,
        "actual_pair_order": actual_order,
        "pair_order_pass": order_pass,
        "baseline_pair_order_pass": baseline_order_pass,
        "baseline_summary": str(baseline_path) if baseline_path is not None else None,
        "full_baseline_comparison_performed": baseline_used,
        "balance_regression_fields": list(BALANCE_REGRESSION_FIELDS),
        "numeric_regression_pass": numeric_pass,
        "required_output_files": required_files,
        "failed_comparison_count": int((~diff["passed"]).sum()),
        "comparison_row_count": int(len(diff)),
    }
    report_path = run_dir / "reliability_regression_report.json"
    report_path.write_text(
        json.dumps(json_safe(report), indent=2, ensure_ascii=False), encoding="utf-8"
    )

    qc_path = run_dir / "reliability_qc.json"
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    qc["six_pair_balance_regression"] = report
    qc_fixed_path = run_dir / "reliability_qc_fixed.json"
    qc_fixed_path.write_text(
        json.dumps(json_safe(qc), indent=2, ensure_ascii=False), encoding="utf-8"
    )

    metadata_path = run_dir / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["six_pair_reliability_regression"] = report
    output_files = metadata.get("reliability_layer", {}).get("output_files", [])
    for name in (
        fixed_summary.name,
        diff_path.name,
        report_path.name,
        qc_fixed_path.name,
    ):
        if name not in output_files:
            output_files.append(name)
    metadata["reliability_layer"]["output_files"] = output_files
    metadata_path.write_text(
        json.dumps(json_safe(metadata), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the preserved six-pair BTC regression set.")
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--context", action="store_true", help="Enable the existing H/V context programs.")
    parser.add_argument(
        "--reliability",
        action="store_true",
        help="Enable the additive v0.65 Reliability Layer sidecars (continuous mode).",
    )
    parser.add_argument(
        "--baseline-reliability-summary",
        type=Path,
        help="Optional preserved v0.65 reliability_pair_summary.csv for full Balance regression.",
    )
    args = parser.parse_args()

    data_root = args.data_root.expanduser().resolve()
    pairs = [
        PairSpec(
            f"sample_{pair_id}",
            _one_match(data_root, f"*_sample_{pair_id}_pre.h5ad"),
            _one_match(data_root, f"*_sample_{pair_id}_post.h5ad"),
        )
        for pair_id in PAIR_IDS
    ]
    config = ComparativeConfig(
        mode="pairwise",
        reference="Pre",
        target="Post",
        c_genes=C_GENES,
        s_genes=S_GENES,
        c_q=0.8,
        s_q=0.8,
        g_q=0.6,
        scoring_options={
            "smoothing_mode": "none",
            "smoothing_k": 6,
            "gaussian_sigma": 0.0,
            "normalization_mode": "raw_mean",
            "perturbation_check": False,
            "c_q_list": [0.75, 0.80, 0.85],
            "s_q_list": [0.75, 0.80, 0.85],
            "g_q_list": [0.50, 0.60, 0.70],
            "parameter_log_export": True,
            "dense_warning_gb": 4.0,
        },
        graph_settings={"method": "knn", "k": 6, "weighting": "binary"},
        enable_h_expr=args.context,
        enable_v_expr=args.context,
        context_smoothing="none",
        seed=42,
        bootstrap_iterations=1000,
        use_cache=False,
    )
    result = run_multi_pair_analysis(
        pairs,
        config,
        args.output_root,
        run_tag=args.run_tag,
        progress=print,
        reliability_config=ReliabilityConfig(
            enabled=args.reliability,
            score_domain="nonnegative",
            epsilon=1.0e-9,
            classification_enabled=False,
            strict_cross_exclusivity=True,
            dependence_qc=True,
            bootstrap_iterations=1000,
            permutation_iterations=1000,
            fdr_method="benjamini-hochberg",
            seed=42,
            minimum_valid_spots=30,
            minimum_valid_fraction=0.80,
            warning_valid_fraction=0.50,
        ),
    )
    if args.reliability:
        baseline = (
            args.baseline_reliability_summary.expanduser().resolve()
            if args.baseline_reliability_summary is not None
            else None
        )
        report = _write_reliability_regression_artifacts(
            result.run_dir,
            result.reliability_pair_summary,
            baseline,
        )
        print(f"Six-pair Balance regression: {report['status']}")
        if report["status"] != "PASS":
            raise RuntimeError(
                "Six-pair Balance regression failed; the fixed reliability result is not approved."
            )
    print(f"SpatialTX {__version__}: {result.run_dir}")
    print(result.pair_results[[
        "pair_label", "regime_pre", "regime_post", "comparability",
        "pre_interface_fraction", "post_interface_fraction",
        "pre_diffuse_fraction", "post_diffuse_fraction",
        "pre_transition_burden", "post_transition_burden",
    ]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
