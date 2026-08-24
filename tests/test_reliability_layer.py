from __future__ import annotations

import unittest
import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from spatialtx_desktop.reliability.audit import (
    CrossExclusivityError,
    audit_cross_exclusivity,
    gene_coverage_audit,
)
from spatialtx_desktop.reliability.core import compute_axis_reliability, compute_reliability_axes
from spatialtx_desktop.reliability.dependence import bh_adjust, compute_axis_dependence
from spatialtx_desktop.reliability.exports import build_metric_qc_table, build_pair_summary
from spatialtx_desktop.reliability.models import ReliabilityConfig
from spatialtx_desktop.comparative.models import ComparativeConfig
from spatialtx_desktop.comparative.multi_pair import PairSpec, run_multi_pair_analysis
from tests.comparative_fixtures import C_GENES, S_GENES, write_comparative_h5ad


class ReliabilityCoreTests(unittest.TestCase):
    def test_zero_zero_is_explicitly_undefined_for_direction_and_fraction(self) -> None:
        result = compute_axis_reliability([0.0], [0.0])
        self.assertEqual(result.activity_A[0], 0.0)
        self.assertEqual(result.balance_B[0], 0.0)
        self.assertEqual(result.ca_strength[0], 0.0)
        self.assertTrue(np.isnan(result.direction_D[0]))
        self.assertTrue(np.isnan(result.ca_fraction[0]))
        self.assertEqual(result.status[0], "undefined_direction_zero_activity")

    def test_equal_positive_scores_are_fully_coactive(self) -> None:
        result = compute_axis_reliability([2.0], [2.0], ReliabilityConfig(epsilon=1.0e-12))
        self.assertEqual(result.balance_B[0], 0.0)
        self.assertAlmostEqual(result.direction_D[0], 0.0)
        self.assertEqual(result.ca_strength[0], result.activity_A[0])
        self.assertAlmostEqual(result.ca_fraction[0], 1.0)

    def test_single_pole_activity_has_zero_coactivation(self) -> None:
        result = compute_axis_reliability([3.0, 0.0], [0.0, 4.0], ReliabilityConfig(epsilon=1.0e-12))
        self.assertAlmostEqual(result.direction_D[0], 1.0, places=10)
        self.assertAlmostEqual(result.direction_D[1], -1.0, places=10)
        np.testing.assert_allclose(result.ca_strength, [0.0, 0.0])

    def test_swapping_poles_preserves_activity_and_coactivation(self) -> None:
        C = np.array([0.1, 0.5, 3.0, 8.0])
        S = np.array([0.9, 0.5, 2.0, 1.0])
        forward = compute_axis_reliability(C, S)
        reverse = compute_axis_reliability(S, C)
        np.testing.assert_allclose(forward.activity_A, reverse.activity_A)
        np.testing.assert_allclose(forward.ca_strength, reverse.ca_strength)
        np.testing.assert_allclose(forward.ca_fraction, reverse.ca_fraction)
        np.testing.assert_allclose(forward.balance_B, -reverse.balance_B)
        np.testing.assert_allclose(forward.direction_D, -reverse.direction_D)

    def test_nonnegative_property_identities_hold(self) -> None:
        rng = np.random.default_rng(42)
        C = rng.uniform(0.001, 20.0, size=500)
        S = rng.uniform(0.001, 20.0, size=500)
        result = compute_axis_reliability(C, S, ReliabilityConfig(epsilon=1.0e-12))
        self.assertTrue(np.all(result.activity_A >= np.abs(result.balance_B)))
        np.testing.assert_allclose(
            result.ca_strength,
            result.activity_A - np.abs(result.balance_B),
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(result.ca_strength, 2.0 * np.minimum(C, S))
        self.assertTrue(np.nanmin(result.direction_D) >= -1.0)
        self.assertTrue(np.nanmax(result.direction_D) <= 1.0)

    def test_nan_inf_and_negative_inputs_are_not_replaced(self) -> None:
        result = compute_axis_reliability(
            [np.nan, np.inf, -1.0, 2.0],
            [1.0, 1.0, 2.0, -3.0],
        )
        self.assertEqual(result.status.tolist(), [
            "invalid_nonfinite_input",
            "invalid_nonfinite_input",
            "invalid_negative_score",
            "invalid_negative_score",
        ])
        self.assertTrue(np.isnan(result.activity_A).all())
        self.assertTrue(np.isnan(result.direction_D).all())
        self.assertTrue(np.isnan(result.ca_strength).all())
        self.assertTrue(np.isnan(result.ca_fraction).all())
        self.assertTrue(np.isnan(result.balance_B[0]))
        self.assertTrue(np.isnan(result.balance_B[1]))
        self.assertEqual(result.balance_B[2], -3.0)
        self.assertEqual(result.balance_B[3], 5.0)

    def test_signed_legacy_balance_is_separate_from_nonnegative_activity(self) -> None:
        legacy_c = np.array([-2.0, 0.5, 3.0])
        legacy_s = np.array([1.0, -0.5, 2.0])
        activity_c = np.array([0.2, 0.5, 3.0])
        activity_s = np.array([1.0, 0.5, 2.0])
        result = compute_axis_reliability(
            legacy_c,
            legacy_s,
            activity_C=activity_c,
            activity_S=activity_s,
        )
        np.testing.assert_allclose(result.balance_B, legacy_c - legacy_s)
        np.testing.assert_allclose(result.activity_A, activity_c + activity_s)
        np.testing.assert_allclose(result.activity_balance, activity_c - activity_s)
        self.assertTrue(result.valid_input.all())
        self.assertFalse(np.array_equal(result.activity_C, result.C))

    def test_classification_requires_explicit_thresholds_and_uses_only_four_states(self) -> None:
        with self.assertRaisesRegex(ValueError, "requires explicit"):
            compute_axis_reliability([1.0], [1.0], ReliabilityConfig(classification_enabled=True))
        config = ReliabilityConfig(
            classification_enabled=True,
            activity_threshold=0.5,
            direction_threshold=0.5,
        )
        result = compute_axis_reliability(
            [0.0, 3.0, 0.1, 1.0],
            [0.0, 0.1, 3.0, 1.0],
            config,
        )
        self.assertEqual(result.reliability_state.tolist(), [
            "low_activity",
            "c_dominant_active",
            "s_dominant_active",
            "active_coactivation_candidate",
        ])

    def test_multiple_axes_use_the_same_core_without_transforming_dependence(self) -> None:
        first_C = np.array([1.0, 2.0, 3.0])
        first_S = np.array([0.5, 1.0, 1.5])
        second_C = first_C * 2.0
        second_S = first_S * 2.0
        results = compute_reliability_axes({
            "axis_1": (first_C, first_S),
            "axis_2": (second_C, second_S),
        })
        np.testing.assert_allclose(results["axis_1"].C, first_C)
        np.testing.assert_allclose(results["axis_2"].C, second_C)
        np.testing.assert_allclose(results["axis_2"].balance_B, 2.0 * results["axis_1"].balance_B)

    def test_pair_summary_retains_values_but_excludes_qc_fail_from_conclusion(self) -> None:
        pre = compute_axis_reliability([1.0, -1.0, -2.0], [1.0, 2.0, 3.0])
        post = compute_axis_reliability([2.0, -1.0, -2.0], [1.0, 2.0, 3.0])
        table = build_pair_summary(
            "pair",
            pre,
            post,
            ReliabilityConfig(
                minimum_valid_spots=3,
                bootstrap_iterations=5,
                permutation_iterations=5,
            ),
        )
        row = table.iloc[0]
        self.assertTrue(np.isfinite(row["delta_B"]))
        self.assertTrue(np.isfinite(row["delta_A"]))
        self.assertTrue(np.isnan(row["delta_A_permutation_p_value"]))
        self.assertEqual(
            row["pair_score_validity"],
            "qc_fail_insufficient_valid_spots_and_fraction",
        )
        self.assertFalse(row["activity_summary_included_in_conclusion"])

    def test_fraction_gate_uses_count_and_fraction(self) -> None:
        def result(valid_count: int, total: int = 100):
            activity_c = np.r_[np.ones(valid_count), -np.ones(total - valid_count)]
            activity_s = np.ones(total)
            return compute_axis_reliability(
                np.linspace(-1.0, 1.0, total),
                np.linspace(1.0, -1.0, total),
                activity_C=activity_c,
                activity_S=activity_s,
            )

        config = ReliabilityConfig(
            minimum_valid_spots=30,
            minimum_valid_fraction=0.80,
            warning_valid_fraction=0.50,
            bootstrap_iterations=5,
            permutation_iterations=5,
        )
        fail = build_pair_summary("fail", result(30), result(30), config).iloc[0]
        warning = build_pair_summary("warning", result(60), result(60), config).iloc[0]
        valid = build_pair_summary("valid", result(80), result(80), config).iloc[0]
        count_fail = build_pair_summary(
            "count_fail", result(19, total=20), result(19, total=20), config
        ).iloc[0]
        self.assertEqual(fail["pair_score_validity"], "qc_fail_insufficient_valid_fraction")
        self.assertEqual(warning["pair_score_validity"], "warning_low_valid_fraction")
        self.assertEqual(valid["pair_score_validity"], "valid")
        self.assertEqual(count_fail["pair_score_validity"], "qc_fail_insufficient_valid_spots")


class ReliabilityAuditTests(unittest.TestCase):
    def test_strict_cross_exclusivity_blocks_any_canonical_duplicate(self) -> None:
        programs = {
            "axis_1": {"C": [" VEGFA "], "S": ["COL1A1"]},
            "axis_2": {"C": ["vegfa"], "S": ["PECAM1"]},
        }
        with self.assertRaises(CrossExclusivityError) as caught:
            audit_cross_exclusivity(programs)
        row = caught.exception.audit.loc[
            caught.exception.audit["canonical_gene"].eq("VEGFA")
        ].iloc[0]
        self.assertEqual(row["overlap_type"], "cross_axis")
        self.assertEqual(row["severity"], "hard_error")
        self.assertEqual(row["action"], "analysis_blocked")

    def test_ensembl_versions_and_aliases_are_canonicalized(self) -> None:
        config = ReliabilityConfig(
            canonical_aliases={"OLD1": "GENE1"},
            canonicalization_source="test-map",
            canonicalization_version="1",
        )
        programs = {
            "axis_1": {"C": ["ENSG000001.4", "OLD1"], "S": ["GENE2"]},
            "axis_2": {"C": ["ENSG000001"], "S": ["GENE1"]},
        }
        audit = audit_cross_exclusivity(programs, config, raise_on_error=False)
        self.assertEqual(audit.loc[audit["canonical_gene"].eq("ENSG000001"), "severity"].item(), "hard_error")
        self.assertEqual(audit.loc[audit["canonical_gene"].eq("GENE1"), "severity"].item(), "hard_error")
        self.assertTrue(audit["normalization_rule"].str.contains("test-map:1", regex=False).all())

    def test_paralogs_are_not_inferred_as_duplicates(self) -> None:
        audit = audit_cross_exclusivity(
            {"axis": {"C": ["COL1A1"], "S": ["COL1A2"]}}
        )
        self.assertTrue(audit["overlap_type"].eq("none").all())

    def test_gene_coverage_reports_missing_and_validity_without_rescoring(self) -> None:
        table = gene_coverage_audit(
            {"axis": {"C": ["A", "B"], "S": ["C", "D"]}},
            ["A", "B", "C"],
        )
        c_row = table.loc[table["pole"].eq("C")].iloc[0]
        s_row = table.loc[table["pole"].eq("S")].iloc[0]
        self.assertEqual(c_row["gene_coverage_fraction"], 1.0)
        self.assertEqual(c_row["score_validity"], "valid")
        self.assertEqual(s_row["missing_genes"], "D")
        self.assertEqual(s_row["coverage_status"], "caution")


class ReliabilityMetricQCGateTests(unittest.TestCase):
    @staticmethod
    def _result(defined_n: int, total: int = 100):
        activity_c = np.r_[np.ones(defined_n), np.zeros(total - defined_n)]
        activity_s = np.zeros(total)
        return compute_axis_reliability(
            np.linspace(-1.0, 1.0, total),
            np.linspace(1.0, -1.0, total),
            activity_C=activity_c,
            activity_S=activity_s,
        )

    @staticmethod
    def _config() -> ReliabilityConfig:
        return ReliabilityConfig(
            minimum_valid_spots=30,
            minimum_valid_fraction=0.80,
            warning_valid_fraction=0.50,
            bootstrap_iterations=7,
            permutation_iterations=7,
            seed=42,
        )

    def _summary(self, defined_n: int, total: int = 100) -> pd.Series:
        result = self._result(defined_n, total)
        return build_pair_summary("pair", result, result, self._config()).iloc[0]

    def _assert_inference_missing(self, row: pd.Series) -> None:
        for metric in ("D", "CA_fraction"):
            self.assertTrue(np.isnan(row[f"delta_{metric}_bootstrap_ci_low"]))
            self.assertTrue(np.isnan(row[f"delta_{metric}_bootstrap_ci_high"]))
            self.assertTrue(np.isnan(row[f"delta_{metric}_permutation_p_value"]))
            self.assertTrue(np.isnan(row[f"delta_{metric}_bh_fdr"]))

    def test_a_two_defined_spots_keep_pair_valid_but_block_direction_and_ca_inference(self) -> None:
        row = self._summary(2)
        self.assertEqual(row["pair_score_validity"], "valid")
        self.assertEqual(row["direction_qc_status"], "FAIL")
        self.assertEqual(row["ca_qc_status"], "FAIL")
        self.assertFalse(row["direction_inference_eligible"])
        self.assertFalse(row["ca_inference_eligible"])
        self.assertTrue(np.isfinite(row["pre_D"]))
        self.assertTrue(np.isfinite(row["pre_CA_fraction"]))
        self._assert_inference_missing(row)

    def test_b_twenty_nine_defined_spots_fail_even_above_eighty_percent(self) -> None:
        row = self._summary(29, total=30)
        self.assertAlmostEqual(row["pre_direction_defined_fraction"], 29 / 30)
        self.assertEqual(row["direction_qc_status"], "FAIL")
        self.assertEqual(row["ca_qc_status"], "FAIL")
        self.assertIn("defined_spot_count_below_minimum", row["direction_qc_reason"])
        self._assert_inference_missing(row)

    def test_c_count_and_fraction_pass_boundaries_are_inclusive(self) -> None:
        # Integer spot counts cannot simultaneously be exactly n=30 and 80%
        # (30 / 0.80 = 37.5), so the two inclusive boundaries are tested
        # separately: exactly 30 spots, then exactly 80% support.
        exactly_thirty = self._summary(30, total=36)
        exactly_eighty_percent = self._summary(32, total=40)
        for row in (exactly_thirty, exactly_eighty_percent):
            self.assertEqual(row["direction_qc_status"], "PASS")
            self.assertEqual(row["ca_qc_status"], "PASS")
            self.assertTrue(row["direction_inference_eligible"])
            self.assertTrue(row["ca_inference_eligible"])
            self.assertTrue(np.isfinite(row["delta_D_permutation_p_value"]))
            self.assertTrue(np.isfinite(row["delta_D_bh_fdr"]))
            self.assertTrue(np.isfinite(row["delta_CA_fraction_permutation_p_value"]))
            self.assertTrue(np.isfinite(row["delta_CA_fraction_bh_fdr"]))
        self.assertEqual(exactly_thirty["pre_direction_defined_n"], 30)
        self.assertEqual(exactly_eighty_percent["pre_direction_defined_fraction"], 0.80)

    def test_d_seventy_nine_percent_is_caution_and_has_no_inference(self) -> None:
        row = self._summary(79)
        self.assertEqual(row["direction_qc_status"], "CAUTION")
        self.assertEqual(row["ca_qc_status"], "CAUTION")
        self.assertIn("defined_fraction_below_pass_threshold", row["ca_qc_reason"])
        self._assert_inference_missing(row)

    def test_e_exactly_fifty_percent_is_caution(self) -> None:
        row = self._summary(50)
        self.assertEqual(row["direction_qc_status"], "CAUTION")
        self.assertEqual(row["ca_qc_status"], "CAUTION")
        self.assertFalse(row["direction_inference_eligible"])
        self.assertFalse(row["ca_inference_eligible"])
        self._assert_inference_missing(row)

    def test_f_below_fifty_percent_fails(self) -> None:
        row = self._summary(49)
        self.assertEqual(row["direction_qc_status"], "FAIL")
        self.assertEqual(row["ca_qc_status"], "FAIL")
        self.assertIn("defined_fraction_below_minimum", row["direction_qc_reason"])
        self._assert_inference_missing(row)

    def test_g_zero_defined_spots_are_explicitly_missing_and_fail(self) -> None:
        row = self._summary(0)
        self.assertTrue(np.isnan(row["pre_D"]))
        self.assertTrue(np.isnan(row["pre_CA_fraction"]))
        self.assertEqual(row["direction_qc_status"], "FAIL")
        self.assertEqual(row["ca_qc_status"], "FAIL")
        self.assertIn("no_defined_spots", row["direction_qc_reason"])
        self._assert_inference_missing(row)

    def test_h_bh_excludes_missing_ineligible_metrics(self) -> None:
        adjusted = bh_adjust([0.01, 0.04, np.nan, np.nan])
        np.testing.assert_allclose(adjusted[:2], [0.02, 0.04])
        self.assertTrue(np.isnan(adjusted[2:]).all())

        row = self._summary(2)
        self.assertTrue(np.isfinite(row["delta_A_permutation_p_value"]))
        self.assertTrue(np.isfinite(row["delta_A_bh_fdr"]))
        self.assertTrue(np.isnan(row["delta_D_permutation_p_value"]))
        self.assertTrue(np.isnan(row["delta_D_bh_fdr"]))
        self.assertTrue(np.isnan(row["delta_CA_fraction_permutation_p_value"]))
        self.assertTrue(np.isnan(row["delta_CA_fraction_bh_fdr"]))

    def test_metric_qc_table_exports_explicit_denominators_masks_and_reasons(self) -> None:
        result = self._result(50)
        table = build_metric_qc_table("pair", result, result, self._config())
        self.assertEqual(len(table), 2)
        for column in (
            "direction_defined_n",
            "direction_valid_input_n",
            "direction_defined_fraction",
            "direction_qc_status",
            "direction_inference_eligible",
            "direction_qc_reason",
            "ca_defined_n",
            "ca_valid_input_n",
            "ca_defined_fraction",
            "ca_qc_status",
            "ca_inference_eligible",
            "ca_qc_reason",
        ):
            self.assertIn(column, table.columns)
        self.assertTrue(table["direction_ca_share_defined_mask"].all())
        self.assertTrue(table["defined_fraction_denominator"].eq("valid_input_n").all())


class ReliabilityDependenceTests(unittest.TestCase):
    def test_high_dependence_is_warned_but_axes_are_not_changed(self) -> None:
        rng = np.random.default_rng(42)
        C = rng.uniform(1.0, 4.0, size=80)
        S = rng.uniform(0.1, 0.8, size=80)
        axis_1 = compute_axis_reliability(C, S, axis="axis_1")
        axis_2 = compute_axis_reliability(C * 2.0, S * 2.0, axis="axis_2")
        original = axis_2.balance_B.copy()
        config = ReliabilityConfig(permutation_iterations=25, minimum_valid_spots=10)
        first = compute_axis_dependence(
            {"axis_1": axis_1, "axis_2": axis_2}, config, sample_id="pre"
        )
        second = compute_axis_dependence(
            {"axis_1": axis_1, "axis_2": axis_2}, config, sample_id="pre"
        )
        self.assertEqual(len(first), 5)
        self.assertTrue(first["qc_status"].eq("warning_high_dependence").any())
        np.testing.assert_allclose(first["permutation_p_value"], second["permutation_p_value"])
        np.testing.assert_allclose(axis_2.balance_B, original)
        self.assertTrue(first["bh_fdr"].dropna().between(0, 1).all())

    def test_single_axis_dependence_is_explicitly_empty(self) -> None:
        axis = compute_axis_reliability([1.0, 2.0, 3.0], [0.5, 0.5, 0.5])
        table = compute_axis_dependence({"axis": axis})
        self.assertTrue(table.empty)
        self.assertIn("bh_fdr", table.columns)

    def test_direction_dependence_inference_uses_metric_support_gate(self) -> None:
        total = 100
        sparse = compute_axis_reliability(
            np.linspace(0.0, 1.0, total),
            np.linspace(1.0, 0.0, total),
            activity_C=np.r_[np.ones(2), np.zeros(total - 2)],
            activity_S=np.zeros(total),
            axis="sparse",
        )
        dense = compute_axis_reliability(
            np.linspace(0.0, 1.0, total),
            np.linspace(1.0, 0.0, total),
            activity_C=np.linspace(0.1, 1.0, total),
            activity_S=np.linspace(1.0, 0.1, total),
            axis="dense",
        )
        table = compute_axis_dependence(
            {"sparse": sparse, "dense": dense},
            ReliabilityConfig(permutation_iterations=7),
        )
        direction_row = table.loc[table["dependence_type"].eq("direction_dependence")].iloc[0]
        self.assertFalse(direction_row["metric_inference_eligible"])
        self.assertEqual(direction_row["qc_status"], "insufficient_direction_defined_support")
        self.assertTrue(np.isnan(direction_row["permutation_p_value"]))
        self.assertTrue(np.isnan(direction_row["bh_fdr"]))


class ReliabilityIntegrationTests(unittest.TestCase):
    @staticmethod
    def _comparative_config(**overrides) -> ComparativeConfig:
        values = {
            "mode": "pairwise",
            "reference": "Pre",
            "target": "Post",
            "c_genes": C_GENES,
            "s_genes": S_GENES,
            "enable_h_expr": False,
            "enable_v_expr": False,
            "bootstrap_iterations": 100,
        }
        values.update(overrides)
        return ComparativeConfig(**values)

    def test_feature_flag_off_preserves_v06_files_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre.h5ad", pattern="localized", seed=31)
            post = write_comparative_h5ad(root / "post.h5ad", pattern="diffuse", seed=32)
            result = run_multi_pair_analysis(
                [PairSpec("off", pre, post)],
                self._comparative_config(),
                root / "results",
                run_tag="off",
            )
            self.assertFalse((result.run_dir / "reliability_spot_results.csv").exists())
            metadata = json.loads((result.run_dir / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertNotIn("reliability_layer", metadata)
            self.assertTrue(result.reliability_spot_results.empty)

    def test_enabled_layer_exports_only_additive_sidecars(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre.h5ad", pattern="localized", seed=41)
            post = write_comparative_h5ad(root / "post.h5ad", pattern="diffuse", seed=42)
            pair = PairSpec("enabled", pre, post)
            baseline = run_multi_pair_analysis(
                [pair], self._comparative_config(), root / "results", run_tag="baseline"
            )
            reliability = ReliabilityConfig(
                enabled=True,
                bootstrap_iterations=10,
                permutation_iterations=10,
                minimum_valid_spots=3,
            )
            result = run_multi_pair_analysis(
                [pair],
                self._comparative_config(),
                root / "results",
                run_tag="reliability",
                reliability_config=reliability,
            )
            for column in baseline.pair_results.columns:
                left, right = baseline.pair_results.iloc[0][column], result.pair_results.iloc[0][column]
                if isinstance(left, str):
                    self.assertEqual(left, right, msg=column)
                elif pd.isna(left) and pd.isna(right):
                    continue
                else:
                    self.assertAlmostEqual(float(left), float(right), places=12, msg=column)
            for name in (
                "reliability_spot_results.csv",
                "reliability_pair_summary.csv",
                "reliability_metric_qc.csv",
                "reliability_gene_coverage.csv",
                "reliability_qc.json",
                "cross_exclusivity_audit.csv",
                "axis_dependence_long.csv",
                "axis_dependence_matrix.csv",
                "axis_dependence_heatmap.png",
                "reliability_valid_fraction_qc.png",
                "reliability_score_domain_diagnostic.csv",
                "reliability_score_domain_diagnostic.json",
            ):
                self.assertTrue((result.run_dir / name).is_file(), msg=name)
            self.assertIn("delta_B", result.reliability_pair_summary.columns)
            self.assertIn("delta_CA_strength", result.reliability_pair_summary.columns)
            self.assertTrue(result.reliability_spot_results["reliability_status"].str.len().gt(0).all())
            self.assertTrue(result.reliability_spot_results["activity_C_input"].ge(0).all())
            self.assertTrue(result.reliability_spot_results["activity_S_input"].ge(0).all())
            self.assertTrue(result.reliability_spot_results["C_input"].lt(0).any())
            summary_row = result.reliability_pair_summary.iloc[0]
            self.assertEqual(summary_row["balance_score_source"], "legacy_signed_cs")
            self.assertEqual(summary_row["balance_score_domain"], "signed")
            self.assertEqual(summary_row["activity_score_domain"], "nonnegative")
            self.assertEqual(summary_row["inference_level"], "spot_distribution_descriptive")
            self.assertFalse(summary_row["treatment_effect_claim_allowed"])
            metadata = json.loads((result.run_dir / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(
                metadata["reliability_layer"]["reliability_schema_version"],
                "v0.65-reliability-v3-metric-qc",
            )
            self.assertFalse(metadata["reliability_layer"]["H_V_reinterpreted_as_paired_pole_axes"])
            self.assertEqual(metadata["reliability_layer"]["balance_score_source"], "legacy_signed_cs")
            self.assertFalse(metadata["treatment_effect_claim_allowed"])
            self.assertFalse(metadata["spotwise_subtraction_performed"])

    def test_strict_overlap_writes_audit_and_block_reason_before_stopping(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = self._comparative_config(
                enable_h_expr=True,
                h_genes=[C_GENES[0]],
            )
            with self.assertRaises(CrossExclusivityError):
                run_multi_pair_analysis(
                    [PairSpec("blocked", root / "pre.h5ad", root / "post.h5ad")],
                    config,
                    root / "results",
                    run_tag="blocked",
                    reliability_config=ReliabilityConfig(enabled=True),
                )
            metadata_files = list((root / "results").rglob("run_metadata.json"))
            self.assertEqual(len(metadata_files), 1)
            run_dir = metadata_files[0].parent
            self.assertTrue((run_dir / "cross_exclusivity_audit.csv").is_file())
            metadata = json.loads((run_dir / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["analysis_status"], "blocked_cross_exclusivity")
            self.assertIn(C_GENES[0], metadata["block_reason"])


if __name__ == "__main__":
    unittest.main()
