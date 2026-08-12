from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import pandas as pd

from spatialtx_desktop.comparative.models import ComparativeConfig
from spatialtx_desktop.comparative.multi_pair import (
    ComparabilityConfig,
    PairInterpretationConfig,
    PairSpec,
    _apply_within_pair_context_thresholds,
    _context_gene_audit_row,
    _pair_interpretation_table,
    direction_symbol,
    evaluate_comparability,
    infer_pair_id_from_filename,
    run_multi_pair_analysis,
    validate_pair_identity,
    validate_pair_specs,
)
from spatialtx_desktop.comparative.multiaxial import (
    SITE_SHIFT_WARNING,
    _context_plot_layers,
    qc_aware_interpretation,
)
from spatialtx_desktop.comparative.metrics import COMPARATIVE_METRIC_LAYER_SCHEMA
from spatialtx_desktop.graph.context import (
    DEFAULT_HYPOXIA_GENES,
    DEFAULT_VASCULAR_PROXY_GENES,
)
from tests.comparative_fixtures import C_GENES, S_GENES, write_comparative_h5ad


def _config() -> ComparativeConfig:
    return ComparativeConfig(
        mode="pairwise",
        reference="Pre",
        target="Post",
        c_genes=C_GENES,
        s_genes=S_GENES,
        enable_h_expr=False,
        enable_v_expr=False,
        bootstrap_iterations=100,
    )


def _context_config() -> ComparativeConfig:
    return ComparativeConfig(
        mode="pairwise",
        reference="Pre",
        target="Post",
        c_genes=C_GENES,
        s_genes=S_GENES,
        enable_h_expr=True,
        enable_v_expr=True,
        h_genes=["CA9", "VEGFA"],
        v_genes=["PECAM1", "VWF"],
        context_smoothing="none",
        bootstrap_iterations=100,
    )


class MultiPairComparativeTests(unittest.TestCase):
    def test_one_pair_exports_pre_post_delta_qc_metadata_and_figures(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre.h5ad", pattern="localized", include_context=False)
            post = write_comparative_h5ad(root / "post.h5ad", pattern="diffuse", include_context=False, seed=2)
            result = run_multi_pair_analysis(
                [PairSpec("Patient_01", pre, post)],
                _config(),
                root / "results",
                run_tag="one_pair",
            )

            self.assertEqual(len(result.pair_results), 1)
            self.assertEqual(result.pair_results.iloc[0]["status"], "PASS")
            self.assertEqual(result.pair_results.iloc[0]["percent_change_status_C"], "not_meaningful")
            self.assertIn("pre_C_mean", result.pair_results.columns)
            for column in (
                "pre_C", "post_C", "delta_C", "pre_S", "post_S", "delta_S",
                "pre_R", "post_R", "delta_R", "pre_interface_fraction",
                "post_interface_fraction", "delta_interface_fraction", "regime_transition",
                "comparability", "comparability_reasons",
            ):
                self.assertIn(column, result.pair_results.columns)
            for name in (
                "pair_results.csv", "balance_changes.csv", "spatial_organization_changes.csv",
                "specimen_reliability.csv", "comparative_overview.csv", "comparability_qc.csv",
                "pair_interpretation_summary.csv", "comparability_details.csv",
                "overview_interpretation.csv", "cohort_summary.csv", "run_metadata.json",
                "context_changes.csv", "multiaxial_pair_summary.csv",
                "comparative_qc_summary.csv", "context_gene_audit.csv",
            ):
                self.assertTrue((result.run_dir / name).is_file(), msg=name)
            self.assertEqual(result.balance_changes.iloc[0]["result_layer"], "Balance change")
            self.assertIn("delta_R", result.balance_changes.columns)
            self.assertNotIn("delta_interface_fraction", result.balance_changes.columns)
            self.assertEqual(
                result.spatial_organization_changes.iloc[0]["result_layer"],
                "Spatial organization change",
            )
            self.assertIn("delta_interface_fraction", result.spatial_organization_changes.columns)
            self.assertNotIn("delta_R", result.spatial_organization_changes.columns)
            self.assertEqual(result.specimen_reliability.iloc[0]["result_layer"], "Specimen reliability")
            self.assertIn("balance_change_summary", result.comparative_overview.columns)
            self.assertIn("spatial_organization_summary", result.comparative_overview.columns)
            self.assertIn("specimen_reliability", result.comparative_overview.columns)
            self.assertEqual(len(result.pair_interpretation_summary), 1)
            self.assertIn(
                result.pair_interpretation_summary.iloc[0]["balance_change_class"],
                ("Minimal", "Moderate", "Large"),
            )
            self.assertIn("spot_count_ratio", result.comparability_details.columns)
            self.assertIn("technical_mismatch_reason", result.comparability_details.columns)
            self.assertIn("interpretive_flag", result.overview_interpretation.columns)
            self.assertTrue((result.run_dir / "figures" / "multi_pair_comparative_overview.png").is_file())
            self.assertTrue((result.run_dir / "figures" / "multiaxial_pair_overview.png").is_file())
            self.assertTrue((root / "results" / ".spatialtx_comparative_cache").is_dir())
            self.assertFalse((result.run_dir / ".spatialtx_comparative_cache").exists())
            metadata = json.loads((result.run_dir / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["delta_definition"], "Post - Pre")
            self.assertEqual(metadata["pair_count_completed"], 1)
            self.assertIn("comparability_configuration", metadata)
            self.assertEqual(metadata["result_layers"]["layer_1"], "Balance change: C, S, and R = C - S")
            self.assertIn("not combined", metadata["result_layers"]["combination_rule"])
            self.assertIn("pair_interpretation_configuration", metadata)
            self.assertFalse(metadata["clinical_interpretation_performed"])
            self.assertEqual(len(result.context_gene_audit), 4)
            self.assertTrue(result.context_gene_audit["context_status"].eq("not_requested").all())

    def test_existing_context_axes_are_parallel_and_do_not_change_core_results(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre.h5ad", pattern="localized", seed=7)
            post = write_comparative_h5ad(root / "post.h5ad", pattern="diffuse", seed=8)
            pair = PairSpec("ContextPair", pre, post, "same_site")

            core = run_multi_pair_analysis(
                [pair], _config(), root / "results", run_tag="core_only"
            ).pair_results.iloc[0]
            context_result = run_multi_pair_analysis(
                [pair], _context_config(), root / "results", run_tag="with_context"
            )
            context = context_result.pair_results.iloc[0]

            for column in (
                "regime_pre", "regime_post", "regime_transition",
                "pre_interface_fraction", "post_interface_fraction",
                "pre_diffuse_fraction", "post_diffuse_fraction",
                "pre_transition_burden", "post_transition_burden",
                "comparability",
            ):
                if isinstance(core[column], str):
                    self.assertEqual(context[column], core[column], msg=column)
                else:
                    self.assertAlmostEqual(float(context[column]), float(core[column]), places=12, msg=column)
            for column in ("pre_H", "post_H", "delta_H", "pre_V", "post_V", "delta_V"):
                self.assertTrue(pd.notna(context[column]), msg=column)
            for axis in ("H", "V"):
                self.assertEqual(context[f"pre_{axis}_context_status"], "available")
                self.assertEqual(context[f"post_{axis}_context_status"], "available")
                self.assertAlmostEqual(
                    float(context[f"{axis}_pair_pooled_q90"]),
                    float(context_result.context_changes.iloc[0][f"{axis}_pair_pooled_q90"]),
                )
            summary = context_result.multiaxial_pair_summary.iloc[0]
            for column in (
                "cs_balance_pre", "cs_balance_post", "delta_cs_balance",
                "interface_pre", "interface_post", "delta_interface",
                "diffuse_pre", "diffuse_post", "delta_diffuse",
                "burden_pre", "burden_post", "delta_burden",
                "H_pre", "H_post", "delta_H", "V_pre", "V_post", "delta_V",
            ):
                self.assertIn(column, context_result.multiaxial_pair_summary.columns)
                self.assertTrue(pd.notna(summary[column]), msg=column)
            for column in (
                "H_q90_pre", "H_q90_post", "delta_H_q90",
                "H_high_fraction_pre", "H_high_fraction_post", "delta_H_high_fraction",
                "H_local_fraction_pre", "H_local_fraction_post", "delta_H_local_fraction",
                "V_q90_pre", "V_q90_post", "delta_V_q90",
                "V_high_fraction_pre", "V_high_fraction_post", "delta_V_high_fraction",
                "V_local_fraction_pre", "V_local_fraction_post", "delta_V_local_fraction",
                "H_gene_coverage_pre", "H_gene_coverage_post",
                "V_gene_coverage_pre", "V_gene_coverage_post",
                "H_context_status_pre", "H_context_status_post",
                "V_context_status_pre", "V_context_status_post",
            ):
                self.assertIn(column, context_result.multiaxial_pair_summary.columns)
            audit = context_result.context_gene_audit
            self.assertEqual(len(audit), 4)
            self.assertEqual(set(audit["axis"]), {"H", "V"})
            self.assertEqual(set(audit["sample_role"]), {"Pre", "Post"})
            self.assertTrue(audit["context_status"].eq("available").all())
            metadata = json.loads((context_result.run_dir / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["comparative_metric_layer_schema"], COMPARATIVE_METRIC_LAYER_SCHEMA)
            self.assertEqual(metadata["pair_pooled_high_quantile"], 0.90)
            self.assertEqual(metadata["pair_pooled_threshold_scope"], "within_pair_pre_plus_post")
            self.assertEqual(
                metadata["context_warning_provenance"]["single_sample_context_warning"],
                "legacy_within_sample_centered_context_q80",
            )
            self.assertEqual(
                metadata["context_warning_provenance"]["pair_pooled_context_warning"],
                "within_pair_pre_plus_post_raw_context_q90",
            )

    def test_default_programs_are_audited_and_low_coverage_is_nan_not_zero(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre.h5ad", include_context=True)
            post = write_comparative_h5ad(root / "post.h5ad", include_context=True, seed=12)
            result = run_multi_pair_analysis(
                [PairSpec("DefaultAudit", pre, post)],
                ComparativeConfig(
                    mode="pairwise",
                    reference="Pre",
                    target="Post",
                    c_genes=C_GENES,
                    s_genes=S_GENES,
                    enable_h_expr=True,
                    enable_v_expr=True,
                    h_genes=None,
                    v_genes=None,
                    bootstrap_iterations=100,
                ),
                root / "results",
                run_tag="default_audit",
            )

            row = result.pair_results.iloc[0]
            for axis in ("H", "V"):
                self.assertEqual(row[f"pre_{axis}_context_status"], "insufficient_gene_coverage")
                self.assertEqual(row[f"post_{axis}_context_status"], "insufficient_gene_coverage")
                self.assertTrue(pd.isna(row[f"pre_{axis}"]))
                self.assertTrue(pd.isna(row[f"post_{axis}"]))
            audit = result.context_gene_audit
            self.assertEqual(len(audit), 4)
            h = audit.loc[audit["axis"].eq("H")].iloc[0]
            v = audit.loc[audit["axis"].eq("V")].iloc[0]
            self.assertEqual(h["requested_genes"].split(";"), DEFAULT_HYPOXIA_GENES)
            self.assertEqual(v["requested_genes"].split(";"), DEFAULT_VASCULAR_PROXY_GENES)
            self.assertEqual(int(h["requested_gene_count"]), 12)
            self.assertEqual(int(v["requested_gene_count"]), 12)
            self.assertAlmostEqual(float(h["coverage_fraction"]), 2 / 12)
            self.assertAlmostEqual(float(v["coverage_fraction"]), 2 / 12)
            metadata = json.loads((result.run_dir / "run_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["effective_context_programs"]["H"]["source"], "default")
            self.assertEqual(metadata["effective_context_programs"]["H"]["requested_genes"], DEFAULT_HYPOXIA_GENES)
            self.assertEqual(metadata["effective_context_programs"]["V"]["requested_genes"], DEFAULT_VASCULAR_PROXY_GENES)

    def test_pair_pooled_threshold_is_shared_within_pair_and_isolated_between_pairs(self) -> None:
        def run_threshold(prefix: str, scale: float) -> tuple[dict, dict]:
            pre_metrics = {"sample_id": f"{prefix}_pre", "group": "Pre", "V_raw_median": 0.0}
            post_metrics = {"sample_id": f"{prefix}_post", "group": "Post", "V_raw_median": 0.0}
            pre_fields = {
                "H_expr_raw": pd.Series([0.0, 0.0, 1.0, 2.0]).to_numpy() * scale,
                "V_expr_raw": pd.Series([0.0, 0.0, 0.0, 1.0]).to_numpy() * scale,
                "context_edge_i": pd.Series([0, 1, 2]).to_numpy(),
                "context_edge_j": pd.Series([1, 2, 3]).to_numpy(),
            }
            post_fields = {
                "H_expr_raw": pd.Series([0.0, 1.0, 2.0, 3.0]).to_numpy() * scale,
                "V_expr_raw": pd.Series([0.0, 0.0, 0.0, 2.0]).to_numpy() * scale,
                "context_edge_i": pd.Series([0, 1, 2]).to_numpy(),
                "context_edge_j": pd.Series([1, 2, 3]).to_numpy(),
            }
            return _apply_within_pair_context_thresholds(
                pre_metrics, post_metrics, pre_fields, post_fields
            )

        first_pre, first_post = run_threshold("first", 1.0)
        second_pre, second_post = run_threshold("second", 10.0)
        for axis in ("H", "V"):
            self.assertEqual(
                first_pre[f"{axis}_pooled_high_threshold"],
                first_post[f"{axis}_pooled_high_threshold"],
            )
            self.assertEqual(
                second_pre[f"{axis}_pooled_high_threshold"],
                second_post[f"{axis}_pooled_high_threshold"],
            )
            self.assertNotEqual(
                first_pre[f"{axis}_pooled_high_threshold"],
                second_pre[f"{axis}_pooled_high_threshold"],
            )
        self.assertEqual(first_pre["V_raw_median"], 0.0)
        self.assertEqual(first_post["V_raw_median"], 0.0)
        self.assertGreater(float(first_post["V_high_fraction"]), 0.0)

    def test_context_audit_separates_legacy_and_pair_pooled_warning_provenance(self) -> None:
        metrics = {
            "V_context_audit": {
                "gene_set_name": "test_v_program",
                "source": "user_supplied",
                "requested_gene_count": 2,
                "matched_gene_count": 2,
                "coverage_fraction": 1.0,
                "requested_genes": ["PECAM1", "VWF"],
                "matched_genes": ["PECAM1", "VWF"],
                "missing_genes": [],
                "genes_expressed_above_min_spot_fraction": ["PECAM1", "VWF"],
                "expressed_gene_count": 2,
                "expressed_gene_fraction": 1.0,
                "expression_scale_guess": "log1p_normalized",
                "detection_source": "X",
                "raw_normalization_method": "existing_nonnegative_expression_then_program_mean",
            },
            "V_context_status": "available",
            "V_context_warning": "V_expr high-state fraction is unexpectedly small or large: 100.0%",
            "V_expr_high_state_fraction": 1.0,
            "V_high_fraction": 0.0504,
            "V_pooled_high_threshold": 0.1733,
        }
        row = _context_gene_audit_row(
            PairSpec("Pair 5", Path("pre.h5ad"), Path("post.h5ad")),
            "Pre",
            Path("pre.h5ad"),
            "V",
            metrics,
            _context_config(),
        )

        self.assertEqual(row["single_sample_context_high_fraction"], 1.0)
        self.assertEqual(row["pair_pooled_context_high_fraction"], 0.0504)
        self.assertIn("100.0%", row["single_sample_context_warning"])
        self.assertEqual(row["pair_pooled_context_warning"], "")
        self.assertEqual(
            row["single_sample_context_warning_provenance"],
            "legacy_within_sample_centered_context_q80",
        )
        self.assertEqual(
            row["pair_pooled_context_warning_provenance"],
            "within_pair_pre_plus_post_raw_context_q90",
        )
        self.assertIn("legacy within-sample centered-context q80", row["context_warning"])
        self.assertNotIn("Pair-pooled context QC", row["context_warning"])

    def test_missing_hv_is_non_blocking_and_exported_as_nan(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre.h5ad", include_context=False)
            post = write_comparative_h5ad(root / "post.h5ad", include_context=False, seed=5)
            result = run_multi_pair_analysis(
                [PairSpec("NoContext", pre, post)],
                _context_config(),
                root / "results",
                run_tag="missing_context",
            )

            pair = result.pair_results.iloc[0]
            self.assertEqual(pair["status"], "PASS")
            self.assertTrue(pd.isna(pair["pre_H"]))
            self.assertTrue(pd.isna(pair["post_H"]))
            self.assertTrue(pd.isna(pair["pre_V"]))
            self.assertTrue(pd.isna(pair["post_V"]))
            self.assertEqual(pair["pre_H_context_status"], "no_matched_genes")
            self.assertEqual(pair["post_H_context_status"], "no_matched_genes")
            self.assertEqual(pair["pre_V_context_status"], "no_matched_genes")
            self.assertEqual(pair["post_V_context_status"], "no_matched_genes")
            note = result.multiaxial_pair_summary.iloc[0]["interpretation_note"]
            self.assertIn("H: not available", note)
            self.assertIn("V: not available", note)

    def test_unsupported_expression_scale_is_explicit_and_non_blocking(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre_source.h5ad", include_context=True)
            post = write_comparative_h5ad(root / "post_source.h5ad", include_context=True, seed=13)
            for source, destination in ((pre, root / "pre.h5ad"), (post, root / "post.h5ad")):
                adata = ad.read_h5ad(source)
                adata.X = adata.X - 20.0
                adata.write_h5ad(destination)
            result = run_multi_pair_analysis(
                [PairSpec("UnsupportedScale", root / "pre.h5ad", root / "post.h5ad")],
                _context_config(),
                root / "results",
                run_tag="unsupported_scale",
            )
            pair = result.pair_results.iloc[0]
            self.assertEqual(pair["status"], "PASS")
            for axis in ("H", "V"):
                self.assertEqual(pair[f"pre_{axis}_context_status"], "unsupported_expression_scale")
                self.assertEqual(pair[f"post_{axis}_context_status"], "unsupported_expression_scale")
                self.assertTrue(pd.isna(pair[f"pre_{axis}"]))
                self.assertTrue(pd.isna(pair[f"post_{axis}"]))

    def test_different_site_is_visible_warning_not_exclusion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre.h5ad", include_context=False)
            post = write_comparative_h5ad(root / "post.h5ad", include_context=False, seed=9)
            result = run_multi_pair_analysis(
                [PairSpec("SiteShift", pre, post, "different_site")],
                _config(),
                root / "results",
                run_tag="site_shift",
            )

            pair = result.pair_results.iloc[0]
            self.assertEqual(pair["status"], "PASS")
            self.assertEqual(pair["site_comparability"], "different_site")
            qc_row = result.comparability_qc.loc[
                result.comparability_qc["category"].eq("site_metadata")
            ].iloc[0]
            self.assertFalse(bool(qc_row["primary_for_classification"]))
            summary = result.multiaxial_pair_summary.iloc[0]
            self.assertIn(summary["interpretation_confidence"], ("CAUTION", "LOW"))
            self.assertIn("SITE-SHIFT WARNING", summary["interpretation_note"])
            self.assertEqual(
                result.comparative_qc_summary.iloc[0]["site_shift_warning"],
                SITE_SHIFT_WARNING,
            )

    def test_qc_aware_wording_covers_low_and_spatial_redistribution(self) -> None:
        confidence, note = qc_aware_interpretation({
            "status": "PASS",
            "comparability": "Low",
            "site_comparability": "same_site",
            "regime_pre": "Type_A_candidate",
            "regime_post": "Type_A_candidate",
            "delta_interface_fraction": -0.08,
            "delta_diffuse_fraction": 0.12,
            "delta_H": float("nan"),
            "delta_V": float("nan"),
        })
        self.assertEqual(confidence, "LOW")
        self.assertIn("redistribution", note)
        self.assertIn("without a regime-level change", note)
        self.assertIn("low technical/spatial comparability", note)
        self.assertNotIn("treatment response", note.casefold())
        _confidence, focal_note = qc_aware_interpretation({
            "status": "PASS",
            "comparability": "Good",
            "site_comparability": "same_site",
            "regime_pre": "Type_A_candidate",
            "regime_post": "Type_A_candidate",
            "delta_interface_fraction": 0.0,
            "delta_diffuse_fraction": 0.0,
            "delta_H": 0.0,
            "delta_H_q90": 0.08,
            "delta_H_high_fraction": 0.10,
            "pre_H_context_status": "available",
            "post_H_context_status": "available",
            "delta_V": 0.0,
            "delta_V_q90": 0.06,
            "delta_V_high_fraction": 0.12,
            "pre_V_context_status": "available",
            "post_V_context_status": "available",
        })
        self.assertIn("median remained broadly stable", focal_note)
        self.assertIn("upper-tail high-context fraction increased", focal_note)
        self.assertNotIn("angiogenesis increased", focal_note.casefold())

    def test_multiaxial_context_plot_uses_separate_median_and_pair_pooled_layers(self) -> None:
        layers = _context_plot_layers(pd.DataFrame({
            "delta_H": [0.0, 0.1493],
            "delta_V": [0.0, 0.05776],
            "delta_H_high_fraction": [0.12, 0.13772],
            "delta_V_high_fraction": [-0.03, 0.10194],
        }))

        self.assertEqual(set(layers), {"median", "pair_pooled_high_fraction"})
        self.assertEqual(layers["median"][0].tolist(), [0.0, 0.1493])
        self.assertEqual(layers["median"][1].tolist(), [0.0, 0.05776])
        self.assertEqual(layers["pair_pooled_high_fraction"][0].tolist(), [0.12, 0.13772])
        self.assertEqual(layers["pair_pooled_high_fraction"][1].tolist(), [-0.03, 0.10194])

    def test_six_pairs_continue_when_one_pair_is_corrupted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pairs: list[PairSpec] = []
            for index in range(1, 7):
                pre = write_comparative_h5ad(root / f"pre_{index}.h5ad", include_context=False, seed=index)
                post = write_comparative_h5ad(
                    root / f"post_{index}.h5ad",
                    pattern="diffuse" if index % 2 else "localized",
                    include_context=False,
                    seed=20 + index,
                )
                pairs.append(PairSpec(f"P{index:02d}", pre, post))
            pairs[2].post_path.write_text("not an h5ad", encoding="utf-8")

            result = run_multi_pair_analysis(pairs, _config(), root / "results", run_tag="six_pairs")

            self.assertEqual(len(result.pair_results), 6)
            self.assertEqual(int(result.pair_results["status"].eq("PASS").sum()), 5)
            self.assertEqual(len(result.pair_interpretation_summary), 6)
            failed = result.pair_results.loc[result.pair_results["pair_label"].eq("P03")].iloc[0]
            self.assertEqual(failed["status"], "ERROR")
            self.assertIn("corrupted", failed["error"])
            self.assertEqual(set(result.comparative_overview["pair_label"]), {f"P{i:02d}" for i in range(1, 7)})

    def test_pair_limit_accepts_six_and_rejects_seven(self) -> None:
        pairs = [PairSpec(f"P{index}", Path(f"pre_{index}.h5ad"), Path(f"post_{index}.h5ad")) for index in range(1, 7)]
        self.assertEqual(len(validate_pair_specs(pairs)), 6)
        with self.assertRaisesRegex(ValueError, "between 1 and 6"):
            validate_pair_specs([
                *pairs,
                PairSpec("P7", Path("pre_7.h5ad"), Path("post_7.h5ad")),
            ])

    def test_intentional_sampling_mismatch_is_low_with_auditable_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pre = write_comparative_h5ad(root / "pre.h5ad", include_context=False)
            full_post = write_comparative_h5ad(root / "post_full.h5ad", include_context=False, seed=4)
            post = root / "post_small.h5ad"
            adata = ad.read_h5ad(full_post)
            adata[:4].copy().write_h5ad(post)

            result = run_multi_pair_analysis(
                [PairSpec("Mismatch", pre, post)],
                _config(),
                root / "results",
                run_tag="mismatch",
            )

            row = result.pair_results.iloc[0]
            self.assertEqual(row["comparability"], "Low")
            self.assertIn("n_spots mismatch", row["comparability_reasons"])
            qc = result.comparability_qc.loc[result.comparability_qc["qc_metric"].eq("n_spots")].iloc[0]
            self.assertEqual(qc["severity"], "low")
            self.assertTrue(bool(qc["primary_for_classification"]))

    def test_composition_proxy_alone_cannot_produce_low(self) -> None:
        technical = {
            "n_spots": 100,
            "n_features": 1000,
            "median_detected_genes_per_spot": 500,
            "median_total_counts": 2000,
            "in_tissue_fraction": 0.9,
            "low_quality_fraction": 0.05,
            "spatial_coordinates_valid": True,
            "C_gene_coverage": 1.0,
            "S_gene_coverage": 1.0,
        }
        pre_fields = {"C": [2.0] * 10, "S": [0.0] * 10}
        post_fields = {"C": [0.0] * 10, "S": [2.0] * 10}
        classification, reasons, table = evaluate_comparability(
            "ProxyOnly",
            technical,
            technical,
            ComparabilityConfig(),
            pre_fields=pre_fields,
            post_fields=post_fields,
        )
        self.assertEqual(classification, "Caution")
        proxy = table.loc[table["category"].eq("composition_proxy_secondary")]
        self.assertTrue(proxy["severity"].eq("low").all())
        self.assertFalse(proxy["primary_for_classification"].any())
        self.assertTrue(any("composition_proxy" in reason for reason in reasons))

    def test_direction_symbols_use_metric_specific_tolerance(self) -> None:
        self.assertEqual(direction_symbol(0.009, 0.01), "→")
        self.assertEqual(direction_symbol(0.011, 0.01), "↑")
        self.assertEqual(direction_symbol(-0.011, 0.01), "↓")

    def test_balance_shift_with_preserved_structure_is_reported_separately_from_low_reliability(self) -> None:
        pair_results = pd.DataFrame([{
            "pair_label": "Patient_42",
            "status": "PASS",
            "error": "",
            "comparability": "Low",
            "delta_C": 0.30,
            "delta_S": 0.10,
            "delta_R": 0.40,
            "delta_interface_fraction": 0.01,
            "delta_diffuse_fraction": 0.01,
            "delta_transition_burden": 0.05,
            "regime_pre": "Type_A_candidate",
            "regime_post": "Type_A_candidate",
            "regime_transition": "Type_A_candidate → Type_A_candidate",
            "pair_id_validation": "matched",
            "pair_id_warning": "",
        }])
        details = pd.DataFrame([{
            "pair_label": "Patient_42",
            "technical_mismatch_reason": "n_spots mismatch (6.97-fold)",
            "sampling_mismatch_reason": "",
            "composition_proxy_reason": "",
        }])
        interpretation = _pair_interpretation_table(
            pair_results,
            details,
            PairInterpretationConfig(),
        ).iloc[0]
        self.assertEqual(interpretation["balance_change_class"], "Moderate")
        self.assertEqual(interpretation["spatial_change_class"], "Minimal")
        self.assertEqual(interpretation["regime_preserved"], "yes")
        self.assertEqual(interpretation["structure_preserved"], "probably")
        self.assertIn("balance shift with preserved structure", interpretation["interpretive_flag"])
        self.assertIn("low comparability", interpretation["interpretive_flag"])
        self.assertIn("CAUTION:", interpretation["caution_message"])
        self.assertIn("Broad spatial regime is preserved", interpretation["note_message"])

    def test_filename_pair_id_validation_is_conservative_and_non_blocking(self) -> None:
        self.assertEqual(
            infer_pair_id_from_filename("GSM9452692_sample_43_pre.h5ad"),
            "sample_43",
        )
        matched = validate_pair_identity(
            "GSM9452692_sample_43_pre.h5ad",
            "GSM9452693_sample_43_post.h5ad",
        )
        self.assertEqual(matched["pair_id_validation"], "matched")
        mismatch = validate_pair_identity("Patient_42_pre.h5ad", "Patient_43_post.h5ad")
        self.assertEqual(mismatch["pair_id_validation"], "warning")
        self.assertIn("Possible pair-ID mismatch", mismatch["pair_id_warning"])
        accession_only = validate_pair_identity("GSM111_pre.h5ad", "GSM222_post.h5ad")
        self.assertEqual(accession_only["pair_id_validation"], "not_available")

    def test_interpretation_thresholds_are_auditable(self) -> None:
        with self.assertRaisesRegex(ValueError, "moderate < large"):
            PairInterpretationConfig(
                balance_moderate_abs_delta=1.0,
                balance_large_abs_delta=0.5,
            ).validate()


if __name__ == "__main__":
    unittest.main()
