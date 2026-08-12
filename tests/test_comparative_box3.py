from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from spatialtx_desktop.comparative.comparative_normalization import (
    CENTERED_HV_WARNING,
    SCALE_WARNING_BANNER,
    add_normalized_topology_metrics,
    apply_pooled_hv_thresholds,
    build_hv_summary,
    build_metric_change_table,
    build_scale_warnings,
    noncentered_context_values,
    percent_change,
    raw_context_summary,
    sample_scale_metrics,
    symmetric_percent_change,
)
from spatialtx_desktop.comparative.metric_registry import (
    METRIC_REGISTRY,
    metric_definition,
    metrics_for_plot_group,
)
from spatialtx_desktop.comparative.metrics import compute_delta_metrics, group_mean_deltas
from spatialtx_desktop.comparative.models import ComparativeConfig
from spatialtx_desktop.comparative.plotting import (
    PRIMARY_SPATIAL_STATE_METRICS,
    TOPOLOGY_COMPONENT_COMPLEXITY_METRICS,
    plot_metric_change,
    plot_regime_transitions,
    plot_standardized_metric_change,
)
from spatialtx_desktop.comparative.reporting import rules_based_summary


def _base_samples() -> pd.DataFrame:
    rows = []
    for sample_id, valid, diffuse in (("reference", 1000, 100), ("target", 500, 50)):
        row = {definition.internal_name: np.nan for definition in METRIC_REGISTRY}
        row.update({
            "sample_id": sample_id,
            "group": "A" if sample_id == "reference" else "B",
            "C_mean": 1.0 if sample_id == "reference" else 1.5,
            "S_mean": 2.0,
            "R_mean": -1.0 if sample_id == "reference" else -0.5,
            "R_std": 0.5,
            "gradient_mean": 0.2,
            "gradient_q90": 0.8,
            "localized_interface_fraction": 0.1,
            "diffuse_fraction": 0.2,
            "transition_burden_score": 0.3,
            "R_crossing_fraction": 0.15,
            "adj_same_fraction": 0.6,
            "adj_zero_fraction": 0.2,
            "adj_opposite_fraction": 0.2,
            "n_diffuse_components": diffuse,
            "n_small_components": diffuse // 2,
            "n_interface_components": 10 if sample_id == "reference" else 5,
            "interface_fragmentation_index": 0.2,
            "small_component_fraction": 0.5,
            "largest_diffuse_component_ratio": 0.3,
            "diffuse_components_per_1000_valid_spots": 100.0,
            "diffuse_components_per_1000_in_tissue_spots": 111.111,
            "diffuse_components_per_tissue_component": 50.0,
            "small_components_per_1000_valid_spots": 50.0,
            "transition_components_per_1000_transition_spots": 500.0,
            "interface_segments_per_1000_valid_spots": 10.0,
            "normalized_fragmentation_score": 110.0,
            "n_total_spots": valid,
            "n_valid_spots": valid,
            "n_in_tissue_spots": int(valid * 0.9),
            "n_transition_spots": int(valid * 0.2),
            "tissue_area_proxy": valid,
            "tissue_component_count": 2 if sample_id == "reference" else 5,
            "spatial_extent_x": 40 if sample_id == "reference" else 20,
            "spatial_extent_y": 20 if sample_id == "reference" else 10,
            "spatial_extent_area_proxy": 800 if sample_id == "reference" else 200,
            "mean_spot_spacing": 1.0,
            "H_expr_mean": 1e-17,
            "V_expr_mean": -1e-18,
            "regime_label": "Type_A_candidate",
            "regime_confidence": 0.4,
            "spatial_qc_status": "PASS",
            "H_expr_available": True,
            "V_expr_available": True,
        })
        rows.append(row)
    return pd.DataFrame(rows)


def _matches() -> pd.DataFrame:
    return pd.DataFrame([{
        "comparison_id": "reference_vs_target",
        "pair_id": "p1",
        "reference_sample_id": "reference",
        "target_sample_id": "target",
    }])


class Box3NormalizationTests(unittest.TestCase):
    def test_raw_delta_is_preserved_while_normalized_density_can_be_unchanged(self) -> None:
        samples = _base_samples()
        deltas = compute_delta_metrics(samples, _matches())
        changes = build_metric_change_table(deltas).set_index("metric_name")
        self.assertEqual(changes.loc["n_diffuse_components", "raw_delta"], -50)
        self.assertEqual(changes.loc["diffuse_components_per_1000_valid_spots", "raw_delta"], 0)
        self.assertEqual(changes.loc["n_diffuse_components", "normalized_delta"], 0)
        warnings = build_scale_warnings(changes.reset_index())
        self.assertTrue(warnings["warning_code"].eq("raw_normalized_change_discordance").any())
        self.assertTrue(warnings["message"].str.contains("raw difference may partly reflect", case=False).any())

    def test_normalized_topology_formulas(self) -> None:
        values = add_normalized_topology_metrics({
            "n_diffuse_components": 4,
            "n_interface_components": 2,
            "small_component_fraction": 0.5,
            "n_valid_spots": 200,
            "n_in_tissue_spots": 160,
            "tissue_component_count": 2,
            "n_transition_spots": 40,
        })
        self.assertEqual(values["n_small_components"], 2)
        self.assertEqual(values["diffuse_components_per_1000_valid_spots"], 20)
        self.assertEqual(values["diffuse_components_per_1000_in_tissue_spots"], 25)
        self.assertEqual(values["diffuse_components_per_tissue_component"], 2)
        self.assertEqual(values["small_components_per_1000_valid_spots"], 10)
        self.assertEqual(values["transition_components_per_1000_transition_spots"], 100)
        self.assertEqual(values["normalized_fragmentation_score"], 30)

    def test_zero_and_near_zero_reference_percent_change_is_na(self) -> None:
        self.assertEqual(percent_change(0.0, 1.0)[1], "unstable_reference")
        self.assertTrue(np.isnan(percent_change(0.0, 1.0)[0]))
        self.assertEqual(percent_change(1e-12, 2.0)[1], "unstable_reference")
        self.assertTrue(np.isnan(percent_change(1e-12, 2.0)[0]))

    def test_symmetric_percent_change(self) -> None:
        self.assertAlmostEqual(symmetric_percent_change(10.0, 20.0), 200.0 / 3.0)
        self.assertAlmostEqual(symmetric_percent_change(20.0, 10.0), -200.0 / 3.0)
        self.assertEqual(symmetric_percent_change(0.0, 0.0), 0.0)

    def test_scale_warning_thresholds_and_missing_in_tissue(self) -> None:
        samples = _base_samples()
        samples["n_in_tissue_spots"] = samples["n_in_tissue_spots"].astype(float)
        samples.loc[:, "n_in_tissue_spots"] = np.nan
        changes = build_metric_change_table(compute_delta_metrics(samples, _matches()))
        warnings = build_scale_warnings(changes)
        self.assertTrue(warnings["message"].eq(SCALE_WARNING_BANNER).any())
        self.assertTrue(warnings["warning_code"].eq("in_tissue_count_unavailable").any())

    def test_sample_scale_uses_spot_proxy_and_separate_physical_area_label(self) -> None:
        coords = np.asarray([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
        adata = ad.AnnData(np.ones((4, 2)), obs=pd.DataFrame({"in_tissue": [1, 1, 0, 1]}))
        adata.obsm["spatial"] = coords
        fields = {name: np.ones(4) for name in ("C", "S", "R", "G")}
        fields["coords"] = coords
        graph = SimpleNamespace(connectivities=sparse.csr_matrix([
            [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]
        ]))
        metrics = sample_scale_metrics(adata, fields, np.asarray([1, 0, 0, 1], dtype=bool), graph)
        self.assertEqual(metrics["n_total_spots"], 4)
        self.assertEqual(metrics["n_valid_spots"], 4)
        self.assertEqual(metrics["n_in_tissue_spots"], 3)
        self.assertEqual(metrics["n_transition_spots"], 2)
        self.assertEqual(metrics["tissue_area_proxy"], 4)
        self.assertEqual(metrics["tissue_component_count"], 2)
        self.assertTrue(np.isnan(metrics["physical_tissue_area_estimate"]))
        self.assertIn("not_physical_area", metrics["spatial_extent_area_proxy_unit"])

    def test_registry_is_central_ordered_and_marks_deprecated_hv_means(self) -> None:
        names = [definition.internal_name for definition in METRIC_REGISTRY]
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual([d.internal_name for d in metrics_for_plot_group("program")], ["C_mean", "S_mean", "R_mean"])
        self.assertEqual(metric_definition("H_expr_mean").deprecated, True)
        self.assertEqual(metric_definition("H_expr_mean").observational_only, True)
        self.assertEqual(metric_definition("n_diffuse_components").scale_sensitive, True)

    def test_group_standardized_delta_uses_pooled_scale_but_pairwise_does_not(self) -> None:
        samples = _base_samples()
        grouped = group_mean_deltas(samples, "A", "B")
        pairwise = compute_delta_metrics(samples, _matches())
        self.assertTrue(grouped.loc[grouped["metric"].eq("C_mean"), "standardized_delta"].notna().all())
        self.assertTrue(pairwise["standardized_delta"].isna().all())
        self.assertTrue(pairwise["standardized_delta_status"].eq("not_computed_for_pairwise").all())

    def test_comparative_overview_uses_two_independent_raw_delta_panels(self) -> None:
        changes = build_metric_change_table(compute_delta_metrics(_base_samples(), _matches()))
        config = ComparativeConfig(mode="pairwise", reference="A", target="B", c_genes=["C1"], s_genes=["S1"])
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "comparative_metric_changes.png"
            plot_metric_change(changes, output, config, effective_mode="pairwise")
            metadata = json.loads(output.with_suffix(".png.metadata.json").read_text(encoding="utf-8"))
        self.assertFalse(metadata["x_axes_shared"])
        self.assertEqual(metadata["display_mode"], "raw_delta")
        self.assertEqual(
            metadata["panel_titles"],
            {
                "panel_a": "Primary spatial-state summary metrics",
                "panel_b": "Topology / component complexity metrics",
            },
        )
        self.assertEqual(metadata["panel_metric_names"]["panel_a"], list(PRIMARY_SPATIAL_STATE_METRICS))
        self.assertEqual(
            metadata["panel_metric_names"]["panel_b"],
            list(TOPOLOGY_COMPONENT_COMPLEXITY_METRICS),
        )
        self.assertTrue(metadata["numeric_value_labels"])
        self.assertIsNone(
            plot_standardized_metric_change(
                changes,
                Path(tempfile.gettempdir()) / "pairwise_standardized_should_not_exist.png",
                config,
                effective_mode="pairwise",
            )
        )

    def test_group_overview_can_export_pooled_scale_standardized_change(self) -> None:
        samples = _base_samples()
        target = samples["group"].eq("B")
        for index, metric in enumerate(
            (*PRIMARY_SPATIAL_STATE_METRICS, *TOPOLOGY_COMPONENT_COMPLEXITY_METRICS),
            start=1,
        ):
            samples.loc[target, metric] = float(samples.loc[target, metric].iloc[0]) + index / 10
        changes = build_metric_change_table(group_mean_deltas(samples, "A", "B"))
        self.assertIn("standardized_delta", changes)
        config = ComparativeConfig(mode="unpaired", reference="A", target="B", c_genes=["C1"], s_genes=["S1"])
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "comparative_metric_changes_standardized.png"
            result = plot_standardized_metric_change(changes, output, config, effective_mode="unpaired")
            self.assertEqual(result, output)
            metadata = json.loads(output.with_suffix(".png.metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["display_mode"], "standardized_delta")
        self.assertFalse(metadata["x_axes_shared"])


class Box3HVTests(unittest.TestCase):
    def test_noncentered_values_and_transition_enrichment(self) -> None:
        adata = ad.AnnData(
            np.asarray([[0, 2], [2, 4], [4, 6], [6, 8]], dtype=float),
            var=pd.DataFrame(index=["H1", "H2"]),
        )
        values, method = noncentered_context_values(adata, ["H1", "H2"], "raw_counts")
        self.assertIsNotNone(values)
        self.assertEqual(method, "log1p_counts_then_program_mean")
        summary = raw_context_summary("H", values, np.asarray([1, 1, 0, 0], dtype=bool), method)
        self.assertTrue(summary["H_raw_available"])
        expected = np.median(values[:2]) - np.median(values[2:])
        self.assertAlmostEqual(summary["H_transition_enrichment"], expected)

    def test_pooled_threshold_does_not_force_each_sample_fraction_constant(self) -> None:
        metrics = pd.DataFrame([
            {"sample_id": "r", "group": "A"}, {"sample_id": "t", "group": "B"}
        ])
        fields = {
            "r": {"H_expr_raw": np.asarray([0.0, 1.0]), "V_expr_raw": np.asarray([0.0, 0.0])},
            "t": {"H_expr_raw": np.asarray([2.0, 3.0]), "V_expr_raw": np.asarray([1.0, 4.0])},
        }
        result = apply_pooled_hv_thresholds(metrics, fields, "A", "B").set_index("sample_id")
        self.assertEqual(result.loc["r", "H_high_fraction"], 0.0)
        self.assertEqual(result.loc["t", "H_high_fraction"], 0.5)
        self.assertEqual(result.loc["r", "V_high_fraction"], 0.0)
        self.assertEqual(result.loc["t", "V_high_fraction"], 0.5)
        self.assertEqual(result.loc["r", "H_high_threshold_method"], "pooled_reference_target_q90")

    def test_missing_noncentered_scores_are_not_fabricated(self) -> None:
        values, method = noncentered_context_values(
            ad.AnnData(np.asarray([[-1.0], [1.0]]), var=pd.DataFrame(index=["H1"])),
            ["H1"],
            "centered_or_scaled",
        )
        self.assertIsNone(values)
        summary = raw_context_summary("H", None, np.asarray([1, 0], dtype=bool), method)
        self.assertFalse(summary["H_raw_available"])
        self.assertTrue(np.isnan(summary["H_raw_median"]))

    def test_centered_means_are_flagged_and_retained_for_compatibility(self) -> None:
        samples = _base_samples()
        changes = build_metric_change_table(compute_delta_metrics(samples, _matches()))
        centered = changes.loc[changes["metric_name"].isin(["H_expr_mean", "V_expr_mean"])]
        self.assertTrue(centered["deprecated"].all())
        self.assertTrue(centered["observational_only"].all())
        self.assertTrue(centered["interpretation_flag"].eq("non_informative_centered_mean").all())
        self.assertTrue(centered["warning"].eq(CENTERED_HV_WARNING).all())
        hv = build_hv_summary(changes, samples)
        self.assertTrue(hv.loc[hv["metric"].isin(["H_expr_mean", "V_expr_mean"]), "warning"].eq(CENTERED_HV_WARNING).all())


class Box3VisualizationAndSummaryTests(unittest.TestCase):
    def _config(self) -> ComparativeConfig:
        return ComparativeConfig(
            mode="pairwise", reference="A", target="B", c_genes=["C1"], s_genes=["S1"],
            enable_h_expr=False, enable_v_expr=False, bootstrap_iterations=100,
        )

    def test_pairwise_regime_uses_card_and_paired_uses_matrix(self) -> None:
        direct = pd.DataFrame([{
            "comparison_id": "p1", "comparison_basis": "matched_sample_operational_change",
            "reference_regime": "Type_A_candidate", "target_regime": "Type_B_candidate",
            "regime_transition": "Type_A_candidate -> Type_B_candidate", "reference_confidence": 0.4,
            "target_confidence": 0.5, "transition_confidence_flag": "adequate_for_descriptive_comparison",
        }])
        with tempfile.TemporaryDirectory() as tmp:
            pair_path = plot_regime_transitions(direct, Path(tmp) / "pair.png", self._config(), effective_mode="pairwise")
            pair_meta = json.loads(pair_path.with_suffix(".png.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(pair_meta["visualization_mode"], "pairwise_transition_card")
            group_path = plot_regime_transitions(
                pd.concat([direct, direct.assign(comparison_id="p2")], ignore_index=True),
                Path(tmp) / "group.png", self._config(), effective_mode="paired",
            )
            group_meta = json.loads(group_path.with_suffix(".png.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(group_meta["visualization_mode"], "matched_group_3x3_transition_matrix")

    def test_summary_prioritizes_normalized_metrics_and_has_no_mechanism_claims(self) -> None:
        samples = _base_samples()
        deltas = compute_delta_metrics(samples, _matches())
        changes = build_metric_change_table(deltas)
        warnings = build_scale_warnings(changes)
        transitions = pd.DataFrame([{
            "comparison_basis": "matched_sample_operational_change",
            "regime_transition": "Type_A_candidate -> Type_A_candidate",
            "reference_regime": "Type_A_candidate", "target_regime": "Type_A_candidate",
            "transition_confidence_flag": "adequate_for_descriptive_comparison",
        }])
        text = rules_based_summary(
            samples, deltas, transitions, self._config(), "pairwise",
            metric_change_table=changes, scale_warnings=warnings, hv_summary=pd.DataFrame(),
        )
        self.assertIn("normalized diffuse-component density", text)
        self.assertIn("raw diffuse-component count", text)
        forbidden = [
            "treatment response confirmed", "tumor architecture improved", "immune activation occurred",
            "mechanism proven", "diagnostic significance", "predictive significance",
        ]
        self.assertFalse(any(term in text.lower() for term in forbidden))
        self.assertIn("no direct spatial registration", text.lower())


if __name__ == "__main__":
    unittest.main()
