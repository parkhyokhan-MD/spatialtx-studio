from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from spatialtx_desktop.comparative.models import ComparativeConfig, SampleRecord
from spatialtx_desktop.comparative.runner import run_comparative_analysis
from tests.comparative_fixtures import C_GENES, S_GENES, write_comparative_h5ad


class ComparativeExportTests(unittest.TestCase):
    def test_required_pairwise_exports_figures_reports_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = write_comparative_h5ad(root / "a.h5ad", pattern="localized", include_context=False)
            b = write_comparative_h5ad(root / "b.h5ad", pattern="diffuse", include_context=False)
            result = run_comparative_analysis(
                [SampleRecord("a", a, "A"), SampleRecord("b", b, "B")],
                ComparativeConfig(
                    mode="pairwise", reference="A", target="B", c_genes=C_GENES, s_genes=S_GENES,
                    enable_h_expr=False, enable_v_expr=False, bootstrap_iterations=100,
                ),
                root / "out",
            )
            required = [
                "comparative_sample_metrics.csv", "comparative_delta_metrics.csv",
                "comparative_group_statistics.csv", "comparative_regime_transitions.csv",
                "comparative_warnings.csv", "comparative_parameters.json",
                "comparative_run_manifest.csv", "comparative_summary_report.html",
                "comparative_summary_report.pdf", "sample_A_summary.csv", "sample_B_summary.csv",
                "pairwise_delta_summary.csv",
                "comparative_sample_scale.csv", "comparative_metric_change_table.csv",
                "comparative_normalized_metrics.csv", "comparative_relative_changes.csv",
                "comparative_scale_warnings.csv", "comparative_HV_summary.csv",
            ]
            self.assertTrue(all((result.run_dir / name).is_file() for name in required))
            self.assertTrue((result.run_dir / "comparative_logs" / "run.log").is_file())
            figures = list((result.run_dir / "comparative_figures").glob("*.png"))
            self.assertGreaterEqual(len(figures), 5)
            expected_figures = {
                "comparative_summary_card.png", "comparative_program_score_changes.png",
                "comparative_transition_changes.png", "comparative_graph_changes.png",
                "comparative_topology_raw_counts.png", "comparative_topology_normalized.png",
                "comparative_sample_scale.png", "comparative_relative_changes.png",
                "comparative_standardized_heatmap.png", "comparative_regime_transitions.png",
                "comparative_side_by_side_maps.png", "comparative_HV_summary.png",
                "comparative_metric_changes.png", "comparative_H_V_context.png",
            }
            self.assertTrue(expected_figures.issubset({path.name for path in figures}))
            self.assertTrue(all(path.with_suffix(path.suffix + ".metadata.json").is_file() for path in figures))
            parameters = json.loads((result.run_dir / "comparative_parameters.json").read_text(encoding="utf-8"))
            self.assertFalse(parameters["spotwise_subtraction_performed"])
            self.assertTrue(parameters["H_V_observational_only"])
            self.assertFalse(parameters["candidate_discovery_added"])
            self.assertIn("software_environment", parameters)
            metrics = pd.read_csv(result.run_dir / "comparative_sample_metrics.csv")
            self.assertTrue(metrics["input_sha256"].str.len().eq(64).all())
            self.assertTrue(metrics["spatialtx_version"].str.len().gt(0).all())
            changes = pd.read_csv(result.run_dir / "comparative_metric_change_table.csv")
            for column in (
                "metric_name", "display_name", "category", "unit", "reference_value", "target_value",
                "raw_delta", "percent_change", "symmetric_percent_change", "normalization_denominator",
                "normalization_status", "scale_sensitive", "observational_only", "interpretation_flag", "warning",
            ):
                self.assertIn(column, changes)
            old_delta = pd.read_csv(result.run_dir / "comparative_delta_metrics.csv")
            old_row = old_delta.loc[old_delta["metric"].eq("transition_burden_score")].iloc[0]
            new_row = changes.loc[changes["metric_name"].eq("transition_burden_score")].iloc[0]
            self.assertEqual(old_row["delta"], new_row["raw_delta"])
            regime_meta = json.loads(
                (result.run_dir / "comparative_figures" / "comparative_regime_transitions.png.metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(regime_meta["visualization_mode"], "pairwise_transition_card")
            heatmap_meta = json.loads(
                (result.run_dir / "comparative_figures" / "comparative_standardized_heatmap.png.metadata.json").read_text(encoding="utf-8")
            )
            self.assertIn("within-run", heatmap_meta["visualization_transform"])
            side_meta = json.loads(
                (result.run_dir / "comparative_figures" / "comparative_side_by_side_maps.png.metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(side_meta["fill_definition"], "R = C - S")
            self.assertIn("transition masks", side_meta["outline_definition"])
            for filename in (
                "comparative_delta_metrics.csv", "comparative_group_statistics.csv",
                "comparative_regime_transitions.csv", "comparative_warnings.csv",
                "comparative_run_manifest.csv", "comparative_sample_scale.csv",
                "comparative_metric_change_table.csv", "comparative_normalized_metrics.csv",
                "comparative_relative_changes.csv", "comparative_scale_warnings.csv",
                "comparative_HV_summary.csv",
            ):
                table = pd.read_csv(result.run_dir / filename)
                for column in (
                    "spatialtx_version", "input_paths_json", "input_hashes_json",
                    "C_gene_program_json", "S_gene_program_json", "graph_settings_json",
                    "thresholds_json", "warnings_json", "failed_samples_json",
                    "software_environment_json",
                ):
                    self.assertIn(column, table.columns, msg=f"{filename}: {column}")

    def test_missing_h_v_continues_with_warning(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = write_comparative_h5ad(root / "a.h5ad", include_context=False)
            b = write_comparative_h5ad(root / "b.h5ad", pattern="diffuse", include_context=False)
            result = run_comparative_analysis(
                [SampleRecord("a", a, "A"), SampleRecord("b", b, "B")],
                ComparativeConfig(
                    mode="pairwise", reference="A", target="B", c_genes=C_GENES, s_genes=S_GENES,
                    bootstrap_iterations=100,
                ),
                root / "out",
            )
            self.assertTrue(result.run_manifest["status"].eq("ok").all())
            self.assertFalse(result.sample_metrics["H_expr_available"].astype(bool).any())
            self.assertFalse(result.sample_metrics["V_expr_available"].astype(bool).any())
            self.assertTrue(result.warnings["message"].str.contains("H_expr|V_expr", regex=True).any())


if __name__ == "__main__":
    unittest.main()
