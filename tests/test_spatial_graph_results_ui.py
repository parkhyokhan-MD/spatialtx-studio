from __future__ import annotations

import tempfile
import tkinter as tk
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, UnidentifiedImageError

from spatialtx_desktop import __version__
from spatialtx_desktop.advanced_analysis_ui import AdvancedAnalysisPanel
from spatialtx_desktop.graph.runner import MODULE_VERSION
from spatialtx_desktop.spatial_graph_results_ui import (
    FIGURE_PRIORITY,
    SpatialGraphResultsPanel,
    classify_figure_filename,
    discover_figures,
    load_context_summary,
    load_image_for_display,
    match_sample_from_figure,
    parse_context_summary_frame,
    parse_manifest,
    parse_warning_value,
    sample_name_from_known_figure,
)


class FigureDiscoveryTests(unittest.TestCase):
    def test_application_version_changes_but_graph_schema_version_does_not(self) -> None:
        self.assertEqual(f"v{__version__}", "v0.6-beta")
        self.assertEqual(MODULE_VERSION, "0.4")

    def test_all_known_suffixes_are_classified_without_short_suffix_collisions(self) -> None:
        expected = {
            "PDAC_P1_fresh_H_expr_map.png": "H_expr hypoxia-associated expression field",
            "PDAC_P1_fresh_H_expr_unsmoothed_map.png": "H_expr unsmoothed",
            "PDAC_P1_fresh_H_expr_smoothed_map.png": "H_expr graph-smoothed",
            "PDAC_P1_fresh_V_expr_map.png": "V_expr endothelial/angiogenic expression proxy",
            "PDAC_P1_fresh_V_expr_unsmoothed_map.png": "V_expr unsmoothed",
            "PDAC_P1_fresh_V_expr_smoothed_map.png": "V_expr graph-smoothed",
            "PDAC_P1_fresh_H_V_association_map.png": "H/V joint high-state map",
            "PDAC_P1_fresh_neighborhood_enrichment_heatmap.png": "Neighborhood enrichment heatmap",
            "GSM4284322_P9_ST_rep1_graph_qc.png": "Graph QC",
        }
        for filename, label in expected.items():
            with self.subTest(filename=filename):
                self.assertEqual(classify_figure_filename(filename)[0], label)

    def test_longest_manifest_sample_prefix_wins(self) -> None:
        samples = ["PDAC_P1", "PDAC_P1_fresh", "GSM4284322_P9_ST_rep1"]
        self.assertEqual(
            match_sample_from_figure("PDAC_P1_fresh_H_expr_map.png", samples),
            "PDAC_P1_fresh",
        )
        self.assertEqual(
            match_sample_from_figure("GSM4284322_P9_ST_rep1_graph_qc.png", samples),
            "GSM4284322_P9_ST_rep1",
        )
        self.assertEqual(
            sample_name_from_known_figure("GSM4284322_P9_ST_rep1_graph_qc.png"),
            "GSM4284322_P9_ST_rep1",
        )
        self.assertIsNone(sample_name_from_known_figure("unknown_plot.png"))

    def test_priority_and_unknown_png_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            figures = run_dir / "figures"
            figures.mkdir()
            names = [
                "PDAC_P1_graph_qc.png",
                "PDAC_P1_custom_diagnostic.png",
                "PDAC_P1_V_expr_map.png",
                "PDAC_P1_H_expr_map.png",
                "PDAC_P1_H_V_association_map.png",
                "PDAC_P1_neighborhood_enrichment_heatmap.png",
            ]
            for name in names:
                Image.new("RGB", (4, 3), "white").save(figures / name)
            records = discover_figures(run_dir, ["PDAC_P1"])["PDAC_P1"]
        self.assertEqual(records[0].suffix, "_H_V_association_map.png")
        self.assertEqual(records[1].suffix, "_H_expr_map.png")
        self.assertEqual(records[2].suffix, "_V_expr_map.png")
        self.assertEqual(records[3].suffix, "_neighborhood_enrichment_heatmap.png")
        self.assertEqual(records[4].suffix, "_graph_qc.png")
        self.assertEqual(records[-1].label, "custom_diagnostic")
        self.assertEqual(FIGURE_PRIORITY["_H_V_association_map.png"], 0)


class ManifestParsingTests(unittest.TestCase):
    def test_success_failure_cancel_and_underscore_samples(self) -> None:
        manifest = pd.DataFrame([
            {"sample": "PDAC_P1_fresh", "status": "ok"},
            {"sample": "GBM_5_2", "status": "error: graph failed"},
            {"sample": "SP4", "status": "cancelled"},
        ])
        parsed = parse_manifest(manifest)
        self.assertEqual(parsed.successful_samples, ("PDAC_P1_fresh",))
        self.assertEqual(parsed.failed_samples[0], ("GBM_5_2", "error: graph failed"))
        self.assertEqual(parsed.failed_samples[1], ("SP4", "cancelled"))
        self.assertEqual(parsed.total_rows, 3)

    def test_missing_status_and_empty_manifest_are_safe(self) -> None:
        missing = parse_manifest(pd.DataFrame([{"sample": "PDAC_P1"}]))
        self.assertEqual(missing.successful_samples, ())
        self.assertEqual(missing.failed_samples, (("PDAC_P1", "status missing"),))
        self.assertIn("no status", missing.warning)
        empty = parse_manifest(pd.DataFrame())
        self.assertEqual(empty.total_rows, 0)
        self.assertIn("empty", empty.warning)


class ContextSummaryTests(unittest.TestCase):
    def test_h_and_v_ok_string_counts_nan_coverage_and_warning_formats(self) -> None:
        table = pd.DataFrame([
            {
                "field": "H_expr", "status": "ok", "matched_gene_count": "10",
                "requested_gene_count": "12", "coverage_fraction": "0.833",
                "score_method": "z_score_mean", "smoothing_method": "none",
                "high_state_fraction": "0.2", "warnings": np.nan,
            },
            {
                "field": "V_expr", "status": "ok", "matched_gene_count": "4",
                "requested_gene_count": "12", "coverage_fraction": np.nan,
                "score_method": "raw_mean", "smoothing_method": "graph_mean",
                "high_state_fraction": "0.25", "warnings": "first; second",
            },
        ])
        summary = parse_context_summary_frame(table)
        self.assertTrue(summary.available)
        self.assertEqual([record.field for record in summary.records], ["H_expr", "V_expr"])
        self.assertEqual(summary.records[0].coverage_fraction, "83.3%")
        self.assertEqual(summary.records[0].warnings, ())
        self.assertEqual(summary.records[1].coverage_fraction, "Unavailable")
        self.assertEqual(summary.records[1].warnings, ("first", "second"))

    def test_skipped_qc_one_field_and_gene_list_count_fallback(self) -> None:
        table = pd.DataFrame([{
            "field": "H_expr", "status": "skipped_qc", "requested_genes": "CA9;VEGFA;LDHA",
            "matched_genes": '["CA9"]', "coverage_fraction": 1 / 3,
            "warnings": '["low H_expr gene coverage"]',
        }])
        summary = parse_context_summary_frame(table)
        self.assertTrue(summary.available)
        self.assertEqual(len(summary.records), 1)
        self.assertEqual(summary.records[0].requested_gene_count, "3")
        self.assertEqual(summary.records[0].matched_gene_count, "1")
        self.assertEqual(summary.records[0].warnings, ("low H_expr gene coverage",))

    def test_v_only_empty_warnings_and_missing_file(self) -> None:
        summary = parse_context_summary_frame(pd.DataFrame([{
            "field": "V_expr", "status": "skipped_qc", "warnings": "",
        }]))
        self.assertEqual(summary.records[0].field, "V_expr")
        self.assertEqual(summary.records[0].warnings, ())
        with tempfile.TemporaryDirectory() as tmp:
            missing = load_context_summary(Path(tmp) / "missing.csv")
        self.assertFalse(missing.available)
        self.assertIn("unavailable", missing.message.lower())

    def test_warning_normalization_handles_list_json_semicolon_and_nan(self) -> None:
        self.assertEqual(parse_warning_value(["a", "b"]), ("a", "b"))
        self.assertEqual(parse_warning_value('["a", "b"]'), ("a", "b"))
        self.assertEqual(parse_warning_value("['a', 'b']"), ("a", "b"))
        self.assertEqual(parse_warning_value("a; b"), ("a", "b"))
        self.assertEqual(parse_warning_value(np.nan), ())


class ImageHandlingTests(unittest.TestCase):
    def test_normal_and_large_images_fit_without_aspect_distortion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            normal = Path(tmp) / "normal.png"
            large = Path(tmp) / "large.png"
            Image.new("RGB", (400, 200), "red").save(normal)
            Image.new("RGB", (4000, 2000), "blue").save(large)
            display, original = load_image_for_display(normal, 200, 200)
            big_display, big_original = load_image_for_display(large, 300, 120)
        self.assertEqual(original, (400, 200))
        self.assertEqual(display.size, (200, 100))
        self.assertEqual(big_original, (4000, 2000))
        self.assertLessEqual(big_display.width, 300)
        self.assertLessEqual(big_display.height, 120)
        self.assertAlmostEqual(big_display.width / big_display.height, 2.0, places=1)

    def test_corrupt_deleted_and_missing_images_raise_safe_read_errors(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            corrupt = Path(tmp) / "corrupt.png"
            corrupt.write_bytes(b"not a png")
            with self.assertRaises((UnidentifiedImageError, OSError)):
                load_image_for_display(corrupt, 100, 100)
            deleted = Path(tmp) / "deleted.png"
            Image.new("RGB", (4, 4)).save(deleted)
            deleted.unlink()
            with self.assertRaises(FileNotFoundError):
                load_image_for_display(deleted, 100, 100)


class SpatialGraphResultsUISmokeTests(unittest.TestCase):
    def test_hidden_panel_load_navigation_sample_change_and_empty_state(self) -> None:
        try:
            root = tk.Tk()
        except tk.TclError as exc:
            self.skipTest(f"Tk display unavailable: {exc}")
        root.withdraw()
        try:
            panel = SpatialGraphResultsPanel(root, lambda _path: None)
            panel.pack(fill="both", expand=True)
            with tempfile.TemporaryDirectory() as tmp:
                run_dir = Path(tmp)
                figures = run_dir / "figures"
                context = run_dir / "context_fields"
                figures.mkdir()
                context.mkdir()
                for filename in (
                    "PDAC_P1_H_expr_map.png",
                    "PDAC_P1_H_V_association_map.png",
                    "PDAC_P1_fresh_graph_qc.png",
                ):
                    Image.new("RGB", (40, 30), "white").save(figures / filename)
                pd.DataFrame([{
                    "field": "H_expr", "status": "ok", "matched_gene_count": 2,
                    "requested_gene_count": 3, "coverage_fraction": .667,
                }]).to_csv(context / "PDAC_P1_context_field_summary.csv", index=False)
                manifest = pd.DataFrame([
                    {"sample": "PDAC_P1", "status": "ok"},
                    {"sample": "PDAC_P1_fresh", "status": "ok"},
                    {"sample": "failed_sample", "status": "error: test"},
                ])
                panel.load_run(run_dir, manifest)
                root.update_idletasks()
                self.assertEqual(panel.sample_var.get(), "PDAC_P1")
                self.assertIn("H/V joint", panel.figure_var.get())
                panel._step_figure(1)
                self.assertIn("H_expr", panel.figure_var.get())
                panel.sample_var.set("PDAC_P1_fresh")
                panel._sample_changed()
                self.assertEqual(len(panel.current_figures), 1)
                self.assertIn("Graph QC", panel.figure_var.get())
                panel.load_run(run_dir, pd.DataFrame([{"sample": "only_failure", "status": "cancelled"}]))
                self.assertEqual(panel.sample_var.get(), "")
                self.assertIn("No successful", panel.status_var.get())
        finally:
            root.destroy()

    def test_graph_done_event_loads_viewer_and_selects_results_tab(self) -> None:
        try:
            root = tk.Tk()
        except tk.TclError as exc:
            self.skipTest(f"Tk display unavailable: {exc}")
        root.withdraw()
        try:
            panel = AdvancedAnalysisPanel(
                root,
                get_samples=lambda: [],
                get_genes=lambda: (["CD8A"], ["COL1A1"]),
                get_quantiles=lambda: (.8, .8, .6),
                get_output=lambda: ".",
            )
            panel.pack(fill="both", expand=True)
            with tempfile.TemporaryDirectory() as tmp:
                run_dir = Path(tmp)
                figures = run_dir / "figures"
                figures.mkdir()
                Image.new("RGB", (40, 30), "white").save(figures / "sample_H_expr_map.png")
                manifest = pd.DataFrame([{"sample": "sample", "status": "ok"}])
                panel.events.put(("graph_done", run_dir, manifest))
                panel._poll()
                root.update_idletasks()
                self.assertEqual(panel.tabs.select(), str(panel.spatial_graph_results_tab))
                self.assertEqual(panel.spatial_graph_results_panel.sample_var.get(), "sample")
                self.assertIn("Loaded 1 figures", panel.status.cget("text"))
        finally:
            root.destroy()


if __name__ == "__main__":
    unittest.main()
