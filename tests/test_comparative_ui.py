from __future__ import annotations

import tkinter as tk
import unittest
from pathlib import Path

import pandas as pd

from spatialtx_desktop.app import SpatialTXDesktop
from spatialtx_desktop.comparative.models import ComparativeRunResult
from spatialtx_desktop.comparative_ui import ComparativeAnalysisPanel, _sample_display_labels


class ComparativeUITests(unittest.TestCase):
    def test_pairwise_sample_labels_keep_filename_visible_and_disambiguate_duplicates(self) -> None:
        paths = [
            Path("C:/a/very/long/directory/reference/sample.h5ad"),
            Path("C:/another/equally/long/directory/target/sample.h5ad"),
            Path("C:/third/long/directory/unique_target.h5ad"),
        ]
        labels, path_by_label = _sample_display_labels(paths)

        self.assertEqual(len(labels), len(set(labels)))
        self.assertTrue(labels[0].startswith("sample.h5ad"))
        self.assertTrue(labels[1].startswith("sample.h5ad"))
        self.assertEqual(labels[2], "unique_target.h5ad")
        self.assertEqual(path_by_label[labels[0]], str(paths[0]))
        self.assertEqual(path_by_label[labels[1]], str(paths[1]))

    def test_desktop_mounts_comparative_tab_without_changing_main_mapper(self) -> None:
        try:
            app = SpatialTXDesktop()
        except tk.TclError as exc:
            self.skipTest(f"Tk display unavailable: {exc}")
        try:
            app.withdraw()
            app.update_idletasks()
            labels = [app.right_tabs.tab(tab_id, "text") for tab_id in app.right_tabs.tabs()]
            self.assertIn("Main Mapper", labels)
            self.assertIn("Comparative Analysis", labels)
            self.assertIsInstance(app.comparative_analysis_panel, ComparativeAnalysisPanel)
            self.assertEqual(app.title(), "SpatialTX Studio Desktop v0.65-dev")
            self.assertIsNotNone(app._app_icon_photo)
            self.assertIn("png", app.app_icon_status)
            self.assertFalse(app.comparative_analysis_panel.busy)
            panel = app.comparative_analysis_panel
            self.assertFalse(panel.multi_pair_panel.busy)
            rules_panel = panel.multi_pair_panel
            rules = rules_panel.rules_text.get("1.0", "end")
            multi_pair_tabs = [
                rules_panel.result_tabs.tab(tab_id, "text") for tab_id in rules_panel.result_tabs.tabs()
            ]
            self.assertIn("1  Balance change", multi_pair_tabs)
            self.assertIn("2  Spatial organization", multi_pair_tabs)
            self.assertIn("3  Specimen reliability", multi_pair_tabs)
            self.assertIn("v0.65 Reliability", multi_pair_tabs)
            self.assertIn("Pair interpretation", multi_pair_tabs)
            self.assertIn("Delta = Post - Pre", rules)
            self.assertIn("THREE SEPARATE RESULT LAYERS", rules)
            self.assertIn("not combined into a single response", rules)
            self.assertIn("TRANSPARENT QUALITATIVE CHANGE LABELS", rules)
            self.assertIn("PAIR-ID SAFETY CHECK", rules)
            self.assertIn("Low: any primary Low reason", rules)
            self.assertIn("Spot-count fold ratio: 1.50 / 2.50", rules)
            self.assertIn("cannot by itself produce Low", rules)
            self.assertIn("do not establish treatment response", rules)
            self.assertIn("OPTIONAL v0.65 RELIABILITY LAYER", rules)
            self.assertFalse(rules_panel._reliability_config().enabled)
            rules_panel.reliability_layer_var.set(True)
            rules_panel.reliability_classification_var.set(True)
            with self.assertRaisesRegex(ValueError, "requires explicit"):
                rules_panel._reliability_config()
            rules_panel.reliability_activity_threshold_var.set("0.5")
            rules_panel.reliability_direction_threshold_var.set("0.4")
            reliability_config = rules_panel._reliability_config()
            self.assertTrue(reliability_config.enabled)
            self.assertTrue(reliability_config.classification_enabled)
            rules_panel._show_rules_tab()
            self.assertEqual(rules_panel.result_tabs.select(), str(rules_panel.rules_tab))
            panel._load_result(ComparativeRunResult(
                run_dir=Path.cwd(),
                sample_metrics=pd.DataFrame([{
                    "sample_id": "a", "group": "A", "QC_flag": "PASS", "regime_label": "Type_A_candidate",
                }]),
                delta_metrics=pd.DataFrame([{
                    "comparison_id": "x", "delta_metric": "delta_R", "reference_value": 0,
                    "target_value": 1, "delta": 1, "status": "ok",
                }]),
                group_statistics=pd.DataFrame(),
                regime_transitions=pd.DataFrame([{
                    "comparison_id": "x", "reference_regime": "Type_A_candidate",
                    "target_regime": "Type_A_candidate", "regime_transition": "stable",
                    "transition_confidence_flag": "adequate_for_descriptive_comparison",
                }]),
                warnings=pd.DataFrame(columns=["sample_id", "message"]),
                run_manifest=pd.DataFrame([{"status": "ok"}]),
                figures=[], summary_text="Exploratory summary", effective_mode="pairwise",
                metric_change_table=pd.DataFrame([{
                    "comparison_id": "x", "display_name": "R balance-field mean", "reference_value": 0,
                    "target_value": 1, "raw_delta": 1, "normalized_delta": float("nan"),
                    "symmetric_percent_change": 200, "scale_sensitive": False,
                    "observational_only": False, "warning": "",
                }]),
                scale_warnings=pd.DataFrame([{
                    "severity": "caution", "message": "Caution: substantial sample-scale difference."
                }]),
            ))
            self.assertIn("Caution", panel.scale_banner_var.get())
            self.assertIn("observational only", panel.hv_notice_var.get())
        finally:
            app.destroy()

    def test_multi_pair_rows_one_two_and_six_run_as_three_pairs(self) -> None:
        try:
            app = SpatialTXDesktop()
        except tk.TclError as exc:
            self.skipTest(f"Tk display unavailable: {exc}")
        try:
            app.withdraw()
            panel = app.comparative_analysis_panel.multi_pair_panel
            self.assertEqual(len(panel.pair_vars), 6)
            for variables in panel.pair_vars:
                variables["label"].set("")
                variables["pre"].set("")
                variables["post"].set("")
            for index in (0, 1, 5):
                panel.pair_vars[index]["label"].set(f"P{index + 1}")
                panel.pair_vars[index]["pre"].set(f"C:/data/pre_{index + 1}.h5ad")
                panel.pair_vars[index]["post"].set(f"C:/data/post_{index + 1}.h5ad")
            pairs = panel._pairs()
            self.assertEqual([pair.label for pair in pairs], ["P1", "P2", "P6"])
            panel.pair_vars[2]["pre"].set("C:/data/incomplete.h5ad")
            with self.assertRaisesRegex(ValueError, "partially filled"):
                panel._pairs()
        finally:
            app.destroy()


if __name__ == "__main__":
    unittest.main()
