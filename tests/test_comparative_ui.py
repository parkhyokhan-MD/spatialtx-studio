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
            self.assertEqual(app.title(), "SpatialTX Studio Desktop v0.5-beta")
            self.assertFalse(app.comparative_analysis_panel.busy)
            panel = app.comparative_analysis_panel
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


if __name__ == "__main__":
    unittest.main()
