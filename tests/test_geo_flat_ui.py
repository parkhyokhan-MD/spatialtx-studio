from __future__ import annotations

import csv
import tempfile
import tkinter as tk
import unittest
from pathlib import Path

from spatialtx_desktop.app import SpatialTXDesktop
from spatialtx_desktop.importers.geo_flat_panel import GeoFlatImportPanel


def descendants(widget):
    for child in widget.winfo_children():
        yield child
        yield from descendants(child)


class GeoFlatUITests(unittest.TestCase):
    def test_geo_panel_mounts_and_confirmed_manifest_handoff_does_not_run(self) -> None:
        try:
            app = SpatialTXDesktop()
        except tk.TclError as exc:
            self.skipTest(f"Tk display unavailable: {exc}")
        app.withdraw()
        try:
            geo_panels = [
                widget for widget in descendants(app.import_convert_panel)
                if isinstance(widget, GeoFlatImportPanel)
            ]
            self.assertEqual(len(geo_panels), 1)
            self.assertFalse(geo_panels[0].recursive_var.get())
            self.assertEqual(geo_panels[0].samples, [])
            self.assertIn("Nothing is converted automatically", geo_panels[0].status_var.get())
            with tempfile.TemporaryDirectory() as tmp:
                manifest = Path(tmp) / "confirmed.csv"
                with manifest.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=[
                        "sample_id", "file_path", "group", "pair_id", "condition", "batch", "notes", "pairing_source",
                    ])
                    writer.writeheader()
                    writer.writerow({"sample_id": "a", "file_path": "a.h5ad", "group": "pre", "pairing_source": "user_confirmed"})
                    writer.writerow({"sample_id": "b", "file_path": "b.h5ad", "group": "post", "pairing_source": "user_confirmed"})
                app._open_comparative_manifest(manifest)
                self.assertEqual(app.comparative_analysis_panel.mode_var.get(), "manifest_batch")
                self.assertEqual(app.comparative_analysis_panel.manifest_var.get(), str(manifest))
                self.assertEqual(app.comparative_analysis_panel.reference_var.get(), "pre")
                self.assertEqual(app.comparative_analysis_panel.target_var.get(), "post")
                self.assertFalse(app.comparative_analysis_panel.busy)
        finally:
            app.destroy()


if __name__ == "__main__":
    unittest.main()
