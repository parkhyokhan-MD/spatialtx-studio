from __future__ import annotations

import tkinter as tk
import unittest

from spatialtx_desktop.advanced import scan_pre_post_pairs
from spatialtx_desktop.advanced_ui import AdvancedToolsPanel


def _widget_texts(widget) -> list[str]:
    texts: list[str] = []
    try:
        text = widget.cget("text")
    except tk.TclError:
        text = ""
    if text:
        texts.append(str(text))
    for child in widget.winfo_children():
        texts.extend(_widget_texts(child))
    return texts


class AdvancedToolsUITests(unittest.TestCase):
    def test_a1_scanner_is_hidden_and_a3_is_explicitly_non_spatial(self) -> None:
        try:
            root = tk.Tk()
        except tk.TclError as exc:
            self.skipTest(f"Tk display unavailable: {exc}")
        try:
            root.withdraw()
            panel = AdvancedToolsPanel(root)
            panel.pack()
            root.update_idletasks()
            texts = _widget_texts(panel)

            self.assertFalse(any("A1" in text or "Scan pre/post pairs" in text for text in texts))
            self.assertTrue(any("non-spatial" in text.lower() for text in texts))
            self.assertIn("Reference H5AD", texts)
            self.assertIn("Target H5AD", texts)
            self.assertTrue(callable(scan_pre_post_pairs))
        finally:
            root.destroy()


if __name__ == "__main__":
    unittest.main()
