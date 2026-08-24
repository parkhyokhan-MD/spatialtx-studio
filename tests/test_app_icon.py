from __future__ import annotations

import unittest

from PIL import Image

from spatialtx_desktop.app import APP_ICON_ICO, APP_ICON_PNG


class ApplicationIconTests(unittest.TestCase):
    def test_packaged_png_has_square_rgba_canvas_and_transparency(self) -> None:
        self.assertTrue(APP_ICON_PNG.is_file())
        with Image.open(APP_ICON_PNG) as image:
            rgba = image.convert("RGBA")
            self.assertEqual(rgba.width, rgba.height)
            alpha_min, alpha_max = rgba.getchannel("A").getextrema()
            self.assertEqual(alpha_min, 0)
            self.assertEqual(alpha_max, 255)

    def test_windows_icon_contains_small_and_large_frames(self) -> None:
        self.assertTrue(APP_ICON_ICO.is_file())
        with Image.open(APP_ICON_ICO) as image:
            sizes = set(image.ico.sizes())
        self.assertTrue({(16, 16), (32, 32), (48, 48), (256, 256)}.issubset(sizes))


if __name__ == "__main__":
    unittest.main()
