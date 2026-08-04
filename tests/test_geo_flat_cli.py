from __future__ import annotations

import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path

import pandas as pd

from geo_flat_cli import main
from tests.geo_flat_fixtures import write_flat_sample


class GeoFlatCLITests(unittest.TestCase):
    def test_list_samples_needs_no_output_and_does_not_convert(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source with spaces"
            write_flat_sample(source, "GSM1_sample_1_pre", include_images=False)
            output = StringIO()
            with redirect_stdout(output):
                status = main(["--input-dir", str(source), "--list-samples"])
            self.assertEqual(status, 0)
            self.assertIn("GSM1_sample_1_pre", output.getvalue())
            self.assertEqual({path.name for path in source.iterdir()}, {
                "GSM1_sample_1_pre_matrix.mtx",
                "GSM1_sample_1_pre_barcodes.tsv",
                "GSM1_sample_1_pre_features.tsv",
                "GSM1_sample_1_pre_tissue_positions.csv",
                "GSM1_sample_1_pre_scalefactors_json.json",
            })

    def test_selected_conversion_and_unconfirmed_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, destination = root / "input", root / "output"
            write_flat_sample(source, "GSM1_sample_1_pre", include_images=False)
            write_flat_sample(source, "GSM2_sample_1_post", include_images=False)
            manifest = root / "draft.csv"
            with redirect_stdout(StringIO()):
                status = main([
                    "--input-dir", str(source),
                    "--samples", "GSM1_sample_1_pre",
                    "--output-dir", str(destination),
                    "--write-comparative-manifest", str(manifest),
                ])
            self.assertEqual(status, 0)
            self.assertTrue((destination / "GSM1_sample_1_pre.h5ad").is_file())
            self.assertFalse((destination / "GSM2_sample_1_post.h5ad").exists())
            table = pd.read_csv(manifest, keep_default_na=False)
            self.assertEqual(table.loc[0, "pairing_source"], "filename_inference_unconfirmed")
            self.assertEqual(table.loc[0, "group"], "")
            self.assertEqual(table.loc[0, "pair_id"], "")

    def test_explicit_cli_confirmation_marks_manifest_user_confirmed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, destination = root / "input", root / "output"
            write_flat_sample(source, "GSM1_sample_1_pre", include_images=False)
            manifest = root / "confirmed.csv"
            with redirect_stdout(StringIO()):
                status = main([
                    "--input-dir", str(source),
                    "--output-dir", str(destination),
                    "--write-comparative-manifest", str(manifest),
                    "--confirm-filename-pairings",
                ])
            self.assertEqual(status, 0)
            table = pd.read_csv(manifest, keep_default_na=False)
            self.assertEqual(table.loc[0, "pairing_source"], "user_confirmed")
            self.assertEqual(table.loc[0, "group"], "pre")
            self.assertEqual(table.loc[0, "pair_id"], "sample_1")

    def test_conversion_requires_output_and_unknown_sample_is_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "input"
            write_flat_sample(source, "GSM1_sample_1_pre", include_images=False)
            errors = StringIO()
            with redirect_stderr(errors):
                self.assertEqual(main(["--input-dir", str(source)]), 1)
                self.assertEqual(main([
                    "--input-dir", str(source), "--output-dir", str(Path(tmp) / "out"),
                    "--samples", "not_detected",
                ]), 1)
            self.assertIn("required", errors.getvalue())
            self.assertIn("not detected", errors.getvalue())


if __name__ == "__main__":
    unittest.main()
