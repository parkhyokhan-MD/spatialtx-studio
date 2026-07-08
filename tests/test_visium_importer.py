from __future__ import annotations

import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse

from spatialtx_desktop.importers.validate_h5ad import validate_h5ad
from spatialtx_desktop.importers.visium_to_h5ad import (
    convert_visium_to_h5ad,
    detect_visium_sample,
)


class VisiumImporterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.sample = self.root / "sample_01"
        self.spatial = self.sample / "spatial"
        self.spatial.mkdir(parents=True)
        self._write_10x_h5(self.sample / "filtered_feature_bc_matrix.h5")
        self._write_spatial_files()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def _write_10x_h5(path: Path) -> None:
        # 10x HDF5 stores a feature x barcode CSC matrix; importer transposes it.
        feature_by_spot = sparse.csc_matrix(np.asarray([
            [1, 0, 3],
            [0, 2, 4],
        ], dtype=np.int32))
        with h5py.File(path, "w") as handle:
            matrix = handle.create_group("matrix")
            matrix.create_dataset("data", data=feature_by_spot.data)
            matrix.create_dataset("indices", data=feature_by_spot.indices)
            matrix.create_dataset("indptr", data=feature_by_spot.indptr)
            matrix.create_dataset("shape", data=np.asarray(feature_by_spot.shape, dtype=np.int64))
            matrix.create_dataset("barcodes", data=np.asarray(["BC1", "BC2", "BC3"], dtype="S"))
            features = matrix.create_group("features")
            features.create_dataset("id", data=np.asarray(["ENSG1", "ENSG2"], dtype="S"))
            features.create_dataset("name", data=np.asarray(["GENE", "GENE"], dtype="S"))
            features.create_dataset("feature_type", data=np.asarray(["Gene Expression", "Gene Expression"], dtype="S"))

    def _write_spatial_files(self) -> None:
        positions = pd.DataFrame({
            "barcode": ["BC1", "BC2", "BC3"],
            "in_tissue": [1, 1, 1],
            "array_row": [0, 1, 2],
            "array_col": [2, 1, 0],
            "pxl_row_in_fullres": [100, 200, 300],
            "pxl_col_in_fullres": [110, 210, 310],
        })
        positions.to_csv(self.spatial / "tissue_positions.csv.gz", index=False, compression="gzip")
        with gzip.open(self.spatial / "scalefactors_json.json.gz", "wt", encoding="utf-8") as handle:
            json.dump({"tissue_hires_scalef": 0.5, "spot_diameter_fullres": 80.0}, handle)
        image_bytes = io.BytesIO()
        Image.new("RGB", (4, 3), color=(20, 40, 60)).save(image_bytes, format="PNG")
        with gzip.open(self.spatial / "tissue_hires_image.png.gz", "wb") as handle:
            handle.write(image_bytes.getvalue())

    def test_detect_and_convert_compressed_spatial_files(self) -> None:
        detection = detect_visium_sample(self.sample)
        self.assertTrue(detection["valid"], detection["errors"])
        self.assertEqual(detection["positions"].suffix, ".gz")
        self.assertEqual(detection["scalefactors"].suffix, ".gz")
        self.assertEqual(detection["hires_image"].suffix, ".gz")

        output, report = convert_visium_to_h5ad(self.sample, self.root / "converted", "sample_01")
        self.assertTrue(output.is_file())
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(report["gene_names_made_unique"])

        adata = ad.read_h5ad(output)
        self.assertEqual(adata.shape, (3, 2))
        np.testing.assert_array_equal(
            adata.X.toarray(),
            np.asarray([[1, 0], [0, 2], [3, 4]], dtype=np.int32),
        )
        self.assertTrue(adata.var_names.is_unique)
        self.assertEqual(list(adata.obs_names), ["BC1", "BC2", "BC3"])
        np.testing.assert_allclose(
            adata.obsm["spatial"],
            np.asarray([[110, 100], [210, 200], [310, 300]], dtype=float),
        )
        self.assertIn("array_row", adata.obs)
        self.assertIn("spatialtx_import", adata.uns)
        self.assertIn("hires", adata.uns["spatial"]["sample_01"]["images"])
        self.assertEqual(adata.uns["spatialtx_import"]["source_folder_name"], "sample_01")
        self.assertTrue(validate_h5ad(output)["valid"])

    def test_detection_fails_when_required_scalefactors_are_missing(self) -> None:
        (self.spatial / "scalefactors_json.json.gz").unlink()
        report = detect_visium_sample(self.sample)
        self.assertFalse(report["valid"])
        self.assertTrue(any("scalefactors" in error for error in report["errors"]))

    def test_detects_geo_prefixed_visium_filenames(self) -> None:
        prefix = "GSM9532669_YUBOISE_"
        for path in list(self.sample.iterdir()):
            if path.is_file():
                path.rename(path.with_name(prefix + path.name))
        for path in list(self.spatial.iterdir()):
            path.rename(self.sample / (prefix + path.name))
        self.spatial.rmdir()

        report = detect_visium_sample(self.sample)
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["suggested_sample_name"], "GSM9532669_YUBOISE")
        self.assertTrue(report["matrix_h5"].name.startswith(prefix))
        self.assertTrue(report["positions"].name.startswith(prefix))

    def test_geo_prefixed_gzip_equivalents_warn_and_prefer_uncompressed(self) -> None:
        prefix = "GSM9532669_YUBOISE_"
        for path in list(self.sample.iterdir()):
            if path.is_file():
                path.rename(path.with_name(prefix + path.name))
        for path in list(self.spatial.iterdir()):
            path.rename(self.sample / (prefix + path.name))
        self.spatial.rmdir()

        for base_name in (
            "tissue_positions.csv",
            "scalefactors_json.json",
            "tissue_hires_image.png",
        ):
            compressed = self.sample / f"{prefix}{base_name}.gz"
            uncompressed = self.sample / f"{prefix}{base_name}"
            mode = "rt" if base_name.endswith((".csv", ".json")) else "rb"
            kwargs = {"encoding": "utf-8"} if mode == "rt" else {}
            with gzip.open(compressed, mode, **kwargs) as source:
                payload = source.read()
            if isinstance(payload, str):
                uncompressed.write_text(payload, encoding="utf-8")
            else:
                uncompressed.write_bytes(payload)

        lowres_payload = io.BytesIO()
        Image.new("RGB", (2, 2), color=(70, 80, 90)).save(lowres_payload, format="PNG")
        lowres_png = self.sample / f"{prefix}tissue_lowres_image.png"
        lowres_png.write_bytes(lowres_payload.getvalue())
        with gzip.open(self.sample / f"{prefix}tissue_lowres_image.png.gz", "wb") as handle:
            handle.write(lowres_payload.getvalue())

        report = detect_visium_sample(self.sample)
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["positions"].name, f"{prefix}tissue_positions.csv")
        self.assertEqual(report["scalefactors"].name, f"{prefix}scalefactors_json.json")
        self.assertEqual(report["hires_image"].name, f"{prefix}tissue_hires_image.png")
        self.assertEqual(report["lowres_image"].name, f"{prefix}tissue_lowres_image.png")
        duplicate_warnings = [
            warning for warning in report["warnings"]
            if "Duplicate-equivalent" in warning
        ]
        self.assertGreaterEqual(len(duplicate_warnings), 4)

    def test_positions_conflict_does_not_skip_scalefactors_and_images(self) -> None:
        with gzip.open(self.spatial / "tissue_positions.csv.gz", "rt", encoding="utf-8") as source:
            payload = source.read()
        (self.spatial / "tissue_positions.csv.gz").unlink()
        (self.sample / "GSM_A_tissue_positions.csv").write_text(payload, encoding="utf-8")
        (self.sample / "GSM_B_tissue_positions.csv").write_text(payload, encoding="utf-8")

        report = detect_visium_sample(self.sample)
        self.assertFalse(report["valid"])
        self.assertTrue(
            any("Multiple conflicting tissue positions" in error for error in report["errors"]),
            report["errors"],
        )
        self.assertIsNotNone(report["scalefactors"])
        self.assertIsNotNone(report["hires_image"])

    def test_validation_rejects_h5ad_without_spatial_coordinates(self) -> None:
        path = self.root / "no_spatial.h5ad"
        ad.AnnData(
            X=np.ones((2, 2)),
            obs=pd.DataFrame(index=["A", "B"]),
            var=pd.DataFrame(index=["G1", "G2"]),
        ).write_h5ad(path)
        report = validate_h5ad(path)
        self.assertFalse(report["valid"])
        self.assertFalse(report["has_spatial"])


if __name__ == "__main__":
    unittest.main()
