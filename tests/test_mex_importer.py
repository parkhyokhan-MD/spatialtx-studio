from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from spatialtx_desktop.importers.mex_to_h5ad import convert_mex_to_h5ad, detect_mex_sample
from spatialtx_desktop.importers.validate_h5ad import validate_h5ad


class MexImporterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.sample = self.root / "mex_sample"
        self.sample.mkdir()
        feature_by_barcode = sparse.coo_matrix(np.asarray([
            [1, 0, 3],
            [0, 2, 4],
        ], dtype=np.int32))
        with gzip.open(self.sample / "matrix.mtx.gz", "wb") as handle:
            mmwrite(handle, feature_by_barcode)
        with gzip.open(self.sample / "features.tsv.gz", "wt", encoding="utf-8") as handle:
            handle.write("ENSG1\tGENE\tGene Expression\nENSG2\tGENE\tGene Expression\n")
        with gzip.open(self.sample / "barcodes.tsv.gz", "wt", encoding="utf-8") as handle:
            handle.write("BC1\nBC2\nBC3\n")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_detect_and_convert_gzipped_mex_to_h5ad(self) -> None:
        detection = detect_mex_sample(self.sample)
        self.assertTrue(detection["valid"], detection["errors"])
        self.assertTrue(any("expression-only scoring" in warning for warning in detection["warnings"]))

        output, report = convert_mex_to_h5ad(self.sample, self.root / "converted", "mex_sample")
        self.assertTrue(output.is_file())
        self.assertTrue(report["valid"], report["errors"])
        self.assertFalse(report["has_spatial"])
        self.assertTrue(report["gene_names_made_unique"])

        adata = ad.read_h5ad(output)
        self.assertEqual(adata.shape, (3, 2))
        np.testing.assert_array_equal(
            adata.X.toarray(),
            np.asarray([[1, 0], [0, 2], [3, 4]], dtype=np.int32),
        )
        self.assertTrue(adata.var_names.is_unique)
        self.assertEqual(adata.uns["spatialtx_import"]["source_format"], "10x MEX/MTX")
        self.assertTrue(validate_h5ad(output, require_spatial=False)["valid"])
        self.assertFalse(validate_h5ad(output, require_spatial=True)["valid"])

    def test_detection_rejects_incomplete_mex_folder(self) -> None:
        (self.sample / "barcodes.tsv.gz").unlink()
        report = detect_mex_sample(self.sample)
        self.assertFalse(report["valid"])
        self.assertTrue(any("barcodes.tsv" in error for error in report["errors"]))


if __name__ == "__main__":
    unittest.main()
