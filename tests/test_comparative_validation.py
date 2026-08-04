from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from spatialtx_desktop.comparative.models import ComparativeConfig, SampleRecord
from spatialtx_desktop.comparative.runner import run_comparative_analysis
from spatialtx_desktop.comparative.validation import load_comparative_manifest, preflight_sample
from tests.comparative_fixtures import C_GENES, S_GENES, write_comparative_h5ad


class ComparativeValidationTests(unittest.TestCase):
    def _config(self, mode: str = "unpaired") -> ComparativeConfig:
        return ComparativeConfig(
            mode=mode, reference="A", target="B", c_genes=C_GENES, s_genes=S_GENES,
            enable_h_expr=False, enable_v_expr=False, bootstrap_iterations=100,
        )

    def test_manifest_required_optional_relative_paths_and_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_comparative_h5ad(root / "a.h5ad")
            write_comparative_h5ad(root / "b.h5ad", pattern="diffuse")
            manifest = root / "manifest.csv"
            pd.DataFrame([
                {"sample_id": "a", "file_path": "a.h5ad", "group": "A"},
                {"sample_id": "b", "file_path": "b.h5ad", "group": "B"},
            ]).to_csv(manifest, index=False)
            records = load_comparative_manifest(manifest)
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0].pair_id, "")
            self.assertTrue(records[0].file_path.is_absolute())
            pd.DataFrame([
                {"sample_id": "same", "file_path": "a.h5ad", "group": "A"},
                {"sample_id": "same", "file_path": "b.h5ad", "group": "B"},
            ]).to_csv(manifest, index=False)
            with self.assertRaisesRegex(ValueError, "Duplicate sample_id"):
                load_comparative_manifest(manifest)

    def test_no_spatial_coordinates_fail_clearly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing = write_comparative_h5ad(root / "missing.h5ad", spatial=False)
            row = preflight_sample(SampleRecord("missing", missing, "A"), self._config())
            self.assertEqual(row["validation_status"], "failed")
            self.assertIn("spatial", row["validation_error"].lower())
            valid = write_comparative_h5ad(root / "valid.h5ad", pattern="diffuse")
            with self.assertRaisesRegex(ValueError, "Fewer than two valid spatial H5AD"):
                run_comparative_analysis(
                    [SampleRecord("missing", missing, "A"), SampleRecord("valid", valid, "B")],
                    self._config(),
                    root / "out",
                )

    def test_mixed_failure_batch_retains_valid_samples_and_logs_invalid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = write_comparative_h5ad(root / "a.h5ad", pattern="localized")
            b = write_comparative_h5ad(root / "b.h5ad", pattern="diffuse")
            bad = write_comparative_h5ad(root / "bad.h5ad", spatial=False)
            result = run_comparative_analysis(
                [SampleRecord("a", a, "A"), SampleRecord("b", b, "B"), SampleRecord("bad", bad, "A")],
                self._config(),
                root / "out",
            )
            status = result.run_manifest.set_index("sample_id")["status"]
            self.assertEqual(status["a"], "ok")
            self.assertEqual(status["b"], "ok")
            self.assertEqual(status["bad"], "failed_validation")
            self.assertTrue(result.warnings["message"].str.contains("spatial", case=False).any())


if __name__ == "__main__":
    unittest.main()
