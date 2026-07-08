from __future__ import annotations

import tempfile
import unittest
import json
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

from spatialtx_desktop.advanced_analysis import calculate_gene_composition, calculate_interface_enrichment
from spatialtx_desktop.workflow import (
    SPATIAL_QC_MESSAGE,
    ScoringOptions,
    inspect_h5ad_memory,
    run_batch,
    score_h5ad,
)


class SpatialQCTests(unittest.TestCase):
    C_GENES = ["C1", "C2"]
    S_GENES = ["S1", "S2"]

    @staticmethod
    def _adata() -> ad.AnnData:
        matrix = np.asarray([
            [9, 7, 1, 0],
            [8, 6, 1, 1],
            [6, 5, 2, 1],
            [2, 1, 6, 7],
            [1, 1, 7, 8],
            [0, 1, 8, 9],
        ], dtype=float)
        return ad.AnnData(matrix, var={"gene": ["C1", "C2", "S1", "S2"]})

    def _write(self, folder: Path, name: str, spatial=None) -> Path:
        adata = self._adata()
        adata.var_names = adata.var["gene"].astype(str)
        if spatial is not None:
            adata.obsm["spatial"] = spatial
        path = folder / f"{name}.h5ad"
        adata.write_h5ad(path)
        return path

    def test_missing_spatial_keeps_expression_scores_and_suppresses_spatial_calls(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp), "missing")
            metrics, fields = score_h5ad(path, self.C_GENES, self.S_GENES)
        self.assertEqual(metrics["status"], "ok")
        self.assertEqual(metrics["analysis_scope"], "expression_only")
        self.assertEqual(metrics["spatial_qc_status"], "WARN")
        self.assertEqual(metrics["regime_label"], "Spatial_QC_incomplete")
        self.assertEqual(metrics["spatial_qc_message"], SPATIAL_QC_MESSAGE)
        self.assertTrue(np.isfinite(metrics["C_mean"]))
        self.assertTrue(np.isfinite(metrics["S_mean"]))
        self.assertTrue(np.isnan(metrics["interface_fraction"]))
        self.assertFalse(fields["spatial_available"])
        self.assertFalse(fields["interface"].any())

    def test_malformed_spatial_arrays_fail_spatial_qc(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            cases = {
                "empty": np.empty((6, 0)),
                "wrong_shape": np.zeros((6, 3)),
                "nonfinite": np.asarray([[0, 0], [1, 0], [2, 0], [np.nan, 1], [1, 1], [2, 1]]),
            }
            for name, coords in cases.items():
                with self.subTest(name=name):
                    path = self._write(folder, name, coords)
                    metrics, _ = score_h5ad(path, self.C_GENES, self.S_GENES)
                    self.assertEqual(metrics["spatial_qc_status"], "FAIL")
                    self.assertEqual(metrics["regime_label"], "Spatial_QC_incomplete")
                    self.assertNotIn(metrics["regime_label"], {"Type_A_candidate", "Type_B_candidate", "Type_C_candidate"})

    def test_batch_writes_expression_only_report_and_no_map(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = self._write(folder, "expression_only")
            run_dir, summary = run_batch([source], folder / "results", self.C_GENES, self.S_GENES)
            sample_dir = next(path for path in run_dir.iterdir() if path.is_dir())
            report = (sample_dir / "analysis_report.txt").read_text(encoding="utf-8")
            self.assertIn("Expression-only results", report)
            self.assertIn("Spatial results", report)
            self.assertIn(SPATIAL_QC_MESSAGE, report)
            self.assertFalse(list(sample_dir.glob("*_spatialtx_maps.png")))
            self.assertEqual(summary.iloc[0]["regime_label"], "Spatial_QC_incomplete")
            self.assertEqual(summary.iloc[0]["spatial_map_png"], "")

    def test_expression_only_advanced_composition_remains_available(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp), "composition")
            table, _ = calculate_gene_composition(path, self.C_GENES, self.S_GENES)
            self.assertEqual(len(table), 4)
            with self.assertRaisesRegex(ValueError, "Spatial interface and transition metrics are not interpretable"):
                calculate_interface_enrichment(path, self.C_GENES, self.S_GENES)

    def test_valid_spatial_coordinates_keep_spatial_workflow_available(self) -> None:
        coords = np.asarray([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp), "valid", coords)
            metrics, fields = score_h5ad(path, self.C_GENES, self.S_GENES)
        self.assertEqual(metrics["spatial_qc_status"], "PASS")
        self.assertEqual(metrics["analysis_scope"], "expression_and_spatial")
        self.assertIn(metrics["regime_label"], {"Type_A_candidate", "Type_B_candidate", "Type_C_candidate"})
        self.assertTrue(fields["spatial_available"])

    def test_optional_smoothing_normalization_and_memory_diagnostics(self) -> None:
        coords = np.asarray([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write(Path(tmp), "robust", coords)
            baseline, _ = score_h5ad(path, self.C_GENES, self.S_GENES)
            options = ScoringOptions(smoothing_mode="knn_mean", smoothing_k=2, normalization_mode="rank_quantile")
            metrics, fields = score_h5ad(path, self.C_GENES, self.S_GENES, options=options)
        self.assertEqual(baseline["normalization_mode"], "raw_mean")
        self.assertEqual(metrics["normalization_mode"], "rank_quantile")
        self.assertEqual(metrics["smoothing_mode"], "knn_mean")
        self.assertTrue(metrics["smoothing_applied"])
        self.assertEqual(metrics["matrix_shape"], "6x4")
        self.assertIn(metrics["matrix_sparse_dense_status"], {"dense", "sparse", "group", "unknown"})
        self.assertFalse(fields["robustness_table"].shape[0])
        self.assertLessEqual(float(np.nanmax(fields["C"])), 1.0)

    def test_perturbation_and_parameter_log_export(self) -> None:
        coords = np.asarray([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = self._write(folder, "perturb", coords)
            options = ScoringOptions(perturbation_check=True)
            run_dir, summary = run_batch([source], folder / "results", self.C_GENES, self.S_GENES, options=options)
            sample_dir = next(path for path in run_dir.iterdir() if path.is_dir())
            parameter_log = sample_dir / "parameter_log.json"
            robustness_csv = sample_dir / "robustness_perturbation.csv"
            parameter_log_exists = parameter_log.is_file()
            robustness_csv_exists = robustness_csv.is_file()
            payload = json.loads(parameter_log.read_text(encoding="utf-8"))
        self.assertTrue(parameter_log_exists)
        self.assertTrue(robustness_csv_exists)
        self.assertEqual(payload["normalization_mode"], "raw_mean")
        self.assertEqual(payload["smoothing"]["mode"], "none")
        self.assertEqual(payload["thresholds"], {"C_q": 0.8, "S_q": 0.8, "G_q": 0.6})
        self.assertEqual(payload["matrix"]["shape"], [6, 4])
        self.assertEqual(int(summary.iloc[0]["robustness_grid_evaluated"]), 27)
        self.assertGreaterEqual(float(summary.iloc[0]["regime_stability"]), 0.0)
        self.assertLessEqual(float(summary.iloc[0]["regime_stability"]), 1.0)

    def test_sparse_memory_preflight_does_not_require_dense_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sparse.h5ad"
            adata = self._adata()
            adata.X = sparse.csr_matrix(adata.X)
            adata.var_names = adata.var["gene"].astype(str)
            adata.obsm["spatial"] = np.asarray([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
            adata.write_h5ad(path)
            memory = inspect_h5ad_memory(path, dense_warning_gb=0.00000001)
            metrics, _ = score_h5ad(path, self.C_GENES, self.S_GENES)
        self.assertEqual(memory["matrix_storage"], "sparse")
        self.assertTrue(memory["dense_conversion_warning"])
        self.assertEqual(metrics["selected_C_S_genes_extracted"], 4)


if __name__ == "__main__":
    unittest.main()
