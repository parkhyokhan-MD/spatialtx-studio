from __future__ import annotations

import tempfile
import unittest
import json
from pathlib import Path
from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from spatialtx_desktop.workflow import DEFAULT_C_GENES, DEFAULT_S_GENES, run_batch, scan_h5ad, score_adata, score_h5ad
from spatialtx_desktop.graph.builder import GraphBuildConfig, build_spatial_graph
from spatialtx_desktop.graph.context import ContextFieldConfig, add_context_field
from spatialtx_desktop.graph.qc import graph_qc
from spatialtx_desktop.graph.runner import SpatialGraphAnalysisConfig, run_spatial_graph_neighborhood_batch
from spatialtx_desktop.graph.metadata import FDR_SCOPE, PERMUTATION_LIMITATION, SMOOTHING_LIMITATION
from spatialtx_desktop.graph.plotting import plot_graph_qc
from spatialtx_studio.io import load_config, packaged_default_config
from spatialtx_studio.frame26 import run_frame26


def _regular_adata(n_rows: int = 4, n_cols: int = 4) -> ad.AnnData:
    genes = DEFAULT_C_GENES + DEFAULT_S_GENES + ["CA9", "VWF"]
    n = n_rows * n_cols
    rng = np.random.default_rng(42)
    matrix = rng.poisson(2.0, size=(n, len(genes))).astype(float)
    obs = pd.DataFrame(
        {
            "array_row": np.repeat(np.arange(n_rows), n_cols),
            "array_col": np.tile(np.arange(n_cols), n_rows) * 2 + np.repeat(np.arange(n_rows) % 2, n_cols),
            "in_tissue": 1,
        },
        index=[f"spot_{i}" for i in range(n)],
    )
    result = ad.AnnData(sparse.csr_matrix(matrix), obs=obs, var=pd.DataFrame(index=genes))
    result.obsm["spatial"] = obs[["array_col", "array_row"]].to_numpy(dtype=float)
    return result


class CanonicalEngineTests(unittest.TestCase):
    def test_score_h5ad_and_score_adata_are_identical(self):
        adata = _regular_adata()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.h5ad"
            adata.write_h5ad(path)
            memory_metrics, memory_fields = score_adata(
                adata,
                DEFAULT_C_GENES,
                DEFAULT_S_GENES,
                source_path=path,
                sample_name="synthetic",
            )
            file_metrics, file_fields = score_h5ad(path, DEFAULT_C_GENES, DEFAULT_S_GENES)
        for key in ("C", "S", "R", "G", "interface", "diffuse"):
            np.testing.assert_allclose(memory_fields[key], file_fields[key], equal_nan=True)
        for key in (
            "regime_label",
            "public_transition_pattern",
            "transition_burden_score",
            "adj_same_fraction",
            "adj_zero_fraction",
            "adj_opposite_fraction",
        ):
            self.assertEqual(memory_metrics[key], file_metrics[key])

    def test_cli_frame26_wrapper_uses_main_mapper_fields(self):
        adata = _regular_adata()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.h5ad"
            output = Path(tmp) / "output"
            adata.write_h5ad(path)
            expected_metrics, expected_fields = score_h5ad(path, DEFAULT_C_GENES, DEFAULT_S_GENES)
            summary, annotated = run_frame26(path, "synthetic", output)
            for key, column in (
                ("C", "frame26_C"),
                ("S", "frame26_B"),
                ("R", "frame26_R"),
                ("G", "frame26_G"),
                ("interface", "frame26_interface_like"),
                ("diffuse", "frame26_diffuse_transition"),
            ):
                np.testing.assert_allclose(np.asarray(annotated.obs[column]), expected_fields[key])
            self.assertEqual(summary.iloc[0]["regime_label"], expected_metrics["regime_label"])
            self.assertEqual(summary.iloc[0]["common_engine"], "spatialtx_desktop.workflow.score_adata")


class GraphStabilizationTests(unittest.TestCase):
    def test_regular_irregular_and_duplicate_datasets_across_all_graphs(self):
        regular = _regular_adata(8, 8)
        irregular = _regular_adata(8, 8)
        irregular.obsm["spatial"] = np.random.default_rng(9).uniform(0, 100, size=(irregular.n_obs, 2))
        irregular.obs = irregular.obs.drop(columns=["array_row", "array_col"])
        duplicate = _regular_adata(8, 8)
        duplicate_coords = np.asarray(duplicate.obsm["spatial"], dtype=float)
        duplicate_coords[1] = duplicate_coords[0]
        duplicate.obsm["spatial"] = duplicate_coords
        for dataset_name, adata in (
            ("regular_visium_like", regular),
            ("irregular_coordinates", irregular),
            ("duplicate_coordinate", duplicate),
        ):
            for method in ("radius", "lattice", "knn"):
                with self.subTest(dataset=dataset_name, method=method):
                    result = build_spatial_graph(
                        adata,
                        GraphBuildConfig(method=method, radius=25.0, k=6),
                    )
                    self.assertTrue(sparse.issparse(result.connectivities))
                    self.assertGreater(result.qc["n_edges"], 0)
                    self.assertEqual(result.qc["n_nodes"], adata.n_obs)

    def test_large_local_graph_is_not_warned_for_low_density_alone(self):
        n = 2000
        row = np.arange(n - 1)
        col = row + 1
        values = np.ones(2 * (n - 1), dtype=float)
        conn = sparse.csr_matrix(
            (values, (np.r_[row, col], np.r_[col, row])), shape=(n, n)
        )
        dist = conn.copy()
        coords = np.column_stack([np.arange(n, dtype=float), np.zeros(n)])
        qc, warnings, _ = graph_qc(
            conn,
            dist,
            coords,
            method="radius",
            isolated_fraction_warning=0.10,
            largest_component_ratio_warning=0.80,
            median_degree_warning=1.0,
        )
        self.assertLess(qc["graph_density"], 0.002)
        self.assertFalse(any("almost empty" in warning for warning in warnings))
        self.assertFalse(any("low density" in warning for warning in warnings))

    def test_inverse_distance_excludes_duplicate_coordinate_edge(self):
        adata = _regular_adata(2, 2)
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        coords[1] = coords[0]
        adata.obsm["spatial"] = coords
        result = build_spatial_graph(
            adata,
            GraphBuildConfig(method="radius", radius=3.0, weighting="inverse_distance"),
        )
        self.assertGreaterEqual(result.qc["duplicate_coordinate_count"], 1)
        self.assertGreaterEqual(result.qc["zero_distance_edges_excluded"], 1)
        self.assertTrue(np.isfinite(result.connectivities.data).all())
        self.assertFalse(np.any(result.distances.data == 0))

    def test_all_zero_distance_inverse_edges_are_invalid(self):
        adata = _regular_adata(1, 3)
        adata.obsm["spatial"] = np.zeros((adata.n_obs, 2), dtype=float)
        result = build_spatial_graph(
            adata,
            GraphBuildConfig(method="radius", radius=1.0, weighting="inverse_distance"),
        )
        self.assertEqual(result.qc["n_edges"], 0)
        self.assertGreater(result.qc["zero_distance_edges_excluded"], 0)
        self.assertTrue(any("no usable spatial edges" in warning for warning in result.warnings))


class StatisticalGuardrailTests(unittest.TestCase):
    def test_centered_context_reports_positive_fraction_not_detection(self):
        adata = _regular_adata(2, 3)
        dense = adata.X.toarray().astype(float)
        dense -= dense.mean(axis=0, keepdims=True)
        adata.X = sparse.csr_matrix(dense)
        table, metadata = add_context_field(
            adata,
            ContextFieldConfig(field="H", genes=DEFAULT_C_GENES[:2], min_coverage=1.0),
        )
        self.assertEqual(metadata["expression_scale_guess"], "centered_or_scaled")
        self.assertEqual(
            metadata["detection_metric_interpretation"],
            "detection_fraction_unavailable_positive_value_fraction_only",
        )
        self.assertTrue(table["detection_fraction"].isna().all())
        self.assertTrue(table["positive_value_fraction"].notna().all())

    def test_counts_layer_is_preferred_for_detection(self):
        adata = _regular_adata(2, 3)
        adata.layers["counts"] = adata.X.copy()
        adata.X = sparse.csr_matrix(adata.X.toarray() - 2.0)
        table, metadata = add_context_field(
            adata,
            ContextFieldConfig(field="H", genes=DEFAULT_C_GENES[:2], min_coverage=1.0),
        )
        self.assertEqual(metadata["detection_source"], "adata.layers['counts']")
        self.assertEqual(metadata["expression_scale_guess"], "raw_counts")
        self.assertTrue(table["detection_fraction"].notna().all())

    def test_smoothing_warning_and_leave_one_out_are_recorded(self):
        adata = _regular_adata(2, 3)
        graph = build_spatial_graph(adata, GraphBuildConfig(method="knn", k=2))
        table, metadata = add_context_field(
            adata,
            ContextFieldConfig(
                field="H",
                genes=DEFAULT_C_GENES[:3],
                min_coverage=1.0,
                smoothing="graph_mean",
            ),
            graph.connectivities,
            active_graph="knn",
        )
        self.assertEqual(metadata["smoothing_warning"], SMOOTHING_LIMITATION)
        self.assertIn(SMOOTHING_LIMITATION, metadata["warnings"])
        self.assertTrue({"field_correlation_without_gene", "mean_absolute_field_change", "variance_change", "rank"}.issubset(table.columns))

    def test_runner_records_derived_source_fdr_scope_and_single_load(self):
        adata = _regular_adata(3, 3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.h5ad"
            adata.write_h5ad(path)
            import spatialtx_desktop.graph.runner as graph_runner

            original_read = graph_runner._read_h5ad
            with mock.patch.object(graph_runner, "_read_h5ad", wraps=original_read) as read_mock:
                run_dir, manifest = run_spatial_graph_neighborhood_batch(
                    [path],
                    Path(tmp) / "results",
                    DEFAULT_C_GENES,
                    DEFAULT_S_GENES,
                    SpatialGraphAnalysisConfig(
                        graph=GraphBuildConfig(method="knn", k=3),
                        enable_h=False,
                        enable_v=False,
                        permutations=5,
                    ),
                )
            self.assertEqual(manifest.iloc[0]["status"], "ok")
            self.assertEqual(read_mock.call_count, 1)
            categorical = pd.read_csv(run_dir / "neighborhood" / "sample_categorical_enrichment.csv")
            self.assertTrue(categorical["label_source"].eq("spatialtx_derived").all())
            self.assertTrue(categorical["analysis_interpretation"].eq("descriptive_spatial_organization").all())
            self.assertTrue(categorical["fdr_scope"].eq(FDR_SCOPE).all())
            parameters = json.loads((run_dir / "neighborhood" / "sample_permutation_parameters.json").read_text(encoding="utf-8"))
            self.assertEqual(parameters["fdr_scope"], FDR_SCOPE)
            self.assertEqual(parameters["permutation_limitation"], PERMUTATION_LIMITATION)

    def test_sparse_input_is_not_replaced_or_densified(self):
        adata = _regular_adata(20, 20)
        original_x = adata.X
        score_adata(adata, DEFAULT_C_GENES, DEFAULT_S_GENES)
        self.assertIs(adata.X, original_x)
        self.assertTrue(sparse.issparse(adata.X))


class PackagingAndPlottingTests(unittest.TestCase):
    def test_packaged_default_config_is_loadable(self):
        resource = packaged_default_config()
        self.assertTrue(resource.is_file())
        config = load_config()
        self.assertEqual(config["analysis"], "frame26")
        self.assertEqual(config["smoothing"]["mode"], "none")

    def test_graph_plot_uses_deterministic_edge_cap_and_metadata(self):
        n = 60
        upper_i, upper_j = np.triu_indices(n, k=1)
        conn = sparse.csr_matrix(
            (
                np.ones(2 * len(upper_i)),
                (np.r_[upper_i, upper_j], np.r_[upper_j, upper_i]),
            ),
            shape=(n, n),
        )
        coords = np.column_stack([np.arange(n, dtype=float), np.sin(np.arange(n))])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "graph.png"
            plot_graph_qc(
                coords,
                conn,
                np.zeros(n, dtype=bool),
                path,
                "test",
                max_plot_edges=50,
                plot_edge_seed=42,
            )
            metadata = json.loads(path.with_suffix(".png.metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["displayed_plot_edges"], 50)
        self.assertEqual(metadata["total_graph_edges"], len(upper_i))
        self.assertTrue(metadata["deterministic_downsampling_applied"])
        self.assertTrue(metadata["statistics_use_full_graph"])


class WindowsPathAndExportTests(unittest.TestCase):
    def test_korean_and_space_path_scan_main_graph_and_exports(self):
        adata = _regular_adata(4, 4)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "한글 경로 with space"
            input_dir = root / "입력 자료"
            output_dir = root / "분석 결과"
            input_dir.mkdir(parents=True)
            path = input_dir / "공간 샘플.h5ad"
            adata.write_h5ad(path)
            self.assertEqual(scan_h5ad(input_dir), [path.resolve()])

            main_run, main_summary = run_batch(
                [path], output_dir / "Main Mapper",
                DEFAULT_C_GENES, DEFAULT_S_GENES,
            )
            self.assertEqual(main_summary.iloc[0]["status"], "ok")
            self.assertTrue(any(main_run.rglob("*.csv")))
            self.assertTrue(any(main_run.rglob("*.json")))
            self.assertTrue(any(main_run.rglob("*.png")))

            graph_run, graph_manifest = run_spatial_graph_neighborhood_batch(
                [path],
                output_dir / "Spatial Graph",
                DEFAULT_C_GENES,
                DEFAULT_S_GENES,
                SpatialGraphAnalysisConfig(
                    graph=GraphBuildConfig(method="knn", k=3),
                    enable_h=False,
                    enable_v=False,
                    permutations=5,
                    write_annotated_h5ad=True,
                    max_plot_edges=100,
                ),
            )
            self.assertEqual(graph_manifest.iloc[0]["status"], "ok")
            for suffix in ("*.csv", "*.json", "*.png", "*.h5ad"):
                self.assertTrue(any(graph_run.rglob(suffix)), suffix)


if __name__ == "__main__":
    unittest.main()
