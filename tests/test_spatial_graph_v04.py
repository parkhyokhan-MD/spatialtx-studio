from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from spatialtx_desktop.graph.builder import GraphBuildConfig, build_spatial_graph, store_graph
from spatialtx_desktop.graph.context import ContextFieldConfig, add_context_field
from spatialtx_desktop.graph.continuous import continuous_edge_statistics
from spatialtx_desktop.graph.enrichment import binary_mask_association, categorical_neighborhood_enrichment
from spatialtx_desktop.graph.runner import SpatialGraphAnalysisConfig, run_spatial_graph_neighborhood_batch


class FakeAdata:
    def __init__(self, X, coords=None, obs=None, var_names=None):
        self.X = np.asarray(X, dtype=float)
        self.n_obs = self.X.shape[0]
        self.n_vars = self.X.shape[1]
        self.obs = pd.DataFrame(index=[f"spot{i}" for i in range(self.n_obs)]) if obs is None else obs.copy()
        self.var_names = pd.Index(var_names or [f"gene{i}" for i in range(self.n_vars)])
        self.obsm = {}
        if coords is not None:
            self.obsm["spatial"] = np.asarray(coords, dtype=float)
        self.obsp = {}
        self.uns = {}


class SpatialGraphV04Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.coords = np.asarray([
            [0.0, 0.0], [1.0, 0.0], [2.0, 0.0],
            [0.0, 1.0], [1.0, 1.0], [2.0, 1.0],
        ])
        self.adata = FakeAdata(
            np.asarray([
                [5, 1, 0, 0],
                [4, 1, 1, 0],
                [0, 1, 4, 2],
                [5, 0, 0, 1],
                [1, 4, 4, 0],
                [0, 3, 5, 1],
            ], dtype=float),
            self.coords,
            var_names=["CA9", "VEGFA", "PECAM1", "VWF"],
        )

    def test_radius_graph_generation_and_sparse_storage(self) -> None:
        result = build_spatial_graph(self.adata, GraphBuildConfig(method="radius", radius=1.01))
        self.assertTrue(sparse.isspmatrix_csr(result.connectivities))
        self.assertTrue(sparse.isspmatrix_csr(result.distances))
        self.assertEqual(result.qc["n_nodes"], 6)
        self.assertGreater(result.qc["n_edges"], 0)
        store_graph(self.adata, result)
        self.assertIn("spatialtx_connectivities_radius", self.adata.obsp)
        self.assertIn("spatialtx_graph", self.adata.uns)
        self.assertEqual(self.adata.uns["spatialtx_graph"]["active_graph"], "radius")

    def test_knn_union_has_at_least_mutual_edges(self) -> None:
        union = build_spatial_graph(self.adata, GraphBuildConfig(method="knn", k=2, symmetrization="union"))
        mutual = build_spatial_graph(self.adata, GraphBuildConfig(method="knn", k=2, symmetrization="mutual"))
        self.assertGreaterEqual(union.qc["n_edges"], mutual.qc["n_edges"])
        self.assertTrue(sparse.issparse(union.connectivities))

    def test_lattice_graph_uses_array_coordinates_when_available(self) -> None:
        obs = pd.DataFrame({
            "array_row": [0, 0, 1, 1],
            "array_col": [0, 2, 1, 3],
        })
        adata = FakeAdata(np.ones((4, 2)), coords=np.asarray([[0, 0], [2, 0], [1, 1], [3, 1]], dtype=float), obs=obs)
        result = build_spatial_graph(adata, GraphBuildConfig(method="lattice"))
        self.assertTrue(result.metadata["lattice_detection_succeeded"])
        self.assertGreater(result.qc["n_edges"], 0)

    def test_lattice_fallback_records_effective_radius_graph(self) -> None:
        result = build_spatial_graph(self.adata, GraphBuildConfig(method="lattice", radius=1.01))
        self.assertEqual(result.method, "radius")
        self.assertEqual(result.metadata["requested_method"], "lattice")
        self.assertEqual(result.metadata["effective_method"], "radius")
        self.assertTrue(result.metadata["fallback_used"])
        self.assertAlmostEqual(result.metadata["radius"], 1.01)
        store_graph(self.adata, result)
        self.assertIn("spatialtx_connectivities_radius", self.adata.obsp)

    def test_duplicate_and_missing_coordinates_warn(self) -> None:
        duplicate = FakeAdata(np.ones((3, 2)), coords=np.asarray([[0, 0], [0, 0], [1, 0]], dtype=float))
        result = build_spatial_graph(duplicate, GraphBuildConfig(method="radius", radius=1.5))
        self.assertGreater(result.qc["duplicate_coordinate_count"], 0)
        missing = FakeAdata(np.ones((3, 2)), coords=None)
        missing_result = build_spatial_graph(missing, GraphBuildConfig(method="radius"))
        self.assertFalse(missing_result.qc["coordinate_valid"])
        self.assertIn("missing spatial coordinates", "; ".join(missing_result.warnings))

    def test_categorical_enrichment_reproducible_and_empirical_p_corrected(self) -> None:
        result = build_spatial_graph(self.adata, GraphBuildConfig(method="radius", radius=1.01))
        labels = np.asarray(["A", "A", "B", "A", "A", "B"], dtype=object)
        table1 = categorical_neighborhood_enrichment(result.connectivities, labels, permutations=49, seed=7)
        table2 = categorical_neighborhood_enrichment(result.connectivities, labels, permutations=49, seed=7)
        pd.testing.assert_frame_equal(table1, table2)
        self.assertTrue(((table1["empirical_p"] >= 1 / 50) & (table1["empirical_p"] <= 1)).all())
        aa = table1[(table1["label_a"] == "A") & (table1["label_b"] == "A")].iloc[0]
        self.assertGreaterEqual(aa["observed_edge_count"], aa["expected_edge_count"])

    def test_binary_mask_association_distinguishes_overlap_and_neighbors(self) -> None:
        result = build_spatial_graph(self.adata, GraphBuildConfig(method="radius", radius=1.01))
        masks = {
            "left": np.asarray([True, True, False, True, False, False]),
            "right": np.asarray([False, False, True, False, True, True]),
        }
        table = binary_mask_association(result.connectivities, masks, pairs=[("left", "right")], permutations=19, seed=3)
        self.assertEqual(set(table["mode"]), {"same_spot_overlap", "neighboring_spot_association"})
        overlap = table[table["mode"].eq("same_spot_overlap")].iloc[0]
        self.assertEqual(overlap["observed_count"], 0)

    def test_h_and_v_context_fields_write_metadata(self) -> None:
        result = build_spatial_graph(self.adata, GraphBuildConfig(method="radius", radius=1.01))
        h_table, h_meta = add_context_field(
            self.adata,
            ContextFieldConfig(field="H", genes=["CA9", "VEGFA", "MISSING"], min_coverage=0.2, smoothing="graph_mean"),
            result.connectivities,
        )
        v_table, v_meta = add_context_field(
            self.adata,
            ContextFieldConfig(field="V", genes=["PECAM1", "VWF"], min_coverage=0.2),
            result.connectivities,
        )
        self.assertIn("spatialtx_H_expr", self.adata.obs)
        self.assertIn("spatialtx_H_expr_smoothed", self.adata.obs)
        self.assertIn("spatialtx_H_high", self.adata.obs)
        self.assertIn("spatialtx_V_expr", self.adata.obs)
        self.assertAlmostEqual(h_meta["coverage_fraction"], 2 / 3)
        self.assertAlmostEqual(v_meta["coverage_fraction"], 1.0)
        self.assertEqual(len(h_table), 3)
        self.assertEqual(len(v_table), 2)
        self.assertIn("spatialtx_context_fields", self.adata.uns)

    def test_continuous_stats_handle_missing_context_fields(self) -> None:
        result = build_spatial_graph(self.adata, GraphBuildConfig(method="radius", radius=1.01))
        fields = {
            "C": np.linspace(0, 1, self.adata.n_obs),
            "S": np.linspace(1, 0, self.adata.n_obs),
            "R": np.linspace(-1, 1, self.adata.n_obs),
        }
        table = continuous_edge_statistics(result.connectivities, fields, permutations=19, seed=5)
        self.assertIn("I_CS", set(table["statistic"]))
        self.assertEqual(table[table["statistic"].eq("I_HV")].iloc[0]["status"], "missing_field")
        ok = table[table["status"].eq("ok")]
        self.assertTrue((ok["empirical_p"] >= 1 / 20).all())

    def test_empty_graph_is_not_interpreted_as_depletion(self) -> None:
        empty = sparse.csr_matrix((self.adata.n_obs, self.adata.n_obs), dtype=float)
        labels = np.asarray(["A", "A", "B", "B", "A", "B"], dtype=object)
        categorical = categorical_neighborhood_enrichment(empty, labels, permutations=9, seed=2)
        self.assertEqual(set(categorical["interpretation"]), {"not_interpretable"})
        fields = {
            "C": np.linspace(0, 1, self.adata.n_obs),
            "S": np.linspace(1, 0, self.adata.n_obs),
            "R": np.linspace(-1, 1, self.adata.n_obs),
        }
        continuous = continuous_edge_statistics(empty, fields, permutations=9, seed=2)
        available = continuous[~continuous["status"].eq("missing_field")]
        self.assertEqual(set(available["status"]), {"insufficient_edges"})
        with self.assertRaisesRegex(ValueError, "permutations"):
            categorical_neighborhood_enrichment(empty, labels, permutations=0, seed=2)

    def test_graph_module_does_not_mutate_existing_cs_columns(self) -> None:
        self.adata.obs["C_score"] = np.arange(self.adata.n_obs, dtype=float)
        self.adata.obs["S_score"] = np.arange(self.adata.n_obs, dtype=float) + 10
        self.adata.obs["R_score"] = self.adata.obs["C_score"] - self.adata.obs["S_score"]
        before = self.adata.obs[["C_score", "S_score", "R_score"]].copy()
        result = build_spatial_graph(self.adata, GraphBuildConfig(method="radius", radius=1.01))
        add_context_field(self.adata, ContextFieldConfig(field="H", genes=["CA9", "VEGFA"], min_coverage=0.2), result.connectivities)
        pd.testing.assert_frame_equal(before, self.adata.obs[["C_score", "S_score", "R_score"]])


class OptionalAnnDataRegressionTests(unittest.TestCase):
    def test_existing_main_mapper_scores_unchanged_when_graph_not_enabled(self) -> None:
        try:
            import anndata as ad
        except ModuleNotFoundError:
            self.skipTest("anndata is not installed in this runtime")
        from spatialtx_desktop.workflow import score_h5ad

        matrix = np.asarray([
            [5, 4, 0, 0],
            [4, 3, 0, 1],
            [0, 1, 5, 4],
            [1, 0, 4, 5],
            [3, 3, 3, 3],
            [2, 2, 2, 2],
        ], dtype=float)
        adata = ad.AnnData(matrix, var=pd.DataFrame(index=["C1", "C2", "S1", "S2"]))
        adata.obsm["spatial"] = np.asarray([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "regression.h5ad"
            adata.write_h5ad(path)
            before, before_fields = score_h5ad(path, ["C1", "C2"], ["S1", "S2"])
            after, after_fields = score_h5ad(path, ["C1", "C2"], ["S1", "S2"])
        np.testing.assert_allclose(before_fields["C"], after_fields["C"])
        np.testing.assert_allclose(before_fields["S"], after_fields["S"])
        np.testing.assert_allclose(before_fields["R"], after_fields["R"])
        np.testing.assert_allclose(before_fields["G"], after_fields["G"])
        np.testing.assert_array_equal(before_fields["interface"], after_fields["interface"])
        np.testing.assert_array_equal(before_fields["diffuse"], after_fields["diffuse"])
        self.assertEqual(before["regime_label"], after["regime_label"])
        self.assertEqual(before["public_transition_pattern"], after["public_transition_pattern"])

    def test_graph_runner_rejects_zero_edge_graph(self) -> None:
        try:
            import anndata as ad
        except ModuleNotFoundError:
            self.skipTest("anndata is not installed in this runtime")
        matrix = np.arange(8 * 4, dtype=float).reshape(8, 4) % 7
        adata = ad.AnnData(matrix, var=pd.DataFrame(index=["C1", "C2", "S1", "S2"]))
        adata.obsm["spatial"] = np.asarray([[x, y] for y in range(2) for x in range(4)], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "no_edges.h5ad"
            adata.write_h5ad(source)
            config = SpatialGraphAnalysisConfig(
                graph=GraphBuildConfig(method="radius", radius=0.01),
                enable_h=False,
                enable_v=False,
                permutations=9,
            )
            _, manifest = run_spatial_graph_neighborhood_batch(
                [source], root / "out", ["C1", "C2"], ["S1", "S2"], config
            )
        self.assertTrue(str(manifest.iloc[0]["status"]).startswith("error:"))
        self.assertIn("no usable spatial edges", str(manifest.iloc[0]["status"]))

    def test_runner_uses_smoothed_context_for_continuous_statistics(self) -> None:
        try:
            import anndata as ad
        except ModuleNotFoundError:
            self.skipTest("anndata is not installed in this runtime")
        from spatialtx_desktop.workflow import score_h5ad

        genes = ["C1", "C2", "S1", "S2", "CA9", "VEGFA"]
        rng = np.random.default_rng(12)
        matrix = rng.poisson(3, size=(16, len(genes))).astype(float)
        matrix[:, 4] += np.repeat([0, 0, 8, 8], 4)
        adata = ad.AnnData(matrix, var=pd.DataFrame(index=genes))
        adata.obsm["spatial"] = np.asarray([[x, y] for y in range(4) for x in range(4)], dtype=float)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "smoothed.h5ad"
            adata.write_h5ad(source)
            config = SpatialGraphAnalysisConfig(
                graph=GraphBuildConfig(method="radius", radius=1.01),
                enable_h=True,
                enable_v=False,
                h_genes=["CA9", "VEGFA"],
                min_gene_coverage=0.5,
                context_smoothing="graph_mean",
                permutations=9,
                seed=4,
                write_annotated_h5ad=True,
            )
            run_dir, manifest = run_spatial_graph_neighborhood_batch(
                [source], root / "out", ["C1", "C2"], ["S1", "S2"], config
            )
            self.assertEqual(manifest.iloc[0]["status"], "ok")
            annotated = ad.read_h5ad(manifest.iloc[0]["annotated_h5ad"])
            _, fields = score_h5ad(source, ["C1", "C2"], ["S1", "S2"])
            expected = continuous_edge_statistics(
                annotated.obsp["spatialtx_connectivities_radius"],
                {
                    "C": fields["C"],
                    "S": fields["S"],
                    "R": fields["R"],
                    "H": np.asarray(annotated.obs["spatialtx_H_expr_smoothed"], dtype=float),
                },
                permutations=9,
                seed=4,
            )
            actual = pd.read_csv(run_dir / "neighborhood" / "smoothed_continuous_edge_statistics.csv")
        expected_value = expected.loc[expected["statistic"].eq("I_RH"), "raw_weighted_statistic"].iloc[0]
        actual_value = actual.loc[actual["statistic"].eq("I_RH"), "raw_weighted_statistic"].iloc[0]
        self.assertAlmostEqual(float(actual_value), float(expected_value))


if __name__ == "__main__":
    unittest.main()
