from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from spatialtx_desktop.graph.audit import audit_input
from spatialtx_desktop.graph.builder import GraphBuildConfig, build_spatial_graph
from spatialtx_desktop.graph.context import ContextFieldConfig, add_context_field
from spatialtx_desktop.graph.continuous import continuous_edge_association
from spatialtx_desktop.graph.permutation import permutation_groups, permute_within_groups
from spatialtx_desktop.graph.runner import SpatialGraphAnalysisConfig, run_spatial_graph_neighborhood_batch
from spatialtx_desktop.graph.semantics import describe_obs_variable


class AuditAndSemanticsTests(unittest.TestCase):
    def test_input_audit_raw_counts_visium_and_sparse(self) -> None:
        adata = ad.AnnData(
            sparse.csr_matrix([[1, 0, 2], [0, 3, 0], [4, 0, 1]]),
            obs=pd.DataFrame({"array_row": [0, 0, 1], "array_col": [0, 2, 1], "in_tissue": [1, 1, 0]}),
            var=pd.DataFrame(index=["A", "B", "C"]),
        )
        adata.obsm["spatial"] = np.asarray([[0, 0], [2, 0], [1, 1]], dtype=float)
        adata.uns["spatial"] = {"sample": {"scalefactors": {"spot_diameter_fullres": 2.0}, "images": {"hires": np.zeros((2, 2, 3))}}}
        audit, table = audit_input(adata, input_filename="sample.h5ad", coordinate_unit="pixel")
        self.assertEqual(audit["expression_state_guess"], "raw_counts")
        self.assertEqual(audit["platform_guess"], "visium")
        self.assertEqual(audit["sparse_or_dense_X"], "sparse")
        self.assertEqual(audit["tissue_spot_count"], 2)
        self.assertTrue(audit["scalefactor_available"])
        self.assertIn("spatialtx_input_audit", adata.uns)
        self.assertEqual(len(table), 1)

    def test_input_audit_log_normalized_and_unknown_states(self) -> None:
        normalized = ad.AnnData(np.asarray([[0.0, 0.2], [1.1, 2.4], [0.3, 1.7]]))
        normalized.obsm["spatial"] = np.asarray([[0, 0], [1, 0], [2, 0]], dtype=float)
        audit, _ = audit_input(normalized)
        self.assertEqual(audit["expression_state_guess"], "log1p_normalized")
        self.assertIn("coordinate unit is unknown", audit["warnings"])
        unknown = ad.AnnData(np.full((3, 2), np.nan))
        unknown_audit, _ = audit_input(unknown)
        self.assertEqual(unknown_audit["expression_state_guess"], "unknown")
        self.assertIn("expression preprocessing state could not be determined", unknown_audit["warnings"])

    def test_four_variable_semantics(self) -> None:
        adata = ad.AnnData(
            np.ones((4, 1)),
            obs=pd.DataFrame({
                "state": ["A", "B", "A", None],
                "mask": [0, 1, 0, 1],
                "score": [-1.0, 0.2, 2.5, 4.0],
                "composition": [0.1, 0.3, 0.6, 0.9],
            }),
        )
        self.assertEqual(describe_obs_variable(adata, "state")["inferred_data_type"], "categorical_state")
        self.assertEqual(describe_obs_variable(adata, "mask")["inferred_data_type"], "binary_mask")
        self.assertEqual(describe_obs_variable(adata, "score")["inferred_data_type"], "continuous_score")
        self.assertEqual(describe_obs_variable(adata, "composition")["inferred_data_type"], "proportion_composition")
        self.assertIn("spatialtx_variable_semantics", adata.uns)

    def test_native_and_calibrated_radius_metadata(self) -> None:
        adata = ad.AnnData(np.ones((4, 1)))
        adata.obsm["spatial"] = np.asarray([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
        native = build_spatial_graph(adata, GraphBuildConfig(method="radius", radius=1.1, coordinate_unit="native"))
        self.assertFalse(native.metadata["physical_calibration_available"])
        self.assertEqual(native.metadata["radius_unit"], "native")
        self.assertIn("not calibrated physical distance", ";".join(native.warnings))
        calibrated = build_spatial_graph(adata, GraphBuildConfig(
            method="radius", radius=10.1, coordinate_unit="pixel", coordinate_scale=10.0, scale_source="micrometers_per_pixel"
        ))
        self.assertTrue(calibrated.metadata["physical_calibration_available"])
        self.assertEqual(calibrated.metadata["source_coordinate_unit"], "pixel")
        self.assertEqual(calibrated.metadata["radius_unit"], "micrometer")
        self.assertAlmostEqual(calibrated.metadata["effective_radius"], 10.1)
        self.assertAlmostEqual(float(calibrated.distances.data.min()), 10.0)


class PermutationAndContinuousTests(unittest.TestCase):
    @staticmethod
    def _two_component_graph() -> sparse.csr_matrix:
        rows = [0, 1, 2, 3]
        cols = [1, 0, 3, 2]
        return sparse.csr_matrix((np.ones(4), (rows, cols)), shape=(4, 4))

    def test_permutation_stays_within_connected_components(self) -> None:
        graph = self._two_component_graph()
        groups, valid, count = permutation_groups(graph, "within_connected_components")
        permuted = permute_within_groups(np.asarray([10, 11, 20, 21]), groups, np.random.default_rng(4), valid)
        self.assertEqual(count, 2)
        self.assertEqual(set(permuted[:2]), {10, 11})
        self.assertEqual(set(permuted[2:]), {20, 21})

    def test_permutation_stays_within_user_strata(self) -> None:
        graph = self._two_component_graph()
        strata = np.asarray(["left", "right", "left", "right"], dtype=object)
        groups, valid, count = permutation_groups(graph, "within_user_strata", strata)
        values = np.asarray([1, 2, 3, 4])
        permuted = permute_within_groups(values, groups, np.random.default_rng(2), valid)
        self.assertEqual(count, 2)
        self.assertEqual(set(permuted[[0, 2]]), {1, 3})
        self.assertEqual(set(permuted[[1, 3]]), {2, 4})

    def test_symmetric_generic_continuous_edge_statistic(self) -> None:
        graph = sparse.csr_matrix(
            (np.ones(4), ([0, 1, 1, 2], [1, 0, 2, 1])), shape=(3, 3)
        )
        kwargs = dict(
            statistic_name="I_XY", x_source="obs/X", y_source="obs/Y", permutations=19, seed=8,
            permutation_scope="whole_slide", active_graph="radius", graph_parameters={"radius": 1.0},
        )
        first = continuous_edge_association(graph, [1, 2, 3], [4, 5, 6], **kwargs)
        second = continuous_edge_association(graph, [1, 2, 3], [4, 5, 6], **kwargs)
        self.assertAlmostEqual(first["raw_weighted_statistic"], 10.0)
        self.assertEqual(first["empirical_p"], second["empirical_p"])
        self.assertGreaterEqual(first["empirical_p"], 1 / 20)
        self.assertTrue(sparse.issparse(graph))


class ContextAndRunnerTests(unittest.TestCase):
    def test_context_qc_warnings_and_low_coverage_gate(self) -> None:
        matrix = np.column_stack([np.arange(1, 9, dtype=float), np.ones(8), np.ones(8)])
        adata = ad.AnnData(matrix, var=pd.DataFrame(index=["CA9", "VEGFA", "OTHER"]))
        coverage, metadata = add_context_field(
            adata,
            ContextFieldConfig(field="H", genes=["CA9", "VEGFA"], min_coverage=0.5, dominant_gene_warning_fraction=0.55, library_correlation_warning=0.7),
        )
        self.assertEqual(len(coverage), 2)
        self.assertIn("one gene dominates", ";".join(metadata["warnings"]))
        self.assertIn("strongly correlated with library size", ";".join(metadata["warnings"]))
        with self.assertRaisesRegex(ValueError, "coverage"):
            add_context_field(adata, ContextFieldConfig(field="H", genes=["CA9", "MISSING1", "MISSING2"], min_coverage=0.8))
        with self.assertRaisesRegex(ValueError, "coverage"):
            add_context_field(adata, ContextFieldConfig(field="V", genes=["OTHER", "MISSING1", "MISSING2"], min_coverage=0.8))

    @staticmethod
    def _write_sample(path: Path, *, valid_spatial: bool = True) -> None:
        genes = ["C1", "C2", "S1", "S2"]
        matrix = np.arange(16 * 4, dtype=float).reshape(16, 4) % 9
        obs = pd.DataFrame({
            "state": ["A", "B"] * 8,
            "mask_a": [0, 1, 1, 0] * 4,
            "mask_b": [1, 0, 1, 0] * 4,
            "x_score": np.linspace(-1, 1, 16),
            "y_composition": np.linspace(0, 1, 16),
            "stratum": ["top"] * 8 + ["bottom"] * 8,
            "in_tissue": [1] * 15 + [0],
        })
        adata = ad.AnnData(matrix, obs=obs, var=pd.DataFrame(index=genes))
        if valid_spatial:
            adata.obsm["spatial"] = np.asarray([[x, y] for y in range(4) for x in range(4)], dtype=float)
        adata.write_h5ad(path)

    def test_runner_exports_audit_semantics_separate_questions_and_robustness(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "extended.h5ad"
            self._write_sample(source)
            config = SpatialGraphAnalysisConfig(
                graph=GraphBuildConfig(method="radius", radius=1.01, coordinate_unit="native"),
                enable_h=False,
                enable_v=False,
                label_source="state",
                label_mode="categorical_state",
                user_mask_a_column="mask_a",
                user_mask_b_column="mask_b",
                continuous_x_column="x_score",
                continuous_y_column="y_composition",
                continuous_x_mode="continuous_score",
                continuous_y_mode="proportion_composition",
                permutation_scope="within_user_strata",
                stratification_column="stratum",
                tissue_only_restriction=True,
                run_graph_robustness=True,
                permutations=9,
                seed=3,
            )
            run_dir, manifest = run_spatial_graph_neighborhood_batch(
                [source], root / "out", ["C1", "C2"], ["S1", "S2"], config
            )
            self.assertEqual(manifest.iloc[0]["status"], "ok")
            expected = [
                run_dir / "input_audit" / "input_audit_extended.json",
                run_dir / "input_audit" / "input_audit_extended.csv",
                run_dir / "variable_semantics" / "variable_semantics_extended.json",
                run_dir / "neighborhood" / "extended_same_spot_overlap.csv",
                run_dir / "neighborhood" / "extended_neighboring_spot_association.csv",
                run_dir / "neighborhood" / "extended_continuous_edge_statistics.csv",
                run_dir / "robustness" / "graph_robustness_summary_extended.csv",
                run_dir / "robustness" / "association_direction_stability_extended.csv",
            ]
            self.assertTrue(all(path.is_file() for path in expected))
            audit = json.loads(expected[0].read_text(encoding="utf-8"))
            self.assertFalse(audit["matrix_altered"])
            continuous = pd.read_csv(run_dir / "neighborhood" / "extended_continuous_edge_statistics.csv")
            self.assertIn("I_x_score_y_composition", set(continuous["statistic"]))

    def test_batch_continues_after_one_sample_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            good, bad = root / "good.h5ad", root / "bad.h5ad"
            self._write_sample(good)
            self._write_sample(bad, valid_spatial=False)
            config = SpatialGraphAnalysisConfig(
                graph=GraphBuildConfig(method="radius", radius=1.01), enable_h=False, enable_v=False, permutations=5
            )
            _, manifest = run_spatial_graph_neighborhood_batch(
                [bad, good], root / "out", ["C1", "C2"], ["S1", "S2"], config
            )
        self.assertEqual(len(manifest), 2)
        self.assertEqual(int(manifest["status"].eq("ok").sum()), 1)
        self.assertEqual(int(manifest["status"].str.startswith("error:").sum()), 1)

    def test_missing_h_and_v_are_skipped_without_losing_cs_statistics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "missing_context.h5ad"
            self._write_sample(source)
            config = SpatialGraphAnalysisConfig(
                graph=GraphBuildConfig(method="radius", radius=1.01),
                enable_h=True,
                enable_v=True,
                permutations=5,
            )
            run_dir, manifest = run_spatial_graph_neighborhood_batch(
                [source], root / "out", ["C1", "C2"], ["S1", "S2"], config
            )
            self.assertEqual(manifest.iloc[0]["status"], "ok")
            context = pd.read_csv(run_dir / "context_fields" / "missing_context_context_field_summary.csv")
            self.assertEqual(set(context["status"]), {"skipped_qc"})
            continuous = pd.read_csv(run_dir / "neighborhood" / "missing_context_continuous_edge_statistics.csv")
            self.assertEqual(continuous.loc[continuous["statistic"].eq("I_RR"), "status"].iloc[0], "ok")
            self.assertEqual(continuous.loc[continuous["statistic"].eq("I_HV"), "status"].iloc[0], "missing_field")


if __name__ == "__main__":
    unittest.main()
