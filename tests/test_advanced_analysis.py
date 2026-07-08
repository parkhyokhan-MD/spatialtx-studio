from __future__ import annotations

import inspect
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from spatialtx_desktop import __version__
from spatialtx_desktop import advanced_analysis as aa
from spatialtx_desktop.advanced_analysis_ui import AdvancedAnalysisPanel
from spatialtx_desktop.workflow import DEFAULT_C_GENES, DEFAULT_S_GENES, score_h5ad


class FakeAdata:
    def __init__(self) -> None:
        self.X = np.asarray([
            [8, 4, 0, 1],
            [6, 3, 1, 2],
            [2, 1, 6, 4],
            [1, 0, 8, 5],
        ], dtype=float)
        self.var_names = np.asarray(["CD8A", "NKG7", "COL1A1", "LUM"])
        self.n_obs, self.n_vars = self.X.shape


class AdvancedAnalysisTests(unittest.TestCase):
    def test_version_and_v01_defaults_are_preserved(self) -> None:
        self.assertEqual(__version__, "0.3-beta")
        self.assertEqual(DEFAULT_C_GENES, ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG"])
        self.assertEqual(DEFAULT_S_GENES, ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN"])
        signature = inspect.signature(score_h5ad)
        self.assertEqual(signature.parameters["c_q"].default, .80)
        self.assertEqual(signature.parameters["s_q"].default, .80)
        self.assertEqual(signature.parameters["g_q"].default, .60)

    def test_gene_composition_lists_missing_genes_and_sums_by_program(self) -> None:
        with patch.object(aa, "_read_h5ad", return_value=FakeAdata()):
            table, metadata = aa.calculate_gene_composition(
                Path("sample.h5ad"), ["CD8A", "NKG7", "MISSING_C"], ["COL1A1", "LUM", "MISSING_S"]
            )
        self.assertEqual(len(table), 6)
        self.assertEqual(set(table.loc[table["status"] == "missing", "gene"]), {"MISSING_C", "MISSING_S"})
        for program in ("Cx", "Sx"):
            self.assertAlmostEqual(table.loc[table["program"] == program, "relative_contribution_percent"].sum(), 100.0)
        self.assertEqual(metadata["Cx_genes_requested"][-1], "MISSING_C")

    def test_bh_adjustment_is_monotone_in_rank(self) -> None:
        raw = np.asarray([0.04, 0.001, np.nan, 0.02, 0.50])
        adjusted = aa._bh_adjust(raw)
        order = np.argsort(raw[np.isfinite(raw)])
        ranked = adjusted[np.flatnonzero(np.isfinite(raw))][order]
        self.assertTrue(np.all(np.diff(ranked) >= -1e-12))
        self.assertTrue(np.isnan(adjusted[2]))

    def test_interaction_metrics_distinguish_coexistence_and_antagonism(self) -> None:
        edges = [(0, 1), (1, 2), (2, 3)]
        coexist = aa._interaction_values(np.asarray([.8, .8, .7, .7]), np.asarray([.8, .7, .8, .7]), edges)
        oppose = aa._interaction_values(np.asarray([1., 1., 0., 0.]), np.asarray([0., 0., 1., 1.]), edges)
        self.assertGreater(coexist["coexistence_index"], oppose["coexistence_index"])
        self.assertLess(coexist["antagonism_index"], oppose["antagonism_index"])
        self.assertGreater(coexist["spatial_overlap_index"], oppose["spatial_overlap_index"])

    def test_seeded_permutation_is_reproducible(self) -> None:
        values = np.linspace(0, 1, 12)
        adjacency = [[1], [0, 2], [1, 3], [2, 4], [3, 5], [4, 6], [5, 7], [6, 8], [7, 9], [8, 10], [9, 11], [10]]
        rng_a, rng_b = np.random.default_rng(20260705), np.random.default_rng(20260705)
        first = aa._local_mean(rng_a.permutation(values), adjacency)
        second = aa._local_mean(rng_b.permutation(values), adjacency)
        np.testing.assert_allclose(first, second)

    def test_dashboard_summaries_cover_all_three_modules(self) -> None:
        panel = AdvancedAnalysisPanel.__new__(AdvancedAnalysisPanel)
        composition = aa.pd.DataFrame([
            {"program": "Cx", "gene": "CD8A", "status": "present", "relative_contribution_percent": 60.0},
            {"program": "Cx", "gene": "NKG7", "status": "present", "relative_contribution_percent": 40.0},
            {"program": "Sx", "gene": "LUM", "status": "present", "relative_contribution_percent": 100.0},
        ])
        enrichment = aa.pd.DataFrame([
            {"status": "present", "gene": "CD8A", "n_interface": 12, "n_noninterface": 88, "hedges_g": 1.25, "significant_fdr_0_05": True},
        ])
        interaction = aa.pd.DataFrame([{
            "coexistence_index": .31, "antagonism_index": .42, "spatial_overlap_index": .37, "balance_index": .55,
        }])
        self.assertIn("CD8A", panel._summary_metrics("composition", composition)[0])
        self.assertIn("Interface 12", panel._summary_metrics("enrichment", enrichment)[0])
        self.assertIn("Coexistence", panel._summary_metrics("interaction", interaction)[0])


if __name__ == "__main__":
    unittest.main()
