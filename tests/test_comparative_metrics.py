from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from spatialtx_desktop.comparative.metrics import compute_delta_metrics, regime_transition_table


class ComparativeMetricTests(unittest.TestCase):
    def test_delta_is_target_minus_reference_for_known_values(self) -> None:
        samples = pd.DataFrame([
            {"sample_id": "A", "C_mean": 1.0, "S_mean": 4.0, "R_mean": -3.0,
             "C_median": .2, "S_median": .8, "R_median": -.6, "R_std": 1.0,
             "gradient_mean": 0.5, "gradient_q90": 1.0, "localized_interface_fraction": .2,
             "diffuse_fraction": .1, "transition_burden_score": .15, "adj_same_fraction": .6,
             "adj_zero_fraction": .2, "adj_opposite_fraction": .2, "R_crossing_fraction": .1,
             "interface_fragmentation_index": .2, "n_diffuse_components": 2, "small_component_fraction": .1,
             "largest_diffuse_component_ratio": .8, "H_expr_mean": .2, "V_expr_mean": .3},
            {"sample_id": "B", "C_mean": 3.0, "S_mean": 2.0, "R_mean": 1.0,
             "C_median": .7, "S_median": .3, "R_median": .4, "R_std": 2.5,
             "gradient_mean": 1.5, "gradient_q90": 3.0, "localized_interface_fraction": .1,
             "diffuse_fraction": .4, "transition_burden_score": .45, "adj_same_fraction": .3,
             "adj_zero_fraction": .4, "adj_opposite_fraction": .3, "R_crossing_fraction": .4,
             "interface_fragmentation_index": .5, "n_diffuse_components": 5, "small_component_fraction": .4,
             "largest_diffuse_component_ratio": .3, "H_expr_mean": .7, "V_expr_mean": .1},
        ])
        matches = pd.DataFrame([{
            "comparison_id": "A_vs_B", "pair_id": "", "reference_sample_id": "A", "target_sample_id": "B"
        }])
        delta = compute_delta_metrics(samples, matches).set_index("delta_metric")
        self.assertEqual(delta.loc["delta_C", "delta"], 2.0)
        self.assertEqual(delta.loc["delta_S", "delta"], -2.0)
        self.assertEqual(delta.loc["delta_R", "delta"], 4.0)
        self.assertAlmostEqual(delta.loc["delta_C_median", "delta"], .5)
        self.assertAlmostEqual(delta.loc["delta_S_median", "delta"], -.5)
        self.assertAlmostEqual(delta.loc["delta_R_median", "delta"], 1.0)
        self.assertAlmostEqual(delta.loc["delta_transition_burden_score", "delta"], .30)
        self.assertAlmostEqual(delta.loc["delta_H_expr", "delta"], .5)

    def test_regime_transition_and_low_confidence_flag(self) -> None:
        samples = pd.DataFrame([
            {"sample_id": "A", "regime_label": "Type_A_candidate", "regime_confidence": .05, "spatial_qc_status": "PASS"},
            {"sample_id": "B", "regime_label": "Type_B_candidate", "regime_confidence": .20, "spatial_qc_status": "PASS"},
        ])
        matches = pd.DataFrame([{
            "comparison_id": "A_vs_B", "pair_id": "", "reference_sample_id": "A", "target_sample_id": "B"
        }])
        table = regime_transition_table(samples, matches, "A", "B", .10)
        self.assertEqual(table.iloc[0]["regime_transition"], "Type_A_candidate → Type_B_candidate")
        self.assertEqual(table.iloc[0]["transition_confidence_flag"], "uncertain")
        self.assertIn("not a validated biological", table.iloc[0]["interpretation"])

    def test_missing_optional_context_delta_is_retained_as_unavailable(self) -> None:
        rows = []
        for sample_id in ("A", "B"):
            row = {"sample_id": sample_id}
            for metric in (
                "C_mean", "S_mean", "R_mean", "R_std", "gradient_mean", "gradient_q90",
                "localized_interface_fraction", "diffuse_fraction", "transition_burden_score",
                "adj_same_fraction", "adj_zero_fraction", "adj_opposite_fraction", "R_crossing_fraction",
                "interface_fragmentation_index", "n_diffuse_components", "small_component_fraction",
                "largest_diffuse_component_ratio",
            ):
                row[metric] = 1.0
            row["H_expr_mean"] = np.nan
            row["V_expr_mean"] = np.nan
            rows.append(row)
        matches = pd.DataFrame([{"comparison_id": "x", "pair_id": "", "reference_sample_id": "A", "target_sample_id": "B"}])
        delta = compute_delta_metrics(pd.DataFrame(rows), matches)
        self.assertEqual(delta.loc[delta["delta_metric"].eq("delta_H_expr"), "status"].iloc[0], "unavailable_optional_metric")


if __name__ == "__main__":
    unittest.main()
