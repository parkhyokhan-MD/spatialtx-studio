from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from spatialtx_desktop.comparative.models import ComparativeConfig, NO_REGISTRATION_NOTICE, SampleRecord
from spatialtx_desktop.comparative.runner import run_comparative_analysis
from tests.comparative_fixtures import C_GENES, S_GENES, write_comparative_h5ad


class ComparativePairwiseTests(unittest.TestCase):
    def test_pairwise_synthetic_delta_direction_and_no_spot_subtraction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference = write_comparative_h5ad(root / "reference.h5ad", pattern="localized", seed=3)
            target = write_comparative_h5ad(root / "target.h5ad", pattern="diffuse", seed=4)
            records = [
                SampleRecord("reference", reference, "A"),
                SampleRecord("target", target, "B"),
            ]
            result = run_comparative_analysis(
                records,
                ComparativeConfig(
                    mode="pairwise",
                    reference="A",
                    target="B",
                    c_genes=C_GENES,
                    s_genes=S_GENES,
                    h_genes=["CA9", "VEGFA"],
                    v_genes=["PECAM1", "VWF"],
                    bootstrap_iterations=100,
                    seed=17,
                ),
                root / "results",
            )
            delta = result.delta_metrics.set_index("metric")
            samples = result.sample_metrics.set_index("sample_id")
            expected = samples.loc["target", "transition_burden_score"] - samples.loc["reference", "transition_burden_score"]
            self.assertAlmostEqual(float(delta.loc["transition_burden_score", "delta"]), float(expected))
            self.assertEqual(delta.loc["transition_burden_score", "direction_definition"], "Target - Reference")
            self.assertTrue((result.sample_metrics["registration_status"] == "not_performed").all())
            self.assertTrue((result.delta_metrics["registration_notice"] == NO_REGISTRATION_NOTICE).all())
            self.assertTrue((result.sample_metrics["H_V_core_effect"].str.startswith("none")).all())
            self.assertTrue((result.run_dir / "pairwise_delta_summary.csv").is_file())
            self.assertTrue((result.run_dir / "comparative_summary_report.pdf").is_file())
            self.assertFalse(np.isnan(float(delta.loc["H_expr_mean", "delta"])))
            self.assertTrue(result.metric_change_table.loc[
                result.metric_change_table["metric_name"].isin(["H_expr_mean", "V_expr_mean"]),
                "interpretation_flag",
            ].eq("non_informative_centered_mean").all())
            self.assertTrue(result.hv_summary["observational_only"].astype(bool).all())
            self.assertTrue(result.hv_summary.loc[
                result.hv_summary["metric"].isin(["H_raw_median", "V_raw_median"]), "raw_delta"
            ].notna().all())
            self.assertTrue(result.sample_metrics["H_high_threshold_method"].eq("pooled_reference_target_q90").all())
            self.assertTrue((result.sample_metrics["H_high_fraction"].between(0, 1)).all())

    def test_same_seed_and_inputs_are_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = write_comparative_h5ad(root / "a.h5ad", pattern="localized", include_context=False, seed=1)
            b = write_comparative_h5ad(root / "b.h5ad", pattern="diffuse", include_context=False, seed=2)
            records = [SampleRecord("a", a, "A"), SampleRecord("b", b, "B")]
            config = ComparativeConfig(
                mode="pairwise", reference="A", target="B", c_genes=C_GENES, s_genes=S_GENES,
                enable_h_expr=False, enable_v_expr=False, bootstrap_iterations=100, seed=42,
            )
            first = run_comparative_analysis(records, config, root / "results")
            second = run_comparative_analysis(records, config, root / "results")
            first_delta = first.delta_metrics[["metric", "delta"]].reset_index(drop=True)
            second_delta = second.delta_metrics[["metric", "delta"]].reset_index(drop=True)
            self.assertTrue(first.sample_metrics["cache_status"].eq("computed").all())
            self.assertTrue(second.sample_metrics["cache_status"].eq("reused").all())
            self.assertEqual(first_delta["metric"].tolist(), second_delta["metric"].tolist())
            np.testing.assert_allclose(first_delta["delta"], second_delta["delta"], equal_nan=True)


if __name__ == "__main__":
    unittest.main()
