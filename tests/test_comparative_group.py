from __future__ import annotations

import json
import unittest
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from spatialtx_desktop.comparative.matching import paired_group_matches
from spatialtx_desktop.comparative.models import ComparativeConfig, SampleRecord
from spatialtx_desktop.comparative.runner import run_comparative_analysis
from spatialtx_desktop.comparative.statistics import comparative_group_statistics
from tests.comparative_fixtures import C_GENES, S_GENES, write_comparative_h5ad


class ComparativeGroupTests(unittest.TestCase):
    @staticmethod
    def _config(mode: str) -> ComparativeConfig:
        return ComparativeConfig(
            mode=mode, reference="pre", target="post", c_genes=["C1"], s_genes=["S1"],
            enable_h_expr=False, enable_v_expr=False, bootstrap_iterations=100, seed=9,
        )

    def test_paired_matching_and_incomplete_pair_rejection(self) -> None:
        records = [
            SampleRecord("p1_pre", Path("p1_pre.h5ad"), "pre", pair_id="p1"),
            SampleRecord("p1_post", Path("p1_post.h5ad"), "post", pair_id="p1"),
            SampleRecord("p2_pre", Path("p2_pre.h5ad"), "pre", pair_id="p2"),
            SampleRecord("p2_post", Path("p2_post.h5ad"), "post", pair_id="p2"),
        ]
        matches = paired_group_matches(records, "pre", "post")
        self.assertEqual(len(matches), 2)
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            paired_group_matches(records[:-1], "pre", "post")

    def test_paired_statistics_uses_matched_values(self) -> None:
        table = pd.DataFrame([
            {"sample_id": "p1_pre", "group": "pre", "R_mean": 1.0},
            {"sample_id": "p1_post", "group": "post", "R_mean": 3.0},
            {"sample_id": "p2_pre", "group": "pre", "R_mean": 2.0},
            {"sample_id": "p2_post", "group": "post", "R_mean": 5.0},
            {"sample_id": "p3_pre", "group": "pre", "R_mean": 4.0},
            {"sample_id": "p3_post", "group": "post", "R_mean": 8.0},
        ])
        matches = pd.DataFrame([
            {"reference_sample_id": f"p{i}_pre", "target_sample_id": f"p{i}_post"} for i in (1, 2, 3)
        ])
        stats = comparative_group_statistics(table, self._config("paired"), matches, effective_mode="paired")
        row = stats.loc[stats["metric"].eq("R_mean")].iloc[0]
        self.assertEqual(row["test"], "wilcoxon_signed_rank")
        self.assertGreater(row["mean_difference_target_minus_reference"], 0)
        self.assertGreater(row["effect_size"], 0)

    def test_unpaired_statistics_effect_size_fdr_and_direction(self) -> None:
        table = pd.DataFrame([
            *[{"sample_id": f"A{i}", "group": "pre", "R_mean": value} for i, value in enumerate([0., 1., 2., 1.])],
            *[{"sample_id": f"B{i}", "group": "post", "R_mean": value} for i, value in enumerate([5., 6., 7., 6.])],
        ])
        stats = comparative_group_statistics(table, self._config("unpaired"), None, effective_mode="unpaired")
        row = stats.loc[stats["metric"].eq("R_mean")].iloc[0]
        self.assertEqual(row["test"], "mann_whitney_u")
        self.assertGreater(row["effect_size"], 0)
        self.assertGreater(row["mean_difference_target_minus_reference"], 0)
        self.assertTrue(np.isfinite(row["adjusted_p_value_bh"]))

    def test_paired_and_unpaired_group_runs_complete_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = []
            for pair_number in (1, 2):
                pre = write_comparative_h5ad(
                    root / f"p{pair_number}_pre.h5ad", pattern="localized", include_context=False, seed=pair_number
                )
                post = write_comparative_h5ad(
                    root / f"p{pair_number}_post.h5ad", pattern="diffuse", include_context=False, seed=pair_number + 10
                )
                records.extend([
                    SampleRecord(f"p{pair_number}_pre", pre, "pre", pair_id=f"p{pair_number}"),
                    SampleRecord(f"p{pair_number}_post", post, "post", pair_id=f"p{pair_number}"),
                ])
            paired = run_comparative_analysis(
                records,
                ComparativeConfig(
                    mode="paired", reference="pre", target="post", c_genes=C_GENES, s_genes=S_GENES,
                    enable_h_expr=False, enable_v_expr=False, bootstrap_iterations=100, seed=11,
                ),
                root / "paired_results",
            )
            self.assertEqual(paired.effective_mode, "paired")
            self.assertEqual(len(paired.regime_transitions), 2)
            self.assertTrue((paired.run_dir / "group_effect_sizes.csv").is_file())
            self.assertTrue(paired.run_manifest["status"].eq("ok").all())
            paired_regime_meta = json.loads(
                (paired.run_dir / "comparative_figures" / "comparative_regime_transitions.png.metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(paired_regime_meta["visualization_mode"], "matched_group_3x3_transition_matrix")

            unpaired_records = [
                SampleRecord(record.sample_id, record.file_path, record.group) for record in records
            ]
            unpaired = run_comparative_analysis(
                unpaired_records,
                ComparativeConfig(
                    mode="unpaired", reference="pre", target="post", c_genes=C_GENES, s_genes=S_GENES,
                    enable_h_expr=False, enable_v_expr=False, bootstrap_iterations=100, seed=11,
                ),
                root / "unpaired_results",
            )
            self.assertEqual(unpaired.effective_mode, "unpaired")
            self.assertTrue(unpaired.regime_transitions["comparison_basis"].eq("unpaired_group_distribution").all())
            self.assertTrue(unpaired.group_statistics["effect_size_method"].astype(str).str.len().gt(0).any())
            self.assertTrue((unpaired.run_dir / "group_fdr_results.csv").is_file())
            unpaired_regime_meta = json.loads(
                (unpaired.run_dir / "comparative_figures" / "comparative_regime_transitions.png.metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(unpaired_regime_meta["visualization_mode"], "unpaired_regime_distribution")


if __name__ == "__main__":
    unittest.main()
