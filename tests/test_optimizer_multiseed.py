from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from spatialtx_desktop.workflow import (
    _optimizer_seed_sequence,
    optimize_genes_multiseed,
    save_multiseed_optimizer_results,
)


class OptimizerMultiSeedTests(unittest.TestCase):
    def test_date_like_seed_sequence_advances_by_calendar_day(self) -> None:
        seeds, method = _optimizer_seed_sequence(20260624, 10)
        self.assertEqual(method, "calendar_day_sequence")
        self.assertEqual(
            seeds,
            [20260624, 20260625, 20260626, 20260627, 20260628, 20260629, 20260630, 20260701, 20260702, 20260703],
        )

    def test_multiseed_consensus_and_outputs(self) -> None:
        genes = ["C1", "C2", "C3", "C4", "S1", "S2"]
        coords = np.asarray([[x, y] for y in range(4) for x in range(4)], dtype=float)
        gradient = np.linspace(0.0, 1.0, len(coords))
        matrix = np.column_stack([
            8 * gradient + 1,
            7 * gradient + np.sin(np.arange(len(coords))),
            6 * gradient + np.cos(np.arange(len(coords))),
            5 * gradient[::-1] + 1,
            8 * gradient[::-1] + 1,
            7 * gradient[::-1] + np.cos(np.arange(len(coords))),
        ])
        adata = ad.AnnData(
            matrix,
            obs=pd.DataFrame(index=[f"spot{i}" for i in range(len(coords))]),
            var=pd.DataFrame(index=genes),
        )
        adata.obsm["spatial"] = coords
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "optimizer.h5ad"
            adata.write_h5ad(source)
            consensus, runs, frequency, overlap, details, summary = optimize_genes_multiseed(
                source,
                "C",
                ["C1", "C2"],
                ["S1", "S2"],
                k=2,
                pool_size=4,
                iterations=40,
                candidate_genes=["C1", "C2", "C3", "C4"],
                seeds=[101, 102, 103],
                consensus_threshold=2 / 3,
            )
            output = save_multiseed_optimizer_results(
                root, source.stem, "C", runs, frequency, overlap, details, summary
            )
            self.assertEqual(len(consensus), 2)
            self.assertEqual(len(runs), 3)
            self.assertEqual(set(runs["random_seed"]), {101, 102, 103})
            self.assertEqual(set(runs["qubo_energy_direction"]), {"lower_is_better"})
            self.assertEqual(len(overlap), 3)
            self.assertTrue(frequency["selection_frequency"].between(0, 1).all())
            self.assertEqual(summary["seed_count"], 3)
            self.assertIn("computational stability", summary["interpretation_limit"])
            self.assertTrue(summary["overlap_constraint_enabled"])
            self.assertEqual(summary["final_overlap_count"], 0)
            self.assertTrue(runs["overlap_constraint_enabled"].all())
            self.assertTrue(runs["final_overlap_count"].eq(0).all())
            for name in (
                "optimizer_multiseed_runs.csv",
                "optimizer_selection_frequency.csv",
                "optimizer_pairwise_overlap.csv",
                "optimizer_stability_summary.json",
                "optimizer_selection_frequency.png",
                "optimizer_energy_stability.png",
            ):
                self.assertTrue((output / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
