from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd

from spatialtx_desktop.advanced_analysis import run_advanced_batch
from spatialtx_desktop.app import SpatialTXDesktop
from spatialtx_desktop.gene_program_validation import (
    GeneProgramOverlapError,
    find_gene_program_overlap,
    normalize_gene_list,
    validate_gene_programs,
)
from spatialtx_desktop.graph.builder import GraphBuildConfig
from spatialtx_desktop.graph.runner import SpatialGraphAnalysisConfig, run_spatial_graph_neighborhood_batch
from spatialtx_desktop.workflow import optimize_genes, run_batch, score_adata
from spatialtx_studio.frame26 import run_frame26
from spatialtx_studio.gene_program import select_gene_programs


def _example_adata() -> ad.AnnData:
    genes = [
        "GZMB", "IFNG", "PRF1", "CD8A", "CD8B", "NKG7", "MS4A1", "IGHA1",
        "IGLC3", "TAGLN", "SPP1", "PDGFRA", "ACTA2", "FN1", "COL5A1",
    ]
    coords = np.asarray([[x, y] for y in range(4) for x in range(4)], dtype=float)
    gradient = np.linspace(0.0, 1.0, len(coords))
    matrix = np.column_stack([
        (index + 1) * gradient + (len(genes) - index) * gradient[::-1] * 0.1 + 1
        for index in range(len(genes))
    ])
    result = ad.AnnData(
        matrix,
        obs=pd.DataFrame(index=[f"spot_{index}" for index in range(len(coords))]),
        var=pd.DataFrame(index=genes),
    )
    result.obsm["spatial"] = coords
    return result


class _FakeText:
    def __init__(self, value: str) -> None:
        self.value = value

    def get(self, *_args) -> str:
        return self.value


class _FakeLabel:
    def __init__(self) -> None:
        self.values: dict = {}

    def configure(self, **kwargs) -> None:
        self.values.update(kwargs)


class GeneProgramValidatorTests(unittest.TestCase):
    def test_case_and_whitespace_overlap_are_canonicalized(self) -> None:
        self.assertEqual(normalize_gene_list([" cd8a ", "CD8A", "nkg7", ""]), ["CD8A", "NKG7"])
        self.assertEqual(find_gene_program_overlap([" IGHA1 "], ["igha1"]), ["IGHA1"])

    def test_same_program_duplicates_are_removed_in_order(self) -> None:
        result = validate_gene_programs(
            [" cd8a ", "CD8A", "NKG7", "nkg7", "PRF1"],
            ["COL1A1", "LUM"],
            mode="custom",
        )
        self.assertEqual(result.normalized_c_genes, ["CD8A", "NKG7", "PRF1"])
        self.assertEqual(result.c_duplicates_removed, ["CD8A", "NKG7"])
        self.assertEqual(result.validation_status, "valid")

    def test_custom_igha1_example_is_blocked_but_iglc3_is_not_overlap(self) -> None:
        c_genes = "GZMB, IFNG, PRF1, CD8A, CD8B, NKG7, MS4A1, IGHA1"
        s_genes = "IGLC3, TAGLN, SPP1, PDGFRA, ACTA2, IGHA1, FN1, COL5A1"
        with self.assertRaises(GeneProgramOverlapError) as caught:
            validate_gene_programs(c_genes, s_genes, mode="custom")
        self.assertEqual(caught.exception.result.overlap_genes, ["IGHA1"])
        self.assertNotIn("IGLC3", caught.exception.result.overlap_genes)
        self.assertIn("R = C - S", str(caught.exception))

    def test_fixed_overlap_raises_development_error(self) -> None:
        with self.assertRaisesRegex(GeneProgramOverlapError, "Fixed gene programs must be mutually exclusive"):
            validate_gene_programs(["CD8A", "IGHA1"], ["COL1A1", "igha1"], mode="fixed")


class GeneProgramIntegrationTests(unittest.TestCase):
    C_VALID = ["GZMB", "IFNG", "PRF1", "CD8A", "CD8B", "NKG7", "MS4A1", "IGHA1"]
    S_VALID = ["IGLC3", "TAGLN", "SPP1", "PDGFRA", "ACTA2", "FN1", "COL5A1"]

    def test_core_score_blocks_overlap(self) -> None:
        with self.assertRaises(GeneProgramOverlapError):
            score_adata(_example_adata(), ["IGHA1"], ["igha1"])

    def test_main_batch_blocks_overlap_before_output_creation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "sample.h5ad"
            _example_adata().write_h5ad(source)
            output = root / "results"
            with self.assertRaises(GeneProgramOverlapError):
                run_batch([source], output, ["IGHA1"], ["igha1"])
            self.assertFalse(output.exists())

    def test_gui_gene_preflight_blocks_overlap(self) -> None:
        fake = SimpleNamespace(
            c_text=_FakeText("CD8A, IGHA1"),
            s_text=_FakeText("COL5A1, igha1"),
            gene_program_warning=_FakeLabel(),
        )
        with self.assertRaises(GeneProgramOverlapError):
            SpatialTXDesktop._genes(fake)
        self.assertIn("IGHA1", fake.gene_program_warning.values["text"])

    def test_gui_shows_nonblocking_s_side_immunoglobulin_warning(self) -> None:
        fake = SimpleNamespace(
            c_text=_FakeText("CD8A, NKG7"),
            s_text=_FakeText("IGLC3, TAGLN"),
            gene_program_warning=_FakeLabel(),
        )
        c_genes, s_genes = SpatialTXDesktop._genes(fake)
        self.assertEqual(c_genes, ["CD8A", "NKG7"])
        self.assertEqual(s_genes, ["IGLC3", "TAGLN"])
        self.assertIn("immune-associated", fake.gene_program_warning.values["text"])

    def test_fixed_config_overlap_is_treated_as_development_error(self) -> None:
        config = {
            "selection": {"topK": 2, "min_genes_per_program": 1},
            "gene_programs": {
                "C_FIXED": ["CD8A", "IGHA1"],
                "B_FIXED": ["FN1", "igha1"],
            },
        }
        with self.assertRaisesRegex(GeneProgramOverlapError, "Fixed gene programs must be mutually exclusive"):
            select_gene_programs(_example_adata(), config, "fixed")

    def test_cli_custom_mode_blocks_overlap(self) -> None:
        config = {
            "selection": {"topK": 2, "min_genes_per_program": 1},
            "gene_programs": {"CUSTOM_C": ["IGHA1"], "CUSTOM_B": ["igha1"]},
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "sample.h5ad"
            _example_adata().write_h5ad(source)
            with self.assertRaises(GeneProgramOverlapError):
                run_frame26(source, "sample", root / "cli", mode="custom", config=config)

    def test_advanced_batch_blocks_overlap_before_output_creation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "sample.h5ad"
            _example_adata().write_h5ad(source)
            output = root / "advanced"
            with self.assertRaises(GeneProgramOverlapError):
                run_advanced_batch("composition", [source], output, ["IGHA1"], ["igha1"])
            self.assertFalse(output.exists())

    def test_spatial_graph_blocks_overlap_before_output_creation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "sample.h5ad"
            _example_adata().write_h5ad(source)
            output = root / "graph"
            config = SpatialGraphAnalysisConfig(
                graph=GraphBuildConfig(method="knn", k=3),
                enable_h=False,
                enable_v=False,
                permutations=1,
            )
            with self.assertRaises(GeneProgramOverlapError):
                run_spatial_graph_neighborhood_batch(
                    [source], output, ["IGHA1"], ["igha1"], config
                )
            self.assertFalse(output.exists())

    def test_gene_composition_csv_uses_validated_disjoint_programs_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "sample.h5ad"
            _example_adata().write_h5ad(source)
            c_requested = [" GZMB ", "gzmb", *self.C_VALID[1:]]
            run_dir, manifest = run_advanced_batch(
                "composition", [source], root / "advanced", c_requested, self.S_VALID
            )
            self.assertEqual(manifest.iloc[0]["status"], "ok")
            table = pd.read_csv(run_dir / "sample" / "gene_composition.csv")
            program_counts = table.groupby("gene")["program"].nunique()
            self.assertTrue((program_counts == 1).all())
            for column in (
                "c_gene_count_requested", "s_gene_count_requested", "c_gene_count_used",
                "s_gene_count_used", "n_overlap_genes", "overlap_genes",
                "overlap_policy", "program_validation_status",
            ):
                self.assertIn(column, table.columns)
            self.assertTrue(table["n_overlap_genes"].eq(0).all())
            self.assertTrue(table["program_validation_status"].eq("valid").all())
            self.assertTrue(table["c_gene_count_requested"].eq(len(c_requested)).all())
            self.assertTrue(table["c_gene_count_used"].eq(len(self.C_VALID)).all())

    def test_provenance_records_requested_used_and_duplicate_removal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "sample.h5ad"
            _example_adata().write_h5ad(source)
            run_dir, summary = run_batch(
                [source],
                root / "results",
                [" CD8A ", "cd8a", "IGHA1"],
                ["IGLC3", "TAGLN"],
            )
            self.assertEqual(summary.iloc[0]["gene_program_validation_status"], "valid")
            run_config = json.loads((run_dir / "run_config.json").read_text(encoding="utf-8"))
            provenance = run_config["gene_program_validation"]
            self.assertEqual(provenance["c_genes_used"], ["CD8A", "IGHA1"])
            self.assertEqual(provenance["c_duplicates_removed"], ["CD8A"])
            self.assertEqual(provenance["n_overlap_genes"], 0)
            parameter = json.loads((run_dir / "sample" / "parameter_log.json").read_text(encoding="utf-8"))
            self.assertEqual(parameter["gene_program_validation"]["validation_status"], "valid")

    def test_adaptive_selection_excludes_opposite_side_gene(self) -> None:
        adata = _example_adata()
        config = {
            "selection": {"topK": 3, "min_genes_per_program": 1},
            "gene_programs": {
                "C_CANDIDATES": ["CD8A", "IGHA1", "NKG7"],
                "B_CANDIDATES": ["igha1", "IGLC3", "TAGLN", "FN1"],
            },
        }
        selected_c, selected_s, note = select_gene_programs(adata, config, "cancer_adaptive")
        self.assertFalse(set(selected_c) & set(selected_s))
        self.assertNotIn("IGHA1", selected_s)
        self.assertIn("genes_excluded_due_to_opposite_side=IGHA1", note)

    def test_qubo_excludes_opposite_side_candidate_and_records_constraint(self) -> None:
        genes = ["C1", "C2", "S1"]
        coords = np.asarray([[x, y] for y in range(4) for x in range(4)], dtype=float)
        gradient = np.linspace(0.0, 1.0, len(coords))
        matrix = np.column_stack([gradient + 1, gradient * 0.8 + 1, gradient[::-1] + 1])
        adata = ad.AnnData(matrix, var=pd.DataFrame(index=genes))
        adata.obsm["spatial"] = coords
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "qubo.h5ad"
            adata.write_h5ad(source)
            selected, detail, summary = optimize_genes(
                source,
                "C",
                ["C1"],
                ["S1"],
                k=2,
                pool_size=3,
                iterations=20,
                candidate_genes=["C1", "S1", "C2"],
                seed=7,
            )
        self.assertEqual(set(selected), {"C1", "C2"})
        self.assertNotIn("S1", set(detail["gene"]))
        self.assertEqual(summary["genes_excluded_due_to_opposite_side"], "S1")
        self.assertTrue(summary["overlap_constraint_enabled"])
        self.assertEqual(summary["final_overlap_count"], 0)

    def test_nonoverlap_canonicalization_preserves_scores_and_masks(self) -> None:
        adata = _example_adata()
        clean_metrics, clean_fields = score_adata(
            adata, ["CD8A", "IGHA1"], ["IGLC3", "TAGLN"]
        )
        normalized_metrics, normalized_fields = score_adata(
            adata,
            [" cd8a ", "CD8A", "igha1"],
            [" iglc3 ", "TAGLN", "tagln"],
        )
        for key in ("C", "S", "R", "G", "interface", "diffuse"):
            np.testing.assert_allclose(clean_fields[key], normalized_fields[key], equal_nan=True)
        for key in ("regime_label", "public_transition_pattern", "transition_burden_score"):
            self.assertEqual(clean_metrics[key], normalized_metrics[key])


if __name__ == "__main__":
    unittest.main()
