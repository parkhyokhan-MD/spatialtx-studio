from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from spatialtx_desktop.comparative.metrics import file_sha256
from spatialtx_desktop.comparative.models import ComparativeConfig, SampleRecord
from spatialtx_desktop.comparative.plotting import plot_side_by_side_maps
from spatialtx_desktop.comparative.runner import run_comparative_analysis
from spatialtx_desktop.workflow import score_h5ad
from tests.comparative_fixtures import C_GENES, S_GENES, write_comparative_h5ad


class ComparativeRegressionTests(unittest.TestCase):
    def test_side_by_side_legend_does_not_overlap_figure_footer(self) -> None:
        captured: dict[str, object] = {}

        def capture_figure(fig, path, config, figure_type, extra=None, **kwargs):
            captured["figure"] = fig
            footer = fig.text(0.01, 0.005, "SpatialTX Studio footer", fontsize=7)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            captured["legend_bottom"] = fig.legends[0].get_window_extent(renderer=renderer).y0
            captured["footer_top"] = footer.get_window_extent(renderer=renderer).y1
            return path

        coords = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        fields = {
            "coords": coords,
            "R": np.asarray([-1.0, -0.2, 0.3, 1.0]),
            "interface": np.asarray([False, True, False, False]),
            "diffuse": np.asarray([False, False, True, False]),
        }
        config = ComparativeConfig(
            mode="pairwise", reference="A", target="B", c_genes=["CD8A"], s_genes=["COL1A1"]
        )
        with tempfile.TemporaryDirectory() as tmp, patch(
            "spatialtx_desktop.comparative.plotting._save", side_effect=capture_figure
        ):
            plot_side_by_side_maps(
                {"sample_a": fields, "sample_b": fields},
                "sample_a",
                "sample_b",
                Path(tmp) / "side_by_side.png",
                config,
            )
        self.assertGreater(float(captured["legend_bottom"]), float(captured["footer_top"]) + 4.0)
        import matplotlib.pyplot as plt

        plt.close(captured["figure"])

    def test_comparative_run_does_not_modify_inputs_or_single_sample_scores(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = write_comparative_h5ad(root / "a.h5ad", pattern="localized", include_context=False)
            b = write_comparative_h5ad(root / "b.h5ad", pattern="diffuse", include_context=False)
            before_hashes = [file_sha256(a), file_sha256(b)]
            before_metrics, before_fields = score_h5ad(a, C_GENES, S_GENES)
            run_comparative_analysis(
                [SampleRecord("a", a, "A"), SampleRecord("b", b, "B")],
                ComparativeConfig(
                    mode="pairwise", reference="A", target="B", c_genes=C_GENES, s_genes=S_GENES,
                    enable_h_expr=False, enable_v_expr=False, bootstrap_iterations=100,
                ),
                root / "out",
            )
            after_metrics, after_fields = score_h5ad(a, C_GENES, S_GENES)
            self.assertEqual(before_hashes, [file_sha256(a), file_sha256(b)])
            for field in ("C", "S", "R", "G"):
                np.testing.assert_allclose(before_fields[field], after_fields[field])
            for mask in ("interface", "diffuse"):
                np.testing.assert_array_equal(before_fields[mask], after_fields[mask])
            for metric in ("regime_label", "transition_burden_score", "interface_fraction", "diffuse_fraction"):
                self.assertEqual(before_metrics[metric], after_metrics[metric])


if __name__ == "__main__":
    unittest.main()
