from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from spatialtx_desktop.comparative.metrics import file_sha256
from spatialtx_desktop.comparative.models import ComparativeConfig, SampleRecord
from spatialtx_desktop.comparative.runner import run_comparative_analysis
from spatialtx_desktop.workflow import score_h5ad
from tests.comparative_fixtures import C_GENES, S_GENES, write_comparative_h5ad


class ComparativeRegressionTests(unittest.TestCase):
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
