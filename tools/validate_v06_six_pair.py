from __future__ import annotations

import argparse
from pathlib import Path

from spatialtx_desktop import __version__
from spatialtx_desktop.comparative.models import ComparativeConfig
from spatialtx_desktop.comparative.multi_pair import PairSpec, run_multi_pair_analysis


C_GENES = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG"]
S_GENES = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN"]
PAIR_IDS = (30, 38, 41, 42, 43, 44)


def _one_match(data_root: Path, pattern: str) -> Path:
    matches = sorted(data_root.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one file for {pattern!r}; found {len(matches)}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the preserved six-pair BTC regression set.")
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--context", action="store_true", help="Enable the existing H/V context programs.")
    args = parser.parse_args()

    data_root = args.data_root.expanduser().resolve()
    pairs = [
        PairSpec(
            f"sample_{pair_id}",
            _one_match(data_root, f"*_sample_{pair_id}_pre.h5ad"),
            _one_match(data_root, f"*_sample_{pair_id}_post.h5ad"),
        )
        for pair_id in PAIR_IDS
    ]
    config = ComparativeConfig(
        mode="pairwise",
        reference="Pre",
        target="Post",
        c_genes=C_GENES,
        s_genes=S_GENES,
        c_q=0.8,
        s_q=0.8,
        g_q=0.6,
        scoring_options={
            "smoothing_mode": "none",
            "smoothing_k": 6,
            "gaussian_sigma": 0.0,
            "normalization_mode": "raw_mean",
            "perturbation_check": False,
            "c_q_list": [0.75, 0.80, 0.85],
            "s_q_list": [0.75, 0.80, 0.85],
            "g_q_list": [0.50, 0.60, 0.70],
            "parameter_log_export": True,
            "dense_warning_gb": 4.0,
        },
        graph_settings={"method": "knn", "k": 6, "weighting": "binary"},
        enable_h_expr=args.context,
        enable_v_expr=args.context,
        context_smoothing="none",
        seed=42,
        bootstrap_iterations=1000,
        use_cache=False,
    )
    result = run_multi_pair_analysis(
        pairs,
        config,
        args.output_root,
        run_tag=args.run_tag,
        progress=print,
    )
    print(f"SpatialTX {__version__}: {result.run_dir}")
    print(result.pair_results[[
        "pair_label", "regime_pre", "regime_post", "comparability",
        "pre_interface_fraction", "post_interface_fraction",
        "pre_diffuse_fraction", "post_diffuse_fraction",
        "pre_transition_burden", "post_transition_burden",
    ]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
