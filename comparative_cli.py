from __future__ import annotations

import argparse
import sys
from pathlib import Path

from spatialtx_desktop import __version__
from spatialtx_desktop.comparative.models import COMPARISON_MODES, ComparativeConfig, SampleRecord
from spatialtx_desktop.comparative.runner import run_comparative_analysis
from spatialtx_desktop.comparative.validation import load_comparative_manifest
from spatialtx_desktop.workflow import DEFAULT_C_GENES, DEFAULT_S_GENES, parse_gene_text


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="spatialtx compare",
        description=f"SpatialTX Studio v{__version__} Comparative Spatial Transition Analysis",
    )
    parser.add_argument("--mode", choices=COMPARISON_MODES, default="pairwise")
    parser.add_argument("--manifest", help="Comparative CSV manifest")
    parser.add_argument("--sample-a", help="Reference H5AD for pairwise mode")
    parser.add_argument("--sample-b", help="Target H5AD for pairwise mode")
    parser.add_argument("--reference", default="A", help="Reference group/condition label")
    parser.add_argument("--target", default="B", help="Target group/condition label")
    parser.add_argument("--outdir", "--output", dest="outdir", required=True, help="Output root")
    parser.add_argument("--c-genes", default=",".join(DEFAULT_C_GENES), help="Comma-separated C-side program")
    parser.add_argument("--s-genes", default=",".join(DEFAULT_S_GENES), help="Comma-separated S-side program")
    parser.add_argument("--c-quantile", type=float, default=.80)
    parser.add_argument("--s-quantile", type=float, default=.80)
    parser.add_argument("--g-quantile", type=float, default=.60)
    parser.add_argument("--normalization", choices=("raw_mean", "z_score", "rank_quantile"), default="raw_mean")
    parser.add_argument("--smoothing", choices=("none", "knn_mean", "gaussian"), default="none")
    parser.add_argument("--smoothing-k", type=int, default=6)
    parser.add_argument("--gaussian-sigma", type=float, default=0.0)
    parser.add_argument("--threshold-robustness", action="store_true")
    parser.add_argument("--graph-method", choices=("knn", "radius", "lattice"), default="knn")
    parser.add_argument("--graph-k", type=int, default=6)
    parser.add_argument("--graph-radius", type=float, default=None)
    parser.add_argument("--disable-h-expr", action="store_true")
    parser.add_argument("--disable-v-expr", action="store_true")
    parser.add_argument("--h-genes", default="")
    parser.add_argument("--v-genes", default="")
    parser.add_argument("--context-smoothing", choices=("none", "graph_mean"), default="none")
    parser.add_argument("--context-min-coverage", type=float, default=.25)
    parser.add_argument("--statistical-test", choices=("auto", "wilcoxon", "mannwhitney", "paired_t", "welch_t"), default="auto")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--no-cache", action="store_true")
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def _records(args: argparse.Namespace) -> list[SampleRecord]:
    if args.manifest:
        if args.sample_a or args.sample_b:
            raise ValueError("Use either --manifest or --sample-a/--sample-b, not both.")
        return load_comparative_manifest(args.manifest)
    if args.mode != "pairwise":
        raise ValueError("Paired, unpaired, and manifest_batch modes require --manifest.")
    if not args.sample_a or not args.sample_b:
        raise ValueError("Pairwise mode requires --sample-a and --sample-b.")
    a = Path(args.sample_a).expanduser().resolve()
    b = Path(args.sample_b).expanduser().resolve()
    if a == b:
        raise ValueError("--sample-a and --sample-b must refer to different H5AD files.")
    a_id, b_id = a.stem, b.stem
    if a_id == b_id:
        a_id, b_id = f"{a_id}_A", f"{b_id}_B"
    return [
        SampleRecord(a_id, a, args.reference, condition=args.reference),
        SampleRecord(b_id, b, args.target, condition=args.target),
    ]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        records = _records(args)
        graph_settings = {"method": args.graph_method, "k": args.graph_k, "weighting": "binary"}
        if args.graph_radius is not None:
            graph_settings["radius"] = args.graph_radius
        config = ComparativeConfig(
            mode=args.mode,
            reference=args.reference,
            target=args.target,
            c_genes=parse_gene_text(args.c_genes),
            s_genes=parse_gene_text(args.s_genes),
            c_q=args.c_quantile,
            s_q=args.s_quantile,
            g_q=args.g_quantile,
            scoring_options={
                "normalization_mode": args.normalization,
                "smoothing_mode": args.smoothing,
                "smoothing_k": args.smoothing_k,
                "gaussian_sigma": args.gaussian_sigma,
                "perturbation_check": args.threshold_robustness,
            },
            graph_settings=graph_settings,
            enable_h_expr=not args.disable_h_expr,
            enable_v_expr=not args.disable_v_expr,
            h_genes=parse_gene_text(args.h_genes) or None,
            v_genes=parse_gene_text(args.v_genes) or None,
            context_smoothing=args.context_smoothing,
            context_min_coverage=args.context_min_coverage,
            statistical_test=args.statistical_test,
            seed=args.seed,
            bootstrap_iterations=args.bootstrap_iterations,
            use_cache=not args.no_cache,
        )
        result = run_comparative_analysis(records, config, args.outdir, progress=print)
    except Exception as exc:
        print(f"Comparative Analysis failed: {exc}", file=sys.stderr)
        return 1
    print(f"Comparative Analysis completed: {result.run_dir}")
    print(result.summary_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
