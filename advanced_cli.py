from __future__ import annotations

import argparse
import sys
from pathlib import Path

from spatialtx_desktop.advanced_analysis import MODULE_LABELS, run_advanced_batch
from spatialtx_desktop.workflow import DEFAULT_C_GENES, DEFAULT_S_GENES, parse_gene_text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SpatialTX Studio Desktop v0.4-beta Advanced Analysis CLI")
    parser.add_argument("--module", required=True, choices=sorted([*MODULE_LABELS, "spatial_graph"]))
    parser.add_argument("--input", required=True, nargs="+", help="One or more local .h5ad files")
    parser.add_argument("--output", required=True, help="Output root")
    parser.add_argument("--c-genes", default=",".join(DEFAULT_C_GENES), help="Comma-separated Cx genes")
    parser.add_argument("--s-genes", default=",".join(DEFAULT_S_GENES), help="Comma-separated Sx genes")
    parser.add_argument("--c-quantile", type=float, default=.80)
    parser.add_argument("--s-quantile", type=float, default=.80)
    parser.add_argument("--g-quantile", type=float, default=.60)
    parser.add_argument("--permutations", type=int, default=499)
    parser.add_argument("--seed", type=int, default=20260705)
    parser.add_argument("--graph-method", choices=("radius", "lattice", "knn"), default="radius")
    parser.add_argument("--radius", type=float, default=None, help="Radius in the active native or calibrated coordinate unit")
    parser.add_argument("--k", type=int, default=6, help="K for symmetric KNN graph")
    parser.add_argument("--weighting", choices=("binary", "inverse_distance", "gaussian"), default="binary")
    parser.add_argument("--symmetrization", choices=("union", "mutual"), default="union")
    parser.add_argument("--disable-h-expr", action="store_true", help="Do not calculate H_expr context field")
    parser.add_argument("--disable-v-expr", action="store_true", help="Do not calculate V_expr context field")
    parser.add_argument("--h-genes", default="", help="Optional comma-separated H_expr genes; defaults are used when omitted")
    parser.add_argument("--v-genes", default="", help="Optional comma-separated V_expr genes; defaults are used when omitted")
    parser.add_argument("--context-smoothing", choices=("none", "graph_mean"), default="none")
    parser.add_argument("--h-quantile", type=float, default=.80)
    parser.add_argument("--v-quantile", type=float, default=.80)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        paths = [Path(value).expanduser().resolve() for value in args.input]
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            raise FileNotFoundError("Input file(s) not found: " + ", ".join(missing))
        c_genes = parse_gene_text(args.c_genes)
        s_genes = parse_gene_text(args.s_genes)
        if args.module == "spatial_graph":
            from spatialtx_desktop.graph.builder import GraphBuildConfig
            from spatialtx_desktop.graph.runner import SpatialGraphAnalysisConfig, run_spatial_graph_neighborhood_batch

            config = SpatialGraphAnalysisConfig(
                graph=GraphBuildConfig(
                    method=args.graph_method,
                    radius=args.radius,
                    k=args.k,
                    weighting=args.weighting,
                    symmetrization=args.symmetrization,
                    random_seed=args.seed,
                ),
                enable_h=not args.disable_h_expr,
                enable_v=not args.disable_v_expr,
                h_genes=parse_gene_text(args.h_genes) or None,
                v_genes=parse_gene_text(args.v_genes) or None,
                context_smoothing=args.context_smoothing,
                h_high_quantile=args.h_quantile,
                v_high_quantile=args.v_quantile,
                permutations=args.permutations,
                seed=args.seed,
                c_q=args.c_quantile,
                s_q=args.s_quantile,
                g_q=args.g_quantile,
            )
            run_dir, manifest = run_spatial_graph_neighborhood_batch(paths, args.output, c_genes, s_genes, config, progress=print)
        else:
            run_dir, manifest = run_advanced_batch(
                args.module,
                paths,
                args.output,
                c_genes,
                s_genes,
                progress=print,
                c_q=args.c_quantile,
                s_q=args.s_quantile,
                g_q=args.g_quantile,
                permutations=args.permutations,
                seed=args.seed,
            )
        failures = int((manifest["status"] != "ok").sum())
        print(f"Advanced Analysis completed: {run_dir}")
        return 1 if failures else 0
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
