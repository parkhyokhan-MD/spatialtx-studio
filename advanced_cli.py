from __future__ import annotations

import argparse
import sys
from pathlib import Path

from spatialtx_desktop.advanced_analysis import MODULE_LABELS, run_advanced_batch
from spatialtx_desktop.workflow import DEFAULT_C_GENES, DEFAULT_S_GENES, parse_gene_text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SpatialTX Studio v0.2-beta Advanced Analysis CLI")
    parser.add_argument("--module", required=True, choices=sorted(MODULE_LABELS))
    parser.add_argument("--input", required=True, nargs="+", help="One or more local .h5ad files")
    parser.add_argument("--output", required=True, help="Output root")
    parser.add_argument("--c-genes", default=",".join(DEFAULT_C_GENES), help="Comma-separated Cx genes")
    parser.add_argument("--s-genes", default=",".join(DEFAULT_S_GENES), help="Comma-separated Sx genes")
    parser.add_argument("--c-quantile", type=float, default=.80)
    parser.add_argument("--s-quantile", type=float, default=.80)
    parser.add_argument("--g-quantile", type=float, default=.60)
    parser.add_argument("--permutations", type=int, default=499)
    parser.add_argument("--seed", type=int, default=20260705)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        paths = [Path(value).expanduser().resolve() for value in args.input]
        missing = [str(path) for path in paths if not path.is_file()]
        if missing:
            raise FileNotFoundError("Input file(s) not found: " + ", ".join(missing))
        run_dir, manifest = run_advanced_batch(
            args.module,
            paths,
            args.output,
            parse_gene_text(args.c_genes),
            parse_gene_text(args.s_genes),
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
