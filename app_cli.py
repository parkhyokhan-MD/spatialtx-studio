from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from spatialtx_studio.runner import run_one


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SpatialTX Studio Desktop exploratory research CLI")
    parser.add_argument("--input", default=None, help="Local .h5ad input path")
    parser.add_argument("--manifest", default=None, help="CSV manifest for batch runs. Required columns: sample,input_path")
    parser.add_argument("--output", required=True, help="Local output directory")
    parser.add_argument("--analysis", default="frame26", choices=["frame26", "istz"], help="Analysis engine. frame26 is the default C/S balance-field workflow.")
    parser.add_argument("--gene-mode", default="fixed", choices=["fixed", "semi_auto", "cancer_adaptive", "custom"], help="Gene program mode")
    parser.add_argument("--config", default=None, help="Optional YAML config override")
    parser.add_argument("--sample", default=None, help="Optional sample name override")
    parser.add_argument("--continue-on-error", action="store_true", help="Continue batch manifest runs after a sample failure")
    return parser.parse_args()


def run_manifest(args: argparse.Namespace) -> list[dict]:
    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest file not found: {manifest_path}")

    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    with manifest_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        required = {"sample", "input_path"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest missing required columns: {','.join(sorted(missing))}")

        for row in reader:
            sample = (row.get("sample") or "").strip()
            input_path = (row.get("input_path") or "").strip()
            if not sample or not input_path:
                raise ValueError("Manifest rows must include non-empty sample and input_path values")
            gene_mode = (row.get("gene_mode") or args.gene_mode).strip()
            analysis = (row.get("analysis") or args.analysis).strip()
            sample_output = output_root / sample
            try:
                metrics = run_one(
                    input_path,
                    sample_output,
                    gene_mode=gene_mode,
                    config_path=args.config,
                    sample=sample,
                    analysis=analysis,
                )
                rows.append({"sample": sample, "status": "ok", "output": str(sample_output), **metrics})
            except Exception as exc:
                rows.append({"sample": sample, "status": "failed", "output": str(sample_output), "error": str(exc)})
                if not args.continue_on_error:
                    break

    if not rows:
        raise ValueError(f"Manifest has no sample rows: {manifest_path}")

    summary_path = output_root / "batch_summary.csv"
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with summary_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> int:
    args = parse_args()
    try:
        if args.manifest:
            rows = run_manifest(args)
            failed = sum(1 for row in rows if row.get("status") == "failed")
            print(f"SpatialTX Studio batch completed: {args.output} ({len(rows) - failed} ok, {failed} failed)")
            return 1 if failed else 0
        if not args.input:
            raise ValueError("Either --input or --manifest is required")
        run_one(args.input, args.output, gene_mode=args.gene_mode, config_path=args.config, sample=args.sample, analysis=args.analysis)
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(f"SpatialTX Studio run completed: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
