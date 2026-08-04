from __future__ import annotations

import argparse
import sys
from pathlib import Path

from spatialtx_desktop import __version__
from spatialtx_desktop.importers.geo_flat_visium import (
    convert_geo_flat_directory,
    inventory_dataframe,
    scan_geo_flat_directory,
    write_comparative_manifest,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="spatialtx import-geo-flat",
        description=(
            f"SpatialTX Studio v{__version__} generic local GEO-style flat standard Visium importer. "
            "The selected source directory is treated as read-only."
        ),
    )
    parser.add_argument("--input-dir", required=True, help="Any readable local or mounted source directory")
    parser.add_argument("--output-dir", help="Any writable output directory; required for conversion")
    parser.add_argument("--list-samples", action="store_true", help="Scan and print inventory without conversion")
    parser.add_argument("--samples", nargs="+", help="Exact full sample prefixes to convert; default is all valid samples")
    parser.add_argument("--recursive", action="store_true", help="Scan subdirectories; off by default")
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing same-name H5AD only when explicitly set")
    parser.add_argument(
        "--write-comparative-manifest",
        help="Write a draft comparative CSV after conversion; the file must not already exist",
    )
    parser.add_argument(
        "--confirm-filename-pairings",
        action="store_true",
        help=(
            "Explicitly confirm filename-inferred group/pair fields in the written manifest. "
            "Without this flag, group and pair_id stay blank and pairing_source remains unconfirmed."
        ),
    )
    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def _print_inventory(samples) -> None:
    table = inventory_dataframe(samples, absolute_paths=False)
    if table.empty:
        print("No supported GEO-style flat Visium sample prefixes were detected.")
        return
    columns = [
        "sample_prefix", "parsed_subject_id", "parsed_condition",
        "validation_status", "orientation_detected", "warnings", "errors",
    ]
    print(table[columns].to_csv(index=False).rstrip())


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        samples = scan_geo_flat_directory(args.input_dir, recursive=args.recursive)
        if args.list_samples:
            _print_inventory(samples)
            return 0
        if not args.output_dir:
            raise ValueError("--output-dir is required unless --list-samples is used.")
        if args.confirm_filename_pairings and not args.write_comparative_manifest:
            raise ValueError("--confirm-filename-pairings requires --write-comparative-manifest.")
        result = convert_geo_flat_directory(
            args.input_dir,
            args.output_dir,
            selected_samples=args.samples,
            recursive=args.recursive,
            overwrite=args.overwrite,
            progress=print,
        )
        if args.write_comparative_manifest:
            manifest = write_comparative_manifest(
                result["samples"],
                result["converted_paths"],
                args.write_comparative_manifest,
                confirmed=args.confirm_filename_pairings,
            )
            print(f"Comparative manifest written: {manifest}")
        summary = result["summary"]
        successes = int(summary["status"].eq("success").sum()) if len(summary) else 0
        failures = int(len(summary) - successes)
        print(f"GEO flat Visium import completed: {successes} successful, {failures} failed")
        print(f"Reports: {result['report_dir']}")
        return 1 if failures else 0
    except Exception as exc:
        print(f"GEO flat Visium import failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
