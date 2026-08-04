from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import tempfile
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from .. import __version__
from .mex_to_h5ad import convert_mex_to_h5ad
from .validate_h5ad import require_valid_h5ad


Progress = Callable[[str], None]
POSITION_COLUMNS = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]
REQUIRED_COMPONENTS = ("matrix", "barcodes", "features", "positions", "scalefactors")
IMAGE_COMPONENTS = ("hires_image", "lowres_image")

# Longest suffix wins. Matching is case-insensitive; source names and paths are preserved.
SUFFIX_COMPONENTS = tuple(sorted((
    ("_tissue_positions_list.csv.gz", "positions"),
    ("_tissue_positions_list.csv", "positions"),
    ("_tissue_positions.csv.gz", "positions"),
    ("_tissue_positions.csv", "positions"),
    ("_scalefactors_json.json.gz", "scalefactors"),
    ("_scalefactors_json.json", "scalefactors"),
    ("_tissue_hires_image.png.gz", "hires_image"),
    ("_tissue_lowres_image.png.gz", "lowres_image"),
    ("_tissue_hires_image.png", "hires_image"),
    ("_tissue_lowres_image.png", "lowres_image"),
    ("_detected_tissue_image.jpeg.gz", "optional_image"),
    ("_detected_tissue_image.jpg.gz", "optional_image"),
    ("_detected_tissue_image.jpeg", "optional_image"),
    ("_detected_tissue_image.jpg", "optional_image"),
    ("_aligned_fiducials.jpeg.gz", "optional_image"),
    ("_aligned_fiducials.jpg.gz", "optional_image"),
    ("_aligned_fiducials.jpeg", "optional_image"),
    ("_aligned_fiducials.jpg", "optional_image"),
    ("_spatial_enrichment.csv.gz", "auxiliary"),
    ("_annotations.csv.gz", "auxiliary"),
    ("_metadata.csv.gz", "auxiliary"),
    ("_spatial_enrichment.csv", "auxiliary"),
    ("_annotations.csv", "auxiliary"),
    ("_metadata.csv", "auxiliary"),
    ("_features.tsv.gz", "features"),
    ("_barcodes.tsv.gz", "barcodes"),
    ("_matrix.mtx.gz", "matrix"),
    ("_features.tsv", "features"),
    ("_barcodes.tsv", "barcodes"),
    ("_matrix.mtx", "matrix"),
    ("_genes.tsv.gz", "features"),
    ("_genes.tsv", "features"),
    ("_image.tiff.gz", "optional_image"),
    ("_image.jpeg.gz", "optional_image"),
    ("_image.tif.gz", "optional_image"),
    ("_image.png.gz", "optional_image"),
    ("_image.jpg.gz", "optional_image"),
    ("_image.tiff", "optional_image"),
    ("_image.jpeg", "optional_image"),
    ("_image.tif", "optional_image"),
    ("_image.png", "optional_image"),
    ("_image.jpg", "optional_image"),
), key=lambda item: len(item[0]), reverse=True))


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _open_binary(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rb")
    return path.open("rb")


def _relative_name(path: Path | None, root: Path) -> str:
    if path is None:
        return ""
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.name


def _file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def _safe_output_stem(value: str) -> str:
    name = re.sub(r'[<>:"/\\|?*]+', "_", str(value).strip()).strip(". ")
    if not name:
        raise ValueError("Sample prefix cannot be converted to a safe output filename.")
    return name


def _component_for_filename(name: str) -> tuple[str, str] | None:
    lower = name.lower()
    for suffix, component in SUFFIX_COMPONENTS:
        if lower.endswith(suffix):
            prefix = name[: -len(suffix)]
            return (prefix, component) if prefix else None
    return None


def _parse_prefix(prefix: str) -> dict[str, str]:
    accession_match = re.search(r"(?i)(?:^|[_-])(GSM\d+)(?:[_-]|$)", prefix)
    accession = accession_match.group(1).upper() if accession_match else ""
    label = prefix
    if accession_match:
        label = (prefix[: accession_match.start(1)] + prefix[accession_match.end(1) :]).strip("_-")
    subject_match = re.search(
        r"(?i)(sample[_-]\d+)[_-](pre|post|before|after|control|treated)$",
        prefix,
    )
    if subject_match:
        subject = subject_match.group(1).replace("-", "_").lower()
        condition = subject_match.group(2).lower()
    else:
        subject = ""
        condition = ""
    return {
        "geo_accession": accession,
        "parsed_subject_id": subject,
        "parsed_condition": condition,
        "parsed_sample_label": label,
    }


@dataclass
class GeoFlatSample:
    sample_prefix: str
    source_directory: Path
    geo_accession: str = ""
    parsed_subject_id: str = ""
    parsed_condition: str = ""
    parsed_sample_label: str = ""
    matrix_file: Path | None = None
    barcodes_file: Path | None = None
    features_file: Path | None = None
    positions_file: Path | None = None
    scalefactors_file: Path | None = None
    hires_image_file: Path | None = None
    lowres_image_file: Path | None = None
    optional_files: list[Path] = field(default_factory=list)
    validation_status: str = "not_validated"
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    matrix_rows: int | None = None
    matrix_columns: int | None = None
    feature_count: int | None = None
    barcode_count: int | None = None
    coordinate_count: int | None = None
    matched_coordinate_count: int | None = None
    missing_coordinate_count: int | None = None
    orientation_detected: str = ""
    orientation_action: str = ""
    nonzero_entries: int | None = None
    components: dict[str, list[Path]] = field(default_factory=dict, repr=False)
    prefix_variants: set[str] = field(default_factory=set, repr=False)

    @property
    def valid(self) -> bool:
        return not self.errors and self.validation_status in {"valid", "warning"}

    @property
    def image_files(self) -> list[Path]:
        return [path for path in (self.hires_image_file, self.lowres_image_file) if path is not None] + [
            path for path in self.optional_files if _component_for_filename(path.name) and _component_for_filename(path.name)[1] == "optional_image"
        ]

    def inventory_record(self, *, absolute_paths: bool = True) -> dict:
        def display(path: Path | None) -> str:
            if path is None:
                return ""
            return str(path) if absolute_paths else _relative_name(path, self.source_directory)

        return {
            "sample_prefix": self.sample_prefix,
            "source_directory": str(self.source_directory) if absolute_paths else self.source_directory.name,
            "geo_accession": self.geo_accession,
            "parsed_subject_id": self.parsed_subject_id,
            "parsed_condition": self.parsed_condition,
            "parsed_sample_label": self.parsed_sample_label,
            "matrix_file": display(self.matrix_file),
            "barcodes_file": display(self.barcodes_file),
            "features_file": display(self.features_file),
            "positions_file": display(self.positions_file),
            "scalefactors_file": display(self.scalefactors_file),
            "hires_image_file": display(self.hires_image_file),
            "lowres_image_file": display(self.lowres_image_file),
            "optional_files": ";".join(display(path) for path in self.optional_files),
            "validation_status": self.validation_status,
            "warnings": "; ".join(self.warnings),
            "errors": "; ".join(self.errors),
            "matrix_rows": self.matrix_rows,
            "matrix_columns": self.matrix_columns,
            "feature_count": self.feature_count,
            "barcode_count": self.barcode_count,
            "coordinate_count": self.coordinate_count,
            "matched_coordinate_count": self.matched_coordinate_count,
            "missing_coordinate_count": self.missing_coordinate_count,
            "orientation_detected": self.orientation_detected,
            "orientation_action": self.orientation_action,
            "image_availability": bool(self.image_files),
        }


def _assign_component(sample: GeoFlatSample, component: str) -> None:
    paths = sample.components.get(component, [])
    if component in REQUIRED_COMPONENTS or component in IMAGE_COMPONENTS:
        if len(paths) > 1:
            conflicts = ", ".join(str(path) for path in paths)
            sample.errors.append(f"Duplicate {component} candidates: {conflicts}")
            return
        value = paths[0] if paths else None
        attribute = {
            "matrix": "matrix_file",
            "barcodes": "barcodes_file",
            "features": "features_file",
            "positions": "positions_file",
            "scalefactors": "scalefactors_file",
            "hires_image": "hires_image_file",
            "lowres_image": "lowres_image_file",
        }[component]
        setattr(sample, attribute, value)
    else:
        sample.optional_files.extend(paths)


def _check_readable(path: Path, sample: GeoFlatSample) -> None:
    try:
        with path.open("rb") as handle:
            handle.read(1)
    except OSError as exc:
        sample.errors.append(f"Unreadable file {path}: {exc}")


def _read_component_table(path: Path, label: str) -> pd.DataFrame:
    try:
        with _open_text(path) as handle:
            table = pd.read_csv(handle, sep="\t", header=None, dtype=str)
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f"Empty {label} file: {path.name}") from exc
    except (OSError, UnicodeError, gzip.BadGzipFile, pd.errors.ParserError, EOFError) as exc:
        raise ValueError(f"Malformed {label} file {path.name}: {exc}") from exc
    if table.empty or table.shape[1] < 1:
        raise ValueError(f"Empty {label} file: {path.name}")
    first = table.iloc[:, 0].fillna("").astype(str).str.strip()
    if first.eq("").any():
        raise ValueError(f"{label.capitalize()} file contains blank identifiers: {path.name}")
    return table


def read_geo_flat_positions(path: Path) -> pd.DataFrame:
    try:
        raw = pd.read_csv(path, compression="infer", header=None, dtype=str)
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f"Tissue-position file contains no rows: {path.name}") from exc
    except (OSError, UnicodeError, gzip.BadGzipFile, pd.errors.ParserError, EOFError) as exc:
        raise ValueError(f"Malformed tissue-position file {path.name}: {exc}") from exc
    if raw.empty or raw.shape[1] < len(POSITION_COLUMNS):
        raise ValueError(f"Tissue-position file must contain at least six columns: {path.name}")
    first_row = [str(value).strip().lower() for value in raw.iloc[0, :].tolist()]
    if "barcode" in first_row:
        table = pd.read_csv(path, compression="infer", header=0, dtype=str)
        if not set(POSITION_COLUMNS).issubset(table.columns):
            raise ValueError(f"Tissue-position header is missing required columns: {path.name}")
        table = table[POSITION_COLUMNS].copy()
    else:
        table = raw.iloc[:, : len(POSITION_COLUMNS)].copy()
        table.columns = POSITION_COLUMNS
    table["barcode"] = table["barcode"].fillna("").astype(str).str.strip()
    if table["barcode"].eq("").any():
        raise ValueError(f"Tissue-position file contains blank barcodes: {path.name}")
    duplicates = table.loc[table["barcode"].duplicated(keep=False), "barcode"].unique()
    if len(duplicates):
        preview = ", ".join(map(str, duplicates[:5]))
        raise ValueError(f"Tissue-position file contains duplicate barcodes ({preview}): {path.name}")
    for column in POSITION_COLUMNS[1:]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    if table[POSITION_COLUMNS[1:]].isna().any().any():
        raise ValueError(f"Tissue-position file contains non-numeric or missing coordinates: {path.name}")
    return table.set_index("barcode", drop=True)


def _read_scalefactors(path: Path) -> dict:
    try:
        with _open_text(path) as handle:
            value = json.load(handle)
    except (OSError, UnicodeError, gzip.BadGzipFile, json.JSONDecodeError, EOFError) as exc:
        raise ValueError(f"Malformed scalefactor JSON {path.name}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Scalefactor JSON must contain an object: {path.name}")
    return value


def validate_geo_flat_sample(sample: GeoFlatSample) -> GeoFlatSample:
    sample.warnings = [warning for warning in sample.warnings if not warning.startswith("Validation: ")]
    sample.errors = [error for error in sample.errors if error.startswith("Duplicate ") or error.startswith("Ambiguous ")]
    sample.optional_files = []
    for component in (*REQUIRED_COMPONENTS, *IMAGE_COMPONENTS, "optional_image", "auxiliary"):
        _assign_component(sample, component)
    for component, attribute in (
        ("matrix", "matrix_file"),
        ("barcodes", "barcodes_file"),
        ("features or genes", "features_file"),
        ("tissue positions", "positions_file"),
        ("scalefactors JSON", "scalefactors_file"),
    ):
        if getattr(sample, attribute) is None and not sample.components.get(attribute.replace("_file", "")):
            sample.errors.append(f"Missing required {component} component for {sample.sample_prefix}.")
    if len(sample.prefix_variants) > 1:
        sample.errors.append(
            "Ambiguous duplicate sample prefix differing only by case: " + ", ".join(sorted(sample.prefix_variants))
        )
    for paths in sample.components.values():
        for path in paths:
            _check_readable(path, sample)

    features = barcodes = positions = None
    if sample.features_file is not None and len(sample.components.get("features", [])) == 1:
        try:
            features = _read_component_table(sample.features_file, "feature")
            sample.feature_count = len(features)
        except ValueError as exc:
            sample.errors.append(str(exc))
    if sample.barcodes_file is not None and len(sample.components.get("barcodes", [])) == 1:
        try:
            barcodes = _read_component_table(sample.barcodes_file, "barcode")
            sample.barcode_count = len(barcodes)
            if barcodes.iloc[:, 0].duplicated().any():
                sample.errors.append(f"Barcode file contains duplicate identifiers: {sample.barcodes_file.name}")
        except ValueError as exc:
            sample.errors.append(str(exc))
    if sample.positions_file is not None and len(sample.components.get("positions", [])) == 1:
        try:
            positions = read_geo_flat_positions(sample.positions_file)
            sample.coordinate_count = len(positions)
        except ValueError as exc:
            sample.errors.append(str(exc))
    if sample.scalefactors_file is not None and len(sample.components.get("scalefactors", [])) == 1:
        try:
            scalefactors = _read_scalefactors(sample.scalefactors_file)
            for key in ("tissue_hires_scalef", "tissue_lowres_scalef", "spot_diameter_fullres"):
                if key not in scalefactors:
                    sample.warnings.append(f"Validation: optional scalefactor key is missing: {key}")
        except ValueError as exc:
            sample.errors.append(str(exc))

    if sample.matrix_file is not None and len(sample.components.get("matrix", [])) == 1:
        try:
            from scipy.io import mminfo

            # A file handle avoids filename-encoding assumptions in SciPy on Windows.
            with _open_binary(sample.matrix_file) as matrix_handle:
                rows, columns, entries, matrix_format, _field, _symmetry = mminfo(matrix_handle)
            sample.matrix_rows = int(rows)
            sample.matrix_columns = int(columns)
            sample.nonzero_entries = int(entries)
            if matrix_format != "coordinate":
                sample.errors.append(f"Unsupported Matrix Market format {matrix_format}; coordinate format is required.")
            if min(int(rows), int(columns)) <= 0 or int(entries) <= 0:
                sample.errors.append("Matrix Market input is empty.")
        except Exception as exc:
            sample.errors.append(f"Malformed Matrix Market file {sample.matrix_file.name}: {exc}")

    if None not in (sample.matrix_rows, sample.matrix_columns, sample.feature_count, sample.barcode_count):
        shape = (int(sample.matrix_rows), int(sample.matrix_columns))
        feature_by_barcode = (int(sample.feature_count), int(sample.barcode_count))
        barcode_by_feature = (int(sample.barcode_count), int(sample.feature_count))
        if feature_by_barcode == barcode_by_feature and shape == feature_by_barcode:
            sample.orientation_detected = "ambiguous_equal_dimensions"
            sample.orientation_action = "fail"
            sample.errors.append("Ambiguous matrix orientation because feature and barcode counts are equal.")
        elif shape == feature_by_barcode:
            sample.orientation_detected = "genes_x_barcodes"
            sample.orientation_action = "transpose_to_barcodes_x_genes"
        elif shape == barcode_by_feature:
            sample.orientation_detected = "barcodes_x_genes"
            sample.orientation_action = "none"
        else:
            sample.orientation_detected = "dimension_mismatch"
            sample.orientation_action = "fail"
            if int(sample.matrix_rows) not in {int(sample.feature_count), int(sample.barcode_count)}:
                sample.errors.append(
                    f"Matrix row count {sample.matrix_rows} matches neither feature count {sample.feature_count} nor barcode count {sample.barcode_count}."
                )
            if int(sample.matrix_columns) not in {int(sample.feature_count), int(sample.barcode_count)}:
                sample.errors.append(
                    f"Matrix column count {sample.matrix_columns} matches neither feature count {sample.feature_count} nor barcode count {sample.barcode_count}."
                )
            sample.errors.append(
                f"Matrix dimensions {shape} do not match features={sample.feature_count}, barcodes={sample.barcode_count}."
            )

    if barcodes is not None and positions is not None:
        barcode_index = pd.Index(barcodes.iloc[:, 0].astype(str).str.strip())
        missing = barcode_index.difference(positions.index)
        extra = positions.index.difference(barcode_index)
        sample.matched_coordinate_count = int(len(barcode_index) - len(missing))
        sample.missing_coordinate_count = int(len(missing))
        if len(missing):
            preview = ", ".join(map(str, missing[:5]))
            sample.errors.append(f"Position table is missing {len(missing)} matrix barcode(s), including: {preview}")
        if len(extra):
            sample.warnings.append(f"Validation: position table contains {len(extra)} extra barcode(s) not present in the matrix.")

    if sample.hires_image_file is None and sample.lowres_image_file is None:
        sample.warnings.append(
            "Validation: no optional hires or lowres tissue image was found; expression conversion is allowed but image overlay is unavailable."
        )
    if not sample.parsed_subject_id and not sample.parsed_condition:
        sample.warnings.append(
            "Validation: filename subject/condition metadata was not confidently inferred; grouping remains based only on the full prefix."
        )
    sample.errors = list(dict.fromkeys(sample.errors))
    sample.warnings = list(dict.fromkeys(sample.warnings))
    sample.validation_status = "invalid" if sample.errors else "warning" if sample.warnings else "valid"
    return sample


def scan_geo_flat_directory(input_dir: str | Path, *, recursive: bool = False) -> list[GeoFlatSample]:
    source = Path(input_dir).expanduser().resolve()
    if not source.is_dir():
        raise ValueError(f"Selected GEO flat source directory does not exist: {source}")
    if not os.access(source, os.R_OK):
        raise PermissionError(f"Selected GEO flat source directory is not readable: {source}")
    try:
        paths = list(source.rglob("*")) if recursive else list(source.iterdir())
    except OSError as exc:
        raise OSError(f"Unable to scan selected source directory {source}: {exc}") from exc
    files = sorted((path for path in paths if path.is_file()), key=lambda path: str(path).casefold())
    grouped: dict[str, GeoFlatSample] = {}
    for path in files:
        match = _component_for_filename(path.name)
        if match is None:
            continue
        prefix, component = match
        key = prefix.casefold()
        if key not in grouped:
            metadata = _parse_prefix(prefix)
            grouped[key] = GeoFlatSample(
                sample_prefix=prefix,
                source_directory=source,
                **metadata,
            )
        sample = grouped[key]
        sample.prefix_variants.add(prefix)
        sample.components.setdefault(component, []).append(path.resolve())
    samples = [validate_geo_flat_sample(sample) for sample in grouped.values()]
    return sorted(samples, key=lambda sample: (sample.sample_prefix.casefold(), sample.sample_prefix))


def validate_output_directory(output_dir: str | Path) -> Path:
    destination = Path(output_dir).expanduser().resolve()
    try:
        destination.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(prefix=".spatialtx_write_test_", dir=destination, delete=True):
            pass
    except OSError as exc:
        raise PermissionError(f"Selected output directory is not writable: {destination}: {exc}") from exc
    return destination


def _reject_source_side_output(source: Path, destination: Path) -> None:
    if destination == source or source in destination.parents:
        raise ValueError(
            "The output directory must be outside the selected read-only source directory."
        )


def _link_or_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, target)
    except OSError:
        shutil.copy2(source, target)


def _canonical_name(path: Path, component: str) -> str:
    compressed = path.name.lower().endswith(".gz")
    if component == "matrix":
        return "matrix.mtx.gz" if compressed else "matrix.mtx"
    if component == "barcodes":
        return "barcodes.tsv.gz" if compressed else "barcodes.tsv"
    if component == "features":
        return "features.tsv.gz" if compressed else "features.tsv"
    if component == "positions":
        old_layout = "tissue_positions_list" in path.name.lower()
        base = "tissue_positions_list.csv" if old_layout else "tissue_positions.csv"
        return base + (".gz" if compressed else "")
    raise ValueError(component)


def _unique_asset_directory(destination: Path, stem: str) -> tuple[Path, str]:
    root = destination / "spatialtx_import_assets"
    candidate = root / stem
    number = 2
    while candidate.exists():
        candidate = root / f"{stem}_{number}"
        number += 1
    return candidate, candidate.relative_to(destination).as_posix()


def convert_geo_flat_sample(
    sample: GeoFlatSample,
    output_dir: str | Path,
    *,
    overwrite: bool = False,
    progress: Progress | None = None,
) -> tuple[Path, dict]:
    import anndata as ad

    validate_geo_flat_sample(sample)
    if not sample.valid:
        raise ValueError(f"Invalid GEO flat Visium sample {sample.sample_prefix}: " + "; ".join(sample.errors))
    destination_candidate = Path(output_dir).expanduser().resolve()
    _reject_source_side_output(sample.source_directory.resolve(), destination_candidate)
    destination = validate_output_directory(destination_candidate)
    stem = _safe_output_stem(sample.sample_prefix)
    output = destination / f"{stem}.h5ad"
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output already exists; no file was overwritten: {output}")
    required_paths = [
        sample.matrix_file,
        sample.barcodes_file,
        sample.features_file,
        sample.positions_file,
        sample.scalefactors_file,
    ]
    assert all(path is not None for path in required_paths)
    hashes = {path.name: _file_sha256(path) for path in required_paths if path is not None}
    started = time.perf_counter()
    created_asset_dir: Path | None = None
    with tempfile.TemporaryDirectory(prefix="spatialtx_geo_flat_") as temporary_name:
        temporary_root = Path(temporary_name)
        canonical = temporary_root / "canonical" / stem
        canonical.mkdir(parents=True)
        for component, source in (
            ("matrix", sample.matrix_file),
            ("barcodes", sample.barcodes_file),
            ("features", sample.features_file),
            ("positions", sample.positions_file),
        ):
            assert source is not None
            _link_or_copy(source, canonical / _canonical_name(source, component))
        if progress:
            progress(f"{sample.sample_prefix}: converting through the existing MEX/MTX engine")
        staged_output_dir = temporary_root / "converted"
        staged_h5ad, _base_validation = convert_mex_to_h5ad(
            canonical,
            staged_output_dir,
            stem,
            overwrite=False,
            progress=progress,
        )
        adata = ad.read_h5ad(staged_h5ad)
        if "spatial" not in adata.obsm:
            raise ValueError("Existing converter did not produce valid spatial coordinates for this flattened sample.")
        scalefactors = _read_scalefactors(sample.scalefactors_file)
        asset_relative_paths: list[str] = []
        staged_assets = temporary_root / "assets"
        for image in sample.image_files:
            staged_asset = staged_assets / image.name
            _link_or_copy(image, staged_asset)
            asset_relative_paths.append(image.name)
        asset_target = None
        asset_relative_root = ""
        if asset_relative_paths:
            asset_target, asset_relative_root = _unique_asset_directory(destination, stem)
        output_image_refs = [f"{asset_relative_root}/{name}" for name in asset_relative_paths]
        relative_originals = {
            "matrix": _relative_name(sample.matrix_file, sample.source_directory),
            "barcodes": _relative_name(sample.barcodes_file, sample.source_directory),
            "features": _relative_name(sample.features_file, sample.source_directory),
            "positions": _relative_name(sample.positions_file, sample.source_directory),
            "scalefactors": _relative_name(sample.scalefactors_file, sample.source_directory),
            "hires_image": _relative_name(sample.hires_image_file, sample.source_directory),
            "lowres_image": _relative_name(sample.lowres_image_file, sample.source_directory),
            "optional": [_relative_name(path, sample.source_directory) for path in sample.optional_files],
        }
        adata.uns["spatial"] = {
            stem: {
                "images": {},
                "scalefactors": scalefactors,
                "metadata": {
                    "source": "GEO-style flat standard Visium directory",
                    "library_id": stem,
                    "image_files_relative_to_h5ad_directory": output_image_refs,
                    "image_overlay_available": bool(output_image_refs),
                },
            }
        }
        adata.uns["spatialtx_import"] = {
            "schema_version": "1.1",
            "source_mode": "geo_flat_visium",
            "source_directory": sample.source_directory.name,
            "sample_prefix": sample.sample_prefix,
            "geo_accession": sample.geo_accession,
            "parsed_subject_id": sample.parsed_subject_id,
            "parsed_condition": sample.parsed_condition,
            "parsed_metadata_status": "filename_inference_unconfirmed",
            "original_files": relative_originals,
            "spatialtx_version": f"v{__version__}",
            "import_timestamp": _now(),
            "matrix_shape": [int(sample.matrix_rows), int(sample.matrix_columns)],
            "barcode_count": int(sample.barcode_count),
            "feature_count": int(sample.feature_count),
            "coordinate_count": int(sample.coordinate_count),
            "matched_coordinate_count": int(sample.matched_coordinate_count or 0),
            "orientation_detected": sample.orientation_detected,
            "orientation_action": sample.orientation_action,
            "input_sha256_by_filename": hashes,
            "warnings": sample.warnings,
            "source_files_modified": False,
        }
        augmented = temporary_root / f"{stem}.augmented.h5ad"
        adata.write_h5ad(augmented, compression="gzip")
        validation = require_valid_h5ad(augmented, require_spatial=True)
        if asset_target is not None:
            asset_target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(staged_assets), str(asset_target))
            created_asset_dir = asset_target
        try:
            if output.exists() and overwrite:
                output.unlink()
            shutil.move(str(augmented), str(output))
        except Exception:
            if created_asset_dir is not None and created_asset_dir.exists():
                shutil.rmtree(created_asset_dir)
            raise
    validation.update({
        "path": str(output),
        "sample_prefix": sample.sample_prefix,
        "warnings": sample.warnings,
        "orientation_detected": sample.orientation_detected,
        "orientation_action": sample.orientation_action,
        "input_hashes": hashes,
        "duration_seconds": round(time.perf_counter() - started, 4),
        "asset_directory": str(created_asset_dir) if created_asset_dir is not None else "",
    })
    if progress:
        progress(f"{sample.sample_prefix}: validated {validation['n_obs']} spots x {validation['n_vars']} genes")
    return output, validation


def inventory_dataframe(samples: Iterable[GeoFlatSample], *, absolute_paths: bool = True) -> pd.DataFrame:
    return pd.DataFrame([sample.inventory_record(absolute_paths=absolute_paths) for sample in samples])


def write_comparative_manifest(
    samples: Iterable[GeoFlatSample],
    converted_paths: dict[str, Path],
    path: str | Path,
    *,
    confirmed: bool = False,
    edits: dict[str, dict[str, str]] | None = None,
) -> Path:
    edits = edits or {}
    rows = []
    for sample in sorted(samples, key=lambda item: item.sample_prefix.casefold()):
        if sample.sample_prefix not in converted_paths:
            continue
        override = edits.get(sample.sample_prefix, {})
        if confirmed:
            group = override.get("group", sample.parsed_condition).strip()
            pair_id = override.get("pair_id", sample.parsed_subject_id).strip()
            condition = override.get("condition", sample.parsed_condition).strip()
            if not group:
                raise ValueError(f"Confirmed manifest requires a group for {sample.sample_prefix}.")
            pairing_source = "user_confirmed"
        else:
            group = ""
            pair_id = ""
            condition = sample.parsed_condition
            pairing_source = "filename_inference_unconfirmed"
        rows.append({
            "sample_id": sample.sample_prefix,
            "file_path": str(converted_paths[sample.sample_prefix].resolve()),
            "group": group,
            "pair_id": pair_id,
            "condition": condition,
            "batch": override.get("batch", "").strip(),
            "notes": override.get("notes", "Filename-derived metadata must be reviewed before comparative analysis.").strip(),
            "pairing_source": pairing_source,
        })
    if not rows:
        raise ValueError("No successfully converted samples are available for a comparative manifest.")
    manifest_path = Path(path).expanduser().resolve()
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        raise FileExistsError(f"Comparative manifest already exists; no file was overwritten: {manifest_path}")
    pd.DataFrame(rows).to_csv(manifest_path, index=False)
    return manifest_path


def convert_geo_flat_directory(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    selected_samples: Iterable[str] | None = None,
    recursive: bool = False,
    overwrite: bool = False,
    progress: Progress | None = None,
) -> dict:
    source = Path(input_dir).expanduser().resolve()
    destination_candidate = Path(output_dir).expanduser().resolve()
    _reject_source_side_output(source, destination_candidate)
    destination = validate_output_directory(destination_candidate)
    samples = scan_geo_flat_directory(source, recursive=recursive)
    requested = set(selected_samples) if selected_samples is not None else {
        sample.sample_prefix for sample in samples if sample.valid
    }
    known = {sample.sample_prefix for sample in samples}
    unknown = sorted(requested.difference(known))
    if unknown:
        raise ValueError("Requested sample prefix was not detected: " + ", ".join(unknown))
    report_dir = destination / "geo_flat_import_logs" / datetime.now().strftime("geo_flat_import_%Y%m%d_%H%M%S_%f")
    report_dir.mkdir(parents=True, exist_ok=False)
    inventory = inventory_dataframe(samples)
    inventory.to_csv(report_dir / "geo_flat_inventory.csv", index=False)
    inventory.to_csv(report_dir / "geo_flat_validation_report.csv", index=False)
    conversions: list[dict] = []
    converted_paths: dict[str, Path] = {}
    for sample in samples:
        if sample.sample_prefix not in requested:
            continue
        row = sample.inventory_record()
        started = time.perf_counter()
        row.update({"output_h5ad": "", "status": "pending", "error_message": "", "duration_seconds": 0.0})
        if not sample.valid:
            row["status"] = "failed_validation"
            row["error_message"] = "; ".join(sample.errors)
            conversions.append(row)
            continue
        try:
            if progress:
                progress(f"Converting {sample.sample_prefix}")
            output, validation = convert_geo_flat_sample(
                sample,
                destination,
                overwrite=overwrite,
                progress=progress,
            )
            converted_paths[sample.sample_prefix] = output
            row.update({
                "output_h5ad": str(output),
                "status": "success",
                "error_message": "",
                "spot_count": validation["n_obs"],
                "gene_count": validation["n_vars"],
                "input_hashes": json.dumps(validation["input_hashes"], sort_keys=True),
            })
        except Exception as exc:
            row["status"] = "failed_conversion"
            row["error_message"] = str(exc)
        row["duration_seconds"] = round(time.perf_counter() - started, 4)
        conversions.append(row)
    summary = pd.DataFrame(conversions)
    summary.to_csv(report_dir / "geo_flat_conversion_summary.csv", index=False)
    log = {
        "spatialtx_version": __version__,
        "source_mode": "geo_flat_visium",
        "source_directory": str(source),
        "output_directory": str(destination),
        "report_directory": str(report_dir),
        "recursive": bool(recursive),
        "requested_samples": sorted(requested),
        "detected_sample_count": len(samples),
        "successful_samples": sorted(converted_paths),
        "failed_samples": summary.loc[~summary.get("status", pd.Series(dtype=str)).eq("success")].to_dict("records") if len(summary) else [],
        "conversion_records": summary.to_dict("records") if len(summary) else [],
        "timestamp": _now(),
        "source_files_modified": False,
    }
    (report_dir / "geo_flat_import_log.json").write_text(
        json.dumps(log, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    return {
        "samples": samples,
        "inventory": inventory,
        "summary": summary,
        "converted_paths": converted_paths,
        "report_dir": report_dir,
        "log": log,
    }
