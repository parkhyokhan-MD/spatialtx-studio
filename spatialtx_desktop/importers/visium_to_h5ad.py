from __future__ import annotations

import gzip
import io
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

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


def _decode(values) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
        for value in np.asarray(values)
    ]


def _relative(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.name


def _read_json(path: Path) -> dict:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path.name}.")
    return value


def _read_positions(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, compression="infer")
    if not set(POSITION_COLUMNS).issubset(table.columns):
        raise ValueError(
            f"{path.name} must contain these columns: {', '.join(POSITION_COLUMNS)}"
        )
    table = table[POSITION_COLUMNS].copy()
    table["barcode"] = table["barcode"].astype(str).str.strip()
    if table["barcode"].eq("").any():
        raise ValueError("Position table contains a blank barcode.")
    if table["barcode"].duplicated().any():
        raise ValueError("Position table contains duplicate barcodes.")
    for column in POSITION_COLUMNS[1:]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    if table[POSITION_COLUMNS[1:]].isna().any().any():
        raise ValueError("Position table contains missing or non-numeric position values.")
    return table.set_index("barcode", drop=True)


def _read_png(path: Path) -> np.ndarray:
    from PIL import Image

    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rb") as handle:
            payload = handle.read()
        image = Image.open(io.BytesIO(payload))
    else:
        image = Image.open(path)
    with image:
        return np.asarray(image.convert("RGB"), dtype=np.uint8)


def _matrix_layout_errors(path: Path) -> list[str]:
    import h5py

    errors: list[str] = []
    try:
        with h5py.File(path, "r") as handle:
            if "matrix" not in handle:
                return ['HDF5 file does not contain the required "matrix" group.']
            matrix = handle["matrix"]
            for name in ("data", "indices", "indptr", "shape", "barcodes", "features"):
                if name not in matrix:
                    errors.append(f'HDF5 matrix group is missing "{name}".')
            if "features" in matrix and not any(name in matrix["features"] for name in ("name", "id")):
                errors.append('HDF5 features group must contain "name" or "id".')
            if "shape" in matrix:
                shape = tuple(int(value) for value in np.asarray(matrix["shape"][:]).ravel())
                if len(shape) != 2 or min(shape) <= 0:
                    errors.append("HDF5 count matrix shape must contain two positive dimensions.")
    except Exception as exc:
        errors.append(f"Unable to inspect filtered_feature_bc_matrix.h5: {exc}")
    return errors


def _strip_gzip_suffix(name: str) -> str:
    return name[:-3] if name.lower().endswith(".gz") else name


def _equivalent_file_key(path: Path) -> tuple[str, str]:
    return (str(path.parent).lower(), _strip_gzip_suffix(path.name).lower())


def _choose_equivalent_file(paths: list[Path], label: str, warnings: list[str]) -> Path:
    ordered = sorted(
        paths,
        key=lambda value: (value.suffix.lower() == ".gz", str(value).lower()),
    )
    selected = ordered[0]
    if len(ordered) > 1:
        alternatives = ", ".join(path.name for path in ordered)
        warnings.append(
            f"Duplicate-equivalent {label} files were found ({alternatives}); "
            f"using {selected.name}."
        )
    return selected


def _find_visium_file(
    root: Path,
    relative_names: tuple[str, ...],
    base_names: tuple[str, ...],
    *,
    label: str,
    warnings: list[str],
) -> Path | None:
    """Find an exact Space Ranger name or one GEO-style prefixed equivalent."""
    exact_by_key: dict[tuple[str, str], list[Path]] = {}
    for name in relative_names:
        path = root / name
        if path.is_file():
            exact_by_key.setdefault(_equivalent_file_key(path), []).append(path)
    if exact_by_key:
        for name in relative_names:
            path = root / name
            if path.is_file():
                return _choose_equivalent_file(exact_by_key[_equivalent_file_key(path)], label, warnings)

    candidates: list[Path] = []
    for directory in (root, root / "spatial"):
        if not directory.is_dir():
            continue
        for path in directory.iterdir():
            lower = path.name.lower()
            if path.is_file() and any(lower.endswith("_" + name.lower()) for name in base_names):
                candidates.append(path)

    grouped: dict[tuple[str, str], list[Path]] = {}
    for path in sorted(set(candidates), key=lambda value: str(value).lower()):
        grouped.setdefault(_equivalent_file_key(path), []).append(path)
    selected = [
        _choose_equivalent_file(paths, label, warnings)
        for paths in grouped.values()
    ]
    if len(selected) > 1:
        names = ", ".join(path.name for path in selected)
        raise ValueError(f"Multiple conflicting {label} files were found: {names}")
    return selected[0] if selected else None


def _suggested_sample_name(matrix_path: Path, fallback: str) -> str:
    suffix = "_filtered_feature_bc_matrix.h5"
    return matrix_path.name[:-len(suffix)] if matrix_path.name.lower().endswith(suffix) else fallback


def detect_visium_sample(sample_folder: str | Path) -> dict:
    """Detect the exact files needed for a 10x/Visium to h5ad conversion."""
    selected = Path(sample_folder).expanduser().resolve()
    report: dict = {
        "valid": False,
        "sample_folder": selected,
        "data_root": None,
        "matrix_h5": None,
        "positions": None,
        "scalefactors": None,
        "hires_image": None,
        "lowres_image": None,
        "suggested_sample_name": selected.name,
        "errors": [],
        "warnings": [],
    }
    if not selected.is_dir():
        report["errors"].append(f"Sample folder does not exist: {selected}")
        return report

    data_root = None
    matrix_h5 = None
    try:
        for root in (selected, selected / "outs"):
            matrix_h5 = _find_visium_file(
                root,
                ("filtered_feature_bc_matrix.h5",),
                ("filtered_feature_bc_matrix.h5",),
                label="filtered feature-barcode matrix",
                warnings=report["warnings"],
            )
            if matrix_h5 is not None:
                data_root = root
                break
    except ValueError as exc:
        report["errors"].append(str(exc))
        return report
    if data_root is None:
        report["errors"].append(
            "Missing filtered_feature_bc_matrix.h5 or *_filtered_feature_bc_matrix.h5 "
            "in the selected folder (or its outs folder)."
        )
        return report
    report["data_root"] = data_root
    report["matrix_h5"] = matrix_h5
    report["suggested_sample_name"] = _suggested_sample_name(matrix_h5, selected.name)
    for key, label, relative_names, base_names in (
        (
            "positions",
            "tissue positions",
            (
                "spatial/tissue_positions.csv",
                "spatial/tissue_positions.csv.gz",
                "tissue_positions.csv",
                "tissue_positions.csv.gz",
            ),
            ("tissue_positions.csv", "tissue_positions.csv.gz"),
        ),
        (
            "scalefactors",
            "scalefactors JSON",
            (
                "spatial/scalefactors_json.json",
                "spatial/scalefactors_json.json.gz",
                "scalefactors_json.json",
                "scalefactors_json.json.gz",
            ),
            ("scalefactors_json.json", "scalefactors_json.json.gz"),
        ),
        (
            "hires_image",
            "high-resolution tissue image",
            (
                "spatial/tissue_hires_image.png",
                "spatial/tissue_hires_image.png.gz",
                "tissue_hires_image.png",
                "tissue_hires_image.png.gz",
            ),
            ("tissue_hires_image.png", "tissue_hires_image.png.gz"),
        ),
        (
            "lowres_image",
            "low-resolution tissue image",
            (
                "spatial/tissue_lowres_image.png",
                "spatial/tissue_lowres_image.png.gz",
                "tissue_lowres_image.png",
                "tissue_lowres_image.png.gz",
            ),
            ("tissue_lowres_image.png", "tissue_lowres_image.png.gz"),
        ),
    ):
        try:
            report[key] = _find_visium_file(
                data_root,
                relative_names,
                base_names,
                label=label,
                warnings=report["warnings"],
            )
        except ValueError as exc:
            report["errors"].append(str(exc))

    if report["positions"] is None:
        report["errors"].append("Missing tissue_positions.csv or tissue_positions.csv.gz.")
    if report["scalefactors"] is None:
        report["errors"].append("Missing scalefactors_json.json or scalefactors_json.json.gz.")
    report["errors"].extend(_matrix_layout_errors(report["matrix_h5"]))
    if report["positions"] is not None:
        try:
            positions = _read_positions(report["positions"])
            if positions.empty:
                report["errors"].append("Position table contains no spots.")
        except Exception as exc:
            report["errors"].append(f"Invalid position table: {exc}")
    if report["scalefactors"] is not None:
        try:
            _read_json(report["scalefactors"])
        except Exception as exc:
            report["errors"].append(f"Invalid scalefactors JSON: {exc}")
    if report["hires_image"] is None and report["lowres_image"] is None:
        report["warnings"].append("No optional tissue image was found; conversion can continue without one.")
    report["valid"] = not report["errors"]
    return report


def _read_10x_h5(path: Path):
    import anndata as ad
    import h5py
    from scipy import sparse

    with h5py.File(path, "r") as handle:
        matrix_group = handle["matrix"]
        shape = tuple(int(value) for value in np.asarray(matrix_group["shape"][:]).ravel())
        matrix = sparse.csc_matrix(
            (
                np.asarray(matrix_group["data"][:]),
                np.asarray(matrix_group["indices"][:]),
                np.asarray(matrix_group["indptr"][:]),
            ),
            shape=shape,
        ).transpose().tocsr()
        barcodes = _decode(matrix_group["barcodes"][:])
        features = matrix_group["features"]
        names = _decode(features["name"][:] if "name" in features else features["id"][:])
        ids = _decode(features["id"][:] if "id" in features else names)
        if matrix.shape != (len(barcodes), len(names)):
            raise ValueError(
                f"10x matrix dimensions do not match barcodes/features: matrix={matrix.shape}, "
                f"barcodes={len(barcodes)}, features={len(names)}"
            )
        if not all(name.strip() for name in names):
            raise ValueError("10x HDF5 feature names contain blank values.")
        names_were_unique = len(names) == len(set(names))
        unique_names: list[str] = []
        used: set[str] = set()
        counts: dict[str, int] = {}
        for original in names:
            candidate = original
            if candidate in used:
                number = counts.get(original, 0) + 1
                candidate = f"{original}-{number}"
                while candidate in used:
                    number += 1
                    candidate = f"{original}-{number}"
                counts[original] = number
            else:
                counts.setdefault(original, 0)
            unique_names.append(candidate)
            used.add(candidate)
        obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
        var = pd.DataFrame({"gene_id": ids}, index=pd.Index(unique_names, name="gene"))
        if not names_were_unique:
            var["gene_name_original"] = names
        for source, target in (("feature_type", "feature_type"), ("genome", "genome")):
            if source in features:
                var[target] = _decode(features[source][:])

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    return adata, names_were_unique


def _safe_sample_name(value: str) -> str:
    name = re.sub(r'[<>:"/\\|?*]+', "_", str(value).strip()).strip(". ")
    if not name:
        raise ValueError("Sample name is empty or invalid.")
    return name


def convert_visium_to_h5ad(
    sample_folder: str | Path,
    output_dir: str | Path,
    sample_name: str | None = None,
    *,
    overwrite: bool = False,
    progress: Progress | None = None,
) -> tuple[Path, dict]:
    """Convert one supported 10x/Visium sample into canonical SpatialTX h5ad."""
    import anndata as ad  # noqa: F401 - import here gives a clear dependency error before conversion starts

    report = detect_visium_sample(sample_folder)
    if not report["valid"]:
        raise ValueError("Invalid raw Visium sample: " + "; ".join(report["errors"]))
    selected: Path = report["sample_folder"]
    data_root: Path = report["data_root"]
    name = _safe_sample_name(sample_name or selected.name)
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / f"{name}.h5ad"
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {output}")
    temporary = destination / f".{name}.tmp.h5ad"
    if temporary.exists():
        temporary.unlink()

    try:
        if progress:
            progress("Reading filtered_feature_bc_matrix.h5")
        adata, gene_names_were_unique = _read_10x_h5(report["matrix_h5"])
        if progress:
            progress("Joining Visium spot positions to filtered barcodes")
        positions = _read_positions(report["positions"])
        missing = adata.obs_names.difference(positions.index)
        if len(missing):
            preview = ", ".join(map(str, missing[:5]))
            raise ValueError(f"Position table is missing {len(missing)} filtered barcode(s), including: {preview}")
        aligned = positions.reindex(adata.obs_names)
        for column in POSITION_COLUMNS[1:]:
            values = aligned[column].to_numpy()
            adata.obs[column] = values.astype(np.int64) if np.allclose(values, np.round(values)) else values
        adata.obsm["spatial"] = aligned[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)

        if progress:
            progress("Reading scalefactors and optional tissue images")
        scalefactors = _read_json(report["scalefactors"])
        images: dict[str, np.ndarray] = {}
        if report["hires_image"] is not None:
            images["hires"] = _read_png(report["hires_image"])
        if report["lowres_image"] is not None:
            images["lowres"] = _read_png(report["lowres_image"])
        adata.uns["spatial"] = {
            name: {
                "images": images,
                "scalefactors": scalefactors,
                "metadata": {"source": "10x/Visium", "library_id": name},
            }
        }
        adata.uns["spatialtx_import"] = {
            "schema_version": "1.0",
            "importer": "SpatialTX raw Visium importer",
            "source_format": "10x/Visium filtered feature-barcode HDF5",
            "sample_name": name,
            "source_folder_name": selected.name,
            "matrix_file": _relative(report["matrix_h5"], data_root),
            "positions_file": _relative(report["positions"], data_root),
            "scalefactors_file": _relative(report["scalefactors"], data_root),
            "image_files": [
                _relative(path, data_root)
                for path in (report["hires_image"], report["lowres_image"])
                if path is not None
            ],
            "imported_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "gene_names_were_unique": gene_names_were_unique,
            "gene_names_made_unique": not gene_names_were_unique,
            "n_spots": int(adata.n_obs),
            "n_features": int(adata.n_vars),
        }

        if progress:
            progress("Writing canonical h5ad")
        adata.write_h5ad(temporary, compression="gzip")
        validation = require_valid_h5ad(temporary)
        temporary.replace(output)
        validation["path"] = str(output)
        validation["warnings"] = list(report["warnings"])
        validation["gene_names_made_unique"] = not gene_names_were_unique
        if progress:
            progress(f"Validated h5ad: {adata.n_obs} spots x {adata.n_vars} genes")
        return output, validation
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
