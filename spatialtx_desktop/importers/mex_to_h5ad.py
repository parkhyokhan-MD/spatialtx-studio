from __future__ import annotations

import gzip
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .validate_h5ad import require_valid_h5ad


Progress = Callable[[str], None]
MEX_NAMES = ("matrix.mtx", "matrix.mtx.gz")
FEATURE_NAMES = ("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")
BARCODE_NAMES = ("barcodes.tsv", "barcodes.tsv.gz")


def _first(directory: Path, names: tuple[str, ...]) -> Path | None:
    return next((directory / name for name in names if (directory / name).is_file()), None)


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix.lower() == ".gz" else path.open("r", encoding="utf-8")


def _mex_files(directory: str | Path) -> tuple[Path, Path, Path]:
    folder = Path(directory).expanduser().resolve()
    matrix = _first(folder, MEX_NAMES)
    features = _first(folder, FEATURE_NAMES)
    barcodes = _first(folder, BARCODE_NAMES)
    if not all((matrix, features, barcodes)):
        raise ValueError(
            "A complete 10x MEX/MTX folder requires matrix.mtx(.gz), "
            "barcodes.tsv(.gz), and features.tsv(.gz) or genes.tsv(.gz)."
        )
    return matrix, features, barcodes


def detect_mex_sample(directory: str | Path) -> dict:
    folder = Path(directory).expanduser().resolve()
    report = {
        "valid": False,
        "input_folder": folder,
        "matrix": None,
        "features": None,
        "barcodes": None,
        "errors": [],
        "warnings": [],
    }
    if not folder.is_dir():
        report["errors"].append(f"Input folder does not exist: {folder}")
        return report
    try:
        report["matrix"], report["features"], report["barcodes"] = _mex_files(folder)
        inspection = inspect_mex(folder)
        if inspection["orientation"] == "mismatch":
            report["errors"].append("MEX matrix dimensions do not match feature and barcode tables.")
        report["inspection"] = inspection
    except Exception as exc:
        report["errors"].append(str(exc))
    if not _position_file(folder):
        report["warnings"].append(
            "No Visium tissue-position table was found. H5AD conversion is allowed for expression-only scoring, but Main Mapper will mark spatial QC incomplete and will not generate spatial maps or metrics."
        )
    report["valid"] = not report["errors"]
    return report


def find_mex_folders(root: str | Path) -> list[Path]:
    base = Path(root).expanduser().resolve()
    if not base.is_dir():
        raise ValueError(f"10x MEX/MTX root folder does not exist: {base}")
    folders = [
        directory.resolve()
        for directory in (base, *(path for path in base.rglob("*") if path.is_dir()))
        if _first(directory, MEX_NAMES) and _first(directory, FEATURE_NAMES) and _first(directory, BARCODE_NAMES)
    ]
    return sorted(set(folders), key=lambda path: str(path).lower())


def inspect_mex(directory: str | Path) -> dict:
    from scipy.io import mminfo

    matrix_path, feature_path, barcode_path = _mex_files(directory)
    rows, columns, entries, _format, _field, _symmetry = mminfo(matrix_path)
    with _open_text(feature_path) as handle:
        features = pd.read_csv(handle, sep="\t", header=None)
    with _open_text(barcode_path) as handle:
        barcodes = pd.read_csv(handle, sep="\t", header=None)
    shape = (int(rows), int(columns))
    orientation = (
        "genes_x_barcodes" if shape == (len(features), len(barcodes))
        else "barcodes_x_genes" if shape == (len(barcodes), len(features))
        else "mismatch"
    )
    return {
        "folder": str(Path(directory).resolve()),
        "matrix": matrix_path.name,
        "features": feature_path.name,
        "barcodes": barcode_path.name,
        "matrix_rows": shape[0],
        "matrix_columns": shape[1],
        "feature_rows": len(features),
        "barcode_rows": len(barcodes),
        "orientation": orientation,
        "nonzero_entries": int(entries),
    }


def _unique_names(values: list[str]) -> tuple[list[str], bool]:
    used: set[str] = set()
    counts: dict[str, int] = {}
    result: list[str] = []
    changed = False
    for raw in values:
        name = str(raw).strip() or "unnamed_gene"
        candidate = name
        if candidate in used:
            changed = True
            number = counts.get(name, 0) + 1
            candidate = f"{name}-{number}"
            while candidate in used:
                number += 1
                candidate = f"{name}-{number}"
            counts[name] = number
        else:
            counts.setdefault(name, 0)
        result.append(candidate)
        used.add(candidate)
    return result, changed


def _position_file(folder: Path) -> Path | None:
    candidates = []
    for directory in (folder / "spatial", folder.parent / "spatial", folder):
        for name in (
            "tissue_positions.csv", "tissue_positions.csv.gz",
            "tissue_positions_list.csv", "tissue_positions_list.csv.gz",
        ):
            candidate = directory / name
            if candidate.is_file():
                candidates.append(candidate)
    return candidates[0] if candidates else None


def _attach_positions(adata, folder: Path) -> str | None:
    path = _position_file(folder)
    if path is None:
        return None
    headerless = path.name.startswith("tissue_positions_list.csv")
    positions = pd.read_csv(
        path,
        compression="infer",
        header=None if headerless else "infer",
        names=(
            ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
            if headerless else None
        ),
    )
    if "barcode" not in positions.columns:
        positions = positions.rename(columns={positions.columns[0]: "barcode"})
    positions["barcode"] = positions["barcode"].astype(str)
    positions = positions.drop_duplicates("barcode", keep="first").set_index("barcode")
    aligned = positions.reindex(adata.obs_names)
    for column in positions.columns:
        adata.obs[column] = aligned[column].to_numpy()
    if {"pxl_col_in_fullres", "pxl_row_in_fullres"}.issubset(aligned.columns):
        coords = aligned[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
        if np.isfinite(coords).all():
            adata.obsm["spatial"] = coords
    elif {"array_col", "array_row"}.issubset(aligned.columns):
        coords = aligned[["array_col", "array_row"]].to_numpy(dtype=float)
        if np.isfinite(coords).all():
            adata.obsm["spatial"] = coords
    return path.name


def _safe_sample_name(value: str) -> str:
    name = re.sub(r'[<>:"/\\|?*]+', "_", str(value).strip()).strip(". ")
    if not name:
        raise ValueError("Sample name is empty or invalid.")
    return name


def convert_mex_to_h5ad(
    input_folder: str | Path,
    output_dir: str | Path,
    sample_name: str | None = None,
    *,
    overwrite: bool = False,
    progress: Progress | None = None,
) -> tuple[Path, dict]:
    import anndata as ad
    from scipy import sparse
    from scipy.io import mmread

    detection = detect_mex_sample(input_folder)
    if not detection["valid"]:
        raise ValueError("Invalid 10x MEX/MTX folder: " + "; ".join(detection["errors"]))
    folder: Path = detection["input_folder"]
    name = _safe_sample_name(sample_name or folder.name)
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / f"{name}.h5ad"
    temporary = destination / f".{name}.tmp.h5ad"
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {output}")
    if temporary.exists():
        temporary.unlink()

    try:
        if progress:
            progress("Reading matrix.mtx, barcodes, and features")
        matrix_path, feature_path, barcode_path = _mex_files(folder)
        matrix = sparse.csr_matrix(mmread(matrix_path))
        with _open_text(feature_path) as handle:
            features = pd.read_csv(handle, sep="\t", header=None)
        with _open_text(barcode_path) as handle:
            barcodes = pd.read_csv(handle, sep="\t", header=None)
        if matrix.shape == (len(features), len(barcodes)):
            matrix = matrix.transpose().tocsr()
        elif matrix.shape != (len(barcodes), len(features)):
            raise ValueError(
                f"MEX dimension mismatch: matrix={matrix.shape}, features={len(features)}, barcodes={len(barcodes)}"
            )

        gene_names = features.iloc[:, 1 if features.shape[1] > 1 else 0].fillna("").astype(str).tolist()
        var_names, names_changed = _unique_names(gene_names)
        obs = pd.DataFrame(index=pd.Index(barcodes.iloc[:, 0].astype(str), name="barcode"))
        var = pd.DataFrame(index=pd.Index(var_names, name="gene"))
        var["gene_id"] = features.iloc[:, 0].astype(str).to_numpy()
        if names_changed:
            var["gene_name_original"] = gene_names
        if features.shape[1] > 2:
            var["feature_type"] = features.iloc[:, 2].astype(str).to_numpy()
        adata = ad.AnnData(X=matrix, obs=obs, var=var)
        if progress:
            progress("Attaching available tissue positions")
        position_source = _attach_positions(adata, folder)
        adata.uns["spatialtx_import"] = {
            "schema_version": "1.0",
            "importer": "SpatialTX Raw 10x MEX/MTX importer",
            "source_format": "10x MEX/MTX",
            "sample_name": name,
            "source_folder_name": folder.name,
            "matrix_file": matrix_path.name,
            "features_file": feature_path.name,
            "barcodes_file": barcode_path.name,
            "positions_file": position_source or "not_found",
            "imported_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "gene_names_made_unique": names_changed,
            "n_spots": int(adata.n_obs),
            "n_features": int(adata.n_vars),
        }
        if progress:
            progress("Writing canonical H5AD")
        adata.write_h5ad(temporary, compression="gzip")
        validation = require_valid_h5ad(temporary, require_spatial=False)
        temporary.replace(output)
        validation["path"] = str(output)
        validation["gene_names_made_unique"] = names_changed
        validation["position_source"] = position_source or "not_found"
        if progress:
            progress(f"Validated H5AD: {adata.n_obs} barcodes x {adata.n_vars} features")
        return output, validation
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
