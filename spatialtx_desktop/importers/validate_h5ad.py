from __future__ import annotations

from pathlib import Path

import numpy as np


def validate_h5ad(path: str | Path, *, require_spatial: bool = True) -> dict:
    """Validate the minimum contract required by the SpatialTX Main Mapper."""
    import anndata as ad

    file_path = Path(path).expanduser().resolve()
    report: dict = {
        "path": str(file_path),
        "valid": False,
        "n_obs": 0,
        "n_vars": 0,
        "has_spatial": False,
        "spatial_shape": None,
        "gene_names_present": False,
        "gene_names_unique": False,
        "errors": [],
        "warnings": [],
    }
    if not file_path.is_file():
        report["errors"].append(f"h5ad file does not exist: {file_path}")
        return report
    if file_path.suffix.lower() != ".h5ad":
        report["errors"].append("Output file must use the .h5ad extension.")
        return report

    adata = None
    try:
        adata = ad.read_h5ad(file_path, backed="r")
        report["n_obs"] = int(adata.n_obs)
        report["n_vars"] = int(adata.n_vars)
        if adata.n_obs <= 0:
            report["errors"].append("AnnData contains no spots/barcodes.")
        if adata.n_vars <= 0:
            report["errors"].append("AnnData contains no genes/features.")

        names = [str(value).strip() for value in adata.var_names]
        report["gene_names_present"] = bool(names) and all(names)
        report["gene_names_unique"] = len(names) == len(set(names))
        if not report["gene_names_present"]:
            report["errors"].append("Gene names are missing or blank.")
        if not report["gene_names_unique"]:
            report["errors"].append("Gene names are not unique.")

        report["has_spatial"] = "spatial" in adata.obsm
        if not report["has_spatial"]:
            message = 'adata.obsm["spatial"] is missing.'
            if require_spatial:
                report["errors"].append(message)
            else:
                report["warnings"].append(
                    message
                    + " Main Mapper can perform expression-only scoring, but spatial QC will be incomplete and spatial maps and metrics will be unavailable."
                )
        else:
            coords = np.asarray(adata.obsm["spatial"], dtype=float)
            report["spatial_shape"] = list(coords.shape)
            if coords.ndim != 2 or coords.shape != (adata.n_obs, 2):
                report["errors"].append("Spatial coordinates must have shape (n_obs, 2).")
            elif not np.isfinite(coords).all():
                report["errors"].append("Spatial coordinates contain non-finite values.")
    except Exception as exc:
        report["errors"].append(f"Unable to read h5ad: {exc}")
    finally:
        if adata is not None and getattr(adata, "file", None) is not None:
            adata.file.close()

    report["valid"] = not report["errors"]
    return report


def require_valid_h5ad(path: str | Path, *, require_spatial: bool = True) -> dict:
    report = validate_h5ad(path, require_spatial=require_spatial)
    if not report["valid"]:
        raise ValueError("Converted h5ad validation failed: " + "; ".join(report["errors"]))
    return report
