from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..gene_program_validation import validate_gene_programs
from .models import ComparativeConfig, SampleRecord


MANIFEST_REQUIRED_COLUMNS = ("sample_id", "file_path", "group")
MANIFEST_OPTIONAL_COLUMNS = ("pair_id", "condition", "batch", "notes")


def load_comparative_manifest(path: str | Path) -> list[SampleRecord]:
    manifest_path = Path(path).expanduser().resolve()
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Comparative manifest not found: {manifest_path}")
    table = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    missing = [column for column in MANIFEST_REQUIRED_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"Comparative manifest missing required columns: {', '.join(missing)}")
    if table.empty:
        raise ValueError("Comparative manifest contains no sample rows.")
    for column in (*MANIFEST_REQUIRED_COLUMNS, *MANIFEST_OPTIONAL_COLUMNS):
        if column not in table:
            table[column] = ""
        table[column] = table[column].astype(str).str.strip()
    empty_required = table[list(MANIFEST_REQUIRED_COLUMNS)].eq("").any(axis=1)
    if empty_required.any():
        rows = ", ".join(str(index + 2) for index in table.index[empty_required])
        raise ValueError(f"Manifest required values are empty on CSV row(s): {rows}")
    duplicate = table["sample_id"].duplicated(keep=False)
    if duplicate.any():
        names = ", ".join(sorted(table.loc[duplicate, "sample_id"].unique()))
        raise ValueError(f"Duplicate sample_id values are not allowed: {names}")
    records: list[SampleRecord] = []
    for _, row in table.iterrows():
        file_path = Path(row["file_path"]).expanduser()
        if not file_path.is_absolute():
            file_path = manifest_path.parent / file_path
        records.append(SampleRecord(
            sample_id=row["sample_id"],
            file_path=file_path.resolve(),
            group=row["group"],
            pair_id=row["pair_id"],
            condition=row["condition"],
            batch=row["batch"],
            notes=row["notes"],
        ))
    return records


def validate_record_structure(records: list[SampleRecord], config: ComparativeConfig) -> None:
    config.validate()
    if len(records) < 2:
        raise ValueError("Comparative analysis requires at least two samples.")
    sample_ids = [record.sample_id for record in records]
    duplicates = sorted({sample for sample in sample_ids if sample_ids.count(sample) > 1})
    if duplicates:
        raise ValueError(f"Duplicate sample_id values are not allowed: {', '.join(duplicates)}")
    validate_gene_programs(config.c_genes, config.s_genes, mode="comparative")
    groups = {record.group for record in records}
    if config.reference not in groups or config.target not in groups:
        raise ValueError(
            f"Reference/target groups must exist in the inputs. Available groups: {', '.join(sorted(groups))}"
        )
    if config.mode == "pairwise":
        reference_count = sum(record.group == config.reference for record in records)
        target_count = sum(record.group == config.target for record in records)
        if len(records) != 2 or reference_count != 1 or target_count != 1:
            raise ValueError("Pairwise mode requires exactly one reference sample and one target sample.")


def preflight_sample(record: SampleRecord, config: ComparativeConfig) -> dict:
    """Read H5AD metadata in backed mode and return a non-silent validation row."""
    row = {
        **record.to_dict(),
        "validation_status": "failed",
        "validation_error": "",
        "n_spots": np.nan,
        "n_genes": np.nan,
        "C_genes_present": 0,
        "S_genes_present": 0,
    }
    path = record.file_path
    if not path.is_file():
        row["validation_error"] = f"File not found: {path}"
        return row
    if path.suffix.lower() != ".h5ad":
        row["validation_error"] = "Canonical comparative input must be an AnnData .h5ad file."
        return row
    try:
        import anndata as ad

        adata = ad.read_h5ad(path, backed="r")
        try:
            row["n_spots"] = int(adata.n_obs)
            row["n_genes"] = int(adata.n_vars)
            if adata.n_obs < 2 or adata.n_vars < 1:
                raise ValueError("H5AD must contain at least two spots and one gene.")
            if "spatial" not in adata.obsm:
                raise ValueError("Valid adata.obsm['spatial'] coordinates are required for comparative analysis.")
            coords = np.asarray(adata.obsm["spatial"], dtype=float)
            if coords.ndim != 2 or coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
                raise ValueError("Spatial coordinates must be finite and shaped (n_obs, 2).")
            names = {str(name).upper() for name in adata.var_names}
            c_present = sum(str(gene).upper() in names for gene in config.c_genes)
            s_present = sum(str(gene).upper() in names for gene in config.s_genes)
            row["C_genes_present"] = c_present
            row["S_genes_present"] = s_present
            if c_present < 1 or s_present < 1:
                raise ValueError(
                    f"Required C/S programs are not represented (C present={c_present}, S present={s_present})."
                )
        finally:
            file_manager = getattr(adata, "file", None)
            if file_manager is not None:
                file_manager.close()
    except Exception as exc:
        row["validation_error"] = str(exc)
        return row
    row["validation_status"] = "ok"
    return row


def preflight_records(records: list[SampleRecord], config: ComparativeConfig) -> pd.DataFrame:
    return pd.DataFrame([preflight_sample(record, config) for record in records])


def valid_records_from_preflight(records: list[SampleRecord], preflight: pd.DataFrame) -> list[SampleRecord]:
    valid_ids = set(preflight.loc[preflight["validation_status"].eq("ok"), "sample_id"].astype(str))
    return [record for record in records if record.sample_id in valid_ids]
