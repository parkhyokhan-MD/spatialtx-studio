from __future__ import annotations

from importlib.resources import files
from pathlib import Path
from typing import Any

import pandas as pd
import anndata as ad
import yaml


def load_h5ad(path: str | Path):
    input_path = Path(path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not input_path.is_file():
        raise ValueError(f"Input path is not a file: {input_path}")
    if input_path.suffix.lower() != ".h5ad":
        raise ValueError(f"Only .h5ad input is supported by the analysis engine: {input_path}")
    return ad.read_h5ad(input_path)


def ensure_output_dirs(output_dir: str | Path) -> dict[str, Path]:
    out = Path(output_dir)
    figures = out / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    return {"output": out, "figures": figures}


def write_yaml(path: str | Path, data: dict[str, Any]) -> None:
    with Path(path).open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True)


def write_metrics(path: str | Path, metrics: dict[str, Any]) -> None:
    pd.DataFrame([metrics]).to_csv(path, index=False)


def write_qc(path: str | Path, qc_rows: list[dict[str, Any]]) -> None:
    pd.DataFrame(qc_rows).to_csv(path, index=False)


def write_selected_genes(path: str | Path, selected_c: list[str], selected_b: list[str]) -> None:
    rows = [{"program": "C", "gene": gene} for gene in selected_c]
    rows.extend({"program": "B", "gene": gene} for gene in selected_b)
    pd.DataFrame(rows).to_csv(path, index=False)


def packaged_default_config():
    """Return the config resource shipped inside source and wheel installs."""
    return files("spatialtx_studio.resources").joinpath("config_default.yaml")


def load_config(default_path: str | Path | None = None, override_path: str | Path | None = None) -> dict[str, Any]:
    resource = Path(default_path) if default_path is not None else packaged_default_config()
    with resource.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}
    if override_path:
        with Path(override_path).open("r", encoding="utf-8") as f:
            override = yaml.safe_load(f) or {}
        config = deep_merge(config, override)
    return config


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged
