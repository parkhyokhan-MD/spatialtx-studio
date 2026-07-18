from __future__ import annotations

import datetime as dt
import json
import math
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from . import __version__
from .gene_program_validation import GeneProgramValidationResult, validate_gene_programs
from .workflow import (
    SPATIAL_QC_MESSAGE,
    _dense,
    _gene_indices,
    _is_count_like,
    _knn,
    _read_h5ad,
    score_h5ad,
)


Progress = Callable[[str], None]
MODULE_LABELS = {
    "composition": "Gene Composition",
    "enrichment": "Interface Enrichment",
    "interaction": "Cx/Sx Interaction",
}


def _prepare(path: str | Path, c_genes: Iterable[str], s_genes: Iterable[str]):
    adata = _read_h5ad(path)
    validation = validate_gene_programs(c_genes, s_genes, mode="custom")
    c_requested = validation.normalized_c_genes
    s_requested = validation.normalized_s_genes
    requested = c_requested + [g for g in s_requested if g.upper() not in {x.upper() for x in c_requested}]
    indices, present, _ = _gene_indices(adata, requested)
    lookup = {gene.upper(): i for i, gene in enumerate(present)}
    matrix = _dense(adata.X[:, indices])
    count_like = _is_count_like(adata.X)
    if count_like:
        matrix = np.log1p(np.maximum(matrix, 0.0))
    if matrix.ndim == 1:
        matrix = matrix[:, None]
    matrix[~np.isfinite(matrix)] = np.nan
    return adata, c_requested, s_requested, present, lookup, matrix, count_like, validation


def _validation_columns(validation: GeneProgramValidationResult) -> dict:
    return {
        "c_gene_count_requested": len(validation.requested_c_genes),
        "s_gene_count_requested": len(validation.requested_s_genes),
        "c_gene_count_used": len(validation.normalized_c_genes),
        "s_gene_count_used": len(validation.normalized_s_genes),
        "n_overlap_genes": validation.n_overlap_genes,
        "overlap_genes": ";".join(validation.overlap_genes),
        "overlap_policy": validation.overlap_policy,
        "program_validation_status": validation.validation_status,
    }


def _program_rows(program: str, requested: list[str], lookup: dict[str, int], present: list[str], matrix: np.ndarray) -> pd.DataFrame:
    rows: list[dict] = []
    available: list[tuple[str, str, np.ndarray]] = []
    for gene in requested:
        column = lookup.get(gene.upper())
        if column is not None:
            available.append((gene, present[column], matrix[:, column]))
    nonnegative = all(np.nanmin(values) >= 0 for _, _, values in available) if available else True
    weights = [float(np.nanmean(values if nonnegative else np.abs(values))) for _, _, values in available]
    total = float(np.nansum(weights))
    weight_by_gene = {requested_gene.upper(): weight for (requested_gene, _, _), weight in zip(available, weights)}
    for order, gene in enumerate(requested, 1):
        column = lookup.get(gene.upper())
        if column is None:
            rows.append({
                "program": program, "gene": gene, "requested_order": order, "status": "missing",
                "mean_expression": np.nan, "median_expression": np.nan, "detection_fraction": np.nan,
                "relative_contribution_percent": 0.0,
                "contribution_basis": "mean_transformed_expression" if nonnegative else "mean_absolute_transformed_expression",
            })
            continue
        values = matrix[:, column]
        weight = weight_by_gene[gene.upper()]
        rows.append({
            "program": program, "gene": present[column], "requested_order": order, "status": "present",
            "mean_expression": float(np.nanmean(values)), "median_expression": float(np.nanmedian(values)),
            "detection_fraction": float(np.nanmean(values > 0)),
            "relative_contribution_percent": 100.0 * weight / total if total > 0 else 0.0,
            "contribution_basis": "mean_transformed_expression" if nonnegative else "mean_absolute_transformed_expression",
        })
    return pd.DataFrame(rows)


def calculate_gene_composition(path: str | Path, c_genes: Iterable[str], s_genes: Iterable[str]) -> tuple[pd.DataFrame, dict]:
    adata, c_requested, s_requested, present, lookup, matrix, count_like, validation = _prepare(path, c_genes, s_genes)
    table = pd.concat([
        _program_rows("Cx", c_requested, lookup, present, matrix),
        _program_rows("Sx", s_requested, lookup, present, matrix),
    ], ignore_index=True)
    for column, value in _validation_columns(validation).items():
        table[column] = value
    metadata = {
        "sample": Path(path).stem, "source_h5ad": str(Path(path).resolve()),
        "n_spots": int(adata.n_obs), "expression_transform": "log1p_count_like" if count_like else "existing_processed_scale",
        "Cx_genes_requested": c_requested, "Sx_genes_requested": s_requested,
        "gene_program_validation": validation.to_provenance(),
        "definition": "Within-program percentage of each gene's mean transformed expression; absolute transformed expression is used when the existing processed scale contains negative values.",
    }
    return table, metadata


def _bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    result = np.full(len(pvalues), np.nan)
    valid = np.flatnonzero(np.isfinite(pvalues))
    if not len(valid):
        return result
    order = valid[np.argsort(pvalues[valid])]
    adjusted = pvalues[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result[order] = np.minimum(adjusted, 1.0)
    return result


def _hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    pooled_df = len(a) + len(b) - 2
    pooled_var = ((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / pooled_df
    if pooled_var <= 0 or not np.isfinite(pooled_var):
        return 0.0 if np.isclose(np.mean(a), np.mean(b)) else float("nan")
    d = (np.mean(a) - np.mean(b)) / math.sqrt(pooled_var)
    correction = 1.0 - 3.0 / (4.0 * (len(a) + len(b)) - 9.0)
    return float(correction * d)


def calculate_interface_enrichment(
    path: str | Path,
    c_genes: Iterable[str],
    s_genes: Iterable[str],
    c_q: float = 0.80,
    s_q: float = 0.80,
    g_q: float = 0.60,
) -> tuple[pd.DataFrame, dict]:
    from scipy.stats import mannwhitneyu

    adata, c_requested, s_requested, present, lookup, matrix, count_like, validation = _prepare(path, c_genes, s_genes)
    metrics, fields = score_h5ad(
        path,
        c_requested,
        s_requested,
        c_q,
        s_q,
        g_q,
        gene_program_mode="custom",
    )
    if not fields.get("spatial_available", False):
        raise ValueError(SPATIAL_QC_MESSAGE + " Interface Enrichment is unavailable.")
    interface = np.asarray(fields["interface"], dtype=bool)
    noninterface = ~interface
    program_map: dict[str, list[str]] = {"Cx": c_requested, "Sx": s_requested}
    rows: list[dict] = []
    for program, genes in program_map.items():
        present_genes = [g for g in genes if g.upper() in lookup]
        int_weights: list[float] = []
        non_weights: list[float] = []
        for gene in present_genes:
            values = matrix[:, lookup[gene.upper()]]
            shift = max(0.0, -float(np.nanmin(values))) if np.isfinite(values).any() else 0.0
            int_weights.append(float(np.nanmean(values[interface] + shift)) if interface.any() else np.nan)
            non_weights.append(float(np.nanmean(values[noninterface] + shift)) if noninterface.any() else np.nan)
        int_total, non_total = float(np.nansum(int_weights)), float(np.nansum(non_weights))
        for order, gene in enumerate(genes, 1):
            column = lookup.get(gene.upper())
            if column is None:
                rows.append({"program": program, "gene": gene, "requested_order": order, "status": "missing"})
                continue
            values = matrix[:, column]
            a, b = values[interface], values[noninterface]
            shift = max(0.0, -float(np.nanmin(values))) if np.isfinite(values).any() else 0.0
            a_shift, b_shift = a + shift, b + shift
            mean_a = float(np.nanmean(a)) if len(a) else np.nan
            mean_b = float(np.nanmean(b)) if len(b) else np.nan
            scale = float(np.nanmedian(np.abs(values[np.isfinite(values)]))) if np.isfinite(values).any() else 0.0
            pseudocount = max(np.finfo(float).eps, scale * 1e-6)
            fold = (float(np.nanmean(a_shift)) + pseudocount) / (float(np.nanmean(b_shift)) + pseudocount) if len(a) and len(b) else np.nan
            pvalue = np.nan
            if np.isfinite(a).sum() >= 2 and np.isfinite(b).sum() >= 2:
                try:
                    pvalue = float(mannwhitneyu(a[np.isfinite(a)], b[np.isfinite(b)], alternative="two-sided").pvalue)
                except ValueError:
                    pvalue = np.nan
            int_weight = float(np.nanmean(a_shift)) if len(a) else np.nan
            non_weight = float(np.nanmean(b_shift)) if len(b) else np.nan
            rows.append({
                "program": program, "gene": present[column], "requested_order": order, "status": "present",
                "n_interface": int(interface.sum()), "n_noninterface": int(noninterface.sum()),
                "mean_interface": mean_a, "mean_noninterface": mean_b,
                "mean_difference": mean_a - mean_b if np.isfinite(mean_a) and np.isfinite(mean_b) else np.nan,
                "fold_enrichment": fold, "log2_fold_enrichment": float(np.log2(fold)) if fold > 0 else np.nan,
                "hedges_g": _hedges_g(a, b), "mannwhitney_p": pvalue,
                "interface_composition_percent": 100 * int_weight / int_total if int_total > 0 else np.nan,
                "noninterface_composition_percent": 100 * non_weight / non_total if non_total > 0 else np.nan,
                "fold_change_basis": "transformed_expression" if shift == 0 else "transformed_expression_shifted_nonnegative",
            })
    table = pd.DataFrame(rows)
    table["mannwhitney_fdr_bh"] = _bh_adjust(table.get("mannwhitney_p", pd.Series(np.nan, index=table.index)).to_numpy())
    table["significant_fdr_0_05"] = table["mannwhitney_fdr_bh"] < 0.05
    metadata = {
        "sample": Path(path).stem, "source_h5ad": str(Path(path).resolve()),
        "n_spots": int(adata.n_obs), "n_interface": int(interface.sum()), "n_noninterface": int(noninterface.sum()),
        "interface_definition": "Unchanged core rule: high Cx AND high Sx AND high local R gradient.",
        "interface_quantiles": {"Cx": c_q, "Sx": s_q, "G": g_q}, "v0_1_regime_label": metrics["regime_label"],
        "statistical_test": "Two-sided Mann-Whitney U with Benjamini-Hochberg correction; reported only when both groups contain at least two finite values.",
        "expression_transform": "log1p_count_like" if count_like else "existing_processed_scale",
        "gene_program_validation": validation.to_provenance(),
    }
    return table, metadata


def _robust_unit(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    low, high = np.nanquantile(values, [0.05, 0.95])
    if not np.isfinite(low + high) or high <= low:
        return np.full(len(values), 0.5)
    return np.clip((values - low) / (high - low), 0.0, 1.0)


def _local_mean(values: np.ndarray, adjacency: list[list[int]]) -> np.ndarray:
    return np.asarray([np.mean(values[[i] + neighbors]) for i, neighbors in enumerate(adjacency)], dtype=float)


def _interaction_values(local_c: np.ndarray, local_s: np.ndarray, edges: list[tuple[int, int]]) -> dict[str, float]:
    epsilon = np.finfo(float).eps
    minimum, maximum = np.minimum(local_c, local_s), np.maximum(local_c, local_s)
    labels = np.full(len(local_c), "neutral", dtype=object)
    active = maximum >= 0.50
    labels[active & (np.abs(local_c - local_s) < 0.15)] = "mixed"
    labels[active & (local_c - local_s >= 0.15)] = "Cx_dominant"
    labels[active & (local_s - local_c >= 0.15)] = "Sx_dominant"
    informative = cross = 0
    for a, b in edges:
        if labels[a] in {"Cx_dominant", "Sx_dominant"} and labels[b] in {"Cx_dominant", "Sx_dominant"}:
            informative += 1
            cross += int(labels[a] != labels[b])
    return {
        "coexistence_index": float(np.mean(minimum)),
        "antagonism_index": float(np.mean(np.abs(local_c - local_s))),
        "balance_index": float(np.mean(1.0 - np.abs(local_c - local_s) / (local_c + local_s + epsilon))),
        "spatial_overlap_index": float(np.sum(minimum) / (np.sum(maximum) + epsilon)),
        "Cx_Sx_edge_mixing_fraction": float(cross / informative) if informative else np.nan,
        "Cx_Sx_cross_edges": int(cross), "informative_dominant_edges": int(informative),
    }


def calculate_spatial_interaction(
    path: str | Path,
    c_genes: Iterable[str],
    s_genes: Iterable[str],
    permutations: int = 499,
    seed: int = 20260705,
    c_q: float = 0.80,
    s_q: float = 0.80,
    g_q: float = 0.60,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    if permutations < 0:
        raise ValueError("Permutations must be zero or greater.")
    validation = validate_gene_programs(c_genes, s_genes, mode="custom")
    metrics, fields = score_h5ad(
        path,
        validation.normalized_c_genes,
        validation.normalized_s_genes,
        c_q,
        s_q,
        g_q,
        gene_program_mode="custom",
    )
    if not fields.get("spatial_available", False):
        raise ValueError(SPATIAL_QC_MESSAGE + " Cx/Sx Interaction is unavailable.")
    adata = _read_h5ad(path)
    spot_ids = np.asarray(adata.obs_names, dtype=str)
    coords = np.asarray(fields["coords"], dtype=float)
    edges, adjacency = _knn(coords, k=6)
    c_unit, s_unit = _robust_unit(fields["C"]), _robust_unit(fields["S"])
    local_c, local_s = _local_mean(c_unit, adjacency), _local_mean(s_unit, adjacency)
    observed = _interaction_values(local_c, local_s, edges)
    null: dict[str, list[float]] = {key: [] for key in (
        "coexistence_index", "antagonism_index", "balance_index", "spatial_overlap_index", "Cx_Sx_edge_mixing_fraction"
    )}
    rng = np.random.default_rng(seed)
    for _ in range(permutations):
        permuted_s = _local_mean(rng.permutation(s_unit), adjacency)
        values = _interaction_values(local_c, permuted_s, edges)
        for key in null:
            if np.isfinite(values[key]):
                null[key].append(values[key])
    row: dict[str, float | int | str] = {
        "sample": metrics["sample"], "n_spots": metrics["n_spots"], "n_edges": len(edges),
        "neighbors_k": 6, "permutations": permutations, "seed": seed, **observed,
    }
    for key, distribution in null.items():
        values = np.asarray(distribution, dtype=float)
        row[f"{key}_null_mean"] = float(np.mean(values)) if len(values) else np.nan
        row[f"{key}_null_sd"] = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
        row[f"{key}_z"] = ((float(observed[key]) - float(np.mean(values))) / float(np.std(values, ddof=1))) if len(values) > 1 and np.std(values, ddof=1) > 0 and np.isfinite(observed[key]) else np.nan
        row[f"{key}_permutation_p_two_sided"] = ((1 + np.sum(np.abs(values - np.mean(values)) >= abs(float(observed[key]) - np.mean(values)))) / (len(values) + 1)) if len(values) and np.isfinite(observed[key]) else np.nan
    epsilon = np.finfo(float).eps
    spot_table = pd.DataFrame({
        "spot": spot_ids, "spot_index": np.arange(len(coords)), "x": coords[:, 0], "y": coords[:, 1],
        "Cx_unit": c_unit, "Sx_unit": s_unit, "local_Cx": local_c, "local_Sx": local_s,
        "local_coexistence": np.minimum(local_c, local_s),
        "local_antagonism": np.abs(local_c - local_s),
        "local_balance": 1.0 - np.abs(local_c - local_s) / (local_c + local_s + epsilon),
        "v0_1_interface": np.asarray(fields["interface"], dtype=bool),
    })
    metadata = {
        "sample": metrics["sample"], "source_h5ad": metrics["source_h5ad"],
        "normalization": "Each core Cx/Sx score is clipped and scaled between its sample 5th and 95th percentiles.",
        "neighborhood": "Undirected 6-nearest-neighbor graph; local fields include each spot and its graph neighbors.",
        "metric_definitions": {
            "coexistence_index": "Mean min(local Cx, local Sx).",
            "antagonism_index": "Mean absolute difference between local Cx and local Sx.",
            "balance_index": "Mean 1-|local Cx-local Sx|/(local Cx+local Sx).",
            "spatial_overlap_index": "Weighted Jaccard: sum min(local Cx,local Sx)/sum max(local Cx,local Sx).",
            "Cx_Sx_edge_mixing_fraction": "Fraction of informative dominance edges joining opposite Cx- and Sx-dominant neighborhoods.",
        },
        "null_model": "Seeded permutation of Sx activation across fixed coordinates, followed by neighborhood recomputation.",
        "gene_program_validation": validation.to_provenance(),
    }
    plot_fields = {"coords": coords, "local_c": local_c, "local_s": local_s, "spot_table": spot_table}
    return pd.DataFrame([row]), spot_table, metadata, plot_fields


def _save_figure(fig, base: Path) -> tuple[Path, Path]:
    png, pdf = base.with_suffix(".png"), base.with_suffix(".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    return png, pdf


def _plot_composition(table: pd.DataFrame, sample: str, base: Path) -> tuple[Path, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, max(4.0, 0.42 * max(6, table.groupby("program").size().max()))))
    colors = {"Cx": "#176B87", "Sx": "#C65D3B"}
    for ax, program in zip(axes, ("Cx", "Sx")):
        data = table[table["program"] == program].sort_values("relative_contribution_percent")
        ax.barh(data["gene"], data["relative_contribution_percent"], color=colors[program])
        ax.set_title(f"{program} gene composition")
        ax.set_xlabel("Relative contribution (%)")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{sample} — Cx/Sx gene composition", fontweight="bold")
    fig.tight_layout()
    paths = _save_figure(fig, base)
    plt.close(fig)
    return paths


def _plot_enrichment(table: pd.DataFrame, sample: str, base: Path) -> tuple[Path, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    data = table[table["status"] == "present"].copy().sort_values("log2_fold_enrichment")
    colors = np.where(data["program"].eq("Cx"), "#176B87", "#C65D3B")
    colors = np.where(data["significant_fdr_0_05"], colors, "#A8B0B8")
    fig, ax = plt.subplots(figsize=(8.2, max(4.5, 0.38 * max(8, len(data)))))
    ax.barh(data["gene"], data["log2_fold_enrichment"], color=colors)
    ax.axvline(0, color="#30343B", lw=0.9)
    ax.set_xlabel("log2 fold enrichment (interface / non-interface)")
    ax.set_title(f"{sample} — interface enrichment")
    ax.legend(
        handles=[Patch(color="#176B87", label="Cx, FDR < 0.05"), Patch(color="#C65D3B", label="Sx, FDR < 0.05"), Patch(color="#A8B0B8", label="FDR ≥ 0.05")],
        frameon=False,
        loc="lower right",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    paths = _save_figure(fig, base)
    plt.close(fig)
    return paths


def _plot_interaction(summary: pd.DataFrame, fields: dict, sample: str, base: Path) -> tuple[Path, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coords, table = fields["coords"], fields["spot_table"]
    panels = [
        (fields["local_c"], "Local Cx", "Blues", 0, 1),
        (fields["local_s"], "Local Sx", "Oranges", 0, 1),
        (table["local_balance"], "Local balance", "viridis", 0, 1),
        (table["local_coexistence"], "Local coexistence", "magma", 0, 1),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.2))
    for ax, (values, title, cmap, vmin, vmax) in zip(axes.ravel(), panels):
        image = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=12, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)
        ax.set_title(title); ax.set_aspect("equal", adjustable="datalim"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(image, ax=ax, fraction=.046, pad=.03)
    row = summary.iloc[0]
    fig.suptitle(
        f"{sample} — Cx/Sx spatial interaction\ncoexistence={row['coexistence_index']:.3f}  "
        f"antagonism={row['antagonism_index']:.3f}  overlap={row['spatial_overlap_index']:.3f}", fontweight="bold"
    )
    fig.tight_layout(rect=[0, 0, 1, .94])
    paths = _save_figure(fig, base)
    plt.close(fig)
    return paths


def _write_metadata(path: Path, module: str, metadata: dict, parameters: dict) -> None:
    payload = {
        "application": "SpatialTX Studio Desktop", "version": __version__, "module": MODULE_LABELS[module],
        "created": dt.datetime.now().isoformat(timespec="seconds"), "parameters": parameters, **metadata,
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def run_advanced_batch(
    module: str,
    paths: list[Path],
    output_root: str | Path,
    c_genes: Iterable[str],
    s_genes: Iterable[str],
    progress: Progress | None = None,
    c_q: float = .80,
    s_q: float = .80,
    g_q: float = .60,
    permutations: int = 499,
    seed: int = 20260705,
) -> tuple[Path, pd.DataFrame]:
    if module not in MODULE_LABELS:
        raise ValueError(f"Unknown advanced analysis module: {module}")
    if not paths:
        raise ValueError("Select at least one h5ad sample.")
    validation = validate_gene_programs(c_genes, s_genes, mode="custom")
    c_genes = validation.normalized_c_genes
    s_genes = validation.normalized_s_genes
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_dir = Path(output_root).expanduser().resolve() / f"advanced_{module}_{stamp}"
    run_dir.mkdir(parents=True, exist_ok=False)
    manifest_rows: list[dict] = []
    for number, source in enumerate(paths, 1):
        sample = source.stem
        sample_dir = run_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        if progress:
            progress(f"[{number}/{len(paths)}] {MODULE_LABELS[module]}: {source.name}")
        try:
            parameters = {
                "Cx_genes": list(c_genes),
                "Sx_genes": list(s_genes),
                "gene_program_validation": validation.to_provenance(),
            }
            if module == "composition":
                table, metadata = calculate_gene_composition(source, c_genes, s_genes)
                for column, value in _validation_columns(validation).items():
                    table[column] = value
                metadata["gene_program_validation"] = validation.to_provenance()
                csv_path = sample_dir / "gene_composition.csv"; table.to_csv(csv_path, index=False)
                png, pdf = _plot_composition(table, sample, sample_dir / "gene_composition")
            elif module == "enrichment":
                table, metadata = calculate_interface_enrichment(source, c_genes, s_genes, c_q, s_q, g_q)
                metadata["gene_program_validation"] = validation.to_provenance()
                csv_path = sample_dir / "interface_enrichment.csv"; table.to_csv(csv_path, index=False)
                png, pdf = _plot_enrichment(table, sample, sample_dir / "interface_enrichment")
                parameters["interface_quantiles"] = {"Cx": c_q, "Sx": s_q, "G": g_q}
            else:
                table, spots, metadata, plot_fields = calculate_spatial_interaction(
                    source, c_genes, s_genes, permutations, seed, c_q, s_q, g_q
                )
                metadata["gene_program_validation"] = validation.to_provenance()
                csv_path = sample_dir / "interaction_summary.csv"; table.to_csv(csv_path, index=False)
                spots.to_csv(sample_dir / "interaction_spot_metrics.csv", index=False)
                png, pdf = _plot_interaction(table, plot_fields, sample, sample_dir / "cx_sx_interaction")
                parameters.update({"permutations": permutations, "seed": seed, "neighbors_k": 6})
            _write_metadata(sample_dir / "analysis_metadata.json", module, metadata, parameters)
            manifest_rows.append({
                "sample": sample, "source_h5ad": str(source.resolve()), "status": "ok", "table_csv": str(csv_path),
                "figure_png": str(png), "figure_pdf": str(pdf), "metadata_json": str(sample_dir / "analysis_metadata.json"),
            })
        except Exception as exc:
            manifest_rows.append({"sample": sample, "source_h5ad": str(source.resolve()), "status": f"error: {exc}"})
            if progress:
                progress(f"  Error: {exc}")
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(run_dir / "advanced_analysis_manifest.csv", index=False)
    if progress:
        progress(f"Completed: {run_dir}")
    return run_dir, manifest
