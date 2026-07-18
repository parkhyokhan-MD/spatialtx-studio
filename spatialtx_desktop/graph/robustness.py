from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from .builder import GraphBuildConfig, build_spatial_graph
from .enrichment import categorical_neighborhood_enrichment


def compare_graph_robustness(
    adata,
    primary_config: GraphBuildConfig,
    labels,
    *,
    permutations: int,
    seed: int,
    permutation_scope: str = "whole_slide",
    strata=None,
    tissue_mask=None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare association stability without replacing the primary graph result."""
    configs = [
        replace(primary_config, method="radius"),
        replace(primary_config, method="lattice"),
        replace(primary_config, method="knn"),
    ]
    graph_rows: list[dict] = []
    association_tables: list[pd.DataFrame] = []
    for config in configs:
        result = build_spatial_graph(adata, config)
        graph_rows.append({
            "requested_graph": config.method,
            "effective_graph": result.method,
            "status": "ok" if int(result.qc.get("n_edges", 0)) > 0 else "insufficient_edges",
            "requested_radius": result.metadata.get("requested_radius"),
            "effective_radius": result.metadata.get("effective_radius"),
            "radius_unit": result.metadata.get("radius_unit"),
            "n_edges": result.qc.get("n_edges"),
            "isolated_fraction": result.qc.get("isolated_fraction"),
            "n_connected_components": result.qc.get("n_connected_components"),
            "largest_component_ratio": result.qc.get("largest_component_ratio"),
            "warnings": "; ".join(result.warnings),
        })
        if int(result.qc.get("n_edges", 0)) < 1:
            continue
        table = categorical_neighborhood_enrichment(
            result.connectivities,
            labels,
            permutations=permutations,
            seed=seed,
            permutation_scope=permutation_scope,
            strata=strata,
            tissue_mask=tissue_mask,
        )
        table.insert(0, "effective_graph", result.method)
        table.insert(0, "requested_graph", config.method)
        association_tables.append(table)
    graphs = pd.DataFrame(graph_rows)
    if len(graphs):
        graphs["isolated_fraction_variation"] = float(graphs["isolated_fraction"].max() - graphs["isolated_fraction"].min())
        graphs["connected_component_variation"] = float(graphs["n_connected_components"].max() - graphs["n_connected_components"].min())
        graphs["valid_edge_variation"] = float(graphs["n_edges"].max() - graphs["n_edges"].min())
    associations = pd.concat(association_tables, ignore_index=True) if association_tables else pd.DataFrame()
    stability_rows: list[dict] = []
    if len(associations):
        for (label_a, label_b), group in associations.groupby(["label_a", "label_b"], dropna=False):
            directions = group["association_direction"].astype(str).tolist()
            significant = (pd.to_numeric(group["fdr_bh"], errors="coerce") <= 0.05).fillna(False)
            ratios = pd.to_numeric(group["observed_expected_ratio"], errors="coerce")
            zscores = pd.to_numeric(group["enrichment_z"], errors="coerce")
            stability_rows.append({
                "label_a": label_a,
                "label_b": label_b,
                "graphs_evaluated": int(len(group)),
                "requested_graphs": ";".join(group["requested_graph"].astype(str)),
                "effective_graphs": ";".join(group["effective_graph"].astype(str)),
                "association_directions": ";".join(directions),
                "direction_stable_across_graphs": len(set(directions)) == 1,
                "significance_states": ";".join("significant" if value else "not_significant" for value in significant),
                "significance_stable_across_graphs": len(set(bool(value) for value in significant)) == 1,
                "observed_expected_ratio_min": float(ratios.min()) if ratios.notna().any() else np.nan,
                "observed_expected_ratio_max": float(ratios.max()) if ratios.notna().any() else np.nan,
                "observed_expected_ratio_variation": float(ratios.max() - ratios.min()) if ratios.notna().any() else np.nan,
                "z_score_min": float(zscores.min()) if zscores.notna().any() else np.nan,
                "z_score_max": float(zscores.max()) if zscores.notna().any() else np.nan,
                "z_score_variation": float(zscores.max() - zscores.min()) if zscores.notna().any() else np.nan,
                "interpretation": "exploratory graph robustness; not biological validation",
            })
    return graphs, pd.DataFrame(stability_rows)
