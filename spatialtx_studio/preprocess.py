from __future__ import annotations

import numpy as np
import scanpy as sc


def preprocess_adata(
    adata,
    min_genes_per_spot: int = 200,
    gene_min_spots: int = 10,
    use_in_tissue: bool = True,
    target_sum: float = 10000.0,
):
    a = adata.copy()
    a.var_names_make_unique()

    if use_in_tissue and "in_tissue" in a.obs.columns:
        before = a.n_obs
        a = a[a.obs["in_tissue"].astype(int) == 1].copy()
        if a.n_obs == 0:
            raise ValueError(f"in_tissue filtering removed all spots; before={before}")

    if "n_genes_by_counts" not in a.obs.columns:
        x = a.X
        detected = np.asarray((x > 0).sum(axis=1)).ravel()
        a.obs["n_genes_by_counts"] = detected

    before = a.n_obs
    a = a[a.obs["n_genes_by_counts"] > min_genes_per_spot].copy()
    if a.n_obs == 0:
        raise ValueError(
            "Spot QC removed all spots: "
            f"n_genes_by_counts>{min_genes_per_spot}, before={before}"
        )

    sc.pp.filter_genes(a, min_cells=gene_min_spots)
    if a.n_vars == 0:
        raise ValueError(f"Gene filtering removed all genes: min_cells={gene_min_spots}")

    sc.pp.normalize_total(a, target_sum=target_sum)
    sc.pp.log1p(a)
    return a
