from __future__ import annotations

import numpy as np

C_FIXED = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG"]
B_FIXED = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN"]

C_CANDIDATES = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG", "GNLY", "CCL5", "CTSW", "TRAC"]
B_CANDIDATES = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN", "POSTN", "ACTA2", "TAGLN", "FAP", "THY1"]


def present_genes(adata, genes: list[str]) -> list[str]:
    return [gene for gene in genes if gene in adata.var_names]


def select_gene_programs(adata, config: dict, mode: str) -> tuple[list[str], list[str], str]:
    mode = mode.lower()
    topk = int(config["selection"]["topK"])
    min_genes = int(config["selection"]["min_genes_per_program"])
    programs = config["gene_programs"]

    if mode == "fixed":
        c_genes = present_genes(adata, programs.get("C_FIXED", C_FIXED))
        b_genes = present_genes(adata, programs.get("B_FIXED", B_FIXED))
        note = "fixed_gene_program"
    elif mode in {"semi_auto", "semi"}:
        c_genes = variance_rank_genes(adata, programs.get("C_CANDIDATES", C_CANDIDATES), topk)
        b_genes = variance_rank_genes(adata, programs.get("B_CANDIDATES", B_CANDIDATES), topk)
        note = "semi_auto_variance_topK"
    elif mode == "cancer_adaptive":
        c_genes = variance_rank_genes(adata, programs.get("C_CANDIDATES", C_CANDIDATES), topk)
        b_genes = variance_rank_genes(adata, programs.get("B_CANDIDATES", B_CANDIDATES), topk)
        note = "cancer_adaptive_v0_1_variance_topK_exploratory"
    elif mode == "custom":
        c_genes = present_genes(adata, programs.get("CUSTOM_C", []))
        b_genes = present_genes(adata, programs.get("CUSTOM_B", []))
        note = "custom_config_gene_program"
    else:
        raise ValueError(f"Unsupported gene mode: {mode}")

    if len(c_genes) < min_genes:
        raise ValueError("insufficient_fixed_gene_availability_immune")
    if len(b_genes) < min_genes:
        raise ValueError("insufficient_fixed_gene_availability_stroma")
    return c_genes, b_genes, note


def variance_rank_genes(adata, genes: list[str], topk: int) -> list[str]:
    present = present_genes(adata, genes)
    scored = []
    for gene in present:
        x = dense_vector(adata, gene)
        scored.append((gene, float(np.nanvar(x))))
    scored.sort(key=lambda item: item[1], reverse=True)
    return [gene for gene, _ in scored[:topk]]


def dense_vector(adata, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).ravel().astype(float)


def mean_program_score(adata, genes: list[str]) -> np.ndarray:
    present = present_genes(adata, genes)
    if not present:
        return np.zeros(adata.n_obs, dtype=float)
    values = [dense_vector(adata, gene) for gene in present]
    return np.nanmean(np.vstack(values), axis=0)
