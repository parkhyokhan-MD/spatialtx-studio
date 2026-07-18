from __future__ import annotations

import numpy as np

from spatialtx_desktop.gene_program_validation import normalize_gene_list, validate_gene_programs

C_FIXED = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG"]
B_FIXED = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN"]

C_CANDIDATES = ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "IFNG", "GNLY", "CCL5", "CTSW", "TRAC"]
B_CANDIDATES = ["COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN", "POSTN", "ACTA2", "TAGLN", "FAP", "THY1"]

_FIXED_PROGRAM_VALIDATION = validate_gene_programs(C_FIXED, B_FIXED, mode="fixed")


def _var_lookup(adata) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for index, gene in enumerate(adata.var_names):
        lookup.setdefault(str(gene).strip().upper(), int(index))
    return lookup


def present_genes(adata, genes: list[str]) -> list[str]:
    lookup = _var_lookup(adata)
    return [gene for gene in normalize_gene_list(genes) if gene in lookup]


def select_gene_programs(adata, config: dict, mode: str) -> tuple[list[str], list[str], str]:
    mode = mode.lower()
    topk = int(config["selection"]["topK"])
    min_genes = int(config["selection"]["min_genes_per_program"])
    programs = config["gene_programs"]

    if mode == "fixed":
        requested = validate_gene_programs(
            programs.get("C_FIXED", C_FIXED),
            programs.get("B_FIXED", B_FIXED),
            mode="fixed",
        )
        c_genes = present_genes(adata, requested.normalized_c_genes)
        b_genes = present_genes(adata, requested.normalized_s_genes)
        note = "fixed_gene_program"
    elif mode in {"semi_auto", "semi"}:
        c_genes = variance_rank_genes(adata, programs.get("C_CANDIDATES", C_CANDIDATES), topk)
        b_candidates = normalize_gene_list(programs.get("B_CANDIDATES", B_CANDIDATES))
        excluded = [gene for gene in b_candidates if gene in set(c_genes)]
        b_genes = variance_rank_genes(adata, b_candidates, topk, exclude_genes=c_genes)
        note = (
            "semi_auto_variance_topK;overlap_constraint_enabled=true;"
            f"genes_excluded_due_to_opposite_side={','.join(excluded)}"
        )
    elif mode == "cancer_adaptive":
        c_genes = variance_rank_genes(adata, programs.get("C_CANDIDATES", C_CANDIDATES), topk)
        b_candidates = normalize_gene_list(programs.get("B_CANDIDATES", B_CANDIDATES))
        excluded = [gene for gene in b_candidates if gene in set(c_genes)]
        b_genes = variance_rank_genes(adata, b_candidates, topk, exclude_genes=c_genes)
        note = (
            "cancer_adaptive_v0_1_variance_topK_exploratory;overlap_constraint_enabled=true;"
            f"genes_excluded_due_to_opposite_side={','.join(excluded)}"
        )
    elif mode == "custom":
        requested = validate_gene_programs(
            programs.get("CUSTOM_C", []),
            programs.get("CUSTOM_B", []),
            mode="custom",
        )
        c_genes = present_genes(adata, requested.normalized_c_genes)
        b_genes = present_genes(adata, requested.normalized_s_genes)
        note = "custom_config_gene_program"
    else:
        raise ValueError(f"Unsupported gene mode: {mode}")

    if len(c_genes) < min_genes:
        raise ValueError("insufficient_fixed_gene_availability_immune")
    if len(b_genes) < min_genes:
        raise ValueError("insufficient_fixed_gene_availability_stroma")
    final_mode = "adaptive" if mode in {"semi_auto", "semi", "cancer_adaptive"} else mode
    validated = validate_gene_programs(c_genes, b_genes, mode=final_mode)
    return validated.normalized_c_genes, validated.normalized_s_genes, note


def variance_rank_genes(
    adata,
    genes: list[str],
    topk: int,
    exclude_genes: list[str] | None = None,
) -> list[str]:
    excluded = set(normalize_gene_list(exclude_genes or []))
    present = [gene for gene in present_genes(adata, genes) if gene not in excluded]
    scored = []
    for gene in present:
        x = dense_vector(adata, gene)
        scored.append((gene, float(np.nanvar(x))))
    scored.sort(key=lambda item: item[1], reverse=True)
    return [gene for gene, _ in scored[:topk]]


def dense_vector(adata, gene: str) -> np.ndarray:
    key = normalize_gene_list([gene])
    if not key or key[0] not in _var_lookup(adata):
        raise KeyError(f"Gene not found in AnnData: {gene}")
    x = adata.X[:, _var_lookup(adata)[key[0]]]
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).ravel().astype(float)


def mean_program_score(adata, genes: list[str]) -> np.ndarray:
    present = present_genes(adata, genes)
    if not present:
        return np.zeros(adata.n_obs, dtype=float)
    values = [dense_vector(adata, gene) for gene in present]
    return np.nanmean(np.vstack(values), axis=0)
