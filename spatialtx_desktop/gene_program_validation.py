from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Iterable


GeneInput = str | Iterable[str]


def _canonical_entries(genes: GeneInput) -> list[str]:
    if isinstance(genes, str):
        values = genes.replace(";", ",").replace("\n", ",").split(",")
    else:
        values = list(genes)
    return [str(value).strip().upper() for value in values if str(value).strip()]


def _normalize_with_duplicates(genes: GeneInput) -> tuple[list[str], list[str], list[str]]:
    requested = _canonical_entries(genes)
    normalized: list[str] = []
    duplicates: list[str] = []
    seen: set[str] = set()
    duplicate_seen: set[str] = set()
    for gene in requested:
        if gene in seen:
            if gene not in duplicate_seen:
                duplicates.append(gene)
                duplicate_seen.add(gene)
            continue
        seen.add(gene)
        normalized.append(gene)
    return requested, normalized, duplicates


def normalize_gene_list(genes: GeneInput) -> list[str]:
    """Return uppercase, trimmed, non-empty, order-preserving unique symbols."""
    return _normalize_with_duplicates(genes)[1]


def find_gene_program_overlap(c_genes: GeneInput, s_genes: GeneInput) -> list[str]:
    """Return canonical overlap symbols in C-side input order."""
    normalized_c = normalize_gene_list(c_genes)
    s_set = set(normalize_gene_list(s_genes))
    return [gene for gene in normalized_c if gene in s_set]


def s_side_biological_warnings(s_genes: GeneInput) -> list[str]:
    genes = normalize_gene_list(s_genes)
    immune_like = [gene for gene in genes if gene.startswith(("IGH", "IGL", "IGK"))]
    if not immune_like:
        return []
    return [
        "Some S-side genes appear immune-associated rather than stromal-associated. "
        "Review the program composition before interpretation. "
        f"Flagged genes: {', '.join(immune_like)}."
    ]


@dataclass
class GeneProgramValidationResult:
    requested_c_genes: list[str]
    requested_s_genes: list[str]
    normalized_c_genes: list[str]
    normalized_s_genes: list[str]
    overlap_genes: list[str]
    n_overlap_genes: int
    c_duplicates_removed: list[str]
    s_duplicates_removed: list[str]
    action_taken: str
    validation_status: str
    mode: str
    overlap_policy: str
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def to_provenance(self) -> dict:
        return {
            "c_genes_requested": list(self.requested_c_genes),
            "s_genes_requested": list(self.requested_s_genes),
            "c_genes_used": list(self.normalized_c_genes),
            "s_genes_used": list(self.normalized_s_genes),
            "overlap_genes": list(self.overlap_genes),
            "n_overlap_genes": int(self.n_overlap_genes),
            "c_duplicates_removed": list(self.c_duplicates_removed),
            "s_duplicates_removed": list(self.s_duplicates_removed),
            "overlap_policy": self.overlap_policy,
            "validation_status": self.validation_status,
            "action_taken": self.action_taken,
            "mode": self.mode,
            "warnings": list(self.warnings),
        }


class GeneProgramOverlapError(ValueError):
    def __init__(self, message: str, result: GeneProgramValidationResult) -> None:
        super().__init__(message)
        self.result = result


def _overlap_message(overlap: list[str], mode: str) -> str:
    joined = ", ".join(overlap)
    if mode == "fixed":
        return (
            f"C/S gene-program overlap detected: {joined}.\n"
            "Fixed gene programs must be mutually exclusive."
        )
    if mode in {"adaptive", "semi_auto", "cancer_adaptive", "qubo"}:
        return (
            f"C/S gene-program overlap detected after {mode} selection: {joined}.\n"
            "This indicates a selection-logic error; the analysis was not run."
        )
    return (
        "The following genes are present in both C-side and S-side programs:\n"
        f"{joined}\n\n"
        "SpatialTX requires mutually exclusive C/S programs because R = C - S.\n"
        "Remove each overlapping gene from one side before continuing."
    )


def validate_gene_programs(
    c_genes: GeneInput,
    s_genes: GeneInput,
    mode: str,
    overlap_policy: str = "error",
) -> GeneProgramValidationResult:
    """Canonicalize C/S lists and enforce mutually exclusive programs.

    ``overlap_policy='report'`` is audit-only and returns an invalid result.
    All analysis entry points use the default hard-error policy.
    """
    policy = str(overlap_policy).strip().lower()
    if policy not in {"error", "report"}:
        raise ValueError("overlap_policy must be 'error' or 'report'.")
    normalized_mode = str(mode or "core").strip().lower()
    requested_c, normalized_c, c_duplicates = _normalize_with_duplicates(c_genes)
    requested_s, normalized_s, s_duplicates = _normalize_with_duplicates(s_genes)
    s_set = set(normalized_s)
    overlap = [gene for gene in normalized_c if gene in s_set]
    if not normalized_c or not normalized_s:
        raise ValueError("Both C-side and S-side gene programs must contain at least one gene.")
    valid = not overlap
    if overlap:
        action = "analysis_blocked_overlap"
    elif c_duplicates or s_duplicates:
        action = "duplicates_removed_and_validated"
    else:
        action = "normalized_and_validated"
    result = GeneProgramValidationResult(
        requested_c_genes=requested_c,
        requested_s_genes=requested_s,
        normalized_c_genes=normalized_c,
        normalized_s_genes=normalized_s,
        overlap_genes=overlap,
        n_overlap_genes=len(overlap),
        c_duplicates_removed=c_duplicates,
        s_duplicates_removed=s_duplicates,
        action_taken=action,
        validation_status="valid" if valid else "invalid_overlap",
        mode=normalized_mode,
        overlap_policy=policy,
        warnings=s_side_biological_warnings(normalized_s),
    )
    if overlap and policy == "error":
        raise GeneProgramOverlapError(_overlap_message(overlap, normalized_mode), result)
    return result
