from __future__ import annotations

import json
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping

import pandas as pd

from .models import ReliabilityConfig


ENSEMBL_VERSION_PATTERN = re.compile(r"^(ENS[A-Z]*G\d+)\.\d+$", re.IGNORECASE)
AUDIT_COLUMNS = [
    "canonical_gene",
    "original_entries",
    "axes",
    "poles",
    "overlap_type",
    "severity",
    "action",
    "normalization_rule",
]


class CrossExclusivityError(ValueError):
    def __init__(self, message: str, audit: pd.DataFrame) -> None:
        super().__init__(message)
        self.audit = audit


def _base_canonical(value: str) -> str:
    normalized = str(value).strip().upper()
    match = ENSEMBL_VERSION_PATTERN.match(normalized)
    return match.group(1).upper() if match else normalized


def canonical_alias_map(aliases: Mapping[str, str] | None) -> dict[str, str]:
    return {
        _base_canonical(source): _base_canonical(target)
        for source, target in dict(aliases or {}).items()
        if _base_canonical(source) and _base_canonical(target)
    }


def canonicalize_gene(value: str, aliases: Mapping[str, str] | None = None) -> str:
    gene = _base_canonical(value)
    alias_lookup = canonical_alias_map(aliases)
    visited: set[str] = set()
    while gene in alias_lookup and gene not in visited:
        visited.add(gene)
        gene = alias_lookup[gene]
    return gene


def normalization_rule(config: ReliabilityConfig) -> str:
    return (
        "trim|uppercase|ensembl_version_removed|"
        f"alias_map:{config.canonicalization_source}:{config.canonicalization_version}"
    )


def audit_cross_exclusivity(
    programs: Mapping[str, Mapping[str, Iterable[str]]],
    config: ReliabilityConfig | Mapping | None = None,
    *,
    raise_on_error: bool = True,
) -> pd.DataFrame:
    """Audit duplicate canonical genes within and across every axis/pole.

    Only exact canonical identity is considered. Paralog and pathway
    relationships are intentionally not inferred.
    """

    cfg = ReliabilityConfig.from_value(config)
    cfg.validate()
    occurrences: dict[str, list[dict]] = defaultdict(list)
    for axis, poles in programs.items():
        for pole, entries in poles.items():
            for original in entries:
                if not str(original).strip():
                    continue
                canonical = canonicalize_gene(str(original), cfg.canonical_aliases)
                if canonical:
                    occurrences[canonical].append({
                        "original": str(original),
                        "axis": str(axis),
                        "pole": str(pole),
                    })

    rows: list[dict] = []
    duplicate_genes: list[str] = []
    for canonical in sorted(occurrences):
        items = occurrences[canonical]
        axes = sorted({item["axis"] for item in items})
        poles = sorted({item["pole"] for item in items})
        if len(items) <= 1:
            overlap_type = "none"
        elif len(axes) > 1:
            overlap_type = "cross_axis"
        elif len(poles) > 1:
            overlap_type = "within_axis_cross_pole"
        else:
            overlap_type = "within_pole"
        if overlap_type != "none":
            duplicate_genes.append(canonical)
        if overlap_type == "none":
            severity, action = "ok", "retained"
        elif cfg.strict_cross_exclusivity:
            severity, action = "hard_error", "analysis_blocked"
        else:
            severity, action = "warning", "reported_not_transformed"
        rows.append({
            "canonical_gene": canonical,
            "original_entries": json.dumps(
                [item["original"] for item in items], ensure_ascii=False
            ),
            "axes": ";".join(axes),
            "poles": ";".join(poles),
            "overlap_type": overlap_type,
            "severity": severity,
            "action": action,
            "normalization_rule": normalization_rule(cfg),
        })
    audit = pd.DataFrame(rows, columns=AUDIT_COLUMNS)
    if duplicate_genes and cfg.strict_cross_exclusivity and raise_on_error:
        raise CrossExclusivityError(
            "Strict cross-exclusivity failed for canonical genes: "
            + ", ".join(duplicate_genes),
            audit,
        )
    return audit


def gene_coverage_audit(
    programs: Mapping[str, Mapping[str, Iterable[str]]],
    available_genes: Iterable[str],
    config: ReliabilityConfig | Mapping | None = None,
    *,
    pair_label: str = "",
    sample_role: str = "",
    sample_file: str = "",
) -> pd.DataFrame:
    cfg = ReliabilityConfig.from_value(config)
    cfg.validate()
    available = {
        canonicalize_gene(gene, cfg.canonical_aliases)
        for gene in available_genes
        if str(gene).strip()
    }
    rows: list[dict] = []
    for axis, poles in programs.items():
        for pole, entries in poles.items():
            requested: list[str] = []
            for entry in entries:
                canonical = canonicalize_gene(entry, cfg.canonical_aliases)
                if canonical and canonical not in requested:
                    requested.append(canonical)
            present = [gene for gene in requested if gene in available]
            missing = [gene for gene in requested if gene not in available]
            fraction = len(present) / len(requested) if requested else 0.0
            if not requested or not present:
                status, validity = "no_genes_present", "invalid"
            elif fraction < cfg.gene_coverage_low:
                status, validity = "low", "invalid"
            elif fraction < cfg.gene_coverage_caution:
                status, validity = "caution", "warning"
            else:
                status, validity = "available", "valid"
            rows.append({
                "pair_label": pair_label,
                "sample_role": sample_role,
                "sample_file": sample_file,
                "axis": str(axis),
                "pole": str(pole),
                "n_genes_requested": len(requested),
                "n_genes_present": len(present),
                "gene_coverage_fraction": fraction,
                "present_genes": ";".join(present),
                "missing_genes": ";".join(missing),
                "coverage_status": status,
                "score_validity": validity,
                "normalization_rule": normalization_rule(cfg),
            })
    return pd.DataFrame(rows)
