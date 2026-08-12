# SpatialTX Studio Desktop v0.6-dev-HV-validation

> H/V computational and audit validation layer. Exploratory research use only; this is not H/V biological validation.

## Scope

This patch makes the existing observational H/V axes auditable and more sensitive to focal upper-tail context without changing the v0.6 C/S + FRAME2.6 result.

## Added

- Effective default or user-supplied H/V program metadata.
- `context_gene_audit.csv` with per-sample matched, missing, expressed, coverage, scale, source, status, and warnings.
- Explicit availability statuses that distinguish valid zero from unavailable `NaN`.
- Raw mean, median, q75, q90, transition enrichment, and coefficient of variation.
- One pooled q90 threshold per axis per Pre/Post pair.
- High-context and graph-local high-context fractions using the shared within-pair threshold.
- Expanded `context_changes.csv`, `multiaxial_pair_summary.csv`, run metadata, and compact H/V UI detail.
- Separate provenance-labelled warnings and high fractions for legacy single-sample centered-context q80 versus pair-pooled raw-context q90 audit layers.
- A two-tier H/V area in `multiaxial_pair_overview.png`: raw-median Delta above and pair-pooled high-context fraction Delta below, each with its own x-axis.
- Comparative cache schema `v0.6-hv-validation-v1`.

## Preserved

The patch does not change C/S programs, `R=C-S`, quantiles, FRAME2.6, Type A/B/C, interface/diffuse/burden, adjacency/fragmentation, comparability QC, site warnings, or any composite/clinical inference boundary.
