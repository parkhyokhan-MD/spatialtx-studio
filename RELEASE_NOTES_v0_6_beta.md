# SpatialTX Studio Desktop v0.6-beta

Release date: 2026-08-12

> Public source beta for exploratory research use only. Not intended for diagnosis, treatment selection, response prediction, or clinical decision-making.

## Release position

v0.6-beta is the first public source release after v0.5-beta. It combines the validated Multi-Pair Pre/Post foundation, multiaxial comparative views, QC-aware paired interpretation, and the H/V computational audit layer while preserving the established Main Mapper and single-pair contracts.

## Added

- One-to-six independent Pre/Post comparisons with shared C/S settings and pair-isolated failures.
- Three separate interpretation layers: C/S balance change, spatial-organization change, and specimen reliability.
- Optional H hypoxia-associated expression context and V endothelial/angiogenic expression proxy as parallel observational axes.
- `same_site`, `different_site`, and `unknown_site` metadata with a visible, non-excluding site-shift warning.
- Rule-based `Minimal`, `Moderate`, and `Large` descriptive change classes alongside `Good`, `Caution`, and `Low` specimen comparability.
- Auditable effective H/V gene programs, coverage, status, raw upper-tail summaries, within-pair pooled q90 high-context fractions, and local high-context fractions.
- A multiaxial overview with independently scaled C/S, spatial, raw H/V median Delta, and pooled high-context-fraction Delta panels. No composite response score is calculated.

## New or expanded outputs

- `balance_changes.csv`
- `spatial_organization_changes.csv`
- `specimen_reliability.csv`
- `pair_interpretation_summary.csv`
- `comparability_details.csv`
- `overview_interpretation.csv`
- `context_gene_audit.csv`
- `context_changes.csv`
- `multiaxial_pair_summary.csv`
- `comparative_qc_summary.csv`
- `figures/multiaxial_pair_overview.png`

Raw Pre, Post, Delta (`Post - Pre`), availability status, and safe percent-change fields remain available. H/V zero values are distinguished from unavailable values.

## Preserved scientific boundaries

- `C(x)`, `S(x)`, `R(x)=C(x)-S(x)`, and `G(x)` are unchanged.
- FRAME2.6 Type A/B/C candidate rules, interface/diffuse masks, transition burden, default thresholds, and single-sample Main Mapper behavior are unchanged.
- H/V do not modify C/S/R, transition detection, candidate regime labels, comparability QC, or site warnings.
- Comparative maps are descriptive side-by-side displays; no spatial registration or spot-wise subtraction is performed.
- The application does not infer treatment efficacy, responder status, drug sensitivity, survival, perfusion, diagnosis, or clinical benefit.

## Compatibility

Existing Main Mapper, Import / Convert, Single Pair Comparative Analysis, QUBO Optimizer, Advanced Analysis, Spatial Graph & Neighborhood, CLI entry points, and v0.5 output contracts remain available. The comparative cache schema is `v0.6-hv-validation-v1` so older cached metrics cannot be mistaken for audited v0.6 output.

## Validation

The complete automated test suite, Python compilation, package metadata, screenshot assets, and release archive are rechecked during packaging. The detailed historical computational audit is recorded in `VALIDATION.md` and `H_V_COMPUTATIONAL_AUDIT_VALIDATION_REPORT.md`.

## Screenshots

Four v0.6 screenshots are included under `docs/screenshots/`: the Multi-Pair workspace, multiaxial change profile, H/V context summary, and an H/V joint-state map. They were captured during final development testing and may display the internal `v0.6-dev-HV-validation` label; the distributed application identifies itself as `v0.6-beta`.
