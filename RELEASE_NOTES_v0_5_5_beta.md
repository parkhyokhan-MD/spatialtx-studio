# SpatialTX Studio Desktop v0.5.5-beta

Release date: 2026-08-10

SpatialTX Studio Desktop v0.5.5-beta is an incremental comparative-analysis upgrade for exploratory research use. It preserves the v0.5 Main Mapper, Import / Convert, Single Pair Comparative Analysis, C/S definitions, `R(x)=C(x)-S(x)`, Type A/B/C candidate rules, QUBO tools, Advanced Analysis, and spatial graph workflows.

## New in this release

- Added **Comparative Analysis → Multi-Pair Pre/Post** for one to six independent labeled pairs.
- Organized Multi-Pair interpretation into three separate layers: **Balance change**, **Spatial organization change**, and **Specimen reliability**. The application does not combine them into an overall response or quality score.
- Added a pair-level interpretation panel with configurable, auditable `Minimal/Moderate/Large` change labels, `regime_preserved`, `structure_preserved`, reliability-qualified status, exact rule basis, and prominent Low-comparability safety messages.
- Added conservative filename pair-ID checking. Explicit patient/sample mismatches warn and request confirmation but do not automatically block a run or override the comparability class; accessions alone are not interpreted as patient IDs.
- Added compact actual comparability ratios/differences and separated technical, sampling/geometry, and secondary composition-proxy reasons.
- Reused one canonical pairwise analysis engine and one shared immutable parameter set across all selected pairs.
- Added separate Pre, Post, Delta, safe percent-change status, and metric-aware direction outputs for C, S, R, localized interface-like fraction, diffuse fraction, transition burden, adjacency, fragmentation, and compatible topology metrics. Primary C/S/R rows use existing field medians; centered means remain separate compatibility columns.
- Added descriptive regime and Type B pattern transitions without clinical-response interpretation.
- Added a rule-based, auditable **Comparability Gate** with `Good`, `Caution`, and `Low` classifications.
- Added technical, sampling, occupancy, gene-coverage, and geometry QC; unavailable optional metrics remain visible as `not_available`.
- Labeled C/S distribution comparisons as secondary `composition proxy` context. They cannot alone produce `Low`.
- Added pair-isolated error handling so one invalid or corrupted pair does not stop the other pairs.
- Retained `pair_results.csv` for compatibility and added layer-specific tables, `pair_interpretation_summary.csv`, `comparability_details.csv`, and `overview_interpretation.csv`, plus detailed QC, cohort-count, metadata, and three-layer figure exports.
- Added **How to read results** and a **Rules & interpretation** tab so users can inspect the exact Good/Caution/Low rules, default thresholds, Delta direction, C/S/R summary basis, and interpretation limits without leaving the program.
- Corrected the side-by-side R-map lower layout so the outline legend no longer overlaps the run footer.
- Split the raw comparative overview into two independent-axis panels: **Primary spatial-state summary metrics** and **Topology / component complexity metrics**. Every plotted bar has a numeric Delta label, and large topology-density values no longer determine the scale of the fraction/score panel.
- Added a separate standardized-change overview when a group comparison has a valid pooled sample scale. Raw Delta remains the default export, and no pairwise z-score is fabricated.

## Output folder

Multi-Pair runs write to:

```text
<output_root>/comparative_multi_pair/<timestamp_or_tag>/
  pair_results.csv
  balance_changes.csv
  spatial_organization_changes.csv
  specimen_reliability.csv
  pair_interpretation_summary.csv
  comparability_details.csv
  overview_interpretation.csv
  comparative_overview.csv
  comparability_qc.csv
  cohort_summary.csv
  run_metadata.json
  figures/
```

Every Delta is `Post - Pre`. Direction arrows are display aids only; exact values remain in the CSV files. Percent change is `NA` for unstable zero or near-zero Pre values.

## Comparability boundary

Specimen reliability describes whether technical and sampling properties are reasonably aligned for cautious review. It is not a biological outcome or a probability that the observed change is true. A `Low` pair remains numerically visible but carries a prominent warning that tissue composition or sampling differences may substantially influence the observed change. Balance, spatial organization, and reliability remain separate outputs.

## Explicit exclusions

This release does not add or activate:

- new H/V/M/P biological axes;
- AI or machine-learning models;
- comparative QUBO or quantum backends;
- tissue registration or spot-wise subtraction;
- composition adjustment or multiple ROI per time point;
- treatment response rate, responder classification, therapeutic efficacy, drug sensitivity, or predictive biomarker accuracy.

## Research-use notice

Outputs are descriptive and exploratory. They are not intended for diagnosis, prognosis, treatment selection, or clinical decision-making and require independent scientific review and validation.
