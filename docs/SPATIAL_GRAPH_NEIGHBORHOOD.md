# Spatial Graph & Neighborhood — Experimental

SpatialTX Studio v0.4-beta includes an optional graph-based analysis module for spatial neighborhood statistics and descriptive context fields. The module is separate from the Main Mapper and remains Experimental.

These analyses provide exploratory spatial association and organization summaries. They do not establish causal, physical, or biological cell-cell interactions.

## Architectural guardrails

- AnnData `.h5ad` remains the canonical internal analysis format.
- The existing C/S framework is unchanged: `C(x)`, `S(x)`, `R(x)=C(x)-S(x)`, and `G(x)` keep their current definitions.
- Existing Type A/B/C classification and Type B internal pattern rules are unchanged.
- `H_expr(x)` and `V_expr(x)` are optional context fields only.
- `H_expr(x)` and `V_expr(x)` do not modify `R(x)`, Type A/B/C labels, or transition masks.
- `V_expr(x)` is an endothelial/angiogenic expression proxy, not true blood-vessel density, perfusion, vascularity, or functional blood supply.

## Graph choices

Radius graph is the default for general spatial analysis. It uses coordinates from `adata.obsm["spatial"]`, stores sparse connectivities and distances separately, and estimates a radius from nearest-neighbor spacing when the user does not provide one. It is described as a calibrated physical-radius graph only when a micrometer unit and scale source are explicitly available; otherwise it remains a native-coordinate radius graph.

Visium lattice graph uses `array_row` and `array_col` when valid Visium lattice metadata are available. If lattice detection fails, the module records `requested_method=lattice`, `effective_method=radius`, `fallback_used=true`, and the actual effective radius. The fallback graph is stored under the radius graph keys.

Symmetric KNN graph is available as a robustness option. It is not treated as the default biological neighborhood definition because KNN can create physically long edges across sparse or irregular tissue regions.

## AnnData storage

SpatialTX graph matrices are stored in `adata.obsp` with explicit names:

- `spatialtx_connectivities_radius`
- `spatialtx_distances_radius`
- `spatialtx_connectivities_knn`
- `spatialtx_distances_knn`
- `spatialtx_connectivities_lattice`
- `spatialtx_distances_lattice`

Graph metadata are stored in `adata.uns["spatialtx_graph"]`. The module does not overwrite third-party keys such as `spatial_connectivities`.

## Graph QC

Each graph records nodes, undirected edges, degree summary, isolated-spot fraction, connected components, largest-component ratio, nearest-neighbor and edge-distance summaries, long-edge fraction, graph density, coordinate validity, duplicate-coordinate counts, and zero-distance edges excluded by inverse-distance weighting.

Default configurable warnings use `isolated_spot_fraction > 0.10`, `largest_component_ratio < 0.80`, `median_degree < 2`, `number_of_edges == 0`, KNN long-edge fraction, and near-complete radius graphs. Graph density is exported for context but low density alone is not a warning or failure because normal local-graph density decreases as sample size grows.

Duplicate coordinates produce a QC warning. With inverse-distance weighting, zero-distance edges are excluded and their count is recorded in graph QC and metadata. If no usable edge remains, the graph is invalid rather than successful.

A graph with zero usable edges is not analyzed. The sample manifest records an error instructing the user to review coordinate scale or choose an automatic/larger radius. Statistics with insufficient valid edges are marked `not_interpretable` or `insufficient_edges`, not depletion or avoidance.

## Neighborhood enrichment

Categorical neighborhood enrichment tests whether graph-neighbor edges connect labels more or less often than expected under seeded permutation. Reported values include observed edge count, expected edge count, null standard deviation, z-score, empirical P-value, FDR, observed/expected ratio, valid spot count, and valid edge count.

The empirical P-value uses:

```text
p = (extreme_count + 1) / (n_permutations + 1)
```

P-values and FDR values should be interpreted as exploratory spatial association statistics, not proof of biological interaction or causality. Permutation P-values rely on exchangeability assumptions and do not fully preserve the original spatial autocorrelation structure. BH-FDR is applied within one sample and one analysis table (`within_sample_within_analysis_table`), not across samples, graph types, analysis families, or robustness runs.

## Overlap versus neighborhood association

Binary-mask analysis reports two distinct modes:

- same-spot overlap: both masks are true in the same spot;
- neighboring-spot association: one mask is true in a spot and the other is true in a graph neighbor.

Avoidance or depletion means the observed count is lower than the permutation expectation.

## Context fields

`H_expr` is a hypoxia-associated expression field from a configurable gene set. `V_expr` is an endothelial/angiogenic expression proxy from a configurable gene set. Both fields store requested genes, matched genes, missing genes, coverage fraction, score method, smoothing method, and high-state quantile threshold.

Low gene coverage generates warnings. High states are quantile-based and should not be interpreted as validated biological cutoffs.

When graph smoothing is enabled, the smoothed H/V field is the active field for context maps and high-state masks. Continuous exports include active compatibility rows plus separate unsmoothed and smoothed R-H/R-V sensitivity rows. Both base and smoothed columns remain available in an optional annotated H5AD for auditability.

Graph-smoothed context fields are intended for visualization and exploratory sensitivity analysis. Association statistics computed on fields smoothed over the same graph may be inflated and should not be interpreted as independent confirmatory evidence.

Context QC also records requested/matched/expressed genes, expression scale and detection source, field range, high-state fraction, experimental legacy dominant-gene contribution, leave-one-gene-out field influence, library-size and detected-gene correlations, C/S/R correlations, active smoothing graph, and explicit warnings. Detection fraction is reported only for verified raw/nonnegative log-like expression or a counts/raw layer; centered/scaled input reports detection as unavailable and provides positive-value fraction instead.

## Input audit and variable semantics

Each sample receives `spatialtx_input_audit` metadata and `input_audit_<sample>.json/.csv`. User-selected `adata.obs` columns are recorded as mutually exclusive categorical state, non-exclusive binary mask, continuous score, or proportion/composition value. Sequencing-based spots are not assumed to represent pure single-cell identities.

## Permutation scope

`whole_slide`, `within_connected_components`, and `within_user_strata` scopes are supported. Missing user strata are excluded, and optional tissue restriction uses `in_tissue=1`. Fixed seeds are reproducible and permutations never cross samples.

## Generic continuous association

User-selected X and Y columns use the same sparse symmetric weighted edge statistic as predefined C/S/R/H/V pairs. Missing observations and invalid edges are excluded and reported. The output includes the null mean/SD, z-score, corrected empirical P-value, BH-FDR, direction, valid observations/edges, graph metadata, and permutation metadata.

SpatialTX-derived masks use `label_source=spatialtx_derived` and `analysis_interpretation=descriptive_spatial_organization`. User annotations use `label_source=user_supplied` and `analysis_interpretation=exploratory_neighbor_association`.

## Graph robustness

Optional radius/lattice/KNN comparison reports association-direction and significance stability, observed/expected and z-score variation, and graph-QC variation. It is supplementary and never replaces the primary graph result.

## Outputs

The module writes structured sample-level outputs under a timestamped `spatial_graph_neighborhood_*` folder:

```text
spatial_graph/
neighborhood/
context_fields/
figures/
annotated/        # optional
```

Multiple samples are analyzed independently. Edges and permutations are never mixed across slides; cohort summaries are combined only after sample-level results are created.

Graph QC figures use all edges for statistics but cap displayed edges at `max_plot_edges=50000` by default using deterministic `plot_edge_seed=42` sampling. A sidecar `*.png.metadata.json` records displayed and total edge counts.
