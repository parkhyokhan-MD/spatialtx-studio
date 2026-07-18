# SpatialTX v0.4 output schema

Spatial Graph and Neighborhood Analysis writes one timestamped result folder:

```text
spatial_graph_neighborhood_<timestamp>/
  run_parameters.json
  spatial_graph_neighborhood_manifest.csv
  combined_cohort_summary.csv
  input_audit/
    input_audit_<sample>.json
    input_audit_<sample>.csv
  variable_semantics/
    variable_semantics_<sample>.json
    variable_semantics_<sample>.csv
  spatial_graph/
    <sample>_graph_metadata.json
    <sample>_graph_qc.csv
    <sample>_graph_degree_distribution.csv
  neighborhood/
    <sample>_categorical_enrichment.csv
    <sample>_binary_mask_association.csv
    <sample>_same_spot_overlap.csv
    <sample>_neighboring_spot_association.csv
    <sample>_continuous_edge_statistics.csv
    <sample>_permutation_parameters.json
  context_fields/
    <sample>_context_field_summary.csv
    <sample>_context_field_gene_coverage.csv
    <sample>_context_field_leave_one_gene_out.csv
    <sample>_H_expr_qc.json
    <sample>_V_expr_qc.json
  robustness/
    graph_robustness_summary_<sample>.csv
    association_direction_stability_<sample>.csv
  figures/
    <sample>_graph_qc.png
    <sample>_graph_qc.png.metadata.json
    <sample>_neighborhood_enrichment_heatmap.png
    <sample>_H_expr_map.png
    <sample>_H_expr_unsmoothed_map.png       # when smoothing is enabled
    <sample>_H_expr_smoothed_map.png         # when smoothing is enabled
    <sample>_V_expr_map.png
    <sample>_V_expr_unsmoothed_map.png       # when smoothing is enabled
    <sample>_V_expr_smoothed_map.png         # when smoothing is enabled
    <sample>_H_V_association_map.png
  annotated/
    <sample>_spatialtx_v0_4_graph_context.h5ad   # optional
```

Every run records the application version, analysis-module version, input-audit state, variable semantics, graph parameters and coordinate calibration, context-field parameters/QC, C/S genes, thresholds, permutation scope/count/seed, input filename, output path, and creation time. Association tables also record label provenance, exploratory interpretation, permutation limitation, and `fdr_scope=within_sample_within_analysis_table`.

Successful Main Mapper, CLI, Advanced Analysis, and Spatial Graph metadata include `gene_program_validation` with requested and used C/S lists, overlap genes/count, duplicate removal, policy, mode, action, warnings, and validation status. Overlap is a hard error, so a successful result has `n_overlap_genes=0` and `validation_status=valid`.

Existing Main Mapper output folders and filenames are unchanged when the v0.4 graph module is not enabled.

Multi-seed optimizer runs write a separate timestamped folder under `optimizer/`:

```text
<sample>_<side>_multiseed_<timestamp>/
  optimizer_multiseed_runs.csv
  optimizer_selection_frequency.csv
  optimizer_pairwise_overlap.csv
  optimizer_multiseed_gene_details.csv
  optimizer_stability_summary.csv
  optimizer_stability_summary.json
  optimizer_selection_frequency.png
  optimizer_energy_stability.png
```

The exact-k consensus is intended for review and is never applied automatically. The summary records that QUBO energy is lower-is-better and that multi-seed consensus measures computational stability only.

Single- and multi-seed optimizer summaries also record `overlap_constraint_enabled`, `genes_excluded_due_to_opposite_side`, `final_overlap_genes`, and `final_overlap_count`. The sequential implementation enforces `x_C,g + x_S,g <= 1` by removing opposite-side genes from the candidate pool before optimization and validating the final program.
