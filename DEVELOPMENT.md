# SpatialTX Studio Desktop developer notes

Current release line: v0.5-beta

Development start date: 2026-08-04

The v0.5-beta source release builds on the v0.4 analysis line while keeping the single-sample scoring contracts unchanged.

## v0.5-beta scope

- Adds a top-level **Comparative Analysis** tab and `spatialtx compare` / `spatialtx-compare` CLI entry points.
- Adds a third Import / Convert input mode for user-selected prefix-grouped GEO flat standard Visium directories without changing the existing folder importers.
- Supports pairwise, paired-group, unpaired-group, and manifest-batch sample-level comparison.
- Calls the existing canonical `score_adata()` engine for every sample instead of duplicating C/S/R/G calculations.
- Keeps `.h5ad` as the canonical analysis input and leaves raw conversion in Import / Convert.
- Uses sample-level summaries by default; no direct coordinate registration or spot-wise subtraction is implemented.
- Keeps `H_expr` and `V_expr` observational only. They cannot modify C/S/R, transition masks, transition burden, or Type A/B/C candidate labels.
- Adds a centralized comparative metric registry and a post-scoring normalization/display layer. This layer preserves canonical raw values and does not modify scoring, graph construction, thresholds, masks, or regimes.
- Treats within-sample centered H/V means as deprecated comparative summaries while preserving the existing spot-level H/V fields and maps.
- Writes new results only to timestamped `comparative_analysis_*` folders and a content-addressed cache under the selected output root.
- Preserves current Main Mapper, Import / Convert, QUBO, Advanced Analysis, Spatial Graph, viewer, and export behavior.

## Explicit exclusions

The comparative module does not add candidate discovery, ligand-receptor analysis, QUBO optimization, AI interpretation, literature search, or multi-axis modeling. Existing tools elsewhere in the application remain present for v0.4 compatibility but are not invoked by Comparative Analysis.

## Development guardrails

- Treat every comparative output as exploratory and non-diagnostic.
- Label delta direction as Target minus Reference.
- Keep invalid and failed samples visible in the run manifest and warnings.
- Do not infer direct spatial correspondence across samples without a future validated registration workflow.
- Do not call operational regime changes validated biological state transitions.
- Keep operational regime changes and candidate genes explicitly exploratory; public-data demonstrations do not constitute biological validation.
