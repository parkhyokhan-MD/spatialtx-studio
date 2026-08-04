# Changelog

All notable public and development changes to SpatialTX Studio Desktop are recorded here.

## v0.5-beta

- Improved pairwise Comparative Analysis sample selection for long directories: filename-first Sample A/B labels, duplicate-name disambiguation, and separate read-only full-path fields now keep the selected H5AD files identifiable in narrow windows.
- Removed the filename-based A1 pre/post pair scanner from the Advanced / Experimental UI while retaining its backend for compatibility; renamed A3 as a non-spatial Reference/Target expression candidate contrast and added direct H5AD browse controls.
- Redesigned Comparative Analysis visualization around a centralized ordered metric registry and category-specific figures so incompatible raw counts, fractions, scores, and ratios no longer share one primary axis.
- Added sample-scale exports and warnings for valid/in-tissue spots, tissue area proxy, tissue graph components, coordinate extent, and nearest-spot spacing.
- Preserved raw topology counts and added per-spot, per-in-tissue-spot, per-tissue-component, transition-spot, and normalized fragmentation metrics.
- Added ordinary and symmetric percent changes with zero/near-zero reference protection; group-only standardized deltas use pooled sample scale, while pairwise z-scores are not invented.
- Replaced the sparse pairwise 3 by 3 regime heatmap with a descriptive transition card; paired mode retains the count and row-percentage matrix.
- Added explicit R-map fill/outline legends, a category-grouped standardized heatmap, sample-scale figure, relative-change figure, and compact exploratory summary card.
- Deprecated centered H/V sample means for comparative interpretation while retaining them for compatibility; added pooled-threshold non-centered distribution, transition-enrichment, and variability summaries when supported by the input expression scale.
- Added six new comparative CSV exports while preserving previous filenames and raw values.

- Added a path-agnostic **GEO Flat Visium Directory** importer mode alongside the two existing standard-folder modes.
- Added exact full-prefix grouping, deterministic inventory, strict duplicate/dimension/orientation/position/scalefactor validation, gzip support, arbitrary input/output paths, read-only source handling, privacy-conscious H5AD provenance, portable output-relative image references, and collision-safe batch conversion.
- Added editable, explicitly confirmed comparative-manifest handoff; unconfirmed filename inference never creates active pair/group assignments.
- Added `spatialtx import-geo-flat` and `spatialtx-import-geo-flat` CLI entry points plus inventory, validation, log, and conversion-summary reports.
- Added a separate Comparative Analysis workspace for pairwise, paired-group, unpaired-group, and manifest-based comparisons.
- Reused the canonical single-sample C/S scoring engine and sample-level metrics without changing Main Mapper behavior.
- Added clearly directed `Target - Reference` deltas, robust group summaries, effect sizes, seeded bootstrap confidence intervals, Wilcoxon/Mann-Whitney defaults, optional explicit t tests, and Benjamini-Hochberg FDR.
- Added descriptive, confidence-flagged Type A/B/C candidate regime transition/distribution summaries.
- Added fresh comparative figures, cautious rules-based summaries, HTML/PDF reports, CSV/JSON exports, input hashes, software-environment metadata, warnings, and failed-sample manifests.
- Added memory preflight, content-addressed sample-summary caching, large-batch warnings, responsive GUI execution, and safe partial export on cancellation.
- Kept H_expr and V_expr optional and observational only; they do not modify C/S/R, masks, transition burden, or regime labels.
- No spatial registration or spot-wise subtraction is performed.
- Candidate discovery, ligand-receptor analysis, comparative QUBO, AI interpretation, literature search, and multi-axis modeling were not added.
- Outputs remain exploratory and are not intended for diagnosis, treatment selection, or clinical decision-making.

## v0.4.1-beta

- Added a read-only in-app Spatial Graph Results viewer for existing generated PNG figures.
- Added automatic post-run navigation, sample/figure selection, previous/next navigation, and fit-to-window image display.
- Added H_expr/V_expr context QC summaries, failed/cancelled sample summaries, and direct access to existing result tables.
- Updated the application version from v0.4-beta to v0.4.1-beta.
- Preserved Main Mapper scoring, C/S/R/G definitions, Type A/B/C rules, H_expr/V_expr calculations, graph construction, permutation statistics, output schema v0.4, existing filenames, and CLI behavior.
- Outputs remain exploratory and are not intended for diagnosis, treatment selection, or clinical decision-making.

## v0.4-beta

- Enforced mutually exclusive C-side and S-side gene programs across Main Mapper, CLI, Advanced Analysis, Spatial Graph, adaptive selection, and QUBO follow-up analysis.
- Added canonical gene normalization (trim, uppercase, empty removal, order-preserving within-program de-duplication) and hard-error overlap validation in both UI preflight and the common scoring engine.
- Added adaptive/QUBO opposite-side candidate exclusion, final overlap revalidation, and machine-readable validation/exclusion provenance.
- Added Gene Composition validation metadata and prevented the same gene from being emitted under both programs.
- Prepared a source-based v0.4-beta stabilization copy from the public v0.3-beta baseline while preserving the original v0.4-dev working copy.
- Preserved the Main Mapper / Import-Convert / Advanced Analysis separation introduced in v0.3-beta.
- Kept Main Mapper centered on AnnData `.h5ad` input.
- Kept raw-format handling in Import / Convert before analysis.
- Added optional Spatial Graph and Neighborhood Analysis under Advanced Analysis.
- Added reusable sparse spatial graph engine with physical radius, Visium lattice, and symmetric KNN graph builders.
- Added graph QC metrics, graph metadata, and stable AnnData `obsp`/`uns` storage for SpatialTX graph outputs.
- Added optional `H_expr` hypoxia-associated expression context field and `V_expr` endothelial/angiogenic expression proxy.
- Added categorical neighborhood enrichment, binary-mask association, and continuous edge-interaction statistics with seeded permutation nulls.
- Added structured v0.4 graph/neighborhood/context/figure exports and documentation.
- Rejects graph runs that contain no usable spatial edges instead of exporting misleading depletion labels.
- Records requested and effective graph methods plus the effective fallback radius when Visium lattice detection falls back to a radius graph.
- Uses the same active H/V context field for maps, high-state masks, and continuous edge statistics when graph smoothing is enabled.
- Added multi-seed QUBO stability analysis with selection frequency, consensus core genes, deterministic exact-k consensus, pairwise overlap, R-field agreement, interface-mask agreement, regime agreement, and objective-stability exports.
- Changed the development optimizer iteration default from 300 to 1,000 while preserving explicit user control and the reproducible single-seed workflow.
- Added Python launcher (`py -3`) fallback to the Windows install and run batch files.
- Added non-mutating H5AD input audit with preprocessing/platform guesses, matrix/coordinate/tissue checks, library summaries, warnings, AnnData metadata, and sample-specific JSON/CSV exports.
- Distinguished native-coordinate and calibrated physical-radius graphs; added coordinate scale, scale source, physical-calibration status, requested/effective radius, radius unit, and platform provenance.
- Added explicit categorical-state, binary-mask, continuous-score, and proportion/composition semantics for user-selected `adata.obs` variables.
- Added generic sparse symmetric X–Y continuous edge association.
- Added whole-slide, connected-component-aware, and user-stratified permutation scopes plus optional tissue-only restriction.
- Expanded H_expr/V_expr QC with detection, dynamic range, dominance, library-size, C/S/R correlation, smoothing-graph, and high-state diagnostics.
- Split same-spot overlap and neighboring-spot association into separate exports while retaining the combined compatibility table.
- Added optional radius/lattice/KNN graph robustness and association-direction stability exports.
- Unified GUI/Main Mapper, CLI FRAME2.6, and Spatial Graph C/S calculations through the canonical `score_adata()` engine.
- Removed duplicate H5AD loading from the Spatial Graph sample workflow.
- Replaced low-density-only graph warnings with configurable degree, isolation, component, and distance QC.
- Excluded zero-distance edges under inverse-distance weighting and recorded the exclusion count.
- Added derived-state provenance, exploratory interpretation, smoothing limitations, permutation limitations, and explicit within-table FDR scope to outputs.
- Added expression-scale-aware context detection metrics and H/V leave-one-gene-out influence QC.
- Bundled the default YAML config as a wheel package resource.
- Added `LineCollection` graph plotting with deterministic display-only edge downsampling.
- No GitHub release or public tag has been created by this stabilization work.
- Outputs remain exploratory and are not intended for diagnosis, treatment selection, or clinical decision-making.

## v0.3-beta

- Added Import / Convert workflow.
- Moved raw input conversion out of Advanced / Experimental.
- Added Raw 10x MEX/MTX → H5AD conversion.
- Added Raw Visium H5 + spatial → H5AD conversion.
- Added GEO-style duplicate-equivalent handling for uncompressed/`.gz` Visium spatial files, preferring uncompressed files with warning-level status logging.
- Added converted H5AD validation workflow.
- Kept Main Mapper centered on AnnData .h5ad input.
- Added clearer separation between file preparation, core mapping, and experimental analysis tools.
- Added strict spatial-coordinate QC with expression-only fallback and no invented spatial geometry.
- Suppressed Type A/B/C regimes, localized interface-like candidates, transition metrics, and maps when spatial QC is incomplete.
- Added optional smoothing, normalization, threshold perturbation, parameter logging, and memory preflight diagnostics with conservative defaults.
- Outputs remain exploratory and are not intended for diagnosis, treatment selection, or clinical decision-making.

## 0.2-beta — 2026-07-05

- Preserved the v0.1-beta Transition Mapper workflow, Cx/Sx definitions, defaults, and output contracts.
- Added an **Advanced Analysis** workspace with Gene Composition, Interface Enrichment, and Cx/Sx Interaction tabs.
- Added per-gene relative-contribution tables and manuscript-ready bar charts.
- Added interface/non-interface composition, fold enrichment, Hedges' g, Mann-Whitney testing, and BH-FDR.
- Added neighborhood-based coexistence, antagonism, balance, weighted spatial overlap, edge mixing, and seeded permutation inference.
- Added 300-dpi PNG, vector PDF, CSV, run-manifest, and JSON provenance outputs.
- Added a separate `advanced_cli.py` entry point; the original `app_cli.py` behavior is unchanged.
- Added an in-app Results Dashboard with a combined module/sample summary, full CSV table viewer, direct CSV/figure actions, and a one-click run-all workflow.
- Expanded Theory & Metrics with the rationale, formulas, statistical assumptions, spatial-permutation reference, and interpretation limitations for all three Advanced Analysis modules.

## 0.1-beta — 2026-07-02

- First public source release.
- Added local desktop workflows for `.h5ad` discovery, spatial scoring, map generation, interpretation, and export.
- Added local C-side and S-side gene-program optimization with a classical simulated-annealing fallback.
- Added opt-in exploratory utilities for MEX conversion, condition comparison, heuristic candidate filtering, and QUBO candidate-pool handoff.
- Added QC flags, gene-coverage reporting, coordinate fallbacks, and research-use guardrails.
- Added source installation guidance, release documentation, citation metadata, and licensing notices.
