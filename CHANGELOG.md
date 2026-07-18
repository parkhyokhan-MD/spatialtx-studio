# Changelog

All notable public and development changes to SpatialTX Studio Desktop are recorded here.

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
