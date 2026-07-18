# SpatialTX Studio Desktop v0.4-beta

Pre-release stabilization notes.

This source-based beta candidate starts from the public v0.3-beta source baseline. Preparing this working copy does not create a GitHub release or public tag.

Current baseline behavior remains exploratory and research-use only:

- Main Mapper analyzes canonical AnnData `.h5ad` inputs.
- Import / Convert prepares supported raw 10x/Visium inputs as `.h5ad` before analysis.
- Advanced Analysis and Advanced / Experimental remain optional hypothesis-generation areas.
- Outputs are not intended for diagnosis, treatment selection, or clinical decision-making.

The Main Mapper C/S engine, Type A/B/C rules, transition masks, and raw import formats remain unchanged unless explicitly noted below.

## C/S gene-program overlap safety fix

- Official analyses now require `C_genes ∩ S_genes = ∅`.
- Gene symbols are trimmed, uppercased, de-duplicated within each program, and checked before analysis and again in the common scoring engine.
- Fixed-program overlap is treated as a development error. Custom overlap is a user-correctable hard error and is never silently assigned to one side.
- Adaptive and QUBO selection exclude genes already selected on the opposite side and validate the final programs again.
- Main Mapper parameter logs, CLI/Advanced/Graph JSON metadata, Gene Composition tables, and QUBO summaries record validation provenance.
- The C, S, `R=C-S`, gradient, interface/diffuse masks, Type A/B/C rules, Type B patterns, and results for pre-existing non-overlapping programs are unchanged.

## Spatial Graph, Neighborhood Statistics, and Context Fields

The v0.4 beta line adds an optional **Spatial Graph & Neighborhood — Experimental** workflow under Advanced Analysis.

- Sparse native/calibrated-radius, Visium-lattice, and symmetric-KNN graph builders.
- Graph QC with isolated-spot, connected-component, degree, duplicate-coordinate, and long-edge warnings.
- Stable SpatialTX graph storage in `adata.obsp` and `adata.uns["spatialtx_graph"]`.
- Optional `H_expr` hypoxia-associated expression field.
- Optional `V_expr` endothelial/angiogenic expression proxy.
- Categorical neighborhood enrichment with empirical permutation P-values and FDR.
- Binary-mask same-spot overlap and neighboring-spot association.
- Continuous edge statistics for C/S/R/H/V fields.
- Structured output folders for graph metadata, neighborhood statistics, context-field coverage, and figures.

These analyses are exploratory spatial association summaries. They do not alter the established C/S fields, `R(x)`, Type A/B/C classifications, Type B patterns, or transition masks.

## Graph correctness and optimizer stability update

- Graph analysis now fails clearly when the selected parameters produce no usable spatial edges.
- Non-interpretable null statistics are no longer labeled as depletion or avoidance.
- Visium lattice fallback records both the requested lattice method and the effective radius method, including the actual radius used.
- Smoothed H/V fields are used consistently by maps, high-state masks, and continuous edge statistics.
- QUBO Optimizer now provides separate single-seed and multi-seed stability runs.
- Multi-seed outputs include seed-level selections and energies, selection frequency, consensus core genes, exact-k consensus genes, pairwise overlap, R-field agreement, interface-mask agreement, regime agreement, and stability figures.
- QUBO energy is explicitly labeled as lower-is-better for otherwise identical optimization problems.
- Multi-seed consensus is a computational stability diagnostic and is not biological validation or a uniquely optimal gene-program claim.

## Input audit and spatial-association refinement

- Added a non-mutating input audit stored in `adata.uns["spatialtx_input_audit"]` and exported per sample.
- Added explicit native versus calibrated coordinate terminology and radius provenance.
- Added four user-variable semantics: spot-level categorical state, binary mask, continuous score, and proportion/composition.
- Added generic user-selected X–Y continuous edge association.
- Added whole-slide, connected-component-aware, and user-stratified permutations.
- Expanded H_expr/V_expr quality control and warnings without changing R, G, transition masks, or regime labels.
- Separated same-spot overlap, neighboring-spot association, and continuous edge association outputs.
- Added optional graph-definition robustness summaries.
- These features remain lightweight, sparse, exploratory, and non-clinical. No quantum SDK, GPU framework, deconvolution, or mandatory Squidpy dependency was added.

## Beta stabilization

- GUI/Main Mapper, CLI FRAME2.6, and Spatial Graph now call the same canonical in-memory C/S scoring engine.
- `score_h5ad()` loads once and delegates to `score_adata()`; Spatial Graph no longer reloads the same H5AD for C/S scoring.
- The retired independent CLI engine is retained only as `legacy_frame26.py` for audit history and is not called by the public CLI.
- Low graph density alone no longer produces an almost-empty warning; configurable degree, isolation, and component criteria drive QC.
- Inverse-distance weighting excludes zero-distance edges and records the excluded count; a graph with no remaining edge is invalid.
- Derived SpatialTX masks are labeled as descriptive spatial-organization summaries rather than independent interaction tests.
- Smoothing, permutation exchangeability, and FDR-scope limitations are recorded in UI text, CSV tables, and JSON metadata.
- Detection fraction is reported only for verified nonnegative raw/log-like expression or a counts/raw layer; centered/scaled inputs receive positive-value fraction instead.
- H/V context QC now includes leave-one-gene-out influence summaries.
- The packaged wheel includes `spatialtx_studio.resources/config_default.yaml` and loads it through `importlib.resources`.
- Graph plotting uses `LineCollection` and deterministic display-only edge downsampling; all statistical calculations continue to use the full graph.
