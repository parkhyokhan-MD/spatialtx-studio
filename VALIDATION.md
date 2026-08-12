# SpatialTX Studio Desktop validation record

## v0.6-beta release packaging

Validation date: 2026-08-12

- Promoted the validated v0.6 multiaxial and H/V audit line to the public `v0.6-beta` version without changing the analysis algorithms or output schemas.
- Python compilation passed for `spatialtx_desktop`, `spatialtx_studio`, `tests`, and `tools` under Python 3.12.13.
- Full release suite: **171 tests passed plus 26 subtests**; 20 existing AnnData implicit-index and SciPy Matrix Market transition warnings were reported.
- Updated two duplicate-coordinate test fixtures to copy AnnData 0.13 read-only coordinate arrays before mutation; application behavior is unchanged.
- Verified the four v0.6 screenshot assets and their README/desktop-guide links.
- Package version is `0.6b0`; application version is `0.6-beta`.

## v0.6-dev-HV-validation

Validation date: 2026-08-11

- Added the H/V computational and audit validation layer in a separate working copy; the v0.6-dev and v0.5.5 checkpoints remain unchanged.
- Full suite: **169 tests passed plus 26 subtests**. The 20 warnings are the existing AnnData implicit-index and SciPy Matrix Market transition warnings.
- Reused the existing default H/V programs and existing pooled-threshold implementation.
- Actual six-pair regression matched v0.5.5 exactly for C/S/R, regimes, interface, diffuse, burden, adjacency, fragmentation, and comparability.
- Regimes remained 5/6 A-to-A and 1/6 A-to-B; comparability remained five Low and one Caution.
- Pair 5 (`sample_43`) had complete 12/12 H and V coverage in both samples and finite median, q90, shared-threshold high fraction, and local high-context fraction outputs.
- Verified an actual zero-median V case (`sample_30`): both medians were zero while q90/high/local context summaries retained focal signal.
- Effective programs, q90 scope, raw method, and cache schema are recorded in `run_metadata.json`; per-sample coverage/status is in `context_gene_audit.csv`.
- See `H_V_COMPUTATIONAL_AUDIT_VALIDATION_REPORT.md` for formulas and exact review values.

## v0.6-dev Multiaxial Comparative Analysis

Validation date: 2026-08-11

- Preserved the v0.5.5 source checkpoint separately and implemented v0.6-dev in a distinct working folder.
- Python compilation passed for `spatialtx_desktop`, `spatialtx_studio`, `tests`, and development validation tools.
- Full suite: **166 tests passed plus 26 subtests**. Existing importer, Main Mapper, Single Pair, GeoFlat, graph, optimizer, Advanced Analysis, spatial-QC, and v0.5.5 Multi-Pair tests remain passing.
- Added focused coverage for available H/V, missing H/V, unchanged C/S/FRAME2.6 outputs with H/V enabled, raw multiaxial exports, different-site warning behavior, Low-QC interpretation, and A-to-A interface-down/diffuse-up redistribution wording.
- Re-ran the actual six-pair GSE316402/BTC working set using the preserved C/S genes, thresholds, scoring options, graph settings, and seed. v0.5.5 and v0.6 produced an exact CSV-string match for regime labels, interface fractions, diffuse fractions, transition burden, and comparability.
- The six-pair result remained **5/6 Type_A_candidate to Type_A_candidate** and **1/6 Type_A_candidate to Type_B_candidate** (`sample_30`). Comparability remained five `Low` and one `Caution` (`sample_41`).
- Enabled existing H/V programs in the v0.6 run and verified that H/V did not alter any listed core result. H/V raw Pre/Post/Delta values were exported independently.
- Verified `multiaxial_pair_summary.csv`, `comparative_qc_summary.csv`, `context_changes.csv`, and `figures/multiaxial_pair_overview.png` on the actual six-pair run.
- Visually reviewed the multiaxial figure. C/S, spatial organization, and H/V use independent axes centered on zero, exact numeric labels, extra label margins, and no composite score.
- No treatment-response, efficacy, responder, survival, perfusion, or clinical-decision inference is generated.

## v0.5.5-beta Multi-Pair Pre/Post

Validation date: 2026-08-11

- Added automated coverage for one-pair and six-pair runs, non-contiguous populated UI rows, the six-pair limit and seven-pair rejection, pair-isolated corrupted input, intentional specimen mismatch, exports, metadata, direction tolerances, and composition-proxy classification limits.
- Full suite: **162 tests passed** on 2026-08-11. Existing importer, Main Mapper, Single Pair, graph, optimizer, Advanced Analysis, and spatial-QC tests remain passing.
- Verified that technical and sampling mismatches produce auditable `Caution` or `Low` reasons and that a C/S composition proxy alone cannot produce `Low`.
- Verified separate Pre, Post, Delta, percent-change, and direction exports for C, S, R, interface, diffuse, transition burden, adjacency, and compatible topology metrics.
- Verified the three-layer UI and output contract: `balance_changes.csv`, `spatial_organization_changes.csv`, and `specimen_reliability.csv` remain separate while `pair_results.csv` is retained for compatibility. The overview contains aligned layer summaries but no composite response or quality score.
- Visually reviewed the generated three-layer overview and pair-level figure for readable layer labels, reliability warning placement, and non-overlapping footer text.
- Re-ran the existing five-pair GSE316402/BTC working set with its preserved v0.5.5 settings: 5/5 pairs passed; all five retained `Type_A_candidate → Type_A_candidate`; four were `Low` and one `Caution`, matching the prior run. All filename IDs matched on the explicit `sample_N` token despite distinct GSM accessions.
- Verified pair-level qualitative classes, exact threshold bases, regime/structure preservation, compact comparability ratios, technical-versus-sampling reason separation, Low-comparability cautions, and all three new interpretation exports.
- Added regression coverage for the required `balance shift with preserved structure, low comparability` case, conservative pair-ID matching/mismatch warnings, accession-only non-inference, and invalid interpretation-threshold rejection.
- Verified that Single Pair Comparative Analysis, Main Mapper, importers, graph tools, Advanced Analysis, C/S definitions, and Type A/B/C rules remain unchanged.
- Multi-Pair outputs are descriptive and do not calculate treatment response rate, therapeutic efficacy, responder status, drug sensitivity, or clinical prediction.

## v0.5-beta Comparative Analysis

Validation date: 2026-08-04

- Prepared a clean public source release separately from the internal development working copy.
- Added comparative validation, matching, metrics, statistics, plotting, reporting, shared runner, CLI, and GUI components.
- Added synthetic GEO-flat importer coverage for arbitrary paths, spaces/non-ASCII names, gzip, exact prefix grouping, required components, duplicates, matrix orientation/dimensions, position/scalefactor validation, source immutability, output collision, reporting, CLI, GUI mounting, and user-confirmed manifest handoff.
- Added Box 3 comparative coverage for sample scale, normalized topology, zero/near-zero percent change, symmetric percent change, raw-versus-normalized warnings, centralized categories, pooled H/V thresholds, robust H/V enrichment, centered-mean deprecation, pairwise regime card, paired matrix preservation, figure metadata/legends, GUI warnings, raw-delta preservation, and prohibited-claim checks.
- Current full suite: **149 passed**. Existing importer, Main Mapper, spatial QC, graph, optimizer, Advanced, and comparative regression tests remain passing.
- Added synthetic tests for delta direction, regime summaries, paired matching/incomplete-pair rejection, unpaired statistics/effect size/FDR, missing H/V, missing spatial coordinates, no-registration behavior, deterministic output, cache reuse, small-n warnings, mixed-failure retention, exports, and unchanged single-sample results.
- A public GEO-style flat Visium import and a pairwise comparative workflow were reviewed as functional demonstrations. These examples do not establish biological validity or clinical utility.
- The public archive is checked for forward-slash paths and excludes H5AD data, raw matrices, local results, caches, logs, wheels, nested archives, and repository metadata.

## v0.4.1-beta Spatial Graph Results patch

Validation date: 2026-07-21

- Added focused tests for figure suffix classification, longest-prefix sample matching, figure priority/fallback, manifest success/failure handling, context-summary edge cases, image fitting/error handling, and hidden-Tk viewer behavior.
- Preserved the Spatial Graph analysis module version and output schema at `0.4`; the new viewer reads existing output files only.
- Full regression-test results for this patch are recorded in `V0_4_1_BETA_SPATIAL_GRAPH_VIEWER_REPORT.md`.

## v0.4-beta stabilization

Validation date: 2026-07-14

- Preserved the public v0.3-beta package and copied v0.4 development into separate `SpatialTX_Studio_Desktop_v0_4_dev_original` and `SpatialTX_Studio_Desktop_v0_4_beta_work` directories.
- Python compile checks passed for desktop, studio, graph, importer, CLI, and test modules.
- All 63 automated tests passed: the original 49 plus 14 beta-stabilization regression tests.
- The canonical `score_adata()` engine produces identical C, S, R, G, interface mask, diffuse mask, regime, Type B pattern, transition burden, and adjacency metrics through `score_h5ad()` and the CLI FRAME2.6 wrapper.
- Spatial Graph uses its already loaded AnnData and no longer reloads the same H5AD for C/S scoring.
- Large local synthetic graphs are not warned solely because density is below 0.002. Degree, isolation, component, and distance QC drive warnings.
- Duplicate-coordinate inverse-distance tests exclude zero-distance edges, record the exclusion count, retain finite weights, and invalidate graphs with no remaining edges.
- Regular Visium-like, irregular-coordinate, and duplicate-coordinate synthetic datasets passed radius, lattice/fallback, and symmetric-KNN construction (nine graph/dataset combinations).
- Smoothing metadata/warnings, SpatialTX-derived label provenance, permutation limitations, within-table FDR scope, context scale-aware detection, counts-layer preference, leave-one-gene-out influence, zero-edge failure, and sparse-input preservation are covered.
- Main Mapper and Spatial Graph succeeded with Korean and space-containing paths and produced CSV, JSON, PNG, and optional annotated H5AD outputs.
- `run_desktop.bat --help`, source GUI startup/capture/clean shutdown, Main Mapper backend, H5AD scanning, and Spatial Graph backend passed on Windows. No Python GUI process remained after the capture harness closed.
- A wheel built as `spatialtx_studio_desktop-0.4b0-py3-none-any.whl`; its packaged `spatialtx_studio/resources/config_default.yaml` loaded through `importlib.resources` in a fresh virtual environment.
- Installed-wheel `spatialtx --help`, `spatialtx-desktop --help`, `spatialtx-advanced --help`, and a minimal synthetic FRAME2.6 H5AD run passed from outside the source directory.
- Mandatory Main Mapper regression confirms unchanged C, S, R, G, localized-interface mask, diffuse-transition mask, regime label, and Type B/public transition pattern when the optional graph/context workflow is not enabled.
- Existing v0.3-beta validation history is retained below for release-line traceability.

## v0.3-beta Import / Convert architecture

Validation date: 2026-07-08

- Twenty-two unit tests passed, including gzipped MEX/MTX conversion, prefixed/gzipped Visium conversion, GEO-style duplicate-equivalent Visium spatial-file selection, H5AD validation, strict missing/malformed/non-finite spatial QC, expression-only reporting, robustness diagnostics, memory preflight, and all existing Advanced Analysis tests.
- The actual GEO-style `GSM9532669_YUBOISE_*` Visium sample converted successfully to a validated 476-spot by 18,085-gene H5AD during development testing.
- MEX conversion code is located under `spatialtx_desktop/importers/` and no longer exists in `advanced.py` or `advanced_ui.py`.
- Advanced UI contains hypothesis-generation and review utilities only; it displays a note directing raw imports to Import / Convert.
- Main Mapper discovery remains restricted to `.h5ad` files.
- UI import smoke checks passed for the desktop app, Advanced tools, and unified Import / Convert panel.
- Robustness and memory-safety defaults preserve the existing Main Mapper workflow: smoothing is off, normalization is raw mean, perturbation checking is off, and parameter log export is on.

## v0.2-beta validation record

Validation date: 2026-07-05

## Automated checks

- All Python source files compiled successfully with Python 3.12.
- Six focused unit tests passed: v0.1 defaults/signature preservation, gene contribution accounting, missing-gene reporting, BH-FDR behavior, spatial metric behavior, seeded reproducibility, and dashboard summaries.
- A synthetic 144-spot `.h5ad` dataset completed Gene Composition, Interface Enrichment, and Cx/Sx Interaction runs.
- The synthetic enrichment run contained 16 interface and 128 non-interface observations and produced finite fold enrichment, Hedges' g, Mann-Whitney p-values, and BH-FDR values.
- The integration run generated CSV, JSON, 300-dpi PNG, and vector PDF artifacts for all modules.
- Generated interaction spot tables retained source observation identifiers.
- A hidden-window Tkinter smoke test loaded all three module manifests into the Results Dashboard, produced three summary rows, and populated the selected detailed table.

Run the focused tests from the source directory with:

```text
python -m unittest discover -s tests -v
```

## v0.1 compatibility audit

The following files are byte-for-byte identical to the supplied v0.1-beta archive:

- `app_cli.py`
- `config_default.yaml`
- `spatialtx_desktop/workflow.py`
- `spatialtx_studio/runner.py`
- `spatialtx_studio/frame26.py`
- `spatialtx_studio/gene_program.py`
- `spatialtx_studio/interface_detection.py`
- `spatialtx_studio/transition_metrics.py`
- `spatialtx_studio/transition_zone.py`

The desktop `app.py` was extended only to mount the new top-level Advanced Analysis tab and update displayed release metadata. Existing tabs and command handlers remain present.

The original `app_cli.py` uses the legacy requirements set, including Scanpy; the new `advanced_cli.py` uses the desktop requirements set. This dependency distinction is inherited from v0.1 and does not change the original CLI.
