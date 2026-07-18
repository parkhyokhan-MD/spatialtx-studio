# SpatialTX Studio Desktop v0.4-beta stabilization report

Validation date: 2026-07-15  
Status: source-based beta candidate; no GitHub release or tag was created.

## Source provenance and revalidation

- User-supplied source archive: `SpatialTX_Studio_Desktop_v0_4_dev_original.zip`
- Source archive SHA-256: `7806B3065B90093DC2AE874ECD8364AC792001041CA8D98EA98F323F055E1A01`
- The archive contained 110 files. A relative-path and per-file SHA-256 comparison against the preserved `SpatialTX_Studio_Desktop_v0_4_dev_original/` copy found zero differences.
- On 2026-07-15, the preserved original again passed 49/49 tests and the beta work copy passed 63/63 tests without weakening or removing tests.
- The existing wheel was installed again into a fresh local virtual environment with `--system-site-packages --no-deps`; all three console entry points, the packaged default config, and a minimal synthetic H5AD FRAME2.6 run passed.

## 1. Modified files

Core and entry points:

- `spatialtx_desktop/workflow.py`
- `spatialtx_studio/frame26.py`
- `spatialtx_studio/legacy_frame26.py` (renamed audit-only legacy implementation)
- `spatialtx_studio/runner.py`
- `spatialtx_studio/io.py`
- `spatialtx_desktop/version.py`
- `spatialtx_desktop/__init__.py`
- `spatialtx_studio/metadata.py`
- `desktop_app.py`
- `advanced_cli.py`
- `run_desktop.bat`

Spatial Graph and context stabilization:

- `spatialtx_desktop/graph/audit.py`
- `spatialtx_desktop/graph/builder.py`
- `spatialtx_desktop/graph/context.py`
- `spatialtx_desktop/graph/continuous.py`
- `spatialtx_desktop/graph/enrichment.py`
- `spatialtx_desktop/graph/metadata.py`
- `spatialtx_desktop/graph/plotting.py`
- `spatialtx_desktop/graph/qc.py`
- `spatialtx_desktop/graph/runner.py`
- `spatialtx_desktop/advanced_analysis_ui.py`
- `spatialtx_desktop/app.py`

Packaging and tests:

- `setup.cfg`
- `config_default.yaml`
- `spatialtx_studio/resources/__init__.py`
- `spatialtx_studio/resources/config_default.yaml`
- `tests/test_beta_stabilization.py`
- `tests/test_advanced_analysis.py`

Documentation and release metadata:

- `README.md`, `README_DESKTOP.md`, `CHANGELOG.md`, `VALIDATION.md`
- `CITATION.cff`, `DEVELOPMENT.md`, `DISCLAIMER.md`, `REFERENCES.md`, `THIRD_PARTY_LICENSES.md`
- `RELEASE_NOTES_v0_4_beta.md` (replaces the beta work copy's former dev-notes filename)
- `docs/SPATIAL_GRAPH_NEIGHBORHOOD.md`
- `docs/INPUT_AUDIT_AND_VARIABLE_SEMANTICS.md`
- `docs/OUTPUT_SCHEMA_v0_4.md`
- three current v0.4-beta screenshots under `docs/screenshots/`

## 2. Engine unification

`spatialtx_desktop.workflow.score_adata()` is now the canonical in-memory C/S and FRAME2.6 engine. `score_h5ad()` performs one file read and delegates to it. Main Mapper continues to use this implementation, Spatial Graph passes its already loaded AnnData directly to it, and the public CLI `spatialtx_studio.frame26.run_frame26()` is a compatibility/export wrapper around it.

The former independent CLI calculation was moved to `spatialtx_studio/legacy_frame26.py`. It is retained for audit history, is not imported by the public runner, and no longer computes a second public result. Existing CLI filenames and legacy B-column aliases are retained, with B documented as the compatibility alias for canonical S.

## 3. Existing versus stabilized results

- Main Mapper canonical C, S, R, G, localized-interface mask, diffuse-transition mask, regime label, Type B/public pattern, transition burden, and adjacency/fragmentation metrics remain unchanged under identical input and settings.
- `score_h5ad()` and `score_adata()` return numerically identical canonical fields and classifications.
- The CLI result intentionally changes from the retired independent grid engine to the unchanged Main Mapper engine. This removes GUI/CLI disagreement rather than changing the Main Mapper definition.
- H_expr and V_expr remain descriptive context-expression fields and are not included in C, S, R, G, Type A/B/C, or transition-mask calculations.
- Spatial Graph remains opt-in and Experimental. Main Mapper output names are unchanged.

## 4. Test results

- Baseline before modification: 49/49 tests passed in the preserved original copy.
- Final beta work copy: 63/63 tests passed.
- New tests cover canonical engine equality, one-load delegation, low-density large graphs, duplicate-coordinate inverse-distance handling, zero-edge invalidation, sparse preservation, detection-scale semantics, counts-layer preference, smoothing warnings, derived-state provenance, FDR scope, permutation limitations, leave-one-gene-out QC, packaged config, deterministic plot-edge sampling, Korean/space paths, and all requested synthetic graph combinations.
- Three synthetic classes were tested: regular Visium-like lattice, irregular coordinates, and duplicate coordinates.
- Radius, Visium lattice/fallback, and symmetric KNN were tested for all three classes (nine combinations).
- Python compile checks passed.

## 5. Wheel installation validation

Final wheel:

```text
dist/spatialtx_studio_desktop-0.4b0-py3-none-any.whl
SHA-256: 2841002598C4D08EBDA7CEFEC84FB651238F9129AA9872866DB5674B801D4FBE
```

The wheel contains `spatialtx_studio/resources/config_default.yaml`. In a fresh virtual environment using the installed system scientific stack, the wheel was installed with `--no-deps` and verified from outside the source directory. `spatialtx --help`, `spatialtx-desktop --help`, `spatialtx-advanced --help`, packaged config loading, and a minimal synthetic H5AD FRAME2.6 run passed. The installed module and config paths resolved inside the virtual environment's `site-packages`, not the source tree.

## 6. Windows GUI and path validation

- `run_desktop.bat --help` passed and now forwards command-line arguments to the launcher.
- The source GUI opened on Windows, displayed `v0.4-beta`, Main Mapper, Import / Convert, Advanced Analysis, and `Spatial Graph & Neighborhood — Experimental`, and closed after screenshot capture.
- Main Mapper scanning and analysis passed with Korean and space-containing input/output paths.
- Spatial Graph analysis passed on the same path class.
- CSV, JSON, PNG, and optional annotated H5AD exports were verified.
- Windows plotting now prefers Malgun Gothic so Korean sample names do not generate missing-glyph warnings.
- The application is Tkinter/ttk, not PySide6; therefore PySide6-specific shutdown checks are not applicable. The capture harness destroyed the Tk root and no Python GUI process remained.

## 7. Remaining statistical limitations

- Graph association results are exploratory and do not establish causal, physical, or biological cell-cell interactions.
- Permutation P-values rely on exchangeability and do not fully preserve the original spatial autocorrelation structure.
- BH-FDR scope is `within_sample_within_analysis_table`; it is not a correction across samples, graph types, analysis families, or robustness runs.
- SpatialTX-derived masks are already spatially constructed and their graph results are descriptive organization summaries, not independent tests.
- Statistics using graph-smoothed H/V fields may be inflated because smoothing and association can use the same graph. They are sensitivity/visualization outputs, not confirmatory inference.
- Input preprocessing and platform states are heuristic guesses unless explicitly documented by the user.
- Coordinate distances are physical only when a verified micrometer scale and source are recorded.

## 8. Publicly presentable features

- H5AD-centered Main Mapper with strict spatial QC and expression-only fallback
- Raw 10x MEX/MTX and Raw Visium H5 + spatial conversion to H5AD
- C/S/R/G mapping and research-use-only reports
- Single-seed QUBO plus multi-seed computational stability summaries
- Input audit, parameter logging, memory preflight, and reproducible source/CLI operation
- Apache-2.0 source distribution with packaged default config

These features remain exploratory and research-use only.

## 9. Features that must remain Experimental

- Spatial Graph & Neighborhood and graph-definition robustness comparison
- H_expr/V_expr context association, especially smoothed-field statistics
- A3-A5 condition comparison, receptor-like filtering, and QUBO candidate handoff
- Neighborhood P-values/FDR and user-label edge association
- QUBO consensus as evidence of computational stability rather than biological validation

## 10. Unresolved or deferred issues

- Full dependency resolution in a completely isolated online environment was not repeated; wheel validation used `--no-deps` with a fresh environment exposing the installed scientific Python stack.
- `install_desktop.bat` was not allowed to redownload all dependencies during this stabilization run; its launcher selection was inspected and `run_desktop.bat` was executed.
- The optional legacy `istz` analysis was not executed because Scanpy was unavailable in the active Anaconda environment. Its import is now lazy so FRAME2.6 and help commands do not fail for that reason.
- No real large cohort or platform-spanning benchmark was run; memory and graph checks used sparse synthetic data.
- Cancellation during the middle of a long single-sample graph computation remains cooperative between samples rather than immediate.
- `CITATION.cff` intentionally omits `date-released` until the actual public release date is known.
- No Git repository exists in this workspace, so the requested small commit sequence could not be created. Work was performed in step-gated changes with tests after each stage.

## Preservation record

```text
Preserved original: SpatialTX_Studio_Desktop_v0_4_dev_original/
Beta stabilization work: SpatialTX_Studio_Desktop_v0_4_beta_work/
```

The public v0.3-beta package and the preserved original v0.4-dev copy were not modified.
