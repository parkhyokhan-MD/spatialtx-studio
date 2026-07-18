# SpatialTX Studio Desktop validation record

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
