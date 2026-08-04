# SpatialTX Studio Desktop v0.4.1-beta Spatial Graph viewer report

## 1. Version update

- Previous application version: `v0.4-beta`
- Updated application version: `v0.4.1-beta`
- Application version sources and visible labels were updated in `spatialtx_desktop/version.py`, `spatialtx_desktop/__init__.py`, `spatialtx_studio/metadata.py`, `spatialtx_desktop/app.py`, `desktop_app.py`, `advanced_cli.py`, package metadata through `setup.cfg`, citation metadata, tests, and current documentation.
- `spatialtx_desktop.graph.runner.MODULE_VERSION` remains `0.4` because it describes the unchanged analysis/output schema rather than the desktop application patch.
- `docs/OUTPUT_SCHEMA_v0_4.md`, output directories, result filenames, and `<sample>_spatialtx_v0_4_graph_context.h5ad` remain unchanged.

## 2. Changed files

New files:

- `spatialtx_desktop/spatial_graph_results_ui.py` — read-only result discovery, parsing, image display, context QC, and file access UI.
- `tests/test_spatial_graph_results_ui.py` — focused parser, image, version, hidden-Tk, and graph-completion integration tests.
- `RELEASE_NOTES_v0_4_1_beta.md` — patch release notes.
- `V0_4_1_BETA_SPATIAL_GRAPH_VIEWER_REPORT.md` — implementation and validation record.

UI and version code updated:

- `spatialtx_desktop/advanced_analysis_ui.py` — adds the new tab and loads/selects it after `graph_done`.
- `spatialtx_desktop/app.py` — About release description and patch edition.
- `spatialtx_desktop/version.py`, `spatialtx_desktop/__init__.py`, `spatialtx_studio/metadata.py` — application version and build date.
- `desktop_app.py`, `advanced_cli.py` — version-aware help/description strings; options and calculations are unchanged.
- `tests/test_advanced_analysis.py` — current-version assertion.

Documentation and metadata updated:

- `README.md`, `README_DESKTOP.md`, `CHANGELOG.md`, `CITATION.cff`, `DISCLAIMER.md`, `THIRD_PARTY_LICENSES.md`, `REFERENCES.md`, `DEVELOPMENT.md`, `VALIDATION.md`, and `docs/SPATIAL_GRAPH_NEIGHBORHOOD.md`.
- The historical `RELEASE_NOTES_v0_4_beta.md` is retained unchanged.

## 3. New UI structure

`SpatialGraphResultsPanel(ttk.Frame)` is a separate UI class. Its `load_run(run_dir, manifest)` method replaces only the Spatial Graph viewer state for a completed run. It does not clear or add records to the existing Advanced Analysis Results Dashboard.

The Advanced Analysis notebook order is now:

1. Gene Composition
2. Interface Enrichment
3. Cx/Sx Interaction
4. Spatial Graph & Neighborhood — Experimental
5. Spatial Graph Results
6. Results Dashboard

On `graph_done`, the main-thread event handler preserves `last_run`, the footer action, and Export preview; loads the run in the new panel; selects the new tab; and reports the successful/total sample count plus the selected sample's figure count.

## 4. Figure discovery

- The viewer scans only `run_dir/figures` for existing PNG files.
- Only samples with manifest status `ok` are offered in the sample selector.
- Manifest sample prefixes are compared longest-first, so names such as `PDAC_P1` and `PDAC_P1_fresh` remain distinct.
- When a successful row lacks a usable sample/source name, a known complete figure suffix can provide a safe fallback sample name.
- Known suffix matching checks the longer smoothed/unsmoothed names before general map names.
- Initial priority is H/V joint map, H_expr map, V_expr map, neighborhood heatmap, graph QC, H smoothed/unsmoothed, V smoothed/unsmoothed, then other PNGs.
- Unknown PNGs are retained under a filename-derived label.

## 5. Context QC display

For the selected sample, the viewer reads the existing `context_fields/<sample>_context_field_summary.csv`. H_expr and V_expr sections show status, matched/requested counts, coverage, method, smoothing, high-state fraction, and normalized warnings. `ok`, `skipped_qc`, one-field-only, missing-file, missing-column, empty/NaN, JSON-list, Python-list text, and semicolon warning forms are handled without warning popups.

## 6. Error handling

- Failed, cancelled, missing-status, and malformed manifest rows do not block successful samples and are listed separately.
- Missing figure directories or PNG files produce a clear empty state and disabled navigation.
- Corrupt PNGs and files deleted after loading produce an in-panel placeholder/status instead of closing the app.
- A deleted or moved run/file disables or safely rejects the related open action.
- A rerun clears only the previous Spatial Graph viewer references and selectors before loading the new run.

## 7. Compatibility

- Main Mapper, Main Mapper Map Viewer, QUBO, Import / Convert, Advanced / Experimental, and the three established Advanced Analysis modules were not changed.
- The complete `spatialtx_desktop/graph` analysis package is byte-identical to the preserved v0.4-beta baseline. H/V field calculations, graph construction, permutation/statistical code, PNG generation, and filenames are therefore unchanged.
- `spatialtx_desktop/workflow.py`, `spatialtx_desktop/advanced.py`, `spatialtx_desktop/advanced_analysis.py`, `spatialtx_desktop/advanced_ui.py`, `app_cli.py`, and importer backends are byte-identical to the baseline.
- Advanced CLI version text changed, but its options, execution behavior, calculations, output paths, and filenames did not.
- The viewer is read-only and introduces no new required output file or schema field.

## 8. Tests

Commands:

```text
python -m pytest -q -p no:cacheprovider tests/test_spatial_graph_results_ui.py
python -m pytest -q -p no:cacheprovider
```

Results:

- Focused viewer tests: 14 passed, 0 failed, 0 skipped.
- Full suite: 94 passed, 0 failed, 0 skipped.
- Python AST parse: 58 files passed.
- Desktop hidden-window smoke: exact `SpatialTX Studio Desktop v0.4.1-beta` title, required six-tab order, and `SpatialGraphResultsPanel` construction passed.
- Launcher checks: `desktop_app.py --help`, `advanced_cli.py --help`, and `app_cli.py --help` exited successfully; the first two display `v0.4.1-beta` as required.
- Remaining warnings: six existing AnnData `ImplicitModificationWarning` messages from synthetic test fixtures that transform indexes to strings. No new viewer warning was emitted.

## 9. Remaining limitations

- The viewer shows generated raster PNGs; it is not an interactive scientific plotting system.
- Very large PNGs must still be decoded in memory before the display copy is resized. Source files are never modified.
- File opening depends on the operating system and its registered default applications.
- The viewer intentionally does not recalculate, edit, or reinterpret graph results.
- H_expr and V_expr remain exploratory expression-derived context fields; V_expr is not a direct measurement of vessel density, perfusion, or functional blood supply.
