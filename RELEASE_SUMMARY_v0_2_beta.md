# SpatialTX Studio Desktop v0.2-beta — Release Summary

Release date: 2026-07-05  
Developer: Hyokhan Park, MD  
Status: Research-use beta; not intended for clinical decision-making

## Release overview

SpatialTX Studio Desktop v0.2-beta extends the publicly released v0.1-beta Transition Mapper into a quantitative Cx/Sx spatial analysis platform. It deepens the existing two-axis framework without introducing additional biological axes.

The definitions of Cx and Sx are unchanged. The original Transition Mapper workflow, default gene programs, thresholds, CLI, configuration, and output contracts are preserved.

## New Advanced Analysis modules

### Gene Composition

- Calculates the relative contribution of every requested gene within the Cx and Sx programs.
- Reports transformed mean and median expression, detection fraction, contribution percentage, and missing-gene status.
- Exports CSV tables, 300-dpi PNG figures, and vector PDF figures.

### Interface Enrichment

- Compares the unchanged v0.1 interface-like regions with non-interface regions.
- Reports group means, composition percentages, fold enrichment, log2 fold enrichment, Hedges' g, two-sided Mann-Whitney p-values, and Benjamini-Hochberg FDR.
- Leaves inferential statistics unavailable when group sizes are insufficient.
- Exports CSV tables, 300-dpi PNG figures, and vector PDF figures.

### Cx/Sx Spatial Interaction

- Uses the six-nearest-neighbor spatial graph rather than relying on simple correlation.
- Quantifies neighborhood coexistence, antagonism, balance, weighted spatial overlap, and Cx/Sx-dominant edge mixing.
- Uses a reproducible seeded spatial permutation null model.
- Exports summary CSV, per-spot CSV, 300-dpi PNG, and vector PDF files.

## In-app Results Dashboard

- Adds **Run All 3 Analyses + Show Dashboard** for all selected samples.
- Displays every module/sample result in one summary table.
- Shows the complete corresponding CSV table when a summary row is selected.
- Provides direct **Open CSV** and **Open figure** actions.
- Individual module runs also update the dashboard automatically.

## Reproducibility and outputs

Each Advanced Analysis run creates a separate timestamped directory containing:

- Analysis CSV tables
- Publication-ready 300-dpi PNG figures
- Vector PDF figures
- `analysis_metadata.json` with metric definitions and parameters
- `advanced_analysis_manifest.csv` with sample-level output locations and status

A separate `advanced_cli.py` entry point supports reproducible command-line execution without changing the original `app_cli.py` workflow.

## Backward compatibility

The following supplied v0.1-beta core files were verified byte-for-byte identical in v0.2-beta:

- `app_cli.py`
- `config_default.yaml`
- `spatialtx_desktop/workflow.py`
- `spatialtx_studio/runner.py`
- `spatialtx_studio/frame26.py`
- `spatialtx_studio/gene_program.py`
- `spatialtx_studio/interface_detection.py`
- `spatialtx_studio/transition_metrics.py`
- `spatialtx_studio/transition_zone.py`

No Hx or Vx axes were introduced.

## Validation summary

- All Python source files compiled successfully.
- Six focused automated tests passed.
- All three modules completed an end-to-end synthetic 144-spot `.h5ad` analysis.
- The enrichment test produced valid interface/non-interface statistics with 16 interface and 128 non-interface observations.
- CSV, JSON, PNG, and PDF artifacts were generated successfully.
- Source spot identifiers were retained in per-spot interaction tables.
- A Tkinter UI smoke test loaded all three module manifests into the Results Dashboard and populated its detailed table.
- The packaged archive was extracted, compiled, and retested successfully.

## License and citation

SpatialTX Studio Desktop is distributed under the Apache License 2.0. Citation metadata is provided in `CITATION.cff`.
