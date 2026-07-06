# v0.2-beta validation record

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
