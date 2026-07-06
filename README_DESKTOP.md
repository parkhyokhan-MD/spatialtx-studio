# SpatialTX Studio Desktop v0.2-beta

The v0.1 Transition Mapper workflow described below is preserved. v0.2 adds an **Advanced Analysis** top-level tab with Gene Composition, Interface Enrichment, and Cx/Sx Interaction modules. These modules consume the current sample selection, Cx/Sx gene text, and C/S/G quantiles and write only to separate timestamped `advanced_*` folders.

Windows desktop research prototype for the main `.h5ad` SpatialTX workflow.

- Creator: **Hyokhan Park, MD**
- Version: **v0.1-beta**
- Release date: **2026-07-02**
- Edition: **First public beta release for Windows**

## Start on Windows

1. Double-click `install_desktop.bat` once to install the Python dependencies.
2. Double-click `run_desktop.bat` to open the application.

If dependencies are already installed, you can start it directly:

```powershell
python desktop_app.py
```

The launcher checks common Miniconda and Anaconda locations before the system Python.

## Workflow

1. Choose a folder and scan recursively for `.h5ad` files.
2. Select one or more samples in the table.
3. Edit the C-side and S-side gene programs if needed.
4. Run scoring to create a summary CSV, per-sample metrics, selected-gene tables, and six-panel PNG maps.
5. Select exactly one sample to optimize either the C-side or S-side program.
6. Apply the optimized genes and recompute/redraw the selected samples.
7. Export the latest result as a folder or ZIP archive.

The application includes seven working tabs and is designed around a Full HD desktop:

- **Analysis** — inputs, C/S programs, thresholds, scoring execution, logs, and export
- **Map Viewer** — displays generated PNG maps inside the application with sample navigation and fit-to-window scaling
- **QUBO Optimizer** — independent C/S optimization and application, per-side or combined restoration of the original fixed gene sets, followed by explicit recompute and map redraw
- **Theory & Metrics** — the C/S/R/G model, interface rules, regimes, metric interpretation, optimizer rationale, and the assumptions and limitations of Advanced Analysis
- **Interpretation** — a sample summary table, automatically generated result explanation, gene-coverage warning, review checklist, and direct access to each PNG map
- **Advanced Tools** — opt-in raw MEX conversion, pre/post candidate workflows, and local ligand/receptor utility exports
- **About & Version** — creator, public version/date, release description, and research-use notice

### QUBO option guide

- **Genes (k)**: fixed number of genes selected for the optimized side. Smaller values produce a more compact program; larger values retain broader signal but may add redundancy. Default: 8.
- **Pool**: maximum candidate genes considered by the optimizer. A larger pool broadens the search but increases computation and potential instability. Default: 40.
- **Iterations**: simulated-annealing swap attempts. More iterations can improve the search at the cost of runtime; this does not change the requested number of selected genes. Default: 300.
- **Seed**: fixes the random search path for reproducible selection with identical input and settings. Default: 20260624.

### How QUBO optimization works

1. Build a bounded candidate gene pool.
2. Score genes for C/S alignment, directional `R`, gradient association, spatial enrichment, detection, and variance.
3. Penalize opposite-side overlap, low detection, and redundant gene pairs.
4. Formulate a binary optimization problem that selects exactly `k` genes.
5. Solve it locally with a classical simulated-annealing heuristic.
6. Apply the selected program, recompute the C/S fields, and redraw the maps.

QUBO does not simply rank genes independently. It selects a complementary combination that explains the requested spatial direction while limiting redundancy.

The scoring implementation uses per-gene z-scores, C/S program means, `R=C-S`, a six-neighbor spatial graph, local balance gradient `G`, and quantile-based localized-interface and diffuse-transition calls. The optimizer uses a side-aware, fixed-cardinality QUBO-inspired objective with a classical simulated-annealing fallback. It is not a quantum backend.

Every scored sample receives a `QC_flag` and machine-readable `QC_notes`. Checks cover C/S gene coverage, coordinate validity, unique feature names, C/S program overlap, and very small spot counts. Missing or non-finite coordinates use a fallback grid and are marked `WARN` rather than being presented as genuine tissue geometry.

## Output layout

Each run creates a timestamped folder under the chosen output root:

```text
spatialtx_run_<timestamp>/
  spatialtx_summary.csv
  run_config.json
  RUN_INFO.txt
  <sample>/
    metrics.csv
    selected_genes.csv
    <sample>_spatialtx_maps.png
```

Optimizer detail and summary CSVs are stored under `optimizer/` in the latest run folder when available.

## Opt-in Advanced tools

Advanced tools are disabled by default and require an explicit enable checkbox. They include:

- complete 10x/Visium MEX folder scanning, inspection, h5ad conversion, and validation
- pre/post h5ad pair scanning and expression candidate comparison
- receptor-like/membrane filtering and QUBO candidate-pool handoff
- sequence-annotation templates and ligand/receptor candidate skeletons
- FASTA/template export when sequence data are supplied
- read-evidence review-plan generation

### A3-A5 hypothesis-generation flow

- **A3 — Pre/Post candidate comparison:** performs exploratory condition-associated comparison using normalized mean-expression contrast and detection-fraction change. Suitable comparisons include pre/post treatment, control/treated, sample A/B, or region A/B.
- **A4 — Receptor-like/membrane filter:** applies lightweight gene-symbol heuristics to prioritize receptor-like, membrane-associated, transporter-like, and surface-like candidates for follow-up review.
- **A5 — Export candidate pool to QUBO:** preserves candidate metadata, writes a bounded QUBO input table, and loads its gene list into downstream C-side or S-side combination selection.

A3-A5 are optional advanced hypothesis-generation utilities. They do not validate drug response, receptor function, ligand-receptor binding, read-level evidence, or clinical biomarkers. A3 candidates should be described as condition-associated or exploratory candidates. A4 results should be described as receptor-like or membrane-associated candidates, not discovered or validated receptors.

The ligand/receptor and sequence utilities are local template/skeleton generators. They do not query or validate against external biological databases.

This software is for exploratory research use only and is not intended for diagnosis, treatment selection, or clinical decision-making.
