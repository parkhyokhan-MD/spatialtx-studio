# SpatialTX Studio Desktop v0.5-beta

> Public source beta. Exploratory research use only.

SpatialTX Studio Desktop is an open-source research workspace for exploratory spatial transcriptomics analysis. It provides a local Python desktop application and command-line workflow for `.h5ad` inputs.

This software is a research prototype. It is not intended for diagnosis, treatment selection, or clinical decision-making. Outputs are exploratory and require independent review and validation.

## What is new in v0.5-beta

v0.5-beta adds **Comparative Spatial Transition Analysis** without changing the established single-sample Main Mapper. It compares sample-level summaries for pairwise, paired-group, unpaired-group, and CSV-manifest designs. AnnData `.h5ad` remains the canonical analysis input.

- Delta direction is always `Target - Reference`; positive values mean Target is higher.
- Default group inference uses Wilcoxon signed-rank for paired data and Mann-Whitney U for unpaired data, with effect sizes, seeded bootstrap confidence intervals, and Benjamini-Hochberg FDR.
- Operational Type A/B/C candidate changes are descriptive, confidence-flagged summaries, not validated biological state transitions.
- Comparative maps are side-by-side displays only. No direct spatial registration or spot-wise subtraction is performed.
- Optional `H_expr` and `V_expr` are observational context only and never alter `C(x)`, `S(x)`, `R(x)`, transition masks, or candidate regime labels.
- Comparative figures are separated by metric category. Existing raw topology counts are preserved alongside sample-scale context, normalized component densities, stable relative-change views, and explicit scale warnings.
- Centered H/V sample means are retained only for compatibility and excluded from primary interpretation; available non-centered H/V summaries use a pooled reference/target threshold and remain observational only.
- Timestamped CSV, JSON, HTML, PDF, logs, fresh figures, input hashes, warnings, failures, and environment metadata are exported.

This release does not add candidate discovery, ligand-receptor analysis, comparative QUBO, AI interpretation, literature search, or multi-axis modeling. See [Comparative Analysis documentation](docs/COMPARATIVE_ANALYSIS.md).

## Preserved v0.4.1 capabilities

- Local Python desktop application and CLI
- Single-sample and manifest-based batch processing
- Fixed, adaptive, and custom C/S gene-program modes
- Spatial C/S balance fields, transition summaries, QC reports, and maps
- Lightweight robustness and memory-safety diagnostics with conservative defaults
- Spot-based distance by default
- Opt-in advanced hypothesis-generation utilities
- A separate **Advanced Analysis** workspace for gene composition, interface enrichment, and local Cx/Sx spatial interaction
- Reproducible CSV tables, 300-dpi PNG figures, vector PDFs, and JSON analysis metadata
- A dedicated **Import / Convert** workspace for Raw 10x MEX/MTX, Raw Visium H5 + spatial, and prefix-grouped GEO flat Visium conversion to canonical H5AD
- An optional **Spatial Graph & Neighborhood — Experimental** workflow for sparse graph QC, context fields, and exploratory neighborhood statistics
- An in-app **Spatial Graph Results** viewer for generated graph/context figures, H_expr/V_expr context QC, failed-sample summaries, and direct result-file access
- Optional multi-seed QUBO stability analysis with selection frequency, consensus core genes, exact-k consensus, objective stability, R-field agreement, and pairwise overlap exports

The established Cx and Sx definitions, Main Mapper scoring workflow, default thresholds, Type A/B/C rules, and existing output contracts are unchanged. The FRAME2.6 CLI now delegates to that same canonical Main Mapper engine instead of maintaining an independent calculation.

## Core definitions

- `C(x)`: C-side gene-program score
- `S(x)`: S-side gene-program score
- `R(x) = C(x) - S(x)`: local C/S balance
- `G(x)`: local spatial gradient of the balance field

Interface-like and transition summaries are operational, exploratory candidates. They are not validated biological subtypes.

### C/S gene-program exclusivity

Official analyses require `C_genes ∩ S_genes = ∅`. Gene symbols are trimmed, converted to uppercase, de-duplicated within each program while preserving order, and checked again inside the canonical scoring engine. A gene present on both sides is a hard error because it can cancel in `R(x)=C(x)-S(x)` and create ambiguous composition results. Custom inputs must be corrected by the user; SpatialTX does not silently choose a side.

Fixed programs are checked as a development invariant. Adaptive and QUBO selection exclude genes already assigned to the opposite side and validate the final programs again. QUBO metadata records the exclusion constraint, excluded genes, and a final overlap count that must be zero. Successful Main Mapper, CLI, Advanced Analysis, and Spatial Graph outputs include gene-program validation provenance.

## Install and start

For the desktop GUI, install `requirements-desktop.txt`. For the legacy CLI workflow, install `requirements.txt`.

On Windows, run:

```text
install_desktop.bat
run_desktop.bat
```

Or install and launch with Python 3.11 or later:

```bash
python -m pip install -r requirements-desktop.txt
python desktop_app.py
```

See [README_DESKTOP.md](README_DESKTOP.md) for the desktop workflow and [README_local_run.md](README_local_run.md) for local CLI examples.

## Comparative Analysis quick start

Select two analyzed H5AD samples in Main Mapper, then open **Comparative Analysis** for a pairwise run. For group designs, load a manifest containing `sample_id`, `file_path`, and `group`; paired designs also require `pair_id`.

Review the sample-scale banner and normalized topology figures before interpreting raw component-count changes. Pairwise regime results use a descriptive transition card; paired cohorts retain a count-and-row-percentage transition matrix. All deltas remain `Target - Reference`, and the original raw metric values remain exported.

```bash
spatialtx compare --sample-a sample_A.h5ad --sample-b sample_B.h5ad --outdir results_pairwise --seed 42
```

```bash
spatialtx compare --manifest examples/comparative_manifest_example.csv --mode paired --reference pre --target post --outdir results_v05 --seed 42
```

Required output files include sample metrics, delta metrics, group statistics, operational regime transitions, warnings, run manifest, parameter JSON, HTML/PDF reports, figures, and logs. Statistical significance alone is not evidence of biological or clinical significance.

## v0.5-beta screenshots

The screenshots were captured during final v0.5 development testing; the published source identifies itself as v0.5-beta.

Comparative Analysis workspace with category-specific metric changes:

![SpatialTX Studio v0.5-beta Comparative Analysis metric changes](docs/screenshots/spatialtx_studio_v0_5_beta_comparative_metric_changes.png)

Side-by-side `R(x) = C(x) - S(x)` maps. These panels are descriptive displays; no registration or spot-wise subtraction is performed:

![SpatialTX Studio v0.5-beta side-by-side R maps](docs/screenshots/spatialtx_studio_v0_5_beta_comparative_side_by_side_r_maps.png)

GEO Flat Visium Directory inventory and conversion workflow:

![SpatialTX Studio v0.5-beta GEO Flat Visium Directory](docs/screenshots/spatialtx_studio_v0_5_beta_geo_flat_visium_directory.png)

## Import / Convert

The Main Mapper remains H5AD-centered. Raw data must first be converted in **Import / Convert**, then opened in the Main Mapper like any other `.h5ad` input.

- **Raw 10x MEX/MTX → H5AD**: `matrix.mtx`, `barcodes.tsv`, and `features.tsv`/`genes.tsv`, including supported `.gz` variants.
- **Raw Visium H5 + spatial → H5AD**: `filtered_feature_bc_matrix.h5`, tissue positions, scalefactors, and optional tissue images, including supported `.gz` spatial files and GEO-style filename prefixes.
- **GEO Flat Visium Directory**: scans any user-selected readable directory, groups standard MEX/spatial files by their exact full filename prefix, validates every detected sample, and converts only explicitly selected valid samples.

The flat-directory importer never groups by accession alone, does not merge similar prefixes such as `sample_3` and `sample_30`, does not auto-pair pre/post names, and writes nothing into the source directory. Filename-derived subject/condition metadata remains unconfirmed until the user reviews and confirms a draft comparative manifest. See [GEO Flat Visium Import](docs/GEO_FLAT_VISIUM_IMPORT.md).

Each section provides runtime input/output selection, conversion and validation feedback, output-folder access, and Main Mapper handoff. Seurat RDS, h5Seurat, parquet, and generic CSV import are not supported.

## Advanced Analysis quick start

In the desktop application, scan and select one or more `.h5ad` files, then open **Advanced Analysis**. The established composition, enrichment, and interaction tabs use the Cx/Sx genes and quantiles currently displayed in the main workspace.

The v0.4 line adds **Spatial Graph & Neighborhood — Experimental** as a separate opt-in tab. It builds sparse radius, Visium-lattice, or symmetric-KNN graphs; distinguishes native-coordinate from calibrated physical radius; reports input/graph/context QC; optionally calculates `H_expr` and `V_expr` context fields; and exports separate same-spot, neighboring-spot, and continuous edge association statistics. These context fields do not modify `R(x)`, Type A/B/C labels, or transition masks.

v0.4.1-beta adds an in-app Spatial Graph Results viewer. After Spatial Graph & Neighborhood analysis completes, generated Graph QC, H_expr, V_expr, H/V joint-state, and neighborhood enrichment figures can be reviewed directly inside the desktop application. The viewer is a read-only UI over the existing generated files: PNG, CSV, and JSON outputs remain in the timestamped output folder, and the output schema remains v0.4. H_expr and V_expr remain descriptive context fields and do not change Main Mapper `R(x)`, Type A/B/C calls, or transition masks.

The Advanced Analysis command-line entry point remains separate from the core CLI:

```bash
python advanced_cli.py --module composition --input sample.h5ad --output results
python advanced_cli.py --module enrichment --input sample.h5ad --output results
python advanced_cli.py --module interaction --input sample.h5ad --output results --permutations 499 --seed 20260705
python advanced_cli.py --module spatial_graph --input sample.h5ad --output results --graph-method radius --permutations 999 --seed 20260713
```

See [RELEASE_NOTES_v0_5_beta.md](RELEASE_NOTES_v0_5_beta.md) for this release and [RELEASE_NOTES_v0_4_beta.md](RELEASE_NOTES_v0_4_beta.md) for the v0.4 beta stabilization history. Earlier public source release notes remain available in the repository.

## CLI quick start

Single sample:

```bash
python app_cli.py --input sample.h5ad --output results/sample1 --gene-mode fixed
```

Batch manifest:

```bash
python app_cli.py --manifest examples/sample_manifest.csv --output results/batch --gene-mode fixed
```

The manifest must contain `sample` and `input_path` columns. Optional per-row columns include `gene_mode` and `analysis`.

Typical outputs include metrics, QC summaries, selected-gene tables, run configuration and logs, and exploratory interface/transition maps. Output folders are generated at run time and are not included in this source release.

The QUBO Optimizer supports a reproducible single-seed run and an optional multi-seed stability run. The default is 1,000 simulated-annealing iterations. Multi-seed mode repeats sequential seeds, reports `QUBO energy (lower is better)`, and separates frequency-based consensus core genes from a deterministic exact-k consensus set. Consensus reflects computational stability for the selected sample, candidate pool, objective weights, and parameters; it is not biological validation or proof of a uniquely optimal gene program.

For C/S exclusivity, the optimizer implements the constraint `x_C,g + x_S,g <= 1` by removing opposite-side genes from the candidate pool before optimization and validating the selected program afterward.

## Research-use guardrails

A3-A5 are optional hypothesis-generation utilities. A3 is a non-spatial Reference/Target expression and detection-fraction contrast; use Comparative Analysis for C/S/R spatial sample or group comparisons. The earlier filename-based A1 pre/post pair scanner is no longer displayed in Advanced / Experimental. These utilities do not discover or validate drug responses, receptor function, membrane localization, ligand-receptor binding, biomarkers, biological subtypes, or clinical effects. See [DISCLAIMER.md](DISCLAIMER.md) and [RELEASE_NOTES_v0_5_beta.md](RELEASE_NOTES_v0_5_beta.md).

## Spatial Graph and Neighborhood Analysis

See [docs/SPATIAL_GRAPH_NEIGHBORHOOD.md](docs/SPATIAL_GRAPH_NEIGHBORHOOD.md), [docs/INPUT_AUDIT_AND_VARIABLE_SEMANTICS.md](docs/INPUT_AUDIT_AND_VARIABLE_SEMANTICS.md), and [docs/OUTPUT_SCHEMA_v0_4.md](docs/OUTPUT_SCHEMA_v0_4.md).

Use cautious interpretation:

- `H_expr` is a hypoxia-associated expression field.
- `V_expr` is an endothelial/angiogenic expression proxy, not vessel density, perfusion, vascularity, or functional blood supply.
- Neighborhood enrichment P-values and FDR values describe exploratory spatial association, not causal biological interaction.
- Graph-smoothed context fields are intended for visualization and exploratory sensitivity analysis. Association statistics computed on fields smoothed over the same graph may be inflated and should not be interpreted as independent confirmatory evidence.
- Permutation P-values are exploratory and rely on exchangeability assumptions. They do not fully preserve the original spatial autocorrelation structure.
- BH-FDR is calculated within one sample and one analysis table (`within_sample_within_analysis_table`), not jointly across samples, graph types, analysis families, or robustness runs.


## Development note

AI-assisted tools were used for documentation support, code organization, and troubleshooting during development. All scientific definitions, software behavior, release decisions, and final review were performed by the author.

## References and related prior archive

Background references are listed in [REFERENCES.md](REFERENCES.md).

Earlier FRAME/ISTZ spatial transcriptomics analysis materials are archived at Zenodo: doi:10.5281/zenodo.19104105. This prior archive is provided for provenance and version lineage only. It is not included as a numbered peer-reviewed reference.

## License and citation

SpatialTX Studio Desktop is released under the Apache License 2.0. See [LICENSE](LICENSE), [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md), and [CITATION.cff](CITATION.cff).
