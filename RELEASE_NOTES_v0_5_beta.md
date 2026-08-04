# SpatialTX Studio Desktop v0.5-beta

Public source beta released on 2026-08-04. Exploratory research use only.

## Highlights

- Added a dedicated Comparative Analysis workspace for pairwise, paired-group, unpaired-group, and manifest-based sample-level comparisons.
- Added a path-agnostic GEO Flat Visium Directory importer with exact filename-prefix grouping, strict validation, selected conversion, read-only source handling, and optional user-reviewed manifest handoff.
- Added category-specific comparative figures, sample-scale context, normalized topology metrics, raw and symmetric relative changes, scale warnings, and descriptive regime-transition summaries.
- Added centered-H/V deprecation guidance and pooled-threshold non-centered H/V summaries when the input expression scale supports them. H/V remains observational only.
- Improved Sample A/B selection in narrow windows with filename-first labels, duplicate-name disambiguation, and separate full-path fields.
- Removed the filename-based A1 pre/post pair scanner from the Advanced / Experimental interface while retaining its backend for compatibility.
- Renamed A3 as a non-spatial Reference/Target expression candidate contrast and added direct H5AD browse controls.

## Preserved architecture

- AnnData `.h5ad` remains the canonical analysis format.
- Raw 10x/Visium inputs are converted to H5AD before Main Mapper analysis.
- Each comparative sample is scored independently by the canonical Main Mapper C/S engine.
- `C(x)`, `S(x)`, `R(x) = C(x) - S(x)`, `G(x)`, Type A/B/C candidate rules, QUBO logic, and existing Advanced algorithms are unchanged.
- Raw comparative topology counts remain available alongside normalized descriptive metrics.

## Interpretation boundaries

- Delta is always `Target - Reference`.
- Comparative maps are side-by-side displays. No direct registration or spot-wise subtraction is performed.
- Operational Type A/B/C labels and transitions are exploratory candidate summaries, not validated biological states.
- A3 is a non-spatial expression/detection contrast and must not be presented as a spatial comparative result.
- Statistical significance alone is not evidence of biological or clinical significance.
- Outputs are not intended for diagnosis, prognosis, treatment selection, or clinical decision-making.

## Validation

- Full automated suite: **149 passed** on the release preparation environment.
- Public-source packaging checks exclude H5AD files, raw matrices, local results, caches, logs, wheels, nested archives, and Git metadata.
- The three included screenshots document the Comparative Analysis workspace, side-by-side R maps, and GEO Flat Visium Directory workflow.

See `README.md`, `README_DESKTOP.md`, `VALIDATION.md`, and `docs/COMPARATIVE_ANALYSIS.md` for usage and limitations.
