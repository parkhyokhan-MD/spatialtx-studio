# SpatialTX Studio Desktop v0.4.1-beta

> Source-based beta patch release. Exploratory research use only.

## Added

- In-app Spatial Graph result viewer
- Automatic navigation after Spatial Graph analysis
- Sample and figure selection
- Previous and next figure navigation
- Fit-to-window image display
- Open image and open results folder actions
- H_expr and V_expr context QC summaries
- Failed sample summary
- Direct access to generated result tables

## Changed

- Application version updated from v0.4-beta to v0.4.1-beta.
- Spatial Graph analysis completion now loads generated figures automatically and opens the dedicated **Spatial Graph Results** tab.

## Unchanged

- Main Mapper scoring
- C/S/R/G definitions
- Type A/B/C rules
- H_expr/V_expr calculations
- Graph construction
- Permutation statistics
- Output schema v0.4
- Existing filenames and timestamped output structure
- CLI options, behavior, and calculations

The new viewer reads existing Spatial Graph output files and does not add a required output artifact. In particular, the Spatial Graph module/schema version remains `0.4`, `OUTPUT_SCHEMA_v0_4.md` remains current, and the annotated H5AD naming contract remains `<sample>_spatialtx_v0_4_graph_context.h5ad`.

## Limitations

- Spatial Graph Results viewer displays existing generated PNG files and does not provide interactive scientific plotting.
- H_expr and V_expr remain exploratory expression-derived context fields.
- V_expr is not a direct measure of vessel density, perfusion, or functional blood supply.
- Graph-smoothed context fields remain exploratory sensitivity outputs.
- Very large PNG files are decoded in memory and resized for display; the source image is never modified.
- Operating-system file opening depends on the configured default application for each file type.

Outputs are exploratory and are not intended for diagnosis, treatment selection, or clinical decision-making.
