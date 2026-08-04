# GEO Flat Visium Directory Import

SpatialTX Studio v0.5-beta can convert flattened standard Visium files downloaded from GEO or similar repositories when several samples are stored in one user-selected directory. AnnData `.h5ad` remains the canonical internal analysis format.

The GSE316402 folder is only a development example. SpatialTX Studio does not depend on this dataset name or location.

## Purpose and path handling

The importer accepts any readable directory selected at runtime through the GUI or supplied through the CLI. Relative and absolute paths, spaces, non-ASCII names, mounted/network directories, external drives, and platform-native Windows, Linux, or macOS paths are handled with `pathlib.Path`. Recursive scanning is off by default.

No source path, user name, drive letter, desktop directory, series accession, or development-machine location is embedded in production code. The output directory is also chosen at runtime and is checked for writability.

## Prefix grouping

One recognized suffix is removed from the end of each filename. The remaining full string is the unique `sample_prefix`. For example:

`GSM9452684_sample_30_pre_matrix.mtx` → `GSM9452684_sample_30_pre`

Only files with exactly the same full prefix are grouped. Grouping never uses a GSM accession alone and never merges `sample_3` with `sample_30` or `pre` with `post`. Recognition is case-insensitive, while source filename case and paths are preserved. Conflicting prefixes that differ only by case and duplicate component candidates are invalid.

## Required suffixes

Each valid sample requires exactly one file in every category:

- Matrix: `_matrix.mtx` or `_matrix.mtx.gz`
- Barcodes: `_barcodes.tsv` or `_barcodes.tsv.gz`
- Features: `_features.tsv`, `_features.tsv.gz`, `_genes.tsv`, or `_genes.tsv.gz`
- Positions: `_tissue_positions.csv`, `_tissue_positions.csv.gz`, `_tissue_positions_list.csv`, or `_tissue_positions_list.csv.gz`
- Scalefactors: `_scalefactors_json.json`; `.json.gz` repository packaging is also accepted

## Optional suffixes

Images:

- `_tissue_hires_image.png`
- `_tissue_lowres_image.png`
- `_detected_tissue_image.jpg` / `.jpeg`
- `_aligned_fiducials.jpg` / `.jpeg`
- `_image.tif`, `.tiff`, `.png`, `.jpg`, or `.jpeg`

GEO repository-level `.gz` packaging of the listed images and auxiliary files is also accepted. Files are not decompressed in the source directory.

Auxiliary files:

- `_spatial_enrichment.csv`
- `_metadata.csv`
- `_annotations.csv`

Images and auxiliary tables are not required for expression conversion. Missing hires/lowres images create a warning and disable image overlay for that sample without blocking H5AD conversion.

## Validation

Before conversion, SpatialTX reports missing or duplicate required components, unreadable or malformed files, broken gzip streams, empty tables or matrices, malformed JSON, missing common optional scalefactor keys, malformed/headerless positions, nonnumeric coordinates, duplicate position barcodes, matrix/feature/barcode dimension mismatches, coordinate mismatches, unsupported or ambiguous matrix orientation, case-ambiguous prefixes, and output collisions.

The matrix orientation is recorded as `genes_x_barcodes` with an explicit transpose action or `barcodes_x_genes` with no action. Equal square feature/barcode dimensions are rejected as ambiguous rather than guessed. All filtered matrix barcodes must have coordinates; extra position rows are retained as a warning.

## Source-file preservation and canonicalization

The source directory is read-only from the importer's perspective. SpatialTX never renames, moves, deletes, rewrites, or decompresses source files and never creates source-side subfolders or caches. An output directory equal to or nested inside the selected source directory is rejected before it is created.

For each selected sample, an operating-system temporary directory receives hard links where safe or copies as a fallback. It exposes canonical `matrix.mtx`, `barcodes.tsv`, `features.tsv`, and tissue-position names to the existing MEX/MTX conversion engine. Temporary files are removed after conversion. Optional images are copied to a sample-specific output asset folder and referenced relative to the H5AD directory; image arrays are not embedded unnecessarily.

## H5AD provenance

`adata.uns["spatialtx_import"]` records source mode, selected directory name, full sample prefix, provisional accession/subject/condition, relative original filenames, SpatialTX version, timestamp, matrix shape, barcode/feature/coordinate counts, orientation/action, input hashes, warnings, and `source_files_modified=false`. Absolute local paths are kept in local run reports instead of portable AnnData provenance.

Scalefactors and memory-safe relative image references are stored under `adata.uns["spatial"]`. Coordinates are stored in `adata.obsm["spatial"]`, and position fields are retained in `adata.obs`.

## GUI workflow

1. Open **Import / Convert → GEO Flat Visium Directory**.
2. Select any readable source directory.
3. Optionally enable recursive scanning; it is off by default.
4. Scan and review every valid/warning/invalid inventory row.
5. Select one or more valid samples explicitly.
6. Choose any writable output directory.
7. Convert and review per-sample status and generated reports.
8. Optionally open outputs in Main Mapper.
9. Optionally review/edit a comparative manifest.

No conversion occurs on folder selection, and no filename-derived pairing is applied automatically.

## Comparative-manifest handoff

Filename patterns ending in values such as `sample_30_pre` or `sample_30_post` may populate provisional subject and condition fields. These values never control file grouping and are not confirmed study metadata.

An unconfirmed draft leaves `group` and `pair_id` blank and records `pairing_source=filename_inference_unconfirmed`. The GUI displays accession, parsed subject, group, pair ID, condition, batch, and notes for editing. Only an explicit confirmation dialog writes `pairing_source=user_confirmed` and offers a handoff to Comparative Analysis. Adjacent GSM accessions are never used as a pairing rule, and no analysis runs automatically.

## CLI examples

List samples without conversion:

```text
spatialtx import-geo-flat --input-dir "ANY_READABLE_DIRECTORY" --list-samples
```

Convert all valid samples:

```text
spatialtx import-geo-flat --input-dir "ANY_READABLE_DIRECTORY" --output-dir "ANY_WRITABLE_DIRECTORY"
```

Convert selected exact prefixes:

```text
spatialtx import-geo-flat --input-dir "ANY_READABLE_DIRECTORY" --samples GSM9452684_sample_30_pre GSM9452685_sample_30_post --output-dir "ANY_WRITABLE_DIRECTORY"
```

Write an unconfirmed draft manifest:

```text
spatialtx import-geo-flat --input-dir "ANY_READABLE_DIRECTORY" --output-dir "ANY_WRITABLE_DIRECTORY" --write-comparative-manifest comparative_manifest_draft.csv
```

Recursive scanning is explicit:

```text
spatialtx import-geo-flat --input-dir "ANY_READABLE_DIRECTORY" --output-dir "ANY_WRITABLE_DIRECTORY" --recursive
```

## Reports

Each run receives a unique `geo_flat_import_logs/geo_flat_import_<timestamp>/` directory under the selected output directory containing:

- `geo_flat_inventory.csv`
- `geo_flat_validation_report.csv`
- `geo_flat_import_log.json`
- `geo_flat_conversion_summary.csv`

Reports include full local paths, detected files, status, warnings/errors, orientation, dimensions, counts, image availability, duration, output H5AD, and feasible input hashes. No report is written to a fixed path.

## Limitations and troubleshooting

- This mode is for flattened standard Visium MEX/spatial files, not Visium HD bins, Xenium, CosMx, MERFISH, RDS, h5Seurat, parquet, or generic CSV expression matrices.
- It does not download GEO data or retrieve online metadata.
- It does not infer confirmed pairs or biological groups.
- A malformed or incomplete sample remains visible but cannot be converted.
- If a destination H5AD exists, conversion stops unless overwrite is explicitly approved.
- Slow or temporarily unavailable network mounts may delay or fail scanning with an explicit path error.
- Candidate discovery, ligand-receptor analysis, QUBO, AI interpretation, literature search, and multi-axis modeling are outside this importer task.
