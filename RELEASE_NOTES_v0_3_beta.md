# SpatialTX Studio Desktop v0.3-beta

Release date: 2026-07-08

This is the public v0.3-beta source release of SpatialTX Studio Desktop.

## Import / Convert architecture

AnnData H5AD remains the canonical analysis format. The Main Mapper continues to scan and analyze `.h5ad` files only.

The dedicated **Import / Convert** tab contains two independent conversion sections:

- **Raw 10x MEX/MTX → H5AD** for `matrix.mtx`, `barcodes.tsv`, and `features.tsv`/`genes.tsv` style folders, including supported `.gz` files.
- **Raw Visium H5 + spatial → H5AD** for filtered-feature HDF5, tissue positions, scalefactors, and optional tissue images, including supported `.gz` files and GEO-style filename prefixes.

Both sections provide folder selection, sample naming, conversion, validation, output-folder access, status logging, and Main Mapper handoff.

## Raw Visium GEO-style duplicate handling

GEO-style duplicated uncompressed and `.gz` spatial files are treated as equivalent alternatives rather than conflicting matches. When both versions are present, the uncompressed file is preferred and the selection is reported as a warning in the status log. This applies to tissue positions, scalefactors JSON, high-resolution tissue PNGs, and low-resolution tissue PNGs.

## Robustness and memory-safety diagnostics

The Main Mapper now includes conservative optional diagnostics:

- smoothing modes: none, kNN mean, and Gaussian spatial smoothing
- normalization modes: raw mean, z-score, and rank/quantile normalization
- optional C/S/G threshold perturbation grid for parameter-sensitivity reporting
- sample-level `parameter_log.json` export
- AnnData memory preflight with dense float32/float64 estimates
- selected C/S gene extraction without full `AnnData.X` dense conversion

## Advanced / Experimental organization

Raw MEX conversion was removed from Advanced tools. Advanced now contains hypothesis-generation comparisons, heuristic candidate filtering, QUBO handoff, and ligand/receptor review utilities only. The underlying C/S scoring, Type A/B/C calls, QUBO implementation, and Advanced Analysis algorithms are unchanged.

## Validation and limitations

- Converted H5AD files are checked for nonzero observations and features, present and unique gene names, and valid spatial coordinates when the source format requires them.
- MEX/MTX inputs without tissue coordinates are allowed with a clear validation warning and expression-only scoring; no fallback geometry, spatial regime, interface metric, or map is generated.
- Missing, empty, malformed, or non-finite `adata.obsm["spatial"]` values produce `Spatial_QC_incomplete`, with expression-only and unavailable spatial results separated in the report.
- Seurat RDS, h5Seurat, parquet, and generic CSV import are not supported.
- This software remains an exploratory research prototype and is not intended for clinical decision-making.

Both supported converter paths have been exercised successfully. Additional representative-dataset review is recommended for future releases.
