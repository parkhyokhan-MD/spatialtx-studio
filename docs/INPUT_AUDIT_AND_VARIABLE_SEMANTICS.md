# SpatialTX v0.4 input audit and variable semantics

## Input preprocessing audit

SpatialTX audits each H5AD before graph or context-field analysis and stores the result in `adata.uns["spatialtx_input_audit"]`. The audit records a cautious preprocessing-state guess, matrix storage, platform guess, coordinate and tissue metadata, duplicate names/coordinates, library-size summaries, detected-gene summaries, zero-expression fraction, and warnings. It never normalizes or changes `AnnData.X`.

`expression_scale_guess` records raw-count, log1p-like, centered/scaled, normalized-unknown, or unknown input. `detection_metric_interpretation` states whether `expression > 0` can be interpreted as detection. Context QC prefers `adata.layers["counts"]`, `adata.layers["raw"]`, or compatible `adata.raw.X` for detection when available.

Preprocessing-state detection is heuristic. Even when two samples receive the same state guess, their score magnitudes may not be directly comparable without harmonized preprocessing.

## Native and calibrated coordinates

A **native-coordinate radius graph** uses the values present in `adata.obsm["spatial"]` when no verified physical scale is available. Its radius must not be described as micrometers.

A **calibrated physical-radius graph** requires an explicit scale source. Coordinates already expressed in micrometers use scale 1. Pixel or native coordinates may be converted using a user-supplied micrometers-per-coordinate-unit scale. SpatialTX does not infer micrometers from platform name alone.

## Spot-level variables

Sequencing-based spatial spots can contain mixtures of cells and expression programs. An `adata.obs` label is therefore treated as a spot-level annotation, not automatically as a pure cell identity.

SpatialTX distinguishes four modes:

- `categorical_state`: mutually exclusive spot-level state;
- `binary_mask`: non-exclusive spot-level membership;
- `continuous_score`: numeric score without a compositional constraint;
- `proportion_composition`: numeric spot-level composition or proportion, typically bounded by 0 and 1.

The inferred type, user-confirmed mode, missing-value handling, valid count, category counts or value range are exported for auditability.

## Permutation scope

- `whole_slide` permutes among all eligible observations in one sample.
- `within_connected_components` permutes only within each connected tissue fragment and is recommended for disconnected graphs.
- `within_user_strata` permutes only within a selected categorical `adata.obs` stratum. Missing stratum values are excluded.

Permutations never cross samples or slides. Tissue-only restriction uses `in_tissue=1` when that column is available.

Permutation P-values are exploratory and rely on exchangeability assumptions. They do not fully preserve the original spatial autocorrelation structure. BH-FDR scope is `within_sample_within_analysis_table`; it does not cover all samples, graph types, analysis families, or robustness runs.

## Three separate spatial questions

Same-spot overlap asks whether two masks are true in the same spot. Neighboring-spot association asks whether a mask in one spot is adjacent to another mask. Continuous edge association evaluates a symmetric weighted X–Y statistic across graph edges. These outputs are exported separately and must not be combined into a claim of cell-cell interaction or causality.

## H_expr and V_expr limitations

`H_expr` is a hypoxia-associated expression field, not an oxygen measurement or a true hypoxia measurement. `V_expr` is an endothelial/angiogenic expression proxy, not vessel density, vascularity, perfusion, or functional blood supply. Coverage, detection, dynamic range, dominant-gene contribution, library-size correlation, C/S/R correlation, smoothing, and high-state fraction must be reviewed before interpretation.

For centered/scaled input, `detection_fraction` is unavailable and `positive_value_fraction` is reported instead. The older absolute-z dominant contribution is retained as experimental; leave-one-gene-out correlation, mean absolute field change, variance change, and rank are exported as the primary influence audit.

## Graph robustness

The optional robustness comparison evaluates radius, lattice, and KNN graphs without replacing the primary result. Direction and significance stability, ratio/z-score variation, isolated fraction, connected components, and valid edges indicate sensitivity to graph definition. KNN remains an auxiliary robustness graph. Stable results remain exploratory and do not establish biological validation.

SpatialTX Studio is a research prototype and is not intended for diagnosis, treatment selection, or clinical decision-making.
