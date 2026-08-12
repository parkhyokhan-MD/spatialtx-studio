# Comparative Spatial Transition Analysis

SpatialTX Studio Desktop v0.6-beta compares established SpatialTX sample summaries across two samples, groups, or one to six independent Pre/Post pairs. It includes an H/V computational and audit validation layer while preserving C/S + FRAME2.6 as the primary framework. It is an exploratory research workflow and is not intended for diagnosis, treatment selection, response prediction, or clinical decision-making.

## Purpose and analysis boundary

The workflow asks what changed in C/S balance, spatial gradient, localized interface-like pattern, diffuse transition pattern, transition burden, adjacency, fragmentation, and optional expression context. Each input is analyzed independently with the canonical Main Mapper engine. Comparisons are then made from sample-level metrics.

No direct spatial registration is performed. Coordinates in different H5AD files are not assumed to correspond, and there is no spot-wise subtraction. Side-by-side maps are visualization only.

Comparative Analysis does not contain candidate discovery, ligand-receptor analysis, QUBO optimization, AI/ML diagnosis, or literature search. Its H/V values reuse existing observational context definitions and are not a new classifier or composite model.

## Multi-Pair Pre/Post workflow

The existing **Single Pair / Existing** workspace remains unchanged. **Multi-Pair Pre/Post** adds one to six independent labeled pairs. Every pair uses the same C/S programs, quantile thresholds, scoring options, and graph settings. Empty rows are ignored, partially filled rows are reported, and Pre and Post must be different H5AD inputs.

For each successful pair, the runner reuses the canonical pairwise sample analyzer and exports separate Pre, Post, Delta (`Post - Pre`), safe percent change, and direction values. C, S, and R are always retained separately. The primary pair-level C/S/R state rows use the existing field medians because gene-wise standardization can make within-sample field means approximately zero by construction; the existing means remain available in separate compatibility columns. Percent change for these signed centered scores is marked `NA/not_meaningful`. Spatial outputs include localized interface-like fraction, diffuse fraction, transition burden, regime and pattern transition, adjacency, fragmentation, component ratios, and normalized component density where available.

### Three separate result layers

Multi-Pair results are presented in three layers that must not be collapsed into a single response, quality, or clinical score:

1. **Balance change** — C, S, and `R=C-S`, including exact Pre, Post, Delta, and direction values.
2. **Spatial organization change** — coordinate-dependent interface, diffuse, adjacency, fragmentation, transition-burden, and compatible topology measurements.
3. **Specimen reliability** — comparability plus technical, sampling, coverage, occupancy, and geometry QC that qualifies interpretation of Layers 1 and 2.

Specimen reliability is not a biological outcome or a calibrated probability that an interpretation is true. Type A/B/C and pattern transitions remain derived exploratory descriptors; they are not reliability, treatment-response, or efficacy classifications. The UI exposes the layers as **1 Balance change**, **2 Spatial organization**, and **3 Specimen reliability**, with a non-composite three-layer overview and in-app interpretation rules.

### Rule-based pair interpretation

The **Pair interpretation** panel keeps the three layers separate while providing concise descriptive labels. `PairInterpretationConfig` centralizes all thresholds:

| Change class input | Moderate | Large |
|---|---:|---:|
| max absolute Delta across C, S, R | 0.25 | 0.75 |
| absolute interface-fraction Delta | 0.03 | 0.10 |
| absolute diffuse-fraction Delta | 0.05 | 0.15 |
| absolute transition-burden Delta | 0.10 | 0.30 |
| absolute adjacency-fraction Delta | 0.05 | 0.15 |
| absolute fragmentation Delta | 0.10 | 0.30 |

Balance change is `Large` when the largest absolute C/S/R Delta reaches the Large threshold, `Moderate` when it reaches only the Moderate threshold, and `Minimal` otherwise. Spatial organization change is `Large` when any configured spatial metric reaches its Large threshold, `Moderate` when none is Large but at least one reaches its Moderate threshold, and `Minimal` otherwise. The exact triggering measurements are exported next to every label.

`regime_preserved=yes` means only that the broad candidate label is unchanged. `structure_preserved` combines regime preservation and spatial-change class: an unchanged regime with Minimal spatial change is `yes`, but Low comparability downgrades it to `probably`; Moderate spatial change is `probably`; a changed regime or Large spatial change is `no`. This does not imply absence of biological change.

The high-level `pair_status_flag` aligns structure, balance class, and reliability, for example `Stable structure / low-confidence pair; balance change: Moderate`. It is an exploratory reading aid, not a response, efficacy, prediction, or clinical-confidence score. A Low pair displays the explicit caution that specimen selection, sampling, or technical differences may substantially influence the observed comparison. An unchanged regime with a quantitative balance shift and Minimal spatial change displays `balance shift with preserved structure` plus a regime-preservation note.

The runner isolates failures. An unreadable, corrupted, non-spatial, gene-inadequate, or otherwise failed pair is recorded as `ERROR`; it does not terminate the remaining valid pairs.

### Parallel H/V context and QC-aware interpretation

v0.6 optionally adds H and V beside the preserved C/S and spatial results. H is the existing hypoxia-associated expression-context axis. V is the existing endothelial/angiogenic expression proxy; it is not direct perfusion, vessel density, measured vascularity, or functional blood supply. Non-centered H/V program medians are used for Pre/Post comparison because within-sample centered means are non-informative for cross-sample direction.

H/V are computed after the core result and never alter C, S, R, gradients, interface/diffuse masks, transition burden, or Type A/B/C labels. If either context program cannot be computed, its Pre/Post/Delta fields remain `NaN` and the C/S/FRAME2.6 pair continues normally. No arbitrary replacement genes are substituted.

The effective existing program is resolved before each sample audit. When the configuration value is `None`, H uses `CA9, VEGFA, SLC2A1, LDHA, ENO1, PGK1, ALDOA, BNIP3, NDRG1, ADM, EGLN3, P4HA1`; V uses `PECAM1, VWF, KDR, FLT1, ESAM, ENG, CDH5, RAMP2, EMCN, CLDN5, ANGPT2, PLVAP`. These are existing SpatialTX defaults, not new definitions.

For supported raw-count input, the non-centered spot value is the mean of `log1p` counts across matched genes. For supported nonnegative log1p input, the stored values are averaged across matched genes. Sample raw mean, median, q75, q90, transition enrichment, and coefficient of variation are retained. Low coverage, no matched genes, unsupported expression scale, unavailable graph, disabled context, and calculation failure remain `NaN` with explicit status rather than being converted to zero.

For each axis within each pair only:

1. Concatenate finite Pre and Post non-centered spot values.
2. Set the common threshold to the pooled q90.
3. Calculate each sample's high-context fraction as the fraction of spots at or above that common threshold.
4. Mark a local high-context spot when it is above the common threshold and at least one spatial-graph neighbor is also above it; export the fraction of all spots meeting both conditions.

No threshold is pooled across different pairs. Raw medians remain backward-compatible H/V aliases, but a stable zero median does not override nonzero q90, high-context fraction, or local high-context fraction.

The rule-based interpretation displays `Comparability: GOOD/CAUTION/LOW` and `Interpretation confidence: GOOD/CAUTION/LOW` next to the regime and metric changes. A Low pair remains visible but receives an explicit warning that observed differences should not be attributed directly to treatment without independent validation. This confidence label is an interpretation gate derived from specimen QC, not a probability of biological truth.

### Anatomical-site metadata

Each Multi-Pair row accepts `same_site`, `different_site`, or `unknown_site`. A `different_site` pair receives `SITE-SHIFT WARNING`: site-specific tissue composition or microenvironmental variation may contribute to observed differences. The warning does not automatically exclude the pair and does not change its comparability classification or FRAME2.6 result. Unknown metadata remains explicitly `unknown_site`.

### Specimen Comparability Gate

Comparability is evaluated before pair-level interpretation and is kept separate from the observed spatial change. The output is `Good`, `Caution`, or `Low`, with every contributing measurement and reason exported.

The compact reliability table reports the actual spot-count ratio, detected-gene ratio, observed median-count ratio, spatial-extent ratio, tissue-component count difference/ratio, and tissue-occupancy difference where available. Technical reasons are separated from sampling/geometry reasons and secondary `composition proxy` reasons. Without direct histologic annotation the software uses cautious proxy wording and does not claim a definitive tumor/stroma mismatch.

Primary checks are:

- spatial-coordinate validity;
- C-side and S-side required-gene coverage;
- spot-count fold difference;
- feature-count fold difference;
- median detected genes per spot;
- median observed matrix total per spot, with a preprocessing-compatibility warning;
- in-tissue occupancy difference when available;
- an existing low-quality-spot flag difference when available;
- coordinate bounding-box area and tissue-component fold differences when analysis succeeds.

C/S dominance fractions may be included only as explicitly labeled secondary `composition proxy` context. Because C/S may legitimately change after treatment, this proxy cannot by itself produce `Low`.

All thresholds reside in `ComparabilityConfig` in `spatialtx_desktop/comparative/multi_pair.py`. A primary `low` reason produces `Low`; two or more primary `caution` reasons also produce `Low`; one primary caution or a secondary composition-proxy warning produces `Caution`; otherwise the pair is `Good`. Missing optional QC fields are exported as `not_available` and do not fail the entire batch.

Default thresholds are:

| QC comparison | Caution | Low |
|---|---:|---:|
| spot-count fold difference | 1.50 | 2.50 |
| feature-count fold difference | 1.25 | 1.75 |
| detected-genes fold difference | 1.50 | 2.50 |
| observed library-size fold difference | 2.00 | 4.00 |
| in-tissue fraction absolute difference | 0.15 | 0.30 |
| existing low-quality fraction absolute difference | 0.15 | 0.30 |
| valid-analysis-spot fraction absolute difference | 0.10 | 0.25 |
| spatial extent-area fold difference | 2.00 | 4.00 |
| tissue-component fold difference | 2.00 | 4.00 |
| minimum C/S gene coverage | below 0.80 | below 0.50 |
| secondary C/S composition-proxy absolute difference | 0.25 | 0.45 |

The last row remains secondary even at its `low` threshold. Comparability is not a calibrated probability and does not validate a biological interpretation.

### Pair-ID validation

Before a GUI run, explicit patient/sample-like identifiers in Pre/Post filenames are compared conservatively. A mismatch produces `Possible pair-ID mismatch` and asks the user to confirm that the comparison is intended. The warning is exported but does not automatically stop the run or change the comparability class. GEO/SRA accessions alone are not treated as patient identities, so files such as different GSM accessions sharing `sample_43` can still match on the explicit sample identifier.

### Multi-Pair exports

Each run creates `comparative_multi_pair/<timestamp_or_tag>/` containing:

- `pair_results.csv`: one row per requested pair with exact Pre/Post/Delta/percent/direction values, regime/pattern transitions, comparability, status, and error;
- `balance_changes.csv`: Layer 1 C/S/R results, including compatibility mean columns, without spatial-organization or reliability fields;
- `spatial_organization_changes.csv`: Layer 2 interface, diffuse, adjacency, fragmentation, topology, and derived pattern-transition results;
- `specimen_reliability.csv`: Layer 3 pair-level Good/Caution/Low classification, reason and warning counts, unavailable-QC count, and interpretation note;
- `pair_interpretation_summary.csv`: pair-level balance/spatial classes, regime and structure preservation, reliability-qualified status, interpretive flag, safety messages, and exact threshold basis;
- `comparability_details.csv`: compact actual QC ratios/differences with technical, sampling, composition-proxy, and pair-ID reasons separated;
- `overview_interpretation.csv`: requested three-layer overview columns with numeric Deltas and interpretive flags;
- `comparability_qc.csv`: long-form technical, sampling, geometry, and secondary composition-proxy measurements with availability, severity, and reasons;
- `comparative_overview.csv`: aligned but non-composite three-layer summary plus exact deltas and directional symbols;
- `context_changes.csv`: optional H/V raw Pre, Post, Delta, and direction fields;
- `context_gene_audit.csv`: one row per sample per axis with effective genes, matched/missing/expressed genes, coverage, expression scale/source, raw method, availability status, and separately sourced single-sample centered-context q80 versus pair-pooled raw-context q90 warnings and high fractions;
- `multiaxial_pair_summary.csv`: raw C/S balance, interface, diffuse, burden, H/V, regime, comparability, site, and interpretation fields for every requested pair;
- `comparative_qc_summary.csv`: compact Good/Caution/Low result, up to three primary mismatch summaries, site warning, and interpretation confidence;
- `cohort_summary.csv`: descriptive counts only;
- `run_metadata.json`: version, immutable shared settings, thresholds, file paths, SHA-256 values where readable, pair labels, timestamp, metric registry, and safety flags;
- `figures/multi_pair_comparative_overview.png`, `figures/multiaxial_pair_overview.png`, and successful pair-level metric panels. The multiaxial figure uses independent raw-Delta panels and does not normalize or overwrite exported scientific values. Its H/V column is split into separately scaled raw-median Delta and pair-pooled high-context fraction Delta tiers; the latter uses each pair's shared raw-context q90 threshold.

Direction symbols use metric-specific tolerances; exact numeric values are never replaced. Ordinary percent change is `NA` when the Pre value is zero or near zero.

## Input requirements

Every sample must be a readable AnnData `.h5ad` file with:

- a non-empty expression matrix;
- finite spatial coordinates shaped `(n_obs, 2)` or wider in `adata.obsm["spatial"]`;
- at least one represented gene from each active C-side and S-side program;
- a unique `sample_id` in manifest workflows.

Required manifest columns are `sample_id`, `file_path`, and `group`. Optional columns are `pair_id`, `condition`, `batch`, and `notes`. Relative file paths are resolved relative to the manifest file. Paired mode requires every included `pair_id` to have exactly one reference and one target sample. Incomplete or duplicate pairs are rejected.

Invalid and failed samples remain visible in `comparative_run_manifest.csv` and `comparative_warnings.csv`. Pairwise analysis stops if either member is invalid. Batch analysis continues when at least two valid samples and both comparison groups remain.

## Comparison modes

- `pairwise`: one reference H5AD and one target H5AD.
- `paired`: matched groups such as pre/post, joined by `pair_id`.
- `unpaired`: independent reference and target groups.
- `manifest_batch`: infers paired design only when all non-empty pair IDs form complete pairs; otherwise uses an unpaired design.

Reference and target labels must be distinct. Every delta is labeled explicitly:

`delta_metric = metric_target - metric_reference`

A positive delta means Target is higher than Reference. For unpaired groups, the delta table uses target-group mean minus reference-group mean.

## Core definitions and collected metrics

For spot `x`, the unchanged Main Mapper definitions are:

- `C(x)`: C-side gene-program score.
- `S(x)`: S-side gene-program score.
- `R(x) = C(x) - S(x)`: local balance field.
- `G(x)`: local spatial gradient magnitude of `R(x)`.

Per-sample outputs include means and medians for C, S, and R; R standard deviation; gradient mean and upper quantile; localized interface-like fraction and spot count; diffuse fraction; transition burden; adjacency and fragmentation metrics; diffuse component counts and sizes; operational Type A/B/C candidate label and cautious confidence; QC, robustness, warnings, input hash, and analysis status.

### Central metric registry

Comparative display and export metadata come from one ordered registry. Each metric has an internal name, human-readable name, category, unit, scale-sensitivity flag, normalization denominator, plot group, interpretation priority, observational-only flag, and deprecation status. Plotting functions do not independently redefine metric categories. Internal metric names remain in CSV files while figures use the human-readable labels.

### Sample scale and QC context

Every independently analyzed sample now reports total spots, valid analysis spots, in-tissue spots when `obs["in_tissue"]` is available, transition-candidate spots, a tissue-area proxy, tissue graph-component count, X/Y coordinate extent, bounding-box area proxy, and mean nearest-spot spacing when feasible.

`tissue_area_proxy` is explicitly defined as `n_valid_spots`. It is a spot-count proxy, not physical area. `spatial_extent_area_proxy` is in squared coordinate units and is also not physical area. A physical-area field remains unavailable unless reliable physical calibration exists; physical and spot-count estimates are never placed under the same name.

A comparison warning is generated when the valid-spot or area-proxy ratio exceeds 1.5, the graph-component count differs by more than 2, or the spatial-extent area ratio exceeds 1.5. These checks describe sampling and tissue-scale differences; they do not correct or register sections.

### Raw and normalized topology

Existing raw component counts remain unchanged and exported. Additional metrics include:

- `diffuse_components_per_1000_valid_spots`;
- `diffuse_components_per_1000_in_tissue_spots`, when in-tissue counts exist;
- `diffuse_components_per_tissue_component`;
- `small_components_per_1000_valid_spots`;
- `transition_components_per_1000_transition_spots`;
- `interface_segments_per_1000_valid_spots`, using the existing interface-component count;
- `normalized_fragmentation_score`, defined as `1000 * (n_diffuse_components + n_interface_components) / n_valid_spots`.

Raw topology counts should not be interpreted without reviewing valid spot count, tissue extent, tissue-component count, and normalized topology metrics.

When a raw component count changes strongly but its normalized density changes only modestly, the report states that sampling scale may partly explain the raw difference. Normalization is descriptive and does not change the existing masks, component definitions, or raw values.

### Raw delta and relative change

`raw_delta` remains `Target - Reference`. Ordinary percent change is `100 * (Target - Reference) / abs(Reference)` and is omitted when the reference is zero or within the documented near-zero stability threshold. No infinite value is shown.

Symmetric percent change is:

`200 * (Target - Reference) / (abs(Target) + abs(Reference))`

It remains bounded for ordinary finite inputs and is suitable for a common relative-change display. For group analyses only, a standardized delta may be computed from the pooled sample scale. A pairwise z-score is never invented from two observations. Raw values and raw deltas remain available even when figures use standardized display values.

Operational terms such as `Type_A_candidate`, `Type_B_candidate`, `Type_C_candidate`, localized interface-like pattern, and diffuse transition pattern are exploratory. They are not validated biological subtypes.

## Statistical methods

Paired mode uses the Wilcoxon signed-rank test by default. Unpaired mode uses the Mann-Whitney U test by default. A paired t test or Welch t test is used only when explicitly selected.

Outputs report group `n`, mean, median, standard deviation, interquartile range, mean Target-minus-Reference difference, an effect size, and a seeded nonparametric bootstrap 95% confidence interval when feasible. Paired effect size is paired standardized mean change for the t-test path or paired rank-biserial direction for the Wilcoxon path. Unpaired effect size is rank-biserial for Mann-Whitney or Hedges' g for Welch's t test. Benjamini-Hochberg adjustment is applied across the comparative metric family.

Small groups (`n < 3`) are flagged. A low p-value does not establish biological, mechanistic, or clinical significance. Confidence intervals and effects are unstable at small sample size.

## Operational regime summaries

Pairwise and paired runs report reference regime, target regime, combined transition label, both cautious confidence values, and a transition confidence flag. Low confidence or failed spatial QC makes the transition uncertain. Unpaired runs report group-level candidate-regime distributions and do not infer direct transitions.

Pairwise mode uses a transition card with the two operational labels, confidence values, and `stable`, `changed`, or `uncertain` status. Paired/multi-sample mode retains the 3 by 3 matrix and shows counts plus reference-row percentages, total pairs, uncertain comparisons, and missing classifications. Unpaired group mode retains a distribution comparison. These are descriptive operational classifications, not inferred biological state transitions.

## H_expr and V_expr context

`H_expr` and `V_expr` are optional observational expression-derived context fields. Their spot-level centered or standardized maps are preserved because within-section contrast may be useful for visualization.

H/V centered sample means are expected to approach zero after within-sample centering or z-standardization and should not be used as comparative biological summaries.

The legacy `H_expr_mean` and `V_expr_mean` values remain available for backward compatibility, but they are marked `deprecated` and `non_informative_centered_mean` and are excluded from primary comparative figures. Values near `1e-17` are floating-point centering artifacts, not interpretable sample differences.

When a non-centered expression scale is available, the comparison independently summarizes raw/log1p program mean, median, 75th and 90th percentiles, variance, median absolute deviation, and coefficient of variation only when mathematically appropriate. These comparison summaries do not replace the existing spatial H/V display field.

For pairwise and group displays, `H_high_fraction` and `V_high_fraction` use the 90th percentile of pooled reference and target non-centered scores. A separate percentile is not calculated inside each sample, so the fraction is not forced to be constant. The threshold value, method (`pooled_reference_target_q90`), and score normalization method are exported.

Transition enrichment is defined robustly as:

`median(score inside the existing transition mask) - median(score outside the transition mask)`

H/V enrichment is observational only and does not influence transition detection. Variance and pooled-threshold local-hotspot fractions are also descriptive. Moran's I is not added in this step because it was not already available in the comparative module and no new spatial-statistics dependency was introduced.

H/V are computed after the core C/S result. They do not alter `R(x)`, gradients, masks, transition burden, or candidate regime labels. Missing H/V genes or unusable context graphs produce warnings and unavailable values without failing the core comparison.

## GUI workflow

1. Scan and select H5AD samples in Main Mapper.
2. Confirm C/S programs, thresholds, normalization, smoothing, and robustness settings.
3. Open **Comparative Analysis**.
4. Choose pairwise mode and two selected samples, or choose a group mode and load a manifest.
5. Set reference, target, context graph, optional H/V, statistical test, and seed.
6. Validate inputs, review all validation rows, then run.
7. Review sample metrics, deltas, operational regimes, warnings, fresh figures, and the timestamped output folder.

The Metric changes view shows reference, target, raw delta, pooled-scale standardized delta where valid, normalized counterpart where defined, symmetric percent change, scale sensitivity, observational status, and warnings. A visible banner directs the user to normalized topology metrics when sample scale differs. The default `comparative_metric_changes.png` raw-Delta overview has two non-shared x-axes: Panel A, **Primary spatial-state summary metrics**, contains localized interface-like fraction, diffuse transition fraction, transition burden score, R zero-crossing fraction, largest diffuse-component ratio, and small-component fraction; Panel B, **Topology / component complexity metrics**, contains the supported component-density measures. Every displayed bar has a numeric value label. The existing per-tissue denominator is tissue graph-component count, not physical tissue area, and is labeled accordingly.

For group analyses with a finite pooled sample scale, `comparative_metric_changes_standardized.png` provides the same two-panel organization using standardized Delta. It is an optional secondary view. Raw Delta remains the default, and pairwise standardized change is not computed from only two observations. The Figures selector also retains program, transition, graph, raw topology, normalized topology, sample scale, relative change, standardized heatmap, regime, side-by-side map, and H/V observation layers.

The standardized heatmap groups columns by the registry category and inserts category separators. Its colors are metric-wise within-run z-scores for visualization only; the transform and deterministic metric order are written to the JSON sidecar, and raw CSV values are not replaced.

The side-by-side map legend states that fill color is `R = C - S` and that the outline is the union of the existing localized interface-like and diffuse transition masks. The map still performs no registration or spot-wise subtraction.

The analysis runs in a worker thread. Cancellation is checked safely between sample computations; partial manifest, sample metrics, and log files are retained.

## CLI examples

Pairwise:

```text
spatialtx compare --sample-a sample_A.h5ad --sample-b sample_B.h5ad --outdir results_pairwise --seed 42
```

Paired manifest:

```text
spatialtx compare --manifest examples/comparative_manifest_example.csv --mode paired --reference pre --target post --outdir results_v05 --seed 42
```

Unpaired manifest without H/V context:

```text
spatialtx-compare --manifest cohort.csv --mode unpaired --reference control --target treated --outdir results_group --disable-h-expr --disable-v-expr --seed 42
```

## Output descriptions

Each run creates `comparative_analysis_<timestamp>/` and never overwrites an older run.

- `comparative_sample_metrics.csv`: canonical per-sample summaries, QC, context availability, hashes, and status.
- `comparative_delta_metrics.csv`: directed deltas and availability status.
- `comparative_group_statistics.csv`: descriptive statistics, tests, effects, confidence intervals, p-values, and FDR.
- `comparative_regime_transitions.csv`: matched operational transitions or unpaired group distributions.
- `comparative_warnings.csv`: run, validation, sample, context, small-n, and matching warnings.
- `comparative_parameters.json`: configuration, provenance, inputs, hashes, environment, warnings, failures, and explicit safety/exclusion flags.
- `comparative_run_manifest.csv`: every requested sample and its validation/analysis/cache status.
- `comparative_sample_scale.csv`: per-sample spot counts, tissue graph-component count, coordinate extents, area proxies, and spacing.
- `comparative_metric_change_table.csv`: registry metadata, reference/target values, raw delta, ordinary and symmetric percent change, normalization status, flags, and warnings.
- `comparative_normalized_metrics.csv`: per-sample normalized topology and fragmentation metrics while retaining the original raw metrics elsewhere.
- `comparative_relative_changes.csv`: relative-change view with undefined/unstable values marked or omitted from plotting.
- `comparative_scale_warnings.csv`: sampling/tissue-scale warnings and raw-versus-normalized discordance notices.
- `comparative_HV_summary.csv`: observational H/V distribution, pooled-threshold, transition-enrichment, and variability comparisons; legacy centered means are flagged.
- `comparative_summary_report.html` and `.pdf`: cautious rules-based report and figures.
- `comparative_figures/`: category-specific program, transition, graph, raw-topology, normalized-topology, sample-scale, relative-change, standardized-heatmap, regime, side-by-side R-map, summary-card, and H/V figures plus JSON sidecars. Legacy figure filenames remain as compatibility outputs.
- `comparative_logs/run.log`: progress and provenance log.
- Pairwise runs add `sample_A_summary.csv`, `sample_B_summary.csv`, and `pairwise_delta_summary.csv`.
- Group runs add `group_A_summary.csv`, `group_B_summary.csv`, `group_effect_sizes.csv`, and `group_fdr_results.csv`.

The selected output root also contains `.spatialtx_comparative_cache/`. Cache entries are keyed by input SHA-256, SpatialTX version, and all analysis settings. Delete that cache only when forced recomputation is desired; disabling cache from the CLI is available with `--no-cache`.

## Limitations

- Samples are compared after independent within-sample scoring; quantile-based masks may reflect relative distributions as well as biology.
- Unmatched tissue sections may differ in capture area, tissue completeness, orientation, fragmentation, and spot count. Normalization reduces scale domination but cannot make sections anatomically registered or causally comparable.
- No batch correction, registration, cell-type deconvolution, causal modeling, or spot matching is performed.
- The first successful reference and target sample are used for representative side-by-side maps in group reports.
- Effect sizes and bootstrap intervals can be unstable for small cohorts.
- Operational regime confidence is a cautious descriptive heuristic, not a calibrated biological probability.
- H/V fields depend on gene availability, expression scale, and graph settings and remain observational.
- Mixed-failure batches can be descriptively useful, but changing group composition may bias inference and must be reviewed.

## Troubleshooting

- **No spatial coordinates:** convert or repair the H5AD so `obsm["spatial"]` contains finite two-dimensional coordinates.
- **Insufficient C/S coverage:** choose biologically appropriate programs represented in every sample; do not silently substitute genes.
- **Incomplete pair:** ensure each `pair_id` contains exactly one reference and one target row.
- **H/V unavailable:** inspect warnings and gene coverage. Core comparison remains valid if C/S and spatial QC pass.
- **Memory preflight warning:** reduce simultaneous batch size, close other applications, or run subsets using the same fixed settings.
- **Unexpected cache reuse:** confirm the input hash and settings in `comparative_parameters.json`, or rerun with `--no-cache`.

## Reproducibility guidance

Keep input H5AD files immutable, preserve manifests, fix the seed, record the exact C/S programs and thresholds, and archive the complete timestamped run directory. Compare input SHA-256 values and SpatialTX version before interpreting reruns. Use the same reference/target direction and statistical design. Review warnings and failures before using any aggregate result.

## v0.5.5 completion status

Release validation covers one- and six-pair Multi-Pair runs, the six-pair limit and seven-pair rejection, non-contiguous populated rows, pair-isolated failure, deliberate specimen mismatch, comparability reasons, exports, metadata, direction tolerances, pairwise, paired, unpaired, manifest, H/V isolation, deterministic behavior, and unchanged existing workflows. Public-data examples remain demonstrations rather than biological validation.
