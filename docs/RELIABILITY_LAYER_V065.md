# SpatialTX Studio v0.65 Multi-axis Reliability Layer

Status: development implementation. Research use only; not for clinical decisions.

## Compatibility contract

The Reliability Layer is an additive Multi-Pair Pre/Post sidecar with two explicit score contracts. Legacy Balance consumes the existing v0.6 signed `score_adata` C/S arrays. Activity, Direction, and Co-activation consume a separately preserved pre-z-score program abundance built from the same selected genes and missing-gene policy. Neither source is clipped, shifted, min-max scaled, or substituted for the other, and established C, S, R, masks, spatial metrics, H/V context, comparability QC, and Type A/B/C candidate labels remain unchanged.

`ReliabilityConfig.enabled` is `false` by default. When disabled, no reliability file or reliability metadata section is produced. The existing v0.6 output tables retain their schemas and values.

## Continuous quantities

For each explicit paired-pole axis:

```text
B = C_legacy_signed - S_legacy_signed
activity_balance = C_activity - S_activity
A = C_activity + S_activity
D = activity_balance / (A + epsilon)
CA_strength = A - abs(activity_balance) = 2 * min(C_activity, S_activity)
CA_fraction = CA_strength / (A + epsilon)
```

Balance is retained for finite signed legacy inputs. The Activity source is the mean of the selected program genes after the established count-like `log1p` input handling but before per-gene z-scoring, configured field normalization, or smoothing. Activity, Direction, and Co-activation are valid only where both Activity inputs are finite and non-negative. When `A <= epsilon`, Direction and `CA_fraction` are undefined (`NaN`). Negative, `NaN`, and infinite inputs are never converted to zero or a normal state.

The first negative values in the legacy path arise in `workflow.score_adata` at `_zscore_columns(expression)`: subtracting each selected gene's across-spot mean makes below-mean expression negative by design. The `normalization_mode=raw_mean` setting is applied after this z-score and therefore means an unmodified mean of signed z-scores, not raw abundance.

Validity uses both `minimum_valid_spots=30` and `minimum_valid_fraction=0.80`, with `warning_valid_fraction=0.50`. A sample below 50% is a QC failure even when it has at least 30 valid spots; 50–80% is a warning; at least 80% and at least 30 spots is valid. Activity values remain visible with their QC status, but QC-failed pairs are excluded from summary conclusions.

Direction and CA_fraction have an additional metric-level validity gate, separate from pair input validity. For each metric, `defined_n` counts finite metric values on its explicit defined mask and `defined_fraction = defined_n / valid_input_n`. PASS requires at least 30 defined spots and at least 80% support; 50–79.9% is CAUTION and below 50% or fewer than 30 spots is FAIL. Both Pre and Post must be PASS before bootstrap CI, permutation p, or BH-FDR is calculated. Descriptive medians/deltas remain available when defined, but unsupported inference remains `NaN` in files and `N/A` in the GUI. See [the metric-level QC audit](RELIABILITY_METRIC_LEVEL_QC_V065.md).

## Optional classification

Continuous mode is primary. Classified mode runs only when both thresholds are explicitly supplied and uses exactly four states:

- `low_activity`
- `c_dominant_active`
- `s_dominant_active`
- `active_coactivation_candidate`

No dataset-specific threshold is learned automatically. The word “transition” is not used for this classification.

## Gene and dependence QC

Strict cross-exclusivity canonicalizes entries by trimming whitespace, uppercasing, removing Ensembl version suffixes, and applying an explicitly supplied alias map. The same canonical gene appearing more than once within a pole, across poles, or across axes blocks the reliability run after the audit and block reason are written. Paralog and pathway relationships are not inferred as duplicates.

Gene coverage is reported per sample, axis, and pole without rescoring. H and V are included in exclusivity and coverage audit but are not reinterpreted as paired-pole reliability axes because the current H/V programs do not define opposing poles.

Dependence QC reports Pearson, Spearman, valid spot count, undefined fraction, seeded permutation p-value, and Benjamini-Hochberg FDR when at least two explicit paired-pole axes exist. Direction-based dependence rows use the same metric-level support gate; an unsupported correlation can remain descriptive while its p/FDR stays missing. High dependence is a warning only. Axes are never PCA-transformed, orthogonalized, removed, or reweighted. The pan-cancer/default axis weight remains 1.0.

## Additive outputs

- `reliability_spot_results.csv`
- `reliability_pair_summary.csv`
- `reliability_metric_qc.csv`
- `reliability_gene_coverage.csv`
- `reliability_qc.json`
- `cross_exclusivity_audit.csv`
- `axis_dependence_long.csv`
- `axis_dependence_matrix.csv`
- `axis_dependence_heatmap.png`
- `reliability_score_domain_diagnostic.csv`
- `reliability_score_domain_diagnostic.json`
- `reliability_valid_fraction_qc.png`

Spot-distribution bootstrap intervals and label-permutation p-values are descriptive for unregistered slides. Every summary and metadata payload states: “Descriptive spot-distribution comparison of unregistered slides. Not specimen-level inference and not evidence of treatment effect.” The machine-readable flags are `inference_level=spot_distribution_descriptive`, `registered_spots=false`, `biological_replicate_inference=false`, and `treatment_effect_claim_allowed=false`.

## Desktop use

Open Comparative Analysis → Multi-Pair Pre/Post and select **Enable v0.65 Reliability Layer (additive sidecars)**. Continuous mode needs only an explicit epsilon. Classified mode additionally requires Activity and Direction thresholds.

The **v0.65 Reliability** result tab explicitly labels **Legacy signed Balance** and **Nonnegative Activity/Co-activation**. It shows Pre/Post summaries, metric defined counts/fractions, QC reasons and eligibility, gene audit/coverage, axis dependence, and a Score domain diagnostic tab. Unsupported inference is shown as `p: N/A — <reason>` rather than as a numerical result. The Figure tab exposes reliability QC figures.

## Six-pair validation

The validation tool accepts `--reliability`:

```powershell
python -X utf8 -m tools.validate_v06_six_pair `
  --data-root <folder-containing-the-12-h5ad-files> `
  --output-root <validation-output-folder> `
  --run-tag v065_reliability `
  --context `
  --reliability `
  --baseline-reliability-summary <preserved-v1-reliability_pair_summary.csv>
```

Development validation used samples 30, 38, 41, 42, 43, and 44 (Pre/Post; 12 H5AD files). The metric-level gate was revalidated on 2026-08-24 against the preserved score-domain-fixed Balance summary.

- All established v0.6 Balance values and the preserved v1 Balance statistics passed 72 regression comparisons at `rtol=atol=1e-12`.
- Eighteen Reliability CSV/JSON/PNG artifacts, including the new metric-QC outputs, were byte-identical in two independent seed-42 runs. Timestamp/run-tag metadata was excluded from the byte comparison.
- The run produced 7,391 spot rows, 48 coverage rows, and 36 exclusivity audit rows with zero hard overlaps.
- The pre-z-score Activity source was finite and nonnegative for 100% of spots in all 12 slides, so all six pairs passed the Activity input gate. Direction/CA_fraction remained undefined at zero-Activity spots as required; defined fractions ranged from 48.97% to 97.29%.
- Under the new predefined metric gate, Direction pair QC was PASS 2, CAUTION 2, FAIL 2; CA_fraction was also PASS 2, CAUTION 2, FAIL 2 because the current masks are identical. Only samples 41 and 42 were inferentially eligible in both Pre and Post. Ineligible Direction/CA_fraction CI, p, and FDR values remained missing.
- Sample 43 retained its legacy Balance change (`delta_B=-0.22667093441959016`) and descriptive significance, while the corrected Activity source produced `delta_A=+0.9267605456354898` (`p=0.000999`, BH-FDR `0.001665`).
- Sample 44 now has a calculated Activity comparison (`delta_A=+0.11552453009332421`, descriptive `p=0.000999`, source-correct BH-FDR `0.0024975`). Median co-activation remained zero in all six pairs and must not be interpreted as evidence of absence at every spot.

The pre-z-score source is log1p-transformed count-like abundance, not library-size-normalized expression. Cross-slide Activity differences can therefore remain sensitive to sampling depth and composition. Registration is not performed, there are no biological replicates in this six-pair run, and all p-values remain spot-distribution descriptive only.

The complete post-fix automated suite passed with 203 tests, 1 environment-dependent skip, and 26 passing subtests. The preserved Balance comparison passed 72/72 fields at `rtol=atol=1e-12`.

This is computational and contract validation, not biological validation.

## Current scope

The v0.65 development integration is limited to Multi-Pair Pre/Post plus its validation CLI. Main Mapper, Single Pair, group Comparative Analysis, and Advanced Analysis continue to use the preserved v0.6 behavior. New biological axes, learned weights, registration, multi-slice comparison, spatial axis interaction, and AI/QUBO expansion are out of scope.
