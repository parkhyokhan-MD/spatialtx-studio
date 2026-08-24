# SpatialTX Studio v0.65 score-domain correction report

Date: 2026-08-19  
Status: implementation and six-pair regression completed  
Scope: default-off Multi-Pair Reliability sidecars only

## Root cause and actual score path

The canonical Desktop/Comparative path is:

1. `spatialtx_desktop.workflow.score_adata` selects the present C/S gene columns from `AnnData.X` with `_gene_indices` and `_dense`.
2. `_is_count_like` selects `log1p_count_like` for the 12 validation slides.
3. `_zscore_columns(expression)` centers and scales every selected gene across spots.
4. `C_raw` and `S_raw` are means of the corresponding z-scored columns.
5. `_normalize_program_fields` applies the configured program-field mode. `raw_mean` returns the signed z-score means unchanged; it does not mean raw abundance.
6. `_smooth_program_fields` applies optional KNN/Gaussian smoothing. The six-pair run used `none`.
7. `R = C - S` produces the existing legacy balance field.
8. `spatialtx_desktop.comparative.metrics.analyze_sample` forwards those arrays to Multi-Pair and legacy result layers.

The first negative values arise at step 3, `_zscore_columns(expression)`, because subtracting each gene's across-spot mean makes below-mean expression negative by design. The old Reliability integration then passed these centered signed `scored_fields["C"]` and `scored_fields["S"]` arrays to Activity, creating the metadata/domain mismatch.

The pre-centering program abundance was previously transient in the local `expression` matrix but was not returned or cached. It is now preserved as `C_activity` and `S_activity` immediately before the z-score. For count-like inputs the Activity source is:

```text
AnnData.X
-> same selected present C/S genes
-> log1p_count_like
-> mean_across_present_program_genes
-> C_activity / S_activity
```

No centering, field scaling, smoothing, clipping, absolute value, dataset-specific shift, offset, or cohort-learned correction is applied.

Desktop Main Mapper and Desktop Comparative Analysis both route through `workflow.score_adata`; `score_h5ad` and `run_batch` delegate to it. Multi-Pair and the six-pair validation CLI route through `comparative.metrics.analyze_sample` and therefore use the same source. The older single-sample `app_cli.py` entry point still uses the separate `spatialtx_studio.runner`/`spatial_fields.build_fields` pipeline and does not expose the v0.65 Reliability Layer; it was not silently changed in this correction.

## Separated score contracts

- `balance_score_source=legacy_signed_cs`
- `balance_score_domain=signed`
- `activity_score_source=selected_gene_program_mean_after_log1p_count_like` for the validation inputs
- `activity_score_domain=nonnegative`
- `activity_source_version=v0.65-nonnegative-program-mean-v1`
- `reliability_schema_version=v0.65-reliability-v2`

Legacy `B=C_legacy-S_legacy`, `pre_B`, `post_B`, and `delta_B` are unchanged. Reliability v2 computes:

```text
activity_balance = C_activity - S_activity
A = C_activity + S_activity
D = activity_balance / (A + epsilon)
CA_strength = A - abs(activity_balance)
CA_fraction = CA_strength / (A + epsilon)
```

When `A <= epsilon`, D and CA_fraction remain `NaN`. Values and QC status are retained even when a sample/pair fails validity; QC-failed Activity results are excluded from summary conclusions.

## Validity correction

Defaults are `minimum_valid_spots=30`, `minimum_valid_fraction=0.80`, and `warning_valid_fraction=0.50`.

- `valid`: count >= 30 and fraction >= 0.80
- `warning_low_valid_fraction`: count >= 30 and 0.50 <= fraction < 0.80
- `qc_fail_insufficient_valid_fraction`: fraction < 0.50
- `qc_fail_insufficient_valid_spots`: count < 30
- `qc_fail_insufficient_valid_spots_and_fraction`: both fail

Pre/Post validity, counts, totals, fractions, and reasons are exported separately. A QC failure in either slide makes the pair a QC failure.

## Files and functions changed

- `spatialtx_desktop/workflow.py`: preserve pre-z-score Activity arrays and source provenance.
- `spatialtx_desktop/comparative/metrics.py`: cache/forward Activity arrays with a cache-contract version.
- `spatialtx_desktop/reliability/models.py`: v2 schema and fraction thresholds.
- `spatialtx_desktop/reliability/core.py`: separate legacy Balance inputs from Activity inputs.
- `spatialtx_desktop/reliability/diagnostics.py`: per-source score-domain diagnostics.
- `spatialtx_desktop/reliability/exports.py`: count+fraction validity, inference flags, diagnostics, QC figure, and compatibility/new-family FDR labels.
- `spatialtx_desktop/comparative/multi_pair.py`: explicit source wiring, diagnostic aggregation, metadata, and QC status.
- `spatialtx_desktop/multi_pair_ui.py`: source labels, valid fractions/reasons, descriptive-only warning, and Score domain tab.
- `tools/validate_v06_six_pair.py`: fixed artifacts and strict six-pair Balance regression.
- `tests/test_reliability_layer.py`: signed/nonnegative separation, identities, invalid source, fraction gates, integration, and regression coverage.

## Six-pair regression

The preserved baseline was compared for pair label/order, `pre_B`, `post_B`, `delta_B`, bootstrap CI, permutation p-value, legacy BH-FDR, and C/S gene coverage. All 72 comparisons passed with `rtol=1e-12` and `atol=1e-12`; all six golden `delta_B` values were exact within that contract.

The old `delta_B_bh_fdr` field had been corrected in a mixed B/A/D/CA family. Correcting Activity changes that family even when B is identical. To preserve the existing field, v2 keeps `delta_B_bh_fdr` as an explicitly labeled v1 compatibility statistic and adds `delta_B_reliability_v2_bh_fdr` for the source-correct v2 metric family. Signed C/S are never exported or interpreted as v2 Activity.

## Before/after source validity

| Pair | Old Pre | Old Post | Old status | New Pre | New Post | New status |
|---|---:|---:|---|---:|---:|---|
| sample_30 | 0.232% | 30.345% | QC fail | 100% (431/431) | 100% (145/145) | valid |
| sample_38 | 0.442% | 0% | QC fail | 100% (453/453) | 100% (245/245) | valid |
| sample_41 | 2.985% | 3.750% | QC fail | 100% (536/536) | 100% (1040/1040) | valid |
| sample_42 | 0.699% | 0.502% | QC fail | 100% (143/143) | 100% (996/996) | valid |
| sample_43 | 5.187% | 11.625% | valid under old count-only gate | 100% (694/694) | 100% (886/886) | valid |
| sample_44 | 0.785% | 4.015% | QC fail | 100% (1274/1274) | 100% (548/548) | valid |

All six pairs have calculable Activity and CA-strength. Direction and CA-fraction are defined only where A exceeds epsilon; the observed per-slide direction-defined fractions ranged from 48.97% to 97.29%. Median CA-strength and CA-fraction were zero in this sparse program set, which does not mean every spot lacked co-expression.

## Pair 5 and Pair 6

- Pair 5 (`sample_43`) retained legacy `delta_B=-0.22667093441959016`, descriptive permutation `p=0.000999`, and compatibility BH-FDR `0.001665`. The corrected source gives `delta_A=+0.9267605456354898`, descriptive `p=0.000999`, v2 BH-FDR `0.001665`. The Activity result changes in the opposite direction from B because the two quantities use intentionally different sources and definitions. Median co-activation did not change.
- Pair 6 (`sample_44`) now has a valid Activity comparison: `delta_A=+0.11552453009332421`, descriptive `p=0.000999`, v2 BH-FDR `0.0024975`. Its legacy Balance remains `delta_B=-0.06837231767997905`. Median co-activation did not change.

These are descriptive spot-distribution comparisons of unregistered slides, not specimen-level inference or evidence of treatment effect.

## Outputs and determinism

The final run directory contains the requested files:

- `reliability_pair_summary_fixed.csv`
- `reliability_score_domain_diagnostic.csv`
- `reliability_score_domain_diagnostic.json`
- `reliability_qc_fixed.json`
- `reliability_regression_report.json`
- `reliability_regression_diff.csv`
- updated `run_metadata.json`
- `reliability_valid_fraction_qc.png`

Two independent seed-42 reruns produced byte-identical SHA-256 values for 15 reliability CSV/JSON/PNG artifacts.

## Remaining scientific and statistical limits

- The Activity source for these count-like H5AD files is log1p count abundance and is not library-size normalized. It can remain sensitive to sequencing depth, tissue composition, and spot occupancy.
- Pre/Post slides are unregistered; spot-wise subtraction is not performed.
- Bootstrap and permutation operate on spot distributions, not biological replicates.
- There is no specimen-level treatment-effect inference, causal interpretation, response classification, or clinical claim.
- Sparse mutually exclusive program detection can yield zero median co-activation even when some individual spots are co-expressing.
- H/V remain observational single-pole contexts and are not paired-pole Reliability axes.

## Test result

The full repository suite passed: 191 tests in 121.550 seconds. Tk cleanup and AnnData implicit-index conversion warnings were emitted by existing test fixtures but did not produce test failures.
