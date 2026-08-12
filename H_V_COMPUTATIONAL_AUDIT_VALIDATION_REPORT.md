# H/V computational and audit validation layer

Validation date: 2026-08-11

This report documents computational behavior and auditability only. It does not establish H/V biological validation, treatment response, perfusion, angiogenesis, or clinical utility.

Automated regression: **171 tests passed plus 26 subtests** on 2026-08-11.

## Summary flow

For each axis and sample:

1. Resolve the existing default or user-supplied gene program.
2. Audit requested, matched, missing, and expressed genes plus coverage and expression-source metadata.
3. If coverage and expression scale are supported, calculate the non-centered matched-gene program mean per spot. Raw counts are transformed with `log1p`; supported nonnegative log1p values are used as stored.
4. Retain raw mean, median, q75, q90, transition enrichment, and coefficient of variation.
5. Within each individual pair and axis, calculate `threshold = q90(concat(Pre values, Post values))`.
6. Calculate each sample's high-context fraction using that shared threshold.
7. Calculate the high-context local fraction: a spot is above the shared threshold and has at least one graph neighbor above the same threshold.

No threshold is pooled across different pairs. H/V do not modify C/S/R, transition masks, Type A/B/C, spatial metrics, comparability QC, or site warnings.

The multiaxial overview visualizes H/V in two independently scaled tiers: raw-median Delta and pair-pooled high-context fraction Delta. This presentation does not alter either metric and prevents a zero median from visually hiding focal upper-tail change.

## Six-pair core regression

The preserved v0.5.5 run with H/V disabled was compared with v0.6-dev-HV-validation using default H/V enabled. CSV string values matched exactly for:

- C, S, and R Pre/Post/Delta;
- regime Pre/Post/transition;
- interface fraction, diffuse fraction, and transition burden;
- same-side, near-zero, and opposite-side adjacency fractions;
- interface fragmentation;
- comparability classification.

Observed regimes remained:

- 5/6 `Type_A_candidate -> Type_A_candidate`;
- 1/6 `Type_A_candidate -> Type_B_candidate` (`sample_30`).

Comparability remained five `Low` and one `Caution` (`sample_41`).

## Pair 5 (`sample_43`) H/V review

Both Pre and Post had 12/12 matched H genes and 12/12 matched V genes. Coverage was 1.0 and context status was `available` for all four sample-axis rows.

| Axis | Metric | Pre | Post | Delta |
|---|---|---:|---:|---:|
| H | raw median | 0.000000 | 0.149313 | 0.149313 |
| H | sample q90 | 0.189248 | 0.383764 | 0.194516 |
| H | pair-pooled high fraction | 0.025937 | 0.163657 | 0.137720 |
| H | high-context local fraction | 0.011527 | 0.138826 | 0.127299 |
| V | raw median | 0.000000 | 0.057762 | 0.057762 |
| V | sample q90 | 0.057762 | 0.183102 | 0.125340 |
| V | pair-pooled high fraction | 0.050432 | 0.152370 | 0.101938 |
| V | high-context local fraction | 0.034582 | 0.135440 | 0.100858 |

Shared within-pair thresholds were 0.332415 for H and 0.173287 for V. The shared threshold is intentionally different from either sample's own q90.

Warning-provenance regression confirmed the reported edge case for Pair 5 Pre V: the legacy within-sample centered-context q80 high fraction was 1.0 and retained its `100.0%` warning, while the pair-pooled raw-context q90 high fraction was 0.050432 (5.04%) with an empty pair-pooled warning. These values now occupy separate columns and cannot be confused without disregarding their explicit provenance labels.

## Zero-median V acceptance example

`sample_30` had a valid V raw median of 0.0 in both Pre and Post, but focal context summaries were nonzero:

| Metric | Pre | Post |
|---|---:|---:|
| V raw median | 0.000000 | 0.000000 |
| V sample q90 | 0.057762 | 0.000000 |
| V pair-pooled high fraction | 0.222738 | 0.089655 |
| V high-context local fraction | 0.192575 | 0.075862 |

The within-pair V threshold was 0.057762. This confirms that a zero whole-tissue median no longer hides an upper-tail/local V-associated expression context signal.

## Audit and cache

- Effective default H/V gene lists are stored in `run_metadata.json` under `effective_context_programs`.
- Per-sample details are stored in `context_gene_audit.csv`.
- Legacy single-sample centered-context q80 warnings and pair-pooled raw-context q90 warnings are stored in separate provenance-labelled columns. Their corresponding high fractions are also separate; a legacy `100.0%` warning must not be read as a pair-pooled high fraction.
- Comparative cache schema is `v0.6-hv-validation-v1`; old cache keys cannot silently supply this audited output layer.
- `H_V_core_effect` remains `none; H/V do not alter C/S/R, transition masks, or Type A/B/C`.
