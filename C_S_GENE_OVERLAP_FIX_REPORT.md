# SpatialTX Studio v0.4-beta C/S gene-overlap fix report

Validation date: 2026-07-15  
Status: completed in a separate working copy; no GitHub release or tag was created.

## 1. Source preservation

- User source: `SpatialTX_Studio_Desktop_v0_4_beta_PUBLIC_SOURCE.zip`
- Source SHA-256: `78717478E81ADDD2FB9F71700FB47D06E325A779BB93B5631770A3044500CD1A`
- Preserved copy: `SpatialTX_Studio_Desktop_v0_4_beta_overlap_fix_original/`
- Modified copy: `SpatialTX_Studio_Desktop_v0_4_beta_overlap_fix_work/`
- The original source, public v0.3 package, and prior v0.4-beta working copy were not modified.

## 2. Modified files

Core and selection:

- `spatialtx_desktop/gene_program_validation.py` (new)
- `spatialtx_desktop/workflow.py`
- `spatialtx_studio/gene_program.py`
- `spatialtx_studio/frame26.py`

Desktop and analysis modules:

- `spatialtx_desktop/app.py`
- `spatialtx_desktop/advanced_analysis.py`
- `spatialtx_desktop/graph/runner.py`

Tests and documentation:

- `tests/test_gene_program_overlap.py` (new)
- `tests/test_optimizer_multiseed.py`
- `README.md`
- `README_DESKTOP.md`
- `CHANGELOG.md`
- `RELEASE_NOTES_v0_4_beta.md`
- `docs/OUTPUT_SCHEMA_v0_4.md`
- `C_S_GENE_OVERLAP_FIX_REPORT.md` (new)

## 3. Overlap detection

All gene lists use one shared validator. Before comparison it:

1. trims leading and trailing whitespace;
2. converts symbols to uppercase;
3. removes empty entries;
4. removes duplicate symbols within each program while preserving first-occurrence order;
5. evaluates `C_genes ∩ S_genes` on the canonical symbols.

The result records normalized lists, overlap genes/count, duplicate removal, action, policy, mode, warnings, and validation status. Analysis entry points use `overlap_policy=error`. The audit-only `report` policy can return an invalid result but is not used to continue an analysis.

## 4. Mode policies

- **Fixed:** configured fixed lists are checked before feature availability filtering. Any overlap is a development error and stops execution.
- **Custom/Main Mapper:** the UI performs preflight validation and displays the conflicting genes. `score_adata()` repeats validation so CLI or direct module calls cannot bypass the rule.
- **Adaptive:** C-side selection is completed first; its genes are excluded from the S-side candidate pool. The final programs are validated again. Excluded genes are recorded in the selection note.
- **Advanced Analysis:** Gene Composition, Interface Enrichment, and Cx/Sx Interaction receive canonical validated lists and perform backend validation before creating a result folder.
- **Spatial Graph:** the graph batch validates C/S lists before creating a graph run; the common scoring engine validates them again.
- **QUBO:** opposite-side genes are excluded from the candidate pool before optimization, and the selected program is validated before recomputation or UI application.

S-side genes beginning with `IGH`, `IGL`, or `IGK` produce an optional non-blocking review warning. This heuristic does not remove genes or determine a biological program automatically.

## 5. QUBO hard exclusion

The intended constraint is:

```text
x_C,g + x_S,g <= 1
```

The current side-by-side optimizer implements this constraint by removing every canonical opposite-side symbol from the active candidate pool. This is equivalent for the sequential single-side workflow and avoids silently deleting a selected gene afterward.

Single- and multi-seed metadata record:

- `overlap_constraint_enabled`
- `genes_excluded_due_to_opposite_side`
- `final_overlap_genes`
- `final_overlap_count`
- final `gene_program_validation`

`final_overlap_count` is required to be zero. A non-zero result raises a selection-logic error.

## 6. Gene Composition and provenance

Gene Composition uses the validated canonical C/S lists, so the same symbol cannot be emitted under both programs in a successful analysis. Its CSV includes:

- `c_gene_count_requested`
- `s_gene_count_requested`
- `c_gene_count_used`
- `s_gene_count_used`
- `n_overlap_genes`
- `overlap_genes`
- `overlap_policy`
- `program_validation_status`

Main Mapper `run_config.json` and per-sample `parameter_log.json` record requested/used lists and validation details. CLI engine metadata, Advanced Analysis metadata, Spatial Graph run/graph metadata, and QUBO JSON summaries also retain validation provenance.

## 7. IGHA1 regression example

Input:

```text
C: GZMB, IFNG, PRF1, CD8A, CD8B, NKG7, MS4A1, IGHA1
S: IGLC3, TAGLN, SPP1, PDGFRA, ACTA2, IGHA1, FN1, COL5A1
```

Observed result:

- `overlap_genes = ["IGHA1"]`
- custom analysis is blocked before scoring or result-folder creation;
- `IGLC3` is not treated as an overlap and is not automatically removed;
- after removing `IGHA1` from one side, validation succeeds;
- `IGLC3` on the S-side produces only the optional immune-associated review warning.

## 8. Tests

Baseline before this fix: 63/63 tests passed.  
Final work copy: 80/80 tests passed.

Seventeen overlap-specific tests cover:

1. case-insensitive overlap;
2. whitespace-insensitive overlap;
3. order-preserving within-program duplicate removal;
4. the requested IGHA1/IGLC3 example;
5. validator-level fixed-program development error;
6. configured fixed-program development error;
7. Main Mapper/core blocking;
8. GUI preflight blocking and message state;
9. non-blocking S-side immunoglobulin warning;
10. CLI custom-mode blocking;
11. Advanced Analysis blocking before output creation;
12. Spatial Graph blocking before output creation;
13. Gene Composition CSV exclusivity and metadata;
14. requested/used/duplicate provenance JSON;
15. adaptive final overlap zero;
16. single- and multi-seed QUBO final overlap zero;
17. unchanged C/S/R/G, interface/diffuse masks, regime, and transition burden for equivalent non-overlapping inputs.

Existing tests were not removed and assertions were not weakened.

Windows entry-point validation also passed for `desktop_app.py --help`, `app_cli.py --help`, `advanced_cli.py --help`, and `run_desktop.bat --help`. The Tkinter/ttk GUI instantiated with the new C/S exclusivity label, reported the expected v0.4-beta title, and was explicitly destroyed without an initialization error.

## 9. Existing-result preservation

For non-overlapping programs, canonical clean input and equivalent mixed-case/whitespace/duplicate input produced numerically identical:

- C score
- S score
- R field
- gradient field
- interface-like mask
- diffuse-transition mask
- regime label
- transition burden

The fixed default lists, C/S field equations, Type A/B/C thresholds, Type B patterns, graph algorithms, H/V context definitions, and biological interpretation logic were not changed.

## 10. Remaining limitations

- The QUBO constraint is implemented by sequential opposite-side candidate exclusion, not by a joint two-side quantum/QUBO solve.
- The optional S-side immune-family warning is a simple symbol-prefix heuristic and is not a biological classifier.
- Canonical comparison assumes gene symbols can be compared by trimmed uppercase text; it does not resolve aliases, Ensembl-to-symbol mappings, or species orthologs.
- Existing result folders created by older versions are not rewritten. Re-run them with this build to obtain validation provenance.
