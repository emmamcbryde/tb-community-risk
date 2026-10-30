# Milestone 3D contracts, audit and validation

## Starting state and isolation

The requested directory started **clean** at
`73d2ea32b6f9cc364206b8c8cd2cd50ace079b8b` on
`feature/country-subnational-tb-model-v1`, tracking the identically named origin
branch, ahead/behind0/0. Read-only remote inspection confirmed agreement before
edits. The starting protection fingerprint matched the committed3C protection
manifest. Existing unrelated TB/Summer processes were inspected and left alone.

No2A–3C source, test, schema or runner was changed. The numerical reference
is a new, deliberately restricted diagnostic utility. The3C public inference
status and production UI remain unchanged. New generated artifacts stay under
`abm/output/`; committed data consist only of new licensed validation snapshots.

## New versioned contracts

Machine schemas: [empirical_validation.schema.json](../../schemas/country_tb/v1/empirical_validation.schema.json).
Python validation enforces cross-field rules beyond JSON Schema. Canonical
serialization uses the retained UTF-8/sorted-key/SHA-256 primitives. Sources and
rows are sorted before sealing; missing values stay explicit. New artifacts
declare `production_available=false` and `externally_validated=false`.

| Family | Function/role |
| --- | --- |
|`country_tb.empirical_tbi_validation_inputs.v1`|Complete original records plus separately transformed positivity; counts/weighted estimates distinct |
|`country_tb.empirical_validation_split.v1`|Explicit fit/validation IDs; participant/assay/parent/nativity overlap refused; one survey cannot be declared an independent second survey |
|`country_tb.empirical_tbi_comparison.v1`|DirectionA/B, metrics/unresolved status and complete internal package; no automatic thresholds |
|`country_tb.review_parameter_set.v1`|Explicit conservative/central-candidate/higher set; no fallback or default adoption |
|`country_tb.scientific_evidence_registry.v1`|Search/selection ledger, eleven definition/effect reviews, unknown fields and evidence gaps |
|`country_tb.model_discrepancy.v1`|Thirteen explicit components; magnitude may be missing; no discrepancy fitted intoalpha |
|`country_tb.scientific_reviewer_response.v1`|Pending accept/change/reject template; no automatic approval |
|`country_tb.assay_validation_scenario.v1`|Fixed named IGRA generation/performance/provenance; pediatric transport unavailable |
|`country_tb.validation_snapshot_manifest.v1`|Immutable TBI files, byte checksums and scientific hash |
|`country_tb.mortality_validation_snapshot.v1`|Original ABS age/sex/cell values plus1-qx; no open-cell expansion or country fallback |
|`country_tb.mortality_validation_manifest.v1`|Immutable mortality payload/checksum/source identity |

The existing2C assay,3B natural-history and3C gated inversion contracts are
consumed without incompatible alteration. Field review status means internal
transcription/mathematical diagnostic scope only; it cannot approve a country
configuration or confer external sign-off.

## Exact bounded commands and results

Run from repository root. Startup-inclusive measured times are in the final
column; pytest also reports its own duration. No full unbounded country suite
was run.

| Command | Result | Runtime |
| --- | --- | ---: |
|`python -m pytest tests/country_tb/test_empirical_validation.py tests/country_tb/test_incidence_inversion.py -q --durations=5 -p no:tmpdir`|**152 passed**;52 new3D and100 retained3C|5.79s pytest;6.3422s process |
|`python -m pytest tests/country_tb/test_baseline_sources.py -q --durations=5 -p no:tmpdir -k "no_fallback or overall_known_scale or age_specific_one_scale or direct_profile_original or direct_recent_remote or subgroup_metadata or incomplete_subgroup or open_tail_direct_value or explicit_conditional_recency"`|**12 passed,62 deselected**|4.43s pytest;4.9503s process |
|`python scripts/audit_country_tb_3c_for_3d.py --output abm/output/country_tb_m3d_audit_20261030/numerical`|Synthetic recovery, multi-root refusal, gates, independent reference and refinement completed|37.4472s |
|`python scripts/validate_country_tb_empirical_3d.py --output abm/output/country_tb_m3d_validation`|96 observations,12TBI-anchored cases,9contextual forward comparisons,2blocked country inversions|49.9441s |
|`python scripts/freeze_country_tb_validation_3d.py --primary abm/output/country_tb_m3d_audit_20261030/primary`|Identical immutable re-extraction|.2303s |
|`python scripts/freeze_country_tb_mortality_3d.py --workbook abm/output/country_tb_m3d_audit_20261030/primary/ABS_lifetables_2022_2024.xlsx`|Identical immutable re-extraction,202raw/200finite rows|1.0015s |

Export verification command:

```
python scripts/check_country_tb_validation_3d.py abm/output/country_tb_m3d_validation/empirical3d-4c5d1df1e4d3/validation.json abm/output/country_tb_m3d_audit_20261030/numerical/audit3c-c6fba9e9989b/audit.json data/country_tb_validation/snapshots/tbi-validation-fe05d6d1ab8c/manifest.json data/country_tb_validation/mortality/abs-life2022-2024-6785449a5200/manifest.json
```

Passed all11new families,258scientific-hash occurrences(160unique), file
checksums and offline schema references in.6173seconds. The snapshot loader
also rederives all transformed observations from source records. Scientific
export hashes:

- Empirical package:
  `4c5d1df1e4d34d7cc0f1cc9080b5ebeb9aca9f72054541fa05c46fc7c140381c`.
- Independent audit:
  `c6fba9e9989b55449ae750c46c02b220e2000051a80990228cc667655c6510c6`.
- ABS validation snapshot:
  `6785449a520081a2e2f4c60708a9f69d03ba66a36bb286752d8525640fa1a7b3`.

Quick checks compiled/parsed11Python files and checked92local links across10
documents plus unstaged `git diff --check` in.2602seconds. The final staged
diff is checked again before commit; only authorized3D paths are included.

Tests cover source fields/licensing, bad percentages/counts/missingness, exact
counts and weighted estimates, migration preservation, immutable imports and
tampering, assay generations, leakage/overlap, empirical TBI anchor/held-out
metrics, independent synthetic incidence-to-TBI validation and incompatible
targets, mortality geography/open ages, parameter selection, factor warning,
competing-event conservation, canonical ordering/hashes, discrepancy/reviewer
contracts, no network imports/use, no production UI imports and retained bytes.
The retained3C group adds structural identifiability, gates, children/no-tail,
concentration, public-mode and refinement regressions. The selected2D group
confirms direct/anchored/user/subgroup functionality without incidence forcing.

An early52-test run passed in2.55s. A subsequent combined run had151 passes
and one fixture setup error: Windows denied pytest's private temp-directory ACL.
Using `--basetemp` also encountered that ACL issue. The new immutable-snapshot
test now creates a unique ordinary workspace directory with inherited ACL,
retains it and disables pytest's temp plugin. All152 tests then passed. No
model condition/tolerance was changed; failed temp folders/logs were preserved.
There were no warnings in the accepted focused runs.

The [machine receipt](validation_receipt_3d.json) retains exact command outputs'
checksums, runtimes, tested file hashes, export checksums and integrity checks.
Model/data/schema bytes stayed unchanged through the successful commands. The
intervening test-source edits were the documented Windows fixture location and
removal of one blank line at EOF found by the staged diff check. After that
formatting correction, `python -m pytest tests/country_tb/test_empirical_validation.py -q --durations=3 -p no:tmpdir`
again passed all52 in2.37s(2.8800s process), fingerprinting the final test file.
The unchanged100retained3C and12user-TBI results remain accepted. No expensive
export was rerun for these fixture/formatting-only changes.

## Protected-state recheck

[protected_fingerprints_3d.json](protected_fingerprints_3d.json) verifies four
protected worktrees/repositories, HEAD/branch/status, index/unstaged diff hashes,
all changed/untracked file bytes, protected refs,128frozen files and61retained
2A–3C files. The unfinished explicit-TBI worktree still has **six modified and
two untracked files**, byte-for-byte unchanged. Existing Starsim reports and
all other protected states remain intact.

Frozen SA Health branch and peeled tag remain
`03cc16e4a52a55e10019dab16bd4f599571c4b90`; annotated tag object remains
`e99d36725dfc3ff40edc0ed68127e0b6cd927da5`. Only the authorized feature branch
and its origin tracking ref may advance. No merge, tag, deployment, reset,
clean, restore, stash, amend or deletion is part of this work.

Final checks include Python compilation, local documentation links, schema
validation, `git diff --check`, explicit staged-file inspection, protected-path
isolation and post-push local/remote0/0 agreement with a clean worktree. These
checks establish software/data integrity, **not scientific parameter approval**.
