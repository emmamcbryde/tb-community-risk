# Milestone 2A validation record

Validation date: 2026-10-27. Environment: Python 3.12.4; pytest 7.4.4; jsonschema 4.19.2. Offline loaders use only the standard library plus the existing pure WHO column contract. The maintainer ABS workbook reader uses the already available openpyxl; no frozen production dependencies were changed.

## Focused tests

Final complete run: 139 passed in 116.65 seconds.

    python -m pytest -q tests/country_tb tests/test_general_who_incidence.py tests/test_general_population_profile.py -p no:cacheprovider --tb=short

The Windows test process set tempfile.tempdir to the short repository-local abm/output/m2a_tmp before pytest.main and ran outside sandbox restrictions. This keeps temporary test data under the permitted output directory.

Coverage comprises 57 existing isolated architecture tests, 50 new Milestone 2A tests, 18 existing WHO snapshot/import tests and 14 existing population-profile tests. It includes complete duplicate decomposition, exact duplicate handling, sex/both reconciliation, vintage/status/scenario separation, conflicts/overlaps, atomic idempotence, WHO corruption/schema/year/uncertainty checks, country/territory stability, offline loading, source ranking, baseline/projection years, provenance/hashing, saved-source restoration, executable schema and protected-state checks.

An earlier run had 133 passes and one canonical rate-format failure (2E+1 rather than 20); plain decimal formatting was corrected. The next run had 138 passes and one Windows MAX_PATH failure in a deeply nested temporary publication test. The same test passed separately with the short scratch root, then all 139 passed together. No test was skipped, weakened or reported passed without execution in the final run.

Python compilation of engine/country_tb, tests/country_tb and the developer CLI passed. Local documentation links, preview schemas and profile hashes passed verification. Git whitespace checks passed before staging and are repeated on the staged change.

No multi-hour suite, MATLAB or Starsim run, projection, TBI/calibration, medical-factor generation, intervention or economics workflow was executed.

## Frozen source checksums

| Artifact | SHA-256 |
| --- | --- |
| Untouched legacy population cache | e51101fa3f5054b9aeb11999d38f6b39a3b0ab4625daa00c5c7fd52e77e535d4 |
| WHO transformed incidence | 6997e3b087d26788b21a47d09968eca87b5b03c3cc50644ae4dcdcd860acfd07 |
| Original ABS workbook | 9d436fca6a0261732f7c6d8f51c834d36707c279f3878d9a815dc4014c8d8086 |
| ABS transformed population | 07ca11bcb092f6bd62c1fc16d3229164680bc68033eac9c5c64404d89f31118a |
| Original WPP 2024 archive | 31804a296b663716236cd26c46415271cc9386fb3b2be56aff6467ad32283dc8 |
| WPP transformed population | 8028a7a3e163eaa6146554a9167a2b6e1a3a87c95d7ad8ed006f4421c36c7ddd |

The two demographic manifests include their own hashes and checksums for published-source evidence. Original sources remain preserved. The legacy raw file was not replaced by the quarantine output.

## Generated read-only previews

Successful CLI command:

    python scripts/preview_country_tb_data.py --countries AUS VNM PRK

| Country | Folder under abm/output/country_tb_m2a_preview | Scientific baseline data-profile SHA-256 |
| --- | --- | --- |
| AUS | AUS-f61fc25afb90 | 11c7e058cbe54db4af1d302c0a95044807bacdabf36a2ea5003449b7c8fc8a29 |
| VNM | VNM-cc13e19cc3c1 | f08621c9e3574eeddc01e83629578962bcd7fa4a60bec4eaf7b5d2127ba2f08a |
| PRK | PRK-de1392b6cfc8 | 78fd3c8c8ce0c98faf2842fefcd0d6e2d4c5602312d197def6653d0f0d005954 |

Each folder contains the full source profile/selection diagnostics and CSVs for incidence uncertainty, sex-specific ages, combined ages, display bands and published-both reconciliation, plus a Markdown summary. Export receipt identity includes presentation; the scientific hash excludes display-only content.

## Protected-state verification

The committed protected_fingerprints_2a.json records before/after HEADs, branches, statuses, tracked-diff hashes, untracked-content hashes and protected refs. All comparisons passed before committing. The unfinished explicit-TBI state retained its six modified tracked files and two untracked files. The Starsim worktree retained 39 untracked reports. WHO and generic worktrees remained clean.

Only the requested country-model branch may advance. Final commit/push and the post-push ref/worktree recheck are recorded in the delivery response to avoid a self-referential commit field.
