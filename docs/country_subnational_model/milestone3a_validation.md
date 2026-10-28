# Milestone 3A validation and protection

## Starting state and retained milestones

Work began in `tb-community-risk-country-model`, branch
`feature/country-subnational-tb-model-v1`, HEAD
`5f56715b4f338338abcdbe04faabbf5a46a7471a`, upstream
`origin/feature/country-subnational-tb-model-v1`, ahead/behind 0/0, clean.
Remote branch agreement was checked before changes. The unrelated unfinished
explicit-TBI worktree was preserved; it has six modified and two untracked files.

Review confirmed repeated Poisson infection, most-recent-event recency, the
prior-remote subset without double counting, shared one-scale calibration,
explicit baseline source selection, functional user inputs, population-weighted
mixtures, visible unresolved subgroups/100+, diagnostic-only 0.005/year examples
and PRK blocking. No defect requiring a 2C/2D source change was identified.
All 38 retained implementation/schema/test/preview files remain byte-identical.
No progression model is connected; incidence inference remains deliberately
`unavailable_pending_progression_and_risk_forward_operator`.

## Bounded validation

Final source/schema/test checksums are recorded before the country suite in
`abm/output/country_tb_m3a_audit_20261028/tested_files.json`. Later checks verify
that none changed. Test XML, intermediate failures and final receipts are retained
in that ignored output directory; no interrupted evidence is overwritten or removed.

| Exact command | Result | Runtime |
| --- | --- | ---: |
| `python -m pytest tests/country_tb/test_risk_allocation.py -q -x -p no:cacheprovider --tb=short --junitxml=abm/output/country_tb_m3a_audit_20261028/risk_tests_complete.xml` | 67 passed | 23.94 s |
| `python -m pytest tests/country_tb -q -p no:cacheprovider --basetemp=abm/output/country_tb_m3a_bounded_tests_20261028 --junitxml=abm/output/country_tb_m3a_audit_20261028/country_suite.xml` | 386 passed; 18 temporary-directory setup errors; no assertion failures | 479.618 s (JUnit) |
| `python -m abm.output.country_tb_m3a_audit_20261028.retry_fixture_tests` | Only those 18 tests retried; 18 passed | 0.62 s (pytest console; 0.528 s JUnit) |
| `python -m engine.country_tb.preview_risk_allocation --output abm/output/country_tb_m3a_preview_20261028` | Six invented diagnostics exported | 4.322 s |

The full run completed all 404 cases, including all 67 new tests, but its 18
existing data-layer tests using `tmp_path` could not access the newly created
pytest directory under the Windows sandbox. The session cleanup encountered the
same `WinError 5`. Its JUnit receipt survived. The targeted retry used the needed
filesystem permission and a fresh output directory; no existing temporary files
were deleted. The exact 18 node IDs and arguments are retained in
`fixture_retry_command.json`. Matching the two XML case sets proves **404 unique
tests successfully covered**, with zero assertion failures. The full suite was
not repeated. No source/schema/test file changed between these runs. The
partition verification is recorded in `validation_summary.json`.

Focused coverage includes all eleven catalogue values and definitions, review
separation, missing versus zero, explicit allocation policies and full-precision
pattern scaling, definition mismatch/mapping, child applicability, conditional
independence and reviewed joint overlap, factor counts, max-RR ties, isolated
sensitivities/caps, feasible association endpoints, phi degeneracy, control
precedence, preservation of all factor/TBI/state marginals, structural covariance,
subgroup weighting, unresolved blocking, positive 100+ denominators, deterministic
hashes, row ordering, offline schema validation of all eleven families, invalid
structure, source/provenance preservation and protection checks. Socket creation
is prohibited in offline tests; production directories are checked for imports.

Early focused failures caught a syntax typo and two underconstrained new schema
fields; these were corrected before final validation. An isolation test initially
matched an unrelated legacy function name and was replaced with an AST import
check. These failures did not affect retained files.

The stochastic allocator is deferred, so no stochastic reconciliation is claimed.
There are no active-TB outcomes, interventions, economics or production controls.

## Diagnostic exports and checks

Final export: `abm/output/country_tb_m3a_preview_20261028/risk-diagnostics-e5697e925f30/`.
The package contains immutable canonical JSON, CSV, a diagnostic warning and a
SHA-256 manifest. The [examples](milestone3a_examples.md) report population-weighted
results with presentation rounding only. The exact command
`python -m abm.output.country_tb_m3a_audit_20261028.verify_export` passed in
9.157 s: 38 unique risk scientific hashes, all ten exported contract families,
file checksums, embedded source hashes, exact population/state marginals,
compilation of seven new Python files and 57 local documentation links. The
definition-mapping family is covered separately by the focused offline schema
test, making eleven new families in total.

## Protected-state checks

The [protection fingerprint](protected_fingerprints_3a.json) records the before/
after comparisons. Four protected sibling repositories were checked by HEAD,
branch, status, tracked/index diff hashes, every surviving changed/untracked file
and detection of new untracked files. The explicit-TBI worktree retains exactly
six modified and two untracked files. All 128 frozen artifacts remain identical.
All non-target refs remain unchanged.

Frozen SA Health release and compatibility branch commit:
`03cc16e4a52a55e10019dab16bd4f599571c4b90`.
Annotated tag object: `e99d36725dfc3ff40edc0ed68127e0b6cd927da5`;
peeled tag commit is the same frozen release commit. Remote values are checked
again with the target branch after the authorized push.

Only the new risk modules, one focused test file, one new schema, new 3A review
documents/fingerprint and the review README/status are included. Production UI,
APY, explicit-TBI, Starsim, MATLAB, economics, DALY/QALY, interventions and frozen
release files remain untouched. No merge, tag or deployment is performed.

The 16 included files are:

- `engine/country_tb/risk_common.py`
- `engine/country_tb/risk_factor_catalog.py`
- `engine/country_tb/risk_prevalence.py`
- `engine/country_tb/risk_policies.py`
- `engine/country_tb/risk_allocation.py`
- `engine/country_tb/preview_risk_allocation.py`
- `schemas/country_tb/v1/risk_allocation.schema.json`
- `tests/country_tb/test_risk_allocation.py`
- `docs/country_subnational_model/milestone3a_allocation.md`
- `docs/country_subnational_model/milestone3a_contracts.md`
- `docs/country_subnational_model/milestone3a_examples.md`
- `docs/country_subnational_model/milestone3a_validation.md`
- `docs/country_subnational_model/decision_register_3a.md`
- `docs/country_subnational_model/protected_fingerprints_3a.json`
- `docs/country_subnational_model/README.md`
- `docs/country_subnational_model/MILESTONE_STATUS.md`

`git diff --check` and `git diff --cached --check` passed. Protection was checked
again after testing, including the unchanged checksums of all 46 tested
source/schema/test files. Remote branch and frozen release/tag values matched
their starting state immediately before the milestone commit.

Commit/push receipts are retained separately under the audit output directory.
The milestone commit is identifiable by this report's Git history, avoiding a
self-referential commit hash inside the committed document.
