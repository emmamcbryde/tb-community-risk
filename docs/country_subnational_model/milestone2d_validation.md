# Milestone 2D validation, recovery and delivery audit

## Starting Git state and capacity-interruption recovery

Work is confined to the requested country-model checkout. Milestone 2D began
on `feature/country-subnational-tb-model-v1` at
`ef3269064bbbac32e250bea3953a6a55bce957a8`, tracking
`origin/feature/country-subnational-tb-model-v1`, clean and 0 ahead / 0 behind.
Read-only remote inspection confirmed the same starting tip. All protected
worktrees, refs and frozen files were fingerprinted before implementation.
The starting receipt is retained under
`abm/output/country_tb_m2d_audit_20261028/starting_state.json`.

The response interruption occurred after successful implementation and test
runs, before the validation report, final protection receipt, commit or push.
Recovery found the same HEAD/upstream and 0/0 state, no staged changes, three
tracked unstaged edits, ten untracked implementation/review files and the
surviving ignored JSON schema. Every surviving change belongs to Milestone
2D. Nothing was discarded, reset, cleaned, restored, checked out, reverted,
stashed, amended or deleted.

The three tracked edits were the review README, replaceable milestone status
and the narrow developer-CLI exception in the existing isolation test. The
ten untracked files were five new pure modules, the preview CLI, focused tests
and three review documents. The schema was hidden by the existing `*.json`
ignore rule and is explicitly included in the milestone commit. The new
validation report and post-test protection receipt complete that change set.
Prior outputs and all test receipts remain in place.

Surviving JUnit files confirm 337 and 74 tests with zero failures, errors or
skips. The full regression began at 2026-10-28 17:57:32 Australia/Sydney; the
latest implementation/test edit was at 17:55:22, before that run and before
the dedicated 2D run began at 17:55:27. All 38 inspected country source,
schema, test and developer-script files predate the successful regression's
start. Tracked 2C engine code is unchanged from the starting commit. The three
new review documents were written at 18:01:59, after the dedicated run but
before the full regression finished at 18:04:55. No implementation, schema or
test file changed after either successful run.

The recovery snapshot records current source SHA-256 values and receipt
checksums in `capacity_recovery_state.json`. Only documentation and audit
receipts are completed after recovery. The completed suites are accepted;
they are not repeated because a response was interrupted. No test or model
process remained. The only surviving Python process was the previously
identified VS Code formatter language server; it was not stopped or modified.

## Milestone 2C review

The retained 2C process allows repeated non-homogeneous Poisson infection
events with no acquisition protection from previous infection. Most recent
infection determines recency; recent reinfection overrides remote-only;
prior remote plus recent is an auditable subset of recent, never additional
TBI. Calibration remains one-scale within each subgroup and progression is
unconnected. No defect requiring a 2C engine change was found.

Australia and Viet Nam outputs using the invented 0.005/year scale remain
clearly diagnostic, not country TBI estimates. PRK still has no usable WHO
incidence shape and remains blocked. The 100+ population remains in
reconciliation and unestimated without direct values or reviewed ages.
The dedicated 2C review passed all 105 tests in 107.66 seconds.

## Implemented contracts and behavior

Ten revision-1.0.0 contract families are specified in the
[contract package](milestone2d_contracts.md), with the complete behavior in
[source/subgroup specification](milestone2d_sources_subgroups.md):

- baseline-TBI source and resolved baseline profile;
- population subgroup and population mixture;
- subgroup demography;
- active-TB inference targets;
- direct age/sex TBI profile;
- pure interface specification;
- synthetic scenario set;
- export receipt.

The two primary choices are **Estimate from active TB incidence** and **Use
observed or user-specified TBI prevalence**. Exactly one is required; unknown
or mixed configurations fail and there is no fallback.

Incidence inference returns
`unavailable_pending_progression_and_risk_forward_operator`, null TBI and
calibration, and the entire unresolved denominator. Recent/remote and age-
specific progression, risk-factor prevalence/effects, competing mortality,
migration/imported disease, ascertainment and applicable transmission
assumptions lack an approved forward operator. WHO estimates, notifications,
prevalent disease, screen-detected disease and imported disease remain
separate target types. Counts require denominators and observation periods.
No legacy progression model, conversion multiplier or existing runner is used.

The functional user pathway supports overall anchors, age/sex anchors with
fixed assay operators and a common within-subgroup scale, direct age/sex
profiles, and direct exclusive recent/remote/never states. Source values and
model comparisons remain preserved. Total-only direct profiles leave recency
unavailable by default. The explicitly selected hybrid policy uses the model's
conditional composition; it never overwrites directly supplied recency or
claims an observed recent fraction. Infection-time detail is not invented.

Subgroups distinguish current residence from birth, previous residence,
exposure-history geography and arrival metadata. No migration model is
created. Mixtures require reviewed disjoint membership, unique IDs,
compatible geography/year/grids and reconciled positive populations. Counts
determine weights. Incomplete groups and 100+ cells remain in denominators;
complete results/full-population prevalences are unavailable when relevant
coverage is missing. Per-component covered summaries are explicitly partial.
Direct 100+ values or reviewed terminal-age weights can resolve the tail
without inventing a terminal age.

## Synthetic scenarios: invented values only

None of the following is an Australian or other country estimate. All
exported demonstration profiles are diagnostic and not simulation-ready.

| Scenario | Result |
| --- | --- |
| 1. Overall anchor | Population 1,000; 980 covered, 20 unresolved in 100+; covered TBI matches 20%; fitted alpha 0.004642127544/year; full prevalence null |
| 2. Australian mixture | Invented Australian-born 8,000 at 4%, recent arrivals 1,000 at 30%, long-term migrants 1,000 at 18%; population-weighted total 8% across 10,000 |
| 3. Complete direct override | Supplied age-band values 0.1000, 0.2000 and 0.3000 retained unchanged in source records; total 22.16%; comparison profile preserved; all 1,000 covered including direct 100+ value |
| 4. Total-only versus explicit hybrid | Default recent/remote unavailable. Explicit hybrid gives recent 3.643806979% and remote-only 18.516193021%, retaining total 22.16%; prior-remote subset 0.949345517% stays within recent |
| 5. Active-TB mode | Required blocked status, null TBI estimate and 1,000 unresolved people |
| 6. Incomplete mixture | Australian-born group selects blocked incidence mode: 8,000 unresolved of 10,000; 2,000 covered with 24% partial TBI; complete prevalence/result null |

Scenario 2's recent prevalence is 0.882550454%, remote-only 7.117449546%, and
prior remote plus recent 0.147394847%. Recent plus remote-only equals total
TBI; the subset is not added again. Each subgroup has one scale across its
ages and sexes; no age-specific scales force a fit.

Developer command actually run:

```
python scripts/preview_country_tb_sources.py --output abm/output/country_tb_m2d_preview_20261028
```

It completed in **92.950 seconds** and published
`synthetic-sources-207fe1ab48f7` with scenario JSON, interface JSON, summary
CSV, README and checksum manifest. Outputs stay under `abm/output/` and are
not added to source control. The intended interface remains a pure
specification with incidence explicitly non-operational; no production
application imports it.

## Accepted tests and verification

The following commands completed before the response interruption:

```
python -m pytest tests/country_tb/test_catalytic_calibration.py -q -p no:cacheprovider --tb=short --junitxml=abm/output/country_tb_m2d_audit_20261028/m2c_review.xml

python -m pytest tests/country_tb/test_baseline_sources.py -q -x -p no:cacheprovider --tb=short --junitxml=abm/output/country_tb_m2d_audit_20261028/m2d_tests_final.xml

python -m pytest tests/country_tb -q -p no:cacheprovider --tb=short --basetemp=abm/output/country_tb_m2d_focused_tests_20261028 --junitxml=abm/output/country_tb_m2d_audit_20261028/country_regression.xml
```

| Accepted run | Result | Reported command runtime |
| --- | --- | ---: |
| 2C review | 105 passed | 107.66 seconds |
| Dedicated 2D | 74 passed | 119.34 seconds |
| Bounded country regression | 337 passed | 442.91 seconds |

The JUnit suite times are respectively 107.596, 119.233 and 442.801 seconds;
they exclude small runner overhead included in the command timings. The
regression reported two existing openpyxl `datetime.utcnow()` deprecation
warnings during read-only native-source reconstruction. No tests were skipped.
An earlier interim assertion-type failure is retained in `m2d_tests.xml`;
it was corrected before both accepted final runs.

Coverage includes explicit/unknown modes and no fallback; deliberately
blocked incidence; distinct disease observation types; overall, age-specific,
assay and actual-count calibration; direct source preservation; exclusive
states and prior-remote subset arithmetic; default/explicit recency policies;
subgroup/membership validation; population-weighted and incomplete mixtures;
age/sex and 100+ coverage; canonical serialization, hashes and row order;
counts/provenance; offline operation; no production exposure; and protection.

The accepted export verification checked **all ten new contract families**,
**57 unique scientific hashes**, **eight resolved profiles** and file
checksums across three root JSON artifacts. Its receipt reports
**12.496 seconds**. The recovery runs only quick compilation, schema/export
integrity and documentation-link checks plus `git diff --check`; no expensive
suite is rerun. The recovery schema/export/hash/link recheck passed in **12.900 seconds**,
including all 38 unchanged tested-file hashes, three root artifacts, 57
scientific hashes, file checksums and 52 documentation links. Compilation and
both staged/unstaged whitespace checks passed. Current source hashes are
checked against the recovery snapshot through staging/commit.

## Change inventory and isolation

The milestone contains 16 intended files:

- Five new pure modules: `active_tb_targets.py`, `baseline_tbi_source.py`,
  `composition_common.py`, `direct_tbi_profile.py`, `population_subgroups.py`.
- New developer CLI: `scripts/preview_country_tb_sources.py`.
- New focused tests: `tests/country_tb/test_baseline_sources.py`, plus the
  narrow allowed-CLI addition in `tests/country_tb/test_contracts.py`.
- New schema: `schemas/country_tb/v1/baseline_source_composition.schema.json`.
- Review README/status updates; new source/subgroup specification, contracts,
  decision register, this validation report and protection fingerprint JSON.

No production UI/navigation, APY workflow, explicit-TBI/Starsim work, MATLAB,
economics, DALYs/QALYs, interventions or frozen-release file is modified.
The fixed 2C mathematical implementation and schema are unchanged. MATLAB is
not invoked because this milestone changes no MATLAB code. No merge, tag or
deployment occurs.

## Protected state and delivery

[Fingerprint receipt](protected_fingerprints_2d.json) compares the starting
and post-test/recovery states. All four protected repositories/worktrees match:
Starsim prototype at `5a4a39e31d01510831068ec128c098290acf8e39` with 39 retained
untracked report files; no-transmission at
`f834057144b3fa17a9de6d995a7997a93ec4f08e` clean; explicit-TBI at
`7c000dd80eef8ec9ce797579f2f53736f75224e0` with all six modified and two
untracked files unchanged byte-for-byte; WHO source sibling at
`666088cac1e20dd1e9e52016ea857a9c48e0ba8d` clean.

All **128 frozen file hashes** and every protected local/remote-tracking ref
match. Read-only remote checks confirm SA Health release branch
`03cc16e4a52a55e10019dab16bd4f599571c4b90`, annotated tag object
`e99d36725dfc3ff40edc0ed68127e0b6cd927da5`, and peeled tag commit
`03cc16e4a52a55e10019dab16bd4f599571c4b90` unchanged. Only the requested
country-model branch/upstream may advance. The delivery response and ignored
post-push receipt record final commit, remote agreement and clean 0/0 state
without requiring a commit to contain its own hash.

## Remaining dependencies and next milestone

Review survey representativeness/design/covariance, assay performance,
within-age and terminal-age distributions, disjoint membership evidence,
historical exposure/migration and the inherited optimizer/compatibility
criteria. Approve a progression/risk forward operator before incidence
inference; no temporary pathway is substituted.

Recommended next work is risk-factor population allocation using user-entered
prevalence, reviewed age/sex patterns, the highest applicable reviewed RR
without additive or multiplicative combination, and explicit TBI-risk-factor
association. See [decision register](decision_register_2d.md). This work is not
started automatically; progression, migration, transmission, intervention,
economics and production integration remain deferred.
