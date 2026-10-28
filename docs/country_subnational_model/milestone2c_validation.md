# Milestone 2C recovery, validation and delivery audit

## Recovery state before editing

Recovery date: 2026-10-28. Working directory matched the requested
`tb-community-risk-country-model` checkout. Branch was
`feature/country-subnational-tb-model-v1`; HEAD and upstream both pointed to
`a1690c65e579d7ee375a3852436ac56aec556705`, 0 ahead / 0 behind. No staged diff.
The recent history retained Milestones 2A/2B and the architecture milestone.
Read-only remote inspection confirmed the requested branch tip before writes
to Git refs.

Legitimate surviving changes were classified before editing:

- One tracked modification: `tests/country_tb/test_contracts.py`, authorizing
  only the isolated new developer CLI import in the isolation test.
- Eleven untracked files: six calibration/shape/observation/profile Python
  modules, one preview CLI, one focused test file and three review documents.
- One executable calibration JSON Schema was present but hidden by the
  repository-wide `*.json` ignore rule. It is explicitly included in this
  milestone commit.
- Ignored prior preview folders (synthetic calibration, AUS/VNM shapes and PRK
  unavailable result) and the starting protection receipt survived. They stay
  in place. Prior temporary directories, including an inaccessible pytest
  directory, were not removed or reused.

The surviving implementation already contained useful source preservation,
review/assay validation, smoothing, backcasts, age averaging, immutable
publication and restoration. Its core still restricted infection events and
its calibration objective was a weighted prevalence moment. Country TBI
profiles/scenario propagation, repeated-event timing, requested likelihood/WLS
methods, final review updates and commit/push were unfinished for the corrected
Milestone 2C objective. The retained tests/docs made the prior intent auditable;
completed components were extended rather than discarded.

Copies of all 13 surviving source/schema files, their SHA-256 values, staged
and unstaged diffs, Git state and protection fingerprints were preserved under
`abm/output/country_tb_m2c_recovery_20261028/` before implementation edits.
No reset, clean, checkout, restore, revert, stash, amend or deletion was used.

Process inspection found no surviving pytest or MATLAB/model run. The single
pre-existing Python process was subsequently identified, read-only, as the
VS Code Black formatter language server. Runtime Node processes belonged to
the assistant tooling; none was stopped or modified. The initial Windows
command-line query was denied, then completed through an authorized read-only
escalation.

## Completed change surface

New isolated modules in `engine/country_tb/`:
`calibration_common.py`, `pressure_shape.py`, `catalytic_core.py`,
`most_recent_event.py`, `tbi_observations.py`, `catalytic_calibration.py` and
`tbi_baseline.py`.

Additional files: `scripts/preview_country_tb_calibration.py`,
`schemas/country_tb/v1/catalytic_calibration.schema.json`,
`tests/country_tb/test_catalytic_calibration.py`, the narrow isolation-test
exception, three Milestone 2C review documents, this validation report and
`protected_fingerprints_2c.json`. The review README, status, specification
preface and decision register now describe the corrected scope. Frozen 2A
schemas, loaders and all historical 2B diagnostic code/tests remain unchanged.

[Scientific semantics](milestone2c_calibration.md) cover reinfection,
identifiability, assays, age/tail treatment and overrides.
[Results](milestone2c_examples.md) report Australia/Viet Nam diagnostics,
PRK blocking, 2020-2022 smoothing behavior, backcast/source/parameter
sensitivities and the reinfection comparison. Candidate primary surveys are
linked there without adopting their values.

## Focused validation

Final bounded regression command:

```
python -m pytest tests/country_tb -q -p no:cacheprovider --tb=short --basetemp=abm/output/country_tb_m2c_recovery_tests_20261028 --junitxml=abm/output/country_tb_m2c_recovery_20261028/focused_tests.xml
```

**263 passed in 332.87 seconds**, no skips. Two existing openpyxl
`datetime.utcnow()` deprecation warnings arose in the read-only native-source
reconstruction. The final focused 2C run separately passed **105 tests in
99.43 seconds**. The new event/likelihood subset passed 17 tests in 9.67 seconds.
After tightening the schema probability aliases, the two schema/end-to-end
checks passed again in 31.34 seconds. No unbounded application/model suite
was run.

Coverage includes zero/constant/step/rising/declining pressure; young ages and
window boundaries; recent overriding remote history; remote-only excluding
recent events; prior-remote-plus-recent subset and exact state sums; expected
event count; analytical conditional last-event CDF/inverse/moments; inverse-CDF
sampling and independent Poisson-event stochastic reconciliation; extremely
small/large hazards and rare-group log probabilities; lifetime backcasts;
supplied, aggregate and multiple-anchor known-scale recovery; zero/near-one
and impossible targets; binomial counts and no pseudo-counts; fixed assays;
sex/display/open-tail reconciliation; order/hash invariance; override
validation; malformed data and explicit missing-incidence/PRK behavior;
restoration, offline loading and code isolation.

A retained synthetic preview exposed cancellation when reconstructing recent
hazard as a subtraction of rounded total and remote exposure endpoints. The
fix retains recent endpoints directly before age averaging. Analytical and
end-to-end tests verify this correction. Interim validation failures were
resolved before the passing final regression run.

The complete developer diagnostic command generated synthetic results,
Australia/Viet Nam six-scenario profiles and the PRK blocked artifact in about
339 seconds, including frozen WPP loading/validation. All outputs stay in
`abm/output/`; old previews remain untouched. Final scenario metadata was
published as new immutable `final-*` folders, preserving intermediates.

Export verification passed for **18 root artifacts**, **221 embedded
scientific hashes**, **17 propagated profiles** and **3,400 estimated age-sex
rows**. Checks included machine schema, source/row hashes, exclusive-state
sums, complete population reconciliation and agreement between integrated
last-event density mass and stored any-TBI probability. The recovery receipt
`export_verification.json` records the folders and constant-pressure
reinfection/time-distribution comparison. Python compilation and whitespace
checks passed. Output paths and forbidden integration imports were checked.
MATLAB is installed but was not invoked: no MATLAB code or runner is changed.

## Protected-state verification

[Committed fingerprint receipt](protected_fingerprints_2c.json) compares the
pre-edit and post-test states. Every protected repository/worktree matched:

| Protected checkout | HEAD | Retained state |
| --- | --- | --- |
| Starsim prototype | `5a4a39e31d01510831068ec128c098290acf8e39` | 39 untracked report artifacts unchanged |
| No-transmission model | `f834057144b3fa17a9de6d995a7997a93ec4f08e` | Clean |
| Explicit-TBI model | `7c000dd80eef8ec9ce797579f2f53736f75224e0` | Six modified plus two untracked files unchanged, including exact diff and file hashes |
| WHO source sibling | `666088cac1e20dd1e9e52016ea857a9c48e0ba8d` | Clean |

All **128 frozen files** and every protected local/remote-tracking ref matched.
The remote SA Health release branch stayed at
`03cc16e4a52a55e10019dab16bd4f599571c4b90`; annotated tag object
`e99d36725dfc3ff40edc0ed68127e0b6cd927da5` still peels to that release commit.
No merge, tag or deployment occurred. Only the authorized feature branch and
its matching upstream may advance. The delivery response and ignored
post-push receipt record the final commit, remote tip and clean/ahead-behind
state without requiring a commit to contain its own hash.

## Remaining scientific decisions

Review empirical TBI targets, survey representativeness/year/design and
covariance; assay-specific thresholds/performance; within-year age and open-
tail distributions; historical exposure, migration/origin/residence and
infection-selective survival. Review departures from the central no-protection,
no-clearance assumptions and the provisional clock-reset rule before any
progression model. Total TBI alone cannot identify those mechanisms or event-
time detail independently. Multi-observation optimizer basin coverage and
fit-compatibility thresholds also need scientific review. A disease-forward
operator is required before any claim of disease-calibrated pressure.

Risk factors, interventions, screening/treatment, economics, DALYs/QALYs,
endogenous transmission, production UI and all engine integrations remain
deferred. The model is for planning and sequencing, not for denying care.
