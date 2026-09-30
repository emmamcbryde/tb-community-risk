# No-transmission workstream status

This document is replaced at the end of each milestone. It does not keep a history; `git log` does.

| Field | Value |
| --- | --- |
| Date | 2026-09-30 |
| Branch | `feature/generic-no-transmission-model` |
| Milestone | Terminology correction and catalytic (background-exposure) specification |
| Commit | The commit that adds this file (see `git log -1 -- docs/no_transmission_workstream_status.md`) |
| Previous milestone | Scope documentation (`c485ff7`, `7920889`) |
| Merged, tagged or deployed | No |

## Objective

1. Correct user-facing wording that implied endogenous transmission effects.
2. Specify an optional exogenous-infection-pressure (catalytic) extension, and scaffold its
   configuration schema, without changing any result.

## Scientific decisions made

* **Model boundary.** "A screening and preventive-treatment model with no endogenous transmission
  feedback." Results are direct outcomes among the modelled population, compared with no screening.
* **Future identity, if background exposure is implemented.** "... with optional exogenous
  infection pressure and no endogenous transmission feedback." The workstream keeps the name
  "no-transmission model". Outputs will say "no transmission feedback" and state the exposure
  mode; this is proposed, pending review.
* **Catalytic extension**, specified in `catalytic_infection_pressure_spec.md` and designed, not
  implemented:
  * an external annual infection hazard (infections per person-year) with three modes: `none`
    (the current model), `constant` and `time_series`;
  * no dependence on the model's own infectious people, and no intervention feedback;
  * never inferred from WHO incidence.
* **Distinct mechanisms.** Background exposure is kept separate from the entry of infected people,
  the entry of people with active TB, and endogenous transmission. Imported active TB is never
  represented as force of infection. "Importation" is not a synonym for background exposure.
* **Applicability.** The framework stays advisory. No incidence value (10, 40, 100 per 100,000 or
  any other) is an engine-selection threshold.

## Implementation completed

* `engine/model_scope.py` holds the canonical scope statements and the patterns for prohibited
  claims. The incidence notes in `trend.py` and `country.py` now reuse them, with no "not yet"
  wording.
* `app/general/terminology.py` holds one general-application label per metric, plus definitions.
  The shared, frozen `app/results_page_display.py` labels are unchanged; the Results page relabels
  its rows.
  * Direct active TB cases averted (`cumulative_cases_averted`, `nPreventedActiveTB` and
    `activeTBCasesPrevented`; one quantity, one label).
  * Relative reduction in directly modelled active TB.
  * Active TB without screening (comparator); Active TB with screening.
  * People screened, or treatment starts, per direct active TB case averted.
* The Results, Run analysis and Health economics pages use the same direct-effects caption. The
  Results page adds an "Outcome definitions" expander.
* The `LIMITATIONS.md` export states:
  * the model identity and the direct-effects boundary;
  * that no infection occurs after baseline;
  * that `dynamicComparison` fields are individual-based no-feedback results.
* `engine/profiles/background_exposure.py` is the `background_exposure_v1` schema. It performs no
  calculation and is not connected to the engine, the engine mapping or any page. Profiles without
  the block map to `none`, and the population-profile contract and its hash are unchanged.
* Docs:
  * new: `catalytic_infection_pressure_spec.md`;
  * updated: `no_transmission_model_scope.md`, the `general_app_milestone2.md` planning row and the
    status of its terminology table, and the README note on the legacy dynamic model.

## Validation performed

* New tests:
  * `test_general_scope_terminology.py` (wording, exports, label consistency);
  * `test_background_exposure_schema.py` (round-trip, hashing, legacy migration to `none`,
    validation, isolation from calculations, no Starsim dependency);
  * a catalytic-spec check in `test_no_transmission_scope_docs.py`;
  * two rendered Results and Health economics tests in `test_general_app_interface.py`.
* Broader suites: all general-app tests, frozen-release integrity, frozen reference loader, results
  presentation, legacy static interface, calibration memoisation, health-economics inputs and the
  SA Health reference package.
* **What the frozen-release check measures.**
  `test_release_files_unchanged_except_documented_exceptions` lists files that existed at tag
  `sa-health-apy-he-v1.0.0` and are now modified, deleted, renamed or type-changed
  (`git diff --diff-filter=MDRT`). It requires that list to contain only `.gitignore`, `README.md`
  and `engine/apy/calibration_policy.py`.
  * Files *added* after the release (all general-application code and documentation, including
    this milestone's new and edited files) are outside that check by design.
  * The other checks are: SHA-256 of the eight frozen reference artifacts; the tag and release
    branch pointing at `03cc16e`; and the frozen numerical headline matching the working-default
    engine configuration.
  * The earlier statement "the diff touches only three allowed files" meant *release files
    modified*, not all files changed on the branch.
* Result: 210 passed and 1 failed (376 subtests). The failure is
  `test_sa_health_reference_package.py::...test_rendered_health_economics_widgets_recalculate_without_changing_health`,
  a 90 s AppTest timeout in the legacy `pages/4_Economics.py`. None of the modules changed here are
  imported by that page, and the test fails identically on a clean checkout of the previous commit
  `7920889` (379 s). This is the ARM64-emulation timing sensitivity already documented in
  `environment_and_reproducibility.md`; it is pre-existing and unchanged.
* MATLAB and the full suite were not run: no engine, economic or MATLAB code changed.

## Unresolved issues

* Background-exposure decisions (spec section 4), especially:
  * reinfection policy and partial protection;
  * whether preventive treatment clears infection or reduces progression;
  * competing mortality, which the current simulation does not model;
  * acquisition versus progression risk factors;
  * time-series missing-year, interpolation and extrapolation rules;
  * the provisional plausibility bound (hazard of 1 per person-year).
* The applicability framework, any numerical flag, and the WHO and national citations still need
  scientific review.
* The future of `engine/dynamic/` and `dynamic_model_readiness_spec.md`, given the Starsim
  workstream.
* Reports built outside the general application (the SA Health Word report and package) keep the
  frozen wording. They should adopt the direct-effects wording only in a new, separately versioned
  release.
* Newly applied country or local incidence profiles carry the revised incidence note, so their
  profile hash differs from one built before this milestone. Results and the epidemiological
  configuration hash are unaffected.

## Recommended next milestone

**Analytical core for background exposure, isolated from the engine:**

1. Pure functions for cumulative hazard, infection probability and first-infection time
   sampling (constant, age-banded and piecewise-constant).
2. The section 10 analytical tests against the closed forms.
3. A scientific review of the reinfection, preventive-treatment-mechanism and mortality
   decisions.
4. A zero-exposure identity harness that proves mode `none` consumes no random draws. It should be
   in place before any engine wiring.
