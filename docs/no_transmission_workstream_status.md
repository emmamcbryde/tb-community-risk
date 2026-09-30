# No-transmission workstream status

This document is replaced at the end of each milestone. It does not keep a history; `git log` does.

| Field | Value |
| --- | --- |
| Date | 2026-09-30 |
| Branch | `feature/generic-no-transmission-model` |
| Milestone | Scope documentation |
| Commit | `c485ff7` (milestone work); this status file is committed immediately after it |
| Merged, tagged or deployed | No |

## Objective

Document the scope and intended use of the screening and preventive-treatment model before
further model development. The model assumes no transmission; the separate Starsim workstream
covers settings with transmission.

## Scientific decisions made

* The model estimates **direct outcomes among the modelled population** only. It assumes that
  onward local transmission during the analysis horizon does not materially affect the comparison
  between strategies.
* It never reports secondary infections or cases prevented, transmission reduction, changes in
  force of infection or community incidence, herd effects, elimination progress or outbreak
  reduction as estimated.
* Intended use: low-incidence communities or other clearly defined populations where the
  assumption is defensible. It is **not selected on incidence alone**.
* High-incidence settings, or any setting with meaningful transmission, go to the Starsim
  workstream. Starsim is not a dependency here.
* No incidence value is a validated cut-off. WHO low-incidence terminology is kept distinct from
  the model's applicability boundary.
* A provisional three-tier applicability framework was adopted **pending scientific review**:
  * very low incidence and negligible transmission: this model may be appropriate;
  * intermediate or uncertain settings: explicit assessment, preferably compared with the
    dynamic model;
  * high incidence or sustained transmission: the dynamic model is preferred.
* Documented a second implicit assumption: infection is assigned at baseline, with no new
  infection or reinfection during follow-up. This can bias results in either direction.

## Implementation completed

* `README.md`: "Scope and intended use" section.
* `docs/no_transmission_model_scope.md`: full scope, applicability domains and framework,
  prohibited outputs, the validation needed before any gate, and open decisions.
* Cross-references in `general_app_milestone1.md`, `general_app_milestone2.md` and
  `dynamic_model_readiness_spec.md`. The dynamic-model items planned for Milestone 3 are marked as
  superseded for this branch.
* Terminology review table in `general_app_milestone2.md` (identified, not corrected).
* `tests/test_no_transmission_scope_docs.py`: documentation checks.
* No change to model calculations, epidemiological or economic logic, the SA Health release or
  frozen artifacts.

## Validation performed

* The documentation tests, frozen-release integrity tests and general-app interface tests pass
  (28 tests plus 57 subtests).
* The diff against `sa-health-apy-he-v1.0.0` shows only the allowed files (`.gitignore`,
  `README.md`, `engine/apy/calibration_policy.py`).
* The full suite and MATLAB checks were not run, because no code changed.

## Unresolved issues

* User-facing wording implies transmission effects: "not yet" phrasing, "people screened and
  treated", unqualified "cases averted" and "relative reduction" labels, and "dynamic" results
  rows (see `general_app_milestone2.md`).
* The framework, any numerical flag (for example 10 or 40 per 100,000) and the WHO and national
  citations need scientific review.
* The future of `engine/dynamic/` and the readiness specification, given the Starsim workstream.
* Whether to model ongoing exposure (background infection and reinfection) without transmission
  feedback.

## Recommended next milestone

1. Wording-only fixes for the terminology risks (no calculation changes).
2. A structured, advisory applicability assessment recorded in the population profile and shown
   in outputs.
3. A matched-scenario comparison specification with Starsim, starting with the check that Starsim
   with transmission switched off reproduces this model's results, plus an agreed definition of a
   material difference.
4. A reviewed decision on optional exogenous reinfection risk.
