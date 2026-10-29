# Milestone 3B validation and protection

## Git and review findings

Started in `tb-community-risk-country-model`, branch
`feature/country-subnational-tb-model-v1`, clean at
`3b5005834e02a057208d883e411f763ecd832e7b`, upstream
`origin/feature/country-subnational-tb-model-v1`, ahead/behind 0/0. Remote agreement
and frozen SA Health refs were checked before edits.

Review retained repeated Poisson infection, most-recent-event recency, the
prior-remote subset without double counting, one-scale latent-TBI calibration,
explicit source selection, functional user inputs, count-weighted mixtures,
baseline TBI/risk marginals, highest-RR multimorbidity, child applicability and
visible 100+ denominators. No defect requiring a retained 2C–3A file change was
identified. All 46 retained source/schema/test/preview files remain byte-identical.

The adult source audit did not verify the proposed two-year knot and did not
establish an unqualified portable cause-specific curve. Pediatric post-baseline
curves remain external-required. These are explicit scientific limits, not filled
with compatibility-engine values. See the [evidence register](milestone3b_evidence.md).

## Bounded tests

New suite:

```
python -m pytest tests/country_tb/test_progression_forward.py -q -x -p no:cacheprovider --tb=short --junitxml=abm/output/country_tb_m3b_audit_20261029/progression_complete.xml
```

The initial 63-test run passed in 37.58 seconds. Final review then added two
effect-followup applicability tests and guards; the final combined run below
validates all 65 new tests.

Combined focused regression command:

```
python -m pytest tests/country_tb/test_progression_forward.py tests/country_tb/test_risk_allocation.py tests/country_tb/test_catalytic_calibration.py::test_reinfection_states_and_event_count tests/country_tb/test_catalytic_calibration.py::test_last_event_analytic_cdf_inverse_mean_and_remote_condition tests/country_tb/test_catalytic_calibration.py::test_population_sexes_do_not_infer_different_acquisition_and_tail_is_honest tests/country_tb/test_catalytic_calibration.py::test_real_frozen_who_vnm_and_prk_missingness tests/country_tb/test_baseline_sources.py::test_incidence_mode_never_calls_catalytic_or_progression_functions tests/country_tb/test_baseline_sources.py::test_direct_recent_remote_retained_and_subset_not_added tests/country_tb/test_baseline_sources.py::test_incomplete_subgroup_blocks_full_mixture_but_keeps_denominator tests/country_tb/test_contracts.py::test_no_protected_app_runner_or_engine_imports_country_prototype -q -p no:cacheprovider --junitxml=abm/output/country_tb_m3b_audit_20261029/combined_final.xml
```

**143 passed in 56.57 seconds**: all 65 new tests, all 67 retained 3A tests, and
11 retained reinfection/source/age-tail/isolation checks. No multi-hour suite or
full historical suite was needed. No test/model assertion remains failed.
The earlier 141-test focused run passed in 71.09 seconds; it was repeated only
because the finite-age-scope guard was subsequently added.

Coverage includes risk-to-hazard conversion and all knots/segments; monotonicity;
tail alternatives; child/adult/sex separation; net-risk versus competing CIF;
conditional hazards; RR 1/3.6/25.3/74 and bounded probabilities; mortality and
depletion; baseline active-TB exclusions; uninfected/recent/remote extremes;
future first/reinfection events; reset/no-reset/protection/partial-reset scopes;
actual recency versus sensitivity clock; common-horizon coverage; age/sex/group
decompositions; extreme-risk concentration; explicit-RNG continuous event times;
canonical hashes/context/row-order invariance; source-preserving 3A adaptation;
unresolved direct time support and 100+; all twelve offline schema families;
validation-only comparison; blocked inversion; no production imports; and live
protected-state/frozen-input integrity.

The final `tested_files_final.json` checksum receipt covers 55 source/schema/test files. It is written before
the combined run and rechecked afterwards. Test XML and earlier development
receipts are retained in `abm/output/country_tb_m3b_audit_20261029/`.

Validation corrected a new-adapter floating-point attribution edge by preserving
exact 3A rational ratios. Actual event recency was separated from the no-reset
sensitivity clock, and common-horizon coverage was added so shrinking coverage
cannot be mistaken for decreasing cumulative disease in a fixed cohort. No
retained scientific artifact was repaired or replaced.
Final review also prevents a pediatric or finite-age reviewed RR from continuing
beyond its declared scope; that follow-up remains unresolved rather than falling
back to an inherited effect.

## Exports and numerical/stochastic checks

Command:

```
python -m engine.country_tb.preview_progression --output abm/output/country_tb_m3b_preview_20261029 --mc-samples 4000
```

Final immutable bundle:
`abm/output/country_tb_m3b_preview_20261029/progression-diagnostics-585539a943d7/`;
preview runtime **20.385 seconds**. The [diagnostic report](milestone3b_examples.md)
contains results and limitations. Scientific hash:
`585539a943d7c1aee57d520aa907734abdce2e0112ac46db3711391aaf8188c4`.

With seed 31829, 4,000 continuous-event samples, five-year reinfection and
synthetic competing mortality, the predeclared checks passed:

| Quantity | Deterministic | Monte Carlo | MC SE | Allowed absolute difference |
| --- | ---: | ---: | ---: | ---: |
| First-TB probability | 0.0128609093 | 0.01375 | 0.00178154 | 0.00900769 |
| Competing-death probability | 0.0245522291 | 0.0265 | 0.00244691 | 0.01233454 |
| Disease-free survival | 0.9625868617 | 0.95975 | 0.00300056 | 0.01510279 |
| Mean disease-free years | 4.910445812 | 4.898798499 | 0.00901234 | 0.04556170 |

Step refinement from 0.125 to 0.0625 years passed at every horizon; the largest
absolute risk change is 0.000000155635 versus the declared 0.0001 tolerance.
Numerical refinement is separate from Monte Carlo or source uncertainty.

`python -m abm.output.country_tb_m3b_audit_20261029.verify_export` passed in
**1.520 seconds**: all twelve contract families, 61 unique new forward hashes,
121 total embedded scientific hashes, all file checksums, input-cohort references,
population/depletion/decomposition/common-horizon checks, compilation of eight
Python files and 65 local document links. The detailed receipt is
`export_verification.json`. No national estimate or inverse fit is claimed.

Full-country computational scaling of continuous age/time quadrature is not
benchmarked here. This milestone validates bounded diagnostic populations; it
does not claim production performance or approval of numerical/evidence defaults.

## Protected state and change scope

The [protection fingerprint](protected_fingerprints_3b.json) records all protected
worktree/ref comparisons. Four protected siblings retain their HEAD, branch,
status, tracked/index diffs and changed/untracked bytes. The explicit-TBI worktree
retains exactly six modified and two untracked files. All 128 frozen files and
all non-target refs remain unchanged.

Frozen SA Health branch and peeled-tag commit:
`03cc16e4a52a55e10019dab16bd4f599571c4b90`; annotated-tag object:
`e99d36725dfc3ff40edc0ed68127e0b6cd927da5`. Remote values are checked again after
the authorized feature-branch push. No APY, explicit-TBI, Starsim, MATLAB,
deployed app, intervention, economics, DALY/QALY, production navigation or frozen
release file is changed. No merge, tag or deployment occurs.

Included files (18):

- `engine/country_tb/progression_common.py`
- `engine/country_tb/progression_natural_history.py`
- `engine/country_tb/progression_policies.py`
- `engine/country_tb/progression_baseline.py`
- `engine/country_tb/progression_forward.py`
- `engine/country_tb/progression_validation.py`
- `engine/country_tb/preview_progression.py`
- `tests/country_tb/test_progression_forward.py`
- `schemas/country_tb/v1/progression_forward.schema.json`
- `docs/country_subnational_model/milestone3b_evidence.md`
- `docs/country_subnational_model/milestone3b_forward.md`
- `docs/country_subnational_model/milestone3b_contracts.md`
- `docs/country_subnational_model/milestone3b_examples.md`
- `docs/country_subnational_model/milestone3b_validation.md`
- `docs/country_subnational_model/decision_register_3b.md`
- `docs/country_subnational_model/protected_fingerprints_3b.json`
- `docs/country_subnational_model/README.md`
- `docs/country_subnational_model/MILESTONE_STATUS.md`

Commit/push and post-push receipts are retained under the audit output directory.
The containing commit identifies the report without embedding a self-referential
commit ID inside it.
