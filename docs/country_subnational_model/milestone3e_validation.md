# Milestone 3E validation and recheck

## Git state and scope

Started clean at `eebbe0242cf7787cf0f2f00e0172070bff99b8c6` on
`feature/country-subnational-tb-model-v1`, tracking the same origin branch,0/0.
Read-only remote inspection agreed before edits. Four protected worktrees,
frozen SA Health branch/tag,128 frozen files and 79 retained 2A–3D files matched
their prior fingerprints. There was no surviving partial 3E work to overwrite.

The review entry and helpers are new. Existing deployed entries, root navigation,
APY, explicit-TBI, Starsim/Summer, MATLAB, interventions, economics and frozen artifacts
are unchanged. Only the country review README/status are replaced/updated.
No merge, tag or deployment is authorized or performed.

## Bounded tests

Final repository-pinned Streamlit 1.39.0 run:

```powershell
python -m pytest tests/country_tb/test_review_ui_model.py tests/country_tb/test_review_ui_apptest.py -q --durations=8 -p no:tmpdir
```

**75 passed in 64.90 seconds**; 65.8646 seconds including interpreter startup.
This includes 59 pure model/contract cases and 16 rendered AppTest cases, covering
all required interface boundaries: isolated navigation/actual page discovery,
all user-TBI formats, CSV/source preservation, mixtures and unresolved groups,
no automatic recency, assay distinction, incidence default/direct-navigation gates,
all 15 unresolved-gate refusals, guarded explicit synthetic recovery, PRK, children,
no-tail semantics, highest RR, RR 74 warning, association marginal preservation/phi,
staleness, lazy full-precision exports, rounding, scenario separation and retained
source/artifact integrity. Tiny positive unresolved coverage never becomes a
complete-population rate.

An earlier installed Streamlit 1.32 run passed 71 cases in 44.73 seconds, including
the seven-page workflow and both restore actions. The final new cases add
scenario/tiny-coverage/subgroup-callback/mortality regressions. The 1.39 runtime
was installed **only** into ignored `abm/output/country_tb_m3e_runtime_139`
with `--no-deps --target`; shared Python and repository requirements are unchanged.
For reproducing that pinned test copy in PowerShell:

```powershell
$env:PYTHONPATH=(Resolve-Path abm/output/country_tb_m3e_runtime_139).Path
python -m pytest tests/country_tb/test_review_ui_model.py tests/country_tb/test_review_ui_apptest.py -q --durations=8 -p no:tmpdir
```

The accepted full run reports 360 existing third-party deprecation warnings from
Altair/jsonschema.RefResolver and protobuf datetime handling. They are not model
failures. Inherited Streamlit logging/CORS configuration may emit bare-mode
warnings; the documented local launch explicitly retains CORS/XSRF protection.
No security or shared application configuration is weakened.

The small retained user-TBI group is separately checked:

```powershell
python -m pytest tests/country_tb/test_baseline_sources.py -q --durations=4 -p no:tmpdir -k "no_fallback or overall_known_scale or age_specific_one_scale or direct_profile_original or direct_recent_remote or subgroup_metadata or incomplete_subgroup or open_tail_direct_value or explicit_conditional_recency"
```

Exact timing/result is in [validation_receipt_3 e.json](validation_receipt_3e.json).
No 337-test, multi-hour or full-country suite is run. Scientific-core bytes are
also checked against the starting 79-file fingerprint.

## Rendered verification and performance

```powershell
python scripts/validate_country_tb_review_ui.py --output abm/output/country_tb_m3e_rendered --prospective
```

Final pinned 1.39 run completed in 49.7545 seconds. It saved immutable AppTest
headers, metrics, warnings, table columns/counts and timings for **all seven views**,
plus full synthetic instantaneous/prospective results. Evidence directory:
`abm/output/country_tb_m3e_rendered/review-ui-20261102T020430996404Z/`.
These are rendered-test records, **not browser screenshots**.

| Operation | Final measured wall time |
| --- | ---: |
|First app load, including imports/initial source view|3.2460 s|
|First Viet Nam UN WPP population load/validation|30.4856 s|
|Working-example population load|.0015 s|
|User-supplied TBI catalytic construction|.3005 s|
|Risk allocation|.2815 s|
|Deterministic instantaneous projection|.0290 s|
|Results-page rendering|.3168 s|
|Default incidence-gate refusal|.000072 s|
|Optional one-year forward operator|10.7115 s|
|Optional one-year complete run|12.6890 s|

Cold WPP cost comes from validating the complete frozen 500677-row source via
the retained loader. It is an explicit, cached load; ordinary edits do not eagerly
repeat it. The default ABS source loads much faster. Download preparation is
explicit and tested not to execute during ordinary results rendering.

A focused mortality performance check initially took 60.8865 seconds with all 200
finite ABS intervals. Exact scoping to the cohort's 8 represented annual age/sex
intervals reduced the regression run to 14.73 seconds, with the identical expected
case value to 1 e-12. This removes unused lookup rows, not people, ages within the
represented horizon, or source information. The full immutable snapshot hash,
licence, source count and interval selection remain in provenance. No 3B operator
or source value was modified.

Optional incidence inversion has **no empirical accepted timing** because no
country configuration is approved. Default refusal is timed above. A complete
explicit uploaded synthetic configuration is exercised in the model tests,
including root uniqueness/refinement and known-scale recovery; it is not a
performance claim for national inversion.

## Clearly classified examples

The restored working cohort is **synthetic/user supplied**:1000 adults in explicit
female/male single-year cells at 20 and 60,20%overall TBI and RR 1. Its age/recency
and infection-time profiles are **model derived**, conditional on an explicitly
selected constant diagnostic pressure history. Results:

- Total TBI 20.00%; recent 2.79%; remote-only 17.21%; prior-remote-plus-recent.49%
  is within recent, not added again.
- Instantaneous first-episode rate 48.1849657735134/100000 PY, case flow
  .481849657735134/year at baseline.
- One-year **diagnostic sensitivity**:.3719675681701176 expected cases;
  incidence 37.204350154859085/100000 disease-freePY; 999.6280324318299 surviving
  disease-free people under the no-mortality diagnostic.
- Explicit ABS mortality sensitivity gives.37144663879619993 one-year cases.
  It is not an adopted mortality model or Australian estimate.
- The invented 700/200/100 subgroup mixture with 5%/25%/15%anchors yields 10.00%
  population-weighted total TBI. Missing 200-person subgroup input makes the
  complete mixture unavailable while retaining the full 1000 denominator.

**Observed** WHO incidence is shown separately for context. Its annual national
new/recurrent episode definition is not the custom cohort's instantaneous
first-episode rate.3D held-out RMSEs are displayed as internal empirical model
checks:8.25,2.71,1.18 and 22.14 percentage points. No fictional Viet Nam result is
used as an empirical country estimate.

## Browser availability and limitations

Live browser verification was **not performed**. The browser access host was
unavailable. Automatic approval review rejected binding the only selected
Chrome window because it contained an unrelated document and could expose
private contents. No document interaction or workaround was performed. A
dedicated localhost window was requested; AppTest work continued independently.
The rendered receipt explicitly records `browser_verified=false`.

Other limits:2024 baseline only; bounded 5000-row CSV/300-stratum inverse and
50000-node prospective requests; first cold WPP load cost; adult/pediatric/tail,
effect, assay, mortality and observation transport not scientifically approved;
no stochastic interface or long-horizon country projection. A reviewed direct
100+TBI value does not invent its exact ages or infection times. Production
incidence inference remains blocked.

## Defects corrected and file inventory

The defects were in new UI/adapter code, not the frozen scientific operator:

1. Streamlit 1.32/AppTest selectbox formatting used raw IDs against formatted
   options; display labels are now real option values with an explicit ID map.
2. Restore-button `st.rerun()` could reconsume the pressed button state; callbacks
   now update state before rendering, without an infinite rerun loop.
3. CSV baseline years needed explicit integer normalization, while preserving raw
   source rows. Percentage conversion now preserves full input precision rather
   than inheriting Decimal's default 28-digit rounding.
4. Partial subgroup values, scope/count/risk/association mismatches and tiny
   positive unresolved coverage are fail-closed. Mortality interval lookup is
   scoped exactly without changing the result.

New files:15 Python entry/helper/view files under `country_model_review_ui/`, its
README and 5 synthetic CSV inputs; `scripts/validate_country_tb_review_ui.py`;
2 focused test modules; `schemas/country_tb/v1/review_ui.schema.json`;
this report, [architecture](milestone3e_ui.md), [decisions](decision_register_3e.md),
and protection/validation receipts. Existing country README/status are updated;
status is replaced with the completed 3E state, not appended as a partial summary.

Quick checks passed in 6.55 seconds: 18 Python files compiled, all nine UI
contract families validated offline, 824 complete scientific-contract hash
occurrences checked (64 distinct hashes), and 81 local documentation links
resolved. The additional synthetic inversion recovered 0.01699999999837018
from a known scale of 0.017. Hash-only provenance references were distinguished
from complete contracts when recomputing scientific hashes. All tested source
bytes still match the successful test receipts; only documentation changed.

Python compilation, offline schemas/scientific hashes, local documentation links,
`git diff --check`, explicit staged inventory and protected-state checks are
recorded in the machine receipt. Post-push local/upstream/origin agreement and a
clean worktree complete release verification. Scientific approval remains pending.
