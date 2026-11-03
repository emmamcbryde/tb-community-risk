# Contracts, input states and missing information

## Versioned objects

- `individual_tb.person_record.v1`: ID and all supported observations.
- `individual_tb.parameter.v1`: discriminated by `name`, source and review status;
  each named family is explicitly versioned within this contract.
- `individual_tb.population_prior.v1`: missing-input distributions and scenarios.
- `individual_tb.model_parameters.v1`: complete parameter assembly and policy.
- `individual_tb.prediction.v1`: estimates, source records, explanation and hashes.
- `individual_tb.synthetic_population_specification.v1`: design-only future input.

Parameter families: infection probability, TST interpretation, IGRA interpretation,
baseline progression hazard, recent progression, remote progression, age effects,
radiographic effects, medical effects, competing mortality and prediction horizons.
`ModelParameters.record()` includes canonical SHA-256 provenance and can be loaded
with `ModelParameters.from_dict()`. Hash mismatch is rejected. Model construction
copies and validates inputs so later caller mutation does not change predictions.
`external_parameter_required` means no numerical value is authorized by this package.
`reviewed`, `sensitivity_only` and `synthetic_not_clinical` are distinct statuses.
A synthetic object cannot be reclassified as a reviewed set. All outputs retain
`clinical_use_approved=false`; source status is not independent expert sign-off.

## Person fields

`record_id` preserves original nonempty strings (including leading zeros) and
integer identifiers. It also permits an explicit unknown/absent/excluded state;
no person ID is invented. Unknown IDs remain null with a source-row audit locator.
All other fields materialize as explicit unknown observations when omitted.
Numeric/categorical/boolean inputs can be supplied directly or as:

```python
{"state": "known_present", "value": True,
 "source": {"citation": "reviewer's source", "definition": "exact exposure"}}
```

States are `unknown`, `known_absent`, `known_present`, `not_applicable`, `excluded`
and `invalid`. False is known absent; None is unknown. A measured numerical zero
is a known value. Not-applicable and excluded fields do not become false, and
are not overwritten by missing-data priors. Required excluded fields make the
calculation unavailable. Explicit invalid observations return structured errors.

Core names are `current_age`, `sex`, `country_of_birth`, `current_country`,
`age_at_immigration`, `immigration_year`, `birth_year`, `close_contact`,
`recent_exposure`, `recent_conversion`, `time_since_infection`, `tst_result`,
`tst_mm`, `igra_result`, `igra_product`, `bcg`, `bcg_age`, `cxr`,
`previous_active_tb`, `previous_preventive_treatment`, `active_tb_at_baseline`,
`hiv`, `hiv_treatment`, `diabetes`, `smoking`, `underweight`, `silicosis`,
`haemodialysis`, `severe_renal_disease`, `tnf_inhibitor`, `transplantation`,
`immunosuppression`, `head_neck_carcinoma`, `observation_date`, `horizon_years`,
`tbi_probability` and `pretest_tbi_probability`. Additional factors use
`factors={"stable_factor_id": observation}`; parameter fields reference
`factor:stable_factor_id`. Wrong or misspelled top-level fields are errors.

TST result: positive/negative/not_done. IGRA additionally accepts indeterminate
and borderline, which need an observation model. CXR: normal/granuloma/
fibronodular/other. Sex: female/male/other, or an explicit non-known state.
Dates are ISO YYYY-MM-DD. Known ages, durations, calendar years and immigration
metadata receive consistency checks. A time since infection cannot exceed age.
No immigration history or infection clock is invented from country of birth.

## Marginalisation

Priors are supplied explicitly, never hidden in code. They must be applicable to
the person/population and observed evidence; the library does not select a prior
from country, sex or ethnicity automatically. Each finite support row has a
value and probability, plus object-level source and status. Dependence policy
must be explicit: `explicit_independence` for marginal distributions, or
`explicit_joint` with complete support rows. The joint table is conditioned on
known inputs and normalized; zero-support configurations remain unavailable.
Marginal distributions alone do not identify dependence.

Central absolute risk is the finite sum `sum(w_s * p_s * r_s)`. It is **not**
computed from an average RR. Returned infection probability is `sum(w_s*p_s)`;
conditional progression is `sum(w_s*p_s*r_s)/sum(w_s*p_s)`. This preserves the
product identity even when infection and progression inputs are correlated.
When infection probability is exactly zero, unconditional reactivation risk is
zero under the explicit no-future-infection scope; progression can remain unknown.

Bounds are separately supplied finite states or joint-bound rows. The lower
and upper results are extrema over these prespecified missing-information
scenarios, not confidence intervals or universal continuous-domain guarantees.
Priors do not automatically become bounds. Scientific reviewers must justify
that supports bracket relevant possibilities. Missing distributions/scenarios
remain unavailable; the engine neither substitutes zero nor drops unsupported
scenario mass. Enumeration is exact and budgeted (8192 scenarios by default).
Exceeding the configured budget blocks the result rather than changing dependence.

Returned completeness counts known total and critical fields. Input contributions
retain known values/source, marginal state weights and conditional mean risks.
Factor explanations also report the absolute-risk difference if each included
factor were absent, with all other inputs held fixed. These differences do not
sum to total risk under nonlinear or highest-RR policies. Unused known evidence,
explicit scope exclusions, unsupported histories and unresolved inputs remain visible.

## Line lists and export

`rank` accepts dictionaries, lists of dictionaries or pandas DataFrames.
`read_csv(path)` and `read_json(path)` read local files only; URLs are not fetched.
CSV columns use person names. Empty cells are unknown; JSON values in cells
preserve booleans, structured states and source objects. IDs stay literal text.
`write_records` creates a lossless input pair. Its optional `__present_fields`
CSV metadata preserves omitted keys separately from explicitly unknown cells;
`__identifier_kind` preserves integer/text/unknown identifier types.
nonempty values outside that mask are rejected, not hidden. Duplicate CSV headers reject the
file; malformed/invalid rows and duplicate person IDs remain structured invalid
results, not silently dropped. Unavailable scientific estimates are distinct
from invalid records.

`export_results` writes full nested JSON and a CSV with scalar summaries plus
`full_result_json`. Numeric outputs are binary64, exported with Python's
round-trip precision, not presentation rounding or arbitrary-precision claims.
Original numeric strings and observation sources remain separate. Non-finite
source anomalies are preserved as tagged source markers, never scores. No raw
person records are sent to the network or written except through explicit exports.
