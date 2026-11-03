# Individual TB-risk core: Milestone 1

An independent Python research library for calculating and ranking risk from a
supplied model definition. This is for planning and sequencing, not denying care.
**No clinical parameter set is supplied or approved.** This is not a diagnostic
service, an externally validated model, or an exact reproduction of TSTin3D.
No interface, deployment or connection to existing Streamlit workflows is added.

## Quick start

From the repository root, with Python 3.10 or newer:

```python
from tb_individual_risk import IndividualTBRiskModel
from examples.individual_risk.synthetic_parameters import synthetic_parameters, fictional_people

# Every numerical parameter below is INVENTED, NOT FOR CLINICAL USE.
model = IndividualTBRiskModel(synthetic_parameters())
result = model.predict(fictional_people()[0])
ranked = model.rank(fictional_people())
comparison = model.compare_policies(fictional_people(), policies=[
    "highest_applicable_rr", "additive_excess_rr", "multiplicative_rr",
])
```

`IndividualTBRiskModel()` without parameters returns an explicit unavailable
result, not a zero-risk placeholder. `ModelParameters`, `Parameter`,
`PopulationPrior`, `PersonRecord` and `Observation` are public Python contracts.
The numerical core uses only the standard library. A pandas DataFrame can be
passed to `rank`; pandas is not imported during ordinary dictionary prediction.

```powershell
python -m pytest tests/test_individual_risk_m1.py -q --durations=8
python examples/individual_risk/run_example.py
python scripts/validate_individual_risk_m1.py
```

Generated CSV/JSON, full logs and receipts go to `abm/output/individual_risk_m1/`.
The example contains 12 fictional people, including missing fields, an unresolved
child and a person with prior TB. Every result retains its synthetic label and
`clinical_use_approved=false`. No clinical effects are loaded from the example
by the production API.

See [science and evidence](science.md), the [parameter evidence table](evidence.csv),
[contracts and missing information](contracts.md), [ranking and worked examples](worked_examples.md),
[population design](population_design.md) and [validation](validation.md).

## Files and boundaries

- `tb_individual_risk/`: explicit contracts, infection probability, progression,
  RR policies, missing data, ranking, local I/O and provenance.
- `examples/individual_risk/`: invented demonstration definitions and runner.
- `tests/test_individual_risk_m1.py`: focused numerical, contract, I/O and isolation tests.
- `scripts/validate_individual_risk_m1.py`: bounded reproducible audit.
- This directory: equations, evidence, decisions, starting protection fingerprint
  and final validation receipt.

Existing release files, pages, MATLAB, APY, Starsim, dynamic models and deployed
entries are not edited. There are no treatment strategies, economics or DALYs.
First disease after the observation date is the only event represented. There
is no future infection/reinfection generation, recurrence or treatment-effect
model in Milestone 1. Prior active TB, preventive treatment and current active TB
require a different model and return unavailable when present.

Next: independent scientific review and acquisition of compatible infection,
assay, age/time-specific progression and conditional risk-effect parameters.
Only then assess external cohort calibration, discrimination, fairness and
transportability. Population generation and an interface remain later work.
