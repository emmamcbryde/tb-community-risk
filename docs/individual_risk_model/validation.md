# Milestone 1 validation and protection

## Starting state

Worktree: `tb-individual-risk-platform`; branch
`feature/individual-tb-risk-platform-v1`; HEAD and origin both
`03cc16e4a52a55e10019dab16bd4f599571c4b90`; upstream
`origin/feature/individual-tb-risk-platform-v1`; ahead/behind 0/0; initially clean.
The frozen tag `sa-health-apy-he-v1.0.0` peeled to that starting commit.
Its annotated object remained `e99d36725dfc3ff40edc0ed68127e0b6cd927da5`.
No previous files were changed and no partial work was discarded.

## Exact bounded commands and results

```powershell
python -m pytest tests/test_individual_risk_m1.py -q --durations=8
python examples/individual_risk/run_example.py
python scripts/validate_individual_risk_m1.py
```

Final focused suite: **69 passed in 2.39 seconds**, 3.1406 seconds including
interpreter startup. Example execution: **0.3558 seconds**. The validator compiles
all **16 Python files** without importing or executing old workflows, captures
command output and SHA-256 checksums, checks that tested source bytes did not
change, and verifies example records and output scientific hashes.
See [machine receipt](validation_receipt.json). Generated complete inputs,
rankings under all three policies, command logs and receipts are under
`abm/output/individual_risk_m1/`. No full legacy or MATLAB suite was run.
No third-party dependency changes or interface changes were necessary.

Coverage includes all six input states, unknown versus false, explicit independent
and joint priors, bounds without central priors, absent bounds despite available
central priors, conditioning joint supports on known evidence, no hidden RR1
fallback for external-required effect models, assay distinction/compatibility,
TST interpretation, hazard probabilities, saturation/overflow, adult/child
separation, unsupported tails, treatment context, optional competing death and
first-event accounting, three RR policies and ties, missingness contributions,
ranking/overlap/common horizons, invalid and duplicate rows, unknown/numeric
identifiers, heterogeneous CSV presence, DataFrame missingness, full-precision
exports, source preservation, deterministic hashes, no prediction network and
frozen-file isolation.

## Independent worked arithmetic

Using 50-digit Decimal exponentials, posterior infection 2/3 and interval
cumulative hazard .04 give:

| Policy for effects 2 and 3 | Multiplier | Unconditional absolute risk |
| --- | ---: | ---: |
| Highest applicable | 3 | .075386375521894989648289847486775098480607611895160 |
| Additive excess | 4 | .098570807355859107695768679020958226817219782700380 |
| Multiplicative | 6 | .14224809262229772718727683498957484101786608428309 |

The independent decimal reference agrees with the binary64 forward calculations
within floating-point tolerance. Missing fictional diabetes/smoking marginalisation
gives .04443032971028184 by averaging risks over RR1/.56, RR2/.14 and RR3/.30.
This differs from applying the average RR1.74. The joint-prior test separately
verifies infection-weighted conditional progression, not a product of unweighted
averages. CSV and JSON exports preserve the exact binary64 results and original
source strings; no display rounding is used.

There are 12 fictional people: 10 central results, one child without a pediatric
curve and one prior-TB record explicitly unavailable. All numerical example
parameters and priors are invented. No clinical estimate is implied.

## Defects and contract checks

End-to-end validation exposed heterogeneous CSV rows losing the distinction
between omitted keys and explicitly unknown values. The source-field presence
metadata now preserves both, with invalid masks rejected. Identifier type metadata
preserves numeric, textual and explicitly unknown IDs without inventing IDs.
These cases are regression tested. Overflow paths and evidence/age applicability
are fail-closed, with no numerical cap or automatic adult fallback.

## Protected state and inventory

The [starting fingerprint](protected_start.json) covers **402 existing files**,
all refs and six other worktrees, including their changed/untracked file bytes.
The final bounded validator confirms these remain identical; only the requested
feature local/tracking refs may advance on commit/push. Existing Streamlit entries
and pages have no import of the new package. The frozen release branch and tag
are unchanged. No APY, MATLAB, Starsim, dynamic model, intervention, economics,
DALY/QALY or deployed workflow is modified. No merge, tag or deployment is done.

New files only: 12 modules in `tb_individual_risk/`; two example modules;
`tests/test_individual_risk_m1.py`; `scripts/validate_individual_risk_m1.py`;
and nine review/evidence/protection/validation files in this directory.
Post-push local/upstream/origin agreement and clean status are checked separately.

## Remaining scientific work

Review target populations and assay calibration, age/time-specific progression,
late hazards, risk-effect definitions and effect measures, pediatric evidence,
prior TB/preventive treatment, mortality and population priors. Validate calibration,
discrimination and fairness in independent cohorts before any clinical use.
Missing-information scenarios are not confidence intervals and rankings do not
establish care priority. Population synthesis remains a design contract only.
