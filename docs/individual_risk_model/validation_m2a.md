# Milestone 2A validation and protection

Started clean on `feature/individual-tb-risk-platform-v1` at
`8b533f81a6a1435bed11f10db52c6d8d916224a1`, upstream the same origin branch,
ahead/behind0/0. Read-only remote inspection agreed. The frozen SA Health release
branch and peeled tag remained at `03cc16e4a52a55e10019dab16bd4f599571c4b90`;
annotated tag object `e99d36725dfc3ff40edc0ed68127e0b6cd927da5`.

## Bounded commands

```powershell
python -m pytest tests/test_individual_risk_m1.py tests/test_individual_risk_m2a_audit.py tests/test_individual_risk_m2a_reference.py -q --durations=8
python scripts/reproduce_menzies_2008.py
python scripts/validate_individual_risk_m2a.py
```

Final suite: **109 passed in 3.50 seconds**, 4.2535 seconds including startup. This
contains 69 retained M1 cases, 7 independent audit regressions and 33 historical/
package/reference cases. The reproduction command completed in **2.2465 seconds**.
The validator compiles **23 Python files**, checks unchanged tested source bytes,
all four source-data JSON objects, package pins and output scientific hashes,
and runs protection and Git whitespace checks. It verified 496 complete scientific
hash occurrences, 186 unique hashes. No network is accessed during prediction;
socket-denial tests cover both formulations. No unrelated Streamlit/MATLAB suite
was run. Exact commands/log output, file checksums and timings are in the
[machine receipt](validation_receipt_m2a.json).

The initial audit deliberately exposed five failing M1 cases before correction;
subsequent focused runs passed. A new historical-input validation error applied
probability bounds to incidence, pack/day and percent-weight units; it was corrected
and regression tested. No frozen engine defect was patched or dependency upgraded.

## Reference evidence and limitations

The paper was downloaded only for development inspection into ignored output.
Relevant complete pages499-503 were rendered with an isolated workspace copy of
PyMuPDF1.26.7 because Poppler and PDF readers were absent. The shared Python
installation and requirements are unchanged; PDF tools are not prediction dependencies.
No copyrighted PDF/screenshots or website scripts are committed. Source URLs,
access dates and checksums are in [source receipts](source_receipts_m2a.json).
No clinical website submissions or author messages occurred.

There are26 fictional reference people. All40 Table 3 cells and two Figure 2 disease
outputs are represented. No whole-algorithm reproduction is claimed.26 partial
progression-stage checks use printed PPV as an input;24 fall within one declared
printed-resolution unit and two differ.16 remain unavailable because the general
combination/contact-lifetime rules are not verified. These are reproduction
outcomes, not failing unit tests hidden or forced to match.

The [reproduction report](reproduction_m2a.md) gives every value, unresolved cause
and partial difference; generated full-precision JSON/CSV also retain relative
errors. The [comparison](formulation_comparison_m2a.md) demonstrates invalid raw
linear probabilities and a bounded exponential mathematical sensitivity. Actual
modern fixtures keep the independent synthetic package and are not tuned.

## Scope and integrity

All402 frozen existing-release files match their original fingerprints. Of427 tracked files at the start,418 remain byte-identical; only the nine
explicitly authorized core/README files were edited. All six other worktrees retain HEAD, status,
index/working diffs and changed/untracked file bytes, including explicit-TBI's
six modified and two untracked files. All refs except the requested feature
local/tracking refs remain unchanged. See [starting fingerprint](protected_start_m2a.json).

The commit contains isolated library/API audit corrections, new reference/package/
comparison modules, four immutable source-data JSON objects, fictional fixtures,
two focused test modules, reproduction/validation scripts and review documents.
No existing deployed entry, pages, APY, frozen artifacts, MATLAB, Starsim, country
model, economics or intervention logic is changed. There is no new interface,
country data repository, population generator or transmission model.

Final staged inventory and `git diff --check` are verified before committing.
Post-push local/upstream/origin agreement,0/0 and clean status are recorded
separately to avoid a self-referential commit receipt. No merge, tag or deployment.

## Review decisions

Potentially reusable: the documented single-factor historical arithmetic and
strictly separated modern hazard engine. Historical-only: baseline/risk effects,
contact assumptions and survival-to80 bookkeeping. Unresolved: original lookup
and BCG mappings, concurrent effects, early-period lifetime integration, IGRA,
serial testing, treatment benefit, pediatric and contemporary clinical parameters.
Current public documentation does not establish an exact Version 3.0 implementation.
Author clarification, reuse terms and independent scientific review are next.

Final visual inspection also identified the smoking threshold disagreement
(Table 2 >=1 versus Figure 1 >1 pack/day) and exact two-year recency boundary.
Both are documented and blocked at the ambiguous boundary. These refinements
were made before the reference package's first commit; final package pins and
tests reflect the final source.
