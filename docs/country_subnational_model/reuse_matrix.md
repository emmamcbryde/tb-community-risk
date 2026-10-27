# Existing model code reuse assessment

This is an assessment before reuse. No code was copied or cherry-picked. Other worktrees were read-only; committed source was inspected with Git so local unfinished changes did not enter the architecture.

| Source and commit | Useful component or pattern | Scientific or software mismatch | Decision |
| --- | --- | --- | --- |
| General community 24b6323 | engine/who_incidence contract, snapshot, import and trend modules | Disease trend is descriptive; loader can fall back to fixture; not infection calibration | Candidate offline data-boundary reuse with explicit production checks |
| General community 24b6323 | engine/profiles/country.py and local_incidence.py | Country selection does not build country demography or infection history | Reuse preview/provenance/units concepts after independent contracts |
| General community 24b6323 | population_profile.py missing/value/excluded states and provenance | Age/risk representation and legacy mapper cannot express the new joint distribution | Concepts only; no silent schema migration |
| General community 24b6323 | engine/profiles/engine_mapping.py | Uniform user prevalence across ages; RR, HR and OR all applied as hazard multipliers; descriptive incidence only | Must not import into the new architecture |
| General community 24b6323 | engine/apy/infection_history.py, cohort.py and calibration.py | APY anchors, legacy risk products and community reference; calibration does not identify country history from WHO disease incidence | Audit mathematical lessons only; no reference assumptions/defaults |
| Explicit recent/remote 7c000dd | engine/apy/explicit_recent_remote_tbi.py: mutually exclusive states, feasibility, infection timing and deterministic/stochastic consistency | Constant recent/remote hazards, 5-year window, capped 100-year history and default oldest age 89; TBI targets differ from infection-trend reconstruction | Candidate math/test patterns after rederivation; no import or copy |
| Explicit recent/remote 7c000dd | population_tbi_progression_risk.py: separate validated scenario multiplier | Risk-enriched multiplier 8 is a scenario default for that workstream, not country evidence | Validation pattern only; never inherit multiplier |
| Generic no transmission f834057 | engine/profiles/background_exposure.py and background_exposure_hazard.py | Strict exogenous hazard units, half-open time intervals and missingness are useful; first infection only, no reinfection, age schedules not executable | Candidate boundary/math concepts; not a complete baseline generator |
| Starsim prototype 5a4a39e | engine/starsim_tb parameters, analytics, prototype; feasibility report | Four-state demonstration; active TB absorbing; no mortality, migration, reinfection, HIV or reviewed country calibration | Future adapter feasibility reference only; no integration or parameter reuse |
| Legacy static/dynamic UI | Country text inputs/cache ages and plots | Uniform fallback, legacy incidence-to-ARI/backcasts and incompatible baseline logic | Do not reuse scientific initialization |
| engine/infection_backcast.py | Lifetime/recent/remote formulas | Direct Styblo-style incidence conversion, probability clamping and implicit earliest-value backfill; clearance calculation is not an approved loss-of-infection model | Must not reuse calibration/backcast behavior |
| ui/dynamic_age_distribution.py | Five-year/single-year transformations | Coercion/clipping, uniform fallback, uniform within-band ages and fixed age 100 for 100+ | Concepts only; no validation or tail defaults |
| Frozen SA Health release 03cc16e | Existing release/integrity guard patterns | Frozen APY population, reference values, economics and validated deployment scope | Must remain unchanged; no country defaults or engine wiring |
| MATLAB v9 and frozen validation/reference reports | Reference artifacts for their own model only | New architecture must not inherit scientific compatibility constraints as country evidence | Protected; no edits or new model versions |
| Existing deployed Streamlit entry points and runners | Existing service workflows | New experimental architecture must stay isolated | No imports or edits |

## Worktree observations

The generic worktree was clean. The explicit-TBI worktree already had six modified tracked files and two untracked files; the Starsim worktree already had 39 untracked report files. Those are user state and were not repaired or included. The audits use the commits above. Protection verification compares HEAD, tracked diff hashes, untracked-content hashes and local branch/tag refs, alongside the permitted changed-file allowlist.

## Reuse prerequisites

Before selective reuse, rederive units, infection-event semantics, age/sex denominators, natural-history reference populations and boundary behavior. Porting an isolated formula is a later explicit implementation choice with tests; a shared package import can pull in defaults and therefore needs dependency review. In particular, no APY odds-ratio product, adult medical algorithm in children, fixed 0.1% progression baseline, unreviewed Styblo conversion or demonstration Starsim active prevalence is permitted.
