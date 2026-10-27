# Country and subnational TB milestone status

Updated 2026-10-27. This file is the replaceable current milestone record; later milestones should replace its contents while retaining historical records in Git.

## Current milestone

Milestone 1 architecture, evidence and data boundaries is complete for review. Scientific assumptions are proposed and require review; no country model has been calibrated or approved for analysis.

All ten requested deliverables are indexed in [README.md](README.md). Seven proposed contracts share version 1.0.0 and schema IDs country_tb.geography.v1, country_tb.incidence_history.v1, country_tb.demography.v1, country_tb.tbi_profile.v1, country_tb.risk_factor_prevalence.v1, country_tb.tbi_risk_association.v1 and country_tb.baseline_population.v1.

The only executable additions are isolated validation, canonical hashing, age reconciliation and Fréchet/phi mathematics. They do not construct or project a population. The seven structural examples use invented data and placeholder hashes, with static and dynamic readiness both false.

## Evidence and data decisions

Use the already frozen national disease-incidence snapshot who-gtb2025-incidence-6997e3b087d2, report round 2025, estimate years 2000–2024, downloaded 2026-09-24. It contains 217 countries/areas and is distinct from the eight-code test fixture. The future production loader must refuse fixture fallback and bind the expected identity.

The WHO sibling main branch at 666088cac1e20dd1e9e52016ea857a9c48e0ba8d was audited read-only, including binary incidence, age/sex, population and risk outputs. It is not a fully reproducible or licensed code dependency for this architecture.

Demography needs a new frozen reviewed source boundary. The legacy 2021 cache has 1,212 duplicate country/year/age keys, no sex column and insufficient source-vintage/oldest-age metadata. No demographic default was selected.

TSTin3D's eleven requested values remain candidate effects, sensitivity-only pending contemporary review. Adult values are not applied to children. Its present script sums selected medical values; that behavior is audited, not adopted. Proposed highest-applicable combination needs review and sensitivity policies.

## Validation

- New focused contract/math/isolation tests: 57 passed.
- Existing focused WHO snapshot/import tests: 18 passed.
- Commands: python -m pytest -q tests/country_tb tests/test_general_who_incidence.py -p no:cacheprovider.
- Combined result: 75 passed in 2.92 seconds with Python 3.12.4, pytest 7.4.4 and jsonschema 4.19.2. jsonschema is a test-only dependency; the new engine prototype uses only the standard library and frozen production requirements were not changed.
- Windows sandbox reruns encountered four temporary-directory permission errors in existing WHO tests. The successful combined run used the same tests outside the sandbox with process-local tempfile.tempdir set under abm/output/country_tb_m1_audit/test_tmp_elevated. No tests were weakened or changed to bypass those failures.
- Structural schema and all seven examples validate locally; source/provenance round-trips and canonical hash behavior are covered.
- New Python/JSON syntax, local-document links, selected WHO snapshot identity/coverage/checksum transcription and absence of pre-existing tracked-file changes passed. Before commit, all protected HEADs, tracked diffs, untracked-content hashes and all local branch/tag refs matched the starting receipt.

No full suite, simulation, MATLAB validation run, Starsim run or production Streamlit workflow was started. None is needed for these isolated Python/schema and documentation additions.

## Protected state

No pre-existing tracked file, application entry point, runner, MATLAB source, frozen artifact or WHO source was changed. No code was copied or cherry-picked. Other worktrees retain their starting HEADs and user changes; their initial dirtiness is documented in [git_audit.md](git_audit.md).

The only branch authorized to advance is feature/country-subnational-tb-model-v1. The delivery response records the resulting commit and push/clean-tree status without creating a self-referential commit-hash field in this file.

## Open scientific decisions and next milestone

Resolve pre-WHO backcasts; incidence-trend elasticity and discrepancy; child/adult and recency-dependent natural history; test-positive versus infection-state calibration; reinfection/clearance; historical HIV/medical factors; demographic revision, migration, mortality and oldest ages; medical effect transport/combination; multimorbidity constraints; active/prior-TB initialization; and source-rights questions.

Recommended Milestone 2 is frozen incidence/demography loaders and cross-artifact semantic validation, alongside the scientific calibration design. Screening, projection integration, economics and production pages remain outside that scope.
