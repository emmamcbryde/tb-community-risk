# Milestone 2B validation and Git protection

Date: 2026-10-28. Starting commit
a4fceee5af7c2871943061ad9635084fd18cec89 on
feature/country-subnational-tb-model-v1. Directory, branch, upstream and remote
HEAD agreed before any tracked edits; worktree clean, ahead/behind 0/0.

## Executed checks

Final focused run: **190 passed, no skips, in 361.25 seconds (6:01)**.
This comprises the previous 139 focused country/WHO/profile tests plus 51
Milestone 2B cases. Two openpyxl dependency warnings concern deprecated
datetime.utcnow during read-only workbook loading; no source workbook is saved.

The process set tempfile.tempdir to the short workspace-local
abm/output/m2b_tmp and invoked:

    pytest.main(["-q", "tests/country_tb",
                 "tests/test_general_who_incidence.py",
                 "tests/test_general_population_profile.py",
                 "-p", "no:cacheprovider", "--tb=short"])

The process ran with sandbox escalation to avoid the documented Windows
temporary-directory ACL/path restrictions. No multi-hour full suite was run.
An earlier pure-mathematics subset passed 41 tests in 0.98 seconds; the final
run includes those cases and all source/protection checks.

Coverage includes:

- Independent native ABS/WPP reconstruction, including reselection from the
  full locally preserved WPP original archive; file/manifest and scientific hashes.
- Exact 1,212 excess-row decomposition, full stored-field identity, row-order
  independence, idempotence and rejection of conflicting keys.
- Frozen WHO rows/entities/years, ordered uncertainty, missingness, distinct
  territories, original-snapshot restoration with catalogue lookup forbidden.
- Female/male counts, audit-only published both, exact display totals,
  source estimate/projection labels and existing post-2024 Medium selection.
- Zero/constant/step/rising/declining pressure, fractional and old ages,
  absent/gapped/overlapping history, first-event recency, invalid targets,
  one-scale synthetic anchor recovery and exact decimal reconciliation.
- Explicit RNG categorical checks (100,000 draws, seed 20241028), predeclared
  per-state tolerance 6*sqrt(p*(1-p)/n)+1/n, repeatability and invalid RNG failure.
- Proposed schema examples, canonical serialization/hashing and material
  changes, offline ordinary loading, source/frozen/protected isolation.

Python compilation passed for engine/country_tb, tests/country_tb and the two
new developer scripts. Git whitespace checks and local documentation links
passed; staged checks are repeated before commit.

The first independent checker draft expected a generic missing label; the
source uses not_estimated. That checker assertion failed and was corrected
to the actual source vocabulary, after which the independent reconstruction
and final tests passed. No frozen source was changed to satisfy a check.
No numerical or classification discrepancy from the 2A claims was found;
the small source both-sex rounding differences are explicitly reported.

## Reproducible diagnostic outputs

    python scripts/recheck_country_tb_frozen.py --output abm/output/country_tb_m2b_audit/independent_recheck.json
    python scripts/preview_country_tb_catalytic.py --output abm/output/country_tb_m2b_audit/synthetic_backcasts.csv

The recheck is a source-artifact audit, not a network download or WHO
estimation rerun. The synthetic CSV contains invented histories only. Its
receipt under the same output directory binds the exact script bytes,
explicit input description, canonical rows and scientific hash:

| Diagnostic | SHA-256 |
| --- | --- |
| Synthetic backcast CSV bytes | b2e08354ce76ead47a152304e26e8de168283748d2f3de994464e9faf690baa8 |
| Canonical synthetic scientific receipt | 9dcf9ad51352e1dc9a835417cce32a3f6df9c71cc0a6051fdbe59a8d3fe9e4b6 |

The scientific hash differs from the file-byte checksum by design. Numeric
closed-form examples use Python floating arithmetic and are independently
checked against fine quadrature; the catalytic probability prototype uses
Decimal with the explicitly documented rounding rule. These are software
validation artifacts, not fitted country outputs. The preserved 2A country
preview profiles retain their original scientific hashes.

## Protected content and refs

[Receipt](protected_fingerprints_2b.json) records before/after branch, HEAD,
status, tracked/index diff hashes, every dirty/untracked file hash, all refs
before commit and 128 frozen tracked-file hashes.

| Worktree | Starting and preserved HEAD | Dirty content preserved |
| --- | --- | --- |
| Starsim | 5a4a39e31d01510831068ec128c098290acf8e39 | 39 untracked report files |
| Generic no-transmission | f834057144b3fa17a9de6d995a7997a93ec4f08e | Clean |
| Explicit recent/remote TBI | 7c000dd80eef8ec9ce797579f2f53736f75224e0 | Six modified tracked and two untracked Milestone 3H files |
| WHO sibling | 666088cac1e20dd1e9e52016ea857a9c48e0ba8d | Clean |

SA Health release branch: 03cc16e4a52a55e10019dab16bd4f599571c4b90.
Annotated tag object: e99d36725dfc3ff40edc0ed68127e0b6cd927da5,
peeling to the same release commit. Local and remote refs matched at start.
All other registered worktrees and refs remained unchanged before commit.

Only the requested country-model branch and its remote-tracking ref may
advance for this commit. A final post-push check compares protected content,
all other refs, frozen files and remote release refs again; delivery reports
the resulting commit and clean status without a self-referential commit hash
inside the committed receipt.

No shared 2A loader/schema/snapshot, MATLAB, frozen artifact, deployed entry
point, Starsim code, intervention or economics logic was changed. No country
calibration, TBI estimate, risk-factor allocation, static/dynamic projection,
production UI, merge, tag or deployment was executed.
