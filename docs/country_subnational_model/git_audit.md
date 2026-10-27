# Country model Git audit

## Starting state

Audit date 2026-10-27. Required directory and branch matched before any work.

| Item | Starting observation |
| --- | --- |
| Directory | C:/Users/emmas/Documents/GITHUB/tb-community-risk-country-model |
| Branch | feature/country-subnational-tb-model-v1 |
| HEAD | 24b632390ce75e24c15088309e45960a990e66ad |
| Tracking upstream | origin/feature/country-subnational-tb-model-v1 |
| Ahead / behind | 0 / 0 |
| Working tree | Clean |
| Origin | https://github.com/emmamcbryde/tb-community-risk.git |
| Branch origin evidence | Reflog: Created from origin/feature/general-community-tb-model |

Starting HEAD equals both local feature/general-community-tb-model and origin/feature/general-community-tb-model. Local remote-tracking refs were inspected without fetching or modifying source branches.

## All registered worktrees

All are siblings under the starting directory's parent.

| Worktree basename | Branch | HEAD |
| --- | --- | --- |
| tb-community-risk-country-model | feature/country-subnational-tb-model-v1 | 24b632390ce75e24c15088309e45960a990e66ad |
| tb-community-risk | prototype/starsim-tb-spectrum | 5a4a39e31d01510831068ec128c098290acf8e39 |
| tb-community-risk-no-transmission | feature/generic-no-transmission-model | f834057144b3fa17a9de6d995a7997a93ec4f08e |
| tb-community-risk-tbi-calibration | feature/explicit-recent-remote-tbi-v2 | 7c000dd80eef8ec9ce797579f2f53736f75224e0 |

## Relevant starting branch refs

| Branch | Local HEAD | origin remote-tracking HEAD |
| --- | --- | --- |
| feature/country-subnational-tb-model-v1 | 24b6323 | 24b6323 |
| feature/general-community-tb-model | 24b6323 | 24b6323 |
| feature/explicit-recent-remote-tbi-v2 | 7c000dd | 7c000dd |
| feature/generic-no-transmission-model | f834057 | f834057 |
| prototype/starsim-tb-spectrum | 5a4a39e | 5a4a39e |
| release/sa-health-apy-he-v1.0.0 | 03cc16e | 03cc16e |
| main | eb209e6 | eb209e6 |

## Existing user state and protection verification

WHO and generic worktrees started clean. The explicit-TBI worktree had six tracked modifications and two untracked files. Its initial tracked diff SHA-256 was 0b145ea9ead01345ed8ec274785abb46f08659e5d4cc1dd65c11255f83f0f851. The Starsim worktree had 39 untracked report files and no tracked modifications.

An ignored audit receipt under abm/output/country_tb_m1_audit records HEAD, status, tracked-diff hash, untracked-content aggregate hash and all local branch/tag refs before milestone writes. End checks compare these values and permit only the requested country-model branch to move. No reset, clean, discard, overwrite of pre-existing files, checkout, merge, tag or deployment is authorized or performed.

The precommit comparison passed for all protected worktrees and local branch/tag refs. Read-only ls-remote also confirmed both the target branch and general-community remote branch at the starting commit immediately before the authorized commit/push. It did not fetch or move protected refs.

Final commit and push results are reported in the delivery response so the commit does not need to contain its own hash. MILESTONE_STATUS.md records the final reviewed scope, tests and remaining decisions; it is the file intended to be replaced at later milestone completion.
