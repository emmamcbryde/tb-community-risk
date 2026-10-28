# Country and subnational TB milestone status

Updated 2026-10-28. Replaces the Milestone 2A status; earlier versions remain in Git.

## Current milestone

Milestone 2B independent frozen-data verification and scientific calibration
design is complete, ready for review. No country pressure calibration, retained
TBI estimator, risk-factor population, projection engine, intervention,
economics or production UI is connected.

## Independent data findings

- WHO 2025 snapshot: 5,347 records, 217 countries/areas, 2000-2024.
  Snapshot/manifest hashes, uncertainty and missingness checks agree.
- Legacy duplicates: 1,010 two-row age groups plus 101 three-row VCT groups:
  1,111 groups and 1,212 excess records, identical in every stored field.
  Canonical quarantine has 20,705 rows; original source remains untouched.
  Modern key dimensions and lost upstream lineage are distinguished explicitly.
- Native ABS workbook and native WPP fields independently reconstruct all
  303 and 502,677 transformed records. Selected WPP rows also match the
  preserved full original archive in this session.
- Australia: ABS June 2024 estimate, 27,204,809. Viet Nam: WPP 2024 Medium
  projection, 100,987,682. PRK: WPP 2024 Medium projection, 26,498,835;
  its incidence remains unavailable. Sexes stay separate, both is audit-only.
- Small derived-minus-published WPP differences remain visible: VNM -4,
  PRK +1, inactive Australia WPP -8. No averaging or correction.
- Fine-age/display counts reconcile; frozen 2A scientific profile hashes and
  original-snapshot restoration are unchanged.

## Scientific design

Six candidate incidence/pressure-shape approaches are compared. A reviewed
forward disease/natural-history operator is required before claiming disease
calibration. Incidence smoothing alone is not that operator, and no Styblo
conversion is used.

Pre-2000 constant pressure is a provisional, explicitly selected reference,
not an automatic default; bounded trend, supplied history, anchor-informed
history and coherent historical uncertainty are required comparisons.

An isolated first-infection catalytic diagnostic provides cumulative hazard,
uninfected/recent-first/remote-only probabilities, synthetic one-scale recovery
and explicit-RNG categorical checks. It excludes reinfection, clearance,
disease exits, selective mortality and migration and is not a country TBI fit.

Six proposed calibration contract families use 1.0.0-design.1, status
design_only. Frozen 2A snapshots and executable 1.1.0 profiles are unchanged.
See calibration_specification.md, calibration_contracts.md and
decision_register_2b.md for assumptions, review gates and identifiers.

## Validation and protection

190 focused tests passed in 361.25 seconds, no skips; two read-only openpyxl
dependency deprecation warnings. Compilation, canonical hash checks, offline
loading, source/native reconstruction, source classification, restoration,
synthetic mathematics, protection and whitespace checks passed.
See milestone2b_validation.md and protected_fingerprints_2b.json.

Protected worktrees, dirty explicit-TBI Milestone 3H files, Starsim reports,
WHO sibling, all protected refs and 128 frozen tracked files were unchanged.
Only the authorized country-model branch may advance; the delivery response
records the commit, push and final post-push state.

## Open decisions and Milestone 2C

Review infection target and recency semantics, age-dependent natural history,
assay operators, mortality/migration and open-tail handling. Obtain reviewed
TBI observations and historical exposure/denominators where required.
Legacy upstream cache lineage and existing WHO terms clarification remain open.

Milestone 2C should validate one-scale catalytic calibration with synthetic
anchors, then a reviewed assay/survey pilot and explicit history scenarios.
Review the disease-forward operator before calling a profile incidence-calibrated.
Do not add risk-factor allocation or connect projection engines yet.
