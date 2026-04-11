# APY v9 — Simple mode vs Advanced mode specification

## Purpose
This file defines exactly what “Simple mode” and “Advanced mode” mean for the APY v9 MATLAB model and any future web app built from it.

Codex should not infer these modes from context. It should follow this file literally when proposing or implementing user-facing inputs.

## General rule
- **Simple mode** is for ordinary public health users.
- **Advanced mode** is for technical users, analysts, or model maintainers.
- If there is any doubt about where a field belongs, prefer **Advanced mode**.
- Simple mode should expose only the minimum information needed to run a valid APY scenario.

## Simple mode fields
These fields should be visible to most users.

### Scenario basics
- Scenario name
- Use APY default scenario (yes/no)
- Reset to APY default
- Load scenario
- Save scenario

### Epidemiology inputs
- LTBI prevalence target
- Active TB prevalence target
- Age distribution file or uploaded table
- Risk-factor prevalence inputs:
  - close contact
  - marijuana use
  - renal disease
  - diabetes
  - smoking
  - chronic lung disease
  - alcohol / drugs
  - female prevalence
  - BCG prevalence

Rules:
- For most risk factors, users may enter either:
  - one overall prevalence, or
  - a 3-element vector for 0–4, 5–14, and >=15
- Female and BCG remain overall scalar inputs only.

### Program pathway
- Test choice
  - IGRA
  - TST
- Treatment choice
  - 3HP
  - 4R
  - 3HR
  - 6H
  - 9H

### Cascade of care
- Treatment start probability
- Completion probability
- ADR stop probability
- Full-course efficacy

### Strategy settings
- Screening strategy
  - random
  - ltbi
  - cure
  - prevent
- Screening coverage

### Outputs shown first in Simple mode
- Number screened
- Number starting treatment
- Number completing treatment
- Number of false positives treated
- Number of infections protected / cured
- Number of active TB cases prevented
- NNS to protect one infection
- NNS to prevent one active TB case
- NNT to protect one infection
- NNT to prevent one active TB case
- Do-nothing comparator summary
- Attributable-risk summary

## Advanced mode fields
These fields should be hidden unless the user explicitly expands Advanced mode.

### Calibration controls
- Age OR target for LTBI calibration
- Early/late progression ratio
- Follow-up horizon
- Screen window

### Test-performance overrides
- IGRA sensitivity
- IGRA specificity
- TST sensitivity
- TST specificity in BCG-vaccinated
- TST specificity in non-BCG users

### Disease / progression overrides
- Disease OR overrides:
  - contact
  - marijuana
  - renal
  - diabetes
  - smoking
  - chronic lung disease
  - alcohol / drugs

### Technical / reproducibility controls
- Cohort size N
- Number of replicates
- Random seed
- Age 85+ upper bound
- Age distribution sheet name / number
- Age distribution table override
- Raw file paths

### Debug / analyst outputs
- Example cohort table
- Full raw results table
- Full attributable-risk table
- Targeting profile tables
- Chart export paths

## UI rules
- The default view should open in **Simple mode**.
- Advanced mode should be collapsed by default.
- Advanced controls should not be required to run a valid APY scenario.
- Every field in Simple mode should have a tooltip or help text.
- Any dangerous or easily misunderstood field should live in Advanced mode.

## Validation rules
- A valid APY run must be possible from Simple mode alone.
- If a user provides invalid values, show a validation message before the run starts.
- The app should clearly label whether the run used:
  - APY defaults, or
  - custom inputs.

## Policy rule
This model is for planning and prioritisation only.
It must not be framed or implemented as a tool for refusing screening or treatment.

## Codex instruction
When asked to implement or modify the APY v9 user interface:
- treat this file as the source of truth for field placement
- do not move fields between Simple and Advanced without explicit instruction
- do not expose advanced technical settings in the default user view
