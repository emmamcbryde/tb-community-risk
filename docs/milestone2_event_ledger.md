# Milestone 2 APY Event Ledger

Milestone 2 adds the versioned `ltbi_screening_event_ledger_v1` contract for
direct-effects APY LTBI screening analyses. The ledger is an epidemiological
event contract, not a health-economic outcome engine. DALYs, ICERs, net
monetary benefit, PSA and early stopping rules are not calculated here.

The current static and individual-based analyses estimate direct benefits,
harms and costs for screened and treated individuals. Transmission-mediated
benefits are not yet included.

## Contract Shape

Each event-ledger bundle contains:

- `metadata`: scenario id/version, population preset, model type, backend,
  screening window, follow-up horizon, model version and timing convention.
- `definitions`: stable event machine names, labels, units and definitions.
- `replicateTotals`: one long-format row per event, arm and replicate.
- `annualEvents`: one long-format row per event, arm, replicate and model year.
- `validation`: reusable coherence checks and warnings.
- `summaries`: mean/median/range summaries by arm and event.

`modelType` is `expected_value` for deterministic fractional expectations and
`agent_based` for simulated replicate counts. `valueType` is `expected` or
`simulated_count` so deterministic values cannot be mistaken for observed
integer counts.

## Arms And Pairing

The comparator arm is current practice / no additional systematic LTBI
screening. It still includes ordinary background TB care and is not absence of
all TB care. Under the current direct-effects comparator, systematic screening
and TPT events are zero unless such background screening is explicitly added to
configuration.

The intervention arm is targeted LTBI screening and preventive treatment.
Agent-based comparator and intervention rows are paired by `replicateId` and
`pairedReplicateId` and use the same simulated cohort. The identity
`comparator active_tb_cases = intervention active_tb_cases +
active_tb_cases_prevented` is validated cumulatively and annually.

## Annual Timing

Model year 0 represents `[0,1)`, model year 1 represents `[1,2)`, and so on.
The final interval may be shorter for a non-integer follow-up horizon.
Post-horizon programme events are reported as `modelYear = -1` with
`timeInterval = post_horizon`.

Agent-based timing uses simulated times: screening and testing at `screenTime`,
treatment start at `screenTime`, treatment stop/completion/effective treatment
at `treatmentStopTime`, active TB at untreated onset time, and prevented active
TB at the untreated onset time of the prevented case.

## Core Definitions

Core events cover population, screening, test outcomes, TPT eligibility and
starts, completion, ADR stops, other stops, partial courses, effective
treatment, active TB cases and active TB cases prevented. BCG-specific
false-positive events are retained additively where supported.

`true_positive_latent` is a test-positive person who is infected and latent at
screening. A test-positive person already active at screening is recorded
separately as `test_positive_active`. `active_tb_at_screen` is not assumed to
be screen-detected active TB. False-positive treatment can contribute treatment
use, harms and later costs, but never infection cure or active-TB prevention.
`infection_effectively_treated_total` is the event-ledger name corresponding to
the legacy `nCuredInfection` summary.

## Deterministic Expected-Value APY Model

The deterministic model is implemented under `engine/apy/expected_value.py`.
It is not the legacy `engine/static` mechanistic model and does not route
through `ui/static_ui.py`.

The expected-value model enumerates exact ages, existing binary infection-risk
factors, existing binary progression-risk factors and BCG status. It uses the
same independence assumptions as APY calibration and population generation.
Targeting scores use the same APY definitions: infection probability for LTBI
targeting, infection probability times survival to random screening for cure
targeting, and infection probability times preventable active-TB risk for
prevention targeting.

Random screening allocates coverage proportionally across all strata. Targeted
screening sorts strata by the selected APY score, then by infection risk,
progression multiplier and age as a deterministic tie rule. The boundary
stratum can be fractionally screened.

Screening is integrated deterministically over a uniform screening window. The
model uses the calibrated early and late progression hazards and the existing
regimen library for starts, completion, ADR stop, other stop, full-course
efficacy, partial-course efficacy, treatment duration and dose fractions.

## Validation Identities

Reusable validation checks enforce the specified screening, test, TPT and
effective-treatment identities, non-negative counts, treatment starts no
greater than eligibility, completions no greater than starts, effective
treatments no greater than true-positive starts, prevented TB no greater than
effective treatments, zero systematic screening/TPT in the comparator, annual
reconciliation, paired comparator/intervention TB identity, and no
false-positive epidemiological benefit.

Expected values use a floating-point tolerance. Simulated ledger rows are
constructed from mutually exclusive masks and should satisfy integer identities.

## Expected Value Versus ABM

Matched validation scenarios compare deterministic expectations with the
agent-based mean, not the median or a single replicate. Unit tests keep this
comparison computationally modest and use a broad Monte Carlo tolerance.
Heavier distributional validation should remain outside the ordinary unit-test
suite.

## MATLAB Compatibility

MATLAB remains the APY reference backend. Existing MATLAB cumulative raw fields
can support only partial total-ledger mapping unless the MATLAB simulation
emits the new event masks and event times. Annual MATLAB ledger output is not
fabricated from cumulative summaries; annual availability should remain false
for that backend until reference MATLAB output includes event timing.
