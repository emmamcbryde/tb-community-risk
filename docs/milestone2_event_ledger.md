# Milestone 2 APY Event Ledger

Milestone 2 adds the versioned `ltbi_screening_event_ledger_v3` contract for
direct-effects APY LTBI screening analyses. The ledger is an epidemiological
event contract, not a health-economic outcome engine. DALYs, ICERs, net
monetary benefit, PSA and early stopping rules are not calculated here.

The current static and individual-based analyses estimate direct benefits,
harms and costs for screened and treated individuals. Transmission-mediated
benefits are not yet included.

## Contract Shape

Each event-ledger bundle contains:

- `metadata`: scenario id/version, population preset, model type, backend,
  screening window, early progression period, active-TB calibration horizon,
  follow-up horizon, model version and timing convention.
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

## Programme And Natural-History Timing

The implementation keeps four concepts distinct:

- `screeningWindowYears`: duration over which screening is delivered. The APY
  demonstration default is 3 years.
- `earlyProgressionPeriodYears`: early-to-late progression-hazard breakpoint.
  The default is 2 years.
- `activeTBCalibrationHorizonYears`: horizon for the existing active-TB
  calibration target. The default is 2 years and retains the legacy
  `targetActive2y` numerical target.
- `followUpHorizonYears`: analytical active-TB follow-up horizon. The APY
  demonstration default is 20 years.
- `recentToRemoteTransitionRatePerYear`: Markov transition rate from recent to
  remote LTBI in the explicit LTBI-state pathway. It, not
  `earlyProgressionPeriodYears`, controls recent-state residence time.

Legacy aliases `screenWindow` and `followHorizon` remain compatibility mirrors,
but internal APY calculations use the explicit fields. `nActiveBy2y` remains
active TB by the two-year calibration horizon and does not change meaning when
the screening delivery window changes.

## Annual Timing

Model year 0 represents `[0,1)`, model year 1 represents `[1,2)`, and so on.
The final interval may be shorter for a non-integer follow-up horizon.
Programme events after follow-up retain the actual non-negative
`modelYear = floor(event time)` and set `withinFollowUp = false`. Active-TB
outcomes remain limited to the configured follow-up horizon.

Agent-based timing uses simulated times: screening and testing at `screenTime`,
treatment start at `screenTime`, treatment stop/completion/effective treatment
at `treatmentStopTime`, active TB at untreated onset time, and prevented active
TB at the untreated onset time of the prevented case.

## Core Definitions

Core events cover population, screening, test outcomes, TPT eligibility and
starts, completion, ADR stops, other stops, partial courses, effective
treatment, active TB cases and active TB cases prevented. BCG-specific
false-positive events are retained additively where supported.

Version 3 adds explicit prevalent-infection states. `infected_at_baseline`
equals `recent_ltbi_at_baseline + remote_ltbi_at_baseline`. Recent LTBI follows
the older static/transmission-dynamic fast/slow latent-state architecture. The
current Python APY implementation uses a continuous Markov `recent -> remote`
transition rate of `0.2` per year, giving a mean residence time of five years
in the recent/fast latent state. This is not a deterministic infection-time
model in which every infection acquired within the last five years becomes
remote exactly at year five. Remote LTBI receives the late progression hazard
from model time zero. Recent LTBI can progress, be effectively treated, or
transition to remote LTBI; remote LTBI can progress or be effectively treated.
The APY direct-effects model is a closed cohort: there is no post-baseline
inflow from uninfected people to recent LTBI.

No validated APY-specific baseline recent-LTBI fraction was found in the APY
MATLAB v9 reference engine. The authoritative
`ltbiStateAssumptions.baselineRecentLTBIProportion` field is therefore `null`
unless a scenario supplies a value or an explicitly documented derivation
method. Development compatibility mode is off by default. When explicitly
selected it can supply a numerical placeholder only to keep old development
workflows runnable; those results are marked provisional and emit a validation
warning. Scenario authors should provide and document an explicit value before
treating the APY exemplar as clinician-ready or scientifically final.

At screening, infected people are classified as `active_tb_at_screen`,
`recent_latent_at_screen`, or `remote_latent_at_screen`. Aggregate latent,
true-positive, TPT start/completion, effective-treatment and prevented-TB
events reconcile exactly with their recent-plus-remote components.

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
It represents a weighted-stratum large-population expectation; targeted
allocation is therefore not necessarily the exact finite-N order statistic for
very small cohorts.
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
reconciliation, paired comparator/intervention TB identity cumulatively and
within every model year, and no false-positive epidemiological benefit.

## Eligibility

The ledger records resolved `eligible_population`. The current APY Python
workflow supports all-population eligibility. If a scenario configures
restricted eligibility but does not provide an implemented eligibility-selection
rule, validation fails with a clear error rather than silently applying
coverage to the full population.

Expected values use a floating-point tolerance. Simulated ledger rows are
constructed from mutually exclusive masks and should satisfy integer identities.

## Expected Value Versus ABM

Matched validation scenarios compare deterministic expectations with the
agent-based mean, not the median or a single replicate. Unit tests keep this
comparison computationally modest and use a broad Monte Carlo tolerance.
Heavier distributional validation should remain outside the ordinary unit-test
suite.

## MATLAB Compatibility

MATLAB remains the APY reference backend. Current Python hardening separates
screening delivery, the early progression period, the active-TB calibration
horizon and follow-up. MATLAB compatibility requires inspection/update before
claiming parity if MATLAB still couples the screening window to the early
progression or active-TB calibration horizon. Existing MATLAB cumulative raw
fields can support only partial total-ledger mapping unless the MATLAB
simulation emits the new recent/remote event masks and event times. The APY
MATLAB v9 engine samples a single `infected` state and does not expose baseline
recent/remote LTBI, recent-to-remote transition times, or screening-time
recent/remote latent states. Annual MATLAB ledger output is not fabricated from
cumulative summaries; annual availability remains false for that backend until
reference MATLAB output includes event timing.
