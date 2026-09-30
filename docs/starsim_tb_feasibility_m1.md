# Starsim TB feasibility prototype — Milestone 1

Status: **feasibility prototype, software-verified only.** Branch
`prototype/starsim-tb-spectrum`.

This prototype is **not** calibrated, **not** a clinical model, **not** a
country model, **not** the SA Health model, **not** the Schwalb/TBsim model and
**not** ready for policy use. Passing its tests shows that the code implements
the specification below. It does not show that the model or its parameter values
are scientifically valid.

## Purpose and exclusions

The purpose is to show that Starsim 3.6.1 can host a small, transparent TB
natural-history model with these features:

- four mutually exclusive biological states;
- infection transmitted over a Starsim contact network;
- elevated progression after recent infection and much slower progression from
  remote infection;
- ageing of infection from recent to remote as a competing exit, not a step
  towards disease.

The following are deliberately **not** included. They belong to later
milestones:

- minimal, incipient, subclinical, radiographic or extrapulmonary TB;
- symptoms, CXR, Xpert, screening, treatment or recovery;
- risk factors or HIV;
- births, deaths or migration;
- reinfection;
- health economics, WHO data, calibration or Streamlit controls.

TBsim is neither installed nor imported, and no part of the Schwalb natural
history is reproduced.

## Files

| File | Role |
| --- | --- |
| `engine/starsim_tb/parameters.py` | Parameter contract, validation, canonical JSON and SHA-256 hash. No Starsim import. |
| `engine/starsim_tb/analytics.py` | Independent analytic and discrete-chain helpers. No Starsim import. |
| `engine/starsim_tb/prototype.py` | Starsim disease module, integrity analyzer, sim builder, result tables and metadata. |
| `scripts/run_starsim_tb_prototype.py` | Runner: prints a summary and optionally writes JSON/CSV to `outputs/starsim_tb_m1/` (git-ignored). |
| `tests/test_starsim_tb_prototype.py` | Focused tests. The simulation tests are skipped when Starsim is absent. |
| `requirements-starsim-prototype.txt` | Separate prototype environment (`starsim==3.6.1`). `requirements.txt` is unchanged. |

The package sits beside the existing engines under `engine/` and imports none of
them. Nothing in the APY, SA Health, MATLAB, dynamic-model or Streamlit code
imports it.

## Architecture note (Starsim 3.6.1)

This design comes from reading the installed source in
`C:\Users\emmas\Envs\tb-starsim\Lib\site-packages\starsim` (`diseases.py`,
`modules.py`, `loop.py`, `distributions.py`, `networks.py`, `time.py`,
`analyzers.py`), not from online examples.

- **Base class.** `TBNaturalHistoryM1` subclasses `ss.Infection`, which
  supplies network transmission, the `susceptible` state, `rel_sus`,
  `rel_trans`, `ti_infected` and new/cumulative infection counting. The
  inherited `infected` `BoolState` is replaced by a derived alias (recent, remote
  or active). The `infectious` alias points only at `active_pulmonary_tb`.
- **Agent states.** The four biological states are `ss.BoolState` arrays owned
  by the module and registered with `ss.People` (`define_states`). Event times
  (`ti_infected`, `ti_remote`, `ti_active`) and `active_origin` are
  `ss.FloatArr`. Each `BoolState` automatically produces an `n_<state>` stock
  result.
- **Initialisation.** `init_post()` assigns every agent one state with a single
  `ss.choice` draw over the four initial proportions. `init_prev` is `None`, so
  `ss.Infection` seeds nothing of its own.
- **Transmission.** This is Starsim's own `ss.Infection.infect()`. For each edge
  of an `ss.StaticNet` (a persistent Erdős–Rényi graph from networkx
  `fast_gnp_random_graph` with the requested mean degree), the per-timestep
  probability is `rel_trans[source] × rel_sus[target] × (1 − exp(−β dt))`.
  `rel_trans` is masked by `infectious`, so only active TB transmits, and
  `rel_sus` is masked by `susceptible`. `β` is an `ss.peryear` rate, which
  Starsim converts with `Rate.to_prob(dt) = 1 − exp(−rate·dt)` (confirmed in
  `time.py` and numerically). New cases enter `recent_infection` in
  `set_prognoses()`.
- **Transitions.** `step_state()` runs once per step, before networks and
  transmission (loop order in `loop.py`). It applies the competing recent exits
  and remote progression described below, drawing from three `ss.random` streams
  (`recent_event`, `recent_which`, `remote_event`).
- **Results.** Stocks come from the automatic `n_<state>` results. Flows,
  person-time and active-TB prevalence are `ss.Result`s filled in
  `step_state()` and `update_results()`; cumulative flows are filled in
  `finalize_results()`. A small `ss.Analyzer` (`StateIntegrityAnalyzer`) checks
  state exclusivity and population balance at every step.
- **Seeds and random numbers.** `ss.Sim(rand_seed=...)` sets the base seed.
  Every `ss.Dist` gets its own stream, offset by a hash of its name, and advances
  by timestep (`jump_dt`). Draws are keyed to agent slots, which are Starsim's
  common random numbers (CRN), so the random value an agent receives does not
  depend on array order. `ss.StaticNet` builds its graph from its own seeded
  stream. Runs sharing a seed therefore share random numbers, even across
  different cohorts.
- **Multiple runs.** Replicates are separate `ss.Sim` runs with different
  `rand_seed`. `ss.MultiSim` exists but was not needed for the few replicates
  in this milestone.

**Starsim infrastructure:** people and state arrays, the integration loop, the
network, edge-level transmission and its hazard→probability conversion, random
streams and seeding, and result storage.
**TB-specific science (ours):** the four-state structure, the initial-state
assignment, the competing recent exits, remote progression, active TB as the
only infectious state, the discrete-time convention, and the result
definitions.

## States and transitions

| State | Meaning | Infectious |
| --- | --- | --- |
| `susceptible` | Never infected | no |
| `recent_infection` | Infected, elevated progression hazard | no |
| `remote_infection` | Infection has aged; low but non-zero progression hazard | no |
| `active_pulmonary_tb` | Infectious pulmonary disease; absorbing in M1 | yes (× `relative_infectiousness_active_pulmonary_tb`) |

```text
                    network transmission
  susceptible  ─────────────────────────────▶  recent_infection
                                                 │          │
                         h_f (progression)       │          │  γ (ageing of infection,
                                                 ▼          ▼   NOT progression)
                              active_pulmonary_tb ◀── h_r ── remote_infection
```

## Equations

All hazards are per person-year; `dt` is in years.

- Per-step probability from a hazard: `p = 1 − exp(−h·dt)`. The approximation
  `p ≈ h·dt` is used only in a diagnostic test that shows how it differs from
  the exact value.
- Remote cohort: `P(active by t) = 1 − exp(−h_r t)`.
- Recent cohort, where `H = h_f + γ`:
  - `P(still recent at t) = exp(−H t)`;
  - `P(progressed directly from recent by t) = (h_f/H)(1 − exp(−H t))`;
  - `P(moved to remote by t) = (γ/H)(1 − exp(−H t))`;
  - `P(active via remote by t) = (γ/H)(1 − e^{−Ht}) − γ e^{−h_r t}(1 − e^{−(H−h_r)t})/(H − h_r)`.
    The limit is used when `H = h_r`. A composite-Simpson integral of the same
    quantity cross-checks it (agreement 1.5 × 10⁻¹⁶ at 20 years).
- Mean residence in recent infection: `1/H`.

## Competing-risk implementation

At each step, each agent in `recent_infection` faces two competing exits with
constant hazards `h_f` and `γ`:

1. `u_event ~ U(0,1)`. If `u_event < 1 − exp(−H dt)`, one exit occurs.
2. `u_which ~ U(0,1)` picks the exit: progression if it falls in the first
   `h_f/H` of `[0,1)`, ageing otherwise. The exits are laid out in sorted-name
   order, so the result does not depend on the order in which hazards are
   listed.

This gives exactly `P(k) = (h_k/H)(1 − exp(−H dt))` per step, the
continuous-time competing-risk probability, with no sequential-Bernoulli bias.
Remote progression is applied only to agents who were remote at the start of
the step, so nobody makes two transitions in one step. The only
time-discretisation error is that someone who ages to remote in a step cannot
also progress from remote in that same step. At monthly steps this changes
`P(active by 20 y)` for a recent cohort by 3.8 × 10⁻⁵.

The tests show that permuting the hazard order leaves the per-step
probabilities unchanged (exactly) and the per-agent assignments unchanged
(exactly), and that running the full simulation with the reversed order gives an
identical event history. A diagnostic test shows that naive sequential
Bernoullis *are* order-dependent.

**Discrete-time convention.** Result index `ti` holds the state at time `t_ti`.
Flows at `ti` are the events in `(t_{ti−1}, t_ti]`. Index 0 is the initial
condition only, so no transitions or infections occur there. Within a step,
transitions are applied before transmission (Starsim's loop order), so a new
active case can transmit in the step it arises. Person-time per step uses the
trapezoid rule, `dt·(N_{ti−1} + N_ti)/2`.

## Demonstration parameters

**These are demonstration values for software verification only. They have not
been selected as final scientific estimates.** All of them are fields of
`PrototypeParameters`.

| Field | Value | Unit |
| --- | --- | --- |
| `n_agents` | 10,000 | agents |
| `initial_state_proportions` | S 0.75, recent 0.02, remote 0.225, active 0.005 | fraction (sums to 1, validated) |
| `start_year` / `stop_year` | 2000 / 2020 | calendar year |
| `timestep_years` | 1/12 (monthly) | years |
| `network_mean_degree` | 8 | persistent contacts per agent |
| `transmission_hazard_per_contact_per_year` | 0.1 | per infectious–susceptible pair per year |
| `relative_infectiousness_active_pulmonary_tb` | 1.0 | multiplier |
| `recent_progression_hazard_per_year` (`h_f`) | 0.01 | per person-year |
| `recent_to_remote_hazard_per_year` (`γ`) | 0.19 | per person-year |
| `remote_progression_hazard_per_year` (`h_r`) | 0.001 | per person-year |
| `rand_seed` | 20260930 | — |

The requested 5-year *mean residence time in recent infection* is treated as
`1/(h_f + γ)`, so `γ = 1/5 − 0.01 = 0.19` per year.
`from_recent_mean_residence()` builds parameters this way. If the intended
meaning was instead `1/γ = 5` years (γ = 0.2), that is a one-field change, and
it is listed below as a decision.

## Analytic validation

Closed cohorts, no transmission, monthly steps, 20 years. The tolerance is 4
binomial SE, plus the known discretisation error where one applies. The
tolerance was fixed before the runs.

| Check | Analytic | Simulated | Difference | Tolerance | n |
| --- | --- | --- | --- | --- | --- |
| Remote, `h_r` = 0.001: P(active by 20 y) | 0.01980 | 0.01775 | −0.00205 | ±0.00394 | 20,000 |
| Remote, `h_r` = 0.05: P(active by 20 y) | 0.63212 | 0.62560 | −0.00652 | ±0.01364 | 20,000 |
| Recent, `h_r` = 0: P(active by 20 y) | 0.04908 | 0.04905 | −0.00003 | ±0.00353 | 60,000 (3 seeds) |
| Recent, `h_r` = 0: P(still recent) | 0.01832 | 0.01802 | −0.00030 | ±0.00219 | 60,000 (3 seeds) |
| Recent, `h_r` = 0: P(remote) | 0.93260 | 0.93293 | +0.00033 | ±0.00409 | 60,000 (3 seeds) |
| Recent, full model: P(active by 20 y) | 0.06330 | 0.06107 | −0.00223 | ±0.00491 | 40,000 |
| Recent, 40 y: fraction progressing while recent | 0.04998 | 0.04770 | −0.00228 | ±0.00616 | 20,000 |

The single-seed rows share Starsim's common random numbers, so their small
same-sign deviations are correlated, not independent evidence of bias. Pooled
across 12 independent seeds (remote cohort, `h_r` = 0.05, 10 y), the simulated
mean is 0.39482 against an analytic 0.39347 (z = +1.36), and the per-seed
deviations have mixed signs. A test repeats this check with 8 seeds and also
checks that the replicate spread is binomial.

Person-time in recent infection for a 40,000-person recent cohort was
195,177 person-years, against an expected 196,337 (−0.6%). In the demonstration
run, incidence per person-year was 912 per 100k from recent infection and 86 per
100k from remote infection, against inputs of 1,000 and 100. These rest on 23
and 43 events, so they fall within Poisson error.

Under the demonstration values, about 5% (`h_f/H`) of recent infections ever
progress while recent, and about 6.3% of a recent cohort has active TB after 20
years. More than 80% of the cohort is in remote infection, so the model is not
a conveyor belt from infection to disease.

## Transmission checks

The tests confirm the following:

- **Only active TB transmits.** With no active TB, no progression and a very
  high `β`, there are zero infections.
- **Every infection comes from an active source.** Every `ss.infection_log`
  source is in active TB, and its onset precedes the infection.
- **Nobody is infected twice.** Each target appears once in the log, the
  susceptibility mask prevents overwriting, and initially infected agents keep
  their original infection time. Any attempt to infect a non-susceptible agent
  raises an error.
- **New infections enter `recent_infection`.**
- **Counts reconcile.** `Σ new_infections` equals the log length, equals the
  fall in susceptibles, and equals the rise in recent infection (with exits
  disabled). With exits enabled, stock changes reconcile with flows exactly at
  every step (`ΔS = −inf`, `ΔR = inf − r→m − R→A`, `ΔM = r→m − M→A`,
  `ΔA = R→A + M→A`).
- **Transmission is consistent across timesteps.** Across 4 seeds,
  demonstration runs averaged 533 infections at monthly steps and 524 at
  quarter-monthly steps. Contacts are persistent and `β` is a rate, so the
  per-pair hazard does not depend on `dt`. This was inspected, not asserted in
  a test.

## Outputs

`results_table()` columns are prefixed by kind: `stock_n_*` (counts at `t_ti`),
`flow_*` (events per step), `cumflow_*`, `stockprop_active_pulmonary_tb`
(prevalence proportion) and `persontime_years_*`. `summary()` adds final and
initial stocks, total flows, total person-time, empirical incidence per 100k
person-years next to the input hazards, and the integrity maxima. The runner's
optional JSON (`--json`) carries `format: starsim_tb_m1_prototype_results/1`,
full metadata, the summary and the time series. It is not connected to any
production result contract.

The metadata records Python version, Starsim version, platform (OS, Python
build platform and host machine), git commit and a
tracked-changes flag, the parameter SHA-256 (from sorted-key canonical JSON),
seed, start/stop/timestep, network configuration and runtime.

Demonstration run (seed 20260930):

| Result | Value |
| --- | --- |
| Initial S / recent / remote / active | 7,504 / 204 / 2,259 / 33 |
| Final S / recent / remote / active | 7,117 / 78 / 2,706 / 99 |
| Cumulative infections | 387 |
| Cumulative recent→remote | 490 |
| Cumulative active from recent / from remote | 23 / 43 |
| Person-years recent / remote | 2,522 / 50,053 |
| Final active prevalence | 0.99% |

Initial active TB (33, expected 50) is low for this seed by chance. Across 8
seeds the mean was 49.1.

## Runtime

On the reference machine (Windows 11 on ARM64, running x86-64 Python 3.12.13
`win-amd64` under emulation, Starsim 3.6.1),
10,000 agents over 20 years with monthly steps (240 steps) took **2.3–3.1 s** per
run after numba compilation. The first run in a fresh process took about 3.3 s.
Peak process memory, including the interpreter and all imported libraries, was
about **277 MB**. Quarter-monthly steps (960 steps) took about 8.6 s. The
focused test module runs in under a minute. No large replicate ensembles were
run.

## Known limitations

- Active TB is absorbing (no treatment, recovery or death), so prevalence only
  rises. The dynamics are not realistic.
- There is no demography, ageing, heterogeneity, reinfection or protection from
  prior infection.
- Hazards are constant and exponential, which makes the unknown age of prevalent
  infections at t0 irrelevant. Non-exponential sojourns would need explicit
  infection ages.
- The network is one static Erdős–Rényi graph. Transmission therefore depletes
  around each source, and results depend on network structure. There are no
  households or age mixing.
- `ss.StaticNet` sets edge probability to `mean_degree / n_agents`, so the
  realised mean degree is `8 × (n−1)/n`.
- The inherited `ss.Infection` `prevalence` result counts *any* infection. The
  active-TB proportion is `prevalence_active_pulmonary_tb`.
- Runs sharing a seed share random numbers across cohorts (CRN). This is useful
  for scenario contrasts, but single-seed validation comparisons are correlated.
- The main repository test environment (pytest, Streamlit) is separate from the
  Starsim environment. The new tests are plain `unittest` and skip their
  simulation part when Starsim is absent.

## Decisions required before the next milestone

1. Should the 5-year recent residence mean `1/(h_f + γ)` (implemented) or `1/γ`?
2. Which evidence sources set `h_f`, `γ` and `h_r`? Should `h_f` decline with
   time since infection rather than being constant?
3. Which states make up the disease spectrum, and in what order: minimal,
   subclinical and clinical, and whether regression or self-cure is allowed?
4. What exits does active TB have (treatment, self-cure, death), and do we need
   minimal demography for background mortality?
5. Should reinfection be allowed, with what partial protection, and what does it
   do to recent/remote status?
6. What contact structure should replace the random graph (household, community,
   age-structured), and how should `β` be parameterised?
7. Should we adopt `ss.MultiSim` and a CRN-aware scenario design for replicate
   ensembles?
