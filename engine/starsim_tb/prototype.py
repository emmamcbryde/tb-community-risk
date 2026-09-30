"""Starsim-native TB natural-history feasibility prototype (Milestone 1).

Four mutually exclusive biological states, stored as Starsim ``BoolState`` arrays:

    susceptible -> recent_infection            (transmission on the contact network)
    recent_infection -> active_pulmonary_tb    (hazard h_f)       } competing
    recent_infection -> remote_infection       (hazard gamma)     } exits
    remote_infection -> active_pulmonary_tb    (hazard h_r)

recent -> remote is ageing of infection, not disease progression. Active pulmonary
TB is absorbing in this milestone (no treatment, recovery or death).

Discrete-time convention: result index ``ti`` is the state at calendar time
``t_ti``; flows at ``ti`` are the events in the interval (t_{ti-1}, t_ti]. Index 0
holds the initial condition only, so no transitions or infections occur at ti = 0.

Not calibrated, not a clinical model, not a country model, not the SA Health
model and not the Schwalb/TBsim model. Not for policy use.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
import platform
import subprocess
import sys
import sysconfig
import time
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import starsim as ss

from engine.starsim_tb.analytics import hazard_to_probability
from engine.starsim_tb.parameters import BIOLOGICAL_STATES, DEMONSTRATION_VALUES_NOTICE, PrototypeParameters

MODULE_NAME = "tb_m1"
ANALYZER_NAME = "tb_m1_integrity"

# Codes stored in ``active_origin`` (which infected state an active case came from).
ORIGIN_NONE, ORIGIN_RECENT, ORIGIN_REMOTE, ORIGIN_INITIAL = 0, 1, 2, 3


def assign_competing_outcomes(
    u_event: np.ndarray, u_which: np.ndarray, hazards_per_year: Mapping[str, float], dt_years: float
) -> np.ndarray:
    """Assign at most one exit per agent from competing constant hazards.

    ``u_event`` decides whether any exit occurs (probability 1 - exp(-H dt));
    ``u_which`` picks which exit, with probability h_k / H. Exits are laid out on
    [0, 1) in sorted-name order, so the assignment for each agent is identical
    however the caller orders ``hazards_per_year``. Returns an array of exit names,
    with "" for no exit.
    """
    names = sorted(hazards_per_year)
    rates = np.array([hazards_per_year[k] for k in names], dtype=float)
    if np.any(rates < 0):
        raise ValueError("hazards must be non-negative")
    outcome = np.full(len(u_event), "", dtype=object)
    total = rates.sum()
    if total == 0 or len(u_event) == 0:
        return outcome
    event = np.asarray(u_event) < hazard_to_probability(total, dt_years)
    edges = np.cumsum(rates / total)
    edges[-1] = 1.0  # Guard against rounding leaving a gap below 1
    which = np.minimum(np.searchsorted(edges, np.asarray(u_which), side="right"), len(names) - 1)
    outcome[event] = np.array(names, dtype=object)[which[event]]
    return outcome


class TBNaturalHistoryM1(ss.Infection):
    """Recent/remote infection with competing exits and network transmission."""

    def __init__(self, prototype_pars: PrototypeParameters, **kwargs: Any) -> None:
        super().__init__(name=MODULE_NAME, label="TB natural history (M1 prototype)")
        self.prototype_pars = prototype_pars
        pp = prototype_pars
        self.define_pars(
            beta=ss.peryear(pp.transmission_hazard_per_contact_per_year),
            init_prev=None,  # Initial states are assigned explicitly in init_post()
            initial_state_proportions=dict(pp.initial_state_proportions),
            recent_progression=ss.peryear(pp.recent_progression_hazard_per_year),
            recent_to_remote=ss.peryear(pp.recent_to_remote_hazard_per_year),
            remote_progression=ss.peryear(pp.remote_progression_hazard_per_year),
            rel_trans_active=pp.relative_infectiousness_active_pulmonary_tb,
        )
        self.update_pars(**kwargs)

        # ss.Infection defines susceptible (kept) and an `infected` BoolState (replaced by
        # a derived alias). `infectious` points only at active pulmonary TB.
        self.define_states(
            ss.BoolState("recent_infection", label="Recent infection"),
            ss.BoolState("remote_infection", label="Remote infection"),
            ss.BoolState("active_pulmonary_tb", label="Active pulmonary TB"),
            ss.FloatArr("ti_remote", label="Timestep of recent -> remote"),
            ss.FloatArr("ti_active", label="Timestep of onset of active pulmonary TB"),
            ss.FloatArr("active_origin", default=ORIGIN_NONE, label="Origin of active TB"),
            reset=["infected"],
            infected=_any_infection,
        )
        self.define_aliases(infectious="active_pulmonary_tb", overwrite=True)

        # One CRN random stream per decision; each is drawn once per timestep.
        self.dist_initial_state = ss.choice(a=len(BIOLOGICAL_STATES), p=self._initial_p(), name="initial_state")
        self.dist_recent_event = ss.random(name="recent_event")
        self.dist_recent_which = ss.random(name="recent_which")
        self.dist_remote_event = ss.random(name="remote_event")

    def _initial_p(self) -> list[float]:
        props = self.prototype_pars.initial_state_proportions
        return [float(props[k]) for k in BIOLOGICAL_STATES]

    # Initialisation ---------------------------------------------------------------------
    def init_post(self) -> None:
        super().init_post()  # init_prev is None, so ss.Infection seeds nothing
        uids = self.sim.people.auids
        codes = np.asarray(self.dist_initial_state.rvs(uids))
        ti = self.ti
        self.susceptible[uids[codes != 0]] = False
        recent, remote, active = uids[codes == 1], uids[codes == 2], uids[codes == 3]
        # Constant hazards are memoryless, so the unknown age of prevalent infections is irrelevant.
        self.recent_infection[recent] = True
        self.ti_infected[recent] = ti
        self.remote_infection[remote] = True
        self.ti_infected[remote] = ti
        self.ti_remote[remote] = ti
        self.active_pulmonary_tb[active] = True
        self.ti_active[active] = ti
        self.active_origin[active] = ORIGIN_INITIAL
        self.rel_trans[active] = self.pars.rel_trans_active
        self.initial_counts = {k: int(np.count_nonzero(codes == i)) for i, k in enumerate(BIOLOGICAL_STATES)}

    def init_results(self) -> None:
        super().init_results()  # n_<state> stocks for each BoolState; new/cum_infections; prevalence
        self.define_results(
            ss.Result("new_recent_to_remote", dtype=int, scale=True, label="Flow: recent -> remote (per step)"),
            ss.Result("new_active_from_recent", dtype=int, scale=True, label="Flow: active TB from recent (per step)"),
            ss.Result("new_active_from_remote", dtype=int, scale=True, label="Flow: active TB from remote (per step)"),
            ss.Result("new_active", dtype=int, scale=True, label="Flow: all new active TB (per step)"),
            ss.Result("cum_recent_to_remote", dtype=int, scale=True, label="Cumulative flow: recent -> remote"),
            ss.Result("cum_active_from_recent", dtype=int, scale=True, label="Cumulative flow: active TB from recent"),
            ss.Result("cum_active_from_remote", dtype=int, scale=True, label="Cumulative flow: active TB from remote"),
            ss.Result("cum_active", dtype=int, scale=True, label="Cumulative flow: all new active TB"),
            ss.Result("prevalence_active_pulmonary_tb", dtype=float, scale=False, label="Stock proportion: active pulmonary TB"),
            ss.Result("person_years_recent", dtype=float, scale=True, label="Person-time (years) in recent infection during step"),
            ss.Result("person_years_remote", dtype=float, scale=True, label="Person-time (years) in remote infection during step"),
        )

    # Integration loop -------------------------------------------------------------------
    @property
    def dt_years(self) -> float:
        return float(self.t.dt_year)

    def step_state(self) -> None:
        """Natural-history transitions over the interval ending at this step."""
        recent = self.recent_infection.uids
        remote = self.remote_infection.uids
        # Draw every stream every step (also at ti = 0) so CRN positions never depend on state.
        u_event = self.dist_recent_event.rvs(recent)
        u_which = self.dist_recent_which.rvs(recent)
        u_remote = self.dist_remote_event.rvs(remote)
        if self.ti == 0:
            return
        res, ti, dt = self.results, self.ti, self.dt_years

        outcome = assign_competing_outcomes(
            u_event,
            u_which,
            {"active": self.pars.recent_progression.rate, "remote": self.pars.recent_to_remote.rate},
            dt,
        )
        to_active = recent[outcome == "active"]
        to_remote = recent[outcome == "remote"]
        # Remote progression uses only agents who were remote at the start of the step, so
        # nobody can make two transitions in one step.
        p_remote = hazard_to_probability(self.pars.remote_progression.rate, dt)
        remote_to_active = remote[np.asarray(u_remote) < p_remote]

        self.recent_infection[to_active] = False
        self._make_active(to_active, ORIGIN_RECENT)
        self.recent_infection[to_remote] = False
        self.remote_infection[to_remote] = True
        self.ti_remote[to_remote] = ti
        self.remote_infection[remote_to_active] = False
        self._make_active(remote_to_active, ORIGIN_REMOTE)

        res.new_recent_to_remote[ti] = len(to_remote)
        res.new_active_from_recent[ti] = len(to_active)
        res.new_active_from_remote[ti] = len(remote_to_active)
        res.new_active[ti] = len(to_active) + len(remote_to_active)

    def _make_active(self, uids: ss.uids, origin: int) -> None:
        self.active_pulmonary_tb[uids] = True
        self.ti_active[uids] = self.ti
        self.active_origin[uids] = origin
        self.rel_trans[uids] = self.pars.rel_trans_active

    def step(self):
        if self.ti == 0:
            return ss.uids(), ss.uids(), np.empty(0, dtype=ss.dtypes.int)
        return super().step()

    def set_prognoses(self, uids: ss.uids, sources: Any = None) -> None:
        """New infections: susceptible -> recent. ss.Infection counts them in new_infections."""
        if np.any(~self.susceptible[uids]):
            raise RuntimeError("Transmission targeted a non-susceptible agent; reinfection is out of scope")
        super().set_prognoses(uids, sources)
        self.recent_infection[uids] = True
        self.ti_infected[uids] = self.ti

    def step_die(self, uids: ss.uids) -> None:  # No mortality in M1; kept for demographic safety
        for name in BIOLOGICAL_STATES:
            getattr(self, name)[uids] = False

    def update_results(self) -> None:
        super().update_results()
        res, ti = self.results, self.ti
        n_alive = self.sim.people.n_alive
        res.prevalence_active_pulmonary_tb[ti] = res.n_active_pulmonary_tb[ti] / n_alive
        if ti > 0:
            # Trapezoid over the step: exact for a straight-line change, O(dt^2) otherwise.
            half_dt = 0.5 * self.dt_years
            res.person_years_recent[ti] = half_dt * (res.n_recent_infection[ti - 1] + res.n_recent_infection[ti])
            res.person_years_remote[ti] = half_dt * (res.n_remote_infection[ti - 1] + res.n_remote_infection[ti])

    def finalize_results(self) -> None:
        res = self.results
        for flow in ("recent_to_remote", "active_from_recent", "active_from_remote", "active"):
            res[f"cum_{flow}"][:] = np.cumsum(res[f"new_{flow}"][:])
        super().finalize_results()


def _any_infection(module: TBNaturalHistoryM1):
    """Derived state: any TB infection (recent, remote or active)."""
    return module.recent_infection | module.remote_infection | module.active_pulmonary_tb


class StateIntegrityAnalyzer(ss.Analyzer):
    """Checks every timestep that each living agent is in exactly one biological state."""

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(name=ANALYZER_NAME, **kwargs)

    def init_results(self) -> None:
        super().init_results()
        self.define_results(
            ss.Result("n_alive", dtype=int, scale=False, label="Stock: living agents"),
            ss.Result("n_in_exactly_one_state", dtype=int, scale=False, label="Stock: agents in exactly one state"),
            ss.Result("n_state_violations", dtype=int, scale=False, label="Stock: agents in zero or several states"),
            ss.Result("population_balance", dtype=int, scale=False, label="Sum of state stocks minus n_alive"),
        )

    def step(self) -> None:
        disease = self.sim.diseases[MODULE_NAME]
        alive = self.sim.people.alive
        membership = sum(getattr(disease, name).raw.astype(np.int64) for name in BIOLOGICAL_STATES)
        membership = membership[alive.raw]
        ti, res = self.ti, self.results
        res.n_alive[ti] = int(alive.sum())
        res.n_in_exactly_one_state[ti] = int(np.count_nonzero(membership == 1))
        res.n_state_violations[ti] = int(np.count_nonzero(membership != 1))
        res.population_balance[ti] = int(sum(disease.results[f"n_{name}"][ti] for name in BIOLOGICAL_STATES)) - res.n_alive[ti]


# Sim construction and running ---------------------------------------------------------

def build_sim(pars: PrototypeParameters, analyzers: list | None = None) -> ss.Sim:
    """Build (but do not run) the Starsim simulation for one parameter set."""
    disease = TBNaturalHistoryM1(pars)
    networks = []
    if pars.network_mean_degree > 0:
        networks.append(ss.StaticNet(n_contacts=pars.network_mean_degree, name="contacts"))
    extra = list(analyzers) if analyzers else []
    return ss.Sim(
        n_agents=pars.n_agents,
        start=pars.start_year,
        stop=pars.stop_year,
        dt=ss.years(pars.timestep_years),
        rand_seed=pars.rand_seed,
        diseases=disease,
        networks=networks,
        analyzers=[StateIntegrityAnalyzer(), *extra],
        verbose=0,
    )


@dataclass
class PrototypeRun:
    pars: PrototypeParameters
    sim: ss.Sim
    runtime_seconds: float
    metadata: dict[str, Any]

    @property
    def disease(self) -> TBNaturalHistoryM1:
        return self.sim.diseases[MODULE_NAME]

    @property
    def integrity(self) -> StateIntegrityAnalyzer:
        return self.sim.analyzers[ANALYZER_NAME]


def run_prototype(pars: PrototypeParameters, analyzers: list | None = None) -> PrototypeRun:
    sim = build_sim(pars, analyzers=analyzers)
    started = time.perf_counter()
    sim.run()
    runtime = time.perf_counter() - started
    n_points = len(sim.diseases[MODULE_NAME].results.n_susceptible)
    if n_points != pars.n_timesteps + 1:
        raise RuntimeError(f"Starsim built {n_points - 1} timesteps but the parameters imply {pars.n_timesteps}")
    return PrototypeRun(pars=pars, sim=sim, runtime_seconds=runtime, metadata=run_metadata(pars, runtime))


def _git_commit() -> str | None:
    root = Path(__file__).resolve().parents[2]
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"], check=True, capture_output=True, text=True, timeout=10
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None


def _git_dirty() -> bool | None:
    root = Path(__file__).resolve().parents[2]
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=no"],
            check=True, capture_output=True, text=True, timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return bool(out.stdout.strip())


def run_metadata(pars: PrototypeParameters, runtime_seconds: float | None = None) -> dict[str, Any]:
    return {
        "model": "Starsim TB natural-history feasibility prototype, Milestone 1",
        "status": "Feasibility prototype. Not calibrated; not a clinical, country, SA Health or Schwalb model; not for policy use.",
        "parameter_notice": DEMONSTRATION_VALUES_NOTICE,
        "python_version": platform.python_version(),
        "starsim_version": ss.__version__,
        "platform": platform.platform(),
        "python_build_platform": sysconfig.get_platform(),  # e.g. win-amd64 (may be emulated on ARM64)
        "host_machine": platform.machine(),
        "git_commit": _git_commit(),
        "git_tracked_changes_present": _git_dirty(),
        "parameter_hash_sha256": pars.parameter_hash(),
        "rand_seed": pars.rand_seed,
        "start_year": pars.start_year,
        "stop_year": pars.stop_year,
        "timestep_years": pars.timestep_years,
        "n_timesteps": pars.n_timesteps,
        "network": {
            "type": "ss.StaticNet (Erdos-Renyi via networkx.fast_gnp_random_graph), fixed for the whole run",
            "mean_degree": pars.network_mean_degree,
            "edge_transmission_hazard_per_year": pars.transmission_hazard_per_contact_per_year,
        } if pars.network_mean_degree > 0 else {"type": "none"},
        "parameters": pars.to_canonical_dict(),
        "runtime_seconds": runtime_seconds,
        "process_id": os.getpid(),
        "argv0": Path(sys.argv[0]).name if sys.argv and sys.argv[0] else None,
    }


# Result tables --------------------------------------------------------------------------

STOCK_KEYS = {f"stock_n_{name}": f"n_{name}" for name in BIOLOGICAL_STATES}
FLOW_KEYS = {
    "flow_new_infections": "new_infections",
    "flow_new_recent_to_remote": "new_recent_to_remote",
    "flow_new_active_from_recent": "new_active_from_recent",
    "flow_new_active_from_remote": "new_active_from_remote",
    "flow_new_active": "new_active",
    "cumflow_infections": "cum_infections",
    "cumflow_recent_to_remote": "cum_recent_to_remote",
    "cumflow_active_from_recent": "cum_active_from_recent",
    "cumflow_active_from_remote": "cum_active_from_remote",
    "cumflow_active": "cum_active",
}
PROPORTION_KEYS = {"stockprop_active_pulmonary_tb": "prevalence_active_pulmonary_tb"}
PERSON_TIME_KEYS = {
    "persontime_years_recent": "person_years_recent",
    "persontime_years_remote": "person_years_remote",
}


def results_table(run: PrototypeRun) -> dict[str, list]:
    """Column-oriented results with labels that separate stocks, flows, proportions and person-time."""
    res = run.disease.results
    years = run.pars.start_year + np.arange(len(res.n_susceptible)) * run.pars.timestep_years
    table: dict[str, list] = {"time_index": list(range(len(years))), "time_years": [float(y) for y in years]}
    for mapping, cast in ((STOCK_KEYS, int), (FLOW_KEYS, int), (PROPORTION_KEYS, float), (PERSON_TIME_KEYS, float)):
        for out_key, res_key in mapping.items():
            table[out_key] = [cast(v) for v in np.asarray(res[res_key][:])]
    return table


def summary(run: PrototypeRun) -> dict[str, Any]:
    res = run.disease.results
    pp = run.pars
    last = -1
    return {
        "stock_final": {name: int(res[f"n_{name}"][last]) for name in BIOLOGICAL_STATES},
        "stock_initial": {name: int(res[f"n_{name}"][0]) for name in BIOLOGICAL_STATES},
        "cumflow_total": {
            "infections": int(res.cum_infections[last]),
            "recent_to_remote": int(res.cum_recent_to_remote[last]),
            "active_from_recent": int(res.cum_active_from_recent[last]),
            "active_from_remote": int(res.cum_active_from_remote[last]),
            "active": int(res.cum_active[last]),
        },
        "persontime_years_total": {
            "recent_infection": float(np.sum(res.person_years_recent[:])),
            "remote_infection": float(np.sum(res.person_years_remote[:])),
        },
        "rate_per_100k_person_years": {
            "active_from_recent": _rate(res.cum_active_from_recent[last], np.sum(res.person_years_recent[:])),
            "active_from_remote": _rate(res.cum_active_from_remote[last], np.sum(res.person_years_remote[:])),
        },
        "hazard_input_per_100k_person_years": {
            "recent_progression": 1e5 * pp.recent_progression_hazard_per_year,
            "remote_progression": 1e5 * pp.remote_progression_hazard_per_year,
        },
        "stockprop_active_pulmonary_tb_final": float(res.prevalence_active_pulmonary_tb[last]),
        "integrity_max_state_violations": int(np.max(run.integrity.results.n_state_violations[:])),
        "integrity_max_abs_population_balance": int(np.max(np.abs(run.integrity.results.population_balance[:]))),
    }


def _rate(events: float, person_years: float) -> float | None:
    return None if person_years <= 0 else float(1e5 * events / person_years)
