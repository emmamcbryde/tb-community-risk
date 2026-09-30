"""Parameter contract for the Starsim TB feasibility prototype (Milestone 1).

Every quantity carries its unit in its name:

* ``*_hazard_per_year``   continuous-time hazard, events per person-year;
* ``*_years``             a duration or a calendar time, in years;
* ``*_proportion``        a unitless fraction of the living population at the start;
* ``network_mean_degree`` expected number of persistent contacts per agent (unitless).

Per-timestep probabilities are never stored; they are derived with
``engine.starsim_tb.analytics.hazard_to_probability`` (``p = 1 - exp(-hazard * dt)``).

DEMONSTRATION VALUES ONLY. The defaults below exist to verify software behaviour.
They have not been selected as scientific estimates, are not calibrated, and must
not be used for clinical, programmatic or policy decisions.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import hashlib
import json
import math
from typing import Any, Mapping

SCHEMA_VERSION = "starsim_tb_m1_parameters/1"

BIOLOGICAL_STATES = ("susceptible", "recent_infection", "remote_infection", "active_pulmonary_tb")

DEMONSTRATION_VALUES_NOTICE = (
    "Demonstration values for software verification only; not evidence-based, "
    "not calibrated, and not selected as final scientific estimates."
)

# Illustrative starting mix; must sum to one. Not derived from any population.
DEMONSTRATION_INITIAL_STATE_PROPORTIONS: dict[str, float] = {
    "susceptible": 0.75,
    "recent_infection": 0.02,
    "remote_infection": 0.225,
    "active_pulmonary_tb": 0.005,
}

PROPORTION_SUM_TOLERANCE = 1e-9


@dataclass(frozen=True)
class PrototypeParameters:
    """Complete, explicit input set for one prototype simulation."""

    # Population (agents). No births, deaths or migration in Milestone 1.
    n_agents: int = 10_000
    initial_state_proportions: Mapping[str, float] = field(
        default_factory=lambda: dict(DEMONSTRATION_INITIAL_STATE_PROPORTIONS)
    )

    # Time (calendar years). A monthly step is 1/12 year.
    start_year: float = 2000.0
    stop_year: float = 2020.0
    timestep_years: float = 1.0 / 12.0

    # Contact network: persistent Erdos-Renyi graph (ss.StaticNet), fixed for the run.
    network_mean_degree: float = 8.0

    # Transmission: hazard of infection per susceptible-infectious contact pair,
    # per year, for a source with relative infectiousness 1.
    transmission_hazard_per_contact_per_year: float = 0.1
    # Multiplier on the transmission hazard for active pulmonary TB. Recent and
    # remote infection are structurally non-infectious (not a parameter).
    relative_infectiousness_active_pulmonary_tb: float = 1.0

    # Natural history (continuous-time hazards per person-year).
    recent_progression_hazard_per_year: float = 0.01
    # 0.19 per year gives a mean residence in recent infection of
    # 1 / (0.01 + 0.19) = 5 years under the demonstration progression hazard.
    recent_to_remote_hazard_per_year: float = 0.19
    remote_progression_hazard_per_year: float = 0.001

    rand_seed: int = 20260930

    def __post_init__(self) -> None:
        object.__setattr__(self, "initial_state_proportions", dict(self.initial_state_proportions))
        validate_parameters(self)

    # Derived quantities -------------------------------------------------------------
    @property
    def recent_total_exit_hazard_per_year(self) -> float:
        return self.recent_progression_hazard_per_year + self.recent_to_remote_hazard_per_year

    @property
    def recent_mean_residence_years(self) -> float:
        """Mean time spent in recent infection, 1 / (h_f + gamma); infinite if both are zero."""
        total = self.recent_total_exit_hazard_per_year
        return math.inf if total == 0 else 1.0 / total

    @property
    def duration_years(self) -> float:
        return self.stop_year - self.start_year

    @property
    def n_timesteps(self) -> int:
        """Number of intervals between start and stop."""
        return int(round(self.duration_years / self.timestep_years))

    def with_updates(self, **changes: Any) -> "PrototypeParameters":
        return replace(self, **changes)

    # Serialisation -----------------------------------------------------------------
    def to_canonical_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["initial_state_proportions"] = {k: float(self.initial_state_proportions[k]) for k in BIOLOGICAL_STATES}
        for key, value in list(data.items()):
            if isinstance(value, float):
                data[key] = float(value)
        return {"schema_version": SCHEMA_VERSION, "parameters": data}

    def canonical_json(self) -> str:
        """Deterministic serialisation: sorted keys, no whitespace, shortest float repr."""
        return json.dumps(self.to_canonical_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False)

    def parameter_hash(self) -> str:
        return hashlib.sha256(self.canonical_json().encode("utf-8")).hexdigest()


def from_recent_mean_residence(
    recent_mean_residence_years: float,
    recent_progression_hazard_per_year: float = PrototypeParameters.recent_progression_hazard_per_year,
    **other: Any,
) -> PrototypeParameters:
    """Build parameters from a mean residence time in recent infection.

    Mean residence is 1 / (h_f + gamma), so gamma = 1 / residence - h_f.
    """
    if not (recent_mean_residence_years > 0 and math.isfinite(recent_mean_residence_years)):
        raise ValueError("recent_mean_residence_years must be positive and finite")
    gamma = 1.0 / recent_mean_residence_years - recent_progression_hazard_per_year
    if gamma < 0:
        raise ValueError("Mean residence is shorter than progression alone allows (gamma would be negative)")
    return PrototypeParameters(
        recent_progression_hazard_per_year=recent_progression_hazard_per_year,
        recent_to_remote_hazard_per_year=gamma,
        **other,
    )


def validate_initial_state_proportions(proportions: Mapping[str, float]) -> None:
    keys = set(proportions)
    expected = set(BIOLOGICAL_STATES)
    if keys != expected:
        raise ValueError(f"initial_state_proportions must have exactly {sorted(expected)}; got {sorted(keys)}")
    for key, value in proportions.items():
        if not (isinstance(value, (int, float)) and math.isfinite(value) and 0.0 <= value <= 1.0):
            raise ValueError(f"initial proportion for {key} must be in [0, 1]; got {value!r}")
    total = math.fsum(proportions.values())
    if abs(total - 1.0) > PROPORTION_SUM_TOLERANCE:
        raise ValueError(f"initial_state_proportions must sum to 1; got {total!r}")


def validate_parameters(pars: PrototypeParameters) -> None:
    if not (isinstance(pars.n_agents, int) and pars.n_agents > 0):
        raise ValueError("n_agents must be a positive integer")
    validate_initial_state_proportions(pars.initial_state_proportions)
    for name in ("start_year", "stop_year", "timestep_years"):
        if not math.isfinite(getattr(pars, name)):
            raise ValueError(f"{name} must be finite")
    if pars.stop_year <= pars.start_year:
        raise ValueError("stop_year must be after start_year")
    if not 0 < pars.timestep_years <= pars.duration_years * (1 + 1e-9):
        raise ValueError("timestep_years must be positive and no longer than the simulation")
    steps = pars.duration_years / pars.timestep_years
    if abs(steps - round(steps)) > 1e-6:
        raise ValueError("stop_year - start_year must be a whole number of timesteps")
    if not (math.isfinite(pars.network_mean_degree) and 0 <= pars.network_mean_degree < pars.n_agents):
        raise ValueError("network_mean_degree must be in [0, n_agents)")
    for name in (
        "transmission_hazard_per_contact_per_year",
        "relative_infectiousness_active_pulmonary_tb",
        "recent_progression_hazard_per_year",
        "recent_to_remote_hazard_per_year",
        "remote_progression_hazard_per_year",
    ):
        value = getattr(pars, name)
        if not (isinstance(value, (int, float)) and math.isfinite(value) and value >= 0):
            raise ValueError(f"{name} must be a finite non-negative number; got {value!r}")
    if not isinstance(pars.rand_seed, int):
        raise ValueError("rand_seed must be an integer")


def demonstration_parameters(**changes: Any) -> PrototypeParameters:
    """The labelled demonstration parameter set, optionally with overrides."""
    return PrototypeParameters(**changes)
