"""Independent analytic helpers for the Starsim TB feasibility prototype.

Nothing in this module imports Starsim or the prototype disease module, so the
simulation can be checked against it. All hazards are per person-year and all
times are in years.

Continuous-time model for one infected person (no transmission, no death):

    recent --h_f--> active
    recent --gamma--> remote        (ageing of infection, not disease progression)
    remote --h_r--> active
"""

from __future__ import annotations

import math
from typing import Mapping, Sequence

import numpy as np


# Hazard and probability conversion -----------------------------------------------------

def hazard_to_probability(hazard_per_year: float | np.ndarray, dt_years: float) -> float | np.ndarray:
    """Probability of at least one event in ``dt_years`` under a constant hazard: 1 - exp(-h dt)."""
    hazard = np.asarray(hazard_per_year, dtype=float)
    if np.any(hazard < 0) or dt_years < 0:
        raise ValueError("hazard and dt must be non-negative")
    out = -np.expm1(-hazard * dt_years)
    return float(out) if out.ndim == 0 else out


def probability_to_hazard(probability: float, dt_years: float) -> float:
    """Inverse of ``hazard_to_probability``: h = -ln(1 - p) / dt."""
    if not 0 <= probability < 1 or dt_years <= 0:
        raise ValueError("probability must be in [0, 1) and dt positive")
    return -math.log1p(-probability) / dt_years


def linear_probability_approximation(hazard_per_year: float, dt_years: float) -> float:
    """DIAGNOSTIC ONLY: the first-order approximation p ~ h dt. Not used by the model."""
    return hazard_per_year * dt_years


def competing_step_probabilities(hazards_per_year: Mapping[str, float], dt_years: float) -> dict[str, float]:
    """Exact per-step outcome probabilities for competing constant hazards.

    For exits k with hazards h_k and total H = sum h_k, the probability that the
    first exit within dt is k equals (h_k / H) * (1 - exp(-H dt)). The result
    depends only on the set of hazards, never on the order they are listed.
    """
    if any(h < 0 for h in hazards_per_year.values()):
        raise ValueError("hazards must be non-negative")
    total = math.fsum(hazards_per_year.values())
    if total == 0:
        return {k: 0.0 for k in hazards_per_year}
    p_any = -math.expm1(-total * dt_years)
    return {k: (h / total) * p_any for k, h in hazards_per_year.items()}


# Cohort formulas -------------------------------------------------------------------------

def remote_cumulative_active(h_r: float, t_years: float | np.ndarray) -> float | np.ndarray:
    """P(active by t | remote at 0) = 1 - exp(-h_r t)."""
    t = np.asarray(t_years, dtype=float)
    out = -np.expm1(-h_r * t)
    return float(out) if out.ndim == 0 else out


def recent_still_recent(h_f: float, gamma: float, t_years: float | np.ndarray) -> float | np.ndarray:
    """P(still recent at t | recent at 0) = exp(-(h_f + gamma) t)."""
    t = np.asarray(t_years, dtype=float)
    out = np.exp(-(h_f + gamma) * t)
    return float(out) if out.ndim == 0 else out


def recent_active_while_recent(h_f: float, gamma: float, t_years: float | np.ndarray) -> float | np.ndarray:
    """P(active by t, having progressed directly from recent), ignoring later remote progression.

    h_f / (h_f + gamma) * [1 - exp(-(h_f + gamma) t)]; equals the full answer when h_r = 0.
    """
    t = np.asarray(t_years, dtype=float)
    total = h_f + gamma
    out = np.zeros_like(t) if total == 0 else (h_f / total) * -np.expm1(-total * t)
    return float(out) if out.ndim == 0 else out


def recent_ever_remote(h_f: float, gamma: float, t_years: float | np.ndarray) -> float | np.ndarray:
    """P(moved recent -> remote by t) = gamma / (h_f + gamma) * [1 - exp(-(h_f + gamma) t)]."""
    t = np.asarray(t_years, dtype=float)
    total = h_f + gamma
    out = np.zeros_like(t) if total == 0 else (gamma / total) * -np.expm1(-total * t)
    return float(out) if out.ndim == 0 else out


def recent_active_via_remote(h_f: float, gamma: float, h_r: float, t_years: float | np.ndarray) -> float | np.ndarray:
    """P(active by t after first ageing to remote), closed form.

    Integral_0^t gamma e^{-H s} (1 - e^{-h_r (t - s)}) ds with H = h_f + gamma:
      = (gamma / H)(1 - e^{-H t}) - gamma e^{-h_r t} (1 - e^{-(H - h_r) t}) / (H - h_r)
    with the limit gamma t e^{-h_r t} for the second term when H == h_r.
    """
    t = np.asarray(t_years, dtype=float)
    total = h_f + gamma
    if gamma == 0:
        out = np.zeros_like(t)
    else:
        first = (gamma / total) * -np.expm1(-total * t)
        diff = total - h_r
        if abs(diff) < 1e-12:
            second = gamma * t * np.exp(-h_r * t)
        else:
            second = gamma * np.exp(-h_r * t) * -np.expm1(-diff * t) / diff
        out = first - second
    return float(out) if out.ndim == 0 else out


def recent_cumulative_active(h_f: float, gamma: float, h_r: float, t_years: float | np.ndarray) -> float | np.ndarray:
    """Total P(active by t | recent at 0), allowing progression while recent and later while remote."""
    return recent_active_while_recent(h_f, gamma, t_years) + recent_active_via_remote(h_f, gamma, h_r, t_years)


def recent_cumulative_active_numerical(h_f: float, gamma: float, h_r: float, t_years: float, n: int = 20_000) -> float:
    """Composite Simpson integration of the same quantity, as an independent check on the closed form.

    P = integral_0^t e^{-H s} [h_f + gamma (1 - e^{-h_r (t - s)})] ds.
    """
    if n % 2:
        n += 1
    s = np.linspace(0.0, t_years, n + 1)
    total = h_f + gamma
    f = np.exp(-total * s) * (h_f + gamma * -np.expm1(-h_r * (t_years - s)))
    h = t_years / n
    return float(h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()))


def expected_person_years_recent(h_f: float, gamma: float, t_years: float) -> float:
    """Expected years spent in recent infection over [0, t], per person recent at 0."""
    total = h_f + gamma
    return t_years if total == 0 else -math.expm1(-total * t_years) / total


def expected_person_years_remote(h_f: float, gamma: float, h_r: float, t_years: float, n: int = 20_000) -> float:
    """Expected years spent in remote infection over [0, t], per person recent at 0 (numerical)."""
    s = np.linspace(0.0, t_years, n + 1)
    total = h_f + gamma
    if gamma == 0:
        return 0.0
    if abs(total - h_r) < 1e-12:
        p_remote = gamma * s * np.exp(-h_r * s)
    else:
        p_remote = gamma * (np.exp(-h_r * s) - np.exp(-total * s)) / (total - h_r)
    return float(np.trapezoid(p_remote, s))


# Exact discrete-time chain ----------------------------------------------------------------

INFECTED_STATE_ORDER = ("recent_infection", "remote_infection", "active_pulmonary_tb")


def discrete_transition_matrix(h_f: float, gamma: float, h_r: float, dt_years: float) -> np.ndarray:
    """One-step transition matrix of the prototype's discrete-time scheme (rows: from).

    Recent exits are competing (exact competing-risk split); remote progression is
    1 - exp(-h_r dt); active is absorbing in Milestone 1.
    """
    recent = competing_step_probabilities({"active": h_f, "remote": gamma}, dt_years)
    p_remote_active = hazard_to_probability(h_r, dt_years)
    return np.array(
        [
            [1 - recent["active"] - recent["remote"], recent["remote"], recent["active"]],
            [0.0, 1 - p_remote_active, p_remote_active],
            [0.0, 0.0, 1.0],
        ]
    )


def discrete_chain_distribution(
    start: Sequence[float], h_f: float, gamma: float, h_r: float, dt_years: float, n_steps: int
) -> np.ndarray:
    """State distribution after ``n_steps`` of the discrete scheme, starting from ``start``."""
    matrix = discrete_transition_matrix(h_f, gamma, h_r, dt_years)
    return np.asarray(start, dtype=float) @ np.linalg.matrix_power(matrix, n_steps)


def binomial_tolerance(p: float, n: int, z: float = 4.0) -> float:
    """Half-width z * sqrt(p(1 - p) / n) for comparing a simulated proportion with p."""
    return z * math.sqrt(max(p * (1 - p), 0.0) / n)
