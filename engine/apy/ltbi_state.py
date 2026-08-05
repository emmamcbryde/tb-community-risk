from __future__ import annotations

from typing import Any

import numpy as np


RECENT_LTBI_DEFINITION = "infection acquired within the last 5 years"
RECENT_TO_REMOTE_RATE_PER_YEAR = 1.0 / 5.0
BASELINE_RECENT_FRACTION_UNRESOLVED_WARNING = (
    "No validated APY-specific baseline recent-LTBI fraction was found. "
    "The configured value is a compatibility placeholder and must be reviewed "
    "before the APY exemplar is treated as scientifically final."
)


def resolve_ltbi_state_assumptions(config: dict[str, Any]) -> dict[str, Any]:
    raw = config.get("ltbiStateAssumptions")
    assumptions = raw if isinstance(raw, dict) else {}
    value = assumptions.get(
        "baselineRecentLTBIProportion",
        config.get("baselineRecentLTBIProportion"),
    )
    if value is None or value == []:
        return {
            "baselineRecentLTBIProportion": None,
            "recentToRemoteTransitionRatePerYear": RECENT_TO_REMOTE_RATE_PER_YEAR,
            "recentDefinitionYears": 5.0,
            "source": "Older static/transmission-dynamic architecture",
            "status": "unresolved",
            "warning": BASELINE_RECENT_FRACTION_UNRESOLVED_WARNING,
        }
    p_recent = float(value)
    if p_recent < 0 or p_recent > 1:
        raise ValueError("baselineRecentLTBIProportion must be in [0,1].")
    rate = assumptions.get(
        "recentToRemoteTransitionRatePerYear",
        config.get("recentToRemoteTransitionRatePerYear", RECENT_TO_REMOTE_RATE_PER_YEAR),
    )
    rate = float(rate)
    if rate <= 0:
        raise ValueError("recentToRemoteTransitionRatePerYear must be > 0.")
    status = str(assumptions.get("status", "configured"))
    return {
        "baselineRecentLTBIProportion": p_recent,
        "recentToRemoteTransitionRatePerYear": rate,
        "recentDefinitionYears": float(assumptions.get("recentDefinitionYears", 5.0)),
        "source": assumptions.get("source", "Configured APY scenario field"),
        "status": status,
        "warning": (
            BASELINE_RECENT_FRACTION_UNRESOLVED_WARNING
            if status.startswith("unresolved")
            else None
        ),
    }


def recent_state_survival(time, mult_disease, lambda_early: float, transition_rate: float):
    time_array = np.maximum(np.asarray(time, dtype=float), 0.0)
    rate = float(lambda_early) * np.asarray(mult_disease, dtype=float)
    return np.exp(-(rate + float(transition_rate)) * time_array)


def recent_no_active_survival(
    time,
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    transition_rate: float,
):
    time_array = np.maximum(np.asarray(time, dtype=float), 0.0)
    early_rate = float(lambda_early) * np.asarray(mult_disease, dtype=float)
    late_rate = float(lambda_late) * np.asarray(mult_disease, dtype=float)
    trans = float(transition_rate)
    still_recent = np.exp(-(early_rate + trans) * time_array)
    denom = late_rate - early_rate - trans
    with np.errstate(divide="ignore", invalid="ignore"):
        moved_remote = (
            trans
            * np.exp(-late_rate * time_array)
            * (np.exp(denom * time_array) - 1.0)
            / denom
        )
    equal_case = trans * time_array * np.exp(-late_rate * time_array)
    moved_remote = np.where(np.isclose(denom, 0.0), equal_case, moved_remote)
    return np.clip(still_recent + moved_remote, 0.0, 1.0)


def remote_survival(time, mult_disease, lambda_late: float):
    return np.exp(-float(lambda_late) * np.asarray(mult_disease, dtype=float) * np.maximum(float(time), 0.0))


def recent_to_remote_latent_probability(
    time,
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    transition_rate: float,
):
    return np.maximum(
        recent_no_active_survival(
            time, mult_disease, lambda_early, lambda_late, transition_rate
        )
        - recent_state_survival(time, mult_disease, lambda_early, transition_rate),
        0.0,
    )


def mixed_baseline_survival(
    time,
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    transition_rate: float,
    baseline_recent_proportion: float,
):
    p_recent = float(baseline_recent_proportion)
    return (
        p_recent
        * recent_no_active_survival(
            time, mult_disease, lambda_early, lambda_late, transition_rate
        )
        + (1.0 - p_recent) * remote_survival(time, mult_disease, lambda_late)
    )


def mixed_baseline_event_between(
    start: float,
    end: float,
    mult_disease,
    lambda_early: float,
    lambda_late: float,
    transition_rate: float,
    baseline_recent_proportion: float,
):
    return np.maximum(
        mixed_baseline_survival(
            start,
            mult_disease,
            lambda_early,
            lambda_late,
            transition_rate,
            baseline_recent_proportion,
        )
        - mixed_baseline_survival(
            end,
            mult_disease,
            lambda_early,
            lambda_late,
            transition_rate,
            baseline_recent_proportion,
        ),
        0.0,
    )
