from __future__ import annotations

from copy import deepcopy
from typing import Any

from engine.apy.ltbi_state import (
    apply_ltbi_state_edit,
    enable_development_compatibility_mode,
    resolve_ltbi_state_assumptions,
)


RISK_FACTOR_LABELS = {
    "contact": "Close contact",
    "MJ": "Marijuana use",
    "renal": "Renal disease",
    "diabetes": "Diabetes",
    "smoking": "Smoking",
    "cld": "Chronic lung disease",
    "alcohol": "Harmful alcohol/other drug use",
    "female": "Female sex",
    "BCG": "BCG vaccination",
}

PRINCIPAL_RISK_FACTORS = ["contact", "MJ", "renal", "diabetes", "smoking", "cld", "alcohol"]
ADVANCED_RISK_FACTORS = ["female", "BCG"]
AGE_GROUP_LABELS = ["0-4", "5-14", ">=15"]


def percent_to_fraction(value: float) -> float:
    number = float(value)
    if number <= 0 or number >= 100:
        raise ValueError("Prevalence percentages must be greater than 0% and less than 100%.")
    return number / 100.0


def optional_percent_to_fraction(value: float | None) -> float | None:
    if value is None:
        return None
    return percent_to_fraction(value)


def fraction_to_percent(value: Any, default_percent: float | None = None) -> float | None:
    if value is None:
        return default_percent
    return float(value) * 100.0


def validate_fraction(value: float) -> float:
    number = float(value)
    if number < 0 or number > 1:
        raise ValueError("Prevalence fractions must be between 0 and 1.")
    return number


def risk_override_mode(value: Any) -> str:
    if value is None or value == []:
        return "Use default"
    if isinstance(value, (list, tuple)):
        return "Three age groups" if len(value) == 3 else "Use default"
    return "Single overall"


def risk_override_percent_values(value: Any) -> tuple[float, float, float]:
    if isinstance(value, (list, tuple)) and len(value) == 3:
        return tuple(float(v) * 100.0 for v in value)
    if value is not None and value != [] and not isinstance(value, (list, tuple, dict)):
        pct = float(value) * 100.0
        return (pct, pct, pct)
    return (0.0, 0.0, 0.0)


def risk_override_from_percentages(mode: str, values: list[float] | tuple[float, ...]) -> float | list[float] | None:
    if mode == "Use default":
        return None
    if mode == "Single overall":
        return percent_to_fraction(float(values[0]))
    if mode == "Three age groups":
        if len(values) != 3:
            raise ValueError("Three age-group prevalence values are required.")
        return [percent_to_fraction(float(value)) for value in values]
    raise ValueError(f"Unknown risk-factor override mode: {mode}")


def apply_epidemiology_updates(
    config: dict[str, Any],
    *,
    use_default_ltbi_prevalence: bool,
    ltbi_prevalence_percent: float | None,
    use_default_active_tb_prevalence: bool,
    active_tb_prevalence_percent: float | None,
    risk_prev_updates: dict[str, Any] | None = None,
) -> dict[str, Any]:
    updated = deepcopy(config)
    updated["ltbiPrevalence"] = (
        None
        if use_default_ltbi_prevalence
        else optional_percent_to_fraction(ltbi_prevalence_percent)
    )
    updated["activeTBPrevalence"] = (
        None
        if use_default_active_tb_prevalence
        else optional_percent_to_fraction(active_tb_prevalence_percent)
    )

    risk_prev = deepcopy(updated.get("riskPrev") or {})
    for key in PRINCIPAL_RISK_FACTORS + ADVANCED_RISK_FACTORS:
        risk_prev.setdefault(key, None)
    if risk_prev_updates:
        for key, value in risk_prev_updates.items():
            risk_prev[key] = value
    updated["riskPrev"] = risk_prev
    return updated


def apply_ltbi_state_assumption_update(
    config: dict[str, Any],
    *,
    baseline_recent_percent: float | None,
    transition_rate_per_year: float,
    source: str,
    status: str,
    notes: str = "",
) -> dict[str, Any]:
    baseline_fraction = (
        None
        if baseline_recent_percent is None
        else validate_fraction(float(baseline_recent_percent) / 100.0)
    )
    return apply_ltbi_state_edit(
        config,
        baseline_recent_proportion=baseline_fraction,
        transition_rate_per_year=transition_rate_per_year,
        source=source,
        status=status,
        notes=notes,
    )


def apply_ltbi_state_development_compatibility(config: dict[str, Any]) -> dict[str, Any]:
    return enable_development_compatibility_mode(config)


def ltbi_state_display_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    state = resolve_ltbi_state_assumptions(config)
    return [
        {"Assumption": "LTBI-state model", "Value": state["transitionModel"]},
        {
            "Assumption": "Baseline recent-LTBI proportion",
            "Value": state["baselineRecentLTBIProportion"],
        },
        {
            "Assumption": "Recent-to-remote transition rate",
            "Value": state["recentToRemoteTransitionRatePerYear"],
        },
        {
            "Assumption": "Implied mean residence time",
            "Value": state["impliedMeanResidenceTimeYears"],
        },
        {"Assumption": "State definition", "Value": state["stateDefinition"]},
        {"Assumption": "Source", "Value": state["source"]},
        {"Assumption": "Status", "Value": state["status"]},
        {
            "Assumption": "Baseline recent-LTBI source",
            "Value": state["baselineRecentLTBIProportionSource"],
        },
        {
            "Assumption": "Baseline recent-LTBI status",
            "Value": state["baselineRecentLTBIProportionStatus"],
        },
        {
            "Assumption": "Baseline recent-LTBI derivation method",
            "Value": state["baselineRecentLTBIDerivationMethod"],
        },
        {"Assumption": "Transition model source", "Value": state["transitionModelSource"]},
        {"Assumption": "Transition model status", "Value": state["transitionModelStatus"]},
        {
            "Assumption": "Development compatibility mode",
            "Value": state["developmentCompatibilityMode"],
        },
        {"Assumption": "Provisional result", "Value": state["provisional"]},
        {"Assumption": "Warning", "Value": state["warning"] or ""},
    ]


def prevalence_source_label(value: Any) -> str:
    return "APY default" if value is None or value == [] else "Custom"
