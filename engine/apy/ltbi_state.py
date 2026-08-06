from __future__ import annotations

from copy import deepcopy
from typing import Any

import numpy as np


LTBI_STATE_MODEL = "continuous_markov_recent_remote"
RECENT_LTBI_DEFINITION = (
    "fast/recent latent state with mean residence time of five years before "
    "transition to remote latent infection"
)
RECENT_TO_REMOTE_RATE_PER_YEAR = 1.0 / 5.0
DEVELOPMENT_COMPATIBILITY_RECENT_PROPORTION = 0.0
BASELINE_RECENT_FRACTION_UNRESOLVED_WARNING = (
    "No validated APY-specific baseline recent-LTBI fraction was found. "
    "A clinician-ready/reference run requires an explicit recent-fraction "
    "value or an explicitly selected derivation method with source and status."
)
COMPATIBILITY_MODE_WARNING = (
    "Development compatibility mode is using a numerical placeholder for the "
    "unresolved baseline recent-LTBI fraction. Results are provisional and "
    "must not be described as clinician-ready or scientifically final."
)

BASELINE_RECENT_KEY = "baselineRecentLTBIProportion"
TRANSITION_RATE_KEY = "recentToRemoteTransitionRatePerYear"
LEGACY_KEYS = {BASELINE_RECENT_KEY, TRANSITION_RATE_KEY}
CLINICIAN_READY_LTBI_STATUSES = {"configured_reviewed", "model_derived_reviewed"}
NON_READY_LTBI_STATUSES = {
    "unresolved",
    "unresolved_development_compatibility",
    "provisional",
    "migrated_legacy_unverified",
}


def resolve_ltbi_state_assumptions(config: dict[str, Any]) -> dict[str, Any]:
    assumptions = canonicalise_ltbi_state_assumptions(config)["ltbiStateAssumptions"]
    value = assumptions.get(BASELINE_RECENT_KEY)
    compatibility_mode = bool(assumptions.get("developmentCompatibilityMode"))
    warnings = list(assumptions.get("warnings") or [])
    provisional = bool(assumptions.get("provisional") or compatibility_mode)
    if value is None or value == []:
        if not compatibility_mode:
            return _resolved(
                baseline_recent=None,
                transition_rate=assumptions.get(TRANSITION_RATE_KEY),
                assumptions=assumptions,
                warning=BASELINE_RECENT_FRACTION_UNRESOLVED_WARNING,
                provisional=True,
            )
        value = DEVELOPMENT_COMPATIBILITY_RECENT_PROPORTION
        warnings.append(COMPATIBILITY_MODE_WARNING)
    p_recent = float(value)
    if p_recent < 0 or p_recent > 1:
        raise ValueError("ltbiStateAssumptions.baselineRecentLTBIProportion must be in [0,1].")
    rate = assumptions.get(TRANSITION_RATE_KEY, RECENT_TO_REMOTE_RATE_PER_YEAR)
    rate = float(rate)
    if rate <= 0:
        raise ValueError("ltbiStateAssumptions.recentToRemoteTransitionRatePerYear must be > 0.")
    status = str(assumptions.get("status", "configured"))
    if status in NON_READY_LTBI_STATUSES:
        provisional = True
    warning = None
    if status.startswith("unresolved"):
        warnings.append(BASELINE_RECENT_FRACTION_UNRESOLVED_WARNING)
    if warnings:
        warning = " ".join(dict.fromkeys(str(item) for item in warnings if item))
    return _resolved(
        baseline_recent=p_recent,
        transition_rate=rate,
        assumptions={
            **assumptions,
            "status": status,
            "provisional": provisional or status.startswith("unresolved"),
        },
        warning=warning,
        provisional=provisional or status.startswith("unresolved"),
    )


def canonicalise_ltbi_state_assumptions(config: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(config)
    raw = out.get("ltbiStateAssumptions")
    nested_present = isinstance(raw, dict) and bool(raw)
    nested = _default_nested_assumptions()
    if nested_present:
        nested.update(deepcopy(raw))
        if "baselineRecentLTBIProportionStatus" not in raw and "status" in raw:
            nested["baselineRecentLTBIProportionStatus"] = raw["status"]
        if "baselineRecentLTBIProportionSource" not in raw and "source" in raw:
            nested["baselineRecentLTBIProportionSource"] = raw["source"]

    legacy_recent_present = BASELINE_RECENT_KEY in out
    legacy_rate_present = TRANSITION_RATE_KEY in out
    if legacy_recent_present:
        legacy_recent = out[BASELINE_RECENT_KEY]
        nested_recent = nested.get(BASELINE_RECENT_KEY)
        if nested_present and nested_recent not in (None, [], legacy_recent):
            raise ValueError(
                "Conflicting LTBI-state configuration: nested "
                "ltbiStateAssumptions.baselineRecentLTBIProportion and top-level "
                "baselineRecentLTBIProportion disagree."
            )
        if not nested_present or nested_recent in (None, []):
            nested[BASELINE_RECENT_KEY] = legacy_recent
            if not nested_present:
                nested["status"] = "migrated_legacy_unverified"
                nested["source"] = "Migrated from legacy top-level field"
                nested["provisional"] = True
    if legacy_rate_present:
        legacy_rate = out[TRANSITION_RATE_KEY]
        nested_rate = nested.get(TRANSITION_RATE_KEY)
        if nested_present and nested_rate not in (None, [], legacy_rate):
            raise ValueError(
                "Conflicting LTBI-state configuration: nested "
                "ltbiStateAssumptions.recentToRemoteTransitionRatePerYear and top-level "
                "recentToRemoteTransitionRatePerYear disagree."
            )
        if not nested_present or nested_rate in (None, []):
            nested[TRANSITION_RATE_KEY] = legacy_rate

    nested["transitionModel"] = str(
        nested.get("transitionModel") or LTBI_STATE_MODEL
    )
    nested["stateDefinition"] = str(nested.get("stateDefinition") or RECENT_LTBI_DEFINITION)
    nested["recentDefinitionYears"] = nested.get("recentDefinitionYears", None)
    rate = nested.get(TRANSITION_RATE_KEY, RECENT_TO_REMOTE_RATE_PER_YEAR)
    rate = float(rate)
    if rate <= 0:
        raise ValueError(
            "ltbiStateAssumptions.recentToRemoteTransitionRatePerYear must be > 0."
        )
    nested[TRANSITION_RATE_KEY] = rate
    nested["impliedMeanResidenceTimeYears"] = 1.0 / rate
    nested.setdefault("source", "")
    nested.setdefault("baselineRecentLTBIProportionSource", nested.get("source", ""))
    nested.setdefault("baselineRecentLTBIProportionStatus", nested.get("status", "unresolved"))
    nested.setdefault("baselineRecentLTBIDerivationMethod", "")
    nested.setdefault("transitionModelSource", nested.get("source", ""))
    nested.setdefault("transitionModelStatus", "configured_reviewed")
    if (
        nested.get("baselineRecentLTBIProportionStatus") in (None, "", "unresolved")
        and str(nested.get("status")) in CLINICIAN_READY_LTBI_STATUSES
    ):
        nested["baselineRecentLTBIProportionStatus"] = nested["status"]
    if (
        not str(nested.get("baselineRecentLTBIProportionSource") or "").strip()
        and str(nested.get("status")) in CLINICIAN_READY_LTBI_STATUSES
        and str(nested.get("source") or "").strip()
        and "APY-specific baseline recent fraction unresolved" not in str(nested.get("source"))
    ):
        nested["baselineRecentLTBIProportionSource"] = nested["source"]
    nested.setdefault("notes", "")
    nested.setdefault("status", "unresolved")
    nested.setdefault("developmentCompatibilityMode", False)
    if str(nested.get("baselineRecentLTBIProportionStatus", "")) in NON_READY_LTBI_STATUSES:
        nested["provisional"] = True
    if str(nested.get("status", "")) in NON_READY_LTBI_STATUSES:
        nested["provisional"] = True
    if "provisional" not in nested:
        nested["provisional"] = bool(
            str(nested["status"]).startswith("unresolved")
            or nested["status"] in NON_READY_LTBI_STATUSES
            or nested.get("developmentCompatibilityMode")
        )
    elif (
        str(nested["status"]) in CLINICIAN_READY_LTBI_STATUSES
        and str(nested.get("baselineRecentLTBIProportionStatus")) in CLINICIAN_READY_LTBI_STATUSES
        and not nested.get(
        "developmentCompatibilityMode"
        )
    ):
        nested["provisional"] = False
    if nested.get(BASELINE_RECENT_KEY) in ("", []):
        nested[BASELINE_RECENT_KEY] = None

    out["ltbiStateAssumptions"] = nested
    out[BASELINE_RECENT_KEY] = nested.get(BASELINE_RECENT_KEY)
    out[TRANSITION_RATE_KEY] = nested.get(TRANSITION_RATE_KEY)
    return out


def apply_ltbi_state_edit(
    config: dict[str, Any],
    *,
    baseline_recent_proportion: float | None,
    transition_rate_per_year: float,
    source: str | None = None,
    status: str = "configured",
    notes: str | None = None,
) -> dict[str, Any]:
    updated = canonicalise_ltbi_state_assumptions(config)
    nested = deepcopy(updated["ltbiStateAssumptions"])
    nested[BASELINE_RECENT_KEY] = baseline_recent_proportion
    nested[TRANSITION_RATE_KEY] = transition_rate_per_year
    nested["status"] = status
    nested["baselineRecentLTBIProportionStatus"] = status
    nested["baselineRecentLTBIProportionSource"] = source or nested.get(
        "baselineRecentLTBIProportionSource",
        "",
    )
    nested["source"] = source or nested.get("source", "")
    nested["notes"] = notes if notes is not None else nested.get("notes", "")
    nested["developmentCompatibilityMode"] = False
    nested["provisional"] = status not in CLINICIAN_READY_LTBI_STATUSES
    updated["ltbiStateAssumptions"] = nested
    updated.pop(BASELINE_RECENT_KEY, None)
    updated.pop(TRANSITION_RATE_KEY, None)
    return canonicalise_ltbi_state_assumptions(updated)


def enable_development_compatibility_mode(config: dict[str, Any]) -> dict[str, Any]:
    updated = canonicalise_ltbi_state_assumptions(config)
    nested = deepcopy(updated["ltbiStateAssumptions"])
    nested[BASELINE_RECENT_KEY] = None
    nested["developmentCompatibilityMode"] = True
    nested["status"] = "unresolved_development_compatibility"
    nested["provisional"] = True
    nested.setdefault("source", "")
    nested["notes"] = (
        str(nested.get("notes") or "").strip()
        + " Development compatibility mode explicitly selected; using 0% placeholder."
    ).strip()
    updated["ltbiStateAssumptions"] = nested
    updated.pop(BASELINE_RECENT_KEY, None)
    updated.pop(TRANSITION_RATE_KEY, None)
    return canonicalise_ltbi_state_assumptions(updated)


def is_clinician_ready_ltbi_state(config: dict[str, Any]) -> bool:
    state = resolve_ltbi_state_assumptions(config)
    return (
        state["baselineRecentLTBIProportion"] is not None
        and not state["provisional"]
        and not state["developmentCompatibilityMode"]
        and state["status"] in CLINICIAN_READY_LTBI_STATUSES
        and state["baselineRecentLTBIProportionStatus"] in CLINICIAN_READY_LTBI_STATUSES
        and bool(str(state.get("baselineRecentLTBIProportionSource") or "").strip())
    )


def require_numeric_ltbi_state_assumptions(config: dict[str, Any]) -> dict[str, Any]:
    state = resolve_ltbi_state_assumptions(config)
    if state["baselineRecentLTBIProportion"] is None:
        raise ValueError(
            "APY run requires ltbiStateAssumptions.baselineRecentLTBIProportion "
            "or explicit ltbiStateAssumptions.developmentCompatibilityMode=true."
        )
    return state


def _default_nested_assumptions() -> dict[str, Any]:
    return {
        BASELINE_RECENT_KEY: None,
        TRANSITION_RATE_KEY: RECENT_TO_REMOTE_RATE_PER_YEAR,
        "transitionModel": LTBI_STATE_MODEL,
        "stateDefinition": RECENT_LTBI_DEFINITION,
        "recentDefinitionYears": None,
        "source": (
            "Transition structure from older static/transmission-dynamic "
            "architecture; APY-specific baseline recent fraction unresolved."
        ),
        "baselineRecentLTBIProportionSource": "",
        "baselineRecentLTBIProportionStatus": "unresolved",
        "baselineRecentLTBIDerivationMethod": "",
        "transitionModelSource": (
            "Transition structure from older static/transmission-dynamic architecture."
        ),
        "transitionModelStatus": "configured_reviewed",
        "status": "unresolved",
        "notes": "",
        "developmentCompatibilityMode": False,
        "provisional": True,
        "warnings": [],
    }


def _resolved(
    *,
    baseline_recent: float | None,
    transition_rate: Any,
    assumptions: dict[str, Any],
    warning: str | None,
    provisional: bool,
) -> dict[str, Any]:
    rate = float(
        RECENT_TO_REMOTE_RATE_PER_YEAR
        if transition_rate in (None, [])
        else transition_rate
    )
    return {
        "model": assumptions.get("transitionModel", LTBI_STATE_MODEL),
        "transitionModel": assumptions.get("transitionModel", LTBI_STATE_MODEL),
        "baselineRecentLTBIProportion": baseline_recent,
        "recentToRemoteTransitionRatePerYear": rate,
        "impliedMeanResidenceTimeYears": 1.0 / rate,
        "stateDefinition": assumptions.get("stateDefinition", RECENT_LTBI_DEFINITION),
        "recentDefinitionYears": _optional_float(assumptions.get("recentDefinitionYears")),
        "source": assumptions.get("source", "Configured APY scenario field"),
        "baselineRecentLTBIProportionSource": assumptions.get(
            "baselineRecentLTBIProportionSource",
            "",
        ),
        "baselineRecentLTBIProportionStatus": assumptions.get(
            "baselineRecentLTBIProportionStatus",
            assumptions.get("status", "configured"),
        ),
        "baselineRecentLTBIDerivationMethod": assumptions.get(
            "baselineRecentLTBIDerivationMethod",
            "",
        ),
        "transitionModelSource": assumptions.get("transitionModelSource", ""),
        "transitionModelStatus": assumptions.get("transitionModelStatus", ""),
        "status": assumptions.get("status", "configured"),
        "notes": assumptions.get("notes", ""),
        "developmentCompatibilityMode": bool(
            assumptions.get("developmentCompatibilityMode", False)
        ),
        "provisional": bool(provisional),
        "warning": warning,
    }


def _optional_float(value: Any) -> float | None:
    if value in (None, "", []):
        return None
    return float(value)


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
