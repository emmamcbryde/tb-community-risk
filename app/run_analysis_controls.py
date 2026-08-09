from __future__ import annotations

from copy import deepcopy
from typing import Any

from app.epidemiology_inputs import apply_ltbi_state_development_compatibility
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions


TECHNICAL_DEMONSTRATION_ROUTE = "technical_demonstration"
REVIEWED_ASSUMPTION_ROUTE = "reviewed_assumption"


def recent_ltbi_decision_required(config: dict[str, Any]) -> bool:
    state = resolve_ltbi_state_assumptions(config)
    return state.get("baselineRecentLTBIProportion") is None


def prepare_run_config_for_recent_ltbi_route(
    config: dict[str, Any],
    *,
    selected_route: str | None,
) -> dict[str, Any]:
    if not recent_ltbi_decision_required(config):
        return deepcopy(config)
    if selected_route == TECHNICAL_DEMONSTRATION_ROUTE:
        return apply_ltbi_state_development_compatibility(config)
    raise ValueError(
        "Recent versus remote LTBI is unresolved. Choose a technical demonstration "
        "or review and enter the assumption before running the analysis."
    )


def technical_demonstration_summary(config: dict[str, Any]) -> dict[str, Any]:
    run_config = prepare_run_config_for_recent_ltbi_route(
        config,
        selected_route=TECHNICAL_DEMONSTRATION_ROUTE,
    )
    state = resolve_ltbi_state_assumptions(run_config)
    return {
        "developmentCompatibilityMode": state["developmentCompatibilityMode"],
        "baselineRecentLTBIProportion": state["baselineRecentLTBIProportion"],
        "status": state["status"],
        "baselineRecentLTBIProportionStatus": state[
            "baselineRecentLTBIProportionStatus"
        ],
        "provisional": state["provisional"],
        "warning": state["warning"],
    }
