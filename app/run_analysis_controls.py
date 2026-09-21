from __future__ import annotations

from copy import deepcopy
from typing import Any

from engine.apy.infection_history import configure_compatibility_reference_assumptions
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions


TECHNICAL_DEMONSTRATION_ROUTE = "technical_demonstration"
REVIEWED_ASSUMPTION_ROUTE = "reviewed_assumption"


def recent_ltbi_decision_required(config: dict[str, Any]) -> bool:
    state = resolve_ltbi_state_assumptions(config)
    if state.get("baselineRecentLTBIDerivationMethod") == "infection_history_trajectory":
        return False
    return state.get("baselineRecentLTBIProportion") is None


def prepare_run_config_for_recent_ltbi_route(
    config: dict[str, Any],
    *,
    selected_route: str | None,
) -> dict[str, Any]:
    del selected_route
    return configure_compatibility_reference_assumptions(deepcopy(config))


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
