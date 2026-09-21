from __future__ import annotations

from datetime import datetime, timezone
import os
from pathlib import Path
from typing import Any

import streamlit as st

from adapters.matlab_backend import MatlabBackend
from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend
from engine.apy.infection_history import (
    configure_compatibility_reference_assumptions,
    is_experimental_infection_history_config,
    is_experimental_infection_history_results,
)

SUPPORTED_SA_HEALTH_ANALYSIS_BASIS = "sa_health_matlab_v9_compatibility_reference"
SUPPORTED_SA_HEALTH_NATURAL_HISTORY_SEMANTICS = "matlab_v9_implicit_early_late"


REFERENCE_ONLY_STALE_RESULTS_MESSAGE = (
    "Previous results used an analysis pathway that is not available in this "
    "SA Health version. Please run the analysis again."
)


def init_session_state() -> None:
    """Initialize JSON-like Streamlit app state."""
    defaults: dict[str, Any] = {
        "backend_status": {
            "name": "python_apy",
            "started": False,
            "error": "",
        },
        "apy_backend_name": "python_apy",
        "config": None,
        "economics_config": None,
        "working_default_preset": None,
        "parameter_workspace": None,
        "parameter_workspace_visible": False,
        "parameter_workspace_validation": None,
        "validation_report": None,
        "load_info": None,
        "save_info": None,
        "results_bundle": None,
        "economics_results": None,
        "dirty_config": False,
        "dirty_economics": False,
        "results_stale": False,
        "last_economics_run_at": "",
        "last_run_at": "",
        "last_validated_at": "",
        "messages": [],
        "dynamic_results_bundle": None,
        "dynamic_last_run_at": "",
        "dynamic_results_stale": False,
        "dynamic_compare_rows": None,
        "dynamic_abm_compare_last_run_at": "",
        "dynamic_abm_compare_warnings": [],
        "compare_baseline_config": None,
        "compare_comparator_config": None,
        "compare_baseline_bundle": None,
        "compare_comparator_bundle": None,
        "compare_economics_config": None,
        "compare_baseline_economics_results": None,
        "compare_comparator_economics_results": None,
        "compare_baseline_validation_report": None,
        "compare_comparator_validation_report": None,
        "compare_dirty": False,
        "compare_results_stale": False,
        "compare_economics_stale": False,
        "compare_outputs_cleared": False,
        "compare_last_run_at": "",
        "compare_last_economics_run_at": "",
        "compare_selected_preset": "Custom",
        "compare_manual_status": "Not reviewed",
        "compare_manual_notes": "",
    }
    for key, value in defaults.items():
        st.session_state.setdefault(key, value)


def sanitize_reference_only_state() -> bool:
    """Remove retired infection-history state from the standard SA Health workflow.

    The research implementation remains importable for tests and development, but the
    Streamlit SA Health workflow must only run the MATLAB-v9-compatible natural history.
    Returns True when active state was changed.
    """
    changed = False

    for key in (
        "experimental_infection_history_enabled",
        "infection_history_trajectory_label",
        "infection_history_calibration",
        "infection_history_diagnostics",
        "recent_ltbi_run_route",
        "stochastic_icer_cloud",
        "stochastic_icer_cloud_source",
    ):
        if key in st.session_state:
            st.session_state.pop(key, None)
            changed = True

    config = st.session_state.get("config")
    if is_experimental_infection_history_config(config):
        try:
            from engine.apy.working_defaults import build_unified_working_default_preset

            preset = build_unified_working_default_preset()
            st.session_state["config"] = configure_compatibility_reference_assumptions(preset["config"])
            st.session_state["economics_config"] = preset["economicsConfig"]
            st.session_state["working_default_preset"] = {
                key: preset[key]
                for key in [
                    "contractVersion",
                    "presetId",
                    "presetVersion",
                    "label",
                    "sourceComponentPresets",
                    "workingDefault",
                    "referenceStatus",
                    "configurationHash",
                    "provisionalAssumptions",
                    "unresolvedAssumptions",
                ]
            }
        except Exception:
            st.session_state["config"] = configure_compatibility_reference_assumptions(config)
        st.session_state["reference_only_migration_notice"] = (
            "A configuration requested an analysis basis that is not available in this "
            "SA Health version. APY defaults have been restored."
        )
        changed = True
        try:
            from app.parameter_workspace import build_parameter_workspace

            econ = st.session_state.get("economics_config") or {}
            st.session_state["parameter_workspace"] = build_parameter_workspace(st.session_state["config"], econ)
            st.session_state["parameter_workspace_validation"] = None
        except Exception:
            st.session_state["parameter_workspace"] = None

    if has_unsupported_sa_health_results(st.session_state.get("results_bundle")):
        _clear_active_completed_analysis()
        changed = True

    for bundle_key, econ_key in (
        ("compare_baseline_bundle", "compare_baseline_economics_results"),
        ("compare_comparator_bundle", "compare_comparator_economics_results"),
    ):
        if has_unsupported_sa_health_results(st.session_state.get(bundle_key)):
            st.session_state[bundle_key] = None
            st.session_state[econ_key] = None
            st.session_state["compare_results_stale"] = False
            st.session_state["compare_economics_stale"] = False
            st.session_state["reference_only_migration_notice"] = REFERENCE_ONLY_STALE_RESULTS_MESSAGE
            changed = True

    return changed


def _clear_active_completed_analysis() -> None:
    for key in (
        "results_bundle",
        "economics_results",
        "validation_report",
        "economic_scenario_comparison",
        "decision_scenario_comparison",
        "decision_sensitivity",
        "decision_threshold",
        "decision_early_review",
    ):
        st.session_state[key] = None
    st.session_state["dirty_economics"] = False
    st.session_state["results_stale"] = False
    st.session_state["last_economics_run_at"] = ""
    st.session_state["last_run_at"] = ""
    st.session_state["reference_only_migration_notice"] = REFERENCE_ONLY_STALE_RESULTS_MESSAGE


def _result_metadata_sources(results_bundle: dict[str, Any] | None) -> list[dict[str, Any]]:
    if not isinstance(results_bundle, dict):
        return []
    sources: list[dict[str, Any]] = []
    metadata = results_bundle.get("metadata")
    if isinstance(metadata, dict):
        sources.append(metadata)
    technical = results_bundle.get("technical") if isinstance(results_bundle.get("technical"), dict) else {}
    interface_config = technical.get("interfaceConfig") if isinstance(technical.get("interfaceConfig"), dict) else {}
    if interface_config:
        sources.append(interface_config)
    event_ledger = technical.get("eventLedger") if isinstance(technical.get("eventLedger"), dict) else {}
    ledger_metadata = event_ledger.get("metadata") if isinstance(event_ledger.get("metadata"), dict) else {}
    if ledger_metadata:
        sources.append(ledger_metadata)
    return sources


def _first_nonempty_metadata_value(sources: list[dict[str, Any]], key: str) -> Any:
    for source in sources:
        value = source.get(key)
        if value not in (None, ""):
            return value
    return None


def is_supported_sa_health_analysis_basis(results_bundle: dict[str, Any] | None) -> bool:
    """Return True only for completed results with supported SA Health semantics.

    Accepted natural-history/calibration identifiers for this release are:
    - naturalHistorySemantics: matlab_v9_implicit_early_late
    - analysisBasis: sa_health_matlab_v9_compatibility_reference, when recorded

    Missing, unknown or legacy recent/remote natural-history provenance is rejected
    rather than reinterpreted.
    """
    if not isinstance(results_bundle, dict):
        return False
    if is_experimental_infection_history_results(results_bundle):
        return False
    sources = _result_metadata_sources(results_bundle)
    if not sources:
        return False
    analysis_basis = _first_nonempty_metadata_value(sources, "analysisBasis")
    natural_history = _first_nonempty_metadata_value(sources, "naturalHistorySemantics")
    if natural_history != SUPPORTED_SA_HEALTH_NATURAL_HISTORY_SEMANTICS:
        return False
    return analysis_basis in (None, SUPPORTED_SA_HEALTH_ANALYSIS_BASIS)


def has_unsupported_sa_health_results(results_bundle: dict[str, Any] | None) -> bool:
    """Return True for completed result bundles not supported in the SA Health workflow."""
    return isinstance(results_bundle, dict) and not is_supported_sa_health_analysis_basis(results_bundle)


def has_retired_infection_history_results(results_bundle: dict[str, Any] | None) -> bool:
    """Return True for result bundles unavailable in the reference-only workflow."""
    return has_unsupported_sa_health_results(results_bundle)



def get_backend_name() -> str:
    name = str(st.session_state.get("apy_backend_name", "python_apy"))
    if name == "matlab" and not matlab_backend_enabled():
        st.session_state["apy_backend_name"] = "python_apy"
        return "python_apy"
    return name


def set_backend_name(name: str) -> None:
    if name not in {"matlab", "python_apy"}:
        raise ValueError(f"Unsupported APY backend: {name}")
    if name == "matlab" and not matlab_backend_enabled():
        raise ValueError(
            "MATLAB reference backend is unavailable in this deployment. "
            "Set APY_ENABLE_MATLAB_BACKEND=true in a local validation environment to enable it."
        )
    if st.session_state.get("apy_backend_name") == name:
        return
    st.session_state["apy_backend_name"] = name
    clear_apy_outputs_for_backend_switch()


def clear_apy_outputs_for_backend_switch() -> None:
    st.session_state["validation_report"] = None
    st.session_state["results_bundle"] = None
    st.session_state["economics_results"] = None
    st.session_state["dirty_config"] = False
    st.session_state["dirty_economics"] = False
    st.session_state["results_stale"] = False
    st.session_state["last_economics_run_at"] = ""
    st.session_state["last_run_at"] = ""
    st.session_state["last_validated_at"] = ""
    st.session_state["compare_baseline_bundle"] = None
    st.session_state["compare_comparator_bundle"] = None
    st.session_state["compare_baseline_economics_results"] = None
    st.session_state["compare_comparator_economics_results"] = None
    st.session_state["compare_baseline_validation_report"] = None
    st.session_state["compare_comparator_validation_report"] = None
    st.session_state["compare_dirty"] = False
    st.session_state["compare_results_stale"] = False
    st.session_state["compare_economics_stale"] = False
    st.session_state["compare_outputs_cleared"] = True
    st.session_state["compare_last_run_at"] = ""
    st.session_state["compare_last_economics_run_at"] = ""


@st.cache_resource(show_spinner=False)
def get_matlab_backend(root: str) -> MatlabBackend:
    return MatlabBackend(Path(root))


@st.cache_resource(show_spinner=False)
def get_python_apy_backend(root: str) -> PythonApyBackend:
    return PythonApyBackend(Path(root))


def get_backend() -> MatlabBackend | PythonApyBackend:
    """Return the selected cached APY backend resource."""
    root = str(repo_root())
    if get_backend_name() == "python_apy":
        return get_python_apy_backend(root)
    return get_matlab_backend(root)


def matlab_backend_enabled() -> bool:
    return os.getenv("APY_ENABLE_MATLAB_BACKEND", "").strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }


def mark_run_completed() -> None:
    st.session_state["last_run_at"] = datetime.now(timezone.utc).isoformat()
    st.session_state["dirty_config"] = False
    st.session_state["results_stale"] = False


def mark_validation_completed() -> None:
    st.session_state["last_validated_at"] = datetime.now(timezone.utc).isoformat()
    st.session_state["dirty_config"] = False


def mark_config_changed() -> None:
    st.session_state["dirty_config"] = True
    st.session_state["validation_report"] = None
    if st.session_state.get("results_bundle"):
        st.session_state["results_stale"] = True
    if st.session_state.get("economics_results"):
        st.session_state["dirty_economics"] = True
    if st.session_state.get("compare_baseline_bundle") or st.session_state.get("compare_comparator_bundle"):
        st.session_state["compare_results_stale"] = True
        st.session_state["compare_economics_stale"] = True
        st.session_state["compare_outputs_cleared"] = False


def mark_economics_changed() -> None:
    if st.session_state.get("economics_results"):
        st.session_state["dirty_economics"] = True
    if st.session_state.get("compare_baseline_economics_results") or st.session_state.get("compare_comparator_economics_results"):
        st.session_state["compare_economics_stale"] = True


def mark_economics_completed() -> None:
    st.session_state["last_economics_run_at"] = datetime.now(timezone.utc).isoformat()
    st.session_state["dirty_economics"] = False


def mark_dynamic_run_completed() -> None:
    st.session_state["dynamic_last_run_at"] = datetime.now(timezone.utc).isoformat()
    st.session_state["dynamic_results_stale"] = False
    st.session_state["dynamic_compare_rows"] = None
    st.session_state["dynamic_abm_compare_warnings"] = []


def clear_dynamic_outputs() -> None:
    st.session_state["dynamic_results_bundle"] = None
    st.session_state["dynamic_results_stale"] = False
    st.session_state["dynamic_compare_rows"] = None
    st.session_state["dynamic_abm_compare_last_run_at"] = ""
    st.session_state["dynamic_abm_compare_warnings"] = []


def mark_dynamic_outputs_stale() -> None:
    if st.session_state.get("dynamic_results_bundle"):
        st.session_state["dynamic_results_stale"] = True


def mark_dynamic_abm_compare_completed(rows: list[dict[str, Any]], warnings: list[str]) -> None:
    st.session_state["dynamic_compare_rows"] = rows
    st.session_state["dynamic_abm_compare_warnings"] = warnings
    st.session_state["dynamic_abm_compare_last_run_at"] = datetime.now(timezone.utc).isoformat()


def sync_backend_status(status: dict[str, Any]) -> None:
    st.session_state["backend_status"] = status


def record_message(level: str, text: str) -> None:
    messages = list(st.session_state.get("messages", []))
    messages.append(
        {
            "level": level,
            "text": text,
            "created_at": datetime.now(timezone.utc).isoformat(),
        }
    )
    st.session_state["messages"] = messages[-10:]
