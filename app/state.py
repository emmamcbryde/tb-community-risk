from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

import streamlit as st

from adapters.matlab_backend import MatlabBackend
from adapters.paths import repo_root


def init_session_state() -> None:
    """Initialize JSON-like Streamlit app state."""
    defaults: dict[str, Any] = {
        "backend_status": {
            "name": "matlab",
            "started": False,
            "error": "",
        },
        "config": None,
        "economics_config": None,
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


@st.cache_resource(show_spinner=False)
def get_backend() -> MatlabBackend:
    """Return the cached backend resource for this Streamlit process."""
    return MatlabBackend(repo_root())


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
