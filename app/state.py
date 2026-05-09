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


def mark_economics_completed() -> None:
    st.session_state["last_economics_run_at"] = datetime.now(timezone.utc).isoformat()
    st.session_state["dirty_economics"] = False


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
