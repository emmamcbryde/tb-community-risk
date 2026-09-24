"""Session state for the general application (kept separate from other workflows)."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import streamlit as st

from engine.profiles.demonstration import build_demonstration_profile
from engine.profiles.engine_mapping import default_analysis, default_intervention
from engine.profiles.population_profile import PopulationProfile
from engine.who_incidence.snapshot import IncidenceSnapshot, load_bundled_snapshot


PROFILE_KEY = "general_profile"
INTERVENTION_KEY = "general_intervention"
ANALYSIS_KEY = "general_analysis"
RESULTS_KEY = "general_results_bundle"
RESULTS_CONFIG_KEY = "general_results_config"
ECONOMICS_KEY = "general_economics_results"
STALE_KEY = "general_results_stale"
EDITOR_VERSION_KEY = "general_editor_version"
WIDGET_KEYS = (
    "general_country_select",
    "general_population_input",
    "general_test_type",
    "general_regimen",
    "general_coverage",
    "general_strategy",
    "general_window",
    "general_analysis_type",
    "general_n_sims",
    "general_seed",
)


def init_general_state() -> None:
    if not isinstance(st.session_state.get(PROFILE_KEY), dict):
        st.session_state[PROFILE_KEY] = build_demonstration_profile().to_dict()
    st.session_state.setdefault(INTERVENTION_KEY, default_intervention())
    st.session_state.setdefault(ANALYSIS_KEY, default_analysis())
    st.session_state.setdefault(RESULTS_KEY, None)
    st.session_state.setdefault(RESULTS_CONFIG_KEY, None)
    st.session_state.setdefault(ECONOMICS_KEY, None)
    st.session_state.setdefault(STALE_KEY, False)
    st.session_state.setdefault(EDITOR_VERSION_KEY, 0)


def get_profile() -> PopulationProfile:
    return PopulationProfile.from_dict(st.session_state[PROFILE_KEY])


def set_profile(profile: PopulationProfile) -> None:
    previous = st.session_state.get(PROFILE_KEY)
    payload = profile.to_dict()
    if payload != previous:
        st.session_state[PROFILE_KEY] = payload
        mark_inputs_changed()


def get_intervention() -> dict[str, Any]:
    return dict(st.session_state[INTERVENTION_KEY])


def set_intervention(intervention: dict[str, Any]) -> None:
    if intervention != st.session_state.get(INTERVENTION_KEY):
        st.session_state[INTERVENTION_KEY] = dict(intervention)
        mark_inputs_changed()


def get_analysis() -> dict[str, Any]:
    return dict(st.session_state[ANALYSIS_KEY])


def set_analysis(analysis: dict[str, Any]) -> None:
    if analysis != st.session_state.get(ANALYSIS_KEY):
        st.session_state[ANALYSIS_KEY] = dict(analysis)
        mark_inputs_changed()


def mark_inputs_changed() -> None:
    if st.session_state.get(RESULTS_KEY):
        st.session_state[STALE_KEY] = True
    st.session_state[ECONOMICS_KEY] = None


def restore_demonstration_defaults() -> None:
    """Reset the profile, intervention and analysis settings; clears all overrides."""
    st.session_state[PROFILE_KEY] = build_demonstration_profile().to_dict()
    st.session_state[INTERVENTION_KEY] = default_intervention()
    st.session_state[ANALYSIS_KEY] = default_analysis()
    for key in WIDGET_KEYS:
        st.session_state.pop(key, None)
    st.session_state[EDITOR_VERSION_KEY] = int(st.session_state.get(EDITOR_VERSION_KEY, 0)) + 1
    st.session_state["general_restore_message"] = "Demonstration defaults restored. All user-defined values were cleared."
    mark_inputs_changed()


def store_results(bundle: dict[str, Any], config: dict[str, Any]) -> None:
    st.session_state[RESULTS_KEY] = bundle
    st.session_state[RESULTS_CONFIG_KEY] = config
    st.session_state[ECONOMICS_KEY] = None
    st.session_state[STALE_KEY] = False
    st.session_state["general_last_run_at"] = datetime.now(timezone.utc).isoformat()


def bundle_outcome_records(bundle: dict[str, Any] | None) -> dict[str, Any] | None:
    """Per-person outcome records from an engine results bundle (internal key)."""
    technical = (bundle or {}).get("technical") or {}
    records = technical.get("eventLedger")
    return records if isinstance(records, dict) else None


def bundled_snapshot() -> IncidenceSnapshot | None:
    try:
        return load_bundled_snapshot()
    except Exception as exc:  # the application must remain usable without the snapshot
        st.session_state["general_snapshot_error"] = str(exc)
        return None


@st.cache_resource(show_spinner=False)
def _python_backend(root: str):
    from adapters.python_apy_backend import PythonApyBackend

    return PythonApyBackend(Path(root))


def analysis_backend():
    from adapters.paths import repo_root

    return _python_backend(str(repo_root()))


def page_link(path: str, *, label: str) -> None:
    try:
        st.page_link(path, label=label)
    except Exception:
        st.caption(label)
