"""Session state for the general application (kept separate from other workflows).

Result currency is decided by comparing hashes, not by flags:

* epidemiological results are current when the epidemiological hash of the
  configuration built from the current profile, intervention and analysis settings
  equals the hash recorded with the results (descriptive incidence and trend edits
  do not change it);
* economic results are current when both the epidemiological results and the
  economics configuration are unchanged, so cost-only edits never rerun epidemiology.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

import streamlit as st

from engine.profiles.demonstration import build_demonstration_profile
from engine.profiles.engine_mapping import (
    ProfileMappingError,
    build_engine_config,
    default_analysis,
    default_intervention,
    epidemiological_config_hash,
)
from engine.profiles.population_profile import PopulationProfile
from engine.who_incidence.snapshot import IncidenceSnapshot, SnapshotError, load_bundled_snapshot
from engine.who_incidence.trend import TrendSettings


PROFILE_KEY = "general_profile"
INTERVENTION_KEY = "general_intervention"
ANALYSIS_KEY = "general_analysis"
RESULTS_KEY = "general_results_bundle"
RESULTS_CONFIG_KEY = "general_results_config"
RESULTS_HASH_KEY = "general_results_epi_hash"
RESULTS_CACHE_KEY = "general_results_cache"
ECONOMICS_KEY = "general_economics_results"
ECONOMICS_CONFIG_KEY = "general_economics_config"
ECONOMICS_BASIS_KEY = "general_economics_basis"
TREND_KEY = "general_trend_settings"
CANDIDATE_KEY = "general_country_candidate"
EDITOR_VERSION_KEY = "general_editor_version"
RESULTS_CACHE_SIZE = 3
WIDGET_KEYS = (
    "general_country_candidate",
    "general_country_candidate_memory",
    "general_population_input",
    "general_test_type",
    "general_regimen",
    "general_coverage",
    "general_strategy",
    "general_window",
    "general_analysis_type",
    "general_n_sims",
    "general_seed",
    "general_confirm_long_run",
    "general_trend_method",
    "general_trend_window",
    "general_trend_start",
    "general_trend_end",
    "general_trend_covid",
    "general_trend_excluded",
    "general_conflict_location",
    "general_conflict_incidence",
    "general_upload_conflict",
)


def init_general_state() -> None:
    if not isinstance(st.session_state.get(PROFILE_KEY), dict):
        st.session_state[PROFILE_KEY] = build_demonstration_profile().to_dict()
    st.session_state.setdefault(INTERVENTION_KEY, default_intervention())
    st.session_state.setdefault(ANALYSIS_KEY, default_analysis())
    st.session_state.setdefault(RESULTS_KEY, None)
    st.session_state.setdefault(RESULTS_CONFIG_KEY, None)
    st.session_state.setdefault(RESULTS_HASH_KEY, None)
    st.session_state.setdefault(RESULTS_CACHE_KEY, {})
    st.session_state.setdefault(ECONOMICS_KEY, None)
    st.session_state.setdefault(ECONOMICS_BASIS_KEY, None)
    st.session_state.setdefault(TREND_KEY, TrendSettings().to_dict())
    st.session_state.setdefault(EDITOR_VERSION_KEY, 0)
    if not isinstance(st.session_state.get(ECONOMICS_CONFIG_KEY), dict):
        from app.general.economics import default_economics_config

        st.session_state[ECONOMICS_CONFIG_KEY] = default_economics_config()


def get_profile() -> PopulationProfile:
    return PopulationProfile.from_dict(st.session_state[PROFILE_KEY])


def set_profile(profile: PopulationProfile) -> None:
    st.session_state[PROFILE_KEY] = profile.to_dict()


def get_intervention() -> dict[str, Any]:
    return dict(st.session_state[INTERVENTION_KEY])


def set_intervention(intervention: dict[str, Any]) -> None:
    st.session_state[INTERVENTION_KEY] = dict(intervention)


def get_analysis() -> dict[str, Any]:
    return dict(st.session_state[ANALYSIS_KEY])


def set_analysis(analysis: dict[str, Any]) -> None:
    st.session_state[ANALYSIS_KEY] = dict(analysis)


def get_trend_settings() -> TrendSettings:
    return TrendSettings.from_dict(st.session_state[TREND_KEY])


def set_trend_settings(settings: TrendSettings) -> None:
    st.session_state[TREND_KEY] = settings.to_dict()


def current_engine_config() -> tuple[dict[str, Any] | None, str | None]:
    """Configuration for the current inputs, or (None, reason) when not runnable."""
    key = _inputs_key()
    cached = st.session_state.get("general_config_cache")
    if isinstance(cached, dict) and cached.get("key") == key:
        return cached["config"], cached["error"]
    snapshot = bundled_snapshot()
    try:
        config = build_engine_config(
            get_profile(),
            intervention=get_intervention(),
            analysis=get_analysis(),
            snapshot_manifest=snapshot.manifest if snapshot else None,
        )
        error = None
    except ProfileMappingError as exc:
        config, error = None, str(exc)
    st.session_state["general_config_cache"] = {"key": key, "config": config, "error": error}
    return config, error


def results_status() -> str:
    """'none', 'current' or 'stale' for the epidemiological results."""
    if not st.session_state.get(RESULTS_KEY):
        return "none"
    config, _ = current_engine_config()
    if config is None:
        return "stale"
    return "current" if epidemiological_config_hash(config) == st.session_state.get(RESULTS_HASH_KEY) else "stale"


def cached_results_for(config: dict[str, Any]) -> dict[str, Any] | None:
    return (st.session_state.get(RESULTS_CACHE_KEY) or {}).get(epidemiological_config_hash(config))


def store_results(bundle: dict[str, Any], config: dict[str, Any]) -> None:
    epi_hash = epidemiological_config_hash(config)
    st.session_state[RESULTS_KEY] = bundle
    st.session_state[RESULTS_CONFIG_KEY] = config
    st.session_state[RESULTS_HASH_KEY] = epi_hash
    st.session_state[ECONOMICS_KEY] = None
    st.session_state[ECONOMICS_BASIS_KEY] = None
    st.session_state["general_last_run_at"] = datetime.now(timezone.utc).isoformat()
    cache = dict(st.session_state.get(RESULTS_CACHE_KEY) or {})
    cache[epi_hash] = {"bundle": bundle, "config": config}
    while len(cache) > RESULTS_CACHE_SIZE:
        cache.pop(next(iter(cache)))
    st.session_state[RESULTS_CACHE_KEY] = cache


def economics_status() -> str:
    """'none', 'current' or 'stale' for the economic results."""
    if not st.session_state.get(ECONOMICS_KEY):
        return "none"
    from app.general.economics import economics_config_hash

    basis = st.session_state.get(ECONOMICS_BASIS_KEY) or {}
    current = {
        "epi": st.session_state.get(RESULTS_HASH_KEY),
        "econ": economics_config_hash(st.session_state[ECONOMICS_CONFIG_KEY]),
    }
    return "current" if basis == current and results_status() == "current" else "stale"


def store_economics(results: dict[str, Any]) -> None:
    from app.general.economics import economics_config_hash

    st.session_state[ECONOMICS_KEY] = results
    st.session_state[ECONOMICS_BASIS_KEY] = {
        "epi": st.session_state.get(RESULTS_HASH_KEY),
        "econ": economics_config_hash(st.session_state[ECONOMICS_CONFIG_KEY]),
    }


def restore_demonstration_defaults() -> None:
    """Reset profile, country data, intervention, analysis, trend and cost settings."""
    from app.general.economics import default_economics_config

    st.session_state[PROFILE_KEY] = build_demonstration_profile().to_dict()
    st.session_state[INTERVENTION_KEY] = default_intervention()
    st.session_state[ANALYSIS_KEY] = default_analysis()
    st.session_state[TREND_KEY] = TrendSettings().to_dict()
    st.session_state[ECONOMICS_CONFIG_KEY] = default_economics_config()
    for key in WIDGET_KEYS:
        st.session_state.pop(key, None)
    st.session_state.pop("general_config_cache", None)
    st.session_state[EDITOR_VERSION_KEY] = int(st.session_state.get(EDITOR_VERSION_KEY, 0)) + 1
    st.session_state["general_restore_message"] = (
        "Demonstration defaults restored. Country data and all user-defined values were cleared."
    )


def bundle_outcome_records(bundle: dict[str, Any] | None) -> dict[str, Any] | None:
    """Per-person outcome records from an engine results bundle (internal key)."""
    technical = (bundle or {}).get("technical") or {}
    records = technical.get("eventLedger")
    return records if isinstance(records, dict) else None


def bundled_snapshot() -> IncidenceSnapshot | None:
    """Installed snapshot, or None with a maintainer-facing message recorded."""
    try:
        snapshot = load_bundled_snapshot()
    except SnapshotError as exc:
        st.session_state["general_snapshot_error"] = str(exc)
        return None
    st.session_state.pop("general_snapshot_error", None)
    return snapshot


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


def _inputs_key() -> str:
    payload = {
        "profile": st.session_state.get(PROFILE_KEY),
        "intervention": st.session_state.get(INTERVENTION_KEY),
        "analysis": st.session_state.get(ANALYSIS_KEY),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()
