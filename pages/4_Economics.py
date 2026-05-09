from __future__ import annotations

from copy import deepcopy

import streamlit as st

from app.display import arrow_safe_dataframe
from app.state import (
    get_backend,
    init_session_state,
    mark_economics_changed,
    mark_economics_completed,
    record_message,
    sync_backend_status,
)


init_session_state()
backend = get_backend()

st.title("Economics")
st.caption("First economics workflow backed by MATLAB.")

ECONOMICS_WIDGET_KEYS = [
    "econ_currency_code",
    "econ_price_year",
    "econ_location_label",
    "econ_active_tb_cost",
    "econ_setup_total",
    "econ_running_total",
]

ECONOMICS_WIDGET_FIELDS = {
    "econ_currency_code": ("metadata", "currencyCode"),
    "econ_price_year": ("metadata", "priceYear"),
    "econ_location_label": ("metadata", "locationLabel"),
    "econ_active_tb_cost": ("costs", "activeTBDiseasePerCase"),
    "econ_setup_total": ("costs", "programSetupTotal"),
    "econ_running_total": ("costs", "programRunningTotal"),
}


def ensure_nested(config: dict, *path: str) -> dict:
    current = config
    for key in path:
        value = current.get(key)
        if not isinstance(value, dict):
            value = {}
            current[key] = value
        current = value
    return current


def optional_number_input(label: str, key: str) -> None:
    st.text_input(label, key=key)


def widget_value(value: object) -> str:
    if isinstance(value, list) and not value:
        return ""
    if value is None:
        return ""
    return str(value)


def sync_econ_widgets_from_config(config: dict) -> None:
    for key, (section, field) in ECONOMICS_WIDGET_FIELDS.items():
        value = config.get(section, {}).get(field)
        st.session_state[key] = widget_value(value)


def sync_econ_widgets_if_missing(config: dict) -> None:
    if any(key not in st.session_state for key in ECONOMICS_WIDGET_KEYS):
        sync_econ_widgets_from_config(config)


def parse_optional_number(key: str, integer: bool = False) -> int | float | list:
    raw = str(st.session_state.get(key, "")).strip()
    if raw == "":
        return []
    value = float(raw)
    if integer and value.is_integer():
        return int(value)
    return value


def economics_config_from_widgets(base_config: dict) -> dict:
    config = deepcopy(base_config)
    ensure_nested(config, "metadata")
    ensure_nested(config, "costs")
    config["metadata"]["currencyCode"] = str(
        st.session_state.get("econ_currency_code", "")
    ).strip()
    config["metadata"]["priceYear"] = parse_optional_number(
        "econ_price_year",
        integer=True,
    )
    config["metadata"]["locationLabel"] = str(
        st.session_state.get("econ_location_label", "")
    ).strip()
    config["costs"]["activeTBDiseasePerCase"] = parse_optional_number(
        "econ_active_tb_cost"
    )
    config["costs"]["programSetupTotal"] = parse_optional_number("econ_setup_total")
    config["costs"]["programRunningTotal"] = parse_optional_number("econ_running_total")
    return config


def economics_overview_rows(config: dict) -> list[dict[str, object]]:
    metadata = config.get("metadata", {})
    costs = config.get("costs", {})
    return [
        {"field": "currencyCode", "value": metadata.get("currencyCode")},
        {"field": "priceYear", "value": metadata.get("priceYear")},
        {"field": "locationLabel", "value": metadata.get("locationLabel")},
        {
            "field": "activeTBDiseasePerCase",
            "value": costs.get("activeTBDiseasePerCase"),
        },
        {"field": "programSetupTotal", "value": costs.get("programSetupTotal")},
        {"field": "programRunningTotal", "value": costs.get("programRunningTotal")},
    ]


cols = st.columns(2)
if cols[0].button("Load economics defaults", type="primary"):
    try:
        econ_config = backend.default_economics_config()
        st.session_state["economics_config"] = econ_config
        st.session_state["economics_results"] = None
        st.session_state["dirty_economics"] = False
        sync_backend_status(backend.status())
        sync_econ_widgets_from_config(econ_config)
        st.rerun()
    except Exception as exc:
        message = f"Could not load economics defaults: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if cols[1].button("Load KWAB150 preset"):
    try:
        econ_config = backend.economics_preset_kwab150()
        st.session_state["economics_config"] = econ_config
        st.session_state["economics_results"] = None
        st.session_state["dirty_economics"] = False
        sync_backend_status(backend.status())
        sync_econ_widgets_from_config(econ_config)
        st.rerun()
    except Exception as exc:
        message = f"Could not load KWAB150 preset: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

config = st.session_state.get("config")
results_bundle = st.session_state.get("results_bundle")
econ_config = st.session_state.get("economics_config")

if not econ_config:
    st.info("Load economics defaults or the KWAB150 preset to begin.")
    st.stop()

sync_econ_widgets_if_missing(econ_config)

if st.session_state.get("dirty_economics") and st.session_state.get("economics_results"):
    st.warning("Economics results are stale because model results or economics inputs changed.")

st.subheader("Edit Economics Inputs")
metadata = ensure_nested(econ_config, "metadata")
costs = ensure_nested(econ_config, "costs")

with st.form("economics_edits"):
    st.text_input(
        "Currency code",
        key="econ_currency_code",
    )
    optional_number_input(
        "Price year",
        "econ_price_year",
    )
    st.text_input(
        "Location label",
        key="econ_location_label",
    )
    optional_number_input(
        "Active TB disease cost per case",
        "econ_active_tb_cost",
    )
    optional_number_input(
        "Program setup total",
        "econ_setup_total",
    )
    optional_number_input(
        "Program running total",
        "econ_running_total",
    )
    submitted = st.form_submit_button("Apply economics edits")

if submitted:
    try:
        updated_config = economics_config_from_widgets(econ_config)
        changed = updated_config != econ_config
        st.session_state["economics_config"] = updated_config
        econ_config = updated_config
        if changed:
            mark_economics_changed()
            st.success("Economics edits applied.")
        else:
            st.info("No economics fields changed.")
    except ValueError as exc:
        st.error(f"Invalid economics number: {exc}")

st.subheader("Current Economics Config")
st.dataframe(
    arrow_safe_dataframe(economics_overview_rows(econ_config)),
    width="content",
    hide_index=True,
)

can_run = bool(config and results_bundle and not st.session_state.get("results_stale"))
if not config:
    st.info("Load a scenario before running economics.")
elif not results_bundle:
    st.info("Run the model before running economics.")
elif st.session_state.get("results_stale"):
    st.warning("Rerun the model before running economics so the economics inputs match current results.")

if st.button("Run economics", type="primary", disabled=not can_run):
    try:
        econ = backend.run_economics_for_config(config, econ_config)
        st.session_state["economics_results"] = econ
        mark_economics_completed()
        sync_backend_status(backend.status())
        st.success("Economics run completed.")
    except Exception as exc:
        message = f"Economics run failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

econ_results = st.session_state.get("economics_results")
if econ_results:
    st.subheader("Economics Summary")
    summary_rows = econ_results.get("summaryRows") or []
    if summary_rows:
        st.dataframe(arrow_safe_dataframe(summary_rows), width="stretch")
    else:
        st.info("No economics summary rows were returned.")

    st.subheader("Economics Status")
    status = econ_results.get("status", {})
    st.json(
        {
            "last_economics_run_at": st.session_state.get("last_economics_run_at"),
            "isComplete": status.get("isComplete"),
            "missingInputs": status.get("missingInputs"),
            "notCalculated": status.get("notCalculated"),
            "messages": status.get("messages"),
            "partialCalculations": status.get("partialCalculations"),
        },
        expanded=False,
    )
