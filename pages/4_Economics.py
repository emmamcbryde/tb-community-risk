from __future__ import annotations

from copy import deepcopy

import streamlit as st

from app.display import (
    arrow_safe_dataframe,
    economics_assumptions_json,
    economics_summary_csv,
    safe_download_stem,
)
from app.state import (
    get_backend,
    get_backend_name,
    init_session_state,
    mark_economics_changed,
    mark_economics_completed,
    record_message,
    sync_backend_status,
)


init_session_state()
backend = get_backend()
backend_name = get_backend_name()

st.title("Economics")
st.caption("First economics workflow. Full economics execution is currently backed by MATLAB.")

if backend_name == "python_apy":
    st.warning(
        "The Python APY backend can load and validate economics assumptions, "
        "but does not yet include full economics execution."
    )
    st.info("Switch to the MATLAB v9 reference backend on the Scenario page to run economics.")

ECONOMICS_WIDGET_KEYS = [
    "econ_currency_code",
    "econ_price_year",
    "econ_location_label",
    "econ_test_igra",
    "econ_test_tst",
    "econ_regimen_3hp",
    "econ_regimen_4r",
    "econ_regimen_3hr",
    "econ_regimen_6h",
    "econ_regimen_9h",
    "econ_false_positive_incremental",
    "econ_active_tb_cost",
    "econ_setup_total",
    "econ_running_total",
]

ECONOMICS_WIDGET_FIELDS = {
    "econ_currency_code": ("metadata", "currencyCode"),
    "econ_price_year": ("metadata", "priceYear"),
    "econ_location_label": ("metadata", "locationLabel"),
    "econ_test_igra": ("costs", "test", "IGRA"),
    "econ_test_tst": ("costs", "test", "TST"),
    "econ_regimen_3hp": ("costs", "regimen", "x3HP"),
    "econ_regimen_4r": ("costs", "regimen", "x4R"),
    "econ_regimen_3hr": ("costs", "regimen", "x3HR"),
    "econ_regimen_6h": ("costs", "regimen", "x6H"),
    "econ_regimen_9h": ("costs", "regimen", "x9H"),
    "econ_false_positive_incremental": (
        "costs",
        "falsePositiveIncrementalPerPerson",
    ),
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
    for key, path in ECONOMICS_WIDGET_FIELDS.items():
        value = nested_get(config, path)
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
    ensure_nested(config, "costs", "test")
    ensure_nested(config, "costs", "regimen")
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
    config["costs"]["test"]["IGRA"] = parse_optional_number("econ_test_igra")
    config["costs"]["test"]["TST"] = parse_optional_number("econ_test_tst")
    config["costs"]["regimen"]["x3HP"] = parse_optional_number("econ_regimen_3hp")
    config["costs"]["regimen"]["x4R"] = parse_optional_number("econ_regimen_4r")
    config["costs"]["regimen"]["x3HR"] = parse_optional_number("econ_regimen_3hr")
    config["costs"]["regimen"]["x6H"] = parse_optional_number("econ_regimen_6h")
    config["costs"]["regimen"]["x9H"] = parse_optional_number("econ_regimen_9h")
    config["costs"]["falsePositiveIncrementalPerPerson"] = parse_optional_number(
        "econ_false_positive_incremental"
    )
    config["costs"]["activeTBDiseasePerCase"] = parse_optional_number(
        "econ_active_tb_cost"
    )
    config["costs"]["programSetupTotal"] = parse_optional_number("econ_setup_total")
    config["costs"]["programRunningTotal"] = parse_optional_number("econ_running_total")
    return config


def nested_get(config: dict, path: tuple[str, ...]) -> object:
    current: object = config
    for key in path:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def economics_overview_rows(config: dict) -> list[dict[str, object]]:
    metadata = config.get("metadata", {})
    costs = config.get("costs", {})
    test_costs = costs.get("test", {})
    regimen_costs = costs.get("regimen", {})
    return [
        {"field": "currencyCode", "value": metadata.get("currencyCode")},
        {"field": "priceYear", "value": metadata.get("priceYear")},
        {"field": "locationLabel", "value": metadata.get("locationLabel")},
        {"field": "test.IGRA", "value": test_costs.get("IGRA")},
        {"field": "test.TST", "value": test_costs.get("TST")},
        {"field": "regimen.3HP", "value": regimen_costs.get("x3HP")},
        {"field": "regimen.4R", "value": regimen_costs.get("x4R")},
        {"field": "regimen.3HR", "value": regimen_costs.get("x3HR")},
        {"field": "regimen.6H", "value": regimen_costs.get("x6H")},
        {"field": "regimen.9H", "value": regimen_costs.get("x9H")},
        {
            "field": "falsePositiveIncrementalPerPerson",
            "value": costs.get("falsePositiveIncrementalPerPerson"),
        },
        {
            "field": "activeTBDiseasePerCase",
            "value": costs.get("activeTBDiseasePerCase"),
        },
        {"field": "programSetupTotal", "value": costs.get("programSetupTotal")},
        {"field": "programRunningTotal", "value": costs.get("programRunningTotal")},
    ]


def validate_economics_config(config: dict) -> dict | None:
    validator = getattr(backend, "validate_economics_config", None)
    if not callable(validator):
        return None
    return validator(config)


def display_economics_config_issues(report: dict | None) -> bool:
    if not report:
        return True
    errors = report.get("errors") or []
    warnings = report.get("warnings") or []
    if errors:
        st.error("Economics assumptions have validation errors.")
        st.dataframe(arrow_safe_dataframe(errors), width="stretch", hide_index=True)
    if warnings:
        st.warning("Economics assumptions have validation warnings.")
        st.dataframe(arrow_safe_dataframe(warnings), width="stretch", hide_index=True)
    return bool(report.get("isValid", not errors))


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
scenario_label = (results_bundle or {}).get("metadata", {}).get("scenarioLabel")

if not econ_config:
    st.info("Load economics defaults or the KWAB150 preset to begin.")
    st.stop()

sync_econ_widgets_if_missing(econ_config)

if st.session_state.get("dirty_economics") and st.session_state.get("economics_results"):
    st.warning("Economics results are stale because model results or economics inputs changed.")

st.subheader("Edit Economics Inputs")
metadata = ensure_nested(econ_config, "metadata")
costs = ensure_nested(econ_config, "costs")
ensure_nested(econ_config, "costs", "test")
ensure_nested(econ_config, "costs", "regimen")

with st.form("economics_edits"):
    st.markdown("Metadata")
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
    st.markdown("Test costs")
    optional_number_input(
        "IGRA test cost",
        "econ_test_igra",
    )
    optional_number_input(
        "TST cost",
        "econ_test_tst",
    )
    st.markdown("Regimen costs")
    optional_number_input(
        "3HP regimen cost",
        "econ_regimen_3hp",
    )
    optional_number_input(
        "4R regimen cost",
        "econ_regimen_4r",
    )
    optional_number_input(
        "3HR regimen cost",
        "econ_regimen_3hr",
    )
    optional_number_input(
        "6H regimen cost",
        "econ_regimen_6h",
    )
    optional_number_input(
        "9H regimen cost",
        "econ_regimen_9h",
    )
    st.markdown("Program and disease costs")
    optional_number_input(
        "False-positive incremental cost per person",
        "econ_false_positive_incremental",
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

validation_report = validate_economics_config(econ_config)
config_is_valid = display_economics_config_issues(validation_report)

st.subheader("Current Assumptions")
st.dataframe(
    arrow_safe_dataframe(economics_overview_rows(econ_config)),
    width="content",
    hide_index=True,
)
st.download_button(
    "Download economics assumptions JSON",
    data=economics_assumptions_json(econ_config),
    file_name=f"{safe_download_stem(scenario_label, 'economics_assumptions')}.json",
    mime="application/json",
)

can_run = bool(config_is_valid and config and results_bundle and not st.session_state.get("results_stale"))
if not config:
    st.info("Load a scenario before running economics.")
elif not results_bundle:
    st.info("Run the model before running economics.")
elif st.session_state.get("results_stale"):
    st.warning("Rerun the model before running economics so the economics inputs match current results.")
elif not config_is_valid:
    st.warning("Fix economics assumption validation errors before running economics.")

if st.button("Run economics", type="primary", disabled=not can_run):
    try:
        econ = backend.run_economics_for_config(config, econ_config)
        st.session_state["economics_results"] = econ
        mark_economics_completed()
        sync_backend_status(backend.status())
        st.success("Economics run completed.")
    except NotImplementedError as exc:
        message = str(exc)
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)
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
        st.markdown("Downloads")
        st.download_button(
            "Download economics summary CSV",
            data=economics_summary_csv(econ_results),
            file_name=f"{safe_download_stem(scenario_label, 'economics_summary')}.csv",
            mime="text/csv",
        )
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
else:
    st.info("Run economics to enable the economics summary CSV download.")
