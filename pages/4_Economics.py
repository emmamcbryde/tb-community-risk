from __future__ import annotations

from copy import deepcopy

import streamlit as st

from app.display import (
    arrow_safe_dataframe,
    economics_assumptions_json,
    economics_summary_csv,
    safe_download_stem,
)
from app.epidemiology_inputs import (
    apply_ltbi_state_assumption_update,
    fraction_to_percent,
)
from app.state import (
    get_backend,
    init_session_state,
    mark_economics_changed,
    mark_economics_completed,
    record_message,
    sync_backend_status,
)
from engine.apy.costing import normalise_cost_table
from engine.apy.evidence import assess_apy_reference_readiness
from engine.apy.economics import update_cost_item_original_values_from_legacy_fields
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions


init_session_state()
st.session_state["apy_backend_name"] = "python_apy"
backend = get_backend()

st.title("Evidence & Assumptions")
st.caption("Review economic perspective, cost inputs, health outcome assumptions and readiness warnings.")

ECONOMICS_WIDGET_KEYS = [
    "econ_currency_code",
    "econ_price_year",
    "econ_location_label",
    "econ_perspective",
    "econ_target_currency",
    "econ_target_price_year",
    "econ_discount_rate",
    "econ_threshold_value",
    "econ_threshold_currency",
    "econ_threshold_year",
    "econ_threshold_source",
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
    "econ_perspective": ("metadata", "perspective"),
    "econ_target_currency": ("metadata", "targetCurrency"),
    "econ_target_price_year": ("metadata", "targetPriceYear"),
    "econ_discount_rate": ("discounting", "selectedAnnualRate"),
    "econ_threshold_value": ("threshold", "value"),
    "econ_threshold_currency": ("threshold", "currency"),
    "econ_threshold_year": ("threshold", "referenceYear"),
    "econ_threshold_source": ("threshold", "source"),
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
    ensure_nested(config, "discounting")
    ensure_nested(config, "threshold")
    ensure_nested(config, "costs")
    ensure_nested(config, "costs", "test")
    ensure_nested(config, "costs", "regimen")
    config["metadata"]["currencyCode"] = str(
        st.session_state.get("econ_currency_code", "")
    ).strip()
    config["metadata"]["priceYear"] = str(
        st.session_state.get("econ_price_year", "")
    ).strip()
    config["metadata"]["locationLabel"] = str(
        st.session_state.get("econ_location_label", "")
    ).strip()
    config["metadata"]["perspective"] = str(
        st.session_state.get("econ_perspective", "")
    ).strip()
    config["metadata"]["targetCurrency"] = str(
        st.session_state.get("econ_target_currency", "")
    ).strip()
    config["metadata"]["targetPriceYear"] = str(
        st.session_state.get("econ_target_price_year", "")
    ).strip()
    config["discounting"]["selectedAnnualRate"] = parse_optional_number(
        "econ_discount_rate"
    )
    config["threshold"]["value"] = parse_optional_number("econ_threshold_value")
    config["threshold"]["currency"] = str(
        st.session_state.get("econ_threshold_currency", "")
    ).strip()
    config["threshold"]["referenceYear"] = parse_optional_number(
        "econ_threshold_year",
        integer=True,
    )
    config["threshold"]["source"] = str(
        st.session_state.get("econ_threshold_source", "")
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
    return update_cost_item_original_values_from_legacy_fields(config)


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
        {"field": "perspective", "value": metadata.get("perspective")},
        {"field": "targetCurrency", "value": metadata.get("targetCurrency")},
        {"field": "targetPriceYear", "value": metadata.get("targetPriceYear")},
        {"field": "discountRate", "value": config.get("discounting", {}).get("selectedAnnualRate")},
        {"field": "primaryHealthOutcome", "value": config.get("healthOutcome", {}).get("primary")},
        {"field": "threshold.value", "value": config.get("threshold", {}).get("value")},
        {"field": "threshold.currency", "value": config.get("threshold", {}).get("currency")},
        {"field": "threshold.referenceYear", "value": config.get("threshold", {}).get("referenceYear")},
        {"field": "threshold.source", "value": config.get("threshold", {}).get("source")},
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

if config:
    st.subheader("Recent versus remote LTBI assumption")
    ltbi_state = resolve_ltbi_state_assumptions(config)
    current_recent = ltbi_state.get("baselineRecentLTBIProportion")
    current_review_label = (
        "Reviewed model-derived estimate"
        if ltbi_state.get("baselineRecentLTBIProportionStatus") == "model_derived_reviewed"
        else "Reviewed direct evidence"
    )
    st.dataframe(
        arrow_safe_dataframe(
            [
                {
                    "Assumption": "Baseline recent-LTBI proportion",
                    "Value": current_recent,
                },
                {
                    "Assumption": "Source",
                    "Value": ltbi_state.get("baselineRecentLTBIProportionSource") or "",
                },
                {
                    "Assumption": "Provisional",
                    "Value": ltbi_state.get("provisional"),
                },
                {
                    "Assumption": "Warning",
                    "Value": ltbi_state.get("warning") or "",
                },
            ]
        ),
        width="stretch",
        hide_index=True,
    )
    with st.form("ltbi_recent_reviewed_assumption"):
        st.write(
            "Enter a reviewed recent-LTBI proportion or a documented model-derived "
            "estimate. This updates the authoritative LTBI-state assumption."
        )
        recent_percent = st.number_input(
            "Baseline recent-LTBI proportion (%)",
            min_value=0.0,
            max_value=100.0,
            value=float(fraction_to_percent(current_recent, 0.0) or 0.0),
            step=1.0,
            format="%.2f",
        )
        review_method = st.selectbox(
            "How was this value established?",
            ["Reviewed direct evidence", "Reviewed model-derived estimate"],
            index=0
            if current_review_label == "Reviewed direct evidence"
            else 1,
        )
        recent_source = st.text_input(
            "Source or derivation reference",
            value=str(ltbi_state.get("baselineRecentLTBIProportionSource") or ""),
        )
        recent_notes = st.text_area(
            "Notes",
            value=str(ltbi_state.get("notes") or ""),
        )
        save_recent = st.form_submit_button("Save recent-LTBI assumption")
    if save_recent:
        if not str(recent_source).strip():
            st.error("Provide a source or documented derivation before saving a reviewed assumption.")
        else:
            status = (
                "model_derived_reviewed"
                if review_method == "Reviewed model-derived estimate"
                else "configured_reviewed"
            )
            updated = apply_ltbi_state_assumption_update(
                config,
                baseline_recent_percent=float(recent_percent),
                transition_rate_per_year=float(
                    ltbi_state["recentToRemoteTransitionRatePerYear"]
                ),
                source=str(recent_source).strip(),
                status=status,
                notes=str(recent_notes).strip(),
            )
            st.session_state["config"] = updated
            st.session_state.pop("recent_ltbi_run_route", None)
            st.success("Recent-LTBI assumption saved.")
            st.rerun()

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
ensure_nested(econ_config, "discounting")
ensure_nested(econ_config, "threshold")

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
    st.text_input(
        "Economic perspective",
        key="econ_perspective",
        disabled=True,
    )
    st.text_input(
        "Target currency",
        key="econ_target_currency",
    )
    st.text_input(
        "Target price year",
        key="econ_target_price_year",
    )
    st.selectbox(
        "Discount rate",
        ["0.03", "0.0"],
        index=0 if str(st.session_state.get("econ_discount_rate", "0.03")) != "0.0" else 1,
        key="econ_discount_rate",
    )
    st.markdown("Primary outcome and benchmark")
    st.info("Default health outcome: DALYs averted.")
    optional_number_input(
        "Illustrative GDP-per-capita threshold value",
        "econ_threshold_value",
    )
    st.text_input(
        "Threshold currency",
        key="econ_threshold_currency",
    )
    optional_number_input(
        "Threshold reference year",
        "econ_threshold_year",
    )
    st.text_input(
        "Threshold source",
        key="econ_threshold_source",
    )
    st.caption("The GDP-per-capita benchmark is illustrative, not an official Australian funding threshold.")
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

st.subheader("Current Assumptions")
warnings = []
ledger = (results_bundle or {}).get("technical", {}).get("eventLedger", {}) if isinstance(results_bundle, dict) else {}
if results_bundle:
    if ledger and ledger.get("validation", {}).get("isValid") is True:
        st.success(
            f"Event ledger available: {ledger.get('metadata', {}).get('contractVersion')} "
            f"({ledger.get('metadata', {}).get('modelType')})."
        )
    else:
        st.warning("A valid APY event ledger is required for authoritative economics.")
for item in normalise_cost_table(econ_config.get("costItems") or []):
    if item.get("conversionStatus") != "valid":
        warnings.append(f"{item.get('costItemId')}: {item.get('conversionStatus')}")
if econ_config.get("threshold", {}).get("value") in (None, "", []):
    warnings.append("GDP-per-capita threshold value is unresolved.")
if warnings:
    st.warning("Unresolved assumptions: " + "; ".join(warnings))
st.dataframe(
    arrow_safe_dataframe(economics_overview_rows(econ_config)),
    width="content",
    hide_index=True,
)

readiness = assess_apy_reference_readiness(config or {}, econ_config)
st.subheader("APY Reference Evidence Readiness")
status_rows = [
    {"category": "epidemiology", "ready": readiness["epidemiologyReady"]},
    {"category": "cost", "ready": readiness["costReady"]},
    {"category": "DALY", "ready": readiness["dalyReady"]},
    {"category": "threshold", "ready": readiness["thresholdReady"]},
    {"category": "overall clinician-ready", "ready": readiness["overallClinicianReady"]},
]
st.dataframe(arrow_safe_dataframe(status_rows), width="content", hide_index=True)
unresolved_rows = [
    row
    for row in readiness["readinessRows"]
    if not row.get("ready")
][:12]
if unresolved_rows:
    st.warning("APY reference evidence remains unresolved or provisional.")
    st.dataframe(
        arrow_safe_dataframe(
            [
                {
                    "category": row.get("category"),
                    "assumptionId": row.get("assumptionId"),
                    "status": row.get("reviewStatus"),
                    "source": row.get("sourceCitation"),
                    "unresolvedReason": row.get("unresolvedReason"),
                }
                for row in unresolved_rows
            ]
        ),
        width="content",
        hide_index=True,
    )
st.download_button(
    "Download economics assumptions JSON",
    data=economics_assumptions_json(econ_config),
    file_name=f"{safe_download_stem(scenario_label, 'economics_assumptions')}.json",
    mime="application/json",
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
        econ = backend.run_economics(results_bundle, econ_config)
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
    metadata = econ_results.get("metadata", {})
    if metadata.get("isProvisional"):
        st.warning("Economic outputs are provisional. Do not interpret them as clinician-ready cost-effectiveness conclusions.")
    if econ_results.get("warnings"):
        st.warning("; ".join(str(item) for item in econ_results.get("warnings") or []))
    st.subheader("Economics Summary")
    summary_rows = econ_results.get("summaryRows") or []
    if summary_rows:
        status = econ_results.get("status", {})
        validation = econ_results.get("validation", {})
        complete = bool(status.get("isComplete"))
        st.dataframe(arrow_safe_dataframe(summary_rows), width="stretch")
        primary = [
            row for row in summary_rows
            if row.get("metric") in {"incrementalCost", "dalysAverted", "primaryICER_ratioOfMeans"}
            and float(row.get("discountRate", 0.03)) in {0.0, 0.03}
        ]
        if primary and complete:
            st.caption("Primary ICER uses mean paired incremental cost divided by mean paired DALYs averted; replicate ICERs are diagnostic only.")
            st.dataframe(arrow_safe_dataframe(primary), width="stretch")
        elif primary:
            st.markdown("Partial calculations - not a complete economic result")
            st.dataframe(arrow_safe_dataframe(primary), width="stretch")
        if validation:
            st.caption(
                "Complete pairs: "
                f"{validation.get('completePairedReplicates')} of "
                f"{validation.get('totalPairedReplicates')}."
            )
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
