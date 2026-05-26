from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import json
import math

import altair as alt
import pandas as pd
import streamlit as st

from app.display import arrow_safe_dataframe, safe_download_stem
from app.state import get_backend, init_session_state, record_message, sync_backend_status
from engine.apy.regimen import default_regimen_library, get_regimen_from_library


init_session_state()
backend = get_backend()

COMPARE_FIELDS = [
    "testType",
    "regimen",
    "screeningStrategy",
    "screenCoverage",
    "pStartTPT",
    "regimenPComplete",
    "regimenADRstop",
    "regimenEffFull",
]

TEST_OPTIONS = ["IGRA", "TST"]
REGIMEN_OPTIONS = ["3HP", "4R", "3HR", "6H", "9H"]
STRATEGY_OPTIONS = ["random", "ltbi", "cure", "prevent"]
DEFAULT_RESOLVED_FIELDS = {
    "pStartTPT",
    "regimenPComplete",
    "regimenADRstop",
    "regimenEffFull",
}
STRATEGY_FIELD_MAP = {
    "pStartTPT": "pStartTPT",
    "regimenPComplete": "pComplete",
    "regimenADRstop": "pADRstop",
    "regimenEffFull": "effFull",
}
STRESS_TEST_PRESETS = {
    "Custom": {
        "baseline": {},
        "comparator": {},
        "note": "Manual comparison setup.",
    },
    "IGRA vs TST": {
        "baseline": {"testType": "IGRA"},
        "comparator": {"testType": "TST"},
        "note": "Check test-specific yield, false positives, treatment starts, and testing cost.",
    },
    "3HP vs 4R": {
        "baseline": {"regimen": "3HP"},
        "comparator": {"regimen": "4R"},
        "note": "Check regimen-driven completion, benefit, adverse-stop behavior, and treatment cost.",
    },
    "Low vs high coverage": {
        "baseline": {"screenCoverage": 0.2},
        "comparator": {"screenCoverage": 0.8},
        "note": "Check scaling: screened, starts, completions, benefits, and costs should generally increase.",
    },
    "Cascade improvement": {
        "baseline": {},
        "comparator": {
            "pStartTPT": 0.9,
            "regimenPComplete": 0.9,
            "regimenADRstop": 0.02,
            "regimenEffFull": 0.9,
        },
        "note": "Check improved starts, completions, infection cured, and active TB prevented.",
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def optional_text(value: object) -> str:
    if value is None or (isinstance(value, list) and not value):
        return ""
    return str(value)


def parse_optional_float(raw: str) -> float | list:
    raw = raw.strip()
    if raw == "":
        return []
    return float(raw)


def choice_index(value: object, options: list[str], fallback: int = 0) -> int:
    try:
        return options.index(str(value))
    except ValueError:
        return fallback


def clone_config(config: dict | None) -> dict | None:
    if config is None:
        return None
    return deepcopy(config)


def apply_overrides(config: dict, overrides: dict[str, object]) -> dict:
    updated = clone_config(config) or {}
    updated.update(overrides)
    return updated


def preset_payload(
    assumptions_rows: list[dict[str, object]],
    outcome_diff: list[dict[str, object]],
    outcome_warnings: list[str],
    econ_diff: list[dict[str, object]],
    econ_warnings: list[str],
) -> dict[str, object]:
    return {
        "exportedAt": utc_now(),
        "selectedPreset": st.session_state.get("compare_selected_preset", "Custom"),
        "manualStatus": st.session_state.get("compare_manual_status", "Not reviewed"),
        "manualNotes": st.session_state.get("compare_manual_notes", ""),
        "stale": {
            "compareResults": bool(st.session_state.get("compare_results_stale")),
            "compareEconomics": bool(st.session_state.get("compare_economics_stale")),
            "compareOutputsCleared": bool(st.session_state.get("compare_outputs_cleared")),
        },
        "timestamps": {
            "comparisonRunAt": st.session_state.get("compare_last_run_at", ""),
            "economicsRunAt": st.session_state.get("compare_last_economics_run_at", ""),
        },
        "baselineConfig": st.session_state.get("compare_baseline_config"),
        "comparatorConfig": st.session_state.get("compare_comparator_config"),
        "compareEconomicsConfig": st.session_state.get("compare_economics_config"),
        "selectedAssumptionsDiff": assumptions_rows,
        "outcomesDiff": outcome_diff,
        "outcomeWarnings": outcome_warnings,
        "economicsDiff": econ_diff,
        "economicsWarnings": econ_warnings,
    }


def json_download(data: dict[str, object]) -> str:
    return json.dumps(json_safe_for_export(data), indent=2, sort_keys=True)


def json_safe_for_export(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): json_safe_for_export(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe_for_export(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def attributable_risk_payload_for_bundle(backend_obj: object, bundle: dict | None) -> dict[str, object]:
    """Return a JSON/CSV-safe attributable-risk payload for one result bundle."""
    run_attributable_risk = getattr(backend_obj, "run_attributable_risk", None)
    if not callable(run_attributable_risk):
        return {
            "status": "unsupported",
            "source": "compare_page_backend_capability_check",
            "calculatedRows": [],
            "missingInputs": [],
            "unsupportedMetrics": [
                {
                    "metric": "attributableRisk",
                    "reason": (
                        "The selected backend does not expose "
                        "run_attributable_risk."
                    ),
                    "backendCapability": "run_attributable_risk",
                }
            ],
            "messages": [
                (
                    "Attributable-risk comparison is unsupported for the "
                    "selected backend because it does not expose "
                    "run_attributable_risk."
                )
            ],
        }

    payload = run_attributable_risk(bundle or {})
    if isinstance(payload, dict):
        return json_safe_for_export(payload)
    return {
        "status": "unsupported",
        "source": "compare_page_backend_capability_check",
        "calculatedRows": [],
        "missingInputs": [],
        "unsupportedMetrics": [
            {
                "metric": "attributableRisk",
                "reason": "Backend run_attributable_risk did not return a dict payload.",
                "backendCapability": "run_attributable_risk",
            }
        ],
        "messages": [
            "Attributable-risk comparison is unsupported because the backend returned an invalid payload."
        ],
    }


def store_attributable_risk_compare_payloads(
    backend_obj: object,
    baseline_bundle: dict | None,
    comparator_bundle: dict | None,
) -> tuple[dict[str, object], dict[str, object]]:
    """Store baseline/comparator attributable-risk payloads for Compare display."""
    baseline_payload = attributable_risk_payload_for_bundle(backend_obj, baseline_bundle)
    comparator_payload = attributable_risk_payload_for_bundle(backend_obj, comparator_bundle)
    st.session_state["compare_baseline_attributable_risk_payload"] = baseline_payload
    st.session_state["compare_comparator_attributable_risk_payload"] = comparator_payload
    return baseline_payload, comparator_payload


def mark_compare_dirty() -> None:
    st.session_state["compare_dirty"] = True
    if has_compare_outputs():
        reset_compare_outputs(outputs_cleared=True)


def has_compare_outputs() -> bool:
    return bool(
        st.session_state.get("compare_baseline_bundle")
        or st.session_state.get("compare_comparator_bundle")
        or st.session_state.get("compare_baseline_economics_results")
        or st.session_state.get("compare_comparator_economics_results")
        or st.session_state.get("compare_baseline_attributable_risk_payload")
        or st.session_state.get("compare_comparator_attributable_risk_payload")
    )


def reset_compare_outputs(outputs_cleared: bool = False) -> None:
    st.session_state["compare_baseline_bundle"] = None
    st.session_state["compare_comparator_bundle"] = None
    st.session_state["compare_economics_config"] = None
    st.session_state["compare_baseline_economics_results"] = None
    st.session_state["compare_comparator_economics_results"] = None
    st.session_state["compare_baseline_attributable_risk_payload"] = None
    st.session_state["compare_comparator_attributable_risk_payload"] = None
    st.session_state["compare_baseline_validation_report"] = None
    st.session_state["compare_comparator_validation_report"] = None
    st.session_state["compare_results_stale"] = False
    st.session_state["compare_economics_stale"] = False
    st.session_state["compare_outputs_cleared"] = outputs_cleared


def bundle_interface_config(bundle: dict | None) -> dict:
    if not bundle:
        return {}
    technical = bundle.get("technical")
    if not isinstance(technical, dict):
        return {}
    interface_config = technical.get("interfaceConfig")
    if isinstance(interface_config, dict):
        return interface_config
    return {}


def bundle_strategy(bundle: dict | None) -> dict:
    if not bundle:
        return {}
    headline = bundle.get("headline")
    if not isinstance(headline, dict):
        return {}
    strategy = headline.get("strategy")
    if isinstance(strategy, dict):
        return strategy
    return {}


def is_blank_default_value(value: object) -> bool:
    return value is None or value == "" or value == []


def default_resolved_compare_value(config: dict, field: str) -> object:
    """Return UI-only default values for blank fields resolved by the backend."""
    if field == "pStartTPT":
        return 0.85
    regimen_field = STRATEGY_FIELD_MAP.get(field)
    if not regimen_field:
        return config.get(field)
    try:
        regimen = get_regimen_from_library(
            str(config.get("regimen", "3HP")),
            default_regimen_library(str(config.get("partialShortCourseMode") or "threshold80")),
        )
    except Exception:
        return config.get(field)
    return regimen.get(regimen_field, config.get(field))


def editable_compare_value(config: dict, field: str) -> object:
    raw_value = config.get(field)
    if field in DEFAULT_RESOLVED_FIELDS and is_blank_default_value(raw_value):
        return default_resolved_compare_value(config, field)
    return raw_value


def populate_editable_compare_defaults(config: dict) -> dict:
    """Fill blank backend-default fields so comparator controls visibly populate."""
    updated = clone_config(config) or {}
    for field in DEFAULT_RESOLVED_FIELDS:
        if is_blank_default_value(updated.get(field)):
            updated[field] = default_resolved_compare_value(updated, field)
    return updated


def value_from_dict_case_insensitive(data: dict, field: str) -> object:
    if field in data:
        return data.get(field)
    lower_field = field.lower()
    for key, value in data.items():
        if str(key).lower() == lower_field:
            return value
    return None


def resolved_assumption_value(field: str, config: dict, bundle: dict | None) -> object:
    interface_config = bundle_interface_config(bundle)
    interface_value = value_from_dict_case_insensitive(interface_config, field)
    if not (field in DEFAULT_RESOLVED_FIELDS and is_blank_default_value(interface_value)):
        if interface_value is not None:
            return interface_value
    strategy_field = STRATEGY_FIELD_MAP.get(field)
    if strategy_field:
        strategy_value = value_from_dict_case_insensitive(bundle_strategy(bundle), strategy_field)
        if strategy_value is not None:
            return strategy_value
    raw_value = config.get(field)
    if field in DEFAULT_RESOLVED_FIELDS and is_blank_default_value(raw_value):
        return "default (resolved at run)"
    return raw_value


def assumptions_diff_rows(
    baseline: dict,
    comparator: dict,
    baseline_bundle: dict | None = None,
    comparator_bundle: dict | None = None,
) -> list[dict[str, object]]:
    fields = ["scenarioLabel", "N", "nReps", "seed", *COMPARE_FIELDS]
    rows: list[dict[str, object]] = []
    for field in fields:
        baseline_value = resolved_assumption_value(field, baseline, baseline_bundle)
        comparator_value = resolved_assumption_value(field, comparator, comparator_bundle)
        rows.append(
            {
                "field": field,
                "baseline": baseline_value,
                "comparator": comparator_value,
                "changed": baseline_value != comparator_value,
            }
        )
    return rows


def coerce_number(value: object) -> float | None:
    try:
        if value in (None, ""):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def bundle_metric_rows(bundle: dict | None) -> list[dict]:
    if not bundle:
        return []
    headline = bundle.get("headline", {})
    return headline.get("keyMetricsRows") or headline.get("summaryRows") or []


def _first_present(row: dict, fields: tuple[str, ...]) -> object:
    for field in fields:
        if field in row:
            return row.get(field)
    return None


def _index_metric_rows(
    rows: list[dict],
    value_columns: tuple[str, ...],
    key_fields: tuple[str, ...] = ("Metric",),
    require_value: bool = True,
) -> tuple[dict[str, dict[str, object]], list[str], list[str], int]:
    values: dict[str, dict[str, object]] = {}
    order: list[str] = []
    duplicates: list[str] = []
    skipped = 0
    for row in rows:
        metric = _first_present(row, key_fields)
        value = _first_present(row, value_columns)
        if metric in (None, "") or (require_value and value is None):
            skipped += 1
            continue
        metric_key = str(metric)
        if metric_key not in values:
            order.append(metric_key)
            values[metric_key] = {
                "value": value,
                "row": row,
            }
        else:
            duplicates.append(metric_key)
    return values, order, duplicates, skipped


def _compare_metric_rows(
    baseline_rows: list[dict],
    comparator_rows: list[dict],
    value_columns: tuple[str, ...],
    key_fields: tuple[str, ...] = ("Metric",),
    require_value: bool = True,
) -> tuple[list[dict[str, object]], list[str]]:
    baseline, baseline_order, baseline_duplicates, baseline_skipped = _index_metric_rows(
        baseline_rows,
        value_columns,
        key_fields,
        require_value,
    )
    comparator, comparator_order, comparator_duplicates, comparator_skipped = _index_metric_rows(
        comparator_rows,
        value_columns,
        key_fields,
        require_value,
    )
    ordered_metrics = baseline_order + [metric for metric in comparator_order if metric not in baseline]
    rows: list[dict[str, object]] = []
    for metric in ordered_metrics:
        baseline_item = baseline.get(metric, {})
        comparator_item = comparator.get(metric, {})
        base_value = baseline_item.get("value")
        comp_value = comparator_item.get("value")
        base_num = coerce_number(base_value)
        comp_num = coerce_number(comp_value)
        abs_diff = None
        rel_diff = None
        if base_num is not None and comp_num is not None:
            abs_diff = comp_num - base_num
            if base_num != 0:
                rel_diff = abs_diff / base_num
        diff_row = {
            "metric": metric,
            "baseline": base_value,
            "comparator": comp_value,
            "absoluteDifference": abs_diff,
            "relativeDifference": rel_diff,
        }
        baseline_source = baseline_item.get("row")
        comparator_source = comparator_item.get("row")
        for metadata_field in ("category", "status", "includedInTotal"):
            if isinstance(baseline_source, dict) and metadata_field in baseline_source:
                diff_row[metadata_field] = baseline_source.get(metadata_field)
            elif isinstance(comparator_source, dict) and metadata_field in comparator_source:
                diff_row[metadata_field] = comparator_source.get(metadata_field)
        rows.append(diff_row)
    warnings: list[str] = []
    duplicate_metrics = sorted(set(baseline_duplicates + comparator_duplicates))
    if duplicate_metrics:
        warnings.append(
            "Duplicate Metric values detected; the first row was used for: "
            + ", ".join(duplicate_metrics)
        )
    if baseline_skipped or comparator_skipped:
        value_column_label = "/".join(value_columns)
        key_label = "/".join(key_fields)
        warnings.append(
            f"Skipped rows missing {key_label} + {value_column_label}: "
            f"baseline={baseline_skipped}, comparator={comparator_skipped}."
        )
    return rows, warnings


def compare_outcome_rows(baseline_rows: list[dict], comparator_rows: list[dict]) -> tuple[list[dict[str, object]], list[str]]:
    return _compare_metric_rows(baseline_rows, comparator_rows, ("Median",))


def compare_economics_rows(baseline_rows: list[dict], comparator_rows: list[dict]) -> tuple[list[dict[str, object]], list[str]]:
    return _compare_metric_rows(
        baseline_rows,
        comparator_rows,
        ("Value", "value"),
        ("Metric", "metric", "component"),
        require_value=False,
    )


def attributable_risk_compare_rows(payload: dict | None) -> list[dict[str, object]]:
    """Flatten one attributable-risk payload into export-friendly Compare rows."""
    if not isinstance(payload, dict):
        payload = {}

    payload_status = payload.get("status")
    rows: list[dict[str, object]] = [
        {
            "payloadField": "status",
            "rowIndex": None,
            "status": payload_status,
            "empty": False,
            "itemJson": None,
        }
    ]

    for field in ("calculatedRows", "missingInputs", "unsupportedMetrics", "messages"):
        items = payload.get(field)
        if not isinstance(items, list):
            items = []
        if not items:
            rows.append(
                {
                    "payloadField": field,
                    "rowIndex": None,
                    "status": payload_status,
                    "empty": True,
                    "itemJson": None,
                }
            )
            continue
        for index, item in enumerate(items):
            row = {
                "payloadField": field,
                "rowIndex": index,
                "status": payload_status,
                "empty": False,
                "itemJson": json.dumps(json_safe_for_export(item), sort_keys=True),
            }
            if isinstance(item, dict):
                row.update(json_safe_for_export(item))
            else:
                row["message"] = item
            rows.append(row)

    return rows


def economics_rows(economics: dict | None) -> list[dict]:
    if not economics:
        return []
    return economics.get("summaryRows") or []


def csv_download(rows: list[dict]) -> bytes:
    return pd.DataFrame(rows).to_csv(index=False).encode("utf-8")


def chart_rows(diff_rows: list[dict[str, object]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in diff_rows[:8]:
        base = coerce_number(row.get("baseline"))
        comp = coerce_number(row.get("comparator"))
        if base is None or comp is None:
            continue
        rows.append({"metric": row["metric"], "scenario": "Baseline", "value": base})
        rows.append({"metric": row["metric"], "scenario": "Comparator", "value": comp})
    return pd.DataFrame(rows)


st.title("APY Compare")
st.caption("Compare one baseline APY scenario with one comparator APY scenario.")

cols = st.columns(4)
if cols[0].button("Load default baseline", type="primary"):
    try:
        had_outputs = has_compare_outputs()
        st.session_state["compare_baseline_config"] = backend.default_config()
        reset_compare_outputs(outputs_cleared=had_outputs)
        st.session_state["compare_dirty"] = False
        sync_backend_status(backend.status())
        st.success("Default baseline loaded.")
    except Exception as exc:
        message = f"Could not load baseline defaults: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if cols[1].button("Use current scenario"):
    current = clone_config(st.session_state.get("config"))
    if current:
        had_outputs = has_compare_outputs()
        st.session_state["compare_baseline_config"] = current
        reset_compare_outputs(outputs_cleared=had_outputs)
        st.session_state["compare_dirty"] = False
        st.success("Current scenario copied to baseline.")
    else:
        st.info("No current APY scenario is loaded.")

if cols[2].button("Copy baseline to comparator"):
    baseline = clone_config(st.session_state.get("compare_baseline_config"))
    if baseline:
        baseline = populate_editable_compare_defaults(baseline)
        baseline["scenarioLabel"] = f"{baseline.get('scenarioLabel', 'APY scenario')} comparator"
        st.session_state["compare_comparator_config"] = baseline
        mark_compare_dirty()
        st.success("Baseline copied to comparator.")
    else:
        st.info("Load a baseline first.")

if cols[3].button("Load default comparator"):
    try:
        st.session_state["compare_comparator_config"] = backend.default_config()
        mark_compare_dirty()
        sync_backend_status(backend.status())
        st.success("Default comparator loaded.")
    except Exception as exc:
        message = f"Could not load comparator defaults: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

baseline_config = st.session_state.get("compare_baseline_config")
comparator_config = st.session_state.get("compare_comparator_config")
baseline_bundle = st.session_state.get("compare_baseline_bundle")
comparator_bundle = st.session_state.get("compare_comparator_bundle")

if not baseline_config or not comparator_config:
    st.info("Load a baseline and comparator to begin.")
    st.stop()

st.subheader("Named Stress Test")
preset_names = list(STRESS_TEST_PRESETS)
selected_preset = st.selectbox(
    "Preset",
    preset_names,
    index=choice_index(st.session_state.get("compare_selected_preset"), preset_names),
)
st.session_state["compare_selected_preset"] = selected_preset
st.caption(STRESS_TEST_PRESETS[selected_preset]["note"])
if st.button("Apply preset to baseline/comparator"):
    preset = STRESS_TEST_PRESETS[selected_preset]
    baseline = apply_overrides(baseline_config, preset["baseline"])
    comparator = apply_overrides(comparator_config, preset["comparator"])
    baseline["scenarioLabel"] = f"{selected_preset} baseline"
    comparator["scenarioLabel"] = f"{selected_preset} comparator"
    st.session_state["compare_baseline_config"] = baseline
    st.session_state["compare_comparator_config"] = comparator
    st.session_state["compare_selected_preset"] = selected_preset
    mark_compare_dirty()
    baseline_config = baseline
    comparator_config = comparator
    baseline_bundle = st.session_state.get("compare_baseline_bundle")
    comparator_bundle = st.session_state.get("compare_comparator_bundle")
    st.success("Stress-test preset applied.")

if st.session_state.get("compare_results_stale"):
    st.warning("Comparison results are stale because assumptions changed. Outputs are hidden until comparison is rerun.")
if st.session_state.get("compare_outputs_cleared"):
    st.warning("Previous compare outputs were cleared because assumptions changed. Rerun comparison to generate current outputs.")

st.subheader("Comparator Assumptions")
with st.form("compare_comparator_edits"):
    edited = clone_config(comparator_config) or {}
    edited["testType"] = st.selectbox(
        "Test type",
        TEST_OPTIONS,
        index=choice_index(edited.get("testType"), TEST_OPTIONS),
    )
    edited["regimen"] = st.selectbox(
        "Regimen",
        REGIMEN_OPTIONS,
        index=choice_index(edited.get("regimen"), REGIMEN_OPTIONS),
    )
    edited["screeningStrategy"] = st.selectbox(
        "Screening strategy",
        STRATEGY_OPTIONS,
        index=choice_index(edited.get("screeningStrategy"), STRATEGY_OPTIONS, 3),
    )
    edited["screenCoverage"] = st.number_input(
        "Screen coverage",
        min_value=0.0,
        max_value=1.0,
        value=float(edited.get("screenCoverage", 0.0)),
        step=0.05,
    )
    p_start = st.text_input(
        "Treatment start probability",
        value=optional_text(editable_compare_value(edited, "pStartTPT")),
    )
    p_complete = st.text_input(
        "Regimen completion probability",
        value=optional_text(editable_compare_value(edited, "regimenPComplete")),
    )
    p_adr = st.text_input(
        "ADR stop probability",
        value=optional_text(editable_compare_value(edited, "regimenADRstop")),
    )
    p_eff = st.text_input(
        "Full-course efficacy",
        value=optional_text(editable_compare_value(edited, "regimenEffFull")),
    )
    submitted = st.form_submit_button("Apply comparator edits")

if submitted:
    try:
        edited["pStartTPT"] = parse_optional_float(p_start)
        edited["regimenPComplete"] = parse_optional_float(p_complete)
        edited["regimenADRstop"] = parse_optional_float(p_adr)
        edited["regimenEffFull"] = parse_optional_float(p_eff)
        if edited != comparator_config:
            st.session_state["compare_comparator_config"] = edited
            mark_compare_dirty()
            comparator_config = edited
            st.success("Comparator edits applied.")
        else:
            st.info("No comparator fields changed.")
    except ValueError as exc:
        st.error(f"Invalid comparator number: {exc}")

run_cols = st.columns(3)
if run_cols[0].button("Validate both"):
    try:
        st.session_state["compare_baseline_validation_report"] = backend.validate_config(baseline_config)
        st.session_state["compare_comparator_validation_report"] = backend.validate_config(comparator_config)
        sync_backend_status(backend.status())
        st.success("Both scenarios validated.")
    except Exception as exc:
        message = f"Validation failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if run_cols[1].button("Run comparison", type="primary"):
    try:
        baseline_run_bundle = backend.run_scenario_bundle(baseline_config)
        comparator_run_bundle = backend.run_scenario_bundle(comparator_config)
        st.session_state["compare_baseline_bundle"] = baseline_run_bundle
        st.session_state["compare_comparator_bundle"] = comparator_run_bundle
        store_attributable_risk_compare_payloads(backend, baseline_run_bundle, comparator_run_bundle)
        st.session_state["compare_economics_config"] = clone_config(st.session_state.get("economics_config"))
        st.session_state["compare_baseline_economics_results"] = None
        st.session_state["compare_comparator_economics_results"] = None
        st.session_state["compare_dirty"] = False
        st.session_state["compare_results_stale"] = False
        st.session_state["compare_economics_stale"] = False
        st.session_state["compare_outputs_cleared"] = False
        st.session_state["compare_last_run_at"] = utc_now()
        sync_backend_status(backend.status())
        st.success("Comparison run completed.")
    except Exception as exc:
        message = f"Comparison run failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

compare_economics_config = st.session_state.get("compare_economics_config")
can_run_econ = bool(
    compare_economics_config
    and st.session_state.get("compare_baseline_bundle")
    and st.session_state.get("compare_comparator_bundle")
    and not st.session_state.get("compare_results_stale")
)
if run_cols[2].button("Run economics comparison", disabled=not can_run_econ):
    try:
        st.session_state["compare_baseline_economics_results"] = backend.run_economics_for_config(
            baseline_config,
            compare_economics_config,
        )
        st.session_state["compare_comparator_economics_results"] = backend.run_economics_for_config(
            comparator_config,
            compare_economics_config,
        )
        st.session_state["compare_economics_stale"] = False
        st.session_state["compare_last_economics_run_at"] = utc_now()
        sync_backend_status(backend.status())
        st.success("Economics comparison completed.")
    except Exception as exc:
        message = f"Economics comparison failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

baseline_bundle = st.session_state.get("compare_baseline_bundle")
comparator_bundle = st.session_state.get("compare_comparator_bundle")
assumptions_rows = assumptions_diff_rows(
    baseline_config,
    comparator_config,
    baseline_bundle,
    comparator_bundle,
)

st.subheader("Selected assumptions diff")
st.dataframe(arrow_safe_dataframe(assumptions_rows), width="stretch", hide_index=True)

outcome_diff, outcome_warnings = compare_outcome_rows(bundle_metric_rows(baseline_bundle), bundle_metric_rows(comparator_bundle))

st.subheader("Outcomes Diff")
if st.session_state.get("compare_results_stale"):
    st.info("Rerun comparison to show current outcome differences.")
elif outcome_diff:
    for warning in outcome_warnings:
        st.warning(warning)
    st.dataframe(arrow_safe_dataframe(outcome_diff), width="stretch", hide_index=True)
    chart_data = chart_rows(outcome_diff)
    if not chart_data.empty:
        chart = (
            alt.Chart(chart_data)
            .mark_bar()
            .encode(
                x=alt.X("metric:N", sort=None, title="Metric"),
                xOffset=alt.XOffset("scenario:N"),
                y=alt.Y("value:Q", title="Value"),
                color=alt.Color("scenario:N", title="Scenario"),
                tooltip=["metric:N", "scenario:N", "value:Q"],
            )
        )
        st.altair_chart(chart, width="stretch")
else:
    st.info("Run comparison to show outcome differences.")

baseline_attributable_risk = st.session_state.get("compare_baseline_attributable_risk_payload")
comparator_attributable_risk = st.session_state.get("compare_comparator_attributable_risk_payload")

st.subheader("Attributable Risk (partial)")
st.caption("Partial Compare support: unsupported or missing-input payloads are shown without blocking the comparison.")
if st.session_state.get("compare_results_stale"):
    st.info("Rerun comparison to show current attributable-risk payloads.")
elif baseline_attributable_risk or comparator_attributable_risk:
    status_rows = []
    detail_rows = []
    for label, payload in (
        ("Baseline", baseline_attributable_risk),
        ("Comparator", comparator_attributable_risk),
    ):
        if not isinstance(payload, dict):
            continue
        messages = payload.get("messages")
        if not isinstance(messages, list):
            messages = []
        status_rows.append(
            {
                "scenario": label,
                "status": payload.get("status"),
                "messages": " | ".join(str(message) for message in messages),
            }
        )
        detail_rows.extend(
            {"scenario": label, **row}
            for row in attributable_risk_compare_rows(payload)
        )

    st.dataframe(arrow_safe_dataframe(status_rows), width="stretch", hide_index=True)
    with st.expander("Flattened attributable-risk payload rows", expanded=False):
        st.dataframe(arrow_safe_dataframe(detail_rows), width="stretch", hide_index=True)
else:
    st.info("Run comparison to show attributable-risk payload status.")

base_econ = st.session_state.get("compare_baseline_economics_results")
comp_econ = st.session_state.get("compare_comparator_economics_results")
econ_diff, econ_warnings = compare_economics_rows(economics_rows(base_econ), economics_rows(comp_econ))

st.subheader("Economics Diff")
if st.session_state.get("compare_economics_stale") and (base_econ or comp_econ):
    st.warning("Economics comparison is stale. Rerun economics comparison.")
if st.session_state.get("compare_results_stale"):
    st.info("Rerun comparison before showing economics differences.")
elif st.session_state.get("compare_economics_stale") and (base_econ or comp_econ):
    st.info("Stale economics differences are hidden until economics comparison is rerun.")
elif econ_diff:
    for warning in econ_warnings:
        st.warning(warning)
    st.dataframe(arrow_safe_dataframe(econ_diff), width="stretch", hide_index=True)
elif compare_economics_config:
    st.info("Run economics comparison to show cost differences.")
else:
    st.info("Load economics assumptions on the Economics page, then rerun comparison to snapshot them for Compare.")

st.subheader("Manual Test Status")
status_options = ["Not reviewed", "Pass", "Concern", "Fail"]
manual_status = st.selectbox(
    "Review status",
    status_options,
    index=choice_index(st.session_state.get("compare_manual_status"), status_options),
)
st.session_state["compare_manual_status"] = manual_status
manual_notes = st.text_area(
    "Review notes",
    value=st.session_state.get("compare_manual_notes", ""),
    placeholder="Record expected-direction checks, anomalies, or follow-up actions.",
)
st.session_state["compare_manual_notes"] = manual_notes

st.subheader("Downloads")
comparison_payload = preset_payload(
    assumptions_rows,
    outcome_diff,
    outcome_warnings,
    econ_diff,
    econ_warnings,
)
download_cols = st.columns(4)
download_cols[0].download_button(
    "Assumptions diff CSV",
    data=csv_download(assumptions_rows),
    file_name=f"{safe_download_stem('apy_compare', 'assumptions_diff')}.csv",
    mime="text/csv",
)
if outcome_diff and not st.session_state.get("compare_results_stale"):
    download_cols[1].download_button(
        "Outcomes diff CSV",
        data=csv_download(outcome_diff),
        file_name=f"{safe_download_stem('apy_compare', 'outcomes_diff')}.csv",
        mime="text/csv",
    )
if econ_diff and not st.session_state.get("compare_results_stale") and not st.session_state.get("compare_economics_stale"):
    download_cols[2].download_button(
        "Economics diff CSV",
        data=csv_download(econ_diff),
        file_name=f"{safe_download_stem('apy_compare', 'economics_diff')}.csv",
        mime="text/csv",
    )
download_cols[3].download_button(
    "Comparison JSON",
    data=json_download(comparison_payload),
    file_name=f"{safe_download_stem('apy_compare', 'comparison')}.json",
    mime="application/json",
)
