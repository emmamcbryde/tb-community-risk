from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone

import altair as alt
import pandas as pd
import streamlit as st

from app.display import arrow_safe_dataframe, safe_download_stem
from app.state import get_backend, init_session_state, record_message, sync_backend_status


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


def mark_compare_dirty() -> None:
    st.session_state["compare_dirty"] = True
    if (
        st.session_state.get("compare_baseline_bundle")
        or st.session_state.get("compare_comparator_bundle")
        or st.session_state.get("compare_baseline_economics_results")
        or st.session_state.get("compare_comparator_economics_results")
    ):
        reset_compare_outputs()


def reset_compare_outputs() -> None:
    st.session_state["compare_baseline_bundle"] = None
    st.session_state["compare_comparator_bundle"] = None
    st.session_state["compare_economics_config"] = None
    st.session_state["compare_baseline_economics_results"] = None
    st.session_state["compare_comparator_economics_results"] = None
    st.session_state["compare_baseline_validation_report"] = None
    st.session_state["compare_comparator_validation_report"] = None
    st.session_state["compare_results_stale"] = False
    st.session_state["compare_economics_stale"] = False


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


def is_blank_default_value(value: object) -> bool:
    return value is None or value == "" or value == []


def resolved_assumption_value(field: str, config: dict, bundle: dict | None) -> object:
    interface_config = bundle_interface_config(bundle)
    if field in interface_config:
        return interface_config.get(field)
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


def _index_metric_rows(rows: list[dict], value_column: str) -> tuple[dict[str, object], list[str], list[str], int]:
    values: dict[str, object] = {}
    order: list[str] = []
    duplicates: list[str] = []
    skipped = 0
    for row in rows:
        if "Metric" not in row or value_column not in row:
            skipped += 1
            continue
        metric = row.get("Metric")
        if metric in (None, ""):
            skipped += 1
            continue
        metric_key = str(metric)
        if metric_key not in values:
            order.append(metric_key)
            values[metric_key] = row.get(value_column)
        else:
            duplicates.append(metric_key)
    return values, order, duplicates, skipped


def _compare_metric_rows(
    baseline_rows: list[dict],
    comparator_rows: list[dict],
    value_column: str,
) -> tuple[list[dict[str, object]], list[str]]:
    baseline, baseline_order, baseline_duplicates, baseline_skipped = _index_metric_rows(baseline_rows, value_column)
    comparator, comparator_order, comparator_duplicates, comparator_skipped = _index_metric_rows(comparator_rows, value_column)
    ordered_metrics = baseline_order + [metric for metric in comparator_order if metric not in baseline]
    rows: list[dict[str, object]] = []
    for metric in ordered_metrics:
        base_value = baseline.get(metric)
        comp_value = comparator.get(metric)
        base_num = coerce_number(base_value)
        comp_num = coerce_number(comp_value)
        abs_diff = None
        rel_diff = None
        if base_num is not None and comp_num is not None:
            abs_diff = comp_num - base_num
            if base_num != 0:
                rel_diff = abs_diff / base_num
        rows.append(
            {
                "metric": metric,
                "baseline": base_value,
                "comparator": comp_value,
                "absoluteDifference": abs_diff,
                "relativeDifference": rel_diff,
            }
        )
    warnings: list[str] = []
    duplicate_metrics = sorted(set(baseline_duplicates + comparator_duplicates))
    if duplicate_metrics:
        warnings.append(
            "Duplicate Metric values detected; the first row was used for: "
            + ", ".join(duplicate_metrics)
        )
    if baseline_skipped or comparator_skipped:
        warnings.append(
            f"Skipped rows missing Metric + {value_column}: "
            f"baseline={baseline_skipped}, comparator={comparator_skipped}."
        )
    return rows, warnings


def compare_outcome_rows(baseline_rows: list[dict], comparator_rows: list[dict]) -> tuple[list[dict[str, object]], list[str]]:
    return _compare_metric_rows(baseline_rows, comparator_rows, "Median")


def compare_economics_rows(baseline_rows: list[dict], comparator_rows: list[dict]) -> tuple[list[dict[str, object]], list[str]]:
    return _compare_metric_rows(baseline_rows, comparator_rows, "Value")


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
        st.session_state["compare_baseline_config"] = backend.default_config()
        reset_compare_outputs()
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
        st.session_state["compare_baseline_config"] = current
        reset_compare_outputs()
        st.session_state["compare_dirty"] = False
        st.success("Current scenario copied to baseline.")
    else:
        st.info("No current APY scenario is loaded.")

if cols[2].button("Copy baseline to comparator"):
    baseline = clone_config(st.session_state.get("compare_baseline_config"))
    if baseline:
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

if st.session_state.get("compare_results_stale"):
    st.warning("Comparison results are stale because assumptions changed. Outputs are hidden until comparison is rerun.")

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
    p_start = st.text_input("Treatment start probability", value=optional_text(edited.get("pStartTPT")))
    p_complete = st.text_input("Regimen completion probability", value=optional_text(edited.get("regimenPComplete")))
    p_adr = st.text_input("ADR stop probability", value=optional_text(edited.get("regimenADRstop")))
    p_eff = st.text_input("Full-course efficacy", value=optional_text(edited.get("regimenEffFull")))
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
        st.session_state["compare_baseline_bundle"] = backend.run_scenario_bundle(baseline_config)
        st.session_state["compare_comparator_bundle"] = backend.run_scenario_bundle(comparator_config)
        st.session_state["compare_economics_config"] = clone_config(st.session_state.get("economics_config"))
        st.session_state["compare_baseline_economics_results"] = None
        st.session_state["compare_comparator_economics_results"] = None
        st.session_state["compare_dirty"] = False
        st.session_state["compare_results_stale"] = False
        st.session_state["compare_economics_stale"] = False
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

st.subheader("Downloads")
download_cols = st.columns(3)
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
