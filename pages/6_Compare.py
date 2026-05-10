from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone

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
    if st.session_state.get("compare_baseline_bundle") or st.session_state.get("compare_comparator_bundle"):
        st.session_state["compare_results_stale"] = True
    if st.session_state.get("compare_baseline_economics_results") or st.session_state.get("compare_comparator_economics_results"):
        st.session_state["compare_economics_stale"] = True


def reset_compare_outputs() -> None:
    st.session_state["compare_baseline_bundle"] = None
    st.session_state["compare_comparator_bundle"] = None
    st.session_state["compare_baseline_economics_results"] = None
    st.session_state["compare_comparator_economics_results"] = None
    st.session_state["compare_baseline_validation_report"] = None
    st.session_state["compare_comparator_validation_report"] = None
    st.session_state["compare_results_stale"] = False
    st.session_state["compare_economics_stale"] = False


def assumptions_diff_rows(baseline: dict, comparator: dict) -> list[dict[str, object]]:
    fields = ["scenarioLabel", "N", "nReps", "seed", *COMPARE_FIELDS]
    return [
        {
            "field": field,
            "baseline": baseline.get(field),
            "comparator": comparator.get(field),
            "changed": baseline.get(field) != comparator.get(field),
        }
        for field in fields
    ]


def metric_name(row: dict) -> str | None:
    for key in ("Metric", "metric", "Name", "name", "Outcome", "outcome"):
        if row.get(key) not in (None, ""):
            return str(row[key])
    return None


def metric_value(row: dict) -> object:
    for key in ("Median", "median", "Value", "value", "Estimate", "estimate", "Mean", "mean"):
        if key in row:
            return row[key]
    return None


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


def compare_rows(baseline_rows: list[dict], comparator_rows: list[dict]) -> list[dict[str, object]]:
    baseline = {metric_name(row): metric_value(row) for row in baseline_rows if metric_name(row)}
    comparator = {metric_name(row): metric_value(row) for row in comparator_rows if metric_name(row)}
    rows: list[dict[str, object]] = []
    for metric in sorted(set(baseline) | set(comparator)):
        base_value = baseline.get(metric)
        comp_value = comparator.get(metric)
        base_num = coerce_number(base_value)
        comp_num = coerce_number(comp_value)
        abs_diff = None
        pct_diff = None
        if base_num is not None and comp_num is not None:
            abs_diff = comp_num - base_num
            if base_num != 0:
                pct_diff = abs_diff / base_num
        rows.append(
            {
                "metric": metric,
                "baseline": base_value,
                "comparator": comp_value,
                "absoluteDifference": abs_diff,
                "percentDifference": pct_diff,
            }
        )
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
        rows.append({"metric": row["metric"], "baseline": base, "comparator": comp})
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).set_index("metric")


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

if not baseline_config or not comparator_config:
    st.info("Load a baseline and comparator to begin.")
    st.stop()

if st.session_state.get("compare_results_stale"):
    st.warning("Comparison results are stale because assumptions changed.")

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

st.subheader("Assumptions Diff")
assumptions_rows = assumptions_diff_rows(baseline_config, comparator_config)
st.dataframe(arrow_safe_dataframe(assumptions_rows), width="stretch", hide_index=True)

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
        st.session_state["compare_baseline_economics_results"] = None
        st.session_state["compare_comparator_economics_results"] = None
        st.session_state["compare_dirty"] = False
        st.session_state["compare_results_stale"] = False
        st.session_state["compare_economics_stale"] = True
        st.session_state["compare_last_run_at"] = utc_now()
        sync_backend_status(backend.status())
        st.success("Comparison run completed.")
    except Exception as exc:
        message = f"Comparison run failed: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

economics_config = st.session_state.get("economics_config")
can_run_econ = bool(
    economics_config
    and st.session_state.get("compare_baseline_bundle")
    and st.session_state.get("compare_comparator_bundle")
    and not st.session_state.get("compare_results_stale")
)
if run_cols[2].button("Run economics comparison", disabled=not can_run_econ):
    try:
        st.session_state["compare_baseline_economics_results"] = backend.run_economics_for_config(
            baseline_config,
            economics_config,
        )
        st.session_state["compare_comparator_economics_results"] = backend.run_economics_for_config(
            comparator_config,
            economics_config,
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
outcome_diff = compare_rows(bundle_metric_rows(baseline_bundle), bundle_metric_rows(comparator_bundle))

st.subheader("Outcomes Diff")
if outcome_diff:
    st.dataframe(arrow_safe_dataframe(outcome_diff), width="stretch", hide_index=True)
    chart_data = chart_rows(outcome_diff)
    if not chart_data.empty:
        st.bar_chart(chart_data)
else:
    st.info("Run comparison to show outcome differences.")

base_econ = st.session_state.get("compare_baseline_economics_results")
comp_econ = st.session_state.get("compare_comparator_economics_results")
econ_diff = compare_rows(economics_rows(base_econ), economics_rows(comp_econ))

st.subheader("Economics Diff")
if st.session_state.get("compare_economics_stale") and (base_econ or comp_econ):
    st.warning("Economics comparison is stale. Rerun economics comparison.")
if econ_diff:
    st.dataframe(arrow_safe_dataframe(econ_diff), width="stretch", hide_index=True)
elif economics_config:
    st.info("Run economics comparison to show cost differences.")
else:
    st.info("Load economics assumptions on the Economics page to enable economics comparison.")

st.subheader("Downloads")
download_cols = st.columns(3)
download_cols[0].download_button(
    "Assumptions diff CSV",
    data=csv_download(assumptions_rows),
    file_name=f"{safe_download_stem('apy_compare', 'assumptions_diff')}.csv",
    mime="text/csv",
)
if outcome_diff:
    download_cols[1].download_button(
        "Outcomes diff CSV",
        data=csv_download(outcome_diff),
        file_name=f"{safe_download_stem('apy_compare', 'outcomes_diff')}.csv",
        mime="text/csv",
    )
if econ_diff:
    download_cols[2].download_button(
        "Economics diff CSV",
        data=csv_download(econ_diff),
        file_name=f"{safe_download_stem('apy_compare', 'economics_diff')}.csv",
        mime="text/csv",
    )
