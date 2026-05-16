from __future__ import annotations

from pathlib import Path
import sys

import altair as alt
import pandas as pd
import streamlit as st

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from app.display import arrow_safe_dataframe, safe_download_stem
from app.state import (
    get_backend_name,
    init_session_state,
    mark_dynamic_abm_compare_completed,
)
from engine.integration.compare_dynamic_abm_v9 import compare_dynamic_abm_v9, number_or_none


TOTAL_ALIGNED_METRICS = 6


def metric_lookup(bundle: dict) -> dict[str, object]:
    headline = bundle.get("headline") or {}
    rows = headline.get("keyMetricsRows") or headline.get("summaryRows") or []
    lookup: dict[str, object] = {}
    for row in rows if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            continue
        metric = row.get("Metric") or row.get("metric")
        if metric in (None, ""):
            continue
        lookup[str(metric)] = row.get("Value", row.get("Median", row.get("value")))
    return lookup


def headline_rows(bundle: dict) -> list[dict[str, object]]:
    headline = bundle.get("headline") or {}
    rows = headline.get("keyMetricsRows") or headline.get("summaryRows") or []
    return rows if isinstance(rows, list) else []


def bundle_debug_rows(bundle: dict) -> list[dict[str, object]]:
    metadata = bundle.get("metadata") or {}
    technical = bundle.get("technical") or {}
    dynamic_comparison = technical.get("dynamicComparison") if isinstance(technical, dict) else {}
    dynamic_comparison = dynamic_comparison if isinstance(dynamic_comparison, dict) else {}
    projection = bundle.get("projection") or {}
    annual_rows = projection.get("annualRows") or []
    metric_names = [str(row.get("Metric") or row.get("metric")) for row in headline_rows(bundle) if isinstance(row, dict)]
    dynamic_comparison_rows = dynamic_comparison.get("metricRows") or []
    dynamic_comparison_metrics = [
        str(row.get("Metric") or row.get("metric"))
        for row in dynamic_comparison_rows
        if isinstance(row, dict)
    ]
    return [
        {"field": "topLevelKeys", "value": sorted(bundle.keys())},
        {"field": "model", "value": bundle.get("model", "")},
        {"field": "modelVersion", "value": bundle.get("modelVersion") or metadata.get("modelVersion", "")},
        {"field": "contractVersion", "value": bundle.get("contractVersion") or metadata.get("contractVersion", "")},
        {"field": "metadata.modelVersion", "value": metadata.get("modelVersion", "")},
        {"field": "headlineRowCount", "value": len(headline_rows(bundle))},
        {"field": "keyMetrics", "value": metric_names},
        {"field": "projectionAnnualRowsCount", "value": len(annual_rows) if isinstance(annual_rows, list) else 0},
        {"field": "technical.dynamicComparison.available", "value": dynamic_comparison.get("available", "")},
        {"field": "technical.dynamicComparison.metrics", "value": dynamic_comparison_metrics},
    ]


def summary_rows(title: str, bundle: dict) -> list[dict[str, object]]:
    lookup = metric_lookup(bundle)
    metadata = bundle.get("metadata") or {}
    return [
        {"section": title, "field": "model", "value": bundle.get("model") or metadata.get("modelVersion", "")},
        {"section": title, "field": "modelVersion", "value": bundle.get("modelVersion") or metadata.get("modelVersion", "")},
        {"section": title, "field": "scenarioLabel", "value": metadata.get("scenarioLabel", "")},
        {"section": title, "field": "population", "value": metadata.get("population", lookup.get("population", ""))},
        {"section": title, "field": "horizon", "value": lookup.get("projection_horizon", metadata.get("time_horizon", ""))},
    ]


def numeric_comparison_rows(rows: list[dict[str, object]]) -> pd.DataFrame:
    numeric_rows: list[dict[str, object]] = []
    for row in rows:
        if not row.get("comparable"):
            continue
        for source, value_key in (("Dynamic", "dynamic_value"), ("APY ABM", "abm_value")):
            value = row.get(value_key)
            try:
                numeric_rows.append({"metric": row["metric"], "model": source, "value": float(value)})
            except (TypeError, ValueError):
                continue
    return pd.DataFrame(numeric_rows)


def nested_dict(bundle: dict, *keys: str) -> dict:
    value: object = bundle
    for key in keys:
        if not isinstance(value, dict):
            return {}
        value = value.get(key)
    return value if isinstance(value, dict) else {}


def match_status(dynamic_value: object, abm_value: object) -> str:
    dynamic_number = number_or_none(dynamic_value)
    abm_number = number_or_none(abm_value)
    if dynamic_number is None or abm_number is None:
        return "unknown"
    return "yes" if dynamic_number == abm_number else "no"


def comparison_summary_rows(
    dynamic_bundle: dict,
    abm_bundle: dict,
    comparison_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    comparable_count = sum(1 for row in comparison_rows if row.get("comparable"))
    dynamic_comparison = nested_dict(abm_bundle, "technical", "dynamicComparison")
    dynamic_metrics = {row["metric"]: row.get("dynamic_value") for row in comparison_rows}
    abm_metrics = {row["metric"]: row.get("abm_value") for row in comparison_rows}
    return [
        {"field": "APY backend", "value": get_backend_name()},
        {
            "field": "Dynamic bundle present",
            "value": "yes" if dynamic_bundle else "no",
        },
        {
            "field": "APY bundle present",
            "value": "yes" if abm_bundle else "no",
        },
        {
            "field": "Comparable metrics available",
            "value": f"{comparable_count} of {TOTAL_ALIGNED_METRICS}",
        },
        {
            "field": "APY dynamic-comparison available",
            "value": dynamic_comparison.get("available", "unknown"),
        },
        {
            "field": "APY dynamic-comparison source",
            "value": dynamic_comparison.get("source") or "unavailable",
        },
        {
            "field": "Population match",
            "value": match_status(dynamic_metrics.get("population"), abm_metrics.get("population")),
        },
        {
            "field": "Horizon match",
            "value": match_status(dynamic_metrics.get("horizon"), abm_metrics.get("horizon")),
        },
    ]


init_session_state()

st.title("Dynamic + APY ABM Compare")
st.caption(
    "Compares the latest Python dynamic model projection with the latest APY v9 ABM results bundle. "
    "Metrics are shown side by side where comparable; structural differences are flagged rather than hidden."
)

dynamic_bundle = st.session_state.get("dynamic_results_bundle")
abm_bundle = st.session_state.get("results_bundle")

if not dynamic_bundle:
    st.warning("No dynamic model results bundle is available. Open Dynamic Model, calibrate, then simulate.")
if not abm_bundle:
    st.warning("No APY ABM results bundle is available. Open APY Run Model and run the selected APY backend.")
if not dynamic_bundle or not abm_bundle:
    st.stop()

summary = summary_rows("Dynamic model", dynamic_bundle) + summary_rows("APY ABM", abm_bundle)
st.subheader("Latest Results")
st.dataframe(arrow_safe_dataframe(summary), width="stretch", hide_index=True)

with st.expander("Dynamic bundle debug", expanded=False):
    st.dataframe(arrow_safe_dataframe(bundle_debug_rows(dynamic_bundle)), width="stretch", hide_index=True)

with st.expander("APY bundle debug", expanded=False):
    st.dataframe(arrow_safe_dataframe(bundle_debug_rows(abm_bundle)), width="stretch", hide_index=True)

comparison = compare_dynamic_abm_v9(dynamic_bundle, abm_bundle)
rows = comparison.get("comparisonRows") or []
warnings = comparison.get("warnings") or []
mark_dynamic_abm_compare_completed(rows, warnings)

comparable_rows = [row for row in rows if row.get("comparable")]
non_comparable_rows = [row for row in rows if not row.get("comparable")]

st.subheader("Comparison Summary")
st.dataframe(
    arrow_safe_dataframe(comparison_summary_rows(dynamic_bundle, abm_bundle, rows)),
    width="content",
    hide_index=True,
)

st.subheader("Comparable Metrics")
if comparable_rows:
    st.dataframe(arrow_safe_dataframe(comparable_rows), width="stretch", hide_index=True)
    chart_data = numeric_comparison_rows(comparable_rows)
    if not chart_data.empty:
        chart = (
            alt.Chart(chart_data)
            .mark_bar()
            .encode(
                x=alt.X("metric:N", sort=None, title="Metric"),
                xOffset=alt.XOffset("model:N"),
                y=alt.Y("value:Q", title="Value"),
                color=alt.Color("model:N", title="Model"),
                tooltip=["metric:N", "model:N", "value:Q"],
            )
        )
        st.altair_chart(chart, width="stretch")
else:
    st.info("No numeric aligned metrics are currently available in both bundles.")

st.subheader("Warnings And Non-Comparable Metrics")
if non_comparable_rows:
    st.dataframe(arrow_safe_dataframe(non_comparable_rows), width="stretch", hide_index=True)
if comparison.get("missing_dynamic_metrics"):
    st.warning("Missing dynamic metrics: " + ", ".join(comparison["missing_dynamic_metrics"]))
if comparison.get("missing_abm_metrics"):
    st.warning("Missing APY ABM metrics: " + ", ".join(comparison["missing_abm_metrics"]))
if comparison.get("structurally_non_comparable_metrics"):
    with st.expander("Structurally non-comparable metrics", expanded=False):
        st.write(comparison["structurally_non_comparable_metrics"])
for warning in warnings:
    st.warning(warning)

st.subheader("Download")
st.download_button(
    "Comparison CSV",
    data=pd.DataFrame(rows).to_csv(index=False).encode("utf-8"),
    file_name=f"{safe_download_stem('dynamic_abm', 'comparison')}.csv",
    mime="text/csv",
)
