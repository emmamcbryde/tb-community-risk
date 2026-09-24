from __future__ import annotations

from typing import Any

import altair as alt
import pandas as pd
import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.economics import (
    apply_cost_edits,
    changed_cost_edits,
    cost_rows,
    default_economics_config,
    icer_cloud_points,
    run_economics,
)
from app.general.state import (
    ECONOMICS_CONFIG_KEY,
    ECONOMICS_KEY,
    EDITOR_VERSION_KEY,
    RESULTS_KEY,
    economics_status,
    init_general_state,
    page_link,
    results_status,
    store_economics,
)
from app.general.terminology import USER_DEFINED_MARK, display_text


SUMMARY_METRICS = {
    "comparatorCost": "Total cost without screening",
    "interventionCost": "Total cost with screening",
    "incrementalCost": "Incremental cost",
    "comparatorDALYs": "DALYs without screening",
    "interventionDALYs": "DALYs with screening",
    "dalysAverted": "DALYs averted",
    "activeTBCasesPrevented": "Active TB cases prevented",
}


def _fmt(value: Any, digits: int = 0) -> str:
    try:
        return f"{float(value):,.{digits}f}"
    except (TypeError, ValueError):
        return "N/A"


init_general_state()
st.title("Health economics")
st.warning(
    "Cost and DALY inputs are demonstration working defaults taken from a published Australian analysis "
    "(Australian dollars, 2019 prices, health-system perspective). Replace them with local costs before "
    "interpreting economic results."
)

bundle = st.session_state.get(RESULTS_KEY)
if not bundle or results_status() != "current":
    st.info("Run the analysis with the current inputs first; then return here to calculate health economics.")
    page_link("general_pages/3_Run_analysis.py", label="Open Run analysis")
    st.stop()
stochastic = (bundle.get("metadata") or {}).get("modelType") == "agent_based"

st.subheader("Unit costs")
st.caption(
    "Edit a unit cost to explore cost assumptions. Cost changes recalculate economics from the completed "
    "analysis; the epidemiological model is not rerun and health outcomes do not change."
)
config = st.session_state[ECONOMICS_CONFIG_KEY]
rows = cost_rows(config)
edited = st.data_editor(
    pd.DataFrame(rows),
    key=f"general_cost_editor_{st.session_state[EDITOR_VERSION_KEY]}",
    hide_index=True,
    use_container_width=True,
    column_order=["Item", "Unit cost", "Currency and year", "Source"],
    disabled=["Item", "Currency and year", "Source"],
    column_config={"Unit cost": st.column_config.NumberColumn("Unit cost", min_value=0.0, format="%.2f")},
)
edits = changed_cost_edits(rows, edited.to_dict(orient="records"))
if edits:
    st.session_state[ECONOMICS_CONFIG_KEY] = apply_cost_edits(config, edits)
    st.rerun()
if any(row["Source"] == "User-defined" for row in rows):
    st.caption(f"{USER_DEFINED_MARK} unit costs are marked 'User-defined' in the Source column.")
    if st.button("Restore demonstration unit costs"):
        st.session_state[ECONOMICS_CONFIG_KEY] = default_economics_config()
        st.session_state[EDITOR_VERSION_KEY] = int(st.session_state[EDITOR_VERSION_KEY]) + 1
        st.rerun()

if stochastic:
    reps = int((bundle.get("metadata") or {}).get("nReps") or 0)
    st.caption(f"Calculating economics for {reps:,} simulated populations takes about {max(0.18 * reps / 60, 0.1):.0f} minute(s).")
status = economics_status()
if status == "stale":
    st.warning("Cost assumptions changed since the last calculation. Recalculate to update the economic results.")
label = "Recalculate health economics" if status != "none" else "Calculate health economics"
if st.button(label, type="primary", disabled=status == "current"):
    try:
        store_economics(run_economics(bundle, st.session_state[ECONOMICS_CONFIG_KEY]))
        st.rerun()
    except Exception as exc:
        st.session_state["general_last_error"] = repr(exc)
        st.error("Health-economic results could not be calculated for this run.")

results = st.session_state.get(ECONOMICS_KEY)
if not results:
    st.stop()

summary_rows = []
for row in results.get("summaryRows") or []:
    if row.get("discountProfile") != "primary" or row.get("metric") not in SUMMARY_METRICS:
        continue
    if stochastic:
        summary_rows.append(
            {
                "Measure": SUMMARY_METRICS[row["metric"]],
                "Mean": _fmt(row.get("mean"), 1),
                "Median": _fmt(row.get("median"), 1),
                "Low 95% (simulation)": _fmt(row.get("p2_5"), 1),
                "High 95% (simulation)": _fmt(row.get("p97_5"), 1),
            }
        )
    else:
        summary_rows.append(
            {
                "Measure": SUMMARY_METRICS[row["metric"]],
                "Expected value": _fmt(row.get("mean"), 1),
                "Low 95%": "N/A",
                "High 95%": "N/A",
            }
        )
st.subheader("Summary (3% annual discounting)")
if summary_rows:
    st.dataframe(arrow_safe_dataframe(summary_rows), use_container_width=True, hide_index=True)
    if stochastic:
        st.caption("Simulation intervals describe variation across simulated populations; they are not confidence intervals.")
    else:
        st.caption("Deterministic preview: single expected values; simulation intervals are not applicable.")

icer = next(
    (row for row in results.get("summaryRows") or [] if row.get("discountProfile") == "primary" and row.get("metric") == "primaryICER_ratioOfMeans"),
    None,
)
if icer:
    classification = icer.get("classification")
    if classification == "dominant":
        st.markdown("**Cost per DALY averted:** screening is less costly and more effective than no screening (dominant).")
    elif icer.get("mean") is not None:
        st.markdown(f"**Cost per DALY averted:** {_fmt(icer.get('mean'))} (ratio of mean incremental cost to mean DALYs averted).")
    else:
        st.markdown(f"**Cost per DALY averted:** {classification or 'not available'}.")
st.caption("No willingness-to-pay threshold is set, so no cost-effectiveness conclusion is drawn.")

points = icer_cloud_points(results)
if stochastic and points:
    st.subheader("Cost-effectiveness plane")
    frame = pd.DataFrame(points)
    origin = pd.DataFrame([{"dalysAverted": 0.0, "incrementalCost": 0.0, "label": "No screening (comparator)"}])
    cloud = alt.Chart(frame).mark_circle(size=28, opacity=0.5, color="#1f5fa8").encode(
        x=alt.X("dalysAverted:Q", title="DALYs averted vs no screening"),
        y=alt.Y("incrementalCost:Q", title="Incremental cost vs no screening"),
        tooltip=["replicate:Q", alt.Tooltip("dalysAverted:Q", format=".1f"), alt.Tooltip("incrementalCost:Q", format=",.0f")],
    )
    comparator = alt.Chart(origin).mark_point(shape="cross", size=160, color="#111827", filled=True).encode(
        x="dalysAverted:Q", y="incrementalCost:Q", tooltip=["label:N"]
    )
    axes = alt.Chart(pd.DataFrame({"zero": [0.0]}))
    st.altair_chart(
        (axes.mark_rule(color="#9ca3af").encode(y="zero:Q") + axes.mark_rule(color="#9ca3af").encode(x="zero:Q") + cloud + comparator).properties(height=340),
        use_container_width=True,
    )
    st.caption(
        "Each point pairs the health outcome and cost of one simulated population with the same population "
        "without screening. The comparator is at the origin. Cost-only changes move points vertically; "
        "health outcomes are unchanged."
    )

unresolved = [display_text(item.get("message", ""), fallback="An economic input is unresolved.") for item in results.get("unresolvedInputs") or []]
if unresolved:
    st.markdown("**Unresolved economic inputs**")
    for message in unresolved:
        st.write(f"- {message}")
page_link("general_pages/6_Evidence_and_technical_information.py", label="Evidence, downloads and technical information")
