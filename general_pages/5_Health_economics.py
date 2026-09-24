from __future__ import annotations

from typing import Any

import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.state import ECONOMICS_KEY, RESULTS_KEY, STALE_KEY, init_general_state, page_link
from app.general.terminology import display_text


DEMONSTRATION_COST_SOURCE = "Demonstration local-pathway working assumption"
SUMMARY_METRICS = {
    "comparatorCost": "Total cost without screening",
    "interventionCost": "Total cost with screening",
    "incrementalCost": "Incremental cost",
    "comparatorDALYs": "DALYs without screening",
    "interventionDALYs": "DALYs with screening",
    "dalysAverted": "DALYs averted",
    "activeTBCasesPrevented": "Active TB cases prevented",
}


def _economics_config() -> dict[str, Any]:
    from engine.apy.working_defaults import build_unified_working_default_preset

    return build_unified_working_default_preset()["economicsConfig"]


def _source_label(citation: Any) -> str:
    text = str(citation or "")
    if text.startswith("Dale KD"):
        return "Dale et al. 2022, Am J Epidemiol (Australian inputs)"
    return display_text(text, fallback=DEMONSTRATION_COST_SOURCE) or "Not recorded"


def _fmt(value: Any, digits: int = 0) -> str:
    try:
        return f"{float(value):,.{digits}f}"
    except (TypeError, ValueError):
        return ""


init_general_state()
st.title("Health economics")
st.warning(
    "Cost and DALY inputs are demonstration working defaults taken from a published Australian analysis "
    "(Australian dollars, 2019 prices, health-system perspective). Replace them with local costs before "
    "interpreting economic results."
)

bundle = st.session_state.get(RESULTS_KEY)
if not bundle or st.session_state.get(STALE_KEY):
    st.info("Run the analysis with the current inputs first; then return here to calculate health economics.")
    page_link("general_pages/3_Run_analysis.py", label="Open Run analysis")
    st.stop()

econ_config = _economics_config()
if st.button("Calculate health economics", type="primary"):
    from engine.apy.event_ledger_economics import run_event_ledger_health_economics

    try:
        st.session_state[ECONOMICS_KEY] = run_event_ledger_health_economics(bundle, econ_config)
    except Exception as exc:
        st.session_state["general_last_error"] = repr(exc)
        st.error("Health-economic results could not be calculated for this run.")

results = st.session_state.get(ECONOMICS_KEY)
if not results:
    st.caption("Health economics uses the outcomes of the current analysis; the epidemiological model is not rerun.")
    st.stop()

stochastic = (bundle.get("metadata") or {}).get("modelType") == "agent_based"
rows = []
for row in results.get("summaryRows") or []:
    if row.get("discountProfile") != "primary" or row.get("metric") not in SUMMARY_METRICS:
        continue
    if stochastic:
        rows.append(
            {
                "Measure": SUMMARY_METRICS[row["metric"]],
                "Mean": _fmt(row.get("mean"), 1),
                "Median": _fmt(row.get("median"), 1),
                "Low 95%": _fmt(row.get("p2_5"), 1),
                "High 95%": _fmt(row.get("p97_5"), 1),
            }
        )
    else:
        rows.append({"Measure": SUMMARY_METRICS[row["metric"]], "Expected value": _fmt(row.get("mean"), 1)})
icer = next(
    (
        row
        for row in results.get("summaryRows") or []
        if row.get("discountProfile") == "primary" and row.get("metric") == "primaryICER_ratioOfMeans"
    ),
    None,
)
st.subheader("Summary (3% annual discounting)")
if rows:
    st.dataframe(arrow_safe_dataframe(rows), use_container_width=True, hide_index=True)
    if stochastic:
        st.caption("Ranges describe variation across simulated populations; they are not confidence intervals.")
    else:
        st.caption("Deterministic preview: single expected values without simulation variation.")
if icer:
    classification = icer.get("classification")
    if classification == "dominant":
        st.markdown("**Cost per DALY averted:** screening is less costly and more effective than no screening (dominant).")
    elif icer.get("mean") is not None:
        st.markdown(f"**Cost per DALY averted:** {_fmt(icer.get('mean'))} (ratio of mean costs to mean DALYs averted).")
    else:
        st.markdown(f"**Cost per DALY averted:** {classification or 'not available'}.")
st.caption("No willingness-to-pay threshold is set, so no cost-effectiveness conclusion is drawn.")

with st.expander("Unit costs used"):
    cost_rows = [
        {
            "Item": item.get("description"),
            "Unit cost": _fmt(item.get("convertedTargetYearCost"), 2),
            "Currency and year": f"{item.get('targetCurrency', '')} {item.get('targetPriceYear', '')}".strip(),
            "Source": _source_label(item.get("sourceCitation")),
        }
        for item in results.get("costItems") or []
    ]
    st.dataframe(arrow_safe_dataframe(cost_rows), use_container_width=True, hide_index=True)

unresolved = [display_text(item.get("message", ""), fallback="An economic input is unresolved.") for item in results.get("unresolvedInputs") or []]
if unresolved:
    st.markdown("**Unresolved economic inputs**")
    for message in unresolved:
        st.write(f"- {message}")
page_link("general_pages/6_Evidence_and_technical_information.py", label="Evidence and technical information")
