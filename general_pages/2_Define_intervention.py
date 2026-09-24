from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.state import get_intervention, init_general_state, page_link, set_intervention
from app.general.terminology import USER_DEFINED_MARK
from engine.profiles.engine_mapping import default_intervention


TEST_OPTIONS = {"IGRA": "Interferon-gamma release assay (IGRA)", "TST": "Tuberculin skin test (TST)"}
REGIMEN_OPTIONS = {
    "3HP": "3HP - 3 months weekly isoniazid and rifapentine",
    "4R": "4R - 4 months daily rifampicin",
    "3HR": "3HR - 3 months daily isoniazid and rifampicin",
    "6H": "6H - 6 months daily isoniazid",
    "9H": "9H - 9 months daily isoniazid",
}
TARGETING_OPTIONS = {
    "prevent": "Prioritise people most likely to avoid active TB",
    "ltbi": "Prioritise people most likely to have LTBI",
    "cure": "Prioritise people most likely to complete effective treatment",
    "random": "No risk-based prioritisation",
}
FIELD_LABELS = {
    "testType": "Screening test",
    "regimen": "Preventive-treatment regimen",
    "screenCoverage": "Screening coverage",
    "screeningStrategy": "Targeting",
    "screeningWindowYears": "Screening period",
}


def _select(label: str, options: dict[str, str], current: str, key: str, help_text: str = "") -> str:
    codes = list(options)
    if st.session_state.get(key) not in codes:
        st.session_state[key] = current if current in codes else codes[0]
    return st.selectbox(label, codes, format_func=options.get, key=key, help=help_text or None)


init_general_state()
defaults = default_intervention()
intervention = get_intervention()

st.title("Define intervention")
st.write("Choose the screening test, preventive treatment and delivery settings to compare against no screening.")

updated = dict(intervention)
updated["testType"] = _select("Screening test", TEST_OPTIONS, intervention["testType"], "general_test_type")
updated["regimen"] = _select("Preventive-treatment regimen", REGIMEN_OPTIONS, intervention["regimen"], "general_regimen")
if st.session_state.get("general_coverage") is None:
    st.session_state["general_coverage"] = int(round(float(intervention["screenCoverage"]) * 100))
coverage = st.slider("Screening coverage (% of the population screened)", 0, 100, key="general_coverage")
updated["screenCoverage"] = coverage / 100.0
updated["screeningStrategy"] = _select(
    "Targeting",
    TARGETING_OPTIONS,
    intervention["screeningStrategy"],
    "general_strategy",
    "How people are prioritised when coverage is below 100%.",
)
if st.session_state.get("general_window") is None:
    st.session_state["general_window"] = int(intervention["screeningWindowYears"])
updated["screeningWindowYears"] = int(
    st.number_input("Screening period (years)", min_value=1, max_value=10, step=1, key="general_window")
)

if updated != intervention:
    set_intervention(updated)
    st.rerun()

changed = [
    {"Setting": FIELD_LABELS[key], "Demonstration default": defaults[key], "Current": updated[key], "Source": USER_DEFINED_MARK}
    for key in FIELD_LABELS
    if updated[key] != defaults[key]
]
if changed:
    st.markdown(f"**{USER_DEFINED_MARK} intervention settings**")
    st.dataframe(arrow_safe_dataframe(changed), use_container_width=True, hide_index=True)
else:
    st.caption("All intervention settings are demonstration defaults.")
st.caption(
    "Test accuracy, treatment completion, adverse-event and regimen-efficacy assumptions are demonstration "
    "values. See Evidence and technical information."
)
page_link("general_pages/3_Run_analysis.py", label="Continue to Run analysis")
