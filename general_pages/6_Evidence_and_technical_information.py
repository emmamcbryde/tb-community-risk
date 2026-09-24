from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.state import bundled_snapshot, get_profile, init_general_state
from app.general.terminology import PLANNING_PRINCIPLE, display_text
from engine.profiles.engine_mapping import effect_interpretation_rows
from engine.profiles.population_profile import unresolved_inputs
from engine.who_incidence.trend import (
    EPIDEMIOLOGICAL_QUANTITIES,
    IMPLEMENTED_METHODS,
    INCIDENCE_TO_INFECTION_POLICY,
    TrendMethod,
)


QUANTITY_LABELS = {
    "estimated_tb_disease_incidence": "Estimated TB disease incidence",
    "infection_pressure": "Infection pressure",
    "force_of_infection": "Force of infection",
    "progression_to_disease": "Progression from infection to disease",
    "case_detection_and_notification": "Case detection and notification",
}
TREND_LABELS = {
    TrendMethod.NOT_ESTIMATED: "Observed values only (current)",
    TrendMethod.LOG_LINEAR_RECENT: "Log-linear recent trend",
    TrendMethod.PENALISED_SPLINE: "Penalised spline",
    TrendMethod.STATE_SPACE: "State-space trend",
    TrendMethod.USER_OVERRIDE: "User-specified annual percentage change",
}

init_general_state()
profile = get_profile()
snapshot = bundled_snapshot()

st.title("Evidence and technical information")
st.info(PLANNING_PRINCIPLE)

st.subheader("Inputs still needing local evidence")
unresolved = unresolved_inputs(profile)
if unresolved:
    st.dataframe(arrow_safe_dataframe(unresolved), use_container_width=True, hide_index=True)
else:
    st.success("All inputs have been reviewed.")

st.subheader("Incidence data source")
if snapshot is not None:
    summary = snapshot.provenance_summary()
    coverage = summary.get("coverage") or {}
    rows = [
        {"Item": "Dataset", "Value": summary["dataset"]},
        {"Item": "Report round", "Value": f"Global Tuberculosis Report {summary['reportYear']}"},
        {"Item": "Snapshot", "Value": summary["snapshotId"]},
        {"Item": "Snapshot scope", "Value": "Complete dataset" if summary["isCompleteDataset"] else "Offline example subset (not the complete dataset)"},
        {"Item": "Coverage", "Value": f"{coverage.get('countriesAndAreas')} countries and areas, {coverage.get('yearRange')}"},
        {"Item": "Source", "Value": summary["sourceUrl"]},
        {"Item": "Access date", "Value": summary["accessDate"]},
        {"Item": "Cross-check", "Value": f"{summary['crossCheckRepository'] or 'none'} @ {(summary['crossCheckCommit'] or '')[:7]}"},
        {"Item": "Data file SHA-256", "Value": summary["dataSha256"]},
        {"Item": "Importer", "Value": summary["importerVersion"]},
        {"Item": "Citation", "Value": summary["citation"]},
    ]
    st.dataframe(arrow_safe_dataframe(rows), use_container_width=True, hide_index=True)
    st.caption(
        "WHO revises its whole retrospective series in each report round, so estimates from different "
        "report years are not combined."
    )
else:
    st.write("No incidence snapshot is bundled with this installation.")

st.subheader("Epidemiological quantities")
st.dataframe(
    arrow_safe_dataframe(
        [{"Quantity": QUANTITY_LABELS.get(key, key.replace("_", " ").capitalize()), "Meaning": text} for key, text in EPIDEMIOLOGICAL_QUANTITIES.items()]
    ),
    use_container_width=True,
    hide_index=True,
)
st.caption(INCIDENCE_TO_INFECTION_POLICY)

st.subheader("Incidence trend estimation")
st.dataframe(
    arrow_safe_dataframe(
        [
            {"Method": label, "Status": "Available" if method in IMPLEMENTED_METHODS else "Planned"}
            for method, label in TREND_LABELS.items()
        ]
    ),
    use_container_width=True,
    hide_index=True,
)
st.caption(
    "Planned methods will keep observed values and WHO uncertainty bounds alongside fitted values, allow the "
    "fitting period to be chosen, flag COVID-era disruption (2020-2022), and report annual percentage change."
)

st.subheader("How risk-factor effects are applied")
st.write(
    "The analysis engine multiplies the hazard of progression from infection to disease by each effect estimate. "
    "The declared measure type is kept exactly; an odds ratio or risk ratio is not converted to a hazard ratio."
)
rows = effect_interpretation_rows(profile)
if rows:
    st.dataframe(arrow_safe_dataframe(rows), use_container_width=True, hide_index=True)
else:
    st.write("This profile has no risk factors; the analysis runs without risk-factor stratification.")
st.caption(
    "User-entered prevalence is applied uniformly across age groups. Demonstration prevalence values are "
    "age-specific and are applied unchanged."
)

st.subheader("Profile record")
st.write(f"Profile: {profile.name}")
st.write(f"Schema: {profile.schema_version} · Fingerprint: {profile.profile_hash()[:16]}")
st.download_button(
    "Download population profile (JSON)",
    data=profile.to_json(),
    file_name=f"{profile.profile_id}.json",
    mime="application/json",
)
last_error = st.session_state.get("general_last_error")
if last_error:
    with st.expander("Last error"):
        st.code(display_text(last_error, fallback="An internal error occurred; see the application log."))
