from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.provenance_export import build_provenance_package, package_filename
from app.general.state import (
    ECONOMICS_CONFIG_KEY,
    ECONOMICS_KEY,
    RESULTS_KEY,
    bundled_snapshot,
    current_engine_config,
    get_analysis,
    get_profile,
    get_trend_settings,
    init_general_state,
)
from app.general.terminology import PLANNING_PRINCIPLE, display_text
from engine.profiles.effect_measures import ACTIVE_CONVERSION_POLICY, POLICY_STATUS, crosswalk_rows, effect_warnings
from engine.profiles.population_profile import unresolved_inputs
from engine.who_incidence.trend import (
    DESCRIPTIVE_STATEMENT,
    EPIDEMIOLOGICAL_QUANTITIES,
    IMPLEMENTED_METHODS,
    INCIDENCE_TO_INFECTION_POLICY,
    TrendMethod,
    estimate_trend,
)


QUANTITY_LABELS = {
    "estimated_tb_disease_incidence": "Estimated TB disease incidence",
    "notification_rate": "Notification rate",
    "active_tb_prevalence": "Prevalence of active TB",
    "force_of_infection": "Force of infection",
    "annual_risk_of_infection": "Annual risk of infection",
    "recent_infection": "Recent infection",
    "remote_infection": "Remote infection",
    "progression_to_disease": "Progression from infection to disease",
    "detection_and_treatment": "Detection and treatment",
}
TREND_LABELS = {
    TrendMethod.LOG_LINEAR_RECENT: ("Log-linear recent trend", "log(I_t) = a + b(t - mean t); annual change = 100(exp(b) - 1)"),
    TrendMethod.PENALISED_SPLINE: ("Penalised smooth (sensitivity)", "Whittaker-Henderson smoother on log incidence; slope = mean change over the last 3 years"),
    TrendMethod.USER_OVERRIDE: ("User-specified annual change", "Entered directly; not estimated from data"),
    TrendMethod.STATE_SPACE: ("State-space trend", "Designed, not implemented"),
}
UNCERTAINTY_ROWS = [
    {"Type": "WHO incidence-estimate uncertainty", "Where shown": "Incidence chart band; propagated trend interval", "Included in model results": "No"},
    {"Type": "Trend-fitting uncertainty", "Where shown": "Regression interval in trend diagnostics", "Included in model results": "No"},
    {"Type": "Stochastic community variation", "Where shown": "Simulation intervals in Results and Health economics", "Included in model results": "Yes (stochastic analysis only)"},
    {"Type": "Parameter uncertainty", "Where shown": "Not yet quantified", "Included in model results": "No"},
    {"Type": "Structural uncertainty", "Where shown": "Not yet quantified", "Included in model results": "No"},
]

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

st.subheader("Incidence data source and terms")
if snapshot is not None:
    summary = snapshot.provenance_summary()
    coverage = summary.get("coverage") or {}
    terms = summary.get("terms") or {}
    rows = [
        {"Item": "Dataset", "Value": summary["dataset"]},
        {"Item": "Report round", "Value": f"Global Tuberculosis Report {summary['reportYear']} (retrospective series; rounds are never combined)"},
        {"Item": "Snapshot", "Value": summary["snapshotId"]},
        {"Item": "Scope", "Value": "Complete dataset" if summary["isCompleteDataset"] else "Offline example subset (not the complete dataset)"},
        {"Item": "Coverage", "Value": f"{coverage.get('countriesAndAreas')} countries and areas ({coverage.get('withEstimates')} with estimates), {coverage.get('yearRange')}"},
        {"Item": "Source", "Value": summary["sourceUrl"]},
        {"Item": "Access date", "Value": summary["accessDate"]},
        {"Item": "Cross-check", "Value": f"{summary['crossCheckRepository'] or 'none'} at commit {(summary['crossCheckCommit'] or '')[:7]}"},
        {"Item": "Data file SHA-256", "Value": summary["dataSha256"]},
        {"Item": "Importer", "Value": summary["importerVersion"]},
    ]
    st.dataframe(arrow_safe_dataframe(rows), use_container_width=True, hide_index=True)
    st.markdown(f"**Citation:** {summary['citation']}")
    st.caption(f"Terms: {terms.get('summary', '')} {terms.get('noEndorsement', '')}")
else:
    st.error(st.session_state.get("general_snapshot_error") or "No incidence snapshot is installed.")

st.subheader("What the incidence data can and cannot tell the model")
st.write(DESCRIPTIVE_STATEMENT)
st.dataframe(
    arrow_safe_dataframe([{"Quantity": QUANTITY_LABELS.get(key, key), "Meaning": text} for key, text in EPIDEMIOLOGICAL_QUANTITIES.items()]),
    use_container_width=True,
    hide_index=True,
)
st.caption(INCIDENCE_TO_INFECTION_POLICY)

st.subheader("Kinds of uncertainty")
st.dataframe(arrow_safe_dataframe(UNCERTAINTY_ROWS), use_container_width=True, hide_index=True)

st.subheader("Incidence trend methods")
st.dataframe(
    arrow_safe_dataframe(
        [
            {"Method": label, "Formulation": formula, "Status": "Available" if method in IMPLEMENTED_METHODS else "Planned"}
            for method, (label, formula) in TREND_LABELS.items()
        ]
    ),
    use_container_width=True,
    hide_index=True,
)
st.caption(
    "WHO bounds are treated as the 2.5% and 97.5% points of a split-normal distribution on log incidence; each of "
    "1,000 draws (fixed seed) is refitted to give a propagated uncertainty interval for the annual change. "
    "COVID-era years are included by default; exclusion and a segmented temporary level shift are optional."
)

st.subheader("How risk-factor effects are applied")
st.write(
    "The analysis engine multiplies each person's progression hazard by the product of their enabled risk-factor "
    "effect estimates. Declared measure types are kept exactly; odds ratios and risk ratios are not converted. The "
    "baseline hazard is recalibrated so that average risk matches the calibration target."
)
rows = crosswalk_rows(profile)
if rows:
    st.dataframe(arrow_safe_dataframe(rows), use_container_width=True, hide_index=True)
else:
    st.write("This profile has no risk factors; the analysis runs without risk-factor stratification.")
for warning in effect_warnings(profile):
    st.caption(f"Warning: {warning.message}")
st.caption(f"Effect-measure conversion policy: {POLICY_STATUS[ACTIVE_CONVERSION_POLICY]}")
st.caption(
    "Country-specific risk-factor prevalence is not imported: the public WHO dataset publishes attributable cases, "
    "not prevalence or relative risks, and attributable fractions cannot be used as either."
)

st.subheader("Profile record and downloads")
st.write(f"Profile: {profile.name}")
st.write(f"Schema: {profile.schema_version} · Fingerprint: {profile.profile_hash()[:16]}")
st.download_button("Download population profile (JSON)", data=profile.to_json(), file_name=f"{profile.profile_id}.json", mime="application/json")
analysis = get_analysis()
analysis_label = f"stochastic-{analysis['nReps']}" if analysis["analysisMethod"] == "agent_based" else "deterministic"
if st.button("Prepare provenance package"):
    trend = None
    if profile.incidence.series:
        try:
            trend = estimate_trend(profile.incidence.series, get_trend_settings())
        except ValueError:
            trend = None
    config, _ = current_engine_config()
    st.session_state["general_provenance_package"] = {
        "name": package_filename(profile, analysis_label),
        "bytes": build_provenance_package(
            profile=profile,
            trend=trend,
            snapshot_manifest=snapshot.manifest if snapshot else None,
            engine_config=config,
            results_bundle=st.session_state.get(RESULTS_KEY),
            economics_config=st.session_state.get(ECONOMICS_CONFIG_KEY),
            economics_results=st.session_state.get(ECONOMICS_KEY),
        ),
    }
package = st.session_state.get("general_provenance_package")
if package:
    st.download_button("Download provenance package (ZIP)", data=package["bytes"], file_name=package["name"], mime="application/zip")
else:
    st.caption(
        "The package is built only when requested. It contains the profile, incidence and trend data, settings, "
        "manifest, risk factors, override audit, configuration, environment, results and limitations."
    )

last_error = st.session_state.get("general_last_error")
if last_error:
    with st.expander("Last error"):
        st.code(display_text(last_error, fallback="An internal error occurred; see the application log."))
