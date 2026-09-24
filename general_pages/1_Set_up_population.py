from __future__ import annotations

from dataclasses import replace

import altair as alt
import pandas as pd
import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.profile_editing import (
    EFFECT_MEASURE_OPTIONS,
    REVIEW_LABELS,
    RISK_FACTOR_COLUMNS,
    apply_risk_factor_rows,
    population_status_text,
    risk_factor_rows,
    value_status_text,
)
from app.general.state import (
    EDITOR_VERSION_KEY,
    bundled_snapshot,
    get_profile,
    init_general_state,
    page_link,
    restore_demonstration_defaults,
    set_profile,
)
from app.general.terminology import RESTORE_DEFAULTS_LABEL, USER_DEFINED_MARK
from engine.profiles.country import apply_snapshot_country, apply_user_incidence
from engine.profiles.demonstration import (
    DEMONSTRATION_PROFILE_LABEL,
    DEMONSTRATION_WARNING,
    build_demonstration_profile,
)
from engine.profiles.population_profile import (
    Location,
    LocationKind,
    Provenance,
    unresolved_inputs,
    user_override_fields,
    validate_profile,
    with_population_size,
)
from engine.who_incidence.adapters import parse_user_incidence_csv


DEMONSTRATION_OPTION = "Demonstration working defaults (no country selected)"

init_general_state()
snapshot = bundled_snapshot()

st.title("Set up population")
restore_message = st.session_state.pop("general_restore_message", "")
if restore_message:
    st.success(restore_message)

profile = get_profile()
st.markdown(f"**Profile:** {profile.name}")
st.warning(DEMONSTRATION_WARNING)
st.button(RESTORE_DEFAULTS_LABEL, on_click=restore_demonstration_defaults)

# ---------------------------------------------------------------------------
# Country or demonstration profile
# ---------------------------------------------------------------------------
st.subheader("Country or demonstration profile")
options = [DEMONSTRATION_OPTION]
labels_to_iso3: dict[str, str] = {}
if snapshot is not None:
    for country in snapshot.countries():
        label = f"{country['name']} ({country['iso3']})"
        options.append(label)
        labels_to_iso3[label] = country["iso3"]
current_iso3 = profile.location.iso3 if profile.incidence.snapshot_id else None
current_label = next((label for label, code in labels_to_iso3.items() if code == current_iso3), DEMONSTRATION_OPTION)
if st.session_state.get("general_country_select") not in options:
    st.session_state["general_country_select"] = current_label
selection = st.selectbox(
    "Country or profile",
    options,
    key="general_country_select",
    help="Selecting a country attaches its WHO incidence estimates. Other inputs keep their demonstration values until you review them.",
)
if selection != current_label:
    if selection == DEMONSTRATION_OPTION:
        demo = build_demonstration_profile()
        profile = replace(
            profile,
            profile_id=demo.profile_id,
            name=demo.name,
            location=Location(name=DEMONSTRATION_PROFILE_LABEL, kind=LocationKind.DEMONSTRATION),
            incidence=demo.incidence,
            data_vintage=demo.data_vintage,
        )
    else:
        profile = apply_snapshot_country(profile, snapshot, labels_to_iso3[selection])
    set_profile(profile)
    st.rerun()

if snapshot is not None:
    manifest = snapshot.manifest
    scope = (
        "complete dataset"
        if snapshot.is_complete_dataset
        else f"offline example snapshot of {len(snapshot.countries())} countries, not the complete WHO dataset"
    )
    st.caption(
        f"Incidence data: {manifest.get('sourceDataset')} · WHO report year {manifest.get('sourceReportYear')} "
        f"· extracted {manifest.get('extractionDate')} · {scope}."
    )
else:
    st.info("No bundled incidence snapshot is available. The demonstration profile can still be used.")

# ---------------------------------------------------------------------------
# Incidence time series
# ---------------------------------------------------------------------------
st.subheader("Estimated TB incidence")
series = profile.incidence.series
if series:
    frame = pd.DataFrame(
        [{"Year": p.year, "Estimate": p.estimate, "Lower": p.lower, "Upper": p.upper} for p in series]
    )
    base = alt.Chart(frame).encode(x=alt.X("Year:O", title="Year"))
    band = base.mark_area(opacity=0.25).encode(
        y=alt.Y("Lower:Q", title="Incidence per 100,000 per year"),
        y2="Upper:Q",
    )
    line = base.mark_line(point=True).encode(
        y="Estimate:Q",
        tooltip=["Year", alt.Tooltip("Estimate:Q", format=".1f"), alt.Tooltip("Lower:Q", format=".1f"), alt.Tooltip("Upper:Q", format=".1f")],
    )
    st.altair_chart(band + line, use_container_width=True)
    first, last = profile.incidence.data_year_range
    source_note = "User-defined file" if profile.incidence.provenance is Provenance.USER_DEFINED else "WHO estimate"
    st.caption(
        f"{profile.location.name}, {first}-{last}. {source_note}; shaded band shows the uncertainty interval. "
        "This is estimated TB disease incidence, not infection pressure. "
        "In this version the series is shown for review and does not yet change the model's epidemiology."
    )
    with st.expander("Incidence values"):
        st.dataframe(arrow_safe_dataframe(frame.to_dict(orient="records")), use_container_width=True, hide_index=True)
else:
    st.info(
        "No incidence series is linked to this profile. The demonstration profile uses bundled calibration "
        "targets that are not specific to any country."
    )

with st.expander("Use a local or subnational incidence file"):
    st.caption(
        "CSV columns: iso3 (optional for subnational areas), country (area name), year, incidence_per_100k, "
        "incidence_per_100k_lo, incidence_per_100k_hi. Values stay in this session only."
    )
    upload = st.file_uploader("Incidence file", type=["csv"], key=f"general_incidence_upload_{st.session_state[EDITOR_VERSION_KEY]}")
    if upload is not None:
        records, report = parse_user_incidence_csv(upload.getvalue().decode("utf-8-sig"))
        for issue in report.warnings:
            st.warning(issue.message)
        if not report.is_valid:
            st.error("The file could not be used:")
            for issue in report.errors[:10]:
                st.write(f"- {issue.message}")
        elif st.button("Use this incidence file"):
            try:
                set_profile(apply_user_incidence(get_profile(), records, source_label=upload.name))
                st.session_state.pop("general_country_select", None)
                st.rerun()
            except ValueError as exc:
                st.error(str(exc))

# ---------------------------------------------------------------------------
# Population
# ---------------------------------------------------------------------------
st.subheader("Population")
size_value = int(profile.population_size.value or 0)
if st.session_state.get("general_population_input") is None:
    st.session_state["general_population_input"] = size_value
new_size = st.number_input(
    "Population size to simulate",
    min_value=1,
    max_value=1_000_000,
    step=100,
    key="general_population_input",
    help="The number of people in the simulated community, not the national population.",
)
if int(new_size) != size_value:
    set_profile(with_population_size(profile, int(new_size)))
    st.rerun()
st.caption(f"Population size: {population_status_text(profile)}.")

st.markdown("**Age distribution**")
age_rows = [
    {
        "Age group": band.label,
        "Proportion": f"{band.proportion.value * 100:.1f}%" if band.proportion.value is not None else "",
        "Status": value_status_text(band.proportion),
    }
    for band in profile.age_distribution
]
st.dataframe(arrow_safe_dataframe(age_rows), use_container_width=True, hide_index=True)
st.caption(f"Source: {profile.age_distribution_source or 'Not recorded'}. Not specific to the selected country.")

ltbi = profile.ltbi_prevalence
st.markdown(
    f"**LTBI prevalence:** {ltbi.value * 100:.1f}% ({value_status_text(ltbi)})"
    if ltbi.value is not None
    else f"**LTBI prevalence:** {value_status_text(ltbi)}"
)

# ---------------------------------------------------------------------------
# Risk factors
# ---------------------------------------------------------------------------
st.subheader("Risk factors (optional)")
st.caption(
    "Edit a cell to create a user-defined value. Leave a value blank to record it as missing; "
    "enter 0 for a true zero. Untick Enabled to run without stratifying by that factor. "
    "The effect-measure type (RR, HR or OR) is kept exactly as entered. Rows can be added or deleted."
)
editor_rows = risk_factor_rows(profile)
edited = st.data_editor(
    pd.DataFrame(editor_rows, columns=RISK_FACTOR_COLUMNS),
    key=f"general_risk_editor_{st.session_state[EDITOR_VERSION_KEY]}",
    num_rows="dynamic",
    hide_index=True,
    use_container_width=True,
    column_order=["Enabled", "Risk factor", "Status", "Prevalence (%)", "Effect estimate", "Effect measure", "Review status", "Source", "Notes"],
    disabled=["Status", "id"],
    column_config={
        "Enabled": st.column_config.CheckboxColumn("Enabled", default=True),
        "Prevalence (%)": st.column_config.NumberColumn("Prevalence (%)", min_value=0.0, max_value=100.0, format="%.1f"),
        "Effect estimate": st.column_config.NumberColumn("Effect estimate", min_value=0.0, format="%.2f"),
        "Effect measure": st.column_config.SelectboxColumn("Effect measure", options=EFFECT_MEASURE_OPTIONS),
        "Review status": st.column_config.SelectboxColumn("Review status", options=list(REVIEW_LABELS.values())),
        "Status": st.column_config.TextColumn("Status", help=f"{USER_DEFINED_MARK} marks values you have changed."),
    },
)
edited_rows = edited.to_dict(orient="records") if hasattr(edited, "to_dict") else list(edited)
try:
    updated_profile = apply_risk_factor_rows(profile, edited_rows)
except ValueError as exc:
    st.error(str(exc))
    updated_profile = profile
if updated_profile != profile:
    set_profile(updated_profile)
    st.rerun()
if not profile.risk_factors:
    st.info("No risk factors are defined. The analysis will run without risk-factor stratification.")

# ---------------------------------------------------------------------------
# Overrides and unresolved local evidence
# ---------------------------------------------------------------------------
overrides = user_override_fields(profile)
if overrides:
    st.markdown(f"**{USER_DEFINED_MARK} values in this profile:** " + "; ".join(overrides))

issues = validate_profile(profile)
if issues:
    st.error("Some inputs need correction before the analysis can run:")
    for issue in issues:
        st.write(f"- {issue['message']} ({issue['field']})")

unresolved = unresolved_inputs(profile)
st.subheader("Inputs still needing local evidence")
if unresolved:
    st.caption(f"{len(unresolved)} input(s) should be reviewed or replaced with locally applicable evidence.")
    st.dataframe(arrow_safe_dataframe(unresolved), use_container_width=True, hide_index=True)
else:
    st.success("All inputs have been reviewed.")

page_link("general_pages/2_Define_intervention.py", label="Continue to Define intervention")
