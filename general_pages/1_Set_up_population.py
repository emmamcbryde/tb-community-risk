from __future__ import annotations

from dataclasses import replace

import pandas as pd
import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.country_panel import render_country_section, render_upload_section
from app.general.profile_editing import (
    EFFECT_MEASURE_OPTIONS,
    REVIEW_LABELS,
    RISK_FACTOR_COLUMNS,
    TRANSITION_LABELS,
    apply_risk_factor_rows,
    parse_risk_factor_csv,
    population_status_text,
    risk_factor_editing_issues,
    risk_factor_rows,
    risk_factor_template_csv,
    risk_factors_csv,
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
from engine.profiles.demonstration import DEMONSTRATION_WARNING
from engine.profiles.effect_measures import effect_warnings
from engine.profiles.population_profile import (
    LocationKind,
    unresolved_inputs,
    user_override_fields,
    validate_profile,
    with_population_size,
)


init_general_state()
snapshot = bundled_snapshot()

st.title("Set up population")
for key in ("general_restore_message", "general_apply_message"):
    message = st.session_state.pop(key, "")
    if message:
        st.success(message)

profile = get_profile()
st.markdown(f"**Profile:** {profile.name}")
st.warning(DEMONSTRATION_WARNING)
st.button(RESTORE_DEFAULTS_LABEL, on_click=restore_demonstration_defaults)

render_country_section(snapshot)
render_upload_section()
profile = get_profile()

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
if profile.location.national_population and profile.location.kind is not LocationKind.DEMONSTRATION:
    st.caption(
        f"National population of {profile.location.name}: {profile.location.national_population:,.0f} "
        f"({profile.location.national_population_year}; {profile.location.national_population_source}). "
        "Shown for context only; it does not change the simulated population."
    )

st.markdown("**Age distribution**")
st.dataframe(
    arrow_safe_dataframe(
        [
            {
                "Age group": band.label,
                "Proportion": f"{band.proportion.value * 100:.1f}%" if band.proportion.value is not None else "",
                "Status": value_status_text(band.proportion),
            }
            for band in profile.age_distribution
        ]
    ),
    use_container_width=True,
    hide_index=True,
)
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
    "Edit a cell to create a user-defined value (marked " + USER_DEFINED_MARK + "). Leave a value blank to record it as "
    "missing; enter 0 for a true zero. Untick Enabled to run without stratifying by that factor. The effect-measure "
    "type (RR, HR or OR) is kept exactly as entered. Rows can be added or deleted."
)
editor_rows = risk_factor_rows(profile)
edited = st.data_editor(
    pd.DataFrame(editor_rows, columns=RISK_FACTOR_COLUMNS),
    key=f"general_risk_editor_{st.session_state[EDITOR_VERSION_KEY]}",
    num_rows="dynamic",
    hide_index=True,
    use_container_width=True,
    column_order=[column for column in RISK_FACTOR_COLUMNS if column != "id"],
    disabled=["Status", "id"],
    column_config={
        "Enabled": st.column_config.CheckboxColumn("Enabled", default=True),
        "Status": st.column_config.TextColumn("Status", help=f"{USER_DEFINED_MARK} marks values you have changed."),
        "Prevalence (%)": st.column_config.NumberColumn("Prevalence (%)", min_value=0.0, max_value=100.0, format="%.1f"),
        "Prevalence low (%)": st.column_config.NumberColumn("Prevalence low (%)", min_value=0.0, max_value=100.0, format="%.1f"),
        "Prevalence high (%)": st.column_config.NumberColumn("Prevalence high (%)", min_value=0.0, max_value=100.0, format="%.1f"),
        "Effect estimate": st.column_config.NumberColumn("Effect estimate", min_value=0.0, format="%.2f"),
        "Effect low": st.column_config.NumberColumn("Effect low", min_value=0.0, format="%.2f"),
        "Effect high": st.column_config.NumberColumn("Effect high", min_value=0.0, format="%.2f"),
        "Effect measure": st.column_config.SelectboxColumn("Effect measure", options=EFFECT_MEASURE_OPTIONS),
        "Affected transition": st.column_config.SelectboxColumn("Affected transition", options=list(TRANSITION_LABELS.values())),
        "Evidence year": st.column_config.NumberColumn("Evidence year", min_value=1900, max_value=2100, format="%d"),
        "Review status": st.column_config.SelectboxColumn("Review status", options=list(REVIEW_LABELS.values())),
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

cols = st.columns(3)
cols[0].download_button("Export risk factors (CSV)", data=risk_factors_csv(profile), file_name="risk_factors.csv", mime="text/csv")
cols[1].download_button("Blank risk-factor template", data=risk_factor_template_csv(), file_name="risk_factor_template.csv", mime="text/csv")
if profile.risk_factors and cols[2].button("Run without risk-factor stratification"):
    set_profile(replace(profile, risk_factors=()))
    st.session_state[EDITOR_VERSION_KEY] = int(st.session_state[EDITOR_VERSION_KEY]) + 1
    st.rerun()
with st.expander("Import risk factors from CSV"):
    imported = st.file_uploader("Risk-factor CSV", type=["csv"], key=f"general_risk_upload_{st.session_state[EDITOR_VERSION_KEY]}")
    if imported is not None:
        rows, errors = parse_risk_factor_csv(imported.getvalue())
        if errors:
            st.error("The file cannot be used: " + " ".join(errors[:8]))
        else:
            st.dataframe(arrow_safe_dataframe(rows), use_container_width=True, hide_index=True)
            if st.button("Replace the risk-factor table with this file"):
                try:
                    set_profile(apply_risk_factor_rows(profile, rows))
                    st.session_state[EDITOR_VERSION_KEY] = int(st.session_state[EDITOR_VERSION_KEY]) + 1
                    st.rerun()
                except ValueError as exc:
                    st.error(str(exc))
if not profile.risk_factors:
    st.info("No risk factors are defined. The analysis will run without risk-factor stratification.")
for warning in effect_warnings(profile):
    st.caption(f"Note: {warning.message}")
for issue in risk_factor_editing_issues(profile):
    st.error(issue)

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
