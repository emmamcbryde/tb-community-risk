"""Country-data, trend and local-upload panels for the Set up population page."""

from __future__ import annotations

from dataclasses import replace
from typing import Any, Sequence

import altair as alt
import pandas as pd
import streamlit as st

from app.display import arrow_safe_dataframe
from app.general.state import (
    CANDIDATE_KEY,
    EDITOR_VERSION_KEY,
    get_profile,
    get_trend_settings,
    set_profile,
    set_trend_settings,
)
from engine.profiles.country import (
    KEEP_CURRENT,
    UNCHANGED_BY_COUNTRY_DATA,
    USE_NEW,
    ConflictRequiresChoice,
    apply_snapshot_country,
    preview_country_application,
)
from engine.profiles.local_incidence import (
    FIELD_GUIDE,
    TEMPLATE_CSV,
    apply_local_incidence,
    parse_local_incidence,
)
from engine.profiles.population_profile import IncidencePoint, PopulationProfile, Provenance, TrendSpec
from engine.who_incidence.snapshot import IncidenceSnapshot
from engine.who_incidence.trend import (
    DESCRIPTIVE_STATEMENT,
    NO_CLEAR_CHANGE_NOTE,
    CovidHandling,
    TrendMethod,
    TrendResult,
    TrendSettings,
    estimate_trend,
)


PLACEHOLDER = "Select a country or area"
CANDIDATE_MEMORY_KEY = "general_country_candidate_memory"
METHOD_LABELS = {
    TrendMethod.LOG_LINEAR_RECENT: "Log-linear trend (primary)",
    TrendMethod.PENALISED_SPLINE: "Penalised smooth (sensitivity)",
}
WINDOW_LABELS = {5: "Last 5 years", 10: "Last 10 years", 15: "Last 15 years", None: "Full series", "custom": "Custom period"}
COVID_LABELS = {
    "include": "Include all WHO estimates",
    "exclude": "Exclude 2020-2022",
    "segmented": "Segmented disruption adjustment (2020-2022)",
    "custom": "Exclude selected years",
}
CONFLICT_LABELS = {USE_NEW: "Replace with the new data", KEEP_CURRENT: "Keep the current values"}


def render_country_section(snapshot: IncidenceSnapshot | None) -> None:
    st.subheader("Country data")
    if snapshot is None:
        message = st.session_state.get("general_snapshot_error") or "No WHO incidence snapshot is installed."
        st.error(f"Country data are unavailable: {message}")
        st.caption("The demonstration profile and local incidence files can still be used.")
        return
    manifest = snapshot.manifest
    scope = (
        f"{manifest['coverage']['countriesAndAreas']} countries and areas"
        if snapshot.is_complete_dataset
        else f"offline example subset of {len(snapshot.countries())} countries, not the complete dataset"
    )
    st.caption(
        f"WHO TB burden estimates, Global Tuberculosis Report {snapshot.report_year} round · "
        f"accessed {(manifest.get('volatile') or {}).get('accessDate')} · {scope}."
    )
    options = [PLACEHOLDER] + [f"{item['name']} ({item['iso3']})" for item in snapshot.countries()]
    codes = {f"{item['name']} ({item['iso3']})": item["iso3"] for item in snapshot.countries()}
    if st.session_state.get(CANDIDATE_KEY) not in options:
        # Widget state is discarded when the user visits another page; restore the last preview.
        remembered = st.session_state.get(CANDIDATE_MEMORY_KEY)
        st.session_state[CANDIDATE_KEY] = remembered if remembered in options else PLACEHOLDER
    choice = st.selectbox(
        "Country or area",
        options,
        key=CANDIDATE_KEY,
        help="Choosing a country previews its WHO estimates. Nothing changes until you apply it.",
    )
    st.session_state[CANDIDATE_MEMORY_KEY] = choice
    profile = get_profile()
    if choice == PLACEHOLDER:
        if profile.incidence.series:
            st.markdown(f"**Current incidence data:** {profile.location.name} ({_provenance_text(profile)})")
            render_trend_section(profile.incidence.series, source_label=_series_label(profile), store_in_profile=True)
        else:
            st.info("No country selected. The demonstration profile uses bundled calibration targets, not country data.")
        return
    iso3 = codes[choice]
    summary = snapshot.country_summary(iso3)
    series = snapshot.series_for(iso3)
    applied = (
        profile.location.iso3 == iso3
        and profile.incidence.provenance is Provenance.WHO_SNAPSHOT
        and profile.incidence.snapshot_id == snapshot.snapshot_id
    )
    _render_country_summary(summary)
    if not series:
        st.warning(f"WHO publishes no incidence estimates for {summary['name']} in this report round.")
        return
    if len(summary["estimatedYears"]) < summary["years"][1] - summary["years"][0] + 1 or summary["years"][1] < snapshot.report_year - 1:
        st.warning(
            f"Estimates are available for {len(summary['estimatedYears'])} year(s), "
            f"{summary['estimatedYears'][0]}-{summary['estimatedYears'][-1]}; the series is incomplete."
        )
    render_trend_section(series, source_label="WHO", store_in_profile=applied)
    if applied:
        st.success(f"{summary['name']} incidence data are applied to this profile.")
        return
    _render_apply_controls(profile, snapshot, iso3)


def _render_country_summary(summary: dict[str, Any]) -> None:
    latest = summary["latest"]
    cols = st.columns(3)
    cols[0].metric("WHO region", summary["region"])
    cols[1].metric("Years with estimates", len(summary["estimatedYears"]))
    if latest:
        cols[2].metric(f"Estimated incidence {latest['year']}", f"{latest['estimate']:g} per 100,000")
        st.caption(
            f"Uncertainty interval {latest['lower']:g}-{latest['upper']:g} per 100,000; about "
            f"{latest['cases']:,.0f} incident cases in a population of {latest['population']:,.0f} ({latest['year']})."
        )


def _render_apply_controls(profile: PopulationProfile, snapshot: IncidenceSnapshot, iso3: str) -> None:
    changes = preview_country_application(profile, snapshot, iso3)
    st.markdown("**Applying these data will change only:**")
    st.dataframe(arrow_safe_dataframe([change.to_row() for change in changes]), use_container_width=True, hide_index=True)
    st.caption("Not changed by country data: " + "; ".join(UNCHANGED_BY_COUNTRY_DATA) + ".")
    resolutions: dict[str, str] = {}
    for change in changes:
        if not change.conflict:
            continue
        resolutions[change.key] = st.radio(
            f"{change.label}: the current value is user-supplied ({change.current}).",
            list(CONFLICT_LABELS),
            format_func=CONFLICT_LABELS.get,
            index=None,
            key=f"general_conflict_{change.key}",
        )
    missing_choice = any(value is None for value in resolutions.values())
    if missing_choice:
        st.info("Choose what to do with each user-supplied value before applying.")
    if st.button("Apply selected country incidence data", type="primary", disabled=missing_choice):
        try:
            updated = apply_snapshot_country(profile, snapshot, iso3, resolutions={k: v for k, v in resolutions.items() if v})
        except ConflictRequiresChoice as exc:
            st.error(str(exc))
            return
        settings = get_trend_settings()
        updated = replace(updated, trend=TrendSpec(method=settings.method.value, settings=tuple(sorted(_freeze(settings.to_dict()).items()))))
        set_profile(updated)
        st.session_state["general_apply_message"] = (
            f"Applied {updated.location.name} WHO incidence data (snapshot {snapshot.snapshot_id}). "
            "Other inputs keep their current values and sources."
        )
        st.rerun()


def render_trend_section(series: Sequence[IncidencePoint], *, source_label: str, store_in_profile: bool) -> TrendResult | None:
    settings = _trend_controls(series, source_label)
    try:
        result = estimate_trend(series, settings)
    except ValueError as exc:
        st.error(str(exc))
        return None
    if store_in_profile:
        profile = get_profile()
        spec = TrendSpec(method=settings.method.value, settings=tuple(sorted(_freeze(settings.to_dict()).items())))
        if profile.trend != spec:
            set_profile(replace(profile, trend=spec))
    st.altair_chart(incidence_chart(series, result, source_label=source_label), use_container_width=True)
    st.caption(chart_caption(result))
    _render_trend_result(result, source_label)
    return result


def _trend_controls(series: Sequence[IncidencePoint], source_label: str = "WHO") -> TrendSettings:
    settings = get_trend_settings()
    years = sorted(point.year for point in series)
    cols = st.columns(3)
    method = cols[0].radio(
        "Trend method",
        list(METHOD_LABELS),
        format_func=METHOD_LABELS.get,
        index=list(METHOD_LABELS).index(settings.method) if settings.method in METHOD_LABELS else 0,
        key="general_trend_method",
    )
    window_options = [5, 10, 15, None, "custom"]
    current_window = "custom" if settings.start_year is not None else settings.window_years
    window = cols[1].selectbox(
        "Fitting period",
        window_options,
        index=window_options.index(current_window) if current_window in window_options else 1,
        format_func=WINDOW_LABELS.get,
        key="general_trend_window",
    )
    start = end = None
    if window == "custom" and years:
        start = int(cols[1].number_input("First year", min_value=years[0], max_value=years[-1], value=settings.start_year or max(years[0], years[-1] - 9), key="general_trend_start"))
        end = int(cols[1].number_input("Last year", min_value=years[0], max_value=years[-1], value=settings.end_year or years[-1], key="general_trend_end"))
    covid_options = list(COVID_LABELS)
    if method is TrendMethod.PENALISED_SPLINE:
        covid_options.remove("segmented")
    current_covid = "custom" if settings.excluded_years else settings.covid_handling.value
    covid = cols[2].selectbox(
        "COVID-era years",
        covid_options,
        index=covid_options.index(current_covid) if current_covid in covid_options else 0,
        format_func=lambda code: _covid_label(code, source_label),
        key="general_trend_covid",
    )
    excluded: tuple[int, ...] = ()
    if covid == "custom":
        excluded = tuple(sorted(cols[2].multiselect("Years to exclude", years, default=[y for y in settings.excluded_years if y in years], key="general_trend_excluded")))
    updated = TrendSettings(
        method=method,
        window_years=None if window in (None, "custom") else int(window),
        start_year=start,
        end_year=end,
        covid_handling=CovidHandling.INCLUDE if covid == "custom" else CovidHandling(covid),
        excluded_years=excluded,
        seed=settings.seed,
        draws=settings.draws,
    )
    if updated != settings:
        set_trend_settings(updated)
    return updated


def _render_trend_result(result: TrendResult, source_label: str = "WHO") -> None:
    covid_text = _covid_label("custom" if result.settings.excluded_years else result.settings.covid_handling.value, source_label)
    if result.status != "estimated":
        st.warning("Trend not estimated: " + " ".join(result.warnings))
        return
    period = f"{result.period[0]}-{result.period[1]}"
    interval = result.propagated_interval
    top = st.columns(2)
    top[0].metric("Annual change in estimated incidence", f"{result.annual_percent_change:+.1f}% per year")
    top[1].metric(
        "Propagated uncertainty interval" if interval else "Regression interval",
        _interval_text(interval or result.fit_interval),
    )
    st.markdown(
        f"**{result.summary_label}** · fitting period {period} ({result.years_used} years) · "
        f"{METHOD_LABELS[result.method].split(' (')[0].lower()} · COVID-era years: {covid_text[0].lower() + covid_text[1:]} · "
        f"fit: {result.classification_label.lower()}"
    )
    if result.summary == "no_clear_change":
        st.caption(NO_CLEAR_CHANGE_NOTE)
    st.caption(DESCRIPTIVE_STATEMENT)
    if result.warnings:
        st.warning(" ".join(result.warnings))
    with st.expander("Trend diagnostics"):
        diag = result.diagnostics
        rows = [
            {"Diagnostic": "Model", "Value": diag.get("model")},
            {"Diagnostic": "Regression (fitting) interval for annual change", "Value": _interval_text(result.fit_interval)},
            {"Diagnostic": "Propagated uncertainty interval (WHO bounds)", "Value": _interval_text(result.propagated_interval)},
            {"Diagnostic": "R-squared (log scale)", "Value": _fmt(diag.get("rSquared"), 3)},
            {"Diagnostic": "Residual SD (log scale)", "Value": _fmt(diag.get("residual_sd"), 3)},
            {"Diagnostic": "Durbin-Watson", "Value": _fmt(diag.get("durbinWatson"), 2)},
            {"Diagnostic": "Curvature test p-value", "Value": _fmt(diag.get("curvature_p"), 3)},
            {"Diagnostic": "Smoothing parameter", "Value": _fmt(diag.get("smoothingParameter"), 3)},
            {"Diagnostic": "Effective degrees of freedom", "Value": _fmt(diag.get("effectiveDegreesOfFreedom"), 2)},
            {"Diagnostic": "Disruption level shift", "Value": _fmt(diag.get("disruptionLevelShiftPercent"), 1, suffix="%")},
            {"Diagnostic": "Missing years in period", "Value": ", ".join(str(y) for y in diag.get("missingYears") or []) or "None"},
            {"Diagnostic": "Excluded years", "Value": ", ".join(str(y) for y in diag.get("excludedYears") or []) or "None"},
        ]
        sensitivity = diag.get("disruptionSensitivity")
        if sensitivity:
            rows.append(
                {
                    "Diagnostic": "Annual change with / without 2020-2022",
                    "Value": f"{sensitivity['apcIncluded']:+.1f}% / {sensitivity['apcExcluded']:+.1f}%",
                }
            )
        propagation = diag.get("propagation") or {}
        if propagation:
            rows.append({"Diagnostic": "Uncertainty draws and seed", "Value": f"{propagation['draws']} draws, seed {propagation['seed']} ({propagation['yearCorrelation']} years)"})
        st.dataframe(arrow_safe_dataframe([row for row in rows if row["Value"] not in (None, "")]), use_container_width=True, hide_index=True)
        st.caption(
            "The propagated interval reflects only the published WHO uncertainty bounds. It excludes WHO methodological "
            "uncertainty, structural change, model structure and infection-pressure uncertainty."
        )


def incidence_chart(series: Sequence[IncidencePoint], result: TrendResult | None, *, source_label: str) -> alt.LayerChart:
    estimate_name = f"{source_label} estimate"
    band_name = f"{source_label} uncertainty interval"
    fitted_name = "Fitted trend"
    points = {p.year: p for p in (result.points if result else [])}
    frame = pd.DataFrame(
        [
            {
                "Year": p.year,
                "Estimate": p.estimate,
                "Lower": p.lower,
                "Upper": p.upper,
                "Used in fit": "Yes" if points.get(p.year) and points[p.year].used_in_fit else "No",
                "series": estimate_name,
            }
            for p in series
        ]
    )
    domain = [estimate_name, band_name, fitted_name]
    scale = alt.Scale(domain=domain, range=["#1f5fa8", "#9ecae1", "#c2410c"])
    x = alt.X("Year:Q", axis=alt.Axis(format="d", title="Year"))
    layers = []
    band = frame.dropna(subset=["Lower", "Upper"]).assign(series=band_name)
    if not band.empty:
        layers.append(
            alt.Chart(band).mark_area(opacity=0.45).encode(
                x=x, y=alt.Y("Lower:Q", title="Incidence per 100,000 per year"), y2="Upper:Q", color=alt.Color("series:N", scale=scale, legend=alt.Legend(title=None, orient="bottom"))
            )
        )
    if result and result.settings.disruption_years and any(p.year in result.settings.disruption_years for p in series):
        years = [y for y in result.settings.disruption_years if any(p.year == y for p in series)]
        rect = pd.DataFrame([{"start": min(years) - 0.5, "end": max(years) + 0.5}])
        layers.append(alt.Chart(rect).mark_rect(opacity=0.08, color="#6b7280").encode(x="start:Q", x2="end:Q"))
    tooltip = ["Year:Q", alt.Tooltip("Estimate:Q", format=".3g"), alt.Tooltip("Lower:Q", format=".3g"), alt.Tooltip("Upper:Q", format=".3g"), "Used in fit:N"]
    color = alt.Color("series:N", scale=scale, legend=alt.Legend(title=None, orient="bottom"))
    layers.append(alt.Chart(frame).mark_line().encode(x=x, y=alt.Y("Estimate:Q", title="Incidence per 100,000 per year"), color=color))
    used = frame[frame["Used in fit"] == "Yes"] if result and result.status == "estimated" else frame
    unused = frame[frame["Used in fit"] == "No"] if result and result.status == "estimated" else frame.iloc[0:0]
    layers.append(alt.Chart(used).mark_point(filled=True, size=45).encode(x=x, y="Estimate:Q", color=color, tooltip=tooltip))
    if not unused.empty:
        layers.append(alt.Chart(unused).mark_point(filled=False, size=45, strokeWidth=1.5).encode(x=x, y="Estimate:Q", color=color, tooltip=tooltip))
    if result and result.status == "estimated":
        fitted = pd.DataFrame([{"Year": p.year, "Fitted": p.fitted, "series": fitted_name} for p in result.points if p.fitted is not None])
        if not fitted.empty:
            layers.append(
                alt.Chart(fitted).mark_line(strokeDash=[6, 3], strokeWidth=2.5).encode(
                    x=x, y="Fitted:Q", color=alt.Color("series:N", scale=scale, legend=alt.Legend(title=None, orient="bottom"))
                )
            )
    return alt.layer(*layers).properties(height=320)


def chart_caption(result: TrendResult | None) -> str:
    text = (
        "Shaded band: published uncertainty interval. Dashed line: fitted trend over the fitting period. "
        "Filled points were used in the fit; hollow points were not (outside the period or excluded)."
    )
    if result and result.settings.disruption_years:
        text += " Grey background: 2020-2022 (COVID-era years)."
    return text


def render_upload_section() -> None:
    with st.expander("Use a local or subnational incidence file"):
        st.caption(
            "Upload estimated TB disease incidence for one location. The file stays in this session and is recorded "
            "as a local file, never as WHO data. Notification rates are not accepted as incidence."
        )
        st.download_button("Download CSV template", data=TEMPLATE_CSV, file_name="local_incidence_template.csv", mime="text/csv")
        st.dataframe(arrow_safe_dataframe([{"Column": k, "Meaning": v} for k, v in FIELD_GUIDE.items()]), use_container_width=True, hide_index=True)
        upload = st.file_uploader("Incidence file (CSV)", type=["csv"], key=f"general_incidence_upload_{st.session_state[EDITOR_VERSION_KEY]}")
        if upload is None:
            return
        result = parse_local_incidence(upload.getvalue(), filename=upload.name)
        for message in result.warnings:
            st.warning(message)
        if not result.is_valid:
            st.error("The file cannot be used:")
            for message in result.errors[:12]:
                st.write(f"- {message}")
            return
        st.markdown(f"**Preview: {result.location}** ({len(result.points)} years; sources: {'; '.join(result.sources)})")
        st.dataframe(arrow_safe_dataframe(result.preview_rows()), use_container_width=True, hide_index=True)
        profile = get_profile()
        choice = USE_NEW
        if profile.incidence.provenance in {Provenance.WHO_SNAPSHOT, Provenance.LOCAL_UPLOAD}:
            choice = st.radio(
                f"This profile already has incidence data ({_provenance_text(profile)}).",
                list(CONFLICT_LABELS),
                format_func=CONFLICT_LABELS.get,
                index=None,
                key="general_upload_conflict",
            )
        if st.button("Apply local incidence data", disabled=choice is None):
            set_profile(apply_local_incidence(profile, result, resolutions={"incidence": choice}))
            st.session_state["general_apply_message"] = f"Applied local incidence data from {result.filename}."
            st.rerun()


def _covid_label(code: str, source_label: str) -> str:
    if code == "include":
        return "Include all WHO estimates" if source_label == "WHO" else "Include all estimates"
    return COVID_LABELS.get(code, code)


def _provenance_text(profile: PopulationProfile) -> str:
    if profile.incidence.provenance is Provenance.WHO_SNAPSHOT:
        return f"WHO snapshot {profile.incidence.snapshot_id}"
    if profile.incidence.provenance is Provenance.LOCAL_UPLOAD:
        return profile.incidence.source
    return "no incidence data"


def _series_label(profile: PopulationProfile) -> str:
    return "WHO" if profile.incidence.provenance is Provenance.WHO_SNAPSHOT else "Local"


def _interval_text(interval) -> str:
    if not interval:
        return "Not available"
    return f"{interval[0]:+.1f}% to {interval[1]:+.1f}%"


def _fmt(value, digits: int, suffix: str = "") -> str | None:
    if value is None:
        return None
    return f"{value:.{digits}f}{suffix}"


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _freeze(item) for key, item in value.items()}
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value
