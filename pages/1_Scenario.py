from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from app.display import arrow_safe_dataframe, display_string
from app.epidemiology_inputs import (
    ADVANCED_RISK_FACTORS,
    AGE_GROUP_LABELS,
    PRINCIPAL_RISK_FACTORS,
    RISK_FACTOR_LABELS,
    apply_epidemiology_updates,
    fraction_to_percent,
    risk_override_from_percentages,
    risk_override_mode,
    risk_override_percent_values,
)
from app.state import (
    get_backend,
    get_backend_name,
    init_session_state,
    mark_config_changed,
    mark_validation_completed,
    record_message,
    set_backend_name,
    sync_backend_status,
)
from adapters.paths import scenarios_dir


init_session_state()

st.title("Scenario")
st.caption("APY v9 scenario setup for the selected backend.")

BACKEND_LABELS = {
    "matlab": "MATLAB v9 reference",
    "python_apy": "Python APY v9 port (experimental)",
}
BACKEND_NAMES = {label: name for name, label in BACKEND_LABELS.items()}

current_backend_name = get_backend_name()
selected_backend_label = st.selectbox(
    "APY backend",
    list(BACKEND_NAMES),
    index=list(BACKEND_NAMES).index(
        BACKEND_LABELS.get(current_backend_name, BACKEND_LABELS["matlab"])
    ),
)
selected_backend_name = BACKEND_NAMES[selected_backend_label]
if selected_backend_name != current_backend_name:
    set_backend_name(selected_backend_name)
    st.rerun()

backend = get_backend()

st.caption(
    "MATLAB remains the reference backend. The Python APY backend is experimental "
    "and has passed diagnostic validation across the committed scenario suite, "
    "but economics and attributable-risk add-ons are not yet ported."
)


def display_backend_status() -> None:
    status = backend.status()
    sync_backend_status(status)
    cols = st.columns(3)
    cols[0].metric(
        "Backend",
        BACKEND_LABELS.get(status.get("name"), status.get("name", "unknown")),
    )
    if status.get("name") == "matlab":
        cols[1].metric("MATLAB started", "yes" if status.get("started") else "no")
    else:
        cols[1].metric("MATLAB required", "no")
    cols[2].metric("Adapter", "error" if status.get("error") else "ready")
    if status.get("experimental"):
        st.warning("Python APY backend is experimental. MATLAB remains the reference backend.")
    if status.get("abm_path"):
        st.caption(f"ABM path: {status.get('abm_path', '')}")
    if status.get("error"):
        st.error(status["error"])


def config_overview_rows(config: dict) -> list[dict[str, object]]:
    fields = [
        "scenarioLabel",
        "configVersion",
        "modelVersion",
        "N",
        "nReps",
        "seed",
        "screenCoverage",
        "screenWindow",
        "followHorizon",
        "screeningStrategy",
        "testType",
        "regimen",
        "outputDir",
    ]
    return [
        {"field": field, "value": config.get(field)}
        for field in fields
        if field in config
    ]


def validation_rows(report: dict) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for group in ("errors", "warnings", "infos"):
        issues = report.get(group) or []
        if isinstance(issues, dict):
            issues = [issues]
        for issue in issues:
            if isinstance(issue, dict):
                rows.append(
                    {
                        "severity": display_string(issue.get("severity", group[:-1])),
                        "field": display_string(issue.get("fieldLabel", issue.get("field", ""))),
                        "config field": display_string(issue.get("fieldName", "")),
                        "code": display_string(issue.get("code", "")),
                        "message": display_string(issue.get("message", "")),
                    }
                )
    return rows


def choice_index(value: object, options: list[str], fallback: int = 0) -> int:
    try:
        return options.index(str(value))
    except ValueError:
        return fallback


def reset_run_state() -> None:
    st.session_state["validation_report"] = None
    st.session_state["results_bundle"] = None
    st.session_state["results_stale"] = False
    st.session_state["dirty_config"] = False


def scenario_path(filename: str) -> Path:
    base = scenarios_dir()
    base.mkdir(parents=True, exist_ok=True)
    name = filename.strip() or "streamlit_scenario.json"
    path = Path(name)
    if not path.is_absolute():
        path = base / path
    return path


st.subheader("Backend Status")
display_backend_status()

if st.button("Load backend defaults", type="primary"):
    try:
        st.session_state["config"] = backend.default_config()
        reset_run_state()
        sync_backend_status(backend.status())
        st.success("Default scenario loaded.")
    except Exception as exc:
        message = f"Could not load defaults: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

config = st.session_state.get("config")
if config:
    if st.session_state.get("dirty_config"):
        st.warning("Config has unvalidated changes.")
    if st.session_state.get("results_stale"):
        st.warning("Existing results are stale because config changed after the last run.")

    st.subheader("Edit Scenario")
    test_options = ["IGRA", "TST"]
    regimen_options = ["3HP", "4R", "3HR", "6H", "9H"]
    strategy_options = ["random", "ltbi", "cure", "prevent"]

    with st.form("scenario_edits"):
        scenario_label = st.text_input(
            "Scenario label",
            value=str(config.get("scenarioLabel", "")),
        )
        screen_coverage = st.number_input(
            "Screen coverage",
            min_value=0.0,
            max_value=1.0,
            value=float(config.get("screenCoverage", 0.0)),
            step=0.05,
        )
        test_type = st.selectbox(
            "Test type",
            test_options,
            index=choice_index(config.get("testType"), test_options),
        )
        regimen = st.selectbox(
            "Regimen",
            regimen_options,
            index=choice_index(config.get("regimen"), regimen_options),
        )
        screening_strategy = st.selectbox(
            "Screening strategy",
            strategy_options,
            index=choice_index(config.get("screeningStrategy"), strategy_options, 3),
        )
        st.subheader("Epidemiological assumptions")
        st.caption(
            "These controls set calibration targets for the simulated cohort. They are not the observed screening yield."
        )
        use_default_ltbi = st.checkbox(
            "Use default APY MTB infection prevalence",
            value=config.get("ltbiPrevalence") is None,
        )
        ltbi_prevalence_percent = None
        if not use_default_ltbi:
            ltbi_prevalence_percent = st.number_input(
                "Assumed prevalence of MTB infection (LTBI)",
                min_value=0.01,
                max_value=99.99,
                value=float(fraction_to_percent(config.get("ltbiPrevalence"), 47 / 624 * 100)),
                step=0.1,
                format="%.2f",
                help=(
                    "Displayed as a percentage and stored internally as a fraction. "
                    "Changing this recalibrates infection probability and affects expected yield, "
                    "false-positive burden and cases prevented."
                ),
            )
        else:
            st.info("Using the APY default MTB infection prevalence calibration target.")

        use_default_active_tb = st.checkbox(
            "Use default APY active-TB prevalence",
            value=config.get("activeTBPrevalence") is None,
        )
        active_tb_prevalence_percent = None
        if not use_default_active_tb:
            active_tb_prevalence_percent = st.number_input(
                "Assumed active-TB prevalence calibration target",
                min_value=0.001,
                max_value=99.99,
                value=float(fraction_to_percent(config.get("activeTBPrevalence"), 10 / 770 * 100)),
                step=0.01,
                format="%.3f",
                help="Displayed as a percentage and stored internally as a fraction.",
            )
        else:
            st.info("Using the APY default active-TB prevalence calibration target.")

        risk_prev_updates: dict[str, object] = {}
        with st.expander("Risk-factor prevalence overrides", expanded=False):
            st.caption(
                "Default behaviour uses values from default_data.csv. Custom values are percentages and are stored internally as fractions."
            )
            risk_prev = config.get("riskPrev") or {}
            mode_options = ["Use default", "Single overall", "Three age groups"]
            for factor in PRINCIPAL_RISK_FACTORS:
                st.markdown(f"**{RISK_FACTOR_LABELS[factor]}**")
                current_value = risk_prev.get(factor)
                mode = st.radio(
                    f"{RISK_FACTOR_LABELS[factor]} source",
                    mode_options,
                    index=choice_index(risk_override_mode(current_value), mode_options),
                    horizontal=True,
                    key=f"risk_mode_{factor}",
                )
                current_pcts = risk_override_percent_values(current_value)
                if mode == "Single overall":
                    pct = st.number_input(
                        f"{RISK_FACTOR_LABELS[factor]} prevalence (%)",
                        min_value=0.01,
                        max_value=99.99,
                        value=float(current_pcts[0] or 1.0),
                        step=0.1,
                        format="%.2f",
                        key=f"risk_single_{factor}",
                    )
                    risk_prev_updates[factor] = risk_override_from_percentages(mode, [pct])
                elif mode == "Three age groups":
                    cols = st.columns(3)
                    values = []
                    for idx, age_label in enumerate(AGE_GROUP_LABELS):
                        values.append(
                            cols[idx].number_input(
                                f"{RISK_FACTOR_LABELS[factor]} {age_label} (%)",
                                min_value=0.01,
                                max_value=99.99,
                                value=float(current_pcts[idx] or 1.0),
                                step=0.1,
                                format="%.2f",
                                key=f"risk_age_{factor}_{idx}",
                            )
                        )
                    risk_prev_updates[factor] = risk_override_from_percentages(mode, values)
                else:
                    risk_prev_updates[factor] = None

            with st.expander("Advanced overall prevalence overrides"):
                for factor in ADVANCED_RISK_FACTORS:
                    current_value = risk_prev.get(factor)
                    use_default = st.checkbox(
                        f"Use default {RISK_FACTOR_LABELS[factor]} prevalence",
                        value=current_value is None,
                        key=f"adv_default_{factor}",
                    )
                    if use_default:
                        risk_prev_updates[factor] = None
                    else:
                        pct = st.number_input(
                            f"{RISK_FACTOR_LABELS[factor]} prevalence (%)",
                            min_value=0.01,
                            max_value=99.99,
                            value=float(fraction_to_percent(current_value, 1.0)),
                            step=0.1,
                            format="%.2f",
                            key=f"adv_pct_{factor}",
                        )
                        risk_prev_updates[factor] = risk_override_from_percentages("Single overall", [pct])

        with st.expander("Advanced / run controls"):
            st.caption("Technical run controls for stress testing and reproducibility.")
            population_size = st.number_input(
                "N",
                min_value=1,
                value=int(config.get("N", 1)),
                step=100,
            )
            n_reps = st.number_input(
                "Simulation replicates",
                min_value=1,
                value=int(config.get("nReps", 1)),
                step=10,
            )
            seed = st.number_input(
                "Random seed",
                min_value=0,
                value=int(config.get("seed", 1)),
                step=1,
            )
            screen_window = st.number_input(
                "Screen window",
                min_value=0.01,
                value=float(config.get("screenWindow", 2.0)),
                step=0.5,
            )
            follow_horizon = st.number_input(
                "Follow-up horizon",
                min_value=0.01,
                value=float(config.get("followHorizon", 20.0)),
                step=1.0,
            )
        submitted = st.form_submit_button("Apply edits")

    if submitted:
        updates = {
            "scenarioLabel": scenario_label,
            "N": int(population_size),
            "nReps": int(n_reps),
            "seed": int(seed),
            "screenWindow": float(screen_window),
            "followHorizon": float(follow_horizon),
            "screenCoverage": float(screen_coverage),
            "testType": test_type,
            "regimen": regimen,
            "screeningStrategy": screening_strategy,
        }
        updated_config = dict(config)
        updated_config.update(updates)
        updated_config = apply_epidemiology_updates(
            updated_config,
            use_default_ltbi_prevalence=use_default_ltbi,
            ltbi_prevalence_percent=ltbi_prevalence_percent,
            use_default_active_tb_prevalence=use_default_active_tb,
            active_tb_prevalence_percent=active_tb_prevalence_percent,
            risk_prev_updates=risk_prev_updates,
        )
        changed = updated_config != config
        if changed:
            st.session_state["config"] = updated_config
            mark_config_changed()
            st.success("Scenario edits applied.")
        else:
            st.info("No scenario fields changed.")

    st.subheader("Current Config")
    st.dataframe(
        arrow_safe_dataframe(config_overview_rows(config)),
        width="content",
        hide_index=True,
    )

    with st.expander("Full config"):
        st.json(config, expanded=False)

    if st.button("Validate current config"):
        try:
            st.session_state["validation_report"] = backend.validate_config(config)
            sync_backend_status(backend.status())
            mark_validation_completed()
            st.success("Validation completed.")
        except Exception as exc:
            message = f"Validation failed: {exc}"
            sync_backend_status(backend.status())
            record_message("error", message)
            st.error(message)

    report = st.session_state.get("validation_report")
    if report:
        st.subheader("Validation Report")
        valid = report.get("isValid")
        if valid is True:
            st.success("Config is valid.")
        elif valid is False:
            st.error("Config has validation errors.")

        rows = validation_rows(report)
        if rows:
            st.dataframe(
                arrow_safe_dataframe(rows),
                width="stretch",
                hide_index=True,
            )
        else:
            st.json(report, expanded=False)

    st.subheader("Scenario Files")
    file_name = st.text_input("Scenario file", value="streamlit_scenario.json")
    cols = st.columns(2)
    if cols[0].button("Save scenario"):
        try:
            path = scenario_path(file_name)
            st.session_state["save_info"] = backend.save_scenario(
                config,
                str(path),
                st.session_state.get("economics_config"),
            )
            sync_backend_status(backend.status())
            st.success(f"Saved scenario to {path}")
        except Exception as exc:
            message = f"Save failed: {exc}"
            sync_backend_status(backend.status())
            record_message("error", message)
            st.error(message)

    if cols[1].button("Load scenario"):
        try:
            path = scenario_path(file_name)
            payload = json.loads(path.read_text(encoding="utf-8"))
            loaded_config = payload.get("config")
            if not isinstance(loaded_config, dict):
                raise ValueError("Scenario JSON does not contain a config object.")
            st.session_state["config"] = loaded_config
            st.session_state["economics_config"] = payload.get("economics")
            st.session_state["validation_report"] = None
            st.session_state["load_info"] = {
                "filename": str(path),
                "contractVersion": payload.get("contractVersion", ""),
                "scenarioLabel": payload.get("scenarioLabel", ""),
            }
            st.session_state["results_bundle"] = None
            st.session_state["results_stale"] = False
            st.session_state["dirty_config"] = True
            sync_backend_status(backend.status())
            st.success(f"Loaded scenario from {path}")
            st.rerun()
        except Exception as exc:
            message = f"Load failed: {exc}"
            sync_backend_status(backend.status())
            record_message("error", message)
            st.error(message)
else:
    st.info("Load backend defaults to initialize the scenario state.")
