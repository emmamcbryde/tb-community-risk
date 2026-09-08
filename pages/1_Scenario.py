from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from app.demographic_profile import (
    age_distribution_rows,
    broad_age_distribution_rows,
    demographic_summary_rows,
    restore_apy_demographic_defaults,
    risk_factor_rows,
)
from app.display import arrow_safe_dataframe, display_string
from app.epidemiology_inputs import (
    ADVANCED_RISK_FACTORS,
    AGE_GROUP_LABELS,
    PRINCIPAL_RISK_FACTORS,
    RISK_FACTOR_LABELS,
    apply_epidemiology_updates,
    apply_ltbi_state_assumption_update,
    fraction_to_percent,
    ltbi_state_display_rows,
    risk_override_from_percentages,
    risk_override_mode,
    risk_override_percent_values,
)
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions
from app.state import (
    get_backend,
    init_session_state,
    mark_config_changed,
    mark_validation_completed,
    record_message,
    sync_backend_status,
)
from adapters.paths import scenarios_dir
from engine.apy.scenario import (
    DEFAULT_POPULATION_PRESET_ID,
    build_scenario_contract,
    config_updates_from_scenario,
)


init_session_state()
st.session_state["apy_backend_name"] = "python_apy"

backend = get_backend()

st.title("Define Strategy")
st.caption("Set the population, testing, treatment and prioritisation assumptions for the analysis.")

STRATEGY_LABELS = {
    "prevent": "Prioritise people most likely to avoid active TB",
    "ltbi": "Prioritise people most likely to have LTBI",
    "cure": "Prioritise people most likely to complete effective treatment",
    "random": "No risk-based prioritisation",
}
STRATEGY_VALUES = {label: value for value, label in STRATEGY_LABELS.items()}
LTBI_REVIEW_LABELS = {
    "unresolved": "Unresolved",
    "configured_reviewed": "Reviewed direct evidence",
    "model_derived_reviewed": "Reviewed model-derived estimate",
    "migrated_legacy_unverified": "Migrated legacy value - not reviewed",
    "provisional": "Provisional - not reviewed",
}
LTBI_REVIEW_VALUES = {label: value for value, label in LTBI_REVIEW_LABELS.items()}


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


def _render_demographic_profile(config: dict) -> None:
    st.subheader("Demography currently used by the model")
    st.caption(
        "The analysis uses this resolved age and risk-factor profile when generating "
        "the cohort, assigning LTBI risk and ranking people for targeted screening. "
        "Blank optional override fields mean these source-default values remain in use."
    )
    st.dataframe(
        arrow_safe_dataframe(demographic_summary_rows(config)),
        use_container_width=True,
        hide_index=True,
    )
    with st.expander("Age-group distribution and risk factors", expanded=False):
        st.write("Broad age structure used by the model")
        st.dataframe(
            arrow_safe_dataframe(broad_age_distribution_rows(config)),
            use_container_width=True,
            hide_index=True,
        )
        st.write("Source age-group distribution")
        st.dataframe(
            arrow_safe_dataframe(age_distribution_rows(config)),
            use_container_width=True,
            hide_index=True,
        )
        st.write("Risk-factor prevalences after resolving APY defaults and user overrides")
        st.dataframe(
            arrow_safe_dataframe(risk_factor_rows(config)),
            use_container_width=True,
            hide_index=True,
        )
    if st.button("Restore APY demographic defaults"):
        restored = restore_apy_demographic_defaults(config)
        if restored != config:
            st.session_state["config"] = restored
            st.session_state.pop("recent_ltbi_run_route", None)
            mark_config_changed()
            st.success("APY demographic defaults restored. Rerun epidemiology before using previous results.")
            st.rerun()
        st.info("APY demographic defaults are already loaded.")


def scenario_path(filename: str) -> Path:
    base = scenarios_dir()
    base.mkdir(parents=True, exist_ok=True)
    name = filename.strip() or "streamlit_analysis.json"
    path = Path(name)
    if not path.is_absolute():
        path = base / path
    return path


if st.button("Create a new analysis", type="primary"):
    try:
        st.session_state["config"] = backend.default_config()
        reset_run_state()
        st.session_state.pop("recent_ltbi_run_route", None)
        sync_backend_status(backend.status())
        st.success("New analysis created.")
    except Exception as exc:
        message = f"Could not load defaults: {exc}"
        sync_backend_status(backend.status())
        record_message("error", message)
        st.error(message)

if st.button("Load APY demonstration preset"):
    try:
        base_config = backend.default_config()
        scenario = build_scenario_contract(DEFAULT_POPULATION_PRESET_ID)
        base_config.update(config_updates_from_scenario(scenario))
        st.session_state["config"] = base_config
        reset_run_state()
        st.session_state.pop("recent_ltbi_run_route", None)
        sync_backend_status(backend.status())
        st.success("APY demonstration loaded.")
    except Exception as exc:
        message = f"Could not load APY demonstration preset: {exc}"
        record_message("error", message)
        st.error(message)

config = st.session_state.get("config")
if config:
    if st.session_state.get("dirty_config"):
        st.warning("Inputs have unvalidated changes.")
    if st.session_state.get("results_stale"):
        st.warning("Existing results are stale because inputs changed after the last run.")

    _render_demographic_profile(config)

    st.subheader("Edit Strategy")
    test_options = ["IGRA", "TST"]
    regimen_options = ["3HP", "4R", "3HR", "6H", "9H"]
    strategy_options = list(STRATEGY_VALUES)

    with st.form("scenario_edits"):
        scenario_label = st.text_input(
            "Analysis label",
            value=str(config.get("scenarioLabel", "")),
        )
        population_preset_id = st.text_input(
            "Population preset identifier",
            value=str(config.get("populationPresetId", DEFAULT_POPULATION_PRESET_ID)),
        )
        population_size = st.number_input(
            "Population size",
            min_value=1,
            value=int(config.get("N", 1)),
            step=100,
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
        test_characteristics = config.get("testCharacteristics", {})
        igra_defaults = test_characteristics.get("IGRA", {})
        tst_defaults = test_characteristics.get("TST", {})
        igra_sensitivity = st.number_input(
            "IGRA sensitivity",
            min_value=0.0,
            max_value=1.0,
            value=float(config.get("testSensitivity", igra_defaults.get("sensitivity", 0.95))),
            step=0.01,
        )
        igra_specificity = st.number_input(
            "IGRA specificity",
            min_value=0.0,
            max_value=1.0,
            value=float(config.get("testSpecificity", igra_defaults.get("specificity", 0.98))),
            step=0.01,
        )
        tst_sensitivity = st.number_input(
            "TST sensitivity",
            min_value=0.0,
            max_value=1.0,
            value=float(config.get("tstSensitivity", tst_defaults.get("sensitivity", 0.80))),
            step=0.01,
        )
        tst_specificity_no_bcg = st.number_input(
            "TST specificity without BCG",
            min_value=0.0,
            max_value=1.0,
            value=float(config.get("tstSpecificityNoBCG", tst_defaults.get("specificity", 0.97))),
            step=0.01,
        )
        tst_specificity_bcg = st.number_input(
            "TST specificity with prior BCG",
            min_value=0.0,
            max_value=1.0,
            value=float(config.get("tstSpecificityBCG", tst_defaults.get("specificityBCG", 0.55))),
            step=0.01,
        )
        regimen = st.selectbox(
            "Regimen",
            regimen_options,
            index=choice_index(config.get("regimen"), regimen_options),
        )
        screening_strategy = st.selectbox(
            "Prioritisation approach",
            strategy_options,
            index=choice_index(
                STRATEGY_LABELS.get(str(config.get("screeningStrategy")), ""),
                strategy_options,
                0,
            ),
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

        ltbi_state = resolve_ltbi_state_assumptions(config)
        with st.expander("Advanced LTBI-state assumptions", expanded=False):
            st.dataframe(
                arrow_safe_dataframe(ltbi_state_display_rows(config)),
                use_container_width=True,
                hide_index=True,
            )
            if ltbi_state.get("warning"):
                st.warning(str(ltbi_state["warning"]))
            ltbi_unresolved = st.checkbox(
                "Baseline recent-LTBI proportion unresolved",
                value=ltbi_state.get("baselineRecentLTBIProportion") is None
                or str(ltbi_state.get("status", "")).startswith("unresolved"),
            )
            ltbi_recent_percent = None
            if not ltbi_unresolved:
                ltbi_recent_percent = st.number_input(
                    "Baseline recent-LTBI proportion (%)",
                    min_value=0.0,
                    max_value=100.0,
                    value=float(
                        fraction_to_percent(
                            ltbi_state.get("baselineRecentLTBIProportion"),
                            0.0,
                        )
                    ),
                    step=1.0,
                    format="%.2f",
                )
            ltbi_transition_rate = st.number_input(
                "Recent-to-remote transition rate (per year)",
                min_value=0.0001,
                value=float(ltbi_state["recentToRemoteTransitionRatePerYear"]),
                step=0.01,
                format="%.4f",
            )
            st.caption(
                f"Implied mean residence time: {1.0 / float(ltbi_transition_rate):.2f} years. "
                "This is a Markov transition, not a deterministic five-year cutoff."
            )
            ltbi_source = st.text_input(
                "LTBI-state assumption source",
                value=str(ltbi_state.get("source") or ""),
            )
            ltbi_status = st.selectbox(
                "LTBI-state assumption review",
                list(LTBI_REVIEW_VALUES),
                index=choice_index(
                    LTBI_REVIEW_LABELS.get(str(ltbi_state.get("status")), "Unresolved"),
                    list(LTBI_REVIEW_VALUES),
                ),
            )
            ltbi_notes = st.text_area(
                "LTBI-state notes",
                value=str(ltbi_state.get("notes") or ""),
            )

        risk_prev_updates: dict[str, object] = {}
        with st.expander("Risk-factor prevalence overrides", expanded=False):
            st.caption(
                "Blank or default risk-factor override fields mean use source defaults from default_data.csv. "
                "Blank fields are optional overrides, not missing model inputs. "
                "Custom values are percentages and are stored internally as fractions."
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
            "populationPresetId": population_preset_id,
            "N": int(population_size),
            "nReps": int(n_reps),
            "seed": int(seed),
            "screenWindow": float(screen_window),
            "followHorizon": float(follow_horizon),
            "screenCoverage": float(screen_coverage),
            "testType": test_type,
            "testSensitivity": float(igra_sensitivity),
            "testSpecificity": float(igra_specificity),
            "tstSensitivity": float(tst_sensitivity),
            "tstSpecificityNoBCG": float(tst_specificity_no_bcg),
            "tstSpecificityBCG": float(tst_specificity_bcg),
            "regimen": regimen,
            "screeningStrategy": STRATEGY_VALUES[screening_strategy],
        }
        updated_config = dict(config)
        updated_config.update(updates)
        updated_config.setdefault("testCharacteristics", {})
        updated_config["testCharacteristics"].setdefault("IGRA", {})
        updated_config["testCharacteristics"].setdefault("TST", {})
        updated_config["testCharacteristics"]["IGRA"].update(
            {"sensitivity": float(igra_sensitivity), "specificity": float(igra_specificity)}
        )
        updated_config["testCharacteristics"]["TST"].update(
            {
                "sensitivity": float(tst_sensitivity),
                "specificity": float(tst_specificity_no_bcg),
                "specificityBCG": float(tst_specificity_bcg),
            }
        )
        updated_config = apply_epidemiology_updates(
            updated_config,
            use_default_ltbi_prevalence=use_default_ltbi,
            ltbi_prevalence_percent=ltbi_prevalence_percent,
            use_default_active_tb_prevalence=use_default_active_tb,
            active_tb_prevalence_percent=active_tb_prevalence_percent,
            risk_prev_updates=risk_prev_updates,
        )
        updated_config = apply_ltbi_state_assumption_update(
            updated_config,
            baseline_recent_percent=ltbi_recent_percent,
            transition_rate_per_year=float(ltbi_transition_rate),
            source=ltbi_source,
            status="unresolved" if ltbi_unresolved else LTBI_REVIEW_VALUES[ltbi_status],
            notes=ltbi_notes,
        )
        changed = updated_config != config
        if changed:
            st.session_state["config"] = updated_config
            st.session_state.pop("recent_ltbi_run_route", None)
            mark_config_changed()
            st.success("Strategy edits applied.")
        else:
            st.info("No strategy fields changed.")

    with st.expander("Technical information", expanded=False):
        st.dataframe(
            arrow_safe_dataframe(config_overview_rows(config)),
            use_container_width=True,
            hide_index=True,
        )
        st.json(config, expanded=False)

    if st.button("Validate inputs"):
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
        st.subheader("Input Checks")
        valid = report.get("isValid")
        if valid is True:
            st.success("Inputs are valid.")
        elif valid is False:
            st.error("Inputs have validation errors.")

        rows = validation_rows(report)
        if rows:
            st.dataframe(
                arrow_safe_dataframe(rows),
                use_container_width=True,
                hide_index=True,
            )
        else:
            st.json(report, expanded=False)

    st.subheader("Analysis Files")
    file_name = st.text_input("Analysis file", value="streamlit_analysis.json")
    cols = st.columns(2)
    if cols[0].button("Save analysis"):
        try:
            path = scenario_path(file_name)
            st.session_state["save_info"] = backend.save_scenario(
                config,
                str(path),
                st.session_state.get("economics_config"),
            )
            sync_backend_status(backend.status())
            st.success(f"Saved analysis to {path}")
        except Exception as exc:
            message = f"Save failed: {exc}"
            sync_backend_status(backend.status())
            record_message("error", message)
            st.error(message)

    if cols[1].button("Load analysis"):
        try:
            path = scenario_path(file_name)
            payload = json.loads(path.read_text(encoding="utf-8"))
            loaded_config = payload.get("config")
            if not isinstance(loaded_config, dict):
                raise ValueError("Analysis JSON does not contain an input object.")
            st.session_state["config"] = loaded_config
            st.session_state["economics_config"] = payload.get("economics")
            st.session_state.pop("recent_ltbi_run_route", None)
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
            st.success(f"Loaded analysis from {path}")
            st.rerun()
        except Exception as exc:
            message = f"Load failed: {exc}"
            sync_backend_status(backend.status())
            record_message("error", message)
            st.error(message)
else:
    st.info("Create a new analysis or load the APY demonstration to begin.")
