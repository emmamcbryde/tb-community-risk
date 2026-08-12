from __future__ import annotations

from copy import deepcopy
import json
from typing import Any

from engine.apy.economics import sync_legacy_cost_fields_from_cost_items
from engine.apy.working_defaults import (
    HIGHER_BURDEN_DALY_OUTCOME_PRESET,
    PRIMARY_DALY_OUTCOME_PRESET,
    apply_outcome_preset,
    apply_program_setup_preset,
    build_unified_working_default_preset,
)


PARAMETER_WORKSPACE_VERSION = "apy_parameter_workspace_v1"
PARAMETER_GROUPS = [
    "Demography",
    "TB epidemiology",
    "Health-service provision and intervention",
    "Costs and outcomes",
    "Analysis settings",
]
MODEL_METHOD_LABELS = {
    "expected_value": "Expected outcomes",
    "agent_based": "Simulated community variation",
}
MODEL_METHOD_CODES = {value: key for key, value in MODEL_METHOD_LABELS.items()}
SIMULATION_MODE_LABELS = {
    "quick_preview": "Quick preview: 5 simulations",
    "standard": "Standard analysis: use the current validated repository standard",
    "custom": "Custom",
}
PROGRAMME_SETUP_PRESET_LABELS = {
    "existing_service_base": "Existing-service/base scenario",
    "new_programme_implementation": "New-programme implementation scenario",
}
TARGETING_LABELS = {
    "prevent": "Prioritise people most likely to avoid active TB",
    "ltbi": "Prioritise people most likely to have LTBI",
    "cure": "Prioritise people most likely to complete effective treatment",
    "random": "No risk-based prioritisation",
}
TARGETING_CODES = {value: key for key, value in TARGETING_LABELS.items()}


def unified_default_session_state() -> dict[str, Any]:
    preset = build_unified_working_default_preset()
    workspace = build_parameter_workspace(preset["config"], preset["economicsConfig"], preset=preset)
    return {
        "config": preset["config"],
        "economics_config": preset["economicsConfig"],
        "parameter_workspace": workspace,
        "working_default_preset": {
            key: preset[key]
            for key in [
                "contractVersion",
                "presetId",
                "presetVersion",
                "label",
                "sourceComponentPresets",
                "workingDefault",
                "referenceStatus",
                "configurationHash",
                "provisionalAssumptions",
                "unresolvedAssumptions",
            ]
        },
    }


def build_parameter_workspace(
    config: dict[str, Any],
    economics_config: dict[str, Any],
    *,
    preset: dict[str, Any] | None = None,
) -> dict[str, Any]:
    preset = preset or build_unified_working_default_preset()
    default_config = preset["config"]
    default_econ = preset["economicsConfig"]
    rows = [_parameter_row(spec, config, economics_config, default_config, default_econ) for spec in _parameter_specs()]
    return {
        "contractVersion": PARAMETER_WORKSPACE_VERSION,
        "presetId": preset["presetId"],
        "presetVersion": preset["presetVersion"],
        "presetLabel": preset["label"],
        "configurationHash": preset["configurationHash"],
        "groups": PARAMETER_GROUPS,
        "rows": rows,
        "changedCount": sum(1 for row in rows if row["changedFromDefault"]),
        "validation": validate_parameter_workspace(rows),
    }


def validate_parameter_workspace(rows: list[dict[str, Any]]) -> dict[str, Any]:
    messages = []
    method_value = next(
        (
            row.get("currentValue")
            for row in rows
            if row.get("parameterId") == "analysis.method"
        ),
        None,
    )
    for row in rows:
        value = row.get("currentValue")
        parameter_id = row.get("parameterId")
        kind = row.get("editableType")
        if parameter_id == "analysis.n_reps" and method_value in {"Expected outcomes", "expected_value"}:
            continue
        if kind == "select" and row.get("selectOptions") and value not in row.get("selectOptions"):
            messages.append({"parameterId": parameter_id, "message": "Select one of the listed options."})
        if kind == "probability" and value not in (None, ""):
            number = _number_or_none(value)
            if number is None or number < 0 or number > 1:
                messages.append({"parameterId": parameter_id, "message": "Probability must be in [0, 1]."})
        elif kind in {"positive_integer", "nonnegative_number", "years", "money"} and value not in (None, ""):
            number = _number_or_none(value)
            if number is None:
                messages.append({"parameterId": parameter_id, "message": "Value must be numeric."})
            elif kind == "positive_integer" and int(number) <= 0:
                messages.append({"parameterId": parameter_id, "message": "Value must be positive."})
            elif kind in {"nonnegative_number", "years", "money"} and number < 0:
                messages.append({"parameterId": parameter_id, "message": "Value must be non-negative."})
    age_rows = [row for row in rows if str(row.get("parameterId", "")).startswith("demography.age.")]
    age_values = [_number_or_none(row.get("currentValue")) for row in age_rows]
    if any(value is not None for value in age_values):
        if any(value is None for value in age_values):
            messages.append({"parameterId": "demography.age", "message": "All age-distribution bands are required when overriding age structure."})
        elif abs(sum(float(value) for value in age_values) - 1.0) > 1e-6:
            messages.append({"parameterId": "demography.age", "message": "Age-distribution probabilities must sum to 1."})
    return {"isValid": not messages, "messages": messages}


def apply_parameter_workspace(
    config: dict[str, Any],
    economics_config: dict[str, Any],
    rows: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    validation = validate_parameter_workspace(rows)
    if not validation["isValid"]:
        raise ValueError("Parameter workspace contains validation errors.")
    cfg = deepcopy(config)
    econ = deepcopy(economics_config)
    for row in rows:
        value = row.get("currentValue")
        parameter_id = row.get("parameterId")
        if parameter_id == "analysis.method":
            cfg["analysisMethod"] = MODEL_METHOD_CODES.get(str(value), str(value))
            cfg["analysisMethodLabel"] = MODEL_METHOD_LABELS.get(cfg["analysisMethod"], str(value))
        elif parameter_id == "analysis.simulation_mode":
            _apply_simulation_mode(cfg, str(value))
        elif parameter_id == "analysis.n_reps" and cfg.get("simulationMode") == "quick_preview":
            continue
        elif parameter_id == "cost.program_setup_preset":
            code = _programme_setup_code(str(value))
            econ = apply_program_setup_preset(econ, code)
        elif parameter_id == "outcome.daly_preset":
            econ = apply_outcome_preset(econ, str(value))
        else:
            spec = next((item for item in _parameter_specs() if item["parameterId"] == parameter_id), None)
            if spec:
                _set_value(cfg, econ, spec, value)
    econ = sync_legacy_cost_fields_from_cost_items(econ)
    return cfg, econ


def reset_parameter_group(rows: list[dict[str, Any]], group: str) -> list[dict[str, Any]]:
    out = deepcopy(rows)
    for row in out:
        if row.get("group") == group:
            row["currentValue"] = row.get("defaultValue")
            row["changedFromDefault"] = False
    return out


def reset_all_parameters(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = deepcopy(rows)
    for row in out:
        row["currentValue"] = row.get("defaultValue")
        row["changedFromDefault"] = False
    return out


def parameter_summary(config: dict[str, Any], economics_config: dict[str, Any]) -> list[dict[str, Any]]:
    metadata = economics_config.get("metadata") or {}
    method = MODEL_METHOD_LABELS.get(str(config.get("analysisMethod") or "expected_value"), "Expected outcomes")
    discounting = economics_config.get("discounting") or {}
    profiles = discounting.get("profiles") or {}
    primary = profiles.get("primary") or {}
    comparison = profiles.get("comparison") or {}
    primary_cost = primary.get("costRate", discounting.get("primaryCostRate", discounting.get("primaryDisplayedRate", 0.03)))
    primary_health = primary.get("healthRate", discounting.get("primaryHealthRate", discounting.get("healthOutcomeAnnualRate", primary_cost)))
    comparison_cost = comparison.get("costRate", discounting.get("comparisonCostRate", discounting.get("comparisonRate", 0.0)))
    comparison_health = comparison.get("healthRate", discounting.get("comparisonHealthRate", comparison_cost))
    return [
        {"Item": "Population", "Value": config.get("populationPresetId")},
        {"Item": "Selected test", "Value": config.get("testType")},
        {"Item": "Selected regimen", "Value": config.get("regimen")},
        {"Item": "Screening coverage", "Value": config.get("screenCoverage")},
        {"Item": "Intervention duration", "Value": f"{config.get('screeningWindowYears')} years"},
        {"Item": "Follow-up horizon", "Value": f"{config.get('followUpHorizonYears')} years"},
        {"Item": "Cost currency/year", "Value": f"{metadata.get('targetCurrency')} {metadata.get('targetPriceYear')}"},
        {
            "Item": "Discount rates",
            "Value": (
                f"Primary cost {_percent(primary_cost)}, primary health {_percent(primary_health)}; "
                f"comparison cost {_percent(comparison_cost)}, comparison health {_percent(comparison_health)}"
            ),
        },
        {"Item": "Health outcome", "Value": "DALYs"},
        {"Item": "Analysis method", "Value": method},
        {"Item": "Number of simulations", "Value": "Not applicable" if method == "Expected outcomes" else config.get("nReps")},
    ]


def changed_parameter_count(workspace: dict[str, Any] | None) -> int:
    return sum(1 for row in (workspace or {}).get("rows") or [] if row.get("changedFromDefault"))


def _parameter_row(
    spec: dict[str, Any],
    config: dict[str, Any],
    econ: dict[str, Any],
    default_config: dict[str, Any],
    default_econ: dict[str, Any],
) -> dict[str, Any]:
    current = _get_value(config, econ, spec)
    default = _get_value(default_config, default_econ, spec)
    return {
        **spec,
        "currentValue": current,
        "defaultValue": default,
        "changedFromDefault": _normalise_compare(current) != _normalise_compare(default),
    }


def _parameter_specs() -> list[dict[str, Any]]:
    return [
        _spec("demography.population_preset", "Demography", "Population preset", "config", ["populationPresetId"], "", "APY demonstration preset", "unreviewed_repository_input", True, "text", operational_status="authoritative_model_input"),
        _spec("demography.population_size", "Demography", "Population size", "config", ["N"], "people", "APY working default", "configured_reviewed", False, "positive_integer"),
        _spec("demography.eligible_population", "Demography", "Eligible population", "config", ["scenario", "eligible", "number"], "people", "Scenario contract", "unreviewed_repository_input", True, "nonnegative_number"),
        _spec("demography.age.0_4", "Demography", "Age distribution: 0-4 years", "ageDistributionBand", ["0-4"], "probability", "APY age-distribution input", "unreviewed_repository_input", True, "probability", advanced=True, operational_status="authoritative_model_input"),
        _spec("demography.age.5_14", "Demography", "Age distribution: 5-14 years", "ageDistributionBand", ["5-14"], "probability", "APY age-distribution input", "unreviewed_repository_input", True, "probability", advanced=True, operational_status="authoritative_model_input"),
        _spec("demography.age.15_plus", "Demography", "Age distribution: 15+ years", "ageDistributionBand", ["15+"], "probability", "APY age-distribution input", "unreviewed_repository_input", True, "probability", advanced=True, operational_status="authoritative_model_input"),
        *[
            _spec(f"demography.risk.{key}", "Demography", label, "config", ["riskPrev", key], "prevalence", "APY data inputs", "unreviewed_repository_input", True, "probability")
            for key, label in [
                ("female", "Sex distribution - female"),
                ("contact", "Contact history"),
                ("smoking", "Smoking"),
                ("diabetes", "Diabetes"),
                ("renal", "Renal impairment"),
                ("cld", "Chronic lung disease"),
                ("alcohol", "Alcohol/drug exposure"),
                ("MJ", "Other current risk-factor prevalence"),
            ]
        ],
        _spec("tb.ltbi_prevalence", "TB epidemiology", "Baseline LTBI prevalence", "config", ["ltbiPrevalence"], "probability", "APY calibration input", "unreviewed_repository_input", True, "probability"),
        _spec("tb.age_pattern", "TB epidemiology", "Age pattern of infection", "config", ["targetAgeOR"], "odds ratio", "APY calibration input", "unreviewed_repository_input", True, "nonnegative_number"),
        _spec("tb.recent_fraction", "TB epidemiology", "Baseline recent-LTBI proportion", "config", ["ltbiStateAssumptions", "baselineRecentLTBIProportion"], "probability", "No APY-specific source established", "unresolved", True, "probability"),
        _spec("tb.transition_model", "TB epidemiology", "Recent-to-remote transition model", "config", ["ltbiStateAssumptions", "transitionModel"], "", "Inherited model structure", "configured_reviewed", False, "select"),
        _spec("tb.transition_rate", "TB epidemiology", "Recent-to-remote transition rate", "config", ["ltbiStateAssumptions", "recentToRemoteTransitionRatePerYear"], "per year", "Inherited model structure", "configured_reviewed", False, "nonnegative_number"),
        _spec("tb.active_tb_target", "TB epidemiology", "Active-TB calibration target", "config", ["targetActive2y"], "probability", "APY calibration input", "unreviewed_repository_input", True, "probability", advanced=True),
        _spec("tb.active_tb_horizon", "TB epidemiology", "Active-TB calibration horizon", "config", ["activeTBCalibrationHorizonYears"], "years", "APY calibration setting", "configured_reviewed", False, "years", advanced=True),
        _spec("tb.direct_scope", "TB epidemiology", "Direct-effects-only scope statement", "static", ["directScope"], "", "Model scope", "configured_reviewed", False, "read_only", operational_status="descriptive_metadata"),
        _spec("service.comparator", "Health-service provision and intervention", "Comparator", "static", ["comparator"], "", "Scenario contract", "configured_reviewed", False, "read_only", operational_status="descriptive_metadata"),
        _spec("service.targeting", "Health-service provision and intervention", "Targeting method", "config", ["screeningStrategy"], "", "APY intervention setting", "unreviewed_repository_input", True, "select"),
        _spec("service.coverage", "Health-service provision and intervention", "Screening coverage", "config", ["screenCoverage"], "proportion eligible", "APY intervention setting", "configured_reviewed", False, "probability"),
        _spec("service.test", "Health-service provision and intervention", "Screening test", "config", ["testType"], "", "APY intervention setting", "configured_reviewed", False, "select"),
        _spec("service.igra_sensitivity", "Health-service provision and intervention", "IGRA sensitivity", "config", ["testSensitivity"], "probability", "APY test assumption", "unreviewed_repository_input", True, "probability"),
        _spec("service.igra_specificity", "Health-service provision and intervention", "IGRA specificity", "config", ["testSpecificity"], "probability", "APY test assumption", "unreviewed_repository_input", True, "probability"),
        _spec("service.tst_sensitivity", "Health-service provision and intervention", "TST sensitivity", "config", ["tstSensitivity"], "probability", "APY test assumption", "unreviewed_repository_input", True, "probability"),
        _spec("service.tst_specificity_bcg", "Health-service provision and intervention", "TST specificity with BCG", "config", ["tstSpecificityBCG"], "probability", "APY test assumption", "unreviewed_repository_input", True, "probability"),
        _spec("service.tst_specificity_no_bcg", "Health-service provision and intervention", "TST specificity without BCG", "config", ["tstSpecificityNoBCG"], "probability", "APY test assumption", "unreviewed_repository_input", True, "probability"),
        _spec("service.regimen", "Health-service provision and intervention", "Preventive regimen", "config", ["regimen"], "", "APY regimen library", "configured_reviewed", False, "select"),
        _spec("service.tpt_start", "Health-service provision and intervention", "Treatment initiation", "config", ["pStartTPT"], "probability", "APY treatment cascade", "unreviewed_repository_input", True, "probability"),
        _spec("service.completion", "Health-service provision and intervention", "Treatment completion", "config", ["regimenPComplete"], "probability", "APY regimen library", "unreviewed_repository_input", True, "probability"),
        _spec("service.adr_stop", "Health-service provision and intervention", "ADR stopping", "config", ["regimenADRstop"], "probability", "APY regimen library", "unreviewed_repository_input", True, "probability"),
        _spec("service.other_stop", "Health-service provision and intervention", "Other stopping", "config", ["regimenOtherStop"], "probability", "APY regimen library", "unreviewed_repository_input", True, "probability"),
        _spec("service.full_efficacy", "Health-service provision and intervention", "Full-course efficacy", "config", ["regimenEffFull"], "probability", "APY regimen library", "unreviewed_repository_input", True, "probability"),
        _spec("service.partial_efficacy", "Health-service provision and intervention", "Partial-course efficacy", "config", ["partialShortCourseMode"], "", "APY regimen library", "unreviewed_repository_input", True, "text"),
        _spec("service.duration", "Health-service provision and intervention", "Intervention duration", "config", ["screeningWindowYears"], "years", "APY working default", "configured_reviewed", False, "years"),
        *[_cost_spec(item_id, label) for item_id, label in [
            ("test_igra", "IGRA cost"),
            ("test_tst", "TST cost"),
            ("regimen_3hp", "3HP regimen cost"),
            ("regimen_4r", "4R regimen cost"),
            ("regimen_3hr", "3HR regimen cost"),
            ("regimen_6h", "6H regimen cost"),
            ("regimen_9h", "9H regimen cost"),
            ("active_tb_disease", "Active-TB disease care cost"),
            ("tpt_adr_management", "ADR-management cost"),
        ]],
        _cost_spec("return_for_results", "Return-for-results cost"),
        _cost_spec("clinical_review", "Clinical-review cost"),
        _cost_spec("active_tb_exclusion_workup", "Active-TB exclusion / CXR / laboratory cost"),
        _cost_spec("travel_outreach_staff_support", "Travel/outreach/staff-support cost"),
        _spec("cost.program_setup_preset", "Costs and outcomes", "Programme setup preset", "econ", ["localSAHealthPathwayCosts", "selectedProgramSetupPreset"], "", "SA Health local working assumption", "configured_reviewed", False, "select"),
        _spec("cost.program_running", "Costs and outcomes", "Programme running cost", "econ", ["localSAHealthPathwayCosts", "programRunning", "value"], "AUD", "Not yet locally costed", "local_working_assumption", True, "money"),
        _spec("outcome.daly_preset", "Costs and outcomes", "Outcome preset", "econ", ["dalyAssumptions", "outcomePreset"], "", "APY outcome working defaults", "configured_reviewed", False, "select"),
        _daly_spec("activeTBDisabilityWeight", "Active-TB disability weight", "disability weight"),
        _daly_spec("activeTBDurationYears", "Active-TB disability duration", "years"),
        _daly_spec("tbCaseFatalityRisk", "TB case-fatality risk", "probability"),
        _daly_spec("yllPerTBDeath", "Years of life lost per TB death", "years"),
        _spec("analysis.intervention_duration", "Analysis settings", "Intervention duration", "config", ["screeningWindowYears"], "years", "APY working default", "configured_reviewed", False, "years"),
        _spec("analysis.followup", "Analysis settings", "Epidemiological follow-up horizon", "config", ["followUpHorizonYears"], "years", "APY working default", "configured_reviewed", False, "years"),
        _spec("analysis.economic_horizon", "Analysis settings", "Economic horizon", "econ", ["metadata", "economicHorizonYears"], "years", "APY working default", "configured_reviewed", False, "years"),
        _spec("analysis.cost_discount", "Analysis settings", "Cost discount rate", "econ", ["discounting", "primaryDisplayedRate"], "annual rate", "APY working default", "configured_reviewed", False, "probability"),
        _spec("analysis.health_discount", "Analysis settings", "Health-outcome discount rate", "econ", ["discounting", "healthOutcomeAnnualRate"], "annual rate", "APY working default", "configured_reviewed", False, "probability"),
        _spec("analysis.comparison_discount", "Analysis settings", "Comparison discount rate", "econ", ["discounting", "comparisonRate"], "annual rate", "APY working default", "configured_reviewed", False, "probability"),
        _spec("analysis.method", "Analysis settings", "Analysis method", "config", ["analysisMethod"], "", "Workflow setting", "configured_reviewed", False, "select"),
        _spec("analysis.simulation_mode", "Analysis settings", "Simulation count option", "config", ["simulationMode"], "", "Workflow setting", "configured_reviewed", False, "select"),
        _spec("analysis.n_reps", "Analysis settings", "Number of simulations", "config", ["nReps"], "simulations", "Repository standard; five simulations is preview only", "configured_reviewed", False, "positive_integer"),
        _spec("analysis.seed", "Analysis settings", "Random seed", "config", ["seed"], "", "Workflow setting", "configured_reviewed", False, "positive_integer"),
    ]


def _spec(
    parameter_id: str,
    group: str,
    label: str,
    source_object: str,
    path: list[str],
    unit: str,
    source: str,
    status: str,
    provisional: bool,
    editable_type: str,
    *,
    advanced: bool = False,
    operational_status: str | None = None,
) -> dict[str, Any]:
    select_options = {
        "service.targeting": list(TARGETING_LABELS.values()),
        "service.test": ["IGRA", "TST"],
        "service.regimen": ["3HP", "4R", "3HR", "6H", "9H"],
        "analysis.method": list(MODEL_METHOD_LABELS.values()),
        "analysis.simulation_mode": list(SIMULATION_MODE_LABELS.values()),
        "cost.program_setup_preset": list(PROGRAMME_SETUP_PRESET_LABELS.values()),
        "outcome.daly_preset": [PRIMARY_DALY_OUTCOME_PRESET, HIGHER_BURDEN_DALY_OUTCOME_PRESET],
        "tb.transition_model": ["continuous_markov_recent_remote"],
    }.get(parameter_id, [])
    if operational_status is None:
        operational_status = {
            "config": "authoritative_model_input",
            "econ": "authoritative_economic_input",
            "costItem": "authoritative_economic_input",
            "dalyAssumption": "authoritative_economic_input",
            "static": "descriptive_metadata",
        }.get(source_object, "descriptive_metadata")
    return {
        "parameterId": parameter_id,
        "group": group,
        "label": label,
        "sourceObject": source_object,
        "path": path,
        "unit": unit,
        "source": source,
        "reviewStatus": status,
        "provisional": provisional,
        "editableType": editable_type,
        "validation": "",
        "notes": "",
        "advanced": advanced,
        "selectOptions": select_options,
        "operationalStatus": operational_status,
    }


def _cost_spec(item_id: str, label: str) -> dict[str, Any]:
    return _spec(
        f"cost.{item_id}",
        "Costs and outcomes",
        label,
        "costItem",
        [item_id, "originalCost"],
        "AUD 2019",
        "Dale 2019 AUD working defaults",
        "configured_reviewed",
        False,
        "money",
    )


def _daly_spec(field: str, label: str, unit: str) -> dict[str, Any]:
    return _spec(
        f"outcome.{field}",
        "Costs and outcomes",
        label,
        "dalyAssumption",
        [field, "value"],
        unit,
        "Dale/APY working outcome assumptions",
        "configured_reviewed",
        False,
        "nonnegative_number" if field not in {"activeTBDisabilityWeight", "tbCaseFatalityRisk"} else "probability",
    )


def _get_value(config: dict[str, Any], econ: dict[str, Any], spec: dict[str, Any]) -> Any:
    source = spec["sourceObject"]
    if source == "config":
        value = _nested_get(config, spec["path"])
        if spec["parameterId"] == "analysis.method":
            return MODEL_METHOD_LABELS.get(str(value), value)
        if spec["parameterId"] == "analysis.simulation_mode":
            return SIMULATION_MODE_LABELS.get(str(value), value)
        if spec["parameterId"] == "service.targeting":
            return TARGETING_LABELS.get(str(value), value)
        return value
    if source == "econ":
        discount_value = _discount_profile_value(econ, spec["parameterId"])
        if discount_value is not None:
            return discount_value
        value = _nested_get(econ, spec["path"])
        if spec["parameterId"] == "cost.program_setup_preset":
            return PROGRAMME_SETUP_PRESET_LABELS.get(str(value), value)
        return value
    if source == "costItem":
        item = _cost_item(econ, spec["path"][0])
        return item.get(spec["path"][1]) if item else None
    if source == "dalyAssumption":
        rec = (econ.get("dalyAssumptions") or {}).get(spec["path"][0])
        return (rec or {}).get("value") if isinstance(rec, dict) else rec
    if source == "ageDistributionBand":
        return _age_band_value(config, spec["path"][0])
    if source == "static" and spec["path"] == ["directScope"]:
        return "Direct benefits, harms and costs only; transmission-mediated benefits are not included."
    if source == "static" and spec["path"] == ["comparator"]:
        return "Current practice / no additional systematic LTBI screening"
    return None


def _set_value(config: dict[str, Any], econ: dict[str, Any], spec: dict[str, Any], value: Any) -> None:
    source = spec["sourceObject"]
    converted = _coerce_value(value, spec["editableType"])
    if source == "config":
        if spec["parameterId"] == "service.targeting":
            converted = TARGETING_CODES.get(str(value), str(value))
        _nested_set(config, spec["path"], converted)
        if spec["path"] == ["screeningWindowYears"]:
            config["screenWindow"] = converted
        if spec["path"] == ["followUpHorizonYears"]:
            config["followHorizon"] = converted
    elif source == "econ":
        _nested_set(econ, spec["path"], converted)
        _sync_discount_profile_value(econ, spec["parameterId"], converted)
    elif source == "costItem":
        item = _cost_item(econ, spec["path"][0])
        if item:
            item[spec["path"][1]] = converted
            item["conversionStatus"] = "not_converted"
            item["conversionApplied"] = False
            item["costRecordType"] = "source"
            item["warnings"] = []
            _sync_local_pathway_value(econ, spec["path"][0], converted)
    elif source == "dalyAssumption":
        rec = (econ.setdefault("dalyAssumptions", {}).setdefault(spec["path"][0], {}))
        if isinstance(rec, dict):
            rec["value"] = converted
    elif source == "ageDistributionBand":
        _set_age_band_value(config, spec["path"][0], converted)


def _apply_simulation_mode(config: dict[str, Any], value: str) -> None:
    code = next((key for key, label in SIMULATION_MODE_LABELS.items() if label == value), value)
    config["simulationMode"] = code
    config["simulationModeLabel"] = SIMULATION_MODE_LABELS.get(code, value)
    if code == "quick_preview":
        config["nReps"] = 5


def _programme_setup_code(label_or_code: str) -> str:
    for code, label in PROGRAMME_SETUP_PRESET_LABELS.items():
        if label_or_code == label:
            return code
    return label_or_code


def _cost_item(econ: dict[str, Any], item_id: str) -> dict[str, Any] | None:
    for item in econ.get("costItems") or []:
        if item.get("costItemId") == item_id:
            return item
    return None


def _sync_local_pathway_value(econ: dict[str, Any], cost_item_id: str, value: Any) -> None:
    mapping = {
        "return_for_results": "returnForResults",
        "clinical_review": "clinicalReview",
        "active_tb_exclusion_workup": "activeTBExclusionWorkup",
        "travel_outreach_staff_support": "travelOutreachStaffSupport",
    }
    key = mapping.get(cost_item_id)
    if key:
        econ.setdefault("localSAHealthPathwayCosts", {}).setdefault(key, {})["value"] = value


def _discount_profile_value(econ: dict[str, Any], parameter_id: str) -> Any:
    profiles = (econ.get("discounting") or {}).get("profiles") or {}
    if parameter_id == "analysis.cost_discount":
        return ((profiles.get("primary") or {}).get("costRate"))
    if parameter_id == "analysis.health_discount":
        return ((profiles.get("primary") or {}).get("healthRate"))
    if parameter_id == "analysis.comparison_discount":
        comparison = profiles.get("comparison") or {}
        cost = comparison.get("costRate")
        health = comparison.get("healthRate")
        return cost if cost == health or health is None else cost
    return None


def _sync_discount_profile_value(econ: dict[str, Any], parameter_id: str, value: Any) -> None:
    if parameter_id not in {"analysis.cost_discount", "analysis.health_discount", "analysis.comparison_discount"}:
        return
    discounting = econ.setdefault("discounting", {})
    profiles = discounting.setdefault("profiles", {})
    primary = profiles.setdefault("primary", {})
    comparison = profiles.setdefault("comparison", {})
    if parameter_id == "analysis.cost_discount":
        primary["costRate"] = value
        discounting["primaryCostRate"] = value
        discounting["primaryDisplayedRate"] = value
    elif parameter_id == "analysis.health_discount":
        primary["healthRate"] = value
        discounting["primaryHealthRate"] = value
        discounting["healthOutcomeAnnualRate"] = value
    elif parameter_id == "analysis.comparison_discount":
        comparison["costRate"] = value
        comparison["healthRate"] = value
        discounting["comparisonCostRate"] = value
        discounting["comparisonHealthRate"] = value
        discounting["comparisonRate"] = value


def _age_band_value(config: dict[str, Any], band: str) -> Any:
    for row in config.get("ageDistributionTable") or []:
        if str(row.get("age") or row.get("Age") or "").strip() == band:
            return row.get("proportion", row.get("Proportion"))
    return None


def _set_age_band_value(config: dict[str, Any], band: str, value: Any) -> None:
    table = list(config.get("ageDistributionTable") or [])
    by_band = {
        str(row.get("age") or row.get("Age") or "").strip(): dict(row)
        for row in table
        if str(row.get("age") or row.get("Age") or "").strip()
    }
    by_band[band] = {"age": band, "proportion": value}
    ordered = ["0-4", "5-14", "15+"]
    config["ageDistributionTable"] = [by_band[item] for item in ordered if item in by_band]
    config["ageDistributionFile"] = ""


def _nested_get(root: dict[str, Any], path: list[str]) -> Any:
    current: Any = root
    for key in path:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def _nested_set(root: dict[str, Any], path: list[str], value: Any) -> None:
    current = root
    for key in path[:-1]:
        current = current.setdefault(key, {})
    current[path[-1]] = value


def _coerce_value(value: Any, editable_type: str) -> Any:
    if editable_type in {"positive_integer"}:
        return int(float(value))
    if editable_type in {"probability", "nonnegative_number", "years", "money"}:
        if value in (None, ""):
            return None
        return float(value)
    if editable_type == "select" and value in MODEL_METHOD_CODES:
        return MODEL_METHOD_CODES[value]
    return value


def _number_or_none(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _normalise_compare(value: Any) -> str:
    return json.dumps(value, sort_keys=True, default=str)


def _percent(value: Any) -> str:
    try:
        return f"{float(value) * 100:g}%"
    except (TypeError, ValueError):
        return "not set"
