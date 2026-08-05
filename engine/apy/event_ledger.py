from __future__ import annotations

from copy import deepcopy
from typing import Any

import numpy as np
import pandas as pd

from engine.apy.ltbi_state import resolve_ltbi_state_assumptions
from engine.apy.scenario import (
    DEFAULT_COMPARATOR,
    DEFAULT_INTERVENTION,
    DIRECT_EFFECTS_SCOPE_STATEMENT,
)


EVENT_LEDGER_CONTRACT_VERSION = "ltbi_screening_event_ledger_v3"
YEAR_BIN_CONVENTION = (
    "model year 0 is [0,1), model year 1 is [1,2), and the final interval "
    "may be shorter for a non-integer follow-up horizon; programme events after "
    "follow-up retain their actual non-negative model year with withinFollowUp=false."
)
FLOAT_TOLERANCE = 1e-7


EVENT_DEFINITIONS = [
    ("population", "Population", "people", "People represented in the model arm."),
    ("eligible_population", "Eligible population", "people", "People eligible for systematic LTBI screening under the scenario."),
    ("screened", "Screened", "people", "People receiving the systematic LTBI screening test."),
    ("infected_screened", "Infected screened", "people", "Screened people with MTB infection at baseline."),
    ("uninfected_screened", "Uninfected screened", "people", "Screened people without MTB infection."),
    ("infected_at_baseline", "Infected at baseline", "people", "People with prevalent MTB infection at intervention start."),
    ("recent_ltbi_at_baseline", "Recent LTBI at baseline", "people", "Prevalent infection in the recent compartment at intervention start."),
    ("remote_ltbi_at_baseline", "Remote LTBI at baseline", "people", "Prevalent infection in the remote compartment at intervention start."),
    ("latent_infected_at_screen", "Latent infected at screen", "people", "Screened infected people who have not developed active TB by screening time."),
    ("recent_latent_at_screen", "Recent latent at screen", "people", "Screened infected people still latent and in the recent compartment at screening time."),
    ("remote_latent_at_screen", "Remote latent at screen", "people", "Screened infected people still latent and in the remote compartment at screening time."),
    ("active_tb_at_screen", "Active TB at screen", "people", "Screened infected people whose untreated active-TB onset time is at or before screening time; this is not assumed to be screen-detected active TB."),
    ("true_positive_latent", "True-positive latent", "people", "Test-positive person who is infected and latent at the time of screening."),
    ("true_positive_recent", "True-positive recent", "people", "Test-positive latent person in the recent compartment at screening time."),
    ("true_positive_remote", "True-positive remote", "people", "Test-positive latent person in the remote compartment at screening time."),
    ("test_positive_active", "Test-positive active", "people", "Test-positive person who already has active TB at the time of screening."),
    ("false_positive", "False positive", "people", "Positive LTBI test in an uninfected person."),
    ("false_negative_latent", "False-negative latent", "people", "Test-negative person who is infected and latent at the time of screening."),
    ("test_negative_active", "Test-negative active", "people", "Test-negative person who already has active TB at the time of screening."),
    ("true_negative", "True negative", "people", "Test-negative uninfected screened person."),
    ("test_positive_total", "Test-positive total", "people", "All positive LTBI test results."),
    ("test_negative_total", "Test-negative total", "people", "All negative LTBI test results."),
    ("tpt_eligible", "TPT eligible", "people", "People eligible for preventive treatment after LTBI testing: true-positive latent plus false-positive results."),
    ("tpt_started_true_positive", "TPT started, true positive", "people", "True-positive latent people starting preventive treatment."),
    ("tpt_started_recent", "TPT started, recent", "people", "True-positive recent latent people starting preventive treatment."),
    ("tpt_started_remote", "TPT started, remote", "people", "True-positive remote latent people starting preventive treatment."),
    ("tpt_started_false_positive", "TPT started, false positive", "people", "False-positive people starting preventive treatment."),
    ("tpt_started_total", "TPT started total", "people", "All preventive treatment starts."),
    ("tpt_completed_true_positive", "TPT completed, true positive", "people", "True-positive latent people completing preventive treatment."),
    ("tpt_completed_recent", "TPT completed, recent", "people", "True-positive recent latent people completing preventive treatment."),
    ("tpt_completed_remote", "TPT completed, remote", "people", "True-positive remote latent people completing preventive treatment."),
    ("tpt_completed_false_positive", "TPT completed, false positive", "people", "False-positive people completing preventive treatment."),
    ("tpt_completed_total", "TPT completed total", "people", "All preventive treatment completions."),
    ("tpt_adr_stop_true_positive", "TPT ADR stop, true positive", "people", "True-positive latent people stopping preventive treatment because of adverse events."),
    ("tpt_adr_stop_false_positive", "TPT ADR stop, false positive", "people", "False-positive people stopping preventive treatment because of adverse events."),
    ("tpt_adr_stop_total", "TPT ADR stop total", "people", "All adverse-event preventive treatment stops."),
    ("tpt_other_stop_true_positive", "TPT other stop, true positive", "people", "True-positive latent people stopping preventive treatment for non-ADR reasons."),
    ("tpt_other_stop_false_positive", "TPT other stop, false positive", "people", "False-positive people stopping preventive treatment for non-ADR reasons."),
    ("tpt_other_stop_total", "TPT other stop total", "people", "All non-ADR preventive treatment stops."),
    ("tpt_partial_course_true_positive", "Partial TPT course, true positive", "people", "True-positive latent people starting but not completing preventive treatment."),
    ("tpt_partial_course_false_positive", "Partial TPT course, false positive", "people", "False-positive people starting but not completing preventive treatment."),
    ("tpt_partial_course_total", "Partial TPT course total", "people", "All started but incomplete preventive treatment courses."),
    ("infection_effectively_treated_full", "Infection effectively treated, full", "people", "True infections effectively treated after full-course preventive treatment."),
    ("infection_effectively_treated_partial", "Infection effectively treated, partial", "people", "True infections effectively treated after partial-course preventive treatment."),
    ("infection_effectively_treated_total", "Infection effectively treated total", "people", "True infections effectively treated by preventive treatment; compatibility alias for nCuredInfection."),
    ("infection_effectively_treated_recent", "Infection effectively treated, recent", "people", "True infections effectively treated among people recent latent at screening."),
    ("infection_effectively_treated_remote", "Infection effectively treated, remote", "people", "True infections effectively treated among people remote latent at screening."),
    ("active_tb_cases", "Active TB cases", "cases", "Untreated active TB cases in the comparator arm or remaining active TB cases in the intervention arm."),
    ("active_tb_cases_prevented", "Active TB cases prevented", "cases", "Active TB cases prevented by effective treatment, timed at untreated onset."),
    ("active_tb_cases_prevented_recent", "Active TB cases prevented, recent", "cases", "Prevented active TB among people recent latent at screening."),
    ("active_tb_cases_prevented_remote", "Active TB cases prevented, remote", "cases", "Prevented active TB among people remote latent at screening."),
    ("false_positive_bcg", "False positive, BCG", "people", "False-positive LTBI tests among BCG-vaccinated screened uninfected people."),
    ("false_positive_no_bcg", "False positive, no BCG", "people", "False-positive LTBI tests among non-BCG screened uninfected people."),
    ("false_positive_due_to_bcg", "False positive due to BCG", "people", "Additive BCG-attributable false-positive tests where supported."),
]


EVENT_NAMES = [row[0] for row in EVENT_DEFINITIONS]
CORE_ZERO_COMPARATOR_EVENTS = [
    name
    for name in EVENT_NAMES
    if name
    not in {
        "population",
        "eligible_population",
        "infected_at_baseline",
        "recent_ltbi_at_baseline",
        "remote_ltbi_at_baseline",
        "active_tb_cases",
    }
]


def event_definitions_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "eventName": name,
                "label": label,
                "units": units,
                "definition": definition,
            }
            for name, label, units, definition in EVENT_DEFINITIONS
        ]
    )


def metadata_from_config(
    config: dict[str, Any],
    *,
    model_type: str,
    backend: str,
    model_version: str,
) -> dict[str, Any]:
    scenario = config.get("scenario") if isinstance(config.get("scenario"), dict) else {}
    ltbi_state = resolve_ltbi_state_assumptions(config)
    return {
        "contractVersion": EVENT_LEDGER_CONTRACT_VERSION,
        "scenarioId": config.get("scenarioLabel") or scenario.get("scenarioName"),
        "scenarioVersion": scenario.get("scenarioVersion", config.get("configVersion")),
        "populationPresetId": config.get("populationPresetId"),
        "modelType": model_type,
        "backend": backend,
        "screeningWindow": config.get("screenWindow"),
        "screeningWindowYears": config.get("screeningWindowYears", config.get("screenWindow")),
        "earlyProgressionPeriodYears": config.get("earlyProgressionPeriodYears"),
        "activeTBCalibrationHorizonYears": config.get("activeTBCalibrationHorizonYears"),
        "followUpHorizon": config.get("followHorizon"),
        "followUpHorizonYears": config.get("followUpHorizonYears", config.get("followHorizon")),
        "ltbiStateModel": ltbi_state.get("transitionModel"),
        "baselineRecentLTBIProportion": ltbi_state.get("baselineRecentLTBIProportion"),
        "recentToRemoteTransitionRatePerYear": ltbi_state.get("recentToRemoteTransitionRatePerYear"),
        "recentStateImpliedMeanResidenceTimeYears": ltbi_state.get("impliedMeanResidenceTimeYears"),
        "ltbiStateDefinition": ltbi_state.get("stateDefinition"),
        "ltbiStateAssumptionStatus": ltbi_state.get("status"),
        "ltbiStateAssumptionSource": ltbi_state.get("source"),
        "ltbiStateAssumptionNotes": ltbi_state.get("notes"),
        "ltbiStateWarning": ltbi_state.get("warning"),
        "ltbiStateProvisional": ltbi_state.get("provisional"),
        "modelVersion": model_version,
        "scopeStatement": DIRECT_EFFECTS_SCOPE_STATEMENT,
        "yearBinConvention": YEAR_BIN_CONVENTION,
        "comparator": DEFAULT_COMPARATOR,
        "intervention": DEFAULT_INTERVENTION,
    }


def make_bundle(
    *,
    metadata: dict[str, Any],
    replicate_totals_wide: pd.DataFrame,
    annual_events: pd.DataFrame,
    validation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if validation is None:
        validation = validate_event_ledger(replicate_totals_wide, annual_events)
    replicate_totals = wide_to_long_totals(replicate_totals_wide)
    return {
        "metadata": deepcopy(metadata),
        "definitions": event_definitions_frame(),
        "replicateTotals": replicate_totals.reset_index(drop=True),
        "annualEvents": annual_events.reset_index(drop=True),
        "validation": validation,
        "summaries": summarise_event_ledger(replicate_totals),
    }


def validate_event_ledger(
    replicate_totals: pd.DataFrame,
    annual_events: pd.DataFrame,
    *,
    tolerance: float = FLOAT_TOLERANCE,
) -> dict[str, Any]:
    errors: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    if replicate_totals.empty:
        errors.append(_issue("replicateTotals", "Event ledger has no replicate totals."))
        return _validation_report(errors, warnings)

    totals_by_key = {
        _row_key(row): row
        for _, row in replicate_totals.iterrows()
    }
    for _, row in replicate_totals.iterrows():
        _validate_identities(row, errors, tolerance)
        if row.get("arm") == "comparator":
            for event in CORE_ZERO_COMPARATOR_EVENTS:
                if abs(_value(row, event)) > tolerance:
                    errors.append(_issue(event, f"Comparator {event} must be zero.", row))
        if row.get("arm") == "intervention":
            comp_key = (
                row.get("modelType"),
                row.get("replicateId"),
                row.get("pairedReplicateId"),
                "comparator",
            )
            comp = totals_by_key.get(comp_key)
            if comp is not None:
                lhs = _value(comp, "active_tb_cases")
                rhs = _value(row, "active_tb_cases") + _value(row, "active_tb_cases_prevented")
                if abs(lhs - rhs) > tolerance:
                    errors.append(_issue("active_tb_identity", "Comparator active TB must equal intervention active TB plus prevented cases.", row))

    if not annual_events.empty:
        grouped = annual_events.groupby(["modelType", "replicateId", "pairedReplicateId", "arm", "eventName"], dropna=False)["value"].sum()
        for _, row in replicate_totals.iterrows():
            key_base = (row.get("modelType"), row.get("replicateId"), row.get("pairedReplicateId"), row.get("arm"))
            for event in EVENT_NAMES:
                annual = float(grouped.get((*key_base, event), 0.0))
                total = _value(row, event)
                if abs(annual - total) > tolerance:
                    errors.append(_issue("annual_reconciliation", f"Annual {event} sums to {annual}, expected {total}.", row))
        _validate_annual_cross_arm_tb_identity(annual_events, errors, tolerance)
    return _validation_report(errors, warnings)


def summarise_event_ledger(replicate_totals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (arm, event), group in replicate_totals.groupby(["arm", "eventName"], dropna=False):
        values = pd.Series(group["value"], dtype="float64")
        rows.append(
            {
                "arm": arm,
                "eventName": event,
                "mean": float(values.mean()) if not values.empty else np.nan,
                "median": float(values.median()) if not values.empty else np.nan,
                "min": float(values.min()) if not values.empty else np.nan,
                "max": float(values.max()) if not values.empty else np.nan,
                "n": int(values.count()),
            }
        )
    return pd.DataFrame(rows)


def wide_to_long_totals(wide: pd.DataFrame) -> pd.DataFrame:
    id_cols = [
        "contractVersion",
        "scenarioId",
        "scenarioVersion",
        "populationPresetId",
        "modelType",
        "backend",
        "modelVersion",
        "arm",
        "comparator",
        "intervention",
        "replicateId",
        "pairedReplicateId",
        "replicateSeed",
        "valueType",
        "screeningWindow",
        "followUpHorizon",
        "screeningWindowYears",
        "earlyProgressionPeriodYears",
        "activeTBCalibrationHorizonYears",
        "followUpHorizonYears",
    ]
    present_id_cols = [col for col in id_cols if col in wide.columns]
    rows = []
    definitions = event_definitions_frame().set_index("eventName").to_dict(orient="index")
    for _, source in wide.iterrows():
        base = {col: source.get(col) for col in present_id_cols}
        for event in EVENT_NAMES:
            meta = definitions[event]
            rows.append(
                {
                    **base,
                    "eventName": event,
                    "label": meta["label"],
                    "units": meta["units"],
                    "value": float(source.get(event, 0.0)),
                }
            )
    return pd.DataFrame(rows)


def annual_records_from_event_times(
    *,
    base: dict[str, Any],
    event_times: dict[str, np.ndarray],
    follow_horizon: float,
) -> pd.DataFrame:
    rows = []
    for event_name in EVENT_NAMES:
        values = np.asarray(event_times.get(event_name, []), dtype=float)
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue
        for (model_year, within_follow_up), count in _bin_times(values, follow_horizon).items():
            rows.append(
                {
                    **base,
                    "modelYear": model_year,
                    "timeInterval": _time_interval_label(model_year, follow_horizon),
                    "withinFollowUp": within_follow_up,
                    "eventName": event_name,
                    "value": float(count),
                }
            )
    return pd.DataFrame(rows)


def annual_records_from_weighted_events(
    *,
    base: dict[str, Any],
    event_values: dict[str, list[tuple[float, float]]],
    follow_horizon: float,
) -> pd.DataFrame:
    rows = []
    for event_name in EVENT_NAMES:
        bins: dict[tuple[int, bool], float] = {}
        for time, value in event_values.get(event_name, []):
            if value == 0 or not np.isfinite(time):
                continue
            model_year = _model_year(float(time), follow_horizon)
            within_follow_up = float(time) < float(follow_horizon)
            key = (model_year, within_follow_up)
            bins[key] = bins.get(key, 0.0) + float(value)
        for (model_year, within_follow_up), value in bins.items():
            rows.append(
                {
                    **base,
                    "modelYear": model_year,
                    "timeInterval": _time_interval_label(model_year, follow_horizon),
                    "withinFollowUp": within_follow_up,
                    "eventName": event_name,
                    "value": float(value),
                }
            )
    return pd.DataFrame(rows)


def zero_comparator_wide(
    base: dict[str, Any],
    population: float,
    active_tb_cases: float,
    eligible_population: float | None = None,
) -> dict[str, Any]:
    row = {**base, "arm": "comparator"}
    for event in EVENT_NAMES:
        row[event] = 0.0
    row["population"] = float(population)
    row["eligible_population"] = float(population if eligible_population is None else eligible_population)
    row["active_tb_cases"] = float(active_tb_cases)
    return row


def matlab_total_ledger_from_raw_rows(
    *,
    config: dict[str, Any],
    raw_rows: pd.DataFrame | list[dict[str, Any]],
    model_version: str = "matlab_apy_v9",
) -> dict[str, Any]:
    raw = pd.DataFrame(raw_rows)
    required = {
        "nScreened",
        "nLatentAtScreen",
        "nActiveAtScreen",
        "nTruePositiveLatent",
        "nTestPositiveActive",
        "nFalseNegativeLatent",
        "nTestNegativeActive",
        "nTrueNegative",
        "nFalsePositiveTests",
        "nTPTEligible",
        "nTPTStartedTruePositive",
        "nTPTStartedFalsePositive",
        "nTPTCompletedTruePositive",
        "nTPTCompletedFalsePositive",
        "nTPTADRstopTruePositive",
        "nTPTADRstopFalsePositive",
        "nTPTStoppedOtherTruePositive",
        "nTPTStoppedOtherFalsePositive",
        "nCuredInfectionFull",
        "nCuredInfectionPartial",
        "nCuredInfection",
        "nPreventedActiveTB",
        "nActiveBy20y",
    }
    missing = sorted(required.difference(raw.columns))
    metadata = metadata_from_config(
        config,
        model_type="agent_based",
        backend="matlab",
        model_version=model_version,
    )
    if missing:
        return {
            "metadata": metadata,
            "available": False,
            "annualAvailable": False,
            "missingFields": missing,
            "notes": (
                "MATLAB raw output cannot be adapted to the common total-event "
                "contract without these unambiguous cumulative fields. Annual "
                "ledger output is not fabricated from cumulative summaries."
            ),
        }

    total_rows = []
    for idx, source in raw.iterrows():
        replicate = int(source.get("rep", idx + 1))
        base = {
            "contractVersion": EVENT_LEDGER_CONTRACT_VERSION,
            "scenarioId": config.get("scenarioLabel"),
            "scenarioVersion": config.get("scenario", {}).get("scenarioVersion", config.get("configVersion"))
            if isinstance(config.get("scenario"), dict)
            else config.get("configVersion"),
            "populationPresetId": config.get("populationPresetId"),
            "modelType": "agent_based",
            "backend": "matlab",
            "modelVersion": model_version,
            "arm": "",
            "comparator": DEFAULT_COMPARATOR,
            "intervention": DEFAULT_INTERVENTION,
            "replicateId": replicate,
            "pairedReplicateId": replicate,
            "replicateSeed": source.get("seed"),
            "valueType": "simulated_count",
            "screeningWindow": float(config.get("screenWindow")),
            "followUpHorizon": float(config.get("followHorizon")),
        }
        comparator_active = float(source["nActiveBy20y"])
        total_rows.append(
            zero_comparator_wide(
                base,
                population=float(config.get("N") or 0.0),
                active_tb_cases=comparator_active,
            )
        )
        total_rows.append(
            _matlab_intervention_row(
                base,
                source,
                comparator_active,
                population=float(config.get("N")),
            )
        )
    validation = validate_event_ledger(pd.DataFrame(total_rows), pd.DataFrame())
    validation["warnings"].append(
        {
            "field": "annualEvents",
            "message": "Annual MATLAB event ledger is unavailable because current MATLAB raw output does not include event timing.",
        }
    )
    validation["hasWarnings"] = True
    bundle = make_bundle(
        metadata=metadata,
        replicate_totals_wide=pd.DataFrame(total_rows),
        annual_events=pd.DataFrame(),
        validation=validation,
    )
    bundle["annualAvailable"] = False
    return bundle


def _matlab_intervention_row(
    base: dict[str, Any],
    source: pd.Series,
    comparator_active: float,
    *,
    population: float,
) -> dict[str, Any]:
    row = {**base, "arm": "intervention"}
    for event in EVENT_NAMES:
        row[event] = 0.0
    row.update(
        {
            "population": population,
            "eligible_population": population,
            "screened": float(source["nScreened"]),
            "latent_infected_at_screen": float(source["nLatentAtScreen"]),
            "active_tb_at_screen": float(source["nActiveAtScreen"]),
            "true_positive_latent": float(source["nTruePositiveLatent"]),
            "test_positive_active": float(source["nTestPositiveActive"]),
            "false_positive": float(source["nFalsePositiveTests"]),
            "false_negative_latent": float(source["nFalseNegativeLatent"]),
            "test_negative_active": float(source["nTestNegativeActive"]),
            "true_negative": float(source["nTrueNegative"]),
            "tpt_eligible": float(source["nTPTEligible"]),
            "tpt_started_true_positive": float(source["nTPTStartedTruePositive"]),
            "tpt_started_false_positive": float(source["nTPTStartedFalsePositive"]),
            "tpt_completed_true_positive": float(source["nTPTCompletedTruePositive"]),
            "tpt_completed_false_positive": float(source["nTPTCompletedFalsePositive"]),
            "tpt_adr_stop_true_positive": float(source["nTPTADRstopTruePositive"]),
            "tpt_adr_stop_false_positive": float(source["nTPTADRstopFalsePositive"]),
            "tpt_other_stop_true_positive": float(source["nTPTStoppedOtherTruePositive"]),
            "tpt_other_stop_false_positive": float(source["nTPTStoppedOtherFalsePositive"]),
            "infection_effectively_treated_full": float(source["nCuredInfectionFull"]),
            "infection_effectively_treated_partial": float(source["nCuredInfectionPartial"]),
            "infection_effectively_treated_total": float(source["nCuredInfection"]),
            "active_tb_cases_prevented": float(source["nPreventedActiveTB"]),
        }
    )
    row["infected_screened"] = row["latent_infected_at_screen"] + row["active_tb_at_screen"]
    row["uninfected_screened"] = row["screened"] - row["infected_screened"]
    row["test_positive_total"] = row["true_positive_latent"] + row["test_positive_active"] + row["false_positive"]
    row["test_negative_total"] = row["false_negative_latent"] + row["test_negative_active"] + row["true_negative"]
    row["tpt_started_total"] = row["tpt_started_true_positive"] + row["tpt_started_false_positive"]
    row["tpt_completed_total"] = row["tpt_completed_true_positive"] + row["tpt_completed_false_positive"]
    row["tpt_adr_stop_total"] = row["tpt_adr_stop_true_positive"] + row["tpt_adr_stop_false_positive"]
    row["tpt_other_stop_total"] = row["tpt_other_stop_true_positive"] + row["tpt_other_stop_false_positive"]
    row["tpt_partial_course_true_positive"] = row["tpt_adr_stop_true_positive"] + row["tpt_other_stop_true_positive"]
    row["tpt_partial_course_false_positive"] = row["tpt_adr_stop_false_positive"] + row["tpt_other_stop_false_positive"]
    row["tpt_partial_course_total"] = row["tpt_partial_course_true_positive"] + row["tpt_partial_course_false_positive"]
    row["active_tb_cases"] = comparator_active - row["active_tb_cases_prevented"]
    for source_field, event in [
        ("nFalsePositiveTestsBCG", "false_positive_bcg"),
        ("nFalsePositiveTestsNoBCG", "false_positive_no_bcg"),
        ("nExcessFalsePositiveTestsDueToBCG", "false_positive_due_to_bcg"),
    ]:
        if source_field in source:
            row[event] = float(source[source_field])
    return row


def _bin_times(times: np.ndarray, follow_horizon: float) -> dict[tuple[int, bool], float]:
    out: dict[tuple[int, bool], float] = {}
    for time in times:
        year = _model_year(float(time), follow_horizon)
        key = (year, float(time) < float(follow_horizon))
        out[key] = out.get(key, 0.0) + 1.0
    return out


def _model_year(time: float, follow_horizon: float) -> int:
    if time < 0:
        return 0
    return int(np.floor(time))


def _time_interval_label(model_year: int, follow_horizon: float) -> str:
    start = float(model_year)
    end = start + 1.0
    return f"[{start:g},{end:g})"


def _validate_annual_cross_arm_tb_identity(
    annual_events: pd.DataFrame,
    errors: list[dict[str, Any]],
    tolerance: float,
) -> None:
    subset = annual_events[
        annual_events["eventName"].isin(["active_tb_cases", "active_tb_cases_prevented"])
    ]
    if subset.empty:
        return
    grouped = subset.groupby(
        ["modelType", "replicateId", "pairedReplicateId", "modelYear", "arm", "eventName"],
        dropna=False,
    )["value"].sum()
    years = subset[["modelType", "replicateId", "pairedReplicateId", "modelYear"]].drop_duplicates()
    for _, row in years.iterrows():
        base = (
            row["modelType"],
            row["replicateId"],
            row["pairedReplicateId"],
            row["modelYear"],
        )
        comparator = float(grouped.get((*base, "comparator", "active_tb_cases"), 0.0))
        intervention = float(grouped.get((*base, "intervention", "active_tb_cases"), 0.0))
        prevented = float(grouped.get((*base, "intervention", "active_tb_cases_prevented"), 0.0))
        if abs(comparator - (intervention + prevented)) > tolerance:
            errors.append(
                {
                    "field": "annual_active_tb_identity",
                    "message": (
                        "Annual comparator active TB must equal intervention active TB "
                        "plus prevented cases in the same model year."
                    ),
                    "modelType": row["modelType"],
                    "replicateId": row["replicateId"],
                    "modelYear": row["modelYear"],
                }
            )


def _validate_identities(row: pd.Series, errors: list[dict[str, Any]], tolerance: float) -> None:
    checks = [
        ("baseline_ltbi_split", _value(row, "infected_at_baseline"), _value(row, "recent_ltbi_at_baseline") + _value(row, "remote_ltbi_at_baseline")),
        ("screened_split", _value(row, "screened"), _value(row, "infected_screened") + _value(row, "uninfected_screened")),
        ("infected_screened_split", _value(row, "infected_screened"), _value(row, "latent_infected_at_screen") + _value(row, "active_tb_at_screen")),
        ("latent_recent_remote_split", _value(row, "latent_infected_at_screen"), _value(row, "recent_latent_at_screen") + _value(row, "remote_latent_at_screen")),
        ("true_positive_recent_remote_split", _value(row, "true_positive_latent"), _value(row, "true_positive_recent") + _value(row, "true_positive_remote")),
        (
            "screened_test_split",
            _value(row, "screened"),
            sum(_value(row, name) for name in ["true_positive_latent", "false_negative_latent", "test_positive_active", "test_negative_active", "false_positive", "true_negative"]),
        ),
        ("test_positive_total", _value(row, "test_positive_total"), _value(row, "true_positive_latent") + _value(row, "test_positive_active") + _value(row, "false_positive")),
        ("test_negative_total", _value(row, "test_negative_total"), _value(row, "false_negative_latent") + _value(row, "test_negative_active") + _value(row, "true_negative")),
        ("tpt_eligible", _value(row, "tpt_eligible"), _value(row, "true_positive_latent") + _value(row, "false_positive")),
        ("tpt_started_total", _value(row, "tpt_started_total"), _value(row, "tpt_started_true_positive") + _value(row, "tpt_started_false_positive")),
        ("tpt_started_recent_remote_split", _value(row, "tpt_started_true_positive"), _value(row, "tpt_started_recent") + _value(row, "tpt_started_remote")),
        ("tpt_completed_total", _value(row, "tpt_completed_total"), _value(row, "tpt_completed_true_positive") + _value(row, "tpt_completed_false_positive")),
        ("tpt_completed_recent_remote_split", _value(row, "tpt_completed_true_positive"), _value(row, "tpt_completed_recent") + _value(row, "tpt_completed_remote")),
        ("tpt_adr_stop_total", _value(row, "tpt_adr_stop_total"), _value(row, "tpt_adr_stop_true_positive") + _value(row, "tpt_adr_stop_false_positive")),
        ("tpt_other_stop_total", _value(row, "tpt_other_stop_total"), _value(row, "tpt_other_stop_true_positive") + _value(row, "tpt_other_stop_false_positive")),
        ("tpt_started_outcomes", _value(row, "tpt_started_total"), _value(row, "tpt_completed_total") + _value(row, "tpt_adr_stop_total") + _value(row, "tpt_other_stop_total")),
        ("effective_total", _value(row, "infection_effectively_treated_total"), _value(row, "infection_effectively_treated_full") + _value(row, "infection_effectively_treated_partial")),
        ("effective_recent_remote_split", _value(row, "infection_effectively_treated_total"), _value(row, "infection_effectively_treated_recent") + _value(row, "infection_effectively_treated_remote")),
        ("prevented_recent_remote_split", _value(row, "active_tb_cases_prevented"), _value(row, "active_tb_cases_prevented_recent") + _value(row, "active_tb_cases_prevented_remote")),
    ]
    for name, lhs, rhs in checks:
        if abs(lhs - rhs) > tolerance:
            errors.append(_issue(name, f"{name} identity failed: {lhs} != {rhs}", row))
    inequalities = [
        ("starts_le_eligible", _value(row, "tpt_started_total"), _value(row, "tpt_eligible")),
        ("completions_le_starts", _value(row, "tpt_completed_total"), _value(row, "tpt_started_total")),
        ("effective_le_true_positive_starts", _value(row, "infection_effectively_treated_total"), _value(row, "tpt_started_true_positive")),
        ("prevented_le_effective", _value(row, "active_tb_cases_prevented"), _value(row, "infection_effectively_treated_total")),
    ]
    for name, lhs, rhs in inequalities:
        if lhs - rhs > tolerance:
            errors.append(_issue(name, f"{name} failed: {lhs} > {rhs}", row))
    if _value(row, "tpt_started_false_positive") > tolerance and _value(row, "infection_effectively_treated_total") - _value(row, "tpt_started_true_positive") > tolerance:
        errors.append(_issue("false_positive_benefit", "False-positive treatment cannot produce epidemiological benefit.", row))
    for event in EVENT_NAMES:
        if _value(row, event) < -tolerance:
            errors.append(_issue("negative_count", f"{event} is materially negative.", row))


def _row_key(row: pd.Series) -> tuple[Any, Any, Any, str]:
    return (row.get("modelType"), row.get("replicateId"), row.get("pairedReplicateId"), row.get("arm"))


def _value(row: pd.Series, name: str) -> float:
    value = row.get(name, 0.0)
    if value is None or pd.isna(value):
        return 0.0
    return float(value)


def _issue(field: str, message: str, row: pd.Series | None = None) -> dict[str, Any]:
    out = {"field": field, "message": message}
    if row is not None:
        out.update(
            {
                "modelType": row.get("modelType"),
                "replicateId": row.get("replicateId"),
                "arm": row.get("arm"),
            }
        )
    return out


def _validation_report(errors: list[dict[str, Any]], warnings: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "isValid": not errors,
        "errors": errors,
        "warnings": warnings,
        "tolerance": FLOAT_TOLERANCE,
    }
