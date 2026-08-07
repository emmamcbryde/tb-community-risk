from __future__ import annotations

from itertools import product
import json
from typing import Any

import numpy as np
import pandas as pd

from engine.apy.age_distribution import broad_age_group_from_years
from engine.apy.calibration import bern_prob
from engine.apy.cohort import (
    disease_multiplier,
    get_counterfactual_no_bcg_specificity,
    get_test_performance,
    infection_probability,
)
from engine.apy.config import normalise_config
from engine.apy.eligibility import resolve_eligibility, screening_coverage_of_population
from engine.apy.event_ledger import (
    EVENT_LEDGER_CONTRACT_VERSION,
    EVENT_NAMES,
    annual_records_from_weighted_events,
    make_bundle as make_event_ledger_bundle,
    metadata_from_config as event_metadata_from_config,
    zero_comparator_wide,
)
from engine.apy.regimen import (
    apply_regimen_overrides,
    default_regimen_library,
    get_regimen_from_library,
    regimen_partial_efficacy,
    validate_regimen,
)
from engine.apy.runner import (
    BACKEND,
    EXPECTED_VALUE_MODEL_VERSION,
    _cached_calibration_from_config,
    build_strategy_metadata,
)
from engine.apy.scenario import DEFAULT_COMPARATOR, DEFAULT_INTERVENTION
from engine.apy.simulation import _simulation_options
from engine.apy.ltbi_state import (
    mixed_baseline_event_between,
    mixed_baseline_survival,
    recent_no_active_survival,
    recent_state_survival,
    recent_to_remote_latent_probability,
    remote_survival,
)
from engine.apy.validation import validate_config


QUADRATURE_POINTS = 24
_STRATA_CACHE: dict[str, pd.DataFrame] = {}


def run_expected_value(
    config: dict[str, Any] | None = None,
    *,
    quadrature_points: int = QUADRATURE_POINTS,
) -> dict[str, Any]:
    cfg = validate_config(normalise_config(config))
    calibration = _cached_calibration_from_config(cfg)
    pars = calibration["parameters"]
    reg = _build_regimen(cfg)
    opts = _simulation_options(cfg)
    opts["N"] = float(cfg["N"])
    opts["quadraturePoints"] = int(quadrature_points)
    strata = _cached_strata(pars, calibration, opts)
    strata = _allocate_screening(strata, opts)
    totals, annual_values, comparator_active_by_year = _expected_events(
        strata, opts, reg, calibration
    )

    base = _ledger_base(cfg)
    comparator = zero_comparator_wide(
        base,
        population=float(cfg["N"]),
        eligible_population=float(opts["eligiblePopulation"]),
        active_tb_cases=sum(comparator_active_by_year.values()),
    )
    for event in [
        "infected_at_baseline",
        "recent_ltbi_at_baseline",
        "remote_ltbi_at_baseline",
    ]:
        comparator[event] = float(totals.get(event, 0.0))
    intervention = {**base, "arm": "intervention"}
    intervention.update(totals)
    totals_wide = pd.DataFrame([comparator, intervention])

    comparator_annual = {
        "population": [(0.0, float(cfg["N"]))],
        "eligible_population": [(0.0, float(cfg["N"]))],
        "infected_at_baseline": [(0.0, float(totals.get("infected_at_baseline", 0.0)))],
        "recent_ltbi_at_baseline": [(0.0, float(totals.get("recent_ltbi_at_baseline", 0.0)))],
        "remote_ltbi_at_baseline": [(0.0, float(totals.get("remote_ltbi_at_baseline", 0.0)))],
        "active_tb_cases": [
            (_representative_time(year), value)
            for year, value in comparator_active_by_year.items()
        ],
    }
    intervention_active_by_year = {
        year: comparator_active_by_year.get(year, 0.0)
        - _sum_event_year(annual_values.get("active_tb_cases_prevented", []), year, opts["followHorizon"])
        for year in comparator_active_by_year
    }
    intervention_annual = {
        **annual_values,
        "population": [(0.0, float(cfg["N"]))],
        "eligible_population": [(0.0, float(cfg["N"]))],
        "active_tb_cases": [
            (_representative_time(year), value)
            for year, value in intervention_active_by_year.items()
        ],
    }
    annual = pd.concat(
        [
            annual_records_from_weighted_events(
                base={**base, "arm": "comparator"},
                event_values=comparator_annual,
                follow_horizon=opts["followHorizon"],
            ),
            annual_records_from_weighted_events(
                base={**base, "arm": "intervention"},
                event_values=intervention_annual,
                follow_horizon=opts["followHorizon"],
            ),
        ],
        ignore_index=True,
    )
    event_ledger = make_event_ledger_bundle(
        metadata=event_metadata_from_config(
            cfg,
            model_type="expected_value",
            backend=BACKEND,
            model_version=EXPECTED_VALUE_MODEL_VERSION,
        ),
        replicate_totals_wide=totals_wide,
        annual_events=annual,
    )
    return {
        "eventLedger": event_ledger,
        "strata": strata,
        "parameters": pars,
        "calibration": {k: v for k, v in calibration.items() if k != "parameters"},
        "strategy": build_strategy_metadata(cfg, reg),
        "interfaceConfig": cfg,
        "modelVersion": EXPECTED_VALUE_MODEL_VERSION,
        "backend": BACKEND,
    }


def _build_strata(pars: dict[str, Any], calibration: dict[str, Any], opts: dict[str, Any]) -> pd.DataFrame:
    rows = []
    factor_names = ["MJ", "contact", "renal", "diabetes", "smoking", "cld", "alcohol", "BCG"]
    for age, age_prob in zip(pars["exactAgeValues"], pars["exactAgeProb"]):
        age_group = broad_age_group_from_years([age])[0]
        age_idx = age_group - 1
        probs = {
            "MJ": pars["mjPrevByAge"][age_idx],
            "contact": pars["contactPrevByAge"][age_idx],
            "renal": pars["renalPrevByAge"][age_idx],
            "diabetes": pars["diabetesPrevByAge"][age_idx],
            "smoking": pars["smokingPrevByAge"][age_idx],
            "cld": pars["cldPrevByAge"][age_idx],
            "alcohol": pars["alcoholPrevByAge"][age_idx],
            "BCG": pars["totalBCGPrev"],
        }
        for flags in product([0, 1], repeat=len(factor_names)):
            flag = dict(zip(factor_names, flags))
            weight_probability = float(age_prob)
            for name in factor_names:
                weight_probability *= bern_prob(flag[name], probs[name])
            mult = float(
                disease_multiplier(
                    pars,
                    flag["MJ"],
                    flag["contact"],
                    flag["renal"],
                    flag["diabetes"],
                    flag["smoking"],
                    flag["cld"],
                    flag["alcohol"],
                )
            )
            if calibration.get("zeroInfectionPrevalence"):
                p_inf = 0.0
            else:
                p_inf = float(
                    infection_probability(
                        [age],
                        pars,
                        calibration["ageInfLogLambda"],
                        calibration["ageInfGamma"],
                        [flag["MJ"]],
                        [flag["contact"]],
                        [flag["renal"]],
                    )[0]
                )
            p_recent = float(opts["baselineRecentLTBIProportion"])
            sens, spec = get_test_performance([bool(flag["BCG"])], opts)
            no_bcg_spec = get_counterfactual_no_bcg_specificity([bool(flag["BCG"])], opts)
            rows.append(
                {
                    "weight": weight_probability,
                    "ageYears": float(age),
                    "ageGroup": int(age_group),
                    **{name: bool(flag[name]) for name in factor_names},
                    "pInfection": p_inf,
                    "pRecentLTBIAtBaseline": p_recent,
                    "pRemoteLTBIAtBaseline": 1.0 - p_recent,
                    "diseaseMultiplier": mult,
                    "testSensitivity": float(sens[0]),
                    "testSpecificity": float(spec[0]),
                    "testSpecificityNoBCG": float(no_bcg_spec[0]),
                }
            )
    out = pd.DataFrame(rows)
    out["populationWeight"] = out["weight"] * float(opts.get("N", 1.0))
    mult = out["diseaseMultiplier"].to_numpy(dtype=float)
    times = (np.arange(32, dtype=float) + 0.5) / 32 * opts["screenWindow"]
    survival_values = [
        mixed_baseline_survival(
            t,
            mult,
            calibration["lambdaEarly"],
            calibration["lambdaLate"],
            opts["recentToRemoteTransitionRatePerYear"],
            float(opts["baselineRecentLTBIProportion"]),
        )
        for t in times
    ]
    avg_latent = np.mean(survival_values, axis=0)
    horizon_survival = mixed_baseline_survival(
        opts["followHorizon"],
        mult,
        calibration["lambdaEarly"],
        calibration["lambdaLate"],
        opts["recentToRemoteTransitionRatePerYear"],
        float(opts["baselineRecentLTBIProportion"]),
    )
    preventable = np.maximum(avg_latent - horizon_survival, 0.0)
    out["ltbiRiskScore"] = out["pInfection"]
    out["cureTargetScore"] = out["pInfection"] * avg_latent
    out["preventTargetScore"] = out["pInfection"] * preventable
    return out


def _cached_strata(pars: dict[str, Any], calibration: dict[str, Any], opts: dict[str, Any]) -> pd.DataFrame:
    key = json.dumps(
        {
            "ageInfLogLambda": calibration["ageInfLogLambda"],
            "ageInfGamma": calibration["ageInfGamma"],
            "lambdaEarly": calibration["lambdaEarly"],
            "lambdaLate": calibration["lambdaLate"],
            "screenWindow": opts["screenWindow"],
            "followHorizon": opts["followHorizon"],
            "earlyProgressionPeriodYears": opts["earlyProgressionPeriodYears"],
            "baselineRecentLTBIProportion": opts["baselineRecentLTBIProportion"],
            "recentToRemoteTransitionRatePerYear": opts[
                "recentToRemoteTransitionRatePerYear"
            ],
            "testType": opts["testType"],
            "testSensitivity": opts["testSensitivity"],
            "testSpecificity": opts["testSpecificity"],
            "tstSensitivity": opts["tstSensitivity"],
            "tstSpecificityBCG": opts["tstSpecificityBCG"],
            "tstSpecificityNoBCG": opts["tstSpecificityNoBCG"],
            "totalBCGPrev": pars.get("totalBCGPrev"),
            "totalFemalePrev": pars.get("totalFemalePrev"),
            "mjPrevByAge": pars.get("mjPrevByAge"),
            "contactPrevByAge": pars.get("contactPrevByAge"),
            "renalPrevByAge": pars.get("renalPrevByAge"),
            "diabetesPrevByAge": pars.get("diabetesPrevByAge"),
            "smokingPrevByAge": pars.get("smokingPrevByAge"),
            "cldPrevByAge": pars.get("cldPrevByAge"),
            "alcoholPrevByAge": pars.get("alcoholPrevByAge"),
            "populationWeight": opts["N"],
        },
        sort_keys=True,
        default=str,
    )
    if key not in _STRATA_CACHE:
        _STRATA_CACHE[key] = _build_strata(pars, calibration, opts)
    return _STRATA_CACHE[key].copy()


def _allocate_screening(strata: pd.DataFrame, opts: dict[str, Any]) -> pd.DataFrame:
    out = strata.copy()
    population = out["populationWeight"].sum()
    requested = min(max(float(opts["screenCoverage"]) * float(opts["eligiblePopulation"]), 0.0), population)
    out["screenedWeight"] = 0.0
    if requested <= 0:
        return out
    strategy = str(opts["screeningStrategy"]).lower()
    if strategy == "random":
        out["screenedWeight"] = out["populationWeight"] * float(opts["screenCoverageOfPopulation"])
        return out
    score_column = {
        "ltbi": "ltbiRiskScore",
        "cure": "cureTargetScore",
        "prevent": "preventTargetScore",
    }[strategy]
    order = out.sort_values(
        [score_column, "pInfection", "diseaseMultiplier", "ageYears"],
        ascending=[False, False, False, False],
        kind="mergesort",
    ).index
    remaining = requested
    for idx in order:
        take = min(float(out.at[idx, "populationWeight"]), remaining)
        out.at[idx, "screenedWeight"] = take
        remaining -= take
        if remaining <= 1e-10:
            break
    return out


def _expected_events(
    strata: pd.DataFrame,
    opts: dict[str, Any],
    reg: dict[str, Any],
    calibration: dict[str, Any],
) -> tuple[dict[str, float], dict[str, list[tuple[float, float]]], dict[int, float]]:
    totals = {name: 0.0 for name in EVENT_NAMES}
    annual_values: dict[str, list[tuple[float, float]]] = {name: [] for name in EVENT_NAMES}
    comparator_active_by_year = _comparator_active_by_year(strata, opts, calibration)
    q_points = int(opts.get("quadraturePoints", QUADRATURE_POINTS))
    q_times = (np.arange(q_points, dtype=float) + 0.5) / q_points * opts["screenWindow"]
    p_complete = float(reg["pComplete"])
    p_adr = float(reg["pADRstop"])
    p_other = max(1.0 - p_complete - p_adr, 0.0)
    duration = float(reg["months"]) / 12.0
    adr_dose = float(opts["partialDoseFractionADR"])
    other_dose = float(opts["partialDoseFractionOther"])
    partial_eff_adr = float(regimen_partial_efficacy(reg, adr_dose))
    partial_eff_other = float(regimen_partial_efficacy(reg, other_dose))
    screened_weight = strata["screenedWeight"].to_numpy(dtype=float)
    selected = screened_weight > 0
    point_weight = screened_weight[selected] / q_points
    p_inf = strata.loc[selected, "pInfection"].to_numpy(dtype=float)
    p_recent = strata.loc[selected, "pRecentLTBIAtBaseline"].to_numpy(dtype=float)
    p_remote = strata.loc[selected, "pRemoteLTBIAtBaseline"].to_numpy(dtype=float)
    sens = strata.loc[selected, "testSensitivity"].to_numpy(dtype=float)
    spec = strata.loc[selected, "testSpecificity"].to_numpy(dtype=float)
    no_bcg_spec = strata.loc[selected, "testSpecificityNoBCG"].to_numpy(dtype=float)
    mult = strata.loc[selected, "diseaseMultiplier"].to_numpy(dtype=float)
    bcg = strata.loc[selected, "BCG"].to_numpy(dtype=bool)
    baseline_infected = (
        strata["populationWeight"].to_numpy(dtype=float)
        * strata["pInfection"].to_numpy(dtype=float)
    )
    baseline_recent = (
        baseline_infected
        * strata["pRecentLTBIAtBaseline"].to_numpy(dtype=float)
    )
    baseline_remote = (
        baseline_infected
        * strata["pRemoteLTBIAtBaseline"].to_numpy(dtype=float)
    )
    _add_many(
        annual_values,
        totals,
        0.0,
        {
            "infected_at_baseline": baseline_infected,
            "recent_ltbi_at_baseline": baseline_recent,
            "remote_ltbi_at_baseline": baseline_remote,
        },
    )
    for screen_time in q_times:
        infected_screened = point_weight * p_inf
        recent_latent = infected_screened * p_recent * _recent_state_survival_array(
            screen_time, mult, calibration, opts
        )
        remote_latent = infected_screened * (
            p_remote * _remote_survival_array(screen_time, mult, calibration, opts)
            + p_recent
            * _recent_to_remote_latent_array(screen_time, mult, calibration, opts)
        )
        latent = recent_latent + remote_latent
        active = infected_screened - latent
        uninfected = point_weight * (1.0 - p_inf)
        tp_recent = recent_latent * sens
        tp_remote = remote_latent * sens
        tp_latent = tp_recent + tp_remote
        tp_active = active * sens
        fn_latent = latent * (1.0 - sens)
        tn_active = active * (1.0 - sens)
        false_positive = uninfected * (1.0 - spec)
        true_negative = uninfected * spec
        false_positive_due_bcg = np.where(
            bcg,
            uninfected * np.maximum((1.0 - spec) - (1.0 - no_bcg_spec), 0.0),
            0.0,
        )
        _add_many(
            annual_values,
            totals,
            screen_time,
            {
                "screened": point_weight,
                "infected_screened": infected_screened,
                "uninfected_screened": uninfected,
                "latent_infected_at_screen": latent,
                "recent_latent_at_screen": recent_latent,
                "remote_latent_at_screen": remote_latent,
                "active_tb_at_screen": active,
                "true_positive_latent": tp_latent,
                "true_positive_recent": tp_recent,
                "true_positive_remote": tp_remote,
                "test_positive_active": tp_active,
                "false_negative_latent": fn_latent,
                "test_negative_active": tn_active,
                "false_positive": false_positive,
                "true_negative": true_negative,
                "test_positive_total": tp_latent + tp_active + false_positive,
                "test_negative_total": fn_latent + tn_active + true_negative,
                "false_positive_bcg": np.where(bcg, false_positive, 0.0),
                "false_positive_no_bcg": np.where(~bcg, false_positive, 0.0),
                "false_positive_due_to_bcg": false_positive_due_bcg,
            },
        )
        tp_started = tp_latent * float(opts["pStartTPT"])
        tp_started_recent = tp_recent * float(opts["pStartTPT"])
        tp_started_remote = tp_remote * float(opts["pStartTPT"])
        fp_started = false_positive * float(opts["pStartTPT"])
        _add_many(
            annual_values,
            totals,
            screen_time,
            {
                "tpt_eligible": tp_latent + false_positive,
                "tpt_started_true_positive": tp_started,
                "tpt_started_recent": tp_started_recent,
                "tpt_started_remote": tp_started_remote,
                "tpt_started_false_positive": fp_started,
                "tpt_started_total": tp_started + fp_started,
            },
        )
        _treatment_outcomes(
            annual_values,
            totals,
            screen_time,
            duration,
            float(tp_started.sum()),
            float(fp_started.sum()),
            p_complete,
            p_adr,
            p_other,
            adr_dose,
            other_dose,
            partial_eff_adr,
            partial_eff_other,
            reg,
            mult,
            calibration,
            opts,
            tp_started_by_stratum=tp_started,
            tp_started_recent_by_stratum=tp_started_recent,
            tp_started_remote_by_stratum=tp_started_remote,
        )
    totals["population"] = float(opts.get("N", 0.0))
    totals["eligible_population"] = float(opts.get("eligiblePopulation", opts.get("N", 0.0)))
    totals["active_tb_cases_prevented"] = totals.pop(
        "infection_effectively_treated_total_prevented_marker", 0.0
    )
    totals["active_tb_cases"] = sum(comparator_active_by_year.values()) - totals["active_tb_cases_prevented"]
    return totals, annual_values, comparator_active_by_year


def _treatment_outcomes(
    annual_values: dict[str, list[tuple[float, float]]],
    totals: dict[str, float],
    screen_time: float,
    duration: float,
    tp_started: float,
    fp_started: float,
    p_complete: float,
    p_adr: float,
    p_other: float,
    adr_dose: float,
    other_dose: float,
    partial_eff_adr: float,
    partial_eff_other: float,
    reg: dict[str, Any],
    mult: np.ndarray,
    calibration: dict[str, Any],
    opts: dict[str, Any],
    *,
    tp_started_by_stratum: np.ndarray,
    tp_started_recent_by_stratum: np.ndarray,
    tp_started_remote_by_stratum: np.ndarray,
) -> None:
    full_time = screen_time + duration
    adr_time = screen_time + duration * adr_dose
    other_time = screen_time + duration * other_dose
    tp_complete = tp_started * p_complete
    fp_complete = fp_started * p_complete
    tp_complete_recent = float(tp_started_recent_by_stratum.sum()) * p_complete
    tp_complete_remote = float(tp_started_remote_by_stratum.sum()) * p_complete
    tp_adr = tp_started * p_adr
    fp_adr = fp_started * p_adr
    tp_other = tp_started * p_other
    fp_other = fp_started * p_other
    _add(annual_values, totals, "tpt_completed_true_positive", full_time, tp_complete)
    _add(annual_values, totals, "tpt_completed_recent", full_time, tp_complete_recent)
    _add(annual_values, totals, "tpt_completed_remote", full_time, tp_complete_remote)
    _add(annual_values, totals, "tpt_completed_false_positive", full_time, fp_complete)
    _add(annual_values, totals, "tpt_completed_total", full_time, tp_complete + fp_complete)
    _add(annual_values, totals, "tpt_adr_stop_true_positive", adr_time, tp_adr)
    _add(annual_values, totals, "tpt_adr_stop_false_positive", adr_time, fp_adr)
    _add(annual_values, totals, "tpt_adr_stop_total", adr_time, tp_adr + fp_adr)
    _add(annual_values, totals, "tpt_other_stop_true_positive", other_time, tp_other)
    _add(annual_values, totals, "tpt_other_stop_false_positive", other_time, fp_other)
    _add(annual_values, totals, "tpt_other_stop_total", other_time, tp_other + fp_other)
    _add(annual_values, totals, "tpt_partial_course_true_positive", adr_time, tp_adr)
    _add(annual_values, totals, "tpt_partial_course_false_positive", adr_time, fp_adr)
    _add(annual_values, totals, "tpt_partial_course_total", adr_time, tp_adr + fp_adr)
    _add(annual_values, totals, "tpt_partial_course_true_positive", other_time, tp_other)
    _add(annual_values, totals, "tpt_partial_course_false_positive", other_time, fp_other)
    _add(annual_values, totals, "tpt_partial_course_total", other_time, tp_other + fp_other)

    full_effective_by_stratum = (
        (
            tp_started_recent_by_stratum
            * _conditional_recent_survival_array(screen_time, full_time, mult, calibration, opts)
            + tp_started_remote_by_stratum
            * _conditional_remote_survival_array(screen_time, full_time, mult, calibration, opts)
        )
        * p_complete
        * float(reg["effFull"])
    )
    adr_effective_by_stratum = (
        (
            tp_started_recent_by_stratum
            * _conditional_recent_survival_array(screen_time, adr_time, mult, calibration, opts)
            + tp_started_remote_by_stratum
            * _conditional_remote_survival_array(screen_time, adr_time, mult, calibration, opts)
        )
        * p_adr
        * partial_eff_adr
    )
    other_effective_by_stratum = (
        (
            tp_started_recent_by_stratum
            * _conditional_recent_survival_array(screen_time, other_time, mult, calibration, opts)
            + tp_started_remote_by_stratum
            * _conditional_remote_survival_array(screen_time, other_time, mult, calibration, opts)
        )
        * p_other
        * partial_eff_other
    )
    full_effective = float(full_effective_by_stratum.sum())
    adr_effective = float(adr_effective_by_stratum.sum())
    other_effective = float(other_effective_by_stratum.sum())
    _add(annual_values, totals, "infection_effectively_treated_full", full_time, full_effective)
    _add(annual_values, totals, "infection_effectively_treated_partial", adr_time, adr_effective)
    _add(annual_values, totals, "infection_effectively_treated_partial", other_time, other_effective)
    _add(annual_values, totals, "infection_effectively_treated_total", full_time, full_effective)
    _add(annual_values, totals, "infection_effectively_treated_total", adr_time, adr_effective)
    _add(annual_values, totals, "infection_effectively_treated_total", other_time, other_effective)
    recent_effective = float(
        (
            tp_started_recent_by_stratum
            * (
                p_complete
                * float(reg["effFull"])
                * _conditional_recent_survival_array(screen_time, full_time, mult, calibration, opts)
                + p_adr
                * partial_eff_adr
                * _conditional_recent_survival_array(screen_time, adr_time, mult, calibration, opts)
                + p_other
                * partial_eff_other
                * _conditional_recent_survival_array(screen_time, other_time, mult, calibration, opts)
            )
        ).sum()
    )
    remote_effective = float(
        (
            tp_started_remote_by_stratum
            * (
                p_complete
                * float(reg["effFull"])
                * _conditional_remote_survival_array(screen_time, full_time, mult, calibration, opts)
                + p_adr
                * partial_eff_adr
                * _conditional_remote_survival_array(screen_time, adr_time, mult, calibration, opts)
                + p_other
                * partial_eff_other
                * _conditional_remote_survival_array(screen_time, other_time, mult, calibration, opts)
            )
        ).sum()
    )
    _add(annual_values, totals, "infection_effectively_treated_recent", full_time, recent_effective)
    _add(annual_values, totals, "infection_effectively_treated_remote", full_time, remote_effective)
    for stop_time, probability, efficacy in [
        (full_time, p_complete, float(reg["effFull"])),
        (adr_time, p_adr, partial_eff_adr),
        (other_time, p_other, partial_eff_other),
    ]:
        horizon = float(opts["followHorizon"])
        for year in range(int(np.ceil(horizon))):
            interval_start = max(float(year), stop_time)
            interval_end = min(float(year + 1), horizon)
            if interval_end <= interval_start:
                continue
            prevented = float(
                (
                    (
                        tp_started_recent_by_stratum
                        * _conditional_recent_interval_event_array(
                            screen_time,
                            interval_start,
                            interval_end,
                            mult,
                            calibration,
                            opts,
                        )
                        + tp_started_remote_by_stratum
                        * _conditional_remote_interval_event_array(
                            screen_time,
                            interval_start,
                            interval_end,
                            mult,
                            calibration,
                            opts,
                        )
                    )
                    * probability
                    * efficacy
                ).sum()
            )
            prevented_recent = float(
                (
                    tp_started_recent_by_stratum
                    * probability
                    * efficacy
                    * _conditional_recent_interval_event_array(
                        screen_time,
                        interval_start,
                        interval_end,
                        mult,
                        calibration,
                        opts,
                    )
                ).sum()
            )
            prevented_remote = float(
                (
                    tp_started_remote_by_stratum
                    * probability
                    * efficacy
                    * _conditional_remote_interval_event_array(
                        screen_time,
                        interval_start,
                        interval_end,
                        mult,
                        calibration,
                        opts,
                    )
                ).sum()
            )
            totals["infection_effectively_treated_total_prevented_marker"] = totals.get("infection_effectively_treated_total_prevented_marker", 0.0) + prevented
            totals["active_tb_cases_prevented_recent"] = totals.get(
                "active_tb_cases_prevented_recent", 0.0
            ) + prevented_recent
            totals["active_tb_cases_prevented_remote"] = totals.get(
                "active_tb_cases_prevented_remote", 0.0
            ) + prevented_remote
            if prevented:
                annual_values["active_tb_cases_prevented"].append(
                    (_representative_time(year), prevented)
                )
                annual_values["active_tb_cases_prevented_recent"].append(
                    (_representative_time(year), prevented_recent)
                )
                annual_values["active_tb_cases_prevented_remote"].append(
                    (_representative_time(year), prevented_remote)
                )


def _add(
    annual_values: dict[str, list[tuple[float, float]]],
    totals: dict[str, float],
    event: str,
    time: float,
    value: float,
) -> None:
    if value == 0:
        return
    totals[event] = totals.get(event, 0.0) + float(value)
    annual_values[event].append((float(time), float(value)))


def _add_many(
    annual_values: dict[str, list[tuple[float, float]]],
    totals: dict[str, float],
    time: float,
    values: dict[str, np.ndarray],
) -> None:
    for event, vector in values.items():
        _add(annual_values, totals, event, time, float(np.asarray(vector, dtype=float).sum()))


def _comparator_active_by_year(
    strata: pd.DataFrame,
    opts: dict[str, Any],
    calibration: dict[str, Any],
) -> dict[int, float]:
    out: dict[int, float] = {}
    horizon = float(opts["followHorizon"])
    n_years = int(np.ceil(horizon))
    for year in range(n_years):
        start = float(year)
        end = min(start + 1.0, horizon)
        if end <= start:
            continue
        value = 0.0
        value = float(
            (
                strata["populationWeight"].to_numpy(dtype=float)
                * strata["pInfection"].to_numpy(dtype=float)
                * _event_between_array(
                    start,
                    end,
                    strata["diseaseMultiplier"].to_numpy(dtype=float),
                    calibration,
                    opts,
                )
            ).sum()
        )
        out[year] = value
    return out


def _event_between(start: float, end: float, mult: float, calibration: dict[str, Any], opts: dict[str, Any]) -> float:
    return max(_survival(start, mult, calibration, opts) - _survival(end, mult, calibration, opts), 0.0)


def _event_between_array(start: float, end: float, mult: np.ndarray, calibration: dict[str, Any], opts: dict[str, Any]) -> np.ndarray:
    return mixed_baseline_event_between(
        start,
        end,
        mult,
        calibration["lambdaEarly"],
        calibration["lambdaLate"],
        opts["recentToRemoteTransitionRatePerYear"],
        opts["baselineRecentLTBIProportion"],
    )


def _conditional_event_between(
    screen_time: float,
    stop_time: float,
    horizon: float,
    mult: float,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> float:
    if horizon <= stop_time:
        return 0.0
    denominator = max(_survival(screen_time, mult, calibration, opts), np.finfo(float).eps)
    return max((_survival(stop_time, mult, calibration, opts) - _survival(horizon, mult, calibration, opts)) / denominator, 0.0)


def _conditional_event_between_array(
    screen_time: float,
    stop_time: float,
    horizon: float,
    mult: np.ndarray,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> np.ndarray:
    if horizon <= stop_time:
        return np.zeros_like(mult, dtype=float)
    denominator = np.maximum(
        _survival_array(screen_time, mult, calibration, opts), np.finfo(float).eps
    )
    return np.maximum(
        (
            _survival_array(stop_time, mult, calibration, opts)
            - _survival_array(horizon, mult, calibration, opts)
        )
        / denominator,
        0.0,
    )


def _conditional_interval_event_array(
    screen_time: float,
    interval_start: float,
    interval_end: float,
    mult: np.ndarray,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> np.ndarray:
    denominator = np.maximum(
        _survival_array(screen_time, mult, calibration, opts), np.finfo(float).eps
    )
    return np.maximum(
        (
            _survival_array(interval_start, mult, calibration, opts)
            - _survival_array(interval_end, mult, calibration, opts)
        )
        / denominator,
        0.0,
    )


def _conditional_survival(
    screen_time: float,
    stop_time: float,
    mult: float,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> float:
    denominator = max(_survival(screen_time, mult, calibration, opts), np.finfo(float).eps)
    return min(max(_survival(stop_time, mult, calibration, opts) / denominator, 0.0), 1.0)


def _conditional_survival_array(
    screen_time: float,
    stop_time: float,
    mult: np.ndarray,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> np.ndarray:
    denominator = np.maximum(
        _survival_array(screen_time, mult, calibration, opts), np.finfo(float).eps
    )
    return np.clip(_survival_array(stop_time, mult, calibration, opts) / denominator, 0.0, 1.0)


def _survival(time: float, mult: float, calibration: dict[str, Any], opts: dict[str, Any]) -> float:
    return float(
        survival_probability(
            time,
            mult,
            calibration["lambdaEarly"],
            calibration["lambdaLate"],
            opts["earlyProgressionPeriodYears"],
        )
    )


def _survival_array(time: float, mult: np.ndarray, calibration: dict[str, Any], opts: dict[str, Any]) -> np.ndarray:
    return mixed_baseline_survival(
        time,
        mult,
        calibration["lambdaEarly"],
        calibration["lambdaLate"],
        opts["recentToRemoteTransitionRatePerYear"],
        opts["baselineRecentLTBIProportion"],
    )


def _recent_state_survival_array(time: float, mult: np.ndarray, calibration: dict[str, Any], opts: dict[str, Any]) -> np.ndarray:
    return recent_state_survival(
        time,
        mult,
        calibration["lambdaEarly"],
        opts["recentToRemoteTransitionRatePerYear"],
    )


def _recent_to_remote_latent_array(time: float, mult: np.ndarray, calibration: dict[str, Any], opts: dict[str, Any]) -> np.ndarray:
    return recent_to_remote_latent_probability(
        time,
        mult,
        calibration["lambdaEarly"],
        calibration["lambdaLate"],
        opts["recentToRemoteTransitionRatePerYear"],
    )


def _remote_survival_array(time: float, mult: np.ndarray, calibration: dict[str, Any], opts: dict[str, Any]) -> np.ndarray:
    return remote_survival(time, mult, calibration["lambdaLate"])


def _conditional_recent_survival_array(
    screen_time: float,
    stop_time: float,
    mult: np.ndarray,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> np.ndarray:
    delta = max(stop_time - screen_time, 0.0)
    return np.clip(
        recent_no_active_survival(
            delta,
            mult,
            calibration["lambdaEarly"],
            calibration["lambdaLate"],
            opts["recentToRemoteTransitionRatePerYear"],
        ),
        0.0,
        1.0,
    )


def _conditional_remote_survival_array(
    screen_time: float,
    stop_time: float,
    mult: np.ndarray,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> np.ndarray:
    delta = max(stop_time - screen_time, 0.0)
    return np.clip(
        remote_survival(delta, mult, calibration["lambdaLate"]),
        0.0,
        1.0,
    )


def _conditional_recent_interval_event_array(
    screen_time: float,
    interval_start: float,
    interval_end: float,
    mult: np.ndarray,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> np.ndarray:
    start_delta = max(interval_start - screen_time, 0.0)
    end_delta = max(interval_end - screen_time, 0.0)
    return np.maximum(
        recent_no_active_survival(
            start_delta,
            mult,
            calibration["lambdaEarly"],
            calibration["lambdaLate"],
            opts["recentToRemoteTransitionRatePerYear"],
        )
        - recent_no_active_survival(
            end_delta,
            mult,
            calibration["lambdaEarly"],
            calibration["lambdaLate"],
            opts["recentToRemoteTransitionRatePerYear"],
        ),
        0.0,
    )


def _conditional_remote_interval_event_array(
    screen_time: float,
    interval_start: float,
    interval_end: float,
    mult: np.ndarray,
    calibration: dict[str, Any],
    opts: dict[str, Any],
) -> np.ndarray:
    start_delta = max(interval_start - screen_time, 0.0)
    end_delta = max(interval_end - screen_time, 0.0)
    return np.maximum(
        remote_survival(start_delta, mult, calibration["lambdaLate"])
        - remote_survival(end_delta, mult, calibration["lambdaLate"]),
        0.0,
    )


def _sum_event_year(events: list[tuple[float, float]], year: int, horizon: float) -> float:
    total = 0.0
    for time, value in events:
        event_year = int(np.floor(max(time, 0.0)))
        if event_year == year:
            total += value
    return total


def _representative_time(year: int) -> float:
    return float(year) + 0.5


def _ledger_base(config: dict[str, Any]) -> dict[str, Any]:
    scenario = config.get("scenario") if isinstance(config.get("scenario"), dict) else {}
    return {
        "contractVersion": EVENT_LEDGER_CONTRACT_VERSION,
        "scenarioId": config.get("scenarioLabel") or scenario.get("scenarioName"),
        "scenarioVersion": scenario.get("scenarioVersion", config.get("configVersion")),
        "populationPresetId": config.get("populationPresetId"),
        "modelType": "expected_value",
        "backend": BACKEND,
        "modelVersion": EXPECTED_VALUE_MODEL_VERSION,
        "arm": "",
        "comparator": DEFAULT_COMPARATOR,
        "intervention": DEFAULT_INTERVENTION,
        "replicateId": 0,
        "pairedReplicateId": 0,
        "replicateSeed": None,
        "valueType": "expected",
        "screeningWindow": float(config["screeningWindowYears"]),
        "screeningWindowYears": float(config["screeningWindowYears"]),
        "earlyProgressionPeriodYears": float(config["earlyProgressionPeriodYears"]),
        "activeTBCalibrationHorizonYears": float(config["activeTBCalibrationHorizonYears"]),
        "followUpHorizon": float(config["followUpHorizonYears"]),
        "followUpHorizonYears": float(config["followUpHorizonYears"]),
    }


def _build_regimen(config: dict[str, Any]) -> dict[str, Any]:
    library = default_regimen_library(config.get("partialShortCourseMode") or "threshold80")
    reg = get_regimen_from_library(config["regimen"], library)
    reg = apply_regimen_overrides(reg, config)
    validation = validate_regimen(reg)
    if not validation["isValid"]:
        raise ValueError("; ".join(validation["errors"]))
    return reg
