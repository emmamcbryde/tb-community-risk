from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from engine.apy.calibration import calibrate_from_config
from engine.apy.cohort import (
    add_targeting_scores,
    draw_base_population,
    draw_untreated_active_times,
    get_counterfactual_no_bcg_specificity as cohort_counterfactual_specificity,
    get_test_performance as cohort_test_performance,
    make_rng,
    select_screened_people,
)
from engine.apy.config import normalise_config
from engine.apy.eligibility import resolve_eligibility, screening_coverage_of_population
from engine.apy.regimen import (
    apply_regimen_overrides,
    default_regimen_library,
    get_regimen_from_library,
    regimen_partial_efficacy,
    validate_regimen,
)
from engine.apy.timing import resolve_time_settings


RAW_FIELDS = [
    "testType",
    "screeningStrategy",
    "regimen",
    "nScreened",
    "nInfected",
    "nLatentAtScreen",
    "nActiveAtScreen",
    "nTruePositiveLatent",
    "nTestPositiveActive",
    "nFalseNegativeLatent",
    "nTestNegativeActive",
    "nTrueNegative",
    "nTPTEligible",
    "nTPTStartedTruePositive",
    "nTPTStartedFalsePositive",
    "nTPTCompletedTruePositive",
    "nTPTCompletedFalsePositive",
    "nTPTADRstopTruePositive",
    "nTPTADRstopFalsePositive",
    "nTPTStoppedOtherTruePositive",
    "nTPTStoppedOtherFalsePositive",
    "nTestPositive",
    "nTestPositiveNonActive",
    "nIGRApos",
    "nFalsePositiveTests",
    "nFalsePositiveTestsBCG",
    "nFalsePositiveTestsNoBCG",
    "nFalsePositiveTreated",
    "nFalsePositiveTreatedBCG",
    "nFalsePositiveTreatedNoBCG",
    "nFalsePositiveCompleted",
    "nFalsePositiveCompletedBCG",
    "nFalsePositiveCompletedNoBCG",
    "nExcessFalsePositiveTestsDueToBCG",
    "nExcessCoursesStartedDueToBCG",
    "nExcessCoursesCompletedDueToBCG",
    "nExcessCoursesDueToBCG",
    "nScreenedBCG",
    "nScreenedNoBCG",
    "nUninfectedScreenedBCG",
    "nUninfectedScreenedNoBCG",
    "nStartTPT",
    "nCompleteTPT",
    "nADRstop",
    "nStoppedOther",
    "nPartialCourses",
    "nTotalCoursesStarted",
    "nTotalCoursesCompleted",
    "nCuredInfection",
    "nCuredInfectionFull",
    "nCuredInfectionPartial",
    "nPreventedActiveTB",
    "nPreventedActiveTBFull",
    "nPreventedActiveTBPartial",
    "nActiveBy2y",
    "nActiveBy20y",
    "completionRateObserved",
    "adrRateObserved",
    "partialCourseRateObserved",
    "propCoursesFalsePositive",
    "falsePositiveRateObservedBCG",
    "falsePositiveRateObservedNoBCG",
    "bcgAttributableFalsePositiveRateObserved",
    "propFalsePositiveTestsDueToBCGAmongBCG",
    "propCoursesDueToBCGAmongFalsePositiveCoursesBCG",
    "NNS_cureInfection",
    "NNS_preventActiveTB",
    "NNS_falsePositiveTreated",
    "NNT_started_cureInfection",
    "NNT_started_preventActiveTB",
]


COHORT_COLUMNS = [
    "id",
    "ageGroup",
    "ageYears",
    "female",
    "BCG",
    "MJ",
    "contact",
    "renal",
    "diabetes",
    "smoking",
    "chronicLungDisease",
    "alcoholDrugs",
    "pInfection",
    "infected",
    "diseaseMultiplier",
    "ltbiRiskScore",
    "cureTargetScore",
    "preventTargetScore",
    "screenPriorityScore",
    "screenPriorityRank",
    "tActiveUntreated",
    "screened",
    "screenTime",
    "activeAtScreen",
    "latentAtScreen",
    "testSensitivityUsed",
    "testSpecificityUsed",
    "testSpecificityNoBCGCounterfactual",
    "testPositive",
    "falsePositiveTest",
    "falsePositiveCounterfactualNoBCG",
    "falsePositiveDueToBCG",
    "eligibleTPT",
    "startedTPT",
    "adrStop",
    "stoppedOther",
    "completedTPT",
    "doseFractionTaken",
    "treatmentStopTime",
    "fullEffAssigned",
    "partialEffAssigned",
    "falsePositiveTreated",
    "falsePositiveTreatedBCG",
    "falsePositiveCompleted",
    "falsePositiveCompletedBCG",
    "excessCourseStartedDueToBCG",
    "excessCourseCompletedDueToBCG",
    "curedInfection",
    "curedInfectionFull",
    "curedInfectionPartial",
    "preventedActiveTB",
    "preventedActiveTBFull",
    "preventedActiveTBPartial",
]


def simulate_one_cohort(
    n: int,
    pars: dict[str, Any],
    reg: dict[str, Any],
    calibration: dict[str, Any],
    config: dict[str, Any],
    rng,
    return_cohort: bool = True,
) -> dict[str, Any]:
    opts = _simulation_options(config)
    population = draw_base_population(
        n,
        pars,
        calibration["ageInfLogLambda"],
        calibration["ageInfGamma"],
        rng,
    )
    population = add_targeting_scores(
        population,
        calibration["lambdaEarly"],
        calibration["lambdaLate"],
        opts["screenWindow"],
        opts["followHorizon"],
        opts["earlyProgressionPeriodYears"],
    )
    t_active = draw_untreated_active_times(
        population["infected"],
        population["diseaseMultiplier"],
        calibration["lambdaEarly"],
        calibration["lambdaLate"],
        opts["screenWindow"],
        rng,
        opts["earlyProgressionPeriodYears"],
    )

    screened, screen_priority_score, screen_priority_rank = select_screened_people(
        population["ltbiRiskScore"],
        population["cureTargetScore"],
        population["preventTargetScore"],
        opts["screenCoverageOfPopulation"],
        opts["screeningStrategy"],
        rng,
    )
    screen_time = np.full(n, np.inf)
    if screened.any():
        screen_time[screened] = rng.random(int(screened.sum())) * opts["screenWindow"]

    infected = population["infected"]
    active_at_screen = screened & infected & (t_active <= screen_time)
    latent_at_screen = screened & infected & (t_active > screen_time)

    test_sensitivity, test_specificity = get_test_performance(population["BCG"], opts)
    no_bcg_specificity = get_counterfactual_no_bcg_specificity(population["BCG"], opts)
    false_positive_rate = 1 - test_specificity
    false_positive_rate_no_bcg = 1 - no_bcg_specificity

    test_positive = np.zeros(n, dtype=bool)
    false_positive_counterfactual_no_bcg = np.zeros(n, dtype=bool)
    inf_screened = screened & infected
    uninf_screened = screened & ~infected
    test_positive[inf_screened] = (
        rng.random(int(inf_screened.sum())) < test_sensitivity[inf_screened]
    )
    if uninf_screened.any():
        u_uninf = rng.random(int(uninf_screened.sum()))
        test_positive[uninf_screened] = u_uninf < false_positive_rate[uninf_screened]
        false_positive_counterfactual_no_bcg[uninf_screened] = (
            u_uninf < false_positive_rate_no_bcg[uninf_screened]
        )

    false_positive_test = screened & ~infected & test_positive
    true_positive_latent = latent_at_screen & test_positive
    test_positive_active = active_at_screen & test_positive
    false_negative_latent = latent_at_screen & ~test_positive
    test_negative_active = active_at_screen & ~test_positive
    true_negative = uninf_screened & ~test_positive
    false_positive_test_bcg = false_positive_test & population["BCG"]
    false_positive_test_no_bcg = false_positive_test & ~population["BCG"]
    false_positive_due_to_bcg = (
        false_positive_test & ~false_positive_counterfactual_no_bcg
    )

    eligible_tpt = screened & test_positive & ~active_at_screen
    started_tpt = eligible_tpt & (rng.random(n) < opts["pStartTPT"])

    adr_stop = np.zeros(n, dtype=bool)
    completed_tpt = np.zeros(n, dtype=bool)
    stopped_other = np.zeros(n, dtype=bool)
    dose_fraction = np.zeros(n)
    course_stop_time = np.full(n, np.inf)
    partial_eff_assigned = np.zeros(n)
    full_eff_assigned = np.zeros(n)

    treat_idx = np.flatnonzero(started_tpt)
    if len(treat_idx) > 0:
        adr_draw = rng.random(len(treat_idx)) < reg["pADRstop"]
        adr_stop[treat_idx[adr_draw]] = True
        remain_idx = treat_idx[~adr_draw]
        if (1 - reg["pADRstop"]) > 0:
            p_complete_given_no_adr = reg["pComplete"] / (1 - reg["pADRstop"])
        else:
            p_complete_given_no_adr = 0
        p_complete_given_no_adr = min(max(p_complete_given_no_adr, 0), 1)

        if len(remain_idx) > 0:
            comp_draw = rng.random(len(remain_idx)) < p_complete_given_no_adr
            completed_tpt[remain_idx[comp_draw]] = True
            stopped_other[remain_idx[~comp_draw]] = True

        dose_fraction[completed_tpt] = 1
        dose_fraction[adr_stop] = opts["partialDoseFractionADR"]
        dose_fraction[stopped_other] = opts["partialDoseFractionOther"]
        course_stop_time[started_tpt] = (
            screen_time[started_tpt] + (reg["months"] / 12) * dose_fraction[started_tpt]
        )
        full_eff_assigned[completed_tpt] = reg["effFull"]
        incomplete = started_tpt & ~completed_tpt
        if incomplete.any():
            partial_eff_assigned[incomplete] = regimen_partial_efficacy(
                reg, dose_fraction[incomplete].tolist()
            )

    false_positive_treated = started_tpt & ~infected
    true_positive_treated = started_tpt & latent_at_screen
    false_positive_completed = completed_tpt & ~infected
    true_positive_completed = completed_tpt & latent_at_screen
    true_positive_adr_stop = adr_stop & latent_at_screen
    false_positive_adr_stop = adr_stop & ~infected
    true_positive_stopped_other = stopped_other & latent_at_screen
    false_positive_stopped_other = stopped_other & ~infected
    false_positive_treated_bcg = false_positive_treated & population["BCG"]
    false_positive_treated_no_bcg = false_positive_treated & ~population["BCG"]
    false_positive_completed_bcg = false_positive_completed & population["BCG"]
    false_positive_completed_no_bcg = false_positive_completed & ~population["BCG"]
    excess_course_started_due_to_bcg = false_positive_treated & false_positive_due_to_bcg
    excess_course_completed_due_to_bcg = (
        false_positive_completed & false_positive_due_to_bcg
    )

    protected_full = (
        completed_tpt
        & infected
        & (t_active > course_stop_time)
        & (rng.random(n) < full_eff_assigned)
    )
    protected_partial = (
        (started_tpt & ~completed_tpt)
        & infected
        & (t_active > course_stop_time)
        & (rng.random(n) < partial_eff_assigned)
    )
    protected_any = protected_full | protected_partial
    cured_infection = protected_any
    prevented_active_tb = protected_any & (t_active <= opts["followHorizon"])
    prevented_active_tb_full = protected_full & (t_active <= opts["followHorizon"])
    prevented_active_tb_partial = protected_partial & (t_active <= opts["followHorizon"])
    active_by_2y = infected & (t_active <= opts["activeTBCalibrationHorizonYears"])
    active_by_20y = infected & (t_active <= opts["followHorizon"])

    raw = _build_raw(
        opts,
        reg,
        population,
        screened,
        uninf_screened,
        active_at_screen,
        latent_at_screen,
        true_positive_latent,
        test_positive_active,
        false_negative_latent,
        test_negative_active,
        true_negative,
        test_positive,
        false_positive_test,
        false_positive_test_bcg,
        false_positive_test_no_bcg,
        false_positive_treated,
        false_positive_treated_bcg,
        false_positive_treated_no_bcg,
        false_positive_completed,
        false_positive_completed_bcg,
        false_positive_completed_no_bcg,
        false_positive_due_to_bcg,
        excess_course_started_due_to_bcg,
        excess_course_completed_due_to_bcg,
        eligible_tpt,
        started_tpt,
        true_positive_treated,
        completed_tpt,
        true_positive_completed,
        adr_stop,
        true_positive_adr_stop,
        false_positive_adr_stop,
        stopped_other,
        true_positive_stopped_other,
        false_positive_stopped_other,
        cured_infection,
        protected_full,
        protected_partial,
        prevented_active_tb,
        prevented_active_tb_full,
        prevented_active_tb_partial,
        active_by_2y,
        active_by_20y,
    )
    event_ledger_data = _build_event_ledger_data(
        opts,
        n,
        population,
        t_active,
        screened,
        screen_time,
        active_at_screen,
        latent_at_screen,
        true_positive_latent,
        test_positive_active,
        false_negative_latent,
        test_negative_active,
        true_negative,
        false_positive_test,
        false_positive_test_bcg,
        false_positive_test_no_bcg,
        false_positive_due_to_bcg,
        eligible_tpt,
        started_tpt,
        true_positive_treated,
        false_positive_treated,
        completed_tpt,
        true_positive_completed,
        false_positive_completed,
        adr_stop,
        true_positive_adr_stop,
        false_positive_adr_stop,
        stopped_other,
        true_positive_stopped_other,
        false_positive_stopped_other,
        cured_infection,
        protected_full,
        protected_partial,
        prevented_active_tb,
        course_stop_time,
    )

    cohort = None
    if return_cohort:
        cohort = _build_cohort_dataframe(
            population,
            screen_priority_score,
            screen_priority_rank,
            t_active,
            screened,
            screen_time,
            active_at_screen,
            latent_at_screen,
            test_sensitivity,
            test_specificity,
            no_bcg_specificity,
            test_positive,
            false_positive_test,
            false_positive_counterfactual_no_bcg,
            false_positive_due_to_bcg,
            eligible_tpt,
            started_tpt,
            adr_stop,
            stopped_other,
            completed_tpt,
            dose_fraction,
            course_stop_time,
            full_eff_assigned,
            partial_eff_assigned,
            false_positive_treated,
            false_positive_treated_bcg,
            false_positive_completed,
            false_positive_completed_bcg,
            excess_course_started_due_to_bcg,
            excess_course_completed_due_to_bcg,
            cured_infection,
            protected_full,
            protected_partial,
            prevented_active_tb,
            prevented_active_tb_full,
            prevented_active_tb_partial,
        )
    return {"raw": raw, "cohort": cohort, "eventLedgerData": event_ledger_data}


def simulate_one_cohort_from_config(
    config: dict[str, Any] | None = None,
    n: int | None = None,
    seed: int | None = None,
    return_cohort: bool = True,
) -> dict[str, Any]:
    cfg = normalise_config(config)
    if n is not None:
        cfg["N"] = int(n)
    if seed is not None:
        cfg["seed"] = int(seed)
    calibration = calibrate_from_config(cfg)
    library = default_regimen_library(cfg.get("partialShortCourseMode") or "threshold80")
    reg = get_regimen_from_library(cfg["regimen"], library)
    reg = apply_regimen_overrides(reg, cfg)
    validation = validate_regimen(reg)
    if not validation["isValid"]:
        raise ValueError("; ".join(validation["errors"]))
    rng = make_rng(int(cfg["seed"]))
    return simulate_one_cohort(
        int(cfg["N"]),
        calibration["parameters"],
        reg,
        calibration,
        cfg,
        rng,
        return_cohort=return_cohort,
    )


def safe_fraction(numerator, denominator) -> float:
    return float("nan") if denominator <= 0 else numerator / denominator


def safe_divide(numerator, denominator) -> float:
    return float("inf") if denominator <= 0 else numerator / denominator


def coerce_probability(value, default=None) -> float:
    if value is None or value == []:
        if default is None:
            raise ValueError("Probability value is missing and no default was supplied.")
        return float(default)
    prob = float(value)
    if prob < 0 or prob > 1:
        raise ValueError(f"Probability must be in [0,1], got {prob}.")
    return prob


def get_test_performance(bcg, config) -> tuple[np.ndarray, np.ndarray]:
    opts = _simulation_options(config)
    return cohort_test_performance(bcg, opts)


def get_counterfactual_no_bcg_specificity(bcg, config) -> np.ndarray:
    opts = _simulation_options(config)
    return cohort_counterfactual_specificity(bcg, opts)


def _simulation_options(config: dict[str, Any]) -> dict[str, Any]:
    cfg = normalise_config(config)
    timing = resolve_time_settings(cfg)
    eligibility = resolve_eligibility(cfg)
    screen_coverage_of_population = screening_coverage_of_population(
        cfg, eligibility["number"]
    )
    tst_no_bcg = cfg.get("tstSpecificityNoBCG")
    test_specificity = coerce_probability(cfg.get("testSpecificity"), 0.98)
    if tst_no_bcg is None or tst_no_bcg == []:
        tst_no_bcg = test_specificity
    return {
        "testType": str(cfg["testType"]).upper(),
        "screeningStrategy": str(cfg["screeningStrategy"]).lower(),
        "regimen": str(cfg["regimen"]).upper(),
        "screenWindow": timing["screeningWindowYears"],
        "screeningWindowYears": timing["screeningWindowYears"],
        "earlyProgressionPeriodYears": timing["earlyProgressionPeriodYears"],
        "activeTBCalibrationHorizonYears": timing["activeTBCalibrationHorizonYears"],
        "followHorizon": timing["followUpHorizonYears"],
        "followUpHorizonYears": timing["followUpHorizonYears"],
        "eligiblePopulation": float(eligibility["number"]),
        "screenCoverage": coerce_probability(cfg["screenCoverage"], 0.30),
        "screenCoverageOfPopulation": screen_coverage_of_population,
        "pStartTPT": coerce_probability(cfg.get("pStartTPT"), 0.85),
        "partialDoseFractionADR": coerce_probability(
            cfg.get("partialDoseFractionADR"), 0.30
        ),
        "partialDoseFractionOther": coerce_probability(
            cfg.get("partialDoseFractionOther"), 0.60
        ),
        "testSensitivity": coerce_probability(cfg.get("testSensitivity"), 0.95),
        "testSpecificity": test_specificity,
        "tstSensitivity": coerce_probability(cfg.get("tstSensitivity"), 0.80),
        "tstSpecificityBCG": coerce_probability(cfg.get("tstSpecificityBCG"), 0.55),
        "tstSpecificityNoBCG": coerce_probability(tst_no_bcg, 0.97),
    }


def _build_raw(opts: dict[str, Any], reg: dict[str, Any], population, *arrays) -> dict[str, Any]:
    (
        screened,
        uninf_screened,
        active_at_screen,
        latent_at_screen,
        true_positive_latent,
        test_positive_active,
        false_negative_latent,
        test_negative_active,
        true_negative,
        test_positive,
        false_positive_test,
        false_positive_test_bcg,
        false_positive_test_no_bcg,
        false_positive_treated,
        false_positive_treated_bcg,
        false_positive_treated_no_bcg,
        false_positive_completed,
        false_positive_completed_bcg,
        false_positive_completed_no_bcg,
        false_positive_due_to_bcg,
        excess_course_started_due_to_bcg,
        excess_course_completed_due_to_bcg,
        eligible_tpt,
        started_tpt,
        true_positive_treated,
        completed_tpt,
        true_positive_completed,
        adr_stop,
        true_positive_adr_stop,
        false_positive_adr_stop,
        stopped_other,
        true_positive_stopped_other,
        false_positive_stopped_other,
        cured_infection,
        protected_full,
        protected_partial,
        prevented_active_tb,
        prevented_active_tb_full,
        prevented_active_tb_partial,
        active_by_2y,
        active_by_20y,
    ) = arrays
    bcg = population["BCG"]
    infected = population["infected"]
    n_screened = int(screened.sum())
    n_started = int(started_tpt.sum())
    n_completed = int(completed_tpt.sum())
    n_adr_stop = int(adr_stop.sum())
    n_stopped_other = int(stopped_other.sum())
    n_cured = int(cured_infection.sum())
    n_prevented = int(prevented_active_tb.sum())
    n_false_positive_treated = int(false_positive_treated.sum())
    n_false_positive_treated_bcg = int(false_positive_treated_bcg.sum())
    n_false_positive_tests_bcg = int(false_positive_test_bcg.sum())
    n_excess_false_positive_due_to_bcg = int(false_positive_due_to_bcg.sum())
    n_excess_courses_started_due_to_bcg = int(excess_course_started_due_to_bcg.sum())

    return {
        "testType": opts["testType"],
        "screeningStrategy": opts["screeningStrategy"],
        "regimen": reg["label"],
        "nScreened": n_screened,
        "nInfected": int(infected.sum()),
        "nLatentAtScreen": int(latent_at_screen.sum()),
        "nActiveAtScreen": int(active_at_screen.sum()),
        "nTruePositiveLatent": int(true_positive_latent.sum()),
        "nTestPositiveActive": int(test_positive_active.sum()),
        "nFalseNegativeLatent": int(false_negative_latent.sum()),
        "nTestNegativeActive": int(test_negative_active.sum()),
        "nTrueNegative": int(true_negative.sum()),
        "nTPTEligible": int(eligible_tpt.sum()),
        "nTPTStartedTruePositive": int(true_positive_treated.sum()),
        "nTPTStartedFalsePositive": int(false_positive_treated.sum()),
        "nTPTCompletedTruePositive": int(true_positive_completed.sum()),
        "nTPTCompletedFalsePositive": int(false_positive_completed.sum()),
        "nTPTADRstopTruePositive": int(true_positive_adr_stop.sum()),
        "nTPTADRstopFalsePositive": int(false_positive_adr_stop.sum()),
        "nTPTStoppedOtherTruePositive": int(true_positive_stopped_other.sum()),
        "nTPTStoppedOtherFalsePositive": int(false_positive_stopped_other.sum()),
        "nTestPositive": int((test_positive & screened).sum()),
        "nTestPositiveNonActive": int((test_positive & screened & ~active_at_screen).sum()),
        "nIGRApos": int((test_positive & screened & ~active_at_screen).sum()),
        "nFalsePositiveTests": int(false_positive_test.sum()),
        "nFalsePositiveTestsBCG": n_false_positive_tests_bcg,
        "nFalsePositiveTestsNoBCG": int(false_positive_test_no_bcg.sum()),
        "nFalsePositiveTreated": n_false_positive_treated,
        "nFalsePositiveTreatedBCG": n_false_positive_treated_bcg,
        "nFalsePositiveTreatedNoBCG": int(false_positive_treated_no_bcg.sum()),
        "nFalsePositiveCompleted": int(false_positive_completed.sum()),
        "nFalsePositiveCompletedBCG": int(false_positive_completed_bcg.sum()),
        "nFalsePositiveCompletedNoBCG": int(false_positive_completed_no_bcg.sum()),
        "nExcessFalsePositiveTestsDueToBCG": n_excess_false_positive_due_to_bcg,
        "nExcessCoursesStartedDueToBCG": n_excess_courses_started_due_to_bcg,
        "nExcessCoursesCompletedDueToBCG": int(excess_course_completed_due_to_bcg.sum()),
        "nExcessCoursesDueToBCG": n_excess_courses_started_due_to_bcg,
        "nScreenedBCG": int((screened & bcg).sum()),
        "nScreenedNoBCG": int((screened & ~bcg).sum()),
        "nUninfectedScreenedBCG": int((uninf_screened & bcg).sum()),
        "nUninfectedScreenedNoBCG": int((uninf_screened & ~bcg).sum()),
        "nStartTPT": n_started,
        "nCompleteTPT": n_completed,
        "nADRstop": n_adr_stop,
        "nStoppedOther": n_stopped_other,
        "nPartialCourses": n_adr_stop + n_stopped_other,
        "nTotalCoursesStarted": n_started,
        "nTotalCoursesCompleted": n_completed,
        "nCuredInfection": n_cured,
        "nCuredInfectionFull": int(protected_full.sum()),
        "nCuredInfectionPartial": int(protected_partial.sum()),
        "nPreventedActiveTB": n_prevented,
        "nPreventedActiveTBFull": int(prevented_active_tb_full.sum()),
        "nPreventedActiveTBPartial": int(prevented_active_tb_partial.sum()),
        "nActiveBy2y": int(active_by_2y.sum()),
        "nActiveBy20y": int(active_by_20y.sum()),
        "completionRateObserved": safe_fraction(n_completed, n_started),
        "adrRateObserved": safe_fraction(n_adr_stop, n_started),
        "partialCourseRateObserved": safe_fraction(n_adr_stop + n_stopped_other, n_started),
        "propCoursesFalsePositive": safe_fraction(n_false_positive_treated, n_started),
        "falsePositiveRateObservedBCG": safe_fraction(
            n_false_positive_tests_bcg, int((uninf_screened & bcg).sum())
        ),
        "falsePositiveRateObservedNoBCG": safe_fraction(
            int(false_positive_test_no_bcg.sum()), int((uninf_screened & ~bcg).sum())
        ),
        "bcgAttributableFalsePositiveRateObserved": safe_fraction(
            n_excess_false_positive_due_to_bcg, int((uninf_screened & bcg).sum())
        ),
        "propFalsePositiveTestsDueToBCGAmongBCG": safe_fraction(
            n_excess_false_positive_due_to_bcg, n_false_positive_tests_bcg
        ),
        "propCoursesDueToBCGAmongFalsePositiveCoursesBCG": safe_fraction(
            n_excess_courses_started_due_to_bcg, n_false_positive_treated_bcg
        ),
        "NNS_cureInfection": safe_divide(n_screened, n_cured),
        "NNS_preventActiveTB": safe_divide(n_screened, n_prevented),
        "NNS_falsePositiveTreated": safe_divide(n_screened, n_false_positive_treated),
        "NNT_started_cureInfection": safe_divide(n_started, n_cured),
        "NNT_started_preventActiveTB": safe_divide(n_started, n_prevented),
    }


def _build_event_ledger_data(
    opts: dict[str, Any],
    n: int,
    population,
    t_active,
    screened,
    screen_time,
    active_at_screen,
    latent_at_screen,
    true_positive_latent,
    test_positive_active,
    false_negative_latent,
    test_negative_active,
    true_negative,
    false_positive_test,
    false_positive_test_bcg,
    false_positive_test_no_bcg,
    false_positive_due_to_bcg,
    eligible_tpt,
    started_tpt,
    true_positive_treated,
    false_positive_treated,
    completed_tpt,
    true_positive_completed,
    false_positive_completed,
    adr_stop,
    true_positive_adr_stop,
    false_positive_adr_stop,
    stopped_other,
    true_positive_stopped_other,
    false_positive_stopped_other,
    cured_infection,
    protected_full,
    protected_partial,
    prevented_active_tb,
    course_stop_time,
) -> dict[str, Any]:
    infected = population["infected"]
    uninfected_screened = screened & ~infected
    active_by_horizon = infected & (t_active <= opts["followHorizon"])
    intervention_active = active_by_horizon & ~prevented_active_tb
    partial_true = true_positive_adr_stop | true_positive_stopped_other
    partial_false = false_positive_adr_stop | false_positive_stopped_other
    totals = {
        "population": int(n),
        "eligible_population": float(opts["eligiblePopulation"]),
        "screened": int(screened.sum()),
        "infected_screened": int((screened & infected).sum()),
        "uninfected_screened": int(uninfected_screened.sum()),
        "latent_infected_at_screen": int(latent_at_screen.sum()),
        "active_tb_at_screen": int(active_at_screen.sum()),
        "true_positive_latent": int(true_positive_latent.sum()),
        "test_positive_active": int(test_positive_active.sum()),
        "false_positive": int(false_positive_test.sum()),
        "false_negative_latent": int(false_negative_latent.sum()),
        "test_negative_active": int(test_negative_active.sum()),
        "true_negative": int(true_negative.sum()),
        "test_positive_total": int((true_positive_latent | test_positive_active | false_positive_test).sum()),
        "test_negative_total": int((false_negative_latent | test_negative_active | true_negative).sum()),
        "tpt_eligible": int(eligible_tpt.sum()),
        "tpt_started_true_positive": int(true_positive_treated.sum()),
        "tpt_started_false_positive": int(false_positive_treated.sum()),
        "tpt_started_total": int(started_tpt.sum()),
        "tpt_completed_true_positive": int(true_positive_completed.sum()),
        "tpt_completed_false_positive": int(false_positive_completed.sum()),
        "tpt_completed_total": int(completed_tpt.sum()),
        "tpt_adr_stop_true_positive": int(true_positive_adr_stop.sum()),
        "tpt_adr_stop_false_positive": int(false_positive_adr_stop.sum()),
        "tpt_adr_stop_total": int(adr_stop.sum()),
        "tpt_other_stop_true_positive": int(true_positive_stopped_other.sum()),
        "tpt_other_stop_false_positive": int(false_positive_stopped_other.sum()),
        "tpt_other_stop_total": int(stopped_other.sum()),
        "tpt_partial_course_true_positive": int(partial_true.sum()),
        "tpt_partial_course_false_positive": int(partial_false.sum()),
        "tpt_partial_course_total": int((partial_true | partial_false).sum()),
        "infection_effectively_treated_full": int(protected_full.sum()),
        "infection_effectively_treated_partial": int(protected_partial.sum()),
        "infection_effectively_treated_total": int(cured_infection.sum()),
        "active_tb_cases": int(intervention_active.sum()),
        "active_tb_cases_prevented": int(prevented_active_tb.sum()),
        "false_positive_bcg": int(false_positive_test_bcg.sum()),
        "false_positive_no_bcg": int(false_positive_test_no_bcg.sum()),
        "false_positive_due_to_bcg": int(false_positive_due_to_bcg.sum()),
    }
    event_times = {
        "screened": screen_time[screened],
        "infected_screened": screen_time[screened & infected],
        "uninfected_screened": screen_time[uninfected_screened],
        "latent_infected_at_screen": screen_time[latent_at_screen],
        "active_tb_at_screen": screen_time[active_at_screen],
        "true_positive_latent": screen_time[true_positive_latent],
        "test_positive_active": screen_time[test_positive_active],
        "false_positive": screen_time[false_positive_test],
        "false_negative_latent": screen_time[false_negative_latent],
        "test_negative_active": screen_time[test_negative_active],
        "true_negative": screen_time[true_negative],
        "test_positive_total": screen_time[true_positive_latent | test_positive_active | false_positive_test],
        "test_negative_total": screen_time[false_negative_latent | test_negative_active | true_negative],
        "tpt_eligible": screen_time[eligible_tpt],
        "tpt_started_true_positive": screen_time[true_positive_treated],
        "tpt_started_false_positive": screen_time[false_positive_treated],
        "tpt_started_total": screen_time[started_tpt],
        "tpt_completed_true_positive": course_stop_time[true_positive_completed],
        "tpt_completed_false_positive": course_stop_time[false_positive_completed],
        "tpt_completed_total": course_stop_time[completed_tpt],
        "tpt_adr_stop_true_positive": course_stop_time[true_positive_adr_stop],
        "tpt_adr_stop_false_positive": course_stop_time[false_positive_adr_stop],
        "tpt_adr_stop_total": course_stop_time[adr_stop],
        "tpt_other_stop_true_positive": course_stop_time[true_positive_stopped_other],
        "tpt_other_stop_false_positive": course_stop_time[false_positive_stopped_other],
        "tpt_other_stop_total": course_stop_time[stopped_other],
        "tpt_partial_course_true_positive": course_stop_time[partial_true],
        "tpt_partial_course_false_positive": course_stop_time[partial_false],
        "tpt_partial_course_total": course_stop_time[partial_true | partial_false],
        "infection_effectively_treated_full": course_stop_time[protected_full],
        "infection_effectively_treated_partial": course_stop_time[protected_partial],
        "infection_effectively_treated_total": course_stop_time[cured_infection],
        "active_tb_cases": t_active[intervention_active],
        "active_tb_cases_prevented": t_active[prevented_active_tb],
        "false_positive_bcg": screen_time[false_positive_test_bcg],
        "false_positive_no_bcg": screen_time[false_positive_test_no_bcg],
        "false_positive_due_to_bcg": screen_time[false_positive_due_to_bcg],
    }
    return {
        "interventionTotals": totals,
        "comparatorActiveTimes": t_active[active_by_horizon],
        "interventionEventTimes": event_times,
    }


def _build_cohort_dataframe(population: dict[str, np.ndarray], *arrays) -> pd.DataFrame:
    (
        screen_priority_score,
        screen_priority_rank,
        t_active,
        screened,
        screen_time,
        active_at_screen,
        latent_at_screen,
        test_sensitivity,
        test_specificity,
        no_bcg_specificity,
        test_positive,
        false_positive_test,
        false_positive_counterfactual_no_bcg,
        false_positive_due_to_bcg,
        eligible_tpt,
        started_tpt,
        adr_stop,
        stopped_other,
        completed_tpt,
        dose_fraction,
        course_stop_time,
        full_eff_assigned,
        partial_eff_assigned,
        false_positive_treated,
        false_positive_treated_bcg,
        false_positive_completed,
        false_positive_completed_bcg,
        excess_course_started_due_to_bcg,
        excess_course_completed_due_to_bcg,
        cured_infection,
        protected_full,
        protected_partial,
        prevented_active_tb,
        prevented_active_tb_full,
        prevented_active_tb_partial,
    ) = arrays
    n = len(population["ageYears"])
    return pd.DataFrame(
        {
            "id": np.arange(1, n + 1),
            "ageGroup": population["ageGroup"],
            "ageYears": population["ageYears"],
            "female": population["female"],
            "BCG": population["BCG"],
            "MJ": population["MJ"],
            "contact": population["contact"],
            "renal": population["renal"],
            "diabetes": population["diabetes"],
            "smoking": population["smoking"],
            "chronicLungDisease": population["cld"],
            "alcoholDrugs": population["alcohol"],
            "pInfection": population["pInfection"],
            "infected": population["infected"],
            "diseaseMultiplier": population["diseaseMultiplier"],
            "ltbiRiskScore": population["ltbiRiskScore"],
            "cureTargetScore": population["cureTargetScore"],
            "preventTargetScore": population["preventTargetScore"],
            "screenPriorityScore": screen_priority_score,
            "screenPriorityRank": screen_priority_rank,
            "tActiveUntreated": t_active,
            "screened": screened,
            "screenTime": screen_time,
            "activeAtScreen": active_at_screen,
            "latentAtScreen": latent_at_screen,
            "testSensitivityUsed": test_sensitivity,
            "testSpecificityUsed": test_specificity,
            "testSpecificityNoBCGCounterfactual": no_bcg_specificity,
            "testPositive": test_positive,
            "falsePositiveTest": false_positive_test,
            "falsePositiveCounterfactualNoBCG": false_positive_counterfactual_no_bcg,
            "falsePositiveDueToBCG": false_positive_due_to_bcg,
            "eligibleTPT": eligible_tpt,
            "startedTPT": started_tpt,
            "adrStop": adr_stop,
            "stoppedOther": stopped_other,
            "completedTPT": completed_tpt,
            "doseFractionTaken": dose_fraction,
            "treatmentStopTime": course_stop_time,
            "fullEffAssigned": full_eff_assigned,
            "partialEffAssigned": partial_eff_assigned,
            "falsePositiveTreated": false_positive_treated,
            "falsePositiveTreatedBCG": false_positive_treated_bcg,
            "falsePositiveCompleted": false_positive_completed,
            "falsePositiveCompletedBCG": false_positive_completed_bcg,
            "excessCourseStartedDueToBCG": excess_course_started_due_to_bcg,
            "excessCourseCompletedDueToBCG": excess_course_completed_due_to_bcg,
            "curedInfection": cured_infection,
            "curedInfectionFull": protected_full,
            "curedInfectionPartial": protected_partial,
            "preventedActiveTB": prevented_active_tb,
            "preventedActiveTBFull": prevented_active_tb_full,
            "preventedActiveTBPartial": prevented_active_tb_partial,
        },
        columns=COHORT_COLUMNS,
    )
