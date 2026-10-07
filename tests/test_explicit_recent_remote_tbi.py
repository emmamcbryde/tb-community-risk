from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
import unittest

import numpy as np

from engine.apy.explicit_recent_remote_tbi import (
    ACTIVE_TB_OBSERVATION_SCHEMA_VERSION,
    AGE85_PLUS_MAX_DEFAULT,
    ANALYSIS_BASIS,
    ASSIGNMENT_CONTRACT_VERSION,
    CALIBRATION_CONTRACT_VERSION,
    CONFIG_CONTRACT_VERSION,
    RECENT_HAZARD_SHAPE,
    REMOTE_HAZARD_SHAPE,
    STATE_RECENT,
    STATE_REMOTE_ONLY,
    STATE_UNINFECTED,
    CalibrationError,
    age_support_calibration_sensitivity,
    assess_calibration_feasibility,
    build_explicit_recent_remote_config,
    calibrate_recent_remote_hazards,
    calibration_result_from_dict,
    deterministic_recent_remote_assignment,
    deterministic_recent_remote_assignment_from_config,
    exposure_durations_for_ages,
    expected_time_since_most_recent_event,
    explicit_recent_remote_config_from_dict,
    explicit_recent_remote_config_hash,
    explicit_recent_remote_config_json,
    infection_time_quantile,
    infection_timing_specification,
    mutually_exclusive_state_probabilities,
    population_weighted_prevalences,
    draw_recent_remote_states_for_ages,
    stochastic_recent_remote_population_assignment,
    stochastic_recent_remote_population_assignment_from_config,
    validate_active_tb_observation,
    validate_active_tb_observations,
    validate_age_distribution,
)
from engine.apy.explicit_recent_remote_progression import (
    OBSERVATION_MODEL_BINOMIAL,
    OBSERVATION_MODEL_POISSON,
    POLICY_BOTH_HAZARDS_SUPPLIED,
    POLICY_EXTERNAL_VALIDATION_ONLY,
    POLICY_EXTERNAL_HAZARDS,
    POLICY_FIXED_EARLY_LATE_RATIO,
    POLICY_FIXED_RATIO_FIT_SCALE,
    POLICY_JOINT_LIKELIHOOD,
    POLICY_JOINT_EARLY_REMOTE_HAZARDS,
    POLICY_ONE_HAZARD_SUPPLIED,
    POLICY_VALIDATION_ONLY,
    PROGRESSION_CONTRACT_VERSION,
    PROGRESSION_CALIBRATION_CONTRACT_VERSION,
    RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC,
    RISK_POLICY_NONE,
    RISK_POLICY_REVIEWED_HAZARD_MULTIPLIERS,
    STATE_BASELINE_ACTIVE_TB,
    TARGET_BASELINE_PREVALENCE,
    TARGET_PROSPECTIVE_INCIDENT,
    TARGET_RETROSPECTIVE_NOTIFICATION,
    TARGET_SCREEN_DETECTED,
    ProgressionCalibrationError,
    active_tb_observation_log_likelihood,
    assess_observation_prospective_calibration_eligibility,
    assess_progression_policy_identifiability,
    ascertainment_progression_scale_confounding_diagnostic,
    baseline_active_tb_state_sequence_specification,
    build_progression_calibration_policy,
    build_risk_factor_application_policy,
    classify_active_tb_observation_target,
    competing_risk_progression_probability,
    evaluate_external_hazard_policy,
    evaluate_validation_only_policy,
    expected_progression_events,
    expected_prospective_incident_cases_for_observation,
    fit_fixed_ratio_progression_scale,
    fixed_ratio_calibration_multiplier_sensitivity,
    progression_calibration_policy_contract_options,
    progression_calibration_policy_hash,
    progression_calibration_policy_json,
    progression_calibration_policy_options,
    progression_cumulative_hazard,
    progression_identifiability_assessment,
    progression_piecewise_hazard,
    progression_probability,
    progression_sampling_specification,
    progression_survival_probability,
    progression_time_quantile,
    risk_factor_multiplier_diagnostics,
    risk_factor_multiplier_for_row,
    remaining_early_risk_years,
    synthetic_natural_history_cost_invariance_example,
    target_progression_sensitivity_vector,
    worked_progression_calibration_examples,
    worked_progression_diagnostic_table,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
FROZEN_COMMIT = "03cc16e4a52a55e10019dab16bd4f599571c4b90"


class ExplicitRecentRemoteConfigTests(unittest.TestCase):
    def test_configuration_serializes_and_hashes_deterministically(self) -> None:
        cfg = build_explicit_recent_remote_config(
            enabled=True,
            recent_tbi_target=0.07,
            remote_only_tbi_target=0.12,
            target_source="unit test source",
            target_reference_year=2026,
            review_status="reviewed_test_fixture",
            notes="configuration round-trip fixture",
        )

        payload = json.loads(explicit_recent_remote_config_json(cfg))
        restored = explicit_recent_remote_config_from_dict(payload)
        shuffled_payload = dict(reversed(list(payload.items())))

        self.assertEqual(payload["configContractVersion"], CONFIG_CONTRACT_VERSION)
        self.assertEqual(restored.as_dict(), cfg.as_dict())
        self.assertEqual(
            explicit_recent_remote_config_hash(payload),
            explicit_recent_remote_config_hash(shuffled_payload),
        )
        self.assertEqual(payload["age85PlusMax"], AGE85_PLUS_MAX_DEFAULT)
        self.assertIn("85+", payload["ageSupportProvenance"])
        self.assertEqual(payload["recentHazardShape"], RECENT_HAZARD_SHAPE)
        self.assertEqual(payload["remoteHazardShape"], REMOTE_HAZARD_SHAPE)

    def test_age85_plus_max_participates_in_configuration_hash(self) -> None:
        base = build_explicit_recent_remote_config(
            enabled=True,
            recent_tbi_target=0.01,
            remote_only_tbi_target=0.02,
            target_source="unit test source",
            review_status="reviewed_test_fixture",
            age85_plus_max=89,
        )
        wider = build_explicit_recent_remote_config(
            enabled=True,
            recent_tbi_target=0.01,
            remote_only_tbi_target=0.02,
            target_source="unit test source",
            review_status="reviewed_test_fixture",
            age85_plus_max=95,
        )

        self.assertNotEqual(
            explicit_recent_remote_config_hash(base),
            explicit_recent_remote_config_hash(wider),
        )
        with self.assertRaisesRegex(ValueError, "age85PlusMax"):
            build_explicit_recent_remote_config(
                enabled=True,
                recent_tbi_target=0.01,
                remote_only_tbi_target=0.02,
                target_source="unit test source",
                review_status="reviewed_test_fixture",
                age85_plus_max=84,
            )

    def test_configuration_does_not_reuse_retired_or_frozen_identifiers(self) -> None:
        cfg = build_explicit_recent_remote_config(
            enabled=True,
            recent_tbi_target=0.01,
            remote_only_tbi_target=0.02,
            target_source="unit test source",
            review_status="reviewed_test_fixture",
        )
        flattened = json.dumps(cfg.as_dict(), sort_keys=True)

        self.assertNotIn("baselineRecentLTBIProportion", flattened)
        self.assertNotIn("continuous_markov_recent_remote", flattened)
        self.assertNotIn("matlab_v9_implicit_early_late", flattened)
        self.assertNotIn("sa_health_matlab_v9_compatibility_reference", flattened)
        with self.assertRaisesRegex(ValueError, "recentHazardShape"):
            explicit_recent_remote_config_from_dict(
                {
                    **cfg.as_dict(),
                    "recentHazardShape": "continuous_markov_recent_remote",
                }
            )

    def test_enabled_configuration_requires_target_source(self) -> None:
        with self.assertRaisesRegex(ValueError, "targetSource"):
            build_explicit_recent_remote_config(
                enabled=True,
                recent_tbi_target=0.01,
                remote_only_tbi_target=0.02,
                target_source="",
                review_status="reviewed_test_fixture",
            )


class ExplicitRecentRemoteProbabilityTests(unittest.TestCase):
    def test_probabilities_sum_to_one(self) -> None:
        probs = mutually_exclusive_state_probabilities([0, 3, 5, 20, 120], 0.02, 0.01)

        for values in zip(probs["recent"], probs["remote_only"], probs["uninfected"]):
            self.assertAlmostEqual(sum(values), 1.0, places=12)
            for value in values:
                self.assertGreaterEqual(value, 0.0)
                self.assertLessEqual(value, 1.0)

    def test_no_remote_infection_is_possible_at_age_five_or_younger(self) -> None:
        durations = exposure_durations_for_ages([0, 4.999, 5.0])

        self.assertEqual(durations.remote_years, (0.0, 0.0, 0.0))
        probs = mutually_exclusive_state_probabilities([0, 4.999, 5.0], 0.0, 10.0)
        self.assertEqual(probs["remote_only"], (0.0, 0.0, 0.0))

    def test_recent_exposure_duration_is_age_for_ages_under_five(self) -> None:
        durations = exposure_durations_for_ages([0, 1.5, 4.9])

        self.assertEqual(durations.recent_years, (0.0, 1.5, 4.9))

    def test_recent_exposure_duration_is_five_years_for_ages_over_five(self) -> None:
        durations = exposure_durations_for_ages([5, 6, 80])

        self.assertEqual(durations.recent_years, (5.0, 5.0, 5.0))

    def test_remote_exposure_accumulates_with_age_up_to_100_year_cap(self) -> None:
        durations = exposure_durations_for_ages([5, 50, 100, 120])

        self.assertEqual(durations.remote_years, (0.0, 45.0, 95.0, 95.0))

    def test_exact_age_window_boundaries(self) -> None:
        durations = exposure_durations_for_ages([0, 4.999, 5, 5.001, 100, 100.1])

        self.assertEqual(durations.recent_years[:3], (0.0, 4.999, 5.0))
        self.assertEqual(durations.remote_years[:3], (0.0, 0.0, 0.0))
        self.assertAlmostEqual(durations.recent_years[3], 5.0)
        self.assertAlmostEqual(durations.remote_years[3], 0.001)
        self.assertAlmostEqual(durations.remote_years[4], 95.0)
        self.assertAlmostEqual(durations.remote_years[5], 95.0)

    def test_recent_event_overrides_remote_in_effective_state(self) -> None:
        probs = mutually_exclusive_state_probabilities([60], 0.08, 0.03)

        recent = probs["recent"][0]
        remote_only = probs["remote_only"][0]
        prior_remote = probs["prior_remote_exposure"][0]
        recent_with_prior_remote = probs["recent_with_prior_remote"][0]
        self.assertGreater(recent_with_prior_remote, 0.0)
        self.assertGreater(prior_remote, remote_only)
        self.assertAlmostEqual(probs["total_tbi"][0], recent + remote_only)
        self.assertLess(probs["total_tbi"][0], recent + prior_remote)

    def test_age_interpretation_distinguishes_population_probability_from_tbi_fraction(self) -> None:
        ages = [10, 80]
        probs = mutually_exclusive_state_probabilities(ages, 0.02, 0.02)
        recent_population_probability = probs["recent"]
        recent_among_tbi = [
            recent / total
            for recent, total in zip(probs["recent"], probs["total_tbi"])
        ]

        self.assertAlmostEqual(
            recent_population_probability[0],
            recent_population_probability[1],
            places=12,
        )
        self.assertGreater(recent_among_tbi[0], recent_among_tbi[1])


class ExplicitRecentRemoteCalibrationTests(unittest.TestCase):
    ages = [2, 10, 30, 80]
    proportions = [0.10, 0.20, 0.40, 0.30]

    def test_zero_recent_and_zero_remote_targets_return_zero_hazards(self) -> None:
        result = calibrate_recent_remote_hazards(self.ages, self.proportions, 0.0, 0.0)

        self.assertEqual(result.fitted_recent_hazard, 0.0)
        self.assertEqual(result.fitted_remote_hazard, 0.0)
        self.assertEqual(result.achieved_total_tbi_prevalence, 0.0)

    def test_recent_only_targets_calibrate_correctly(self) -> None:
        result = calibrate_recent_remote_hazards(self.ages, self.proportions, 0.08, 0.0)

        self.assertAlmostEqual(result.achieved_recent_prevalence, 0.08, places=10)
        self.assertAlmostEqual(result.achieved_remote_only_prevalence, 0.0, places=12)
        self.assertGreater(result.fitted_recent_hazard, 0.0)
        self.assertEqual(result.fitted_remote_hazard, 0.0)

    def test_remote_only_targets_calibrate_correctly(self) -> None:
        result = calibrate_recent_remote_hazards(self.ages, self.proportions, 0.0, 0.10)

        self.assertAlmostEqual(result.achieved_recent_prevalence, 0.0, places=12)
        self.assertAlmostEqual(result.achieved_remote_only_prevalence, 0.10, places=10)
        self.assertEqual(result.fitted_recent_hazard, 0.0)
        self.assertGreater(result.fitted_remote_hazard, 0.0)

    def test_joint_targets_calibrate_correctly(self) -> None:
        result = calibrate_recent_remote_hazards(self.ages, self.proportions, 0.07, 0.12)

        self.assertEqual(result.convergence_status, "converged")
        self.assertAlmostEqual(result.achieved_recent_prevalence, 0.07, places=10)
        self.assertAlmostEqual(result.achieved_remote_only_prevalence, 0.12, places=10)
        self.assertAlmostEqual(
            result.achieved_total_tbi_prevalence,
            result.achieved_recent_prevalence + result.achieved_remote_only_prevalence,
            places=12,
        )
        self.assertAlmostEqual(
            result.achieved_uninfected_prevalence + result.achieved_total_tbi_prevalence,
            1.0,
            places=12,
        )

    def test_achieved_population_prevalences_match_feasible_targets(self) -> None:
        result = calibrate_recent_remote_hazards(self.ages, self.proportions, 0.15, 0.25)
        achieved = population_weighted_prevalences(
            self.ages,
            self.proportions,
            result.fitted_recent_hazard,
            result.fitted_remote_hazard,
        )

        self.assertAlmostEqual(achieved.recent, result.requested_recent_prevalence, places=10)
        self.assertAlmostEqual(achieved.remote_only, result.requested_remote_only_prevalence, places=10)

    def test_impossible_targets_are_rejected_rather_than_clipped(self) -> None:
        assessment = assess_calibration_feasibility([2, 3], [0.5, 0.5], 0.2, 0.01)
        self.assertFalse(assessment["isFeasible"])

        with self.assertRaises(CalibrationError):
            calibrate_recent_remote_hazards([2, 3], [0.5, 0.5], 0.2, 0.01)

    def test_empty_or_malformed_age_distributions_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one"):
            validate_age_distribution([], [])
        with self.assertRaisesRegex(ValueError, "same length"):
            validate_age_distribution([10], [0.5, 0.5])
        with self.assertRaisesRegex(ValueError, "non-negative"):
            validate_age_distribution([-1], [1.0])

    def test_age_proportions_must_be_non_negative_and_sum_to_one(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-negative"):
            validate_age_distribution([10, 20], [1.1, -0.1])
        with self.assertRaisesRegex(ValueError, "sum to one"):
            validate_age_distribution([10, 20], [0.6, 0.5])

    def test_calibration_is_deterministic_and_repeatable(self) -> None:
        first = calibrate_recent_remote_hazards(self.ages, self.proportions, 0.09, 0.11)
        second = calibrate_recent_remote_hazards(self.ages, self.proportions, 0.09, 0.11)

        self.assertEqual(first.as_dict(), second.as_dict())

    def test_calibration_serialization_round_trips(self) -> None:
        result = calibrate_recent_remote_hazards(self.ages, self.proportions, 0.09, 0.11)
        payload = json.loads(json.dumps(result.as_dict()))
        restored = calibration_result_from_dict(payload)

        self.assertEqual(restored.as_dict(), result.as_dict())
        self.assertEqual(payload["calibrationContractVersion"], CALIBRATION_CONTRACT_VERSION)
        self.assertEqual(payload["analysisBasis"], ANALYSIS_BASIS)

    def test_custom_window_parameters_flow_through_calibration(self) -> None:
        result = calibrate_recent_remote_hazards(
            self.ages,
            self.proportions,
            0.04,
            0.08,
            recent_window_years=3.0,
            remote_history_cap_years=50.0,
        )
        achieved = population_weighted_prevalences(
            self.ages,
            self.proportions,
            result.fitted_recent_hazard,
            result.fitted_remote_hazard,
            recent_window_years=3.0,
            remote_history_cap_years=50.0,
        )

        self.assertAlmostEqual(achieved.recent, 0.04, places=10)
        self.assertAlmostEqual(achieved.remote_only, 0.08, places=10)

    def test_age_support_sensitivity_reports_remote_hazard_changes(self) -> None:
        narrow = {"age85PlusMax": 89, "ages": list(range(85, 90)), "proportions": [0.2] * 5}
        wider = {"age85PlusMax": 95, "ages": list(range(85, 96)), "proportions": [1 / 11] * 11}

        rows = age_support_calibration_sensitivity(
            [narrow, wider],
            requested_recent_prevalence=0.05,
            requested_remote_only_prevalence=0.40,
        )

        self.assertEqual([row["age85PlusMax"] for row in rows], [89, 95])
        self.assertNotAlmostEqual(
            rows[0]["fittedRemoteHazard"],
            rows[1]["fittedRemoteHazard"],
            places=6,
        )
        self.assertAlmostEqual(rows[0]["achievedRemoteOnlyPrevalence"], 0.40, places=10)
        self.assertAlmostEqual(rows[1]["achievedRemoteOnlyPrevalence"], 0.40, places=10)

    def test_infection_time_quantile_uses_truncated_window(self) -> None:
        median_recent = infection_time_quantile(
            window_start_years_before_baseline=0.0,
            window_duration_years=5.0,
            hazard_per_year=0.03,
            quantile=0.5,
        )
        median_remote = infection_time_quantile(
            window_start_years_before_baseline=5.0,
            window_duration_years=30.0,
            hazard_per_year=0.03,
            quantile=0.5,
        )

        self.assertGreaterEqual(median_recent, 0.0)
        self.assertLessEqual(median_recent, 5.0)
        self.assertGreaterEqual(median_remote, 5.0)
        self.assertLessEqual(median_remote, 35.0)

    def test_timing_specification_records_unresolved_runner_mechanism(self) -> None:
        spec = infection_timing_specification()

        self.assertIn("most recent infection", spec["effectiveStateRule"])
        self.assertEqual(spec["bothRemoteAndRecent"]["effectiveState"], "recent")
        self.assertIn("does not currently distinguish", spec["unresolvedIntegrationQuestion"])


class ExplicitRecentRemoteAssignmentTests(unittest.TestCase):
    ages = [2, 10, 30, 80]
    proportions = [0.10, 0.20, 0.40, 0.30]

    def test_recent_state_overrides_remote_and_prior_remote_is_auditable(self) -> None:
        assignment = draw_recent_remote_states_for_ages([60], 1000.0, 1000.0, seed=123)

        self.assertEqual(assignment["effectiveStates"], [STATE_RECENT])
        self.assertEqual(assignment["priorRemoteExposure"], [True])
        self.assertEqual(assignment["priorRemotePlusRecent"], [True])
        self.assertLessEqual(assignment["timeSinceMostRecentInfection"][0], 5.0)

    def test_conditional_infection_times_stay_inside_state_windows(self) -> None:
        recent = draw_recent_remote_states_for_ages([3, 8], 1000.0, 0.0, seed=1)
        remote = draw_recent_remote_states_for_ages([10, 150], 0.0, 1000.0, seed=1)

        for idx, (age, time) in enumerate(
            zip(recent["ages"], recent["timeSinceMostRecentInfection"])
        ):
            self.assertEqual(recent["effectiveStates"][idx], STATE_RECENT)
            self.assertGreaterEqual(time, 0.0)
            self.assertLessEqual(time, min(5.0, age))
        for idx, (age, time) in enumerate(
            zip(remote["ages"], remote["timeSinceMostRecentInfection"])
        ):
            self.assertEqual(remote["effectiveStates"][idx], STATE_REMOTE_ONLY)
            self.assertGreaterEqual(time, 5.0)
            self.assertLessEqual(time, min(100.0, age))

    def test_inverse_cdf_sampling_agrees_with_analytic_quantiles(self) -> None:
        hazard = 0.25
        duration = 5.0
        q = 0.8
        event_probability = 1.0 - math.exp(-hazard * duration)
        expected = -math.log1p(-q * event_probability) / hazard

        observed = infection_time_quantile(
            window_start_years_before_baseline=0.0,
            window_duration_years=duration,
            hazard_per_year=hazard,
            quantile=q,
        )

        self.assertAlmostEqual(observed, expected, places=12)
        self.assertAlmostEqual(
            infection_time_quantile(
                window_start_years_before_baseline=5.0,
                window_duration_years=duration,
                hazard_per_year=hazard,
                quantile=0.0,
            ),
            5.0,
            places=12,
        )

    def test_recent_and_remote_time_sampling_handles_extreme_hazards(self) -> None:
        small_recent = infection_time_quantile(
            window_start_years_before_baseline=0.0,
            window_duration_years=5.0,
            hazard_per_year=1e-12,
            quantile=0.75,
        )
        large_recent = infection_time_quantile(
            window_start_years_before_baseline=0.0,
            window_duration_years=5.0,
            hazard_per_year=1e6,
            quantile=0.75,
        )
        small_remote = infection_time_quantile(
            window_start_years_before_baseline=5.0,
            window_duration_years=95.0,
            hazard_per_year=1e-12,
            quantile=0.25,
        )
        large_remote = infection_time_quantile(
            window_start_years_before_baseline=5.0,
            window_duration_years=95.0,
            hazard_per_year=1e6,
            quantile=0.25,
        )

        self.assertAlmostEqual(small_recent, 3.75, places=5)
        self.assertTrue(0.0 <= large_recent <= 5.0)
        self.assertAlmostEqual(small_remote, 5.0 + 23.75, places=4)
        self.assertTrue(5.0 <= large_remote <= 100.0)

    def test_expected_time_uses_truncated_distribution(self) -> None:
        mean_time = expected_time_since_most_recent_event(
            window_start_years_before_baseline=0.0,
            window_duration_years=5.0,
            hazard_per_year=1e-12,
        )

        self.assertAlmostEqual(mean_time, 2.5, places=5)

    def test_remaining_early_risk_duration_is_calculated_by_infection_time(self) -> None:
        assignment = draw_recent_remote_states_for_ages([10], 10.0, 0.0, seed=4)
        time_since = assignment["timeSinceMostRecentInfection"][0]

        self.assertEqual(assignment["effectiveStates"], [STATE_RECENT])
        self.assertAlmostEqual(
            assignment["remainingEarlyRiskYears"][0],
            max(0.0, assignment["recentWindowYears"] - time_since),
            places=12,
        )
        self.assertLess(assignment["remainingEarlyRiskYears"][0], 5.0)

    def test_non_default_recent_window_controls_remaining_risk(self) -> None:
        recent_window = 7.0
        assignment = draw_recent_remote_states_for_ages(
            [12],
            1000.0,
            0.0,
            seed=9,
            recent_window_years=recent_window,
            remote_history_cap_years=100.0,
        )
        time_since = assignment["timeSinceMostRecentInfection"][0]

        self.assertEqual(assignment["effectiveStates"], [STATE_RECENT])
        self.assertLessEqual(time_since, recent_window)
        self.assertAlmostEqual(
            assignment["remainingEarlyRiskYears"][0],
            max(0.0, recent_window - time_since),
            places=12,
        )

    def test_deterministic_non_default_recent_window_controls_remaining_risk(self) -> None:
        assignment = deterministic_recent_remote_assignment(
            [20],
            [1.0],
            1.0,
            0.0,
            recent_window_years=7.0,
            remote_history_cap_years=100.0,
        )
        row = assignment["ageSpecificStateDistributions"][0]

        self.assertEqual(assignment["earlyRiskPeriodYears"], 7.0)
        self.assertAlmostEqual(
            row["expectedRemainingEarlyRiskYearsAmongRecent"],
            7.0 - row["expectedTimeSinceMostRecentInfectionAmongRecent"],
            places=12,
        )

    def test_remote_only_receives_no_remaining_early_risk_duration(self) -> None:
        assignment = draw_recent_remote_states_for_ages([80], 0.0, 10.0, seed=5)

        self.assertEqual(assignment["effectiveStates"], [STATE_REMOTE_ONLY])
        self.assertEqual(assignment["remainingEarlyRiskYears"], [0.0])

    def test_deterministic_state_totals_match_calibration_targets(self) -> None:
        calibration = calibrate_recent_remote_hazards(
            self.ages, self.proportions, 0.07, 0.12
        )
        assignment = deterministic_recent_remote_assignment(
            self.ages,
            self.proportions,
            calibration.fitted_recent_hazard,
            calibration.fitted_remote_hazard,
            population_size=1000,
        )

        self.assertAlmostEqual(assignment["proportions"][STATE_RECENT], 0.07, places=10)
        self.assertAlmostEqual(
            assignment["proportions"][STATE_REMOTE_ONLY], 0.12, places=10
        )
        self.assertAlmostEqual(
            assignment["proportions"][STATE_UNINFECTED]
            + assignment["proportions"]["totalTBI"],
            1.0,
            places=12,
        )
        self.assertEqual(assignment["assignmentContractVersion"], ASSIGNMENT_CONTRACT_VERSION)

    def test_deterministic_age_specific_states_sum_correctly(self) -> None:
        assignment = deterministic_recent_remote_assignment(
            [0, 5, 80], [0.2, 0.3, 0.5], 0.02, 0.03
        )

        for row in assignment["ageSpecificStateDistributions"]:
            self.assertAlmostEqual(sum(row["stateProportions"].values()), 1.0, places=12)
            if row["ageYears"] <= 5.0:
                self.assertEqual(row["stateProportions"][STATE_REMOTE_ONLY], 0.0)

    def test_recent_population_probability_not_tbi_fraction_constant_by_age(self) -> None:
        assignment = deterministic_recent_remote_assignment(
            [10, 80], [0.5, 0.5], 0.02, 0.02
        )
        rows = assignment["ageSpecificStateDistributions"]

        self.assertAlmostEqual(
            rows[0]["recentPrevalenceInTotalAgeGroup"],
            rows[1]["recentPrevalenceInTotalAgeGroup"],
            places=12,
        )
        self.assertGreater(rows[0]["recentFractionAmongTBI"], rows[1]["recentFractionAmongTBI"])

    def test_deterministic_assignment_from_config_reconciles_targets(self) -> None:
        cfg = build_explicit_recent_remote_config(
            enabled=True,
            recent_tbi_target=0.06,
            remote_only_tbi_target=0.14,
            target_source="unit test source",
            review_status="reviewed_test_fixture",
        )

        assignment = deterministic_recent_remote_assignment_from_config(
            cfg, self.ages, self.proportions, population_size=5000
        )

        self.assertAlmostEqual(assignment["proportions"][STATE_RECENT], 0.06, places=10)
        self.assertAlmostEqual(
            assignment["proportions"][STATE_REMOTE_ONLY], 0.14, places=10
        )
        self.assertFalse(assignment["drawsUsed"])

    def test_stochastic_state_frequencies_converge_to_analytic_probabilities(self) -> None:
        calibration = calibrate_recent_remote_hazards(
            self.ages, self.proportions, 0.08, 0.18
        )
        assignment = stochastic_recent_remote_population_assignment(
            self.ages,
            self.proportions,
            calibration.fitted_recent_hazard,
            calibration.fitted_remote_hazard,
            population_size=80000,
            seed=44,
        )

        self.assertAlmostEqual(assignment["proportions"][STATE_RECENT], 0.08, delta=0.006)
        self.assertAlmostEqual(
            assignment["proportions"][STATE_REMOTE_ONLY], 0.18, delta=0.006
        )

    def test_fixed_seeds_reproduce_states_and_infection_times(self) -> None:
        first = stochastic_recent_remote_population_assignment(
            self.ages, self.proportions, 0.03, 0.02, population_size=200, seed=99
        )
        second = stochastic_recent_remote_population_assignment(
            self.ages, self.proportions, 0.03, 0.02, population_size=200, seed=99
        )
        third = stochastic_recent_remote_population_assignment(
            self.ages, self.proportions, 0.03, 0.02, population_size=200, seed=100
        )

        self.assertEqual(first["effectiveStates"], second["effectiveStates"])
        self.assertEqual(
            first["timeSinceMostRecentInfection"],
            second["timeSinceMostRecentInfection"],
        )
        self.assertNotEqual(first["effectiveStates"], third["effectiveStates"])

    def test_stochastic_assignment_does_not_consume_global_rng_state(self) -> None:
        np.random.seed(1234)
        stochastic_recent_remote_population_assignment(
            self.ages, self.proportions, 0.03, 0.02, population_size=100, seed=77
        )
        after = np.random.random(5)
        np.random.seed(1234)
        expected = np.random.random(5)

        np.testing.assert_allclose(after, expected)

    def test_disabled_pathway_causes_no_draws(self) -> None:
        class ExplodingRng:
            def random(self, *args, **kwargs):
                raise AssertionError("disabled pathway should not draw random numbers")

            def choice(self, *args, **kwargs):
                raise AssertionError("disabled pathway should not draw random numbers")

        cfg = build_explicit_recent_remote_config(
            enabled=False,
            recent_tbi_target=0.0,
            remote_only_tbi_target=0.0,
            review_status="not_enabled",
        )

        assignment = stochastic_recent_remote_population_assignment_from_config(
            cfg,
            self.ages,
            self.proportions,
            population_size=10,
            rng=ExplodingRng(),
        )

        self.assertFalse(assignment["enabled"])
        self.assertFalse(assignment["drawsUsed"])


class ExplicitRecentRemoteProgressionTests(unittest.TestCase):
    def valid_row(self, **overrides) -> dict:
        row = {
            "observationId": "obs-prog",
            "startYear": 2026,
            "endYear": 2026,
            "observedActiveTBCaseCount": 5,
            "populationDenominator": 10000,
            "personYears": 10000,
            "denominatorType": "census_population",
            "populationScope": "whole_population",
            "caseClassification": "incident_follow_up",
            "observationWindowMeaning": "follow_up_incident",
            "ascertainmentMethod": "combined",
            "activeTBClassification": "all_active_tb",
            "source": "unit test fixture",
            "reviewStatus": "unreviewed_test_fixture",
            "notes": "progression fixture",
            "uncertainty": {},
        }
        row.update(overrides)
        return row

    def ascertainment_kwargs(self, probability: float = 1.0) -> dict:
        return {
            "ascertainment_probability": probability,
            "ascertainment_source": (
                "unit test complete ascertainment assumption"
                if probability == 1.0
                else "unit test incomplete ascertainment assumption"
            ),
            "ascertainment_review_status": "reviewed_test_fixture",
        }

    def test_recent_cumulative_hazard_is_piecewise_correct(self) -> None:
        self.assertAlmostEqual(
            progression_cumulative_hazard(
                state=STATE_RECENT,
                horizon_years=1.0,
                early_hazard=0.2,
                remote_hazard=0.02,
                multiplier=2.0,
                remaining_early_risk_years=2.0,
            ),
            2.0 * 0.2 * 1.0,
            places=12,
        )
        self.assertAlmostEqual(
            progression_cumulative_hazard(
                state=STATE_RECENT,
                horizon_years=5.0,
                early_hazard=0.2,
                remote_hazard=0.02,
                multiplier=2.0,
                remaining_early_risk_years=2.0,
            ),
            2.0 * (0.2 * 2.0 + 0.02 * 3.0),
            places=12,
        )

    def test_remote_only_and_uninfected_progression_math(self) -> None:
        self.assertAlmostEqual(
            progression_cumulative_hazard(
                state=STATE_REMOTE_ONLY,
                horizon_years=10.0,
                early_hazard=0.2,
                remote_hazard=0.02,
                multiplier=3.0,
            ),
            3.0 * 0.02 * 10.0,
            places=12,
        )
        self.assertEqual(
            progression_probability(
                state=STATE_UNINFECTED,
                horizon_years=20.0,
                early_hazard=1.0,
                remote_hazard=1.0,
                multiplier=100.0,
            ),
            0.0,
        )

    def test_survival_plus_cumulative_incidence_equals_one(self) -> None:
        survival = progression_survival_probability(
            state=STATE_RECENT,
            horizon_years=6.0,
            early_hazard=0.15,
            remote_hazard=0.01,
            multiplier=1.4,
            remaining_early_risk_years=4.0,
        )
        incidence = progression_probability(
            state=STATE_RECENT,
            horizon_years=6.0,
            early_hazard=0.15,
            remote_hazard=0.01,
            multiplier=1.4,
            remaining_early_risk_years=4.0,
        )

        self.assertAlmostEqual(survival + incidence, 1.0, places=12)

    def test_progression_probability_is_monotonic(self) -> None:
        p1 = progression_probability(
            state=STATE_RECENT,
            horizon_years=1.0,
            early_hazard=0.05,
            remote_hazard=0.01,
            remaining_early_risk_years=2.0,
        )
        p2 = progression_probability(
            state=STATE_RECENT,
            horizon_years=2.0,
            early_hazard=0.05,
            remote_hazard=0.01,
            remaining_early_risk_years=2.0,
        )
        p3 = progression_probability(
            state=STATE_RECENT,
            horizon_years=2.0,
            early_hazard=0.10,
            remote_hazard=0.01,
            remaining_early_risk_years=2.0,
        )
        p4 = progression_probability(
            state=STATE_REMOTE_ONLY,
            horizon_years=2.0,
            early_hazard=0.10,
            remote_hazard=0.02,
        )
        p5 = progression_probability(
            state=STATE_REMOTE_ONLY,
            horizon_years=2.0,
            early_hazard=0.10,
            remote_hazard=0.04,
        )

        self.assertLess(p1, p2)
        self.assertLess(p2, p3)
        self.assertLess(p4, p5)

    def test_transition_from_early_to_late_cumulative_hazard_is_continuous(self) -> None:
        boundary = progression_cumulative_hazard(
            state=STATE_RECENT,
            horizon_years=2.0,
            early_hazard=0.2,
            remote_hazard=0.02,
            multiplier=1.5,
            remaining_early_risk_years=2.0,
        )
        left = progression_cumulative_hazard(
            state=STATE_RECENT,
            horizon_years=2.0 - 1e-9,
            early_hazard=0.2,
            remote_hazard=0.02,
            multiplier=1.5,
            remaining_early_risk_years=2.0,
        )
        right = progression_cumulative_hazard(
            state=STATE_RECENT,
            horizon_years=2.0 + 1e-9,
            early_hazard=0.2,
            remote_hazard=0.02,
            multiplier=1.5,
            remaining_early_risk_years=2.0,
        )

        self.assertAlmostEqual(boundary, 1.5 * 0.2 * 2.0, places=12)
        self.assertAlmostEqual(left, boundary, places=8)
        self.assertAlmostEqual(right, boundary, places=8)
        self.assertAlmostEqual(
            progression_piecewise_hazard(
                state=STATE_RECENT,
                time_years=1.999,
                early_hazard=0.2,
                remote_hazard=0.02,
                multiplier=1.5,
                remaining_early_risk_years=2.0,
            ),
            0.3,
            places=12,
        )
        self.assertAlmostEqual(
            progression_piecewise_hazard(
                state=STATE_RECENT,
                time_years=2.0,
                early_hazard=0.2,
                remote_hazard=0.02,
                multiplier=1.5,
                remaining_early_risk_years=2.0,
            ),
            0.03,
            places=12,
        )

    def test_remaining_risk_examples_use_recent_window(self) -> None:
        self.assertAlmostEqual(
            remaining_early_risk_years(4.9, recent_window_years=5.0),
            0.1,
            places=12,
        )
        self.assertAlmostEqual(
            remaining_early_risk_years(0.5, recent_window_years=5.0),
            4.5,
            places=12,
        )
        self.assertAlmostEqual(
            remaining_early_risk_years(0.5, recent_window_years=7.0),
            6.5,
            places=12,
        )

    def test_multipliers_scale_hazard_but_probability_saturates(self) -> None:
        base_hazard = progression_cumulative_hazard(
            state=STATE_REMOTE_ONLY,
            horizon_years=20.0,
            early_hazard=0.2,
            remote_hazard=0.02,
            multiplier=1.0,
        )
        high_hazard = progression_cumulative_hazard(
            state=STATE_REMOTE_ONLY,
            horizon_years=20.0,
            early_hazard=0.2,
            remote_hazard=0.02,
            multiplier=10.0,
        )
        base_probability = progression_probability(
            state=STATE_REMOTE_ONLY,
            horizon_years=20.0,
            early_hazard=0.2,
            remote_hazard=0.02,
            multiplier=1.0,
        )
        high_probability = progression_probability(
            state=STATE_REMOTE_ONLY,
            horizon_years=20.0,
            early_hazard=0.2,
            remote_hazard=0.02,
            multiplier=10.0,
        )

        self.assertAlmostEqual(high_hazard, 10.0 * base_hazard, places=12)
        self.assertGreater(high_probability, base_probability)
        self.assertLess(high_probability, 1.0)
        self.assertLess(high_probability, 10.0 * base_probability)

    def test_invalid_negative_hazards_and_multipliers_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "early_hazard"):
            progression_probability(
                state=STATE_RECENT,
                horizon_years=1,
                early_hazard=-0.1,
                remote_hazard=0.01,
            )
        with self.assertRaisesRegex(ValueError, "multiplier"):
            progression_probability(
                state=STATE_RECENT,
                horizon_years=1,
                early_hazard=0.1,
                remote_hazard=0.01,
                multiplier=-1,
            )

    def test_expected_count_aggregation_and_ascertainment(self) -> None:
        strata = [
            {"state": STATE_RECENT, "weight": 100.0, "remainingEarlyRiskYears": 2.0},
            {"state": STATE_UNINFECTED, "weight": 900.0},
        ]
        full = expected_progression_events(
            strata,
            horizon_years=1.0,
            early_hazard=0.1,
            remote_hazard=0.01,
        )
        half = expected_progression_events(
            strata,
            horizon_years=1.0,
            early_hazard=0.1,
            remote_hazard=0.01,
            default_ascertainment_probability=0.5,
        )

        expected = 100.0 * (1.0 - math.exp(-0.1))
        self.assertAlmostEqual(full["expectedCases"], expected, places=12)
        self.assertAlmostEqual(half["expectedCases"], expected * 0.5, places=12)

    def test_progression_time_quantile_and_sampling_specification(self) -> None:
        q = 1.0 - math.exp(-0.1)
        self.assertAlmostEqual(
            progression_time_quantile(
                state=STATE_RECENT,
                quantile=q,
                early_hazard=0.1,
                remote_hazard=0.01,
                multiplier=1.0,
                remaining_early_risk_years=2.0,
            ),
            1.0,
            places=12,
        )
        spec = progression_sampling_specification()
        self.assertEqual(spec["contractVersion"], PROGRESSION_CONTRACT_VERSION)
        self.assertIn("explicit RNG", spec["rngRequirement"])

    def test_active_tb_target_classification(self) -> None:
        baseline = classify_active_tb_observation_target(
            self.valid_row(
                caseClassification="prevalent_baseline",
                observationWindowMeaning="baseline_prevalent",
            ),
            model_baseline_year=2026,
        )
        screen = classify_active_tb_observation_target(
            self.valid_row(
                caseClassification="screen_detected",
                observationWindowMeaning="screen_detected_prevalent",
                ascertainmentMethod="active_screening",
            ),
            model_baseline_year=2026,
        )
        prospective = classify_active_tb_observation_target(
            self.valid_row(startYear=2026, endYear=2027),
            model_baseline_year=2026,
        )
        retrospective = classify_active_tb_observation_target(
            self.valid_row(startYear=2023, endYear=2024, ascertainmentMethod="passive_notification"),
            model_baseline_year=2026,
        )

        self.assertEqual(baseline["targetType"], TARGET_BASELINE_PREVALENCE)
        self.assertEqual(screen["targetType"], TARGET_SCREEN_DETECTED)
        self.assertEqual(prospective["targetType"], TARGET_PROSPECTIVE_INCIDENT)
        self.assertEqual(retrospective["targetType"], TARGET_RETROSPECTIVE_NOTIFICATION)

    def test_prospective_incident_expected_cases_retain_observation_components(self) -> None:
        result = expected_prospective_incident_cases_for_observation(
            self.valid_row(startYear=2026, endYear=2026, observedActiveTBCaseCount=8),
            [{"state": STATE_REMOTE_ONLY, "weight": 100.0, "multiplier": 1.0}],
            model_baseline_year=2026,
            early_hazard=0.2,
            remote_hazard=0.01,
            **self.ascertainment_kwargs(0.8),
        )

        self.assertEqual(result["observedActiveTBCaseCount"], 8.0)
        self.assertEqual(result["populationDenominator"], 10000.0)
        self.assertEqual(result["sourcePopulationDenominator"], 10000.0)
        self.assertEqual(result["personYears"], 10000.0)
        self.assertAlmostEqual(
            result["expectedActiveTBCaseCount"],
            100.0 * (1.0 - math.exp(-0.01)) * 0.8,
            places=12,
        )

    def test_retrospective_targets_are_not_silently_accepted(self) -> None:
        with self.assertRaisesRegex(ValueError, "prospective incident"):
            expected_prospective_incident_cases_for_observation(
                self.valid_row(
                    startYear=2023,
                    endYear=2024,
                    ascertainmentMethod="passive_notification",
                ),
                [{"state": STATE_REMOTE_ONLY, "weight": 100.0}],
                model_baseline_year=2026,
                early_hazard=0.2,
                remote_hazard=0.01,
                **self.ascertainment_kwargs(),
            )

    def test_one_active_tb_target_is_insufficient_for_two_free_hazards(self) -> None:
        assessment = progression_identifiability_assessment(
            [self.valid_row(startYear=2026, endYear=2026)],
            model_baseline_year=2026,
        )

        self.assertFalse(assessment["isIdentifiedByAggregateTargets"])
        self.assertIn("cannot identify both", assessment["diagnosticMessages"][0])

    def test_baseline_active_tb_sequence_is_separate_from_tbi(self) -> None:
        spec = baseline_active_tb_state_sequence_specification()

        self.assertEqual(
            spec["mutuallyExclusiveBaselineSequence"][0],
            "baseline_prevalent_active_tb",
        )
        self.assertIn("not counted simultaneously", spec["rule"])
        self.assertIn(STATE_RECENT, spec["futureIntegrationStates"])

    def test_policy_options_are_explicit_and_no_hidden_default_is_chosen(self) -> None:
        policies = {item["policyId"] for item in progression_calibration_policy_options()}

        self.assertEqual(
            policies,
            {
                POLICY_BOTH_HAZARDS_SUPPLIED,
                POLICY_ONE_HAZARD_SUPPLIED,
                POLICY_FIXED_EARLY_LATE_RATIO,
                POLICY_JOINT_LIKELIHOOD,
                POLICY_EXTERNAL_VALIDATION_ONLY,
            },
        )

    def test_worked_diagnostics_and_cost_invariance_are_reviewable(self) -> None:
        rows = worked_progression_diagnostic_table(
            early_hazard=0.04,
            remote_hazard=0.004,
            multipliers=(1.0, 4.0),
            horizons=(1.0, 20.0),
        )
        first = next(
            row
            for row in rows
            if row["scenario"] == "recent_s_0_5"
            and row["multiplier"] == 1.0
            and row["horizonYears"] == 1.0
        )
        low_cost = synthetic_natural_history_cost_invariance_example(
            cost_parameter=1.0,
            intervention_parameter=0.0,
            early_hazard=0.04,
            remote_hazard=0.004,
        )
        high_cost = synthetic_natural_history_cost_invariance_example(
            cost_parameter=999999.0,
            intervention_parameter=1.0,
            early_hazard=0.04,
            remote_hazard=0.004,
        )

        self.assertAlmostEqual(
            first["progressionProbability"],
            1.0 - math.exp(-0.04),
            places=12,
        )
        self.assertAlmostEqual(
            low_cost["expectedProgressionEvents"],
            high_cost["expectedProgressionEvents"],
            places=12,
        )


class ExplicitProgressionCalibrationPolicyTests(unittest.TestCase):
    def valid_row(self, **overrides) -> dict:
        row = {
            "observationId": "obs-policy",
            "startYear": 2026,
            "endYear": 2026,
            "observedActiveTBCaseCount": 2,
            "populationDenominator": 1000,
            "personYears": 1000,
            "denominatorType": "census_population",
            "populationScope": "whole_population",
            "caseClassification": "incident_follow_up",
            "observationWindowMeaning": "follow_up_incident",
            "ascertainmentMethod": "combined",
            "activeTBClassification": "all_active_tb",
            "source": "unit test fixture",
            "reviewStatus": "unreviewed_test_fixture",
            "notes": "policy fixture",
            "uncertainty": {},
        }
        row.update(overrides)
        return row

    def ascertainment_kwargs(self, probability: float = 1.0) -> dict:
        return {
            "ascertainment_probability": probability,
            "ascertainment_source": (
                "unit test complete ascertainment assumption"
                if probability == 1.0
                else "unit test incomplete ascertainment assumption"
            ),
            "ascertainment_review_status": "reviewed_test_fixture",
        }

    def strata(self):
        return [
            {"state": STATE_RECENT, "weight": 100.0, "remainingEarlyRiskYears": 4.5},
            {"state": STATE_REMOTE_ONLY, "weight": 200.0, "remainingEarlyRiskYears": 0.0},
            {"state": STATE_UNINFECTED, "weight": 700.0, "remainingEarlyRiskYears": 0.0},
        ]

    def test_progression_calibration_policy_serialization_and_hashing(self) -> None:
        policy = build_progression_calibration_policy(
            policy_id=POLICY_EXTERNAL_HAZARDS,
            early_hazard=0.02,
            remote_hazard=0.002,
            source="unit test source",
            reference_population="unit test population",
            review_status="reviewed_test_fixture",
            notes="external hazard policy",
            observation_model_id=OBSERVATION_MODEL_POISSON,
        )
        payload = json.loads(progression_calibration_policy_json(policy))
        shuffled = dict(reversed(list(payload.items())))

        self.assertEqual(
            payload["contractVersion"],
            PROGRESSION_CALIBRATION_CONTRACT_VERSION,
        )
        self.assertEqual(
            progression_calibration_policy_hash(payload),
            progression_calibration_policy_hash(shuffled),
        )
        self.assertEqual(
            {item["policyId"] for item in progression_calibration_policy_contract_options()},
            {
                POLICY_EXTERNAL_HAZARDS,
                POLICY_FIXED_RATIO_FIT_SCALE,
                POLICY_VALIDATION_ONLY,
                POLICY_JOINT_EARLY_REMOTE_HAZARDS,
            },
        )

    def test_external_hazard_policy_performs_no_fitting(self) -> None:
        policy = build_progression_calibration_policy(
            policy_id=POLICY_EXTERNAL_HAZARDS,
            early_hazard=0.02,
            remote_hazard=0.002,
            source="unit test source",
            reference_population="unit test population",
            review_status="reviewed_test_fixture",
        )
        result = evaluate_external_hazard_policy(
            self.valid_row(),
            self.strata(),
            policy,
            model_baseline_year=2026,
            **self.ascertainment_kwargs(),
        )

        self.assertFalse(result["fittingPerformed"])
        self.assertEqual(result["earlyHazard"], 0.02)
        self.assertEqual(result["remoteHazard"], 0.002)

    def test_fixed_ratio_policy_preserves_ratio_and_reproduces_target(self) -> None:
        result = fit_fixed_ratio_progression_scale(
            self.valid_row(),
            self.strata(),
            model_baseline_year=2026,
            early_to_remote_ratio=10.0,
            ratio_source="unit test ratio",
            ratio_review_status="reviewed_test_fixture",
            ratio_provenance="synthetic fixture",
            **self.ascertainment_kwargs(),
        )

        self.assertEqual(result["policyId"], POLICY_FIXED_RATIO_FIT_SCALE)
        self.assertAlmostEqual(
            result["derivedEarlyHazard"] / result["fittedRemoteHazard"],
            10.0,
            places=12,
        )
        self.assertAlmostEqual(result["achievedExpectedCases"], 2.0, places=8)
        self.assertEqual(result["convergenceStatus"], "converged")

    def test_zero_target_returns_zero_hazards(self) -> None:
        result = fit_fixed_ratio_progression_scale(
            self.valid_row(observedActiveTBCaseCount=0),
            self.strata(),
            model_baseline_year=2026,
            early_to_remote_ratio=10.0,
            ratio_source="unit test ratio",
            ratio_review_status="reviewed_test_fixture",
            ratio_provenance="synthetic fixture",
            **self.ascertainment_kwargs(),
        )

        self.assertEqual(result["fittedRemoteHazard"], 0.0)
        self.assertEqual(result["derivedEarlyHazard"], 0.0)
        self.assertEqual(result["achievedExpectedCases"], 0.0)

    def test_missing_ascertainment_blocks_fixed_ratio_policy(self) -> None:
        with self.assertRaises(ProgressionCalibrationError) as caught:
            fit_fixed_ratio_progression_scale(
                self.valid_row(),
                self.strata(),
                model_baseline_year=2026,
                early_to_remote_ratio=10.0,
                ratio_source="unit test ratio",
                ratio_review_status="reviewed_test_fixture",
                ratio_provenance="synthetic fixture",
            )

        self.assertEqual(
            caught.exception.diagnostics["feasibilityStatus"],
            "ineligible_ascertainment_identifiability",
        )

    def test_impossible_targets_are_rejected_rather_than_clipped(self) -> None:
        with self.assertRaises(ProgressionCalibrationError) as caught:
            fit_fixed_ratio_progression_scale(
                self.valid_row(observedActiveTBCaseCount=500),
                self.strata(),
                model_baseline_year=2026,
                early_to_remote_ratio=10.0,
                ratio_source="unit test ratio",
                ratio_review_status="reviewed_test_fixture",
                ratio_provenance="synthetic fixture",
                **self.ascertainment_kwargs(),
            )

        self.assertEqual(
            caught.exception.diagnostics["feasibilityStatus"],
            "infeasible_above_achievable_range",
        )

    def test_ineligible_targets_cannot_be_used_for_fitting(self) -> None:
        retrospective = self.valid_row(startYear=2023, endYear=2024)
        baseline = self.valid_row(
            caseClassification="prevalent_baseline",
            observationWindowMeaning="baseline_prevalent",
        )

        for target in (retrospective, baseline):
            with self.assertRaises(ProgressionCalibrationError):
                fit_fixed_ratio_progression_scale(
                    target,
                    self.strata(),
                    model_baseline_year=2026,
                    early_to_remote_ratio=10.0,
                    ratio_source="unit test ratio",
                    ratio_review_status="reviewed_test_fixture",
                    ratio_provenance="synthetic fixture",
                    **self.ascertainment_kwargs(),
                )

    def test_observation_eligibility_reports_reasons(self) -> None:
        eligible = assess_observation_prospective_calibration_eligibility(
            self.valid_row(),
            model_baseline_year=2026,
            strata=self.strata(),
            ascertainment_probability=0.8,
            ascertainment_source="unit test incomplete ascertainment assumption",
            ascertainment_review_status="reviewed_test_fixture",
        )
        retrospective = assess_observation_prospective_calibration_eligibility(
            self.valid_row(startYear=2023, endYear=2024),
            model_baseline_year=2026,
            strata=self.strata(),
            ascertainment_probability=1.0,
            ascertainment_source="unit test complete ascertainment assumption",
            ascertainment_review_status="reviewed_test_fixture",
        )

        self.assertTrue(eligible["eligible"])
        self.assertIn("Incomplete ascertainment", eligible["warnings"][0])
        self.assertFalse(retrospective["eligible"])
        self.assertIn("before baseline", " ".join(retrospective["reasons"]))

    def test_source_denominator_is_preserved_and_at_risk_population_is_separate(self) -> None:
        strata = [
            {"state": STATE_BASELINE_ACTIVE_TB, "weight": 2.0},
            {"state": STATE_RECENT, "weight": 100.0, "remainingEarlyRiskYears": 4.5},
            {"state": STATE_REMOTE_ONLY, "weight": 1398.0},
        ]
        result = expected_prospective_incident_cases_for_observation(
            self.valid_row(
                observedActiveTBCaseCount=10,
                populationDenominator=1500,
                personYears=None,
                numeratorIncludesBaselineActiveTB=False,
            ),
            strata,
            model_baseline_year=2026,
            early_hazard=0.02,
            remote_hazard=0.002,
            **self.ascertainment_kwargs(),
        )

        self.assertEqual(result["sourcePopulationDenominator"], 1500.0)
        self.assertEqual(result["populationDenominator"], 1500.0)
        self.assertEqual(result["baselineActiveTBCount"], 2.0)
        self.assertEqual(result["prospectiveAtRiskPopulation"], 1498.0)
        self.assertEqual(result["tbiEligiblePopulation"], 1498.0)
        self.assertEqual(result["personTimeAtRisk"], 1498.0)
        self.assertTrue(result["denominatorSummary"]["sourceDenominatorPreserved"])

    def test_ambiguous_baseline_numerator_composition_blocks_fitting(self) -> None:
        with self.assertRaises(ProgressionCalibrationError) as caught:
            fit_fixed_ratio_progression_scale(
                self.valid_row(populationDenominator=1500, observedActiveTBCaseCount=10),
                [
                    {"state": STATE_BASELINE_ACTIVE_TB, "weight": 2.0},
                    {"state": STATE_REMOTE_ONLY, "weight": 1498.0},
                ],
                model_baseline_year=2026,
                early_to_remote_ratio=10.0,
                ratio_source="unit test ratio",
                ratio_review_status="reviewed_test_fixture",
                ratio_provenance="synthetic fixture",
                **self.ascertainment_kwargs(),
            )

        reasons = " ".join(caught.exception.diagnostics["eligibility"]["reasons"])
        self.assertIn("numerator composition is ambiguous", reasons)

    def test_validation_only_mode_never_modifies_hazards(self) -> None:
        policy = build_progression_calibration_policy(
            policy_id=POLICY_VALIDATION_ONLY,
            source="unit test source",
            reference_population="unit test population",
            review_status="reviewed_test_fixture",
        )
        result = evaluate_validation_only_policy(
            self.valid_row(startYear=2023, endYear=2024),
            self.strata(),
            policy,
            model_baseline_year=2026,
            early_hazard=0.02,
            remote_hazard=0.002,
        )

        self.assertFalse(result["fittingPerformed"])
        self.assertEqual(result["earlyHazard"], 0.02)
        self.assertEqual(result["remoteHazard"], 0.002)
        self.assertEqual(result["targetType"], TARGET_RETROSPECTIVE_NOTIFICATION)

    def test_ascertainment_enters_expected_cases_explicitly(self) -> None:
        full = expected_prospective_incident_cases_for_observation(
            self.valid_row(),
            self.strata(),
            model_baseline_year=2026,
            early_hazard=0.02,
            remote_hazard=0.002,
            **self.ascertainment_kwargs(),
        )
        half = expected_prospective_incident_cases_for_observation(
            self.valid_row(),
            self.strata(),
            model_baseline_year=2026,
            early_hazard=0.02,
            remote_hazard=0.002,
            **self.ascertainment_kwargs(0.5),
        )

        self.assertAlmostEqual(
            half["expectedActiveTBCaseCount"],
            full["expectedActiveTBCaseCount"] * 0.5,
            places=12,
        )

    def test_ascertainment_and_progression_scale_are_confounding(self) -> None:
        full_q = fit_fixed_ratio_progression_scale(
            self.valid_row(),
            self.strata(),
            model_baseline_year=2026,
            early_to_remote_ratio=10.0,
            ratio_source="unit test ratio",
            ratio_review_status="reviewed_test_fixture",
            ratio_provenance="synthetic fixture",
            **self.ascertainment_kwargs(),
        )
        half_q = fit_fixed_ratio_progression_scale(
            self.valid_row(),
            self.strata(),
            model_baseline_year=2026,
            early_to_remote_ratio=10.0,
            ratio_source="unit test ratio",
            ratio_review_status="reviewed_test_fixture",
            ratio_provenance="synthetic fixture",
            **self.ascertainment_kwargs(0.5),
        )
        diagnostic = ascertainment_progression_scale_confounding_diagnostic(
            observed_cases=2,
            expected_true_cases_at_scale=2,
            ascertainment_probability=0.5,
        )

        self.assertGreater(half_q["fittedRemoteHazard"], full_q["fittedRemoteHazard"])
        self.assertAlmostEqual(full_q["achievedExpectedCases"], 2.0, places=8)
        self.assertAlmostEqual(half_q["achievedExpectedCases"], 2.0, places=8)
        self.assertEqual(diagnostic["requiredProgressionScaleMultiplier"], 2.0)

    def test_binomial_and_poisson_likelihoods_are_analytic_and_explicit(self) -> None:
        binomial = active_tb_observation_log_likelihood(
            observed_cases=2,
            expected_cases=10,
            denominator=100,
            observation_model_id=OBSERVATION_MODEL_BINOMIAL,
        )
        poisson = active_tb_observation_log_likelihood(
            observed_cases=2,
            expected_cases=3,
            denominator=100,
            observation_model_id=OBSERVATION_MODEL_POISSON,
        )
        expected_binomial = (
            math.log(math.comb(100, 2))
            + 2 * math.log(0.1)
            + 98 * math.log(0.9)
        )
        expected_poisson = 2 * math.log(3) - 3 - math.lgamma(3)

        self.assertAlmostEqual(binomial["logLikelihood"], expected_binomial, places=12)
        self.assertAlmostEqual(poisson["logLikelihood"], expected_poisson, places=12)
        with self.assertRaisesRegex(ValueError, "observation_model_id"):
            active_tb_observation_log_likelihood(
                observed_cases=2,
                expected_cases=3,
                denominator=100,
                observation_model_id="auto",
            )

    def test_identifiability_enforcement(self) -> None:
        one_target = assess_progression_policy_identifiability(
            POLICY_JOINT_EARLY_REMOTE_HAZARDS,
            [self.valid_row()],
            model_baseline_year=2026,
            strata_by_target=[self.strata()],
            **self.ascertainment_kwargs(),
        )
        fixed_ratio = assess_progression_policy_identifiability(
            POLICY_FIXED_RATIO_FIT_SCALE,
            [self.valid_row()],
            model_baseline_year=2026,
            strata_by_target=[self.strata()],
            early_to_remote_ratio=10.0,
            **self.ascertainment_kwargs(),
        )
        duplicate = assess_progression_policy_identifiability(
            POLICY_JOINT_EARLY_REMOTE_HAZARDS,
            [self.valid_row(observationId="a"), self.valid_row(observationId="b")],
            model_baseline_year=2026,
            strata_by_target=[self.strata(), self.strata()],
            **self.ascertainment_kwargs(),
        )

        self.assertFalse(one_target["identified"])
        self.assertTrue(fixed_ratio["identified"])
        self.assertEqual(fixed_ratio["identifiabilityStatus"], "identified_one_dimensional")
        self.assertFalse(duplicate["identified"])
        self.assertEqual(duplicate["sensitivityRank"], 1)

    def test_near_collinear_targets_fail_practical_identifiability(self) -> None:
        near_a = [{"state": STATE_RECENT, "weight": 100.0, "remainingEarlyRiskYears": 1.0}]
        near_b = [
            {"state": STATE_RECENT, "weight": 100.0, "remainingEarlyRiskYears": 1.0},
            {"state": STATE_REMOTE_ONLY, "weight": 0.1},
        ]
        assessment = assess_progression_policy_identifiability(
            POLICY_JOINT_EARLY_REMOTE_HAZARDS,
            [self.valid_row(observationId="near-a"), self.valid_row(observationId="near-b")],
            model_baseline_year=2026,
            strata_by_target=[near_a, near_b],
            **self.ascertainment_kwargs(),
        )

        self.assertEqual(assessment["sensitivityRank"], 2)
        self.assertFalse(assessment["identified"])
        self.assertLess(assessment["compositionContrast"], 0.05)

    def test_target_sensitivity_vectors_can_support_joint_identifiability(self) -> None:
        recent_heavy = [{"state": STATE_RECENT, "weight": 100.0, "remainingEarlyRiskYears": 1.0}]
        remote_heavy = [{"state": STATE_REMOTE_ONLY, "weight": 100.0, "remainingEarlyRiskYears": 0.0}]

        self.assertEqual(
            target_progression_sensitivity_vector(recent_heavy, horizon_years=1.0),
            (100.0, 0.0),
        )
        self.assertEqual(
            target_progression_sensitivity_vector(remote_heavy, horizon_years=1.0),
            (0.0, 100.0),
        )
        identified = assess_progression_policy_identifiability(
            POLICY_JOINT_EARLY_REMOTE_HAZARDS,
            [self.valid_row(observationId="a"), self.valid_row(observationId="b")],
            model_baseline_year=2026,
            strata_by_target=[recent_heavy, remote_heavy],
            **self.ascertainment_kwargs(),
        )
        self.assertTrue(identified["identified"])
        self.assertGreaterEqual(identified["compositionContrast"], 0.05)

    def test_risk_multiplier_policies(self) -> None:
        none = risk_factor_multiplier_for_row(
            {"diabetes": True, "smoking": True},
            build_risk_factor_application_policy(policy_id=RISK_POLICY_NONE),
        )
        reviewed = risk_factor_multiplier_for_row(
            {"diabetes": True, "smoking": True},
            build_risk_factor_application_policy(
                policy_id=RISK_POLICY_REVIEWED_HAZARD_MULTIPLIERS,
                reviewed_factor_names=("diabetes",),
            ),
        )
        legacy = risk_factor_multiplier_for_row(
            {"contact": True, "renal": True, "diabetes": True},
            build_risk_factor_application_policy(
                policy_id=RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC,
            ),
        )

        self.assertEqual(none["combinedMultiplier"], 1.0)
        self.assertEqual(reviewed["combinedMultiplier"], 3.0)
        self.assertTrue(any("diagnostic only" in item for item in legacy["warnings"]))
        self.assertAlmostEqual(legacy["combinedMultiplier"], 5.0 * 3.6 * 3.0)

    def test_legacy_multiplier_policy_cannot_be_production_default(self) -> None:
        with self.assertRaisesRegex(ValueError, "production/default"):
            build_risk_factor_application_policy(
                policy_id=RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC,
                production_default=True,
            )

    def test_large_joint_multipliers_generate_review_status(self) -> None:
        row = {
            "MJ": True,
            "contact": True,
            "renal": True,
            "diabetes": True,
            "smoking": True,
            "cld": True,
            "alcohol": True,
        }
        result = risk_factor_multiplier_for_row(
            row,
            build_risk_factor_application_policy(
                policy_id=RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC,
                review_warning_threshold=10.0,
                review_block_threshold=100.0,
            ),
        )

        self.assertEqual(result["combinedMultiplier"], 2916.0)
        self.assertEqual(result["reviewStatus"], "blocking_review_required")

    def test_calibration_reports_high_risk_concentration(self) -> None:
        high_risk_strata = [
            {
                "state": STATE_RECENT,
                "weight": 1.0,
                "remainingEarlyRiskYears": 4.5,
                "contact": True,
                "renal": True,
                "diabetes": True,
            },
            {"state": STATE_REMOTE_ONLY, "weight": 99.0, "remainingEarlyRiskYears": 0.0},
        ]
        scenarios = fixed_ratio_calibration_multiplier_sensitivity(
            self.valid_row(populationDenominator=100, observedActiveTBCaseCount=1),
            [
                {
                    "scenarioId": "none",
                    "riskPolicy": build_risk_factor_application_policy(policy_id=RISK_POLICY_NONE),
                    "strata": high_risk_strata,
                },
                {
                    "scenarioId": "legacy",
                    "riskPolicy": build_risk_factor_application_policy(
                        policy_id=RISK_POLICY_LEGACY_OR_AS_HAZARD_DIAGNOSTIC
                    ),
                    "strata": high_risk_strata,
                },
            ],
            model_baseline_year=2026,
            early_to_remote_ratio=10.0,
            ratio_source="unit test ratio",
            ratio_review_status="reviewed_test_fixture",
            ratio_provenance="synthetic fixture",
            **self.ascertainment_kwargs(),
        )

        none, legacy = scenarios
        self.assertLess(legacy["fittedRemoteHazard"], none["fittedRemoteHazard"])
        self.assertGreater(
            legacy["caseConcentration"]["shareExpectedCasesTop1PercentByWeight"],
            none["caseConcentration"]["shareExpectedCasesTop1PercentByWeight"],
        )

    def test_competing_mortality_constant_hazard_formula_is_correct(self) -> None:
        h_tb = 0.02
        h_death = 0.03
        horizon = 10.0
        result = competing_risk_progression_probability(
            state=STATE_REMOTE_ONLY,
            horizon_years=horizon,
            early_hazard=0.2,
            remote_hazard=h_tb,
            competing_mortality_hazard=h_death,
        )
        expected = h_tb / (h_tb + h_death) * (1.0 - math.exp(-(h_tb + h_death) * horizon))

        self.assertAlmostEqual(result["cumulativeIncidence"], expected, places=12)

    def test_piecewise_competing_mortality_integrates_across_early_late_boundary(self) -> None:
        result = competing_risk_progression_probability(
            state=STATE_RECENT,
            horizon_years=5.0,
            early_hazard=0.1,
            remote_hazard=0.01,
            remaining_early_risk_years=2.0,
            competing_mortality_hazard=0.03,
        )
        first = 0.1 / 0.13 * (1.0 - math.exp(-0.13 * 2.0))
        second = math.exp(-0.13 * 2.0) * 0.01 / 0.04 * (1.0 - math.exp(-0.04 * 3.0))

        self.assertAlmostEqual(result["cumulativeIncidence"], first + second, places=12)

    def test_horizon_level_survival_multiplication_is_not_used(self) -> None:
        h_tb = 0.02
        h_death = 0.03
        horizon = 10.0
        result = competing_risk_progression_probability(
            state=STATE_REMOTE_ONLY,
            horizon_years=horizon,
            early_hazard=0.2,
            remote_hazard=h_tb,
            competing_mortality_hazard=h_death,
        )
        horizon_product = (1.0 - math.exp(-h_tb * horizon)) * math.exp(-h_death * horizon)

        self.assertNotAlmostEqual(result["cumulativeIncidence"], horizon_product, places=4)
        with self.assertRaisesRegex(ValueError, "scalar horizon survival"):
            expected_progression_events(
                self.strata(),
                horizon_years=20.0,
                early_hazard=0.02,
                remote_hazard=0.002,
                competing_survival_probability=0.8,
            )

    def test_invalid_survival_functions_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "begin at one"):
            competing_risk_progression_probability(
                state=STATE_REMOTE_ONLY,
                horizon_years=2.0,
                early_hazard=0.1,
                remote_hazard=0.01,
                competing_survival_curve={0.0: 0.9, 2.0: 0.8},
            )
        with self.assertRaisesRegex(ValueError, "non-increasing"):
            competing_risk_progression_probability(
                state=STATE_REMOTE_ONLY,
                horizon_years=2.0,
                early_hazard=0.1,
                remote_hazard=0.01,
                competing_survival_curve={0.0: 1.0, 1.0: 0.8, 2.0: 0.9},
            )

    def test_optional_mortality_lowers_or_preserves_expected_cases(self) -> None:
        no_mortality = expected_progression_events(
            self.strata(),
            horizon_years=20.0,
            early_hazard=0.02,
            remote_hazard=0.002,
        )
        with_mortality = expected_progression_events(
            self.strata(),
            horizon_years=20.0,
            early_hazard=0.02,
            remote_hazard=0.002,
            competing_mortality_hazard=0.02,
        )
        curve_survival = expected_progression_events(
            self.strata(),
            horizon_years=20.0,
            early_hazard=0.02,
            remote_hazard=0.002,
            competing_survival_curve={0.0: 1.0, 10.0: 0.9, 20.0: 0.8},
        )

        self.assertEqual(no_mortality["competingMortality"], "not_modelled")
        self.assertIn("not modelled", no_mortality["diagnosticMessages"][0])
        self.assertLessEqual(with_mortality["expectedCases"], no_mortality["expectedCases"])
        self.assertLessEqual(curve_survival["expectedCases"], no_mortality["expectedCases"])

    def test_baseline_active_tb_is_excluded_from_expected_cases(self) -> None:
        result = expected_progression_events(
            [
                {"state": STATE_BASELINE_ACTIVE_TB, "weight": 900.0},
                {"state": STATE_RECENT, "weight": 100.0, "remainingEarlyRiskYears": 4.5},
            ],
            horizon_years=1.0,
            early_hazard=0.02,
            remote_hazard=0.002,
        )

        self.assertEqual(result["totalWeight"], 1000.0)
        self.assertEqual(result["includedWeight"], 100.0)
        self.assertEqual(result["excludedBaselineActiveTBWeight"], 900.0)
        self.assertAlmostEqual(result["expectedCases"], 100.0 * (1 - math.exp(-0.02)))

    def test_worked_policy_examples_are_reviewable(self) -> None:
        examples = worked_progression_calibration_examples()
        by_id = {example["exampleId"]: example["result"] for example in examples}

        self.assertIn("eligible_external_hazards", by_id)
        self.assertIn("eligible_fixed_ratio_fit", by_id)
        self.assertIn("no_multipliers_vs_moderate_reviewed_multipliers", by_id)
        self.assertIn("legacy_or_as_hazard_diagnostic", by_id)
        self.assertIn("retrospective_validation_only", by_id)
        self.assertIn("baseline_prevalence_validation_only", by_id)
        self.assertIn("impossible_high_target", by_id)
        self.assertIn("zero_case_target", by_id)
        self.assertIn("incomplete_ascertainment", by_id)
        self.assertIn("survival_curve", by_id)
        self.assertIsNotNone(by_id["impossible_high_target"])


class ActiveTBObservationSchemaTests(unittest.TestCase):
    def valid_row(self, **overrides) -> dict:
        row = {
            "observationId": "obs-1",
            "startDate": "2025-01-01",
            "endDate": "2025-12-31",
            "observedActiveTBCaseCount": 7,
            "populationDenominator": 100000,
            "personYears": 99500,
            "denominatorType": "census_population",
            "populationScope": "whole_population",
            "caseClassification": "prevalent_baseline",
            "observationWindowMeaning": "baseline_prevalent",
            "ascertainmentMethod": "passive_notification",
            "activeTBClassification": "all_active_tb",
            "source": "unit test fixture",
            "reviewStatus": "unreviewed_test_fixture",
            "notes": "schema fixture",
            "uncertainty": {"lowCount": 5, "highCount": 9},
        }
        row.update(overrides)
        return row

    def test_active_tb_observation_schema_round_trips(self) -> None:
        validated = validate_active_tb_observation(self.valid_row())
        payload = json.loads(json.dumps(validated))

        self.assertEqual(payload, validated)
        self.assertEqual(payload["schemaVersion"], ACTIVE_TB_OBSERVATION_SCHEMA_VERSION)
        self.assertAlmostEqual(payload["observedRatePer100000Population"], 7.0)

    def test_invalid_active_tb_denominators_and_dates_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "populationDenominator"):
            validate_active_tb_observation(self.valid_row(populationDenominator=0))
        with self.assertRaisesRegex(ValueError, "start"):
            validate_active_tb_observation(
                self.valid_row(startDate="2026-01-01", endDate="2025-01-01")
            )

    def test_baseline_prevalent_and_incident_observations_remain_distinguishable(self) -> None:
        baseline, incident = validate_active_tb_observations(
            [
                self.valid_row(
                    observationId="baseline",
                    caseClassification="prevalent_baseline",
                    observationWindowMeaning="baseline_prevalent",
                ),
                self.valid_row(
                    observationId="incident",
                    startYear=2027,
                    endYear=2027,
                    startDate=None,
                    endDate=None,
                    caseClassification="incident_follow_up",
                    observationWindowMeaning="follow_up_incident",
                    ascertainmentMethod="combined",
                ),
            ]
        )

        self.assertEqual(baseline["caseClassification"], "prevalent_baseline")
        self.assertEqual(incident["caseClassification"], "incident_follow_up")
        self.assertNotEqual(
            baseline["observationWindowMeaning"],
            incident["observationWindowMeaning"],
        )


class ScopeProtectionTests(unittest.TestCase):
    def test_existing_frozen_sa_health_outputs_and_artifacts_remain_unchanged(self) -> None:
        protected_paths = [
            "app/reference_data",
            "reports/sa_health_apy_health_economic_working_report_v1",
            "validation/matlab_reference",
            "abm",
        ]
        changed = subprocess.check_output(
            ["git", "diff", "--name-only", "--", *protected_paths],
            cwd=REPO_ROOT,
            text=True,
        ).splitlines()

        self.assertEqual(changed, [])

    def test_no_frozen_release_branch_or_tag_pointer_changes(self) -> None:
        branch = subprocess.check_output(
            ["git", "rev-parse", "release/sa-health-apy-he-v1.0.0"],
            cwd=REPO_ROOT,
            text=True,
        ).strip()
        tag = subprocess.check_output(
            ["git", "rev-parse", "sa-health-apy-he-v1.0.0^{}"],
            cwd=REPO_ROOT,
            text=True,
        ).strip()

        self.assertEqual(branch, FROZEN_COMMIT)
        self.assertEqual(tag, FROZEN_COMMIT)

    def test_no_streamlit_page_exposes_new_pathway_during_milestone(self) -> None:
        combined = "\n".join(
            path.read_text(encoding="utf-8")
            for path in (REPO_ROOT / "pages").glob("*.py")
        )

        self.assertNotIn("explicit_recent_remote_tbi", combined)
        self.assertNotIn("explicit_recent_remote_progression", combined)
        self.assertNotIn("Recently infected within 5 years", combined)
        self.assertNotIn("Remote infection only", combined)

    def test_no_runner_ledger_economics_daly_or_frozen_loader_consumes_new_assignments(self) -> None:
        protected_files = [
            REPO_ROOT / "engine" / "apy" / "runner.py",
            REPO_ROOT / "engine" / "apy" / "simulation.py",
            REPO_ROOT / "engine" / "apy" / "expected_value.py",
            REPO_ROOT / "engine" / "apy" / "event_ledger.py",
            REPO_ROOT / "engine" / "apy" / "event_ledger_economics.py",
            REPO_ROOT / "engine" / "apy" / "economics.py",
            REPO_ROOT / "engine" / "apy" / "frozen_reference.py",
        ]
        protected_files.extend((REPO_ROOT / "engine" / "dynamic").rglob("*.py"))
        combined = "\n".join(
            path.read_text(encoding="utf-8") for path in protected_files if path.exists()
        )

        self.assertNotIn("explicit_recent_remote_tbi", combined)
        self.assertNotIn("explicit_recent_remote_progression", combined)
        self.assertNotIn("deterministic_recent_remote_assignment", combined)
        self.assertNotIn("stochastic_recent_remote_population_assignment", combined)
        self.assertNotIn("ExplicitRecentRemoteConfig", combined)

    def test_retired_pathway_is_not_reactivated(self) -> None:
        combined = "\n".join(
            path.read_text(encoding="utf-8")
            for path in (REPO_ROOT / "pages").glob("*.py")
        )
        new_module = (REPO_ROOT / "engine" / "apy" / "explicit_recent_remote_tbi.py").read_text(
            encoding="utf-8"
        )

        self.assertNotIn("Enable experimental infection-history analysis", combined)
        self.assertNotIn("Historical TB infection pressure", combined)
        self.assertNotIn("continuous_markov_recent_remote", new_module)
        self.assertNotIn("engine.apy.infection_history", new_module)
        self.assertNotIn("engine.apy.ltbi_state", new_module)


if __name__ == "__main__":
    unittest.main()
