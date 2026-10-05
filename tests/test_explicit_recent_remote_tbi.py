from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
import unittest

import numpy as np

from engine.apy.explicit_recent_remote_tbi import (
    ACTIVE_TB_OBSERVATION_SCHEMA_VERSION,
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
        self.assertEqual(payload["recentHazardShape"], RECENT_HAZARD_SHAPE)
        self.assertEqual(payload["remoteHazardShape"], REMOTE_HAZARD_SHAPE)

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
            max(0.0, 5.0 - time_since),
            places=12,
        )
        self.assertLess(assignment["remainingEarlyRiskYears"][0], 5.0)

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
