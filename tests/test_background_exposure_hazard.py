"""Background-exposure hazard mathematics: analytical, numerical, statistical and isolation tests."""

from __future__ import annotations

import math
from pathlib import Path
import random
import subprocess
import sys
import unittest
from unittest import mock

import numpy as np

from engine.profiles.background_exposure import (
    BACKGROUND_EXPOSURE_SCHEMA_VERSION,
    AgeBandHazard,
    BackgroundExposure,
    ExposureMode,
    HazardYear,
    hazard_value,
    no_background_exposure,
)
from engine.profiles.background_exposure_hazard import (
    AgeSpecificHazardNotExecutableError,
    HazardCoverageError,
    HazardInputError,
    HazardSchedule,
    InfectionOccurs,
    NoInfectionInInterval,
    cumulative_hazard,
    exponential_threshold,
    first_infection_time,
    first_infection_times,
    hazard_at,
    infection_probabilities,
    infection_probability,
    plausibility_warnings,
    probability_from_cumulative_hazard,
    requires_random_draw,
    schedule_from_config,
)


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "engine" / "profiles" / "background_exposure_hazard.py"
NONE = HazardSchedule.none()
# Hand-checked series used throughout: year y applies on [y, y + 1).
SERIES = HazardSchedule.annual_series({2020: 0.02, 2021: 0.05, 2022: 0.10, 2023: 0.01})
# Statistical tolerance, fixed before running: |empirical - analytical| <= 4 standard errors.
SE_MULTIPLIER = 4.0
SAMPLE_SIZE = 200_000
SEED = 20261001


class ModeNoneTests(unittest.TestCase):
    def test_zero_hazard_probability_and_no_event(self) -> None:
        for t0, t1 in ((0.0, 0.0), (0.0, 1.0), (2020.3, 2045.9), (-5.0, 1e6)):
            with self.subTest(interval=(t0, t1)):
                self.assertEqual(hazard_at(NONE, t0), 0.0)
                self.assertEqual(cumulative_hazard(NONE, t0, t1), 0.0)
                self.assertEqual(infection_probability(NONE, t0, t1), 0.0)
                result = first_infection_time(NONE, t0, t1)
                self.assertIsInstance(result, NoInfectionInInterval)
                self.assertFalse(result.infected)
                self.assertIsNone(result.residual_threshold)
        infected, times = first_infection_times(NONE, 0.0, 10.0, np.empty(0))
        self.assertEqual(infected.size, 0)
        self.assertEqual(times.size, 0)

    def test_none_requires_and_accepts_no_random_variate(self) -> None:
        self.assertFalse(requires_random_draw(NONE))
        self.assertTrue(requires_random_draw(HazardSchedule.constant(0.0)))
        with self.assertRaises(HazardInputError):
            first_infection_time(NONE, 0.0, 1.0, uniform=0.5)
        with self.assertRaises(HazardInputError):
            first_infection_time(NONE, 0.0, 1.0, threshold=0.1)
        with self.assertRaises(HazardInputError):
            first_infection_times(NONE, 0.0, 1.0, np.array([0.1]))

    def test_none_consumes_no_random_draws(self) -> None:
        """A protocol-following caller draws nothing for mode none, and the stream state is untouched."""

        class CountingGenerator:
            def __init__(self, seed: int) -> None:
                self.generator = np.random.default_rng(seed)
                self.calls = 0
                self.draws = 0

            def random(self, size: int) -> np.ndarray:
                self.calls += 1
                self.draws += size
                return self.generator.random(size)

        def reference_caller(schedule: HazardSchedule, rng: CountingGenerator, people: int):
            if not requires_random_draw(schedule):
                return first_infection_times(schedule, 2025.0, 2035.0, np.empty(0))
            return first_infection_times(schedule, 2025.0, 2035.0, -np.log1p(-rng.random(people)))

        rng = CountingGenerator(SEED)
        before = rng.generator.bit_generator.state
        reference_caller(NONE, rng, 1000)
        self.assertEqual((rng.calls, rng.draws), (0, 0))
        self.assertEqual(rng.generator.bit_generator.state, before)
        reference_caller(HazardSchedule.constant(0.01), rng, 1000)
        self.assertEqual((rng.calls, rng.draws), (1, 1000))

    def test_module_never_draws_random_numbers(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertNotRegex(source, r"np\.random|numpy\.random|^\s*import random|^\s*from random|default_rng|secrets")

        def forbidden(*args, **kwargs):
            raise AssertionError("random number generator called")

        with mock.patch("numpy.random.default_rng", forbidden), mock.patch("numpy.random.random", forbidden), mock.patch.object(
            random, "random", forbidden
        ):
            numpy_state = np.random.get_state()[1].copy()
            python_state = random.getstate()
            for schedule in (NONE, HazardSchedule.constant(0.1), SERIES):
                cumulative_hazard(schedule, 2020.5, 2023.5)
                infection_probability(schedule, 2020.5, 2023.5)
            first_infection_time(NONE, 2020.5, 2023.5)
            first_infection_time(SERIES, 2020.5, 2023.5, uniform=0.3)
            first_infection_times(SERIES, 2020.5, 2023.5, np.array([0.01, 0.2]))
        np.testing.assert_array_equal(np.random.get_state()[1], numpy_state)
        self.assertEqual(random.getstate(), python_state)

    def test_existing_none_configuration_unchanged(self) -> None:
        config = no_background_exposure()
        self.assertEqual(
            config.to_dict(),
            {
                "schemaVersion": BACKGROUND_EXPOSURE_SCHEMA_VERSION,
                "mode": "none",
                "quantity": "exogenous_infection_hazard",
                "unit": "infections per person-year",
                "constantHazard": None,
                "ageSpecificHazards": [],
                "timeSeries": [],
                "source": "",
                "citation": "",
                "reviewStatus": "not_required",
                "applicabilityNotes": "",
            },
        )
        # Hash computed from the schema at commit 9003146, before this milestone.
        self.assertEqual(config.exposure_hash(), "39cac6a7ff746ea8d6a85f39d5c1217ceaa33f4952b205f856b01328faf91967")
        self.assertEqual(schedule_from_config(config), NONE)
        self.assertEqual(schedule_from_config(BackgroundExposure.from_dict(None)), NONE)


class ConstantHazardTests(unittest.TestCase):
    def test_cumulative_hazard_and_probability(self) -> None:
        for hazard in (0.0, 0.001, 0.01, 0.1, 0.7, 3.0):
            for t0, t1 in ((0.0, 1.0), (0.0, 10.0), (2020.25, 2023.75), (2031.9, 2032.1)):
                with self.subTest(hazard=hazard, interval=(t0, t1)):
                    schedule = HazardSchedule.constant(hazard)
                    expected = hazard * (t1 - t0)
                    self.assertTrue(math.isclose(cumulative_hazard(schedule, t0, t1), expected, rel_tol=1e-15, abs_tol=0.0))
                    self.assertTrue(math.isclose(infection_probability(schedule, t0, t1), 1.0 - math.exp(-expected), rel_tol=1e-12, abs_tol=1e-16))
                    self.assertEqual(hazard_at(schedule, t0), hazard)

    def test_exact_inversion(self) -> None:
        for hazard in (0.002, 0.05, 1.7):
            schedule = HazardSchedule.constant(hazard)
            for uniform in (0.999999, 0.95, 0.5, 0.2):
                with self.subTest(hazard=hazard, uniform=uniform):
                    threshold = -math.log(uniform)
                    result = first_infection_time(schedule, 2025.0, 2025.0 + 1e4, uniform=uniform)
                    self.assertIsInstance(result, InfectionOccurs)
                    self.assertTrue(math.isclose(result.time, 2025.0 + threshold / hazard, rel_tol=1e-15))
                    # Round trip: limited only by float spacing of calendar time near 2025 (about 2e-13).
                    spacing = 4 * math.ulp(2025.0 + 1e4) * hazard
                    self.assertTrue(math.isclose(cumulative_hazard(schedule, 2025.0, result.time), threshold, rel_tol=1e-12, abs_tol=spacing))

    def test_no_event_when_threshold_exceeds_interval_hazard(self) -> None:
        schedule = HazardSchedule.constant(0.01)
        result = first_infection_time(schedule, 0.0, 5.0, threshold=0.06)
        self.assertIsInstance(result, NoInfectionInInterval)
        self.assertAlmostEqual(result.cumulative_hazard, 0.05, places=15)
        self.assertAlmostEqual(result.residual_threshold, 0.01, places=15)
        # E == H is not an infection within [t0, t1): it would occur at t1, in the next interval.
        self.assertFalse(first_infection_time(schedule, 0.0, 5.0, threshold=0.05).infected)
        # Zero hazard never infects, even at E = 0.
        self.assertFalse(first_infection_time(HazardSchedule.constant(0.0), 0.0, 5.0, threshold=0.0).infected)
        # E = 0 with a positive hazard infects at t0.
        self.assertEqual(first_infection_time(schedule, 3.0, 5.0, threshold=0.0).time, 3.0)


class PiecewiseConstantTests(unittest.TestCase):
    def test_hazard_lookup_and_year_boundaries(self) -> None:
        self.assertEqual(hazard_at(SERIES, 2020.0), 0.02)
        self.assertEqual(hazard_at(SERIES, math.nextafter(2021.0, 0.0)), 0.02)
        self.assertEqual(hazard_at(SERIES, 2021.0), 0.05)
        self.assertEqual(hazard_at(SERIES, 2023.999), 0.01)
        for outside in (2019.999, 2024.0):
            with self.subTest(t=outside), self.assertRaises(HazardCoverageError):
                hazard_at(SERIES, outside)

    def test_hand_calculated_cumulative_hazards(self) -> None:
        cases = {
            "within one year": ((2021.25, 2021.75), 0.05 * 0.5),
            "spanning two years": ((2020.5, 2021.25), 0.02 * 0.5 + 0.05 * 0.25),
            "spanning several years": ((2020.75, 2023.5), 0.02 * 0.25 + 0.05 + 0.10 + 0.01 * 0.5),
            "exact year boundaries": ((2021.0, 2023.0), 0.05 + 0.10),
            "whole series": ((2020.0, 2024.0), 0.18),
            "fractional within first year": ((2020.1, 2020.3), 0.02 * 0.2),
        }
        for name, ((t0, t1), expected) in cases.items():
            with self.subTest(case=name):
                self.assertTrue(math.isclose(cumulative_hazard(SERIES, t0, t1), expected, rel_tol=1e-11))
                self.assertTrue(math.isclose(infection_probability(SERIES, t0, t1), -math.expm1(-expected), rel_tol=1e-11))

    def test_exact_boundaries_need_only_touched_years(self) -> None:
        short = HazardSchedule.annual_series({2021: 0.05, 2022: 0.10})
        self.assertAlmostEqual(cumulative_hazard(short, 2021.0, 2023.0), 0.15, places=15)
        with self.assertRaises(HazardCoverageError):
            cumulative_hazard(short, 2021.0, math.nextafter(2023.0, 3000.0))
        with self.assertRaises(HazardCoverageError):
            cumulative_hazard(short, math.nextafter(2021.0, 0.0), 2022.0)

    def test_additivity_and_inversion_across_years(self) -> None:
        whole = cumulative_hazard(SERIES, 2020.5, 2023.5)
        parts = sum(cumulative_hazard(SERIES, a, b) for a, b in ((2020.5, 2021.3), (2021.3, 2022.0), (2022.0, 2023.5)))
        self.assertAlmostEqual(whole, parts, places=14)
        # E = 0.03 from 2020.5: 0.01 accrues by 2021, then 0.02 / 0.05 = 0.4 years.
        self.assertAlmostEqual(first_infection_time(SERIES, 2020.5, 2023.5, threshold=0.03).time, 2021.4, places=12)
        # E = 0.1: 0.06 by 2022, then 0.04 / 0.10 = 0.4 years.
        single = first_infection_time(SERIES, 2020.5, 2023.5, threshold=0.1)
        self.assertAlmostEqual(single.time, 2022.4, places=12)
        first = first_infection_time(SERIES, 2020.5, 2022.0, threshold=0.1)
        self.assertIsInstance(first, NoInfectionInInterval)
        second = first_infection_time(SERIES, 2022.0, 2023.5, threshold=first.residual_threshold)
        self.assertAlmostEqual(second.time, single.time, places=12)

    def test_zero_hazard_year_is_skipped_by_inversion(self) -> None:
        gap = HazardSchedule.annual_series({2020: 0.1, 2021: 0.0, 2022: 0.3})
        self.assertAlmostEqual(cumulative_hazard(gap, 2020.5, 2022.5), 0.2, places=15)
        # tau = inf{s : H(t0, s) > E}: the end of the zero-hazard plateau.
        self.assertEqual(first_infection_time(gap, 2020.5, 2022.5, threshold=0.05).time, 2022.0)
        self.assertAlmostEqual(first_infection_time(gap, 2020.5, 2022.5, threshold=0.06).time, 2022.0 + 0.01 / 0.3, places=12)

    def test_scalar_and_vector_forms_agree(self) -> None:
        thresholds = np.concatenate(([0.0, 0.01, 0.06, 0.1599, 0.16, 0.5], -np.log(np.linspace(0.01, 0.99, 197))))
        infected, times = first_infection_times(SERIES, 2020.75, 2023.5, thresholds)
        for threshold, flag, time in zip(thresholds, infected, times):
            result = first_infection_time(SERIES, 2020.75, 2023.5, threshold=float(threshold))
            self.assertEqual(result.infected, bool(flag))
            if flag:
                self.assertTrue(math.isclose(result.time, time, rel_tol=0.0, abs_tol=1e-9))
            else:
                self.assertTrue(np.isnan(time))
        probs = infection_probabilities(SERIES, [(2020.0, 2021.0), (2020.0, 2024.0)])
        np.testing.assert_allclose(probs, [-math.expm1(-0.02), -math.expm1(-0.18)], rtol=1e-12)


class NumericalBehaviourTests(unittest.TestCase):
    def test_very_small_hazards_avoid_cancellation(self) -> None:
        for cumulative in (1e-18, 1e-12, 1e-8):
            with self.subTest(H=cumulative):
                exact = cumulative - cumulative**2 / 2 + cumulative**3 / 6
                self.assertTrue(math.isclose(probability_from_cumulative_hazard(cumulative), exact, rel_tol=1e-15))
        # The naive formula loses most significant digits at H = 1e-12; expm1 does not.
        naive = 1.0 - math.exp(-1e-12)
        self.assertGreater(abs(naive - 1e-12) / 1e-12, 1e-6)
        schedule = HazardSchedule.constant(1e-12)
        self.assertTrue(math.isclose(infection_probability(schedule, 0.0, 1.0), 1e-12, rel_tol=1e-11))

    def test_large_finite_hazards_and_long_intervals(self) -> None:
        large = HazardSchedule.constant(1e3)
        self.assertEqual(infection_probability(large, 0.0, 50.0), 1.0)
        self.assertAlmostEqual(cumulative_hazard(large, 0.0, 50.0), 5e4, places=8)
        self.assertAlmostEqual(first_infection_time(large, 0.0, 50.0, threshold=2.0).time, 0.002, places=15)
        long_series = HazardSchedule.annual_series({year: 0.001 * (1 + year % 7) for year in range(1950, 2150)})
        expected = math.fsum(0.001 * (1 + year % 7) for year in range(1950, 2150))
        self.assertTrue(math.isclose(cumulative_hazard(long_series, 1950.0, 2150.0), expected, rel_tol=1e-13))
        with self.assertRaises(HazardInputError):
            cumulative_hazard(HazardSchedule.constant(1e308), 0.0, 1e10)

    def test_monotonic_and_bounded(self) -> None:
        hazards = [0.0, 1e-6, 1e-3, 0.01, 0.1, 1.0, 10.0, 100.0]
        stops = [0.0, 0.01, 0.5, 1.0, 5.0, 50.0]
        for t1 in stops:
            values = [infection_probability(HazardSchedule.constant(h), 0.0, t1) for h in hazards]
            self.assertEqual(values, sorted(values))
            self.assertTrue(all(0.0 <= v <= 1.0 for v in values))
        for h in hazards:
            values = [cumulative_hazard(HazardSchedule.constant(h), 0.0, t1) for t1 in stops]
            self.assertEqual(values, sorted(values))
        series_values = [infection_probability(SERIES, 2020.0, 2020.0 + k * 0.25) for k in range(17)]
        self.assertEqual(series_values, sorted(series_values))
        self.assertTrue(all(0.0 <= v <= 1.0 for v in series_values))


class ValidationTests(unittest.TestCase):
    def test_invalid_hazards(self) -> None:
        for bad in (-0.001, float("nan"), float("inf"), float("-inf"), True, "0.1", None):
            with self.subTest(hazard=bad):
                with self.assertRaises(HazardInputError):
                    HazardSchedule.constant(bad)
                with self.assertRaises(HazardInputError):
                    HazardSchedule.annual_series({2020: 0.01, 2021: bad})

    def test_series_structure(self) -> None:
        with self.assertRaises(HazardInputError):
            HazardSchedule.annual_series([(2020, 0.1), (2020, 0.2)])
        with self.assertRaises(HazardInputError):
            HazardSchedule.annual_series({})
        with self.assertRaises(HazardInputError):
            HazardSchedule.annual_series([(2020.5, 0.1)])
        self.assertEqual(HazardSchedule.annual_series([(2021, 0.2), (2020, 0.1)]).years, (2020, 2021))

    def test_missing_years_and_coverage_are_errors_not_fills(self) -> None:
        gappy = HazardSchedule.annual_series({2020: 0.02, 2022: 0.10})
        with self.assertRaisesRegex(HazardCoverageError, r"\[2021\]"):
            cumulative_hazard(gappy, 2020.5, 2022.5)
        for interval in ((2019.5, 2020.5), (2022.5, 2023.5), (2030.0, 2031.0)):
            with self.subTest(interval=interval), self.assertRaises(HazardCoverageError):
                infection_probability(gappy, *interval)
        with self.assertRaises(HazardCoverageError):
            first_infection_time(gappy, 2020.5, 2022.5, threshold=1.0)
        self.assertAlmostEqual(cumulative_hazard(gappy, 2020.5, 2021.0), 0.01, places=15)

    def test_missing_configuration_values_stay_missing(self) -> None:
        config = BackgroundExposure(
            mode=ExposureMode.TIME_SERIES,
            time_series=tuple(
                HazardYear(year, hazard_value(value))
                for year, value in ((2020, 0.02), (2021, None), (2022, 0.1))
            ),
        )
        schedule = schedule_from_config(config)
        self.assertEqual(schedule.years, (2020, 2022))
        self.assertAlmostEqual(cumulative_hazard(schedule, 2022.0, 2022.5), 0.05, places=15)
        with self.assertRaises(HazardCoverageError):
            cumulative_hazard(schedule, 2020.0, 2022.5)

    def test_intervals(self) -> None:
        schedule = HazardSchedule.constant(0.1)
        with self.assertRaises(HazardInputError):
            cumulative_hazard(schedule, 2.0, 1.0)
        for bad in (float("nan"), float("inf"), float("-inf"), "2020", True):
            with self.subTest(time=bad), self.assertRaises(HazardInputError):
                cumulative_hazard(schedule, bad, 2030.0)
        # Zero duration accumulates nothing and needs no coverage.
        self.assertEqual(cumulative_hazard(schedule, 3.0, 3.0), 0.0)
        self.assertEqual(infection_probability(SERIES, 1900.0, 1900.0), 0.0)
        self.assertFalse(first_infection_time(SERIES, 1900.0, 1900.0, threshold=0.0).infected)

    def test_random_variates(self) -> None:
        for bad in (0.0, 1.0, -0.1, 1.1, float("nan"), float("inf"), True, "0.5"):
            with self.subTest(uniform=bad):
                with self.assertRaises(HazardInputError):
                    exponential_threshold(bad)
                with self.assertRaises(HazardInputError):
                    first_infection_time(SERIES, 2020.0, 2021.0, uniform=bad)
        for bad in (-1e-9, float("nan"), float("inf")):
            with self.subTest(threshold=bad):
                with self.assertRaises(HazardInputError):
                    first_infection_time(SERIES, 2020.0, 2021.0, threshold=bad)
                with self.assertRaises(HazardInputError):
                    first_infection_times(SERIES, 2020.0, 2021.0, np.array([0.1, bad]))
        with self.assertRaises(HazardInputError):
            first_infection_time(SERIES, 2020.0, 2021.0)
        with self.assertRaises(HazardInputError):
            first_infection_time(SERIES, 2020.0, 2021.0, uniform=0.5, threshold=0.1)
        self.assertAlmostEqual(exponential_threshold(math.exp(-0.25)), 0.25, places=15)

    def test_high_hazards_are_valid_with_warning(self) -> None:
        self.assertEqual(plausibility_warnings(HazardSchedule.constant(0.5)), ())
        high = HazardSchedule.constant(2.5)
        self.assertEqual(len(plausibility_warnings(high)), 1)
        self.assertAlmostEqual(infection_probability(high, 0.0, 1.0), -math.expm1(-2.5), places=15)
        self.assertEqual(plausibility_warnings(NONE), ())

    def test_configuration_bridge(self) -> None:
        constant = BackgroundExposure(mode=ExposureMode.CONSTANT, constant_hazard=hazard_value(0.004))
        self.assertEqual(schedule_from_config(constant), HazardSchedule.constant(0.004))
        age_specific = BackgroundExposure(
            mode=ExposureMode.CONSTANT,
            age_specific_hazards=(AgeBandHazard(0, 15, hazard_value(0.001)), AgeBandHazard(15, None, hazard_value(0.004))),
        )
        with self.assertRaisesRegex(AgeSpecificHazardNotExecutableError, "not yet executable"):
            schedule_from_config(age_specific)


class StatisticalValidationTests(unittest.TestCase):
    """Fixed-seed Monte Carlo with caller-supplied variates; tolerance fixed in advance (4 SE)."""

    SCENARIOS = {
        "constant 0.01 over 10 years": (HazardSchedule.constant(0.01), 0.0, 10.0),
        "constant 0.1 over fractional years": (HazardSchedule.constant(0.1), 2020.25, 2023.75),
        "constant 2.0 over half a year": (HazardSchedule.constant(2.0), 0.0, 0.5),
        "series with a zero year": (HazardSchedule.annual_series({2020: 0.02, 2021: 0.0, 2022: 0.15, 2023: 0.05}), 2020.4, 2023.6),
        "rising series over 4 years": (SERIES, 2020.0, 2024.0),
    }

    def test_empirical_proportions_and_time_distribution(self) -> None:
        rng = np.random.default_rng(SEED)
        for name, (schedule, t0, t1) in self.SCENARIOS.items():
            uniforms = rng.random(SAMPLE_SIZE)
            self.assertTrue(np.all(uniforms > 0.0))
            infected, times = first_infection_times(schedule, t0, t1, -np.log(uniforms))
            expected = infection_probability(schedule, t0, t1)
            se = math.sqrt(expected * (1.0 - expected) / SAMPLE_SIZE)
            with self.subTest(scenario=name, quantity="proportion infected"):
                self.assertLessEqual(abs(infected.mean() - expected), SE_MULTIPLIER * se)
            for fraction in (0.2, 0.5, 0.8):
                s = t0 + fraction * (t1 - t0)
                cdf = infection_probability(schedule, t0, s)
                empirical = np.mean(infected & (times < s))
                se_s = math.sqrt(cdf * (1.0 - cdf) / SAMPLE_SIZE)
                with self.subTest(scenario=name, quantity=f"P(T < {s:g})"):
                    self.assertLessEqual(abs(empirical - cdf), SE_MULTIPLIER * se_s)
            with self.subTest(scenario=name, quantity="times inside interval"):
                self.assertTrue(np.all((times[infected] >= t0) & (times[infected] < t1)))


class IsolationTests(unittest.TestCase):
    def test_not_imported_by_execution_paths(self) -> None:
        code = (
            "import sys\n"
            "import engine.apy.runner, engine.apy.simulation, engine.apy.expected_value, engine.apy.event_ledger\n"
            "import engine.apy.event_ledger_economics, engine.apy.health_economics, engine.apy.frozen_reference\n"
            "import engine.profiles.engine_mapping, engine.profiles.country, engine.profiles.background_exposure\n"
            "import app.general.state, app.general.economics, app.general.provenance_export, app.general.terminology\n"
            "print('engine.profiles.background_exposure_hazard' in sys.modules)\n"
        )
        result = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, timeout=300)
        self.assertEqual(result.returncode, 0, result.stderr[-2000:])
        self.assertEqual(result.stdout.strip().splitlines()[-1], "False")

    def test_only_tests_reference_the_module(self) -> None:
        scanned = [
            *sorted((ROOT / "engine").rglob("*.py")),
            *sorted((ROOT / "app").rglob("*.py")),
            *sorted((ROOT / "ui").rglob("*.py")),
            *sorted((ROOT / "general_pages").glob("*.py")),
            *sorted((ROOT / "pages").glob("*.py")),
            *sorted((ROOT / "adapters").rglob("*.py")),
            *sorted(ROOT.glob("*.py")),
        ]
        for path in scanned:
            if path == MODULE:
                continue
            text = path.read_text(encoding="utf-8")
            with self.subTest(file=path.relative_to(ROOT).as_posix()):
                self.assertNotRegex(text, r"(?m)^\s*(from|import)\s.*background_exposure_hazard")
                self.assertNotRegex(text, r"[\"']engine\.profiles\.background_exposure_hazard[\"']")

    def test_no_starsim_or_streamlit(self) -> None:
        source = MODULE.read_text(encoding="utf-8")
        self.assertNotRegex(source, r"^\s*(import|from)\s+(starsim|streamlit)")


if __name__ == "__main__":
    unittest.main()
