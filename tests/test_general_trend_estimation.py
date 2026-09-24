from __future__ import annotations

import json
import math
import unittest

import numpy as np

from engine.profiles.population_profile import IncidencePoint
from engine.who_incidence import trend as trend_module
from engine.who_incidence.trend import (
    CLASSIFICATION_LABELS,
    SUMMARY_LABELS,
    CovidHandling,
    TrendMethod,
    TrendSettings,
    estimate_trend,
)


def series(values, start=2010, spread=(0.8, 1.25)):
    low, high = spread
    return [
        IncidencePoint(year=start + i, estimate=v, lower=None if v is None else v * low, upper=None if v is None else v * high)
        for i, v in enumerate(values)
        if v is not None
    ]


def geometric(rate, n=15, start=2010, base=100.0):
    return series([base * (1 + rate) ** i for i in range(n)], start=start)


class LogLinearTrendTests(unittest.TestCase):
    def test_recovers_rising_and_falling_change(self) -> None:
        rising = estimate_trend(geometric(0.05), TrendSettings(window_years=None))
        falling = estimate_trend(geometric(-0.03), TrendSettings(window_years=None))
        self.assertAlmostEqual(rising.annual_percent_change, 5.0, places=8)
        self.assertAlmostEqual(falling.annual_percent_change, -3.0, places=8)
        self.assertEqual(rising.summary, "increasing")
        self.assertEqual(falling.summary, "decreasing")
        self.assertEqual(rising.classification, "adequate_log_linear")

    def test_constant_series(self) -> None:
        narrow = series([100.0] * 15, spread=(0.97, 1.03))
        result = estimate_trend(narrow, TrendSettings(window_years=None))
        self.assertAlmostEqual(result.annual_percent_change, 0.0, places=10)
        self.assertEqual(result.summary, "no_clear_change")
        wide = estimate_trend(geometric(0.0), TrendSettings(window_years=None))
        self.assertAlmostEqual(wide.annual_percent_change, 0.0, places=10)
        self.assertEqual(wide.summary, "uncertain")
        self.assertIn("not evidence that the epidemic is constant", trend_module.NO_CLEAR_CHANGE_NOTE)

    def test_default_window_is_ten_years_ending_at_latest(self) -> None:
        result = estimate_trend(geometric(-0.02, n=25, start=2000))
        self.assertEqual(result.period, (2015, 2024))
        self.assertEqual(result.years_used, 10)
        five = estimate_trend(geometric(-0.02, n=25, start=2000), TrendSettings(window_years=5))
        self.assertEqual(five.period, (2020, 2024))
        custom = estimate_trend(geometric(-0.02, n=25, start=2000), TrendSettings(start_year=2005, end_year=2014))
        self.assertEqual(custom.period, (2005, 2014))

    def test_asymmetric_bounds_stay_asymmetric(self) -> None:
        points = [IncidencePoint(2010 + i, 100.0 * 0.97**i, 100.0 * 0.97**i * 0.9, 100.0 * 0.97**i * 1.6) for i in range(10)]
        result = estimate_trend(points, TrendSettings(window_years=None))
        ratios = result.diagnostics["propagation"]["asymmetry"]
        self.assertTrue(all(ratio > 3 for ratio in ratios))
        symmetric = [IncidencePoint(p.year, p.estimate, p.estimate / 1.25, p.estimate * 1.25) for p in points]
        other = estimate_trend(symmetric, TrendSettings(window_years=None))
        self.assertTrue(all(abs(r - 1.0) < 1e-9 for r in other.diagnostics["propagation"]["asymmetry"]))
        self.assertNotEqual(result.propagated_interval, other.propagated_interval)
        # split-normal sampler: upper tail wider than lower tail on the log scale
        rng = np.random.default_rng(1)
        z = rng.standard_normal(200_000)
        m, s_low, s_high = 0.0, 0.1 / 1.96, 0.5 / 1.96
        draws = m + np.where(z > 0, s_high, s_low) * z
        self.assertGreater(np.percentile(draws, 97.5), -np.percentile(draws, 2.5) * 3)

    def test_missing_years_reported_not_imputed(self) -> None:
        values = [100 * 0.97**i for i in range(12)]
        values[4] = None
        values[5] = None
        result = estimate_trend(series(values), TrendSettings(window_years=None))
        self.assertEqual(result.diagnostics["missingYears"], [2014, 2015])
        self.assertEqual(result.years_used, 10)
        self.assertTrue(any("not imputed" in warning for warning in result.warnings))
        self.assertNotIn(2014, [p.year for p in result.points])

    def test_short_series_flagged_and_very_short_rejected(self) -> None:
        short = estimate_trend(geometric(-0.02, n=4), TrendSettings(window_years=None))
        self.assertEqual(short.classification, "short_series")
        tiny = estimate_trend(geometric(-0.02, n=2), TrendSettings(window_years=None))
        self.assertEqual(tiny.classification, "unsuitable")
        self.assertIsNone(tiny.annual_percent_change)
        self.assertEqual(tiny.summary, "not_estimated")

    def test_zero_and_near_zero_incidence(self) -> None:
        zero = estimate_trend(series([1.0, 0.5, 0.0, 0.2, 0.1, 0.1]), TrendSettings(window_years=None))
        self.assertEqual(zero.classification, "unsuitable")
        near = estimate_trend(series([0.09, 0.08, 0.07, 0.06, 0.05, 0.05]), TrendSettings(window_years=None))
        self.assertTrue(any("near zero" in warning for warning in near.warnings))

    def test_implausible_bounds_rejected(self) -> None:
        points = [IncidencePoint(2010 + i, 100.0, 120.0, 130.0) for i in range(8)]
        result = estimate_trend(points, TrendSettings(window_years=None))
        self.assertEqual(result.classification, "unsuitable")

    def test_missing_bounds_disable_propagation(self) -> None:
        points = [IncidencePoint(2010 + i, 100 * 0.98**i, None, None) for i in range(10)]
        result = estimate_trend(points, TrendSettings(window_years=None))
        self.assertIsNone(result.propagated_interval)
        self.assertIn("bounds_incomplete", result.flags)
        self.assertIsNotNone(result.fit_interval)

    def test_poor_log_linear_fit_warns(self) -> None:
        values = [100 * math.exp(-0.05 * (i - 7) ** 2) + 10 for i in range(15)]
        result = estimate_trend(series(values), TrendSettings(window_years=None))
        self.assertIn("non_linear", result.flags)
        self.assertTrue(any("log-linear description is poor" in warning or "Residual scatter" in warning for warning in result.warnings))


class CovidHandlingTests(unittest.TestCase):
    def setUp(self) -> None:
        values = [100 * 0.97**i for i in range(15)]
        for year in (2020, 2021, 2022):
            values[year - 2010] *= 0.7
        self.points = series(values)

    def test_default_includes_all_years(self) -> None:
        self.assertEqual(TrendSettings().covid_handling, CovidHandling.INCLUDE)
        result = estimate_trend(self.points, TrendSettings(window_years=None))
        self.assertTrue(all(p.used_in_fit for p in result.points))
        self.assertIn("disruption_sensitive", result.flags)

    def test_exclusion_changes_only_intended_years(self) -> None:
        included = estimate_trend(self.points, TrendSettings(window_years=None))
        excluded = estimate_trend(self.points, TrendSettings(window_years=None, covid_handling=CovidHandling.EXCLUDE))
        unused = sorted(p.year for p in excluded.points if not p.used_in_fit)
        self.assertEqual(unused, [2020, 2021, 2022])
        self.assertEqual(excluded.years_used, included.years_used - 3)
        self.assertAlmostEqual(excluded.annual_percent_change, -3.0, places=6)
        custom = estimate_trend(self.points, TrendSettings(window_years=None, excluded_years=(2012,)))
        self.assertEqual(sorted(p.year for p in custom.points if not p.used_in_fit), [2012])

    def test_segmented_disruption_is_reproducible_and_recovers_trend(self) -> None:
        settings = TrendSettings(window_years=None, covid_handling=CovidHandling.SEGMENTED)
        first = estimate_trend(self.points, settings)
        second = estimate_trend(self.points, settings)
        self.assertEqual(json.dumps(first.to_dict(), sort_keys=True), json.dumps(second.to_dict(), sort_keys=True))
        self.assertAlmostEqual(first.annual_percent_change, -3.0, places=6)
        self.assertAlmostEqual(first.diagnostics["disruptionLevelShiftPercent"], -30.0, places=6)
        self.assertIn("not a causal", first.diagnostics["disruptionNote"])

    def test_segmented_not_available_for_smooth(self) -> None:
        with self.assertRaises(ValueError):
            estimate_trend(self.points, TrendSettings(method=TrendMethod.PENALISED_SPLINE, covid_handling=CovidHandling.SEGMENTED))


class SmoothTrendTests(unittest.TestCase):
    def test_spline_deterministic_and_conservative(self) -> None:
        rng = np.random.default_rng(3)
        values = [100 * 0.96**i * math.exp(rng.normal(0, 0.03)) for i in range(20)]
        settings = TrendSettings(method=TrendMethod.PENALISED_SPLINE, window_years=None)
        first = estimate_trend(series(values), settings)
        second = estimate_trend(series(values), settings)
        self.assertEqual(json.dumps(first.to_dict(), sort_keys=True), json.dumps(second.to_dict(), sort_keys=True))
        self.assertGreaterEqual(first.diagnostics["smoothingParameter"], trend_module.SPLINE_LAMBDA_FLOOR)
        self.assertLess(first.diagnostics["effectiveDegreesOfFreedom"], 10)
        self.assertAlmostEqual(first.annual_percent_change, -4.0, delta=1.5)
        self.assertIn("smooth_boundary", first.flags)
        self.assertIsNotNone(first.propagated_interval)

    def test_spline_recovers_exact_geometric_trend(self) -> None:
        result = estimate_trend(geometric(-0.03, n=12), TrendSettings(method=TrendMethod.PENALISED_SPLINE, window_years=None))
        self.assertAlmostEqual(result.annual_percent_change, -3.0, places=6)


class SummaryAndTerminologyTests(unittest.TestCase):
    def test_summaries_agree_with_intervals(self) -> None:
        for rate in (-0.08, -0.02, 0.0, 0.03, 0.1):
            result = estimate_trend(geometric(rate), TrendSettings(window_years=None))
            lower, upper = result.propagated_interval
            if lower > 0:
                self.assertEqual(result.summary, "increasing")
            elif upper < 0:
                self.assertEqual(result.summary, "decreasing")
            else:
                self.assertIn(result.summary, {"no_clear_change", "uncertain"})
            self.assertLessEqual(lower, result.annual_percent_change + 1e-9)
            self.assertGreaterEqual(upper, result.annual_percent_change - 1e-9)

    def test_propagation_is_seeded(self) -> None:
        a = estimate_trend(geometric(-0.02), TrendSettings(seed=7))
        b = estimate_trend(geometric(-0.02), TrendSettings(seed=7))
        c = estimate_trend(geometric(-0.02), TrendSettings(seed=8))
        self.assertEqual(a.propagated_interval, b.propagated_interval)
        self.assertNotEqual(a.propagated_interval, c.propagated_interval)
        self.assertEqual(a.diagnostics["propagation"]["seed"], 7)

    def test_incidence_trend_never_labelled_transmission(self) -> None:
        labels = list(CLASSIFICATION_LABELS.values()) + list(SUMMARY_LABELS.values())
        for label in labels:
            self.assertNotIn("transmission", label.lower())
            self.assertNotIn("infection", label.lower())
        payload = json.dumps(estimate_trend(geometric(-0.02)).to_dict()).lower()
        self.assertNotIn("transmission", payload)
        self.assertEqual(estimate_trend(geometric(-0.02)).quantity, "estimated_tb_disease_incidence")

    def test_settings_round_trip(self) -> None:
        settings = TrendSettings(method=TrendMethod.PENALISED_SPLINE, window_years=15, covid_handling=CovidHandling.EXCLUDE, excluded_years=(2011,))
        self.assertEqual(TrendSettings.from_dict(json.loads(json.dumps(settings.to_dict()))), settings)

    def test_real_snapshot_country(self) -> None:
        from engine.who_incidence.snapshot import load_snapshot, find_manifests, FIXTURE_DIR

        snapshot = load_snapshot(find_manifests(FIXTURE_DIR)[0])
        result = estimate_trend(snapshot.series_for("ZAF"))
        self.assertEqual(result.period, (2015, 2024))
        self.assertLess(result.annual_percent_change, 0)
        self.assertEqual(result.summary, "decreasing")
        empty = estimate_trend(snapshot.series_for("PRK"))
        self.assertEqual(empty.classification, "unsuitable")


if __name__ == "__main__":
    unittest.main()
