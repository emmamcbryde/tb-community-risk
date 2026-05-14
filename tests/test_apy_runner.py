from __future__ import annotations

import unittest

import pandas as pd

from engine.apy.runner import run_replicates, run_scenario


class ApyRunnerTests(unittest.TestCase):
    def test_run_replicates_small_scenario(self) -> None:
        results = run_replicates(n=100, n_reps=5, seed=1)

        self.assertEqual(len(results["raw"]), 5)
        self.assertIsInstance(results["summary"], pd.DataFrame)
        self.assertIsNotNone(results["exampleCohort"])
        self.assertEqual(results["modelVersion"], "python_apy_v9_port")
        self.assertEqual(results["backend"], "python")

    def test_summary_dataframe_has_expected_columns(self) -> None:
        summary = run_replicates(n=100, n_reps=3, seed=2)["summary"]

        self.assertEqual(list(summary.columns), ["Metric", "Median", "Low95", "High95"])

    def test_reproducible_with_same_seed(self) -> None:
        first = run_replicates(n=100, n_reps=3, seed=3)["raw"]
        second = run_replicates(n=100, n_reps=3, seed=3)["raw"]

        pd.testing.assert_frame_equal(first, second)

    def test_different_seed_changes_stochastic_output(self) -> None:
        first = run_replicates(n=100, n_reps=3, seed=4)["raw"]
        second = run_replicates(n=100, n_reps=3, seed=5)["raw"]

        self.assertFalse(first.equals(second))

    def test_all_screening_strategies_run(self) -> None:
        for strategy in ["random", "ltbi", "cure", "prevent"]:
            with self.subTest(strategy=strategy):
                results = run_replicates(
                    {"screeningStrategy": strategy}, n=80, n_reps=2, seed=6
                )
                self.assertEqual(results["strategy"]["screeningStrategy"], strategy)

    def test_igra_and_tst_both_run(self) -> None:
        igra = run_replicates({"testType": "IGRA"}, n=80, n_reps=2, seed=7)
        tst = run_replicates({"testType": "TST"}, n=80, n_reps=2, seed=7)

        self.assertEqual(igra["strategy"]["testType"], "IGRA")
        self.assertEqual(tst["strategy"]["testType"], "TST")

    def test_nreps_one_works(self) -> None:
        results = run_replicates(n=50, n_reps=1, seed=8)

        self.assertEqual(len(results["raw"]), 1)

    def test_invalid_n_or_config_fails_clearly(self) -> None:
        with self.assertRaisesRegex(ValueError, "N must be > 0"):
            run_replicates(n=0, n_reps=1, seed=9)
        with self.assertRaises(ValueError):
            run_replicates({"screenCoverage": 2}, n=50, n_reps=1, seed=9)

    def test_run_scenario_uses_config_nreps(self) -> None:
        results = run_scenario({"N": 50, "nReps": 2, "seed": 10})

        self.assertEqual(len(results["raw"]), 2)

    def test_raw_contains_rep_and_seed(self) -> None:
        raw = run_replicates(n=50, n_reps=2, seed=11)["raw"]

        self.assertEqual(raw["rep"].tolist(), [1, 2])
        self.assertTrue((raw["seed"] >= 0).all())


if __name__ == "__main__":
    unittest.main()
