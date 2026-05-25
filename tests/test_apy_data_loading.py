from __future__ import annotations

import unittest
from pathlib import Path

from engine.apy.age_distribution import (
    broad_age_group_from_years,
    expand_age_distribution_table,
    parse_age_band,
)
from engine.apy.config import REPO_ROOT, build_default_config, resolve_repo_path
from engine.apy.data import (
    apply_disease_or_overrides,
    apply_risk_prevalence_overrides,
    load_age_distribution,
    load_parameters_from_config,
    load_tb_csv,
    resolve_age_distribution_file,
)
from engine.apy.reference_loader import load_reference_scenario_config


class ApyDataLoadingTests(unittest.TestCase):
    def test_stale_windows_absolute_csv_file_resolves_to_repo_local_data(self) -> None:
        stale_csv = r"C:\stale\checkout\tb-community-risk\abm\default_data.csv"

        resolved = resolve_repo_path(stale_csv)
        pars = load_tb_csv(stale_csv)

        self.assertEqual(resolved, REPO_ROOT / "abm" / "default_data.csv")
        self.assertEqual(pars["ageNames"], ["0-4", "5-14", ">=15"])

    def test_unresolved_stale_windows_absolute_csv_file_raises(self) -> None:
        stale_csv = r"C:\stale\checkout\tb-community-risk\abm\definitely_missing.csv"

        with self.assertRaisesRegex(FileNotFoundError, "definitely_missing.csv"):
            load_tb_csv(stale_csv)

    def test_repo_relative_csv_file_still_resolves_to_repo_local_data(self) -> None:
        resolved = resolve_repo_path("abm/default_data.csv")
        pars = load_tb_csv("abm/default_data.csv")

        self.assertEqual(resolved, REPO_ROOT / "abm" / "default_data.csv")
        self.assertEqual(pars["popFrac"], [0.09, 0.17, 0.74])

    def test_default_apy_data_loads_expected_values(self) -> None:
        pars = load_tb_csv("abm/default_data.csv")

        self.assertEqual(pars["ageNames"], ["0-4", "5-14", ">=15"])
        self.assertEqual(pars["popFrac"], [0.09, 0.17, 0.74])
        self.assertEqual(pars["baseInfByAge"], [0.015, 0.025, 0.08])
        self.assertEqual(pars["mjPrevByAge"], [0.0, 0.05, 0.2])
        self.assertEqual(pars["contactPrevByAge"], [0.13, 0.13, 0.13])
        self.assertEqual(pars["renalPrevByAge"], [0.0, 0.01, 0.15])
        self.assertAlmostEqual(pars["totalFemalePrev"], 0.56)
        self.assertAlmostEqual(pars["totalBCGPrev"], 0.68)

    def test_default_age_distribution_resolves_and_expands(self) -> None:
        pars = load_tb_csv("abm/default_data.csv")
        age_file = resolve_age_distribution_file(None, "abm/default_data.csv")

        loaded = load_age_distribution(pars, age_file, age85_plus_max=89)

        self.assertIsNotNone(age_file)
        self.assertIn(85, loaded["exactAgeValues"])
        self.assertIn(89, loaded["exactAgeValues"])
        self.assertNotIn(90, loaded["exactAgeValues"])
        self.assertAlmostEqual(sum(loaded["exactAgeProb"]), 1.0)

    def test_age_band_parsing_and_grouping(self) -> None:
        self.assertEqual(parse_age_band("0-4"), (0, 4))
        self.assertEqual(parse_age_band("85+", age85_plus_max=89), (85, 89))
        self.assertEqual(parse_age_band(">=15"), (15, 15))
        self.assertEqual(broad_age_group_from_years([0, 4, 5, 14, 15]), [1, 1, 2, 2, 3])

    def test_expand_age_distribution_prefers_smoothed_proportion(self) -> None:
        import pandas as pd

        df = pd.DataFrame(
            [
                {"age_group": "0-4", "proportion": 0.9, "smoothed proportion": 0.5},
                {"age_group": "5-9", "proportion": 0.1, "smoothed proportion": 0.5},
            ]
        )
        ages, probs = expand_age_distribution_table(df)

        self.assertEqual(ages, list(range(10)))
        self.assertAlmostEqual(sum(probs), 1.0)
        self.assertAlmostEqual(probs[0], 0.1)

    def test_risk_prevalence_overrides_scalar_and_vector(self) -> None:
        pars = load_tb_csv("abm/default_data.csv")

        updated = apply_risk_prevalence_overrides(
            pars, {"smoking": 0.2, "diabetes": [0.0, 0.1, 0.3]}
        )

        self.assertEqual(updated["smokingPrevByAge"], [0.2, 0.2, 0.2])
        self.assertEqual(updated["diabetesPrevByAge"], [0.0, 0.1, 0.3])

    def test_invalid_risk_prevalence_override_raises(self) -> None:
        pars = load_tb_csv("abm/default_data.csv")

        with self.assertRaises(ValueError):
            apply_risk_prevalence_overrides(pars, {"smoking": [0.1, 0.2]})

    def test_disease_or_overrides_work_and_reject_invalid_values(self) -> None:
        pars = load_tb_csv("abm/default_data.csv")

        updated = apply_disease_or_overrides(pars, {"renal": 4.0})

        self.assertEqual(updated["disOR"]["renal"], 4.0)
        with self.assertRaises(ValueError):
            apply_disease_or_overrides(pars, {"renal": 0})

    def test_load_parameters_from_config_refreshes_totals(self) -> None:
        pars = load_parameters_from_config(build_default_config())

        self.assertAlmostEqual(sum(pars["exactAgeProb"]), 1.0)
        self.assertAlmostEqual(sum(pars["popFrac"]), 1.0)
        self.assertAlmostEqual(pars["totalFemalePrev"], 0.56)
        self.assertAlmostEqual(pars["totalBCGPrev"], 0.68)

    def test_committed_reference_configs_csv_file_can_be_consumed(self) -> None:
        reference_root = Path("validation") / "matlab_reference"
        config_paths = sorted(reference_root.glob("*/scenario_config.json"))
        if not config_paths:
            self.skipTest("MATLAB reference scenario configs are unavailable.")

        for config_path in config_paths:
            with self.subTest(config_path=config_path):
                cfg = load_reference_scenario_config(config_path)
                resolved = resolve_repo_path(cfg["csvFile"])
                pars = load_tb_csv(cfg["csvFile"])

                self.assertEqual(resolved, REPO_ROOT / "abm" / "default_data.csv")
                self.assertEqual(pars["ageNames"], ["0-4", "5-14", ">=15"])


if __name__ == "__main__":
    unittest.main()
