from __future__ import annotations

import unittest

from app.demographic_profile import demographic_profile_hash, restore_apy_demographic_defaults
from app.parameter_workspace import apply_parameter_workspace, unified_default_session_state
from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend
from engine.apy.age_distribution import (
    broad_age_group_from_years,
    expand_age_distribution_table,
    parse_age_band,
)
from engine.apy.calibration import calibrate_from_config
from engine.apy.cohort import draw_base_population, make_rng
from engine.apy.config import build_default_config
from engine.apy.data import (
    apply_disease_or_overrides,
    apply_risk_prevalence_overrides,
    load_age_distribution,
    load_parameters_from_config,
    load_tb_csv,
    resolve_age_distribution_file,
)


class ApyDataLoadingTests(unittest.TestCase):
    def _compat_config(self) -> dict:
        config = build_default_config()
        config["ltbiStateAssumptions"]["developmentCompatibilityMode"] = True
        config["ltbiStateAssumptions"]["baselineRecentLTBIProportion"] = 0.0
        return config

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

    def test_unified_defaults_resolve_complete_apy_demographic_profile(self) -> None:
        config = unified_default_session_state()["config"]

        pars = load_parameters_from_config(config)

        self.assertTrue(str(pars["ageDistributionSource"]).endswith("abm\\default_age_distribution.csv") or str(pars["ageDistributionSource"]).endswith("abm/default_age_distribution.csv"))
        self.assertEqual(len(pars["exactAgeValues"]), 90)
        self.assertEqual(min(pars["exactAgeValues"]), 0)
        self.assertEqual(max(pars["exactAgeValues"]), 89)
        self.assertAlmostEqual(sum(pars["exactAgeProb"]), 1.0)
        self.assertTrue(all(value >= 0 for value in pars["exactAgeProb"]))
        self.assertAlmostEqual(sum(pars["popFrac"]), 1.0)
        self.assertGreater(pars["totalDiabetesPrev"], 0)
        self.assertGreater(pars["totalCurrentSmokerPrev"], 0)

    def test_scenario_export_reload_preserves_resolved_demographic_profile(self) -> None:
        backend = PythonApyBackend(repo_root())
        state = unified_default_session_state()
        before = demographic_profile_hash(state["config"])

        tmp = repo_root() / ".tmp_demographic_roundtrip"
        tmp.mkdir(exist_ok=True)
        path = tmp / "scenario.json"
        try:
            backend.save_scenario(state["config"], path, state["economics_config"])
            loaded, _, info = backend.load_scenario(str(path))
        finally:
            if path.exists():
                path.unlink()
            if tmp.exists():
                tmp.rmdir()

        self.assertIn("economics", info)
        self.assertEqual(demographic_profile_hash(loaded), before)

    def test_supplied_age_distribution_changes_age_sensitive_intermediate(self) -> None:
        base = self._compat_config()
        young = dict(base)
        old = dict(base)
        young["ageDistributionFile"] = ""
        old["ageDistributionFile"] = ""
        young["ageDistributionTable"] = [{"age": "0-4", "proportion": 1.0}]
        old["ageDistributionTable"] = [{"age": "15+", "proportion": 1.0}]
        calib = calibrate_from_config(base)

        young_pop = draw_base_population(
            1000,
            load_parameters_from_config(young),
            calib["ageInfLogLambda"],
            calib["ageInfGamma"],
            make_rng(123),
        )
        old_pop = draw_base_population(
            1000,
            load_parameters_from_config(old),
            calib["ageInfLogLambda"],
            calib["ageInfGamma"],
            make_rng(123),
        )

        self.assertLess(max(young_pop["ageYears"]), 5)
        self.assertGreaterEqual(min(old_pop["ageYears"]), 15)
        self.assertLess(float(young_pop["pInfection"].mean()), float(old_pop["pInfection"].mean()))

    def test_fixed_seed_reproducibly_draws_demographic_profile(self) -> None:
        config = self._compat_config()
        pars = load_parameters_from_config(config)
        calib = calibrate_from_config(config)

        first = draw_base_population(
            200,
            pars,
            calib["ageInfLogLambda"],
            calib["ageInfGamma"],
            make_rng(99),
        )
        second = draw_base_population(
            200,
            pars,
            calib["ageInfLogLambda"],
            calib["ageInfGamma"],
            make_rng(99),
        )

        self.assertEqual(first["ageYears"].tolist(), second["ageYears"].tolist())
        self.assertEqual(first["infected"].tolist(), second["infected"].tolist())

    def test_economic_only_workspace_edit_does_not_change_demographic_profile(self) -> None:
        state = unified_default_session_state()
        before = demographic_profile_hash(state["config"])
        rows = [dict(row) for row in state["parameter_workspace"]["rows"]]
        next(row for row in rows if row["parameterId"] == "cost.return_for_results")["currentValue"] = 75.0

        config, _ = apply_parameter_workspace(state["config"], state["economics_config"], rows)

        self.assertEqual(demographic_profile_hash(config), before)

    def test_restore_apy_demographic_defaults_does_not_reset_strategy(self) -> None:
        config = build_default_config()
        config["testType"] = "TST"
        config["regimen"] = "4R"
        config["screenCoverage"] = 0.45
        config["screeningWindowYears"] = 2
        config["screenWindow"] = 2
        config["ageDistributionFile"] = ""
        config["ageDistributionTable"] = [{"age": "0-4", "proportion": 1.0}]
        config["riskPrev"]["diabetes"] = 0.9

        restored = restore_apy_demographic_defaults(config)

        self.assertEqual(restored["testType"], "TST")
        self.assertEqual(restored["regimen"], "4R")
        self.assertEqual(restored["screenCoverage"], 0.45)
        self.assertEqual(restored["screeningWindowYears"], 2)
        self.assertEqual(restored["ageDistributionFile"], "abm/default_age_distribution.csv")
        self.assertEqual(restored["ageDistributionTable"], [])
        self.assertIsNone(restored["riskPrev"]["diabetes"])


if __name__ == "__main__":
    unittest.main()
