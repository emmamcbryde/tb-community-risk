from __future__ import annotations

import sys
import unittest

from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend
from engine.apy.ltbi_state import enable_development_compatibility_mode


class PythonApyBackendTests(unittest.TestCase):
    def setUp(self) -> None:
        self.backend = PythonApyBackend(repo_root())

    def test_status_reports_python_backend(self) -> None:
        status = self.backend.status()

        self.assertEqual(status["name"], "python_apy")
        self.assertTrue(status["started"])
        self.assertTrue(status["experimental"])
        self.assertFalse(status["matlabRequired"])

    def test_default_config_returns_dict(self) -> None:
        config = self.backend.default_config()

        self.assertIsInstance(config, dict)
        self.assertEqual(config["modelVersion"], "v9")

    def test_validate_default_config_passes(self) -> None:
        report = self.backend.validate_config(self.backend.default_config())

        self.assertTrue(report["isValid"])
        self.assertIn("warnings", report)

    def test_run_scenario_bundle_small_run(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 100, "nReps": 5, "seed": 1})
        config = enable_development_compatibility_mode(config)

        bundle = self.backend.run_scenario_bundle(config)

        self.assertIn("metadata", bundle)
        self.assertIn("headline", bundle)
        self.assertIn("technical", bundle)
        self.assertIn("dynamicComparison", bundle["technical"])
        self.assertEqual(bundle["metadata"]["backend"], "python")
        self.assertEqual(
            bundle["technical"]["dynamicComparison"]["source"],
            "doNothing.derived",
        )
        self.assertEqual(bundle["metadata"]["modelType"], "agent_based")
        self.assertEqual(bundle["metadata"]["nReps"], 5)
        self.assertEqual(bundle["metadata"]["seed"], 1)
        self.assertEqual(bundle["technical"]["tableMetadata"]["rawRows"], 5)

    def test_expected_value_bundle_uses_deterministic_runner_not_replicates(self) -> None:
        config = self.backend.default_config()
        config.update(
            {
                "N": 60,
                "analysisMethod": "expected_value",
                "analysisMethodLabel": "Expected outcomes",
                "nReps": 5000,
                "seed": 999,
            }
        )
        config = enable_development_compatibility_mode(config)

        bundle = self.backend.run_scenario_bundle(config)

        self.assertEqual(bundle["metadata"]["modelType"], "expected_value")
        self.assertIsNone(bundle["metadata"]["nReps"])
        self.assertIsNone(bundle["metadata"]["seed"])
        self.assertEqual(bundle["technical"]["tableMetadata"]["rawRows"], 1)
        self.assertEqual(
            bundle["technical"]["eventLedger"]["metadata"]["modelType"],
            "expected_value",
        )

    def test_stochastic_bundle_executes_selected_repetition_count(self) -> None:
        config = self.backend.default_config()
        config.update(
            {
                "N": 60,
                "analysisMethod": "agent_based",
                "analysisMethodLabel": "Simulated community variation",
                "nReps": 3,
                "seed": 44,
            }
        )
        config = enable_development_compatibility_mode(config)

        bundle = self.backend.run_scenario_bundle(config)

        self.assertEqual(bundle["metadata"]["modelType"], "agent_based")
        self.assertEqual(bundle["metadata"]["nReps"], 3)
        self.assertEqual(bundle["metadata"]["seed"], 44)
        self.assertEqual(bundle["technical"]["tableMetadata"]["rawRows"], 3)

    def test_default_economics_config_returns_dict(self) -> None:
        config = self.backend.default_economics_config()

        self.assertIsInstance(config, dict)
        self.assertIn("metadata", config)
        self.assertIn("costs", config)

    def test_economics_preset_kwab150_returns_dict(self) -> None:
        config = self.backend.economics_preset_kwab150()

        self.assertEqual(config["metadata"]["currencyCode"], "AUD")
        self.assertEqual(config["metadata"]["priceYear"], "2025-26")
        self.assertEqual(
            config["costItems"][0]["conversionStatus"],
            "not_converted",
        )

    def test_economics_preset_dale2019_returns_working_default(self) -> None:
        config = self.backend.economics_preset_dale2019_aud("4R")

        self.assertEqual(config["metadata"]["presetName"], "Dale 2019 AUD working defaults")
        self.assertEqual(config["metadata"]["targetPriceYear"], "2019")
        adr = {item["costItemId"]: item for item in config["costItems"]}["tpt_adr_management"]
        self.assertAlmostEqual(adr["originalCost"], 23.141)

    def test_run_economics_for_config_small_run(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 100, "nReps": 5, "seed": 1})
        config = enable_development_compatibility_mode(config)
        economics_config = self.backend.economics_preset_kwab150()

        economics = self.backend.run_economics_for_config(config, economics_config)

        self.assertIn("summaryRows", economics)
        self.assertIn("status", economics)
        self.assertTrue(economics["summaryRows"])
        self.assertTrue(economics["status"]["validationReport"]["isValid"])

    def test_python_backend_does_not_import_matlab_engine(self) -> None:
        sys.modules.pop("matlab.engine", None)
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 1, "seed": 2})
        config = enable_development_compatibility_mode(config)

        self.backend.run_scenario_bundle(config)

        self.assertNotIn("matlab.engine", sys.modules)

    def test_python_backend_economics_does_not_import_matlab_engine(self) -> None:
        sys.modules.pop("matlab.engine", None)
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 2, "seed": 2})
        config = enable_development_compatibility_mode(config)

        self.backend.run_economics_for_config(
            config,
            self.backend.economics_preset_kwab150(),
        )

        self.assertNotIn("matlab.engine", sys.modules)


if __name__ == "__main__":
    unittest.main()
