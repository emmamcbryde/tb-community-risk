from __future__ import annotations

import sys
import unittest

from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend


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

    def test_run_economics_for_config_small_run(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 100, "nReps": 5, "seed": 1})
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

        self.backend.run_scenario_bundle(config)

        self.assertNotIn("matlab.engine", sys.modules)

    def test_python_backend_economics_does_not_import_matlab_engine(self) -> None:
        sys.modules.pop("matlab.engine", None)
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 2, "seed": 2})

        self.backend.run_economics_for_config(
            config,
            self.backend.economics_preset_kwab150(),
        )

        self.assertNotIn("matlab.engine", sys.modules)


if __name__ == "__main__":
    unittest.main()
