from __future__ import annotations

import sys
import unittest

from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend


class PythonApyBackendEconomicsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.backend = PythonApyBackend(repo_root())

    def test_default_economics_config_returns_dict(self) -> None:
        config = self.backend.default_economics_config()

        self.assertIsInstance(config, dict)
        self.assertIn("metadata", config)
        self.assertIn("costs", config)

    def test_kwab150_preset_returns_dict(self) -> None:
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

        economics = self.backend.run_economics_for_config(
            config,
            self.backend.economics_preset_kwab150(),
        )

        self.assertIn("summaryRows", economics)
        self.assertIn("status", economics)
        self.assertTrue(economics["summaryRows"])

    def test_python_economics_does_not_import_matlab_engine(self) -> None:
        sys.modules.pop("matlab.engine", None)
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 2, "seed": 4})

        self.backend.run_economics_for_config(
            config,
            self.backend.economics_preset_kwab150(),
        )

        self.assertNotIn("matlab.engine", sys.modules)


if __name__ == "__main__":
    unittest.main()
