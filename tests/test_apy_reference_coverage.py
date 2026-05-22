from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

from engine.apy.capabilities import (
    DOCUMENTATION_ONLY,
    PARTIAL,
    UNSUPPORTED,
    capability_rows,
)
from engine.apy.reference_coverage import CONTRACT_VERSION, get_reference_coverage


REFERENCE_ROOT = Path("validation") / "matlab_reference"


def _components_by_name(report: dict) -> dict[str, dict]:
    return {row["component"]: row for row in report["components"]}


class ApyReferenceCoverageTests(unittest.TestCase):
    def test_reference_coverage_report_is_json_serialisable(self) -> None:
        report = get_reference_coverage()

        self.assertEqual(report["contractVersion"], CONTRACT_VERSION)
        self.assertFalse(report["matlabRequired"])
        json.dumps(report, allow_nan=False, sort_keys=True)

    def test_every_capability_component_appears_exactly_once(self) -> None:
        report = get_reference_coverage()
        expected_components = [row["component"] for row in capability_rows()]
        reported_components = [row["component"] for row in report["components"]]

        self.assertEqual(sorted(reported_components), sorted(expected_components))
        self.assertEqual(len(reported_components), len(set(reported_components)))

    def test_committed_matlab_reference_suite_and_fixture_dirs_are_detected(self) -> None:
        report = get_reference_coverage()
        fixture_dirs = set(report["referenceFixtureDirs"])

        self.assertEqual(report["referenceSuite"], str(REFERENCE_ROOT / "scenario_suite_v1.json"))
        self.assertTrue(report["referenceSuiteExists"])
        self.assertGreaterEqual(report["referenceFixtureCount"], 8)
        self.assertIn("default_random_igra_3hp_N1500_seed1", fixture_dirs)
        self.assertIn("random_tst_3hp_N1500_seed1", fixture_dirs)

        matlab_reference_fixtures = _components_by_name(report)["matlab_reference_fixtures"]
        self.assertTrue(matlab_reference_fixtures["matlabReferenceFixturesExist"])
        self.assertTrue(matlab_reference_fixtures["parityTestsExist"])
        self.assertTrue(matlab_reference_fixtures["distributionalTestsExist"])

    def test_partial_runtime_helpers_are_not_marked_complete(self) -> None:
        report = get_reference_coverage()
        by_component = _components_by_name(report)

        partial_components = {
            "attributable_risk",
            "attributable_risk_analysis",
            "chart_rendering_parity",
            "economics",
            "export_display_helpers",
            "exports",
            "health_economics",
            "natural_history",
            "natural_history_addon_reporting",
            "targeting_compare",
            "targeting_strategy_comparison",
        }

        for component in partial_components:
            with self.subTest(component=component):
                self.assertEqual(by_component[component]["migrationStatus"], PARTIAL)
                self.assertEqual(by_component[component]["coverageStatus"], PARTIAL)

    def test_documentation_only_matlab_ui_helpers_are_not_runtime_requirements(self) -> None:
        report = get_reference_coverage()
        by_component = _components_by_name(report)

        for component in {"matlab_app_callbacks", "matlab_app_designer_helpers"}:
            with self.subTest(component=component):
                row = by_component[component]
                self.assertEqual(row["migrationStatus"], DOCUMENTATION_ONLY)
                self.assertEqual(row["coverageStatus"], DOCUMENTATION_ONLY)
                self.assertIsNone(row["pythonEntryPoint"])
                self.assertFalse(row["parityTestsExist"])
                self.assertFalse(row["distributionalTestsExist"])
                self.assertFalse(row["pythonOnlyTestsExist"])

    def test_importing_and_calling_report_does_not_import_matlab_engine(self) -> None:
        code = (
            "import importlib, json, sys; "
            "before = {'matlab': 'matlab' in sys.modules, "
            "'matlab.engine': 'matlab.engine' in sys.modules}; "
            "module = importlib.import_module('engine.apy.reference_coverage'); "
            "module.get_reference_coverage(); "
            "after = {'matlab': 'matlab' in sys.modules, "
            "'matlab.engine': 'matlab.engine' in sys.modules}; "
            "print(json.dumps({'before': before, 'after': after}, sort_keys=True))"
        )
        completed = subprocess.run(
            [sys.executable, "-c", code],
            check=True,
            capture_output=True,
            text=True,
        )
        result = json.loads(completed.stdout)

        if not result["before"]["matlab"]:
            self.assertFalse(result["after"]["matlab"])
        if not result["before"]["matlab.engine"]:
            self.assertFalse(result["after"]["matlab.engine"])

    def test_unsupported_component_is_not_marked_complete(self) -> None:
        unsupported_capability = {
            "component": "unsupported_probe",
            "label": "Unsupported probe",
            "status": UNSUPPORTED,
            "pythonEntryPoint": None,
            "matlabReference": None,
            "notes": "Test-only unsupported component.",
        }

        with patch(
            "engine.apy.reference_coverage.capability_rows",
            return_value=[unsupported_capability],
        ):
            report = get_reference_coverage()

        row = report["components"][0]
        self.assertEqual(row["migrationStatus"], UNSUPPORTED)
        self.assertEqual(row["coverageStatus"], UNSUPPORTED)
        self.assertNotEqual(row["coverageStatus"], "complete")


if __name__ == "__main__":
    unittest.main()
