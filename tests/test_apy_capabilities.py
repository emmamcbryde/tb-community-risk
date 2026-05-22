from __future__ import annotations

import json
import unittest

from engine.apy.capabilities import (
    COMPLETE,
    DOCUMENTATION_ONLY,
    MATLAB_REFERENCE,
    PARTIAL,
    REFERENCE_ONLY,
    UNSUPPORTED,
    UNKNOWN,
    capability_rows,
    get_apy_capabilities,
    get_component_status,
    unsupported_or_partial_components,
)


ALLOWED_CAPABILITY_STATUSES = {
    COMPLETE,
    PARTIAL,
    REFERENCE_ONLY,
    DOCUMENTATION_ONLY,
    UNSUPPORTED,
    MATLAB_REFERENCE,
}

HUMAN_SCOPE_COMPONENTS = {
    "matlab_reference_engine",
    "default_config",
    "config_validation",
    "data_loading",
    "age_distribution",
    "regimen_handling",
    "calibration",
    "cohort_primitives",
    "core_simulation",
    "scenario_runner",
    "scenario_save_load",
    "summary_rows",
    "results_bundle",
    "exports",
    "export_display_helpers",
    "chart_rendering_parity",
    "natural_history",
    "do_nothing_dynamic_comparison",
    "natural_history_addon_reporting",
    "targeting_compare",
    "targeting_strategy_comparison",
    "economics",
    "health_economics",
    "attributable_risk",
    "attributable_risk_analysis",
    "matlab_reference_fixtures",
    "distributional_validation",
    "matlab_user_options_compatibility",
    "matlab_app_callbacks",
    "matlab_app_designer_helpers",
}


class ApyCapabilitiesTests(unittest.TestCase):
    def test_capability_contract_is_json_serialisable(self) -> None:
        capabilities = get_apy_capabilities()

        self.assertEqual(
            capabilities["contractVersion"],
            "apy_python_migration_capabilities_v1",
        )
        self.assertEqual(capabilities["backend"], "python_apy")
        self.assertFalse(capabilities["matlabRequired"])
        self.assertEqual(capabilities["overallStatus"], PARTIAL)
        self.assertEqual(
            capabilities["statusValues"],
            [
                "complete",
                "partial",
                "reference_only",
                "documentation_only",
                "unsupported",
                "matlab_reference",
            ],
        )
        self.assertGreater(len(capabilities["components"]), 0)
        json.dumps(capabilities, allow_nan=False, sort_keys=True)

    def test_component_rows_have_stable_shape_and_known_statuses(self) -> None:
        rows = capability_rows()

        self.assertGreater(len(rows), 0)
        for row in rows:
            self.assertEqual(
                set(row),
                {
                    "component",
                    "label",
                    "status",
                    "pythonEntryPoint",
                    "matlabReference",
                    "notes",
                },
            )
            self.assertIn(row["status"], ALLOWED_CAPABILITY_STATUSES)

        by_component = {row["component"]: row for row in rows}
        self.assertEqual(by_component["matlab_reference_engine"]["status"], MATLAB_REFERENCE)
        self.assertEqual(by_component["default_config"]["status"], COMPLETE)
        self.assertEqual(by_component["core_simulation"]["status"], PARTIAL)
        self.assertEqual(by_component["economics"]["status"], PARTIAL)
        self.assertEqual(by_component["health_economics"]["status"], PARTIAL)
        self.assertEqual(by_component["attributable_risk"]["status"], PARTIAL)
        self.assertEqual(by_component["matlab_app_callbacks"]["status"], DOCUMENTATION_ONLY)

    def test_human_scope_components_are_listed_with_allowed_statuses(self) -> None:
        by_component = {row["component"]: row for row in capability_rows()}

        self.assertEqual(HUMAN_SCOPE_COMPONENTS.difference(by_component), set())
        for component in HUMAN_SCOPE_COMPONENTS:
            self.assertIn(
                by_component[component]["status"],
                ALLOWED_CAPABILITY_STATUSES,
            )

    def test_required_components_have_exact_classifications(self) -> None:
        by_component = {row["component"]: row for row in capability_rows()}

        expected = {
            "matlab_reference_engine": MATLAB_REFERENCE,
            "default_config": COMPLETE,
            "config_validation": COMPLETE,
            "data_loading": COMPLETE,
            "age_distribution": COMPLETE,
            "regimen_handling": COMPLETE,
            "calibration": COMPLETE,
            "cohort_primitives": COMPLETE,
            "core_simulation": PARTIAL,
            "scenario_runner": PARTIAL,
            "scenario_save_load": PARTIAL,
            "summary_rows": COMPLETE,
            "results_bundle": PARTIAL,
            "exports": PARTIAL,
            "export_display_helpers": PARTIAL,
            "chart_rendering_parity": PARTIAL,
            "natural_history": PARTIAL,
            "do_nothing_dynamic_comparison": PARTIAL,
            "natural_history_addon_reporting": PARTIAL,
            "targeting_compare": PARTIAL,
            "targeting_strategy_comparison": PARTIAL,
            "economics": PARTIAL,
            "health_economics": PARTIAL,
            "attributable_risk": PARTIAL,
            "attributable_risk_analysis": PARTIAL,
            "matlab_reference_fixtures": REFERENCE_ONLY,
            "distributional_validation": REFERENCE_ONLY,
            "matlab_user_options_compatibility": REFERENCE_ONLY,
            "matlab_app_callbacks": DOCUMENTATION_ONLY,
            "matlab_app_designer_helpers": DOCUMENTATION_ONLY,
        }

        self.assertEqual(set(expected).difference(by_component), set())
        for component, status in expected.items():
            self.assertEqual(by_component[component]["status"], status)

    def test_get_component_status_returns_unknown_for_missing_component(self) -> None:
        self.assertEqual(get_component_status("default_config"), COMPLETE)
        self.assertEqual(get_component_status("core_simulation"), PARTIAL)
        self.assertEqual(get_component_status("does_not_exist"), UNKNOWN)

    def test_unsupported_or_partial_components_filters_to_named_statuses(self) -> None:
        rows = unsupported_or_partial_components()
        components = {row["component"] for row in rows}

        self.assertIn("core_simulation", components)
        self.assertIn("economics", components)
        self.assertIn("attributable_risk", components)
        self.assertNotIn("default_config", components)
        self.assertNotIn("matlab_reference_engine", components)
        self.assertNotIn("matlab_reference_fixtures", components)
        self.assertNotIn("matlab_app_callbacks", components)
        self.assertNotIn("matlab_app_designer_helpers", components)
        self.assertTrue(
            all(row["status"] in {PARTIAL, UNSUPPORTED} for row in rows)
        )

    def test_helpers_return_copies(self) -> None:
        capabilities = get_apy_capabilities()
        capabilities["components"][0]["status"] = "mutated"
        rows = capability_rows()
        rows[0]["status"] = "mutated"

        self.assertEqual(get_component_status("default_config"), COMPLETE)
        self.assertEqual(
            get_component_status("matlab_reference_engine"),
            MATLAB_REFERENCE,
        )


if __name__ == "__main__":
    unittest.main()
