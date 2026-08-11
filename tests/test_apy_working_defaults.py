from __future__ import annotations

from io import BytesIO
import unittest

from openpyxl import load_workbook

from app.parameter_workspace import (
    HIGHER_BURDEN_DALY_OUTCOME_PRESET,
    PARAMETER_GROUPS,
    PRIMARY_DALY_OUTCOME_PRESET,
    apply_parameter_workspace,
    build_parameter_workspace,
    reset_all_parameters,
    reset_parameter_group,
    unified_default_session_state,
    validate_parameter_workspace,
)
from app.results_workbook import build_results_workbook
from engine.apy.working_defaults import (
    UNIFIED_WORKING_DEFAULT_LABEL,
    UNIFIED_WORKING_DEFAULT_PRESET_ID,
    build_unified_working_default_preset,
)


class APYWorkingDefaultsTests(unittest.TestCase):
    def test_unified_preset_loads_all_components_atomically(self) -> None:
        state = unified_default_session_state()

        config = state["config"]
        econ = state["economics_config"]
        workspace = state["parameter_workspace"]

        self.assertEqual(config["workingDefaultPresetId"], UNIFIED_WORKING_DEFAULT_PRESET_ID)
        self.assertEqual(config["populationPresetId"], "apy_demonstration")
        self.assertEqual(config["testType"], "IGRA")
        self.assertEqual(config["regimen"], "3HP")
        self.assertEqual(config["screeningWindowYears"], 3)
        self.assertEqual(config["followUpHorizonYears"], 20)
        self.assertEqual(config["analysisMethod"], "expected_value")
        self.assertEqual(econ["metadata"]["presetName"], "Dale 2019 AUD working defaults")
        self.assertEqual(econ["metadata"]["targetCurrency"], "AUD")
        self.assertEqual(econ["metadata"]["targetPriceYear"], "2019")
        self.assertEqual(econ["dalyAssumptions"]["outcomePreset"], PRIMARY_DALY_OUTCOME_PRESET)
        self.assertEqual(workspace["presetLabel"], UNIFIED_WORKING_DEFAULT_LABEL)
        self.assertEqual(set(workspace["groups"]), set(PARAMETER_GROUPS))

    def test_parameter_groups_and_changed_from_default_tracking(self) -> None:
        state = unified_default_session_state()
        rows = state["parameter_workspace"]["rows"]
        edited = [dict(row) for row in rows]
        population_row = next(row for row in edited if row["parameterId"] == "demography.population_size")
        population_row["currentValue"] = 12345

        config, econ = apply_parameter_workspace(state["config"], state["economics_config"], edited)
        workspace = build_parameter_workspace(config, econ)

        self.assertEqual(config["N"], 12345)
        changed_ids = {row["parameterId"] for row in workspace["rows"] if row["changedFromDefault"]}
        self.assertIn("demography.population_size", changed_ids)
        self.assertIn("Demography", {row["group"] for row in workspace["rows"]})
        self.assertIn("TB epidemiology", {row["group"] for row in workspace["rows"]})
        self.assertIn("Health-service provision and intervention", {row["group"] for row in workspace["rows"]})
        self.assertIn("Costs and outcomes", {row["group"] for row in workspace["rows"]})
        self.assertIn("Analysis settings", {row["group"] for row in workspace["rows"]})

    def test_group_reset_and_reset_all_restore_defaults(self) -> None:
        state = unified_default_session_state()
        rows = [dict(row) for row in state["parameter_workspace"]["rows"]]
        next(row for row in rows if row["parameterId"] == "demography.population_size")["currentValue"] = 999
        next(row for row in rows if row["parameterId"] == "service.coverage")["currentValue"] = 0.75

        demography_reset = reset_parameter_group(rows, "Demography")
        self.assertEqual(
            next(row for row in demography_reset if row["parameterId"] == "demography.population_size")["currentValue"],
            next(row for row in demography_reset if row["parameterId"] == "demography.population_size")["defaultValue"],
        )
        self.assertEqual(next(row for row in demography_reset if row["parameterId"] == "service.coverage")["currentValue"], 0.75)

        all_reset = reset_all_parameters(rows)
        self.assertTrue(all(row["currentValue"] == row["defaultValue"] for row in all_reset))

    def test_dale_costs_and_local_sa_health_costs_are_populated(self) -> None:
        econ = build_unified_working_default_preset()["economicsConfig"]
        items = {item["costItemId"]: item for item in econ["costItems"]}

        self.assertEqual(items["test_igra"]["originalCost"], 113.48)
        self.assertEqual(items["test_tst"]["originalCost"], 116.07)
        self.assertEqual(items["regimen_3hp"]["originalCost"], 165.5072)
        self.assertEqual(items["regimen_4r"]["originalCost"], 123.3172)
        self.assertEqual(items["regimen_3hr"]["originalCost"], 134.2272)
        self.assertEqual(items["regimen_6h"]["originalCost"], 187.7508)
        self.assertEqual(items["regimen_9h"]["originalCost"], 254.8544)
        self.assertEqual(items["active_tb_disease"]["originalCost"], 19079.60)
        self.assertEqual(items["program_setup"]["originalCost"], 0.0)
        self.assertEqual(items["program_running"]["originalCost"], 0.0)
        self.assertEqual(econ["localSAHealthPathwayCosts"]["returnForResults"]["value"], 50.0)
        self.assertEqual(econ["localSAHealthPathwayCosts"]["clinicalReview"]["value"], 50.0)
        self.assertEqual(econ["localSAHealthPathwayCosts"]["activeTBExclusionWorkup"]["value"], 150.0)

    def test_programme_setup_base_and_new_programme_options_apply(self) -> None:
        state = unified_default_session_state()
        rows = [dict(row) for row in state["parameter_workspace"]["rows"]]
        setup = next(row for row in rows if row["parameterId"] == "cost.program_setup_preset")
        setup["currentValue"] = "New-programme implementation scenario"

        _, econ = apply_parameter_workspace(state["config"], state["economics_config"], rows)
        item = next(item for item in econ["costItems"] if item["costItemId"] == "program_setup")

        self.assertEqual(item["originalCost"], 500000.0)
        self.assertEqual(econ["localSAHealthPathwayCosts"]["selectedProgramSetupPreset"], "new_programme_implementation")

    def test_standard_workspace_uses_dalys_and_hides_qalys(self) -> None:
        workspace = unified_default_session_state()["parameter_workspace"]
        labels = "\n".join(str(row["label"]) for row in workspace["rows"])
        self.assertIn("Outcome preset", labels)
        self.assertIn("Active-TB disability weight", labels)
        self.assertNotIn("QALY", labels)

    def test_standard_workspace_uses_plain_language_targeting_labels(self) -> None:
        workspace = unified_default_session_state()["parameter_workspace"]
        targeting = next(row for row in workspace["rows"] if row["parameterId"] == "service.targeting")

        self.assertEqual(targeting["currentValue"], "Prioritise people most likely to avoid active TB")
        self.assertNotEqual(targeting["currentValue"], "prevent")
        self.assertIn("No risk-based prioritisation", targeting["selectOptions"])

    def test_invalid_select_override_is_rejected(self) -> None:
        workspace = unified_default_session_state()["parameter_workspace"]
        targeting = next(row for row in workspace["rows"] if row["parameterId"] == "service.targeting")
        targeting["currentValue"] = "invented targeting option"

        validation = validate_parameter_workspace(workspace["rows"])

        self.assertFalse(validation["isValid"])
        self.assertIn("service.targeting", {row["parameterId"] for row in validation["messages"]})

    def test_higher_burden_option_adds_post_tb_dalys_without_combining_qalys(self) -> None:
        state = unified_default_session_state()
        rows = [dict(row) for row in state["parameter_workspace"]["rows"]]
        outcome = next(row for row in rows if row["parameterId"] == "outcome.daly_preset")
        outcome["currentValue"] = HIGHER_BURDEN_DALY_OUTCOME_PRESET

        _, econ = apply_parameter_workspace(state["config"], state["economics_config"], rows)
        daly = econ["dalyAssumptions"]

        self.assertTrue(daly["includePostTBSequelae"])
        self.assertEqual(daly["postTBDALYsPerActiveTBCase"]["value"], 5.8)
        self.assertEqual(daly["standardOutcomeMetric"], "DALYs")
        self.assertFalse(daly["qalyStandardInterfaceVisible"])

    def test_expected_outcomes_mode_does_not_require_simulation_count(self) -> None:
        workspace = unified_default_session_state()["parameter_workspace"]
        n_reps = next(row for row in workspace["rows"] if row["parameterId"] == "analysis.n_reps")
        n_reps["currentValue"] = ""

        validation = validate_parameter_workspace(workspace["rows"])
        self.assertTrue(validation["isValid"])

    def test_quick_simulation_is_labelled_preview_and_sets_five_replicates(self) -> None:
        state = unified_default_session_state()
        rows = [dict(row) for row in state["parameter_workspace"]["rows"]]
        next(row for row in rows if row["parameterId"] == "analysis.method")["currentValue"] = "Simulated community variation"
        next(row for row in rows if row["parameterId"] == "analysis.simulation_mode")["currentValue"] = "Quick preview: 5 simulations"

        config, _ = apply_parameter_workspace(state["config"], state["economics_config"], rows)

        self.assertEqual(config["analysisMethod"], "agent_based")
        self.assertEqual(config["simulationMode"], "quick_preview")
        self.assertEqual(config["nReps"], 5)
        self.assertIn("preview", config["simulationModeLabel"].lower())

    def test_workbook_exports_defaults_and_overrides(self) -> None:
        state = unified_default_session_state()
        rows = [dict(row) for row in state["parameter_workspace"]["rows"]]
        next(row for row in rows if row["parameterId"] == "service.coverage")["currentValue"] = 0.4
        config, econ = apply_parameter_workspace(state["config"], state["economics_config"], rows)

        payload = build_results_workbook(
            config=config,
            economics_config=econ,
            bundle={"metadata": {}, "headline": {}, "technical": {}},
        )
        wb = load_workbook(BytesIO(payload), read_only=True)
        self.assertIn("Parameter_workspace", wb.sheetnames)
        self.assertIn("Parameter_defaults", wb.sheetnames)
        self.assertIn("Parameter_overrides", wb.sheetnames)
        rows = list(wb["Parameter_overrides"].iter_rows(values_only=True))
        headers = rows[0]
        data = [dict(zip(headers, row)) for row in rows[1:]]
        coverage = next(row for row in data if row["parameterId"] == "service.coverage")
        self.assertEqual(coverage["currentValue"], 0.4)
        self.assertTrue(coverage["changedFromDefault"])
        wb.close()


if __name__ == "__main__":
    unittest.main()
