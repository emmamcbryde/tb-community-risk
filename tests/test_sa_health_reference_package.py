from __future__ import annotations

from copy import deepcopy
import csv
import json
from pathlib import Path
import unittest

import pyarrow as pa

from app.health_economics_inputs import editable_assumption_rows
from engine.apy.event_ledger_economics import run_event_ledger_health_economics
from engine.apy.sa_health_reference_package import (
    PROGRAMME_NOT_COSTED_LABEL,
    SA_HEALTH_REFERENCE_PACKAGE_ID,
    build_sa_health_reference_package,
    build_same_ledger_economic_scenario_comparison,
    write_sa_health_reference_package,
    _replicate_row,
)
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]


class SAHealthReferencePackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.package = build_sa_health_reference_package(
            n_reps=50,
            include_technical_recent_remote=True,
        )

    def test_manifest_is_complete_and_provisional(self) -> None:
        manifest = self.package["manifest"]

        self.assertEqual(manifest["packageId"], SA_HEALTH_REFERENCE_PACKAGE_ID)
        self.assertEqual(manifest["modelType"], "agent_based")
        self.assertEqual(manifest["epidemiologicalAnchor"]["naturalHistorySemantics"], "matlab_v9_implicit_early_late")
        self.assertFalse(manifest["dynamicTransmissionIncluded"])
        self.assertEqual(manifest["eventLedgerContractVersion"], "ltbi_screening_event_ledger_v3")
        self.assertEqual(manifest["healthEconomicsContractVersion"], "ltbi_health_economics_results_v3")
        self.assertTrue(manifest["configurationHash"])
        self.assertTrue(manifest["economicsConfigurationHash"])
        self.assertTrue(manifest["evidenceRegistryHash"])
        self.assertFalse(manifest["readiness"]["overallClinicianReady"])
        self.assertIn("The inherited 10/770 active-TB calibration target remains unresolved and must not be described as validated incident progression from LTBI.", manifest["interpretationGuardrails"])

    def test_reference_config_matches_sa_health_decisions(self) -> None:
        config = self.package["config"]
        econ = self.package["economicsConfig"]

        self.assertEqual(config["N"], 1500)
        self.assertEqual(config["testType"], "IGRA")
        self.assertEqual(config["regimen"], "3HP")
        self.assertEqual(config["screeningStrategy"], "prevent")
        self.assertEqual(config["screenCoverage"], 0.30)
        self.assertEqual(config["screeningWindowYears"], 2)
        self.assertEqual(config["followUpHorizonYears"], 20)
        self.assertEqual(config["nReps"], 50)
        self.assertEqual(config["seed"], 1)
        self.assertEqual(config["naturalHistorySemantics"], "matlab_v9_implicit_early_late")
        self.assertEqual(econ["metadata"]["perspective"], "Australian health-care system")
        self.assertEqual(econ["metadata"]["targetCurrency"], "AUD")
        self.assertEqual(econ["metadata"]["targetPriceYear"], "2019")
        self.assertIsNone(econ["threshold"]["value"])

    def test_primary_anchor_uses_matlab_v9_compatible_hazards(self) -> None:
        calibration = self.package["epidemiology"]["calibration"]

        self.assertAlmostEqual(calibration["lambdaEarly"], 0.0037975986627, places=10)
        self.assertAlmostEqual(calibration["lambdaLate"], 0.0007595197325, places=10)
        self.assertNotAlmostEqual(calibration["lambdaEarly"], 0.0189879933116, places=5)

    def test_streamlit_cost_workspace_contains_requested_editable_costs(self) -> None:
        rows = editable_assumption_rows(economics_config=self.package["economicsConfig"])
        by_id = {row["assumptionId"]: row for row in rows}
        required = {
            "cost.test_igra",
            "cost.test_tst",
            "cost.regimen_3hp",
            "cost.regimen_4r",
            "cost.regimen_3hr",
            "cost.regimen_6h",
            "cost.regimen_9h",
            "cost.active_tb_disease",
            "cost.tpt_adr_management",
            "cost.false_positive_incremental",
            "cost.return_for_results",
            "cost.clinical_review",
            "cost.active_tb_exclusion_workup",
            "cost.program_setup",
            "cost.program_running",
            "cost.travel_outreach_staff_support",
        }

        self.assertTrue(required.issubset(by_id))
        for assumption_id in required:
            row = by_id[assumption_id]
            self.assertIn("currentValue", row)
            self.assertTrue(row.get("unit") or row.get("costBasis"))
            self.assertIn("originalCurrency", row)
            self.assertIn("originalPriceYear", row)
            self.assertIn("sourceCitation", row)
            self.assertIn("reviewStatus", row)
            self.assertIn("inclusionStatus", row)

    def test_report_tables_reconcile_to_authoritative_outputs(self) -> None:
        tables = self.package["tables"]
        primary = _replicate_row(self.package["economics"], "primary")
        executive = {row["metric"]: row for row in tables["executive_summary"]}

        self.assertAlmostEqual(executive["People screened"]["value"], 450.0)
        self.assertAlmostEqual(
            executive["Comparator active TB cases"]["value"] - executive["Intervention active TB cases"]["value"],
            executive["Active TB cases averted"]["value"],
        )
        self.assertAlmostEqual(executive["Incremental cost"]["value"], primary["incrementalCost"])
        self.assertAlmostEqual(
            executive["Net health-system cost or saving per active TB case averted"]["value"],
            primary["costPerActiveTBCasePrevented"],
        )
        annual = tables["annual_budget_impact"]
        self.assertAlmostEqual(
            annual[-1]["cumulativeIncrementalCostDiscounted"],
            primary["incrementalCost"],
        )

        category_total = next(row for row in tables["cost_categories"] if row["categoryId"] == "total")
        self.assertAlmostEqual(category_total["incrementalDiscountedCost"], primary["incrementalCost"])

    def test_unresolved_programme_costs_are_distinguishable_from_zero(self) -> None:
        categories = {
            row["categoryId"]: row
            for row in self.package["tables"]["cost_categories"]
        }

        for category_id in ["program_setup", "program_running", "travel_outreach_staff"]:
            self.assertEqual(categories[category_id]["interpretation"], PROGRAMME_NOT_COSTED_LABEL)
            self.assertIsNotNone(categories[category_id]["interventionDiscountedCost"])

    def test_missing_wtp_keeps_nmb_unavailable(self) -> None:
        economics = self.package["economics"]
        primary = economics["replicateResults"][economics["replicateResults"]["discountProfile"] == "primary"].iloc[0]
        unresolved = {item["field"] for item in economics["unresolvedInputs"]}

        self.assertIn("threshold.value", unresolved)
        self.assertIsNone(primary["netMonetaryBenefit"])
        self.assertNotIn("probabilityPositiveNMB_fixedParameterSimulation", set(economics["summaries"]["metric"]))

    def test_economic_scenarios_use_same_event_ledger_and_change_costs_only(self) -> None:
        comparison = build_same_ledger_economic_scenario_comparison(
            self.package["epidemiology"],
            self.package["economicsConfig"],
        )
        self.assertFalse(comparison["epidemiologyRerunRequired"])
        scenarios = {row["scenarioId"]: row for row in comparison["rows"]}
        primary = scenarios["primary_working_reference"]
        setup = scenarios["illustrative_500k_setup"]
        bundled = scenarios["pathway_components_bundled"]
        higher = scenarios["higher_burden_post_tb_daly"]

        self.assertTrue(all(row["usesSameEventLedger"] for row in scenarios.values()))
        self.assertAlmostEqual(setup["incrementalCost"] - primary["incrementalCost"], 500000.0)
        self.assertLess(bundled["interventionCost"], primary["interventionCost"])
        self.assertGreater(higher["dalysAverted"], primary["dalysAverted"])

    def test_distributional_agreement_with_stored_matlab_reference(self) -> None:
        rows = {row["metric"]: row for row in self.package["tables"]["matlab_reference_validation"]}

        for metric in [
            "nScreened",
            "nFalsePositiveTests",
            "nStartTPT",
            "nCompleteTPT",
            "nADRstop",
            "nCuredInfection",
            "nPreventedActiveTB",
            "nActiveBy20y",
        ]:
            self.assertTrue(rows[metric]["withinMatlabInterval"], metric)

    def test_technical_recent_remote_scenario_is_not_primary_anchor(self) -> None:
        technical = self.package["technicalRecentRemoteScenario"]

        self.assertFalse(technical["primarySAHealthAnchor"])
        self.assertFalse(technical["suitableAsPrimaryEstimate"])
        self.assertEqual(technical["packageId"], "sa_health_apy_technical_recent_remote_igra_3hp_prevent_30pct")
        self.assertNotEqual(
            technical["epidemiology"]["eventLedger"]["metadata"]["modelType"],
            self.package["eventLedger"]["metadata"]["modelType"],
        )

    def test_gross_and_net_ratio_tables_are_separate(self) -> None:
        gross = {row["ratioId"]: row for row in self.package["tables"]["gross_delivery_ratios"]}
        net = {row["ratioId"]: row for row in self.package["tables"]["net_health_system_ratios"]}

        self.assertGreater(gross["person_screened"]["value"], 0)
        self.assertIn("before active-TB care offsets", gross["person_screened"]["interpretation"])
        self.assertIn("after active-TB care offsets", net["person_screened"]["interpretation"])
        self.assertNotEqual(gross["person_screened"]["value"], net["person_screened"]["value"])

    def test_health_economics_page_uses_workspace_and_same_ledger_scenarios(self) -> None:
        page = (ROOT / "pages" / "4_Economics.py").read_text(encoding="utf-8")

        self.assertIn("Economic result", page)
        self.assertIn("Change cost assumptions", page)
        self.assertLess(
            page.index('st.subheader("Economic result")'),
            page.index('st.expander("Change cost assumptions"'),
        )
        self.assertIn("Recalculate economics", page)
        self.assertIn("Economic changes reuse the current screening outcomes", page)
        self.assertIn("All scenarios reuse the same screening outcomes", page)
        self.assertIn("Gross delivery expenditure", page)
        self.assertIn("delivery_scenario_comparison_rows", page)

    def test_rendered_health_economics_page_leads_with_decision_result(self) -> None:
        app = AppTest.from_file(str(ROOT / "pages" / "4_Economics.py"), default_timeout=60)
        app.session_state["config"] = self.package["config"]
        app.session_state["results_bundle"] = {
            "metadata": {"scenarioLabel": "test_reference"},
            "technical": {"eventLedger": self.package["eventLedger"]},
        }
        app.session_state["economics_config"] = self.package["economicsConfig"]
        app.session_state["economics_results"] = self.package["economics"]
        app.session_state["results_stale"] = False

        app.run()

        self.assertFalse(app.exception)
        self.assertEqual(
            [item.value for item in app.subheader[:3]],
            ["Economic result", "Programme-delivery scenarios", "Cost breakdown and budget impact"],
        )
        self.assertIn("Change cost assumptions", [item.label for item in app.expander])
        self.assertIn("Limitations and methods", [item.label for item in app.expander])
        button_labels = [item.label for item in app.button]
        self.assertIn("Recalculate economics", button_labels)
        self.assertIn("Restore SA Health economic defaults", button_labels)
        self.assertIn(
            "Incremental cost-effectiveness plane",
            [item.value for item in app.markdown],
        )
        self.assertGreaterEqual(sum(1 for item in app if getattr(item, "type", None) == "arrow_vega_lite_chart"), 1)
        headline = app.dataframe[0].value
        self.assertIn("Dominant", str(headline.loc[headline["Result"] == "Economic result", "Value"].iloc[0]))
        self.assertIn("DALYs averted", set(headline["Result"]))
        scenarios = next(
            item.value for item in app.dataframe
            if "Scenario" in set(getattr(item.value, "columns", []))
        )
        self.assertEqual(
            list(scenarios["Scenario"]),
            ["No additional programme overhead entered"],
        )
        self.assertEqual(len(set(scenarios["Active TB averted"])), 1)
        self.assertIn("Quadrant", scenarios.columns)
        self.assertIn("Classification", scenarios.columns)
        self.assertIn("Arithmetic ICER", scenarios.columns)
        self.assertIn("Additional setup cost", scenarios.columns)
        self.assertIn("Programme-cost assumption", scenarios.columns)
        self.assertEqual(
            str(scenarios.loc[scenarios["Scenario"] == "No additional programme overhead entered", "Additional setup cost"].iloc[0]),
            "AUD 0",
        )
        self.assertEqual(
            str(scenarios.loc[scenarios["Scenario"] == "No additional programme overhead entered", "Quadrant"].iloc[0]),
            "Lower right",
        )
        self.assertIn(
            "Dominant",
            str(scenarios.loc[scenarios["Scenario"] == "No additional programme overhead entered", "Classification"].iloc[0]),
        )
        self.assertIn(
            "-AUD",
            str(scenarios.loc[scenarios["Scenario"] == "No additional programme overhead entered", "Arithmetic ICER"].iloc[0]),
        )
        for dataframe in app.dataframe:
            columns = set(getattr(dataframe.value, "columns", []))
            self.assertFalse({"Event ledger", "Analysis basis"}.issubset(columns))
        source = (ROOT / "pages" / "4_Economics.py").read_text(encoding="utf-8")
        self.assertIn("DALYs averted compared with business as usual", source)
        self.assertIn("Incremental cost compared with business as usual (AUD)", source)
        self.assertIn("Business as usual", source)
        self.assertIn("st.altair_chart(_cost_effectiveness_plane_chart", source)
        self.assertNotIn('"Event ledger": row.get("Event ledger"', source)

    def test_health_economics_programme_cost_controls_are_explicit(self) -> None:
        app = AppTest.from_file(str(ROOT / "pages" / "4_Economics.py"), default_timeout=60)
        app.session_state["config"] = self.package["config"]
        app.session_state["results_bundle"] = {
            "metadata": {"scenarioLabel": "test_reference"},
            "technical": {"eventLedger": self.package["eventLedger"]},
        }
        app.session_state["economics_config"] = self.package["economicsConfig"]
        app.session_state["economics_results"] = self.package["economics"]
        app.session_state["results_stale"] = False
        app.session_state["health_econ_delivery_scenarios"] = {
            "includeAdditionalProgramCosts": True,
            "standaloneSetupCost": 500000.0,
            "illustrativeSetupCost": 500000.0,
            "standaloneAnnualRunningCost": 0.0,
            "standaloneRunningYears": 2,
            "annualCostFirstYear": 0,
            "standaloneTravelOutreachCost": 0.0,
            "standaloneStaffSupportCost": 0.0,
            "sharedAttributionShare": 0.50,
        }

        app.run()

        self.assertFalse(app.exception)
        scenarios = next(
            item.value for item in app.dataframe
            if "Scenario" in set(getattr(item.value, "columns", []))
        )
        self.assertEqual(
            list(scenarios["Scenario"]),
            [
                "No additional programme overhead entered",
                "Standalone programme - user-defined additional costs",
                "Shared delivery - user-defined attributable share",
            ],
        )
        base = scenarios.iloc[0]
        standalone = scenarios.iloc[1]
        shared = scenarios.iloc[2]
        self.assertEqual(standalone["Additional setup cost"], "AUD 500,000")
        self.assertEqual(standalone["Attributed share"], "100%")
        self.assertEqual(shared["Additional setup cost"], "AUD 250,000")
        self.assertEqual(shared["Attributed share"], "50%")
        self.assertIn("Common setup AUD 500,000", shared["Programme-cost assumption"])
        self.assertEqual(base["DALYs averted"], standalone["DALYs averted"])
        self.assertEqual(base["DALYs averted"], shared["DALYs averted"])
        self.assertAlmostEqual(
            self._money_display_to_float(standalone["Incremental cost"])
            - self._money_display_to_float(base["Incremental cost"]),
            500000.0,
            delta=2.0,
        )
        self.assertAlmostEqual(
            self._money_display_to_float(shared["Incremental cost"])
            - self._money_display_to_float(base["Incremental cost"]),
            250000.0,
            delta=2.0,
        )
        self.assertTrue(
            any(
                "AUD 500,000 is illustrative and user-changeable" in item.value
                for item in app.caption
            )
        )

    def test_health_economics_icer_plane_shows_vertical_programme_cost_movement(self) -> None:
        setup_controls = self._programme_controls(setup=500000.0, annual=0.0, years=2, first_year=0)
        bundle = self._results_bundle_with_programme_timing(setup_controls)
        applied = run_event_ledger_health_economics(
            bundle,
            self._economics_config_with_programme_costs(setup=500000.0),
        )

        app = self._render_health_economics(
            results_bundle=bundle,
            economics_results=applied,
            controls=setup_controls,
            applied_controls=setup_controls,
        )

        chart_rows = self._icer_chart_rows(app)
        scenarios = set(chart_rows["Scenario"])
        self.assertIn("Current analysis - no additional programme overhead", scenarios)
        self.assertIn("Current analysis - user-defined costs", scenarios)
        self.assertIn("SA Health report reference", scenarios)
        self.assertNotIn("Event ledger", self._scenario_table(app).columns)

        no_overhead = self._chart_row(chart_rows, "Current analysis - no additional programme overhead")
        user_defined = self._chart_row(chart_rows, "Current analysis - user-defined costs")
        reference = self._chart_row(chart_rows, "SA Health report reference")
        x_col = "DALYs averted compared with business as usual"
        y_col = "Incremental cost compared with business as usual (AUD)"

        self.assertAlmostEqual(no_overhead[x_col], user_defined[x_col], places=10)
        self.assertAlmostEqual(
            no_overhead["Active TB averted"],
            user_defined["Active TB averted"],
            places=10,
        )
        self.assertAlmostEqual(user_defined[y_col] - no_overhead[y_col], 500000.0, places=6)
        self.assertAlmostEqual(reference[x_col], 16.173794759755, places=10)
        self.assertAlmostEqual(reference[y_col], -92369.6296, places=2)
        self.assertEqual(no_overhead["Legend label"], "Current: no overhead")
        self.assertEqual(user_defined["Legend label"], "Current: user costs")
        self.assertEqual(reference["Legend label"], "Report reference")
        self.assertTrue(
            any(
                "Economic-only changes move the current analysis vertically on the ICER plane" in item.value
                for item in app.caption
            )
        )

        scenarios_table = self._scenario_table(app)
        base_row = scenarios_table.loc[scenarios_table["Scenario"] == "No additional programme overhead entered"].iloc[0]
        standalone_row = scenarios_table.loc[
            scenarios_table["Scenario"] == "Standalone programme - user-defined additional costs"
        ].iloc[0]
        self.assertEqual(base_row["DALYs averted"], standalone_row["DALYs averted"])
        self.assertEqual(base_row["Active TB averted"], standalone_row["Active TB averted"])

    def test_health_economics_icer_plane_reconciles_annual_programme_cost_timing(self) -> None:
        controls = self._programme_controls(setup=0.0, annual=50000.0, years=3, first_year=0)
        bundle = self._results_bundle_with_programme_timing(controls)
        applied = run_event_ledger_health_economics(
            bundle,
            self._economics_config_with_programme_costs(annual=50000.0),
        )

        app = self._render_health_economics(
            results_bundle=bundle,
            economics_results=applied,
            controls=controls,
            applied_controls=controls,
        )
        chart_rows = self._icer_chart_rows(app)
        no_overhead = self._chart_row(chart_rows, "Current analysis - no additional programme overhead")
        user_defined = self._chart_row(chart_rows, "Current analysis - user-defined costs")
        x_col = "DALYs averted compared with business as usual"
        y_col = "Incremental cost compared with business as usual (AUD)"
        expected_added = sum(50000.0 / (1.03 ** year) for year in range(0, 3))

        self.assertAlmostEqual(no_overhead[x_col], user_defined[x_col], places=10)
        self.assertAlmostEqual(
            no_overhead["Active TB averted"],
            user_defined["Active TB averted"],
            places=10,
        )
        self.assertAlmostEqual(user_defined[y_col] - no_overhead[y_col], expected_added, places=6)

    def test_health_economics_icer_plane_reconciles_combined_programme_costs(self) -> None:
        controls = self._programme_controls(setup=100000.0, annual=20000.0, years=3, first_year=1)
        bundle = self._results_bundle_with_programme_timing(controls)
        applied = run_event_ledger_health_economics(
            bundle,
            self._economics_config_with_programme_costs(setup=100000.0, annual=20000.0),
        )

        app = self._render_health_economics(
            results_bundle=bundle,
            economics_results=applied,
            controls=controls,
            applied_controls=controls,
        )
        chart_rows = self._icer_chart_rows(app)
        no_overhead = self._chart_row(chart_rows, "Current analysis - no additional programme overhead")
        user_defined = self._chart_row(chart_rows, "Current analysis - user-defined costs")
        x_col = "DALYs averted compared with business as usual"
        y_col = "Incremental cost compared with business as usual (AUD)"
        expected_added = 100000.0 + sum(20000.0 / (1.03 ** year) for year in range(1, 4))

        self.assertAlmostEqual(no_overhead[x_col], user_defined[x_col], places=10)
        self.assertAlmostEqual(user_defined[y_col] - no_overhead[y_col], expected_added, places=6)
        summary = next(
            item.value for item in app.dataframe
            if "Discounted total additional programme cost" in set(getattr(item.value, "columns", []))
        )
        self.assertEqual(
            summary["Discounted total additional programme cost"].iloc[0],
            f"AUD {expected_added:,.0f}",
        )

    def test_health_economics_icer_plane_collapses_identical_current_points_after_defaults(self) -> None:
        app = self._render_health_economics(
            results_bundle={
                "metadata": {"scenarioLabel": "test_reference"},
                "technical": {"eventLedger": self.package["eventLedger"]},
            },
            economics_results=self.package["economics"],
            controls=self._programme_controls(),
            applied_controls=self._programme_controls(),
        )
        chart_rows = self._icer_chart_rows(app)
        scenarios = list(chart_rows["Scenario"])
        self.assertIn("Current analysis - no additional programme overhead", scenarios)
        self.assertNotIn("Current analysis - user-defined costs", scenarios)

    def test_rendered_health_economics_widgets_recalculate_without_changing_health(self) -> None:
        app = AppTest.from_file(str(ROOT / "pages" / "4_Economics.py"), default_timeout=60)
        app.session_state["config"] = self.package["config"]
        app.session_state["results_bundle"] = {
            "metadata": {"scenarioLabel": "test_reference"},
            "technical": {"eventLedger": self.package["eventLedger"]},
        }
        app.session_state["economics_config"] = self.package["economicsConfig"]
        app.session_state["economics_results"] = self.package["economics"]
        app.session_state["results_stale"] = False

        app.run()
        self.assertFalse(app.exception)
        self.assertNotIn("Include additional programme costs", [item.label for item in app.toggle])
        self.assertIn("One-off programme setup cost (AUD)", [item.label for item in app.number_input])
        self.assertIn("Annual programme running cost (AUD per year)", [item.label for item in app.number_input])
        self.assertIn("First year annual cost is incurred", [item.label for item in app.number_input])
        before = app.session_state["economics_results"]
        before_cost = self._summary_mean(before, "incrementalCost")
        before_dalys = self._summary_mean(before, "dalysAverted")
        before_tb = self._summary_mean(before, "activeTBCasesPrevented")

        next(item for item in app.number_input if item.label == "One-off programme setup cost (AUD)").set_value(100000.0).run()
        next(item for item in app.number_input if item.label == "IGRA screening test per person - Value used by model").set_value(125.0).run()
        next(item for item in app.button if item.label == "Recalculate economics").click().run(timeout=90)

        self.assertFalse(app.exception)
        after = app.session_state["economics_results"]
        after_cost = self._summary_mean(after, "incrementalCost")
        after_dalys = self._summary_mean(after, "dalysAverted")
        after_tb = self._summary_mean(after, "activeTBCasesPrevented")
        self.assertGreater(after_cost, before_cost + 100000.0)
        self.assertAlmostEqual(after_dalys, before_dalys, places=10)
        self.assertAlmostEqual(after_tb, before_tb, places=10)
        self.assertEqual(
            app.session_state["applied_programme_cost_controls"]["standaloneSetupCost"],
            100000.0,
        )
        applied_rows = app.session_state["health_econ_workspace"]["rows"]
        igra = next(row for row in applied_rows if row["assumptionId"] == "cost.test_igra")
        self.assertEqual(igra["sourceCitation"], "User-defined")
        self.assertEqual(float(igra["currentValue"]), 125.0)

        self.assertIn("Restore SA Health economic defaults", [item.label for item in app.button])

    def _render_health_economics(
        self,
        *,
        results_bundle: dict,
        economics_results: dict,
        controls: dict,
        applied_controls: dict,
    ) -> AppTest:
        app = AppTest.from_file(str(ROOT / "pages" / "4_Economics.py"), default_timeout=60)
        app.session_state["config"] = self.package["config"]
        app.session_state["results_bundle"] = results_bundle
        app.session_state["economics_config"] = self.package["economicsConfig"]
        app.session_state["economics_results"] = economics_results
        app.session_state["results_stale"] = False
        app.session_state["health_econ_delivery_scenarios"] = controls
        app.session_state["applied_programme_cost_controls"] = applied_controls
        app.run(timeout=90)
        self.assertFalse(app.exception)
        return app

    @staticmethod
    def _programme_controls(
        *,
        setup: float = 0.0,
        annual: float = 0.0,
        years: int = 2,
        first_year: int = 0,
        travel: float = 0.0,
        staff: float = 0.0,
        share: float = 0.0,
    ) -> dict:
        return {
            "includeAdditionalProgramCosts": any(value > 0 for value in [setup, annual, travel, staff]),
            "standaloneSetupCost": setup,
            "illustrativeSetupCost": 500000.0,
            "standaloneAnnualRunningCost": annual,
            "standaloneRunningYears": years,
            "annualCostFirstYear": first_year,
            "standaloneTravelOutreachCost": travel,
            "standaloneStaffSupportCost": staff,
            "sharedAttributionShare": share,
        }

    def _results_bundle_with_programme_timing(self, controls: dict) -> dict:
        bundle = {
            "metadata": {"scenarioLabel": "test_reference"},
            "technical": {"eventLedger": deepcopy(self.package["eventLedger"])},
        }
        metadata = bundle["technical"]["eventLedger"].setdefault("metadata", {})
        metadata["programRunningDurationYears"] = int(controls["standaloneRunningYears"])
        metadata["programRunningFirstYear"] = int(controls["annualCostFirstYear"])
        return bundle

    def _economics_config_with_programme_costs(
        self,
        *,
        setup: float = 0.0,
        annual: float = 0.0,
        travel_staff: float = 0.0,
    ) -> dict:
        config = deepcopy(self.package["economicsConfig"])
        cost_values = {
            "program_setup": setup,
            "program_running": annual,
            "travel_outreach_staff_support": travel_staff,
        }
        for item in config.get("costItems") or []:
            item_id = item.get("costItemId")
            if item_id in cost_values:
                item["originalCost"] = float(cost_values[item_id])
                item["convertedTargetYearCost"] = float(cost_values[item_id])
                item.setdefault("resourceUse", {})["costBasis"] = (
                    "annual_during_screening_window"
                    if item_id == "program_running"
                    else "per_person_screened"
                    if item_id == "travel_outreach_staff_support"
                    else "total_once_at_program_start"
                )
        return config

    @staticmethod
    def _scenario_table(app: AppTest):
        return next(
            item.value for item in app.dataframe
            if "Scenario" in set(getattr(item.value, "columns", []))
        )

    @staticmethod
    def _icer_chart_rows(app: AppTest):
        x_col = "DALYs averted compared with business as usual"
        y_col = "Incremental cost compared with business as usual (AUD)"
        for item in app:
            if getattr(item, "type", None) != "arrow_vega_lite_chart":
                continue
            if x_col not in str(item.proto.spec) or y_col not in str(item.proto.spec):
                continue
            for dataset in item.proto.datasets:
                try:
                    frame = pa.ipc.open_stream(dataset.data.data).read_pandas()
                except (pa.ArrowInvalid, OSError):
                    continue
                if {"Scenario", x_col, y_col}.issubset(set(frame.columns)):
                    return frame
        raise AssertionError("Rendered ICER chart rows were not found")

    @staticmethod
    def _chart_row(frame, scenario: str):
        matches = frame.loc[frame["Scenario"] == scenario]
        if matches.empty:
            raise AssertionError(f"Missing chart scenario: {scenario}")
        return matches.iloc[0]

    @staticmethod
    def _summary_mean(economics_results: dict, metric: str) -> float:
        for row in economics_results.get("summaryRows") or []:
            if row.get("metric") == metric and row.get("discountProfile") == "primary":
                return float(row.get("mean"))
        raise AssertionError(f"Missing economics metric: {metric}")

    @staticmethod
    def _money_display_to_float(value: object) -> float:
        text = str(value).replace("-AUD", "-").replace("AUD", "").replace(",", "").strip()
        text = text.replace("- ", "-")
        if text.endswith("saving"):
            text = text.removesuffix("saving").strip()
        return float(text)

    def test_economic_changes_do_not_alter_event_counts(self) -> None:
        ledger = self.package["epidemiology"]
        base_econ = self.package["economicsConfig"]
        edited_econ = deepcopy(base_econ)
        for item in edited_econ["costItems"]:
            if item["costItemId"] == "test_igra":
                item["originalCost"] = item["originalCost"] + 10
        base = run_event_ledger_health_economics(ledger, base_econ)
        edited = run_event_ledger_health_economics(ledger, edited_econ)

        self.assertEqual(
            base["replicateResults"].iloc[0]["activeTBCasesPrevented"],
            edited["replicateResults"].iloc[0]["activeTBCasesPrevented"],
        )
        self.assertNotEqual(
            base["replicateResults"].iloc[0]["incrementalCost"],
            edited["replicateResults"].iloc[0]["incrementalCost"],
        )

    def test_exported_and_reloaded_economics_config_reproduces_results(self) -> None:
        out_dir = ROOT / "outputs" / "test_tmp" / "sa_health_reference_package_test"
        result = write_sa_health_reference_package(out_dir, n_reps=10)
        econ_path = out_dir / "economics_config.json"
        loaded_econ = json.loads(econ_path.read_text(encoding="utf-8"))
        rerun = run_event_ledger_health_economics(result["package"]["epidemiology"], loaded_econ)
        original = result["package"]["economics"]["replicateResults"]
        reloaded = rerun["replicateResults"]

        self.assertAlmostEqual(
            original[original["discountProfile"] == "primary"].iloc[0]["incrementalCost"],
            reloaded[reloaded["discountProfile"] == "primary"].iloc[0]["incrementalCost"],
        )
        self.assertTrue((out_dir / "sa_health_apy_working_reference_outputs.xlsx").exists())
        self.assertTrue((out_dir / "technical_recent_remote_scenario.json").exists())
        self.assertTrue((out_dir / "figures" / "screening_treatment_cascade.svg").exists())
        with (out_dir / "tables" / "executive_summary.csv").open("r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        self.assertTrue(any(row["metric"] == "Analysis status" and "provisional" in row["value"] for row in rows))


if __name__ == "__main__":
    unittest.main()
