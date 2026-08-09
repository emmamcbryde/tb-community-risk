from __future__ import annotations

from copy import deepcopy
from io import BytesIO
import unittest

from openpyxl import load_workbook

from app.health_economics_inputs import (
    assumptions_csv,
    assumptions_workbook,
    apply_assumptions_to_economics_config,
    editable_assumption_rows,
    parse_assumptions_csv,
    validate_editable_assumptions,
)
from app.results_workbook import build_results_workbook
from engine.apy.config import build_default_config
from engine.apy.economics import build_default_economics_config
from engine.apy.evidence import load_apy_evidence_registry
from tests.test_apy_event_ledger_economics import _synthetic_econ
from tests.test_apy_evidence_registry import _minimal_reviewed_registry, _reviewed_epi_config


class HealthEconomicsInputsWorkspaceTests(unittest.TestCase):
    def test_registry_rows_load_into_editable_assumptions_table(self) -> None:
        rows = editable_assumption_rows(economics_config=build_default_economics_config())
        ids = {row["assumptionId"] for row in rows}

        self.assertIn("cost.test_igra", ids)
        self.assertIn("daly.active_tb_disability_weight", ids)
        self.assertIn("threshold.gdp_per_capita", ids)
        cost_row = next(row for row in rows if row["assumptionId"] == "cost.test_igra")
        self.assertEqual(cost_row["costBasis"], "per_person_screened")

    def test_editing_cost_updates_working_copy_and_authoritative_cost_item(self) -> None:
        econ_config = _synthetic_econ()
        original_registry = load_apy_evidence_registry()
        rows = editable_assumption_rows(_minimal_reviewed_registry(), econ_config)
        row = next(item for item in rows if item["assumptionId"] == "cost.test_igra")
        row.update(
            {
                "currentValue": "42.5",
                "originalCurrency": "AUD",
                "originalPriceYear": "2025-26",
                "targetCurrency": "AUD",
                "targetPriceYear": "2025-26",
                "sourceCitation": "Synthetic edited cost source",
                "reviewStatus": "configured_reviewed",
                "provisional": False,
            }
        )

        updated = apply_assumptions_to_economics_config(econ_config, rows, config=_reviewed_epi_config())
        item = next(item for item in updated["costItems"] if item["costItemId"] == "test_igra")

        self.assertEqual(item["originalCost"], 42.5)
        self.assertEqual(item["sourceCitation"], "Synthetic edited cost source")
        self.assertEqual(updated["costs"]["test"]["IGRA"], 42.5)
        self.assertTrue(updated["assumptionEvidenceValidation"]["epidemiologyReady"])
        self.assertEqual(load_apy_evidence_registry(), original_registry)

    def test_incomplete_cost_assumptions_block_icer(self) -> None:
        econ_config = _synthetic_econ()
        rows = editable_assumption_rows(_minimal_reviewed_registry(), econ_config)
        row = next(item for item in rows if item["assumptionId"] == "cost.test_igra")
        row["sourceCitation"] = ""
        row["provisional"] = False

        report = validate_editable_assumptions(rows, econ_config, config=_reviewed_epi_config())

        self.assertFalse(report["summary"]["costReady"])
        self.assertFalse(report["summary"]["icerReady"])
        self.assertIn("cost.test_igra", report["summary"]["remainingBlockers"])

    def test_cost_complete_but_daly_incomplete_allows_cost_consequence_only(self) -> None:
        econ_config = _synthetic_econ()
        rows = editable_assumption_rows(_minimal_reviewed_registry(), econ_config)
        daly_row = next(item for item in rows if item["assumptionId"] == "daly.active_tb_disability_weight")
        daly_row["sourceCitation"] = ""
        daly_row["reviewStatus"] = "unresolved"
        daly_row["provisional"] = True

        report = validate_editable_assumptions(rows, econ_config, config=_reviewed_epi_config())

        self.assertTrue(report["summary"]["costConsequenceReady"])
        self.assertFalse(report["summary"]["dalyReady"])
        self.assertFalse(report["summary"]["icerReady"])

    def test_reviewed_exclusions_require_rationale(self) -> None:
        rows = editable_assumption_rows(_minimal_reviewed_registry(), _synthetic_econ())
        row = next(item for item in rows if item["assumptionId"] == "daly.tpt_health_loss")
        row["inclusionStatus"] = "excluded"
        row["reviewStatus"] = "reviewed_exclusion"
        row["notes"] = ""

        report = validate_editable_assumptions(rows, _synthetic_econ(), config=_reviewed_epi_config())

        self.assertIn("daly.tpt_health_loss", report["rowMessages"])
        self.assertFalse(report["summary"]["dalyReady"])

    def test_typed_numerical_value_without_source_does_not_become_ready(self) -> None:
        rows = editable_assumption_rows(_minimal_reviewed_registry(), _synthetic_econ())
        row = next(item for item in rows if item["assumptionId"] == "daly.active_tb_duration")
        row["currentValue"] = "0.5"
        row["sourceCitation"] = ""
        row["reviewStatus"] = "configured_reviewed"
        row["provisional"] = False

        report = validate_editable_assumptions(rows, _synthetic_econ(), config=_reviewed_epi_config())

        self.assertFalse(report["summary"]["dalyReady"])
        self.assertIn("DALY source citation is missing.", report["rowMessages"]["daly.active_tb_duration"])

    def test_exported_assumptions_preserve_review_and_provisional_status(self) -> None:
        rows = editable_assumption_rows(_minimal_reviewed_registry(), _synthetic_econ())
        row = next(item for item in rows if item["assumptionId"] == "cost.test_igra")
        row["reviewStatus"] = "configured_reviewed"
        row["provisional"] = False

        csv_payload = assumptions_csv(rows)
        workbook_payload = assumptions_workbook(rows, validate_editable_assumptions(rows, _synthetic_econ(), config=_reviewed_epi_config()))
        wb = load_workbook(BytesIO(workbook_payload), read_only=True)
        ws = wb["Edited_assumptions"]
        headers = [cell.value for cell in next(ws.iter_rows(min_row=1, max_row=1))]
        status_idx = headers.index("reviewStatus")
        provisional_idx = headers.index("provisional")
        exported = {
            row_values[0].value: (row_values[status_idx].value, row_values[provisional_idx].value)
            for row_values in ws.iter_rows(min_row=2)
        }
        wb.close()

        self.assertIn("configured_reviewed", csv_payload)
        self.assertEqual(exported["cost.test_igra"], ("configured_reviewed", "false"))

    def test_uploaded_assumptions_validate_before_application(self) -> None:
        rows = editable_assumption_rows(_minimal_reviewed_registry(), _synthetic_econ())
        row = next(item for item in rows if item["assumptionId"] == "cost.test_igra")
        row["costBasis"] = "wrong_basis"

        uploaded = parse_assumptions_csv(assumptions_csv(rows))
        report = validate_editable_assumptions(uploaded, _synthetic_econ(), config=_reviewed_epi_config())

        self.assertFalse(report["isValidForApplication"])
        self.assertIn("cost.test_igra", report["rowMessages"])

    def test_missing_economic_value_is_exported_blank_not_zero(self) -> None:
        payload = build_results_workbook(
            config=build_default_config(),
            bundle={"metadata": {}, "technical": {}, "headline": {}, "downloads": {}},
            economics_config=build_default_economics_config(),
        )
        wb = load_workbook(BytesIO(payload), read_only=True)
        ws = wb["Edited_assumptions"]
        rows = list(ws.iter_rows(values_only=True))
        wb.close()

        self.assertEqual(rows[0][0], "Status")
        self.assertEqual(rows[1][0], "No edited assumptions applied")
        self.assertNotEqual(rows[1][0], 0)

    def test_workbook_includes_original_and_applied_assumption_sheets(self) -> None:
        econ_config = apply_assumptions_to_economics_config(
            _synthetic_econ(),
            editable_assumption_rows(_minimal_reviewed_registry(), _synthetic_econ()),
            config=_reviewed_epi_config(),
        )

        payload = build_results_workbook(
            config=_reviewed_epi_config(),
            bundle={"metadata": {}, "technical": {}, "headline": {}, "downloads": {}},
            economics_config=econ_config,
        )
        wb = load_workbook(BytesIO(payload), read_only=True)
        names = set(wb.sheetnames)
        wb.close()

        self.assertIn("Edited_assumptions", names)
        self.assertIn("Applied_session_assumptions", names)
        self.assertIn("Original_evidence_registry", names)
        self.assertIn("Assumption_validation", names)


if __name__ == "__main__":
    unittest.main()
