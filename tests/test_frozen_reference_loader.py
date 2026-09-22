from __future__ import annotations

from copy import deepcopy
from unittest.mock import patch
import unittest

import pandas as pd

from app.parameter_workspace import unified_default_session_state
from app.state import ensure_frozen_reference_loaded_if_eligible
from app.health_economics_inputs import (
    apply_assumptions_to_economics_config,
    reconcile_workspace_state,
)
from engine.apy.frozen_reference import (
    FROZEN_REFERENCE_PRIMARY_REPLICATES_PATH,
    FROZEN_REFERENCE_REPS,
    is_frozen_sa_health_reference_eligible,
    load_frozen_reference_results,
    recalculate_frozen_reference_economics,
    validate_stochastic_replicates,
)


class _FakeStreamlit:
    def __init__(self) -> None:
        self.session_state: dict = {}


class FrozenReferenceLoaderTests(unittest.TestCase):
    def test_default_report_configuration_is_eligible_and_changed_inputs_are_not(self) -> None:
        state = unified_default_session_state()
        config = state["config"]

        self.assertTrue(is_frozen_sa_health_reference_eligible(config))

        changed = deepcopy(config)
        changed["screenCoverage"] = 0.4
        self.assertFalse(is_frozen_sa_health_reference_eligible(changed))

        deterministic = deepcopy(config)
        deterministic["analysisMethod"] = "expected_value"
        self.assertFalse(is_frozen_sa_health_reference_eligible(deterministic))

    def test_frozen_artifact_contains_two_thousand_complete_paired_replicates(self) -> None:
        payload = load_frozen_reference_results()
        bundle = payload["resultsBundle"]
        reps = payload["referenceEconomics"]["replicateResults"]
        artifact = pd.read_csv(FROZEN_REFERENCE_PRIMARY_REPLICATES_PATH)

        self.assertEqual(len(artifact), FROZEN_REFERENCE_REPS)
        self.assertEqual(artifact["replicateId"].nunique(), FROZEN_REFERENCE_REPS)
        self.assertTrue(
            {
                "replicateId",
                "dalysAverted",
                "incrementalCost",
                "activeTBCasesPrevented",
                "comparatorActiveTBCases",
                "interventionActiveTBCases",
                "configurationHash",
                "economicConfigurationHash",
                "analysisBasis",
                "naturalHistorySemantics",
                "seed",
                "nReps",
                "replicateContractVersion",
            }.issubset(artifact.columns)
        )
        validation = validate_stochastic_replicates(
            reps,
            expected_reps=FROZEN_REFERENCE_REPS,
            configuration_hash=bundle["metadata"]["configurationHash"],
        )
        self.assertTrue(validation["isValid"], validation["errors"])
        self.assertAlmostEqual(float(reps["dalysAverted"].mean()), 16.173794759755, places=10)
        self.assertAlmostEqual(float(reps["incrementalCost"].mean()), -92369.6295674358, places=6)
        self.assertAlmostEqual(float(reps["activeTBCasesPrevented"].mean()), 11.9665, places=6)

    def test_frozen_cost_recalculation_preserves_pairing_and_vertical_setup_shift(self) -> None:
        payload = load_frozen_reference_results()
        base = recalculate_frozen_reference_economics(
            payload["resultsBundle"],
            payload["economicsConfig"],
        )
        edited_config = deepcopy(payload["economicsConfig"])
        setup = next(item for item in edited_config["costItems"] if item.get("costItemId") == "program_setup")
        setup.update(
            {
                "originalCost": 500000.0,
                "convertedTargetYearCost": 500000.0,
                "originalCurrency": "AUD",
                "targetCurrency": "AUD",
                "originalPriceYear": "2019",
                "targetPriceYear": "2019",
            }
        )
        setup.setdefault("resourceUse", {})["costBasis"] = "total_once_at_program_start"
        edited = recalculate_frozen_reference_economics(payload["resultsBundle"], edited_config)

        before = base["replicateResults"].sort_values("replicateId").reset_index(drop=True)
        after = edited["replicateResults"].sort_values("replicateId").reset_index(drop=True)
        self.assertEqual(list(before["replicateId"]), list(after["replicateId"]))
        self.assertTrue(before["dalysAverted"].round(12).equals(after["dalysAverted"].round(12)))
        self.assertTrue(
            before["activeTBCasesPrevented"].round(12).equals(after["activeTBCasesPrevented"].round(12))
        )
        deltas = after["incrementalCost"] - before["incrementalCost"]
        self.assertTrue(((deltas - 500000.0).abs() < 1e-6).all())

    def test_workspace_setup_cost_edit_preserves_reference_cost_components(self) -> None:
        payload = load_frozen_reference_results()
        base = recalculate_frozen_reference_economics(
            payload["resultsBundle"],
            payload["economicsConfig"],
        )
        workspace = reconcile_workspace_state(
            None,
            payload["economicsConfig"],
            registry=payload["economicsConfig"].get("assumptionEvidenceRegistry"),
        )
        edited_config = apply_assumptions_to_economics_config(
            payload["economicsConfig"],
            workspace["rows"],
        )
        setup = next(item for item in edited_config["costItems"] if item.get("costItemId") == "program_setup")
        setup.update(
            {
                "originalCost": 500000.0,
                "originalCurrency": "AUD",
                "originalPriceYear": "2019",
                "targetCurrency": "AUD",
                "targetPriceYear": "2019",
            }
        )
        setup.setdefault("resourceUse", {})["costBasis"] = "total_once_at_program_start"
        edited = recalculate_frozen_reference_economics(payload["resultsBundle"], edited_config)

        self.assertIsNotNone(edited)
        unchanged_components = [
            "testingCost",
            "treatmentCost",
            "adrManagementCost",
            "returnForResultsCost",
            "clinicalReviewCost",
            "activeTBExclusionWorkupCost",
            "baselineTBDiseaseCost",
            "interventionTBDiseaseCost",
        ]
        for component in unchanged_components:
            self.assertAlmostEqual(
                float(edited["costs"][component]),
                float(base["costs"][component]),
                places=6,
                msg=component,
            )
        self.assertAlmostEqual(float(base["costs"]["programSetupCost"]), 0.0, places=6)
        self.assertAlmostEqual(float(edited["costs"]["programSetupCost"]), 500000.0, places=6)
        base_mean = float(
            next(
                row["mean"]
                for row in base["summaryRows"]
                if row.get("metric") == "incrementalCost" and row.get("discountProfile") == "primary"
            )
        )
        edited_mean = float(
            next(
                row["mean"]
                for row in edited["summaryRows"]
                if row.get("metric") == "incrementalCost" and row.get("discountProfile") == "primary"
            )
        )
        self.assertAlmostEqual(base_mean, -92369.62956743593, places=6)
        self.assertAlmostEqual(edited_mean, 407630.37043256406, places=6)

    def test_default_state_loader_loads_frozen_reference_without_runner_state(self) -> None:
        fake_st = _FakeStreamlit()
        with patch("app.state.st", fake_st):
            loaded = ensure_frozen_reference_loaded_if_eligible(force=True)

        self.assertTrue(loaded)
        bundle = fake_st.session_state["results_bundle"]
        economics = fake_st.session_state["economics_results"]
        self.assertEqual(bundle["metadata"]["nReps"], FROZEN_REFERENCE_REPS)
        self.assertEqual(len(economics["replicateResults"]), FROZEN_REFERENCE_REPS)
        self.assertTrue(fake_st.session_state["validation_report"]["loadedFrozenReference"])
        self.assertFalse(fake_st.session_state["results_stale"])
        self.assertFalse(fake_st.session_state["dirty_economics"])

    def test_direct_state_loader_provides_cloud_ready_economics(self) -> None:
        fake_st = _FakeStreamlit()
        with patch("app.state.st", fake_st):
            ensure_frozen_reference_loaded_if_eligible(force=True)

        economics = fake_st.session_state["economics_results"]
        reps = economics["replicateResults"]
        self.assertEqual(len(reps), FROZEN_REFERENCE_REPS)
        self.assertAlmostEqual(float(reps["dalysAverted"].mean()), 16.173794759755, places=10)
        self.assertAlmostEqual(float(reps["incrementalCost"].mean()), -92369.6295674358, places=6)
        self.assertAlmostEqual(float(reps["activeTBCasesPrevented"].mean()), 11.9665, places=6)
        self.assertIn("summaryRows", economics)
        summary_metrics = {row["metric"] for row in economics["summaryRows"]}
        self.assertTrue({"dalysAverted", "incrementalCost", "activeTBCasesPrevented"}.issubset(summary_metrics))


if __name__ == "__main__":
    unittest.main()
