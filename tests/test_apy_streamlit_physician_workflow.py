from __future__ import annotations

from io import BytesIO
import json
import os
import sys
import types
import unittest
from pathlib import Path

from openpyxl import load_workbook

from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend
from app.epidemiology_inputs import (
    apply_epidemiology_updates,
    fraction_to_percent,
    percent_to_fraction,
    risk_override_from_percentages,
)
from app.results_workbook import build_results_workbook
from engine.apy.calibration import calibrate_from_config
from engine.apy.config import build_default_config
from engine.apy.runner import run_scenario_with_do_nothing


class PhysicianWorkflowHelperTests(unittest.TestCase):
    def test_percentage_fraction_conversion(self) -> None:
        self.assertAlmostEqual(percent_to_fraction(7.5), 0.075)
        self.assertAlmostEqual(fraction_to_percent(0.0125), 1.25)
        with self.assertRaises(ValueError):
            percent_to_fraction(0)
        with self.assertRaises(ValueError):
            percent_to_fraction(100)

    def test_default_and_custom_prevalence_updates(self) -> None:
        cfg = build_default_config()
        custom = apply_epidemiology_updates(
            cfg,
            use_default_ltbi_prevalence=False,
            ltbi_prevalence_percent=1.0,
            use_default_active_tb_prevalence=False,
            active_tb_prevalence_percent=0.2,
        )
        self.assertEqual(custom["ltbiPrevalence"], 0.01)
        self.assertEqual(custom["activeTBPrevalence"], 0.002)
        restored = apply_epidemiology_updates(
            custom,
            use_default_ltbi_prevalence=True,
            ltbi_prevalence_percent=9.0,
            use_default_active_tb_prevalence=True,
            active_tb_prevalence_percent=1.0,
        )
        self.assertIsNone(restored["ltbiPrevalence"])
        self.assertIsNone(restored["activeTBPrevalence"])

    def test_risk_factor_scalar_and_age_group_inputs(self) -> None:
        self.assertEqual(risk_override_from_percentages("Single overall", [20]), 0.2)
        self.assertEqual(
            risk_override_from_percentages("Three age groups", [1, 2, 3]),
            [0.01, 0.02, 0.03],
        )
        with self.assertRaises(ValueError):
            risk_override_from_percentages("Three age groups", [1, 2])


class PythonApyPrevalencePathTests(unittest.TestCase):
    def test_custom_ltbi_prevalence_reaches_calibration(self) -> None:
        cfg = build_default_config()
        cfg["ltbiPrevalence"] = 0.01
        cfg["activeTBPrevalence"] = (10 / 770) * 0.01 / (47 / 624)
        calibration = calibrate_from_config(cfg)
        self.assertAlmostEqual(calibration["targetInfPrev"], 0.01)
        self.assertAlmostEqual(calibration["expectedInfPrev"], 0.01, places=4)

    def test_low_prevalence_changes_python_yield_without_matlab(self) -> None:
        self.assertNotIn("matlab.engine", sys.modules)
        default_cfg = build_default_config()
        default_cfg.update({"N": 300, "nReps": 20, "seed": 4, "screenCoverage": 1.0})
        low_cfg = dict(default_cfg)
        low_cfg["ltbiPrevalence"] = 0.01
        low_cfg["activeTBPrevalence"] = (10 / 770) * 0.01 / (47 / 624)
        default_out = run_scenario_with_do_nothing(default_cfg)
        low_out = run_scenario_with_do_nothing(low_cfg)
        default_inf = _median(default_out, "nInfected")
        low_inf = _median(low_out, "nInfected")
        default_positive = _median(default_out, "nTestPositiveNonActive")
        low_positive = _median(low_out, "nTestPositiveNonActive")
        self.assertGreater(default_inf, low_inf)
        self.assertGreater(default_positive, low_positive)
        self.assertNotIn("matlab.engine", sys.modules)

    def test_python_backend_save_load_round_trip_preserves_custom_prevalence(self) -> None:
        backend = PythonApyBackend(repo_root())
        cfg = backend.default_config()
        cfg["ltbiPrevalence"] = 0.02
        cfg["activeTBPrevalence"] = 0.003
        cfg["riskPrev"]["contact"] = [0.01, 0.02, 0.03]
        tmp = repo_root() / ".tmp_test_scenario_roundtrip"
        tmp.mkdir(exist_ok=True)
        path = tmp / "scenario.json"
        try:
            backend.save_scenario(cfg, str(path))
            loaded, _, info = backend.load_scenario(str(path))
        finally:
            if path.exists():
                path.unlink()
            if tmp.exists():
                tmp.rmdir()
        self.assertEqual(loaded["ltbiPrevalence"], 0.02)
        self.assertEqual(loaded["activeTBPrevalence"], 0.003)
        self.assertEqual(loaded["riskPrev"]["contact"], [0.01, 0.02, 0.03])
        self.assertEqual(info["backend"], "python_apy")


class WorkbookExportTests(unittest.TestCase):
    def test_workbook_contains_required_sheets_and_values(self) -> None:
        cfg = build_default_config()
        cfg.update({"N": 80, "nReps": 3, "seed": 3, "ltbiPrevalence": 0.02})
        bundle = PythonApyBackend(repo_root()).run_scenario_bundle(cfg)
        payload = build_results_workbook(
            config=cfg,
            bundle=bundle,
            backend_status={"name": "python_apy", "experimental": True},
            economics_results=None,
            economics_config=None,
            results_stale=False,
        )
        wb = load_workbook(BytesIO(payload), read_only=True, data_only=True)
        expected = {
            "Read_me",
            "Scenario_inputs",
            "Risk_factor_inputs",
            "Headline_results",
            "Summary_results",
            "Natural_history",
            "Technical_metadata",
            "Economics",
        }
        self.assertTrue(expected.issubset(set(wb.sheetnames)))
        scenario_values = list(wb["Scenario_inputs"].iter_rows(values_only=True))
        self.assertIn(("LTBI prevalence source", "Custom", None), scenario_values)
        self.assertIn(("LTBI prevalence input", 0.02, "percentage"), scenario_values)
        economics_values = list(wb["Economics"].iter_rows(values_only=True))
        self.assertIn(("Economics not run", None, "No zero values have been substituted for missing economics outputs."), economics_values)
        wb.close()

    def test_stale_state_function_marks_results_and_economics(self) -> None:
        fake_streamlit = types.SimpleNamespace(
            session_state={
                "results_bundle": {"metadata": {}},
                "economics_results": {"summaryRows": []},
                "compare_baseline_bundle": {"metadata": {}},
            },
            cache_resource=lambda **_: (lambda fn: fn),
        )
        old_streamlit = sys.modules.get("streamlit")
        sys.modules["streamlit"] = fake_streamlit
        try:
            import importlib
            import app.state as state

            importlib.reload(state)
            state.mark_config_changed()
            self.assertTrue(fake_streamlit.session_state["dirty_config"])
            self.assertTrue(fake_streamlit.session_state["results_stale"])
            self.assertTrue(fake_streamlit.session_state["dirty_economics"])
            self.assertTrue(fake_streamlit.session_state["compare_results_stale"])
        finally:
            if old_streamlit is not None:
                sys.modules["streamlit"] = old_streamlit
            else:
                sys.modules.pop("streamlit", None)
            import importlib
            import app.state as state

            importlib.reload(state)


class HostedDeploymentSafetyTests(unittest.TestCase):
    def test_new_session_defaults_to_python_and_blocks_matlab_without_opt_in(self) -> None:
        fake_streamlit = types.SimpleNamespace(
            session_state={},
            cache_resource=lambda **_: (lambda fn: fn),
        )
        old_streamlit = sys.modules.get("streamlit")
        old_env = os.environ.pop("APY_ENABLE_MATLAB_BACKEND", None)
        sys.modules["streamlit"] = fake_streamlit
        try:
            import importlib
            import app.state as state

            importlib.reload(state)
            state.init_session_state()
            self.assertEqual(state.get_backend_name(), "python_apy")
            self.assertFalse(state.matlab_backend_enabled())
            with self.assertRaises(ValueError):
                state.set_backend_name("matlab")
        finally:
            if old_env is not None:
                os.environ["APY_ENABLE_MATLAB_BACKEND"] = old_env
            if old_streamlit is not None:
                sys.modules["streamlit"] = old_streamlit
            else:
                sys.modules.pop("streamlit", None)
            import importlib
            import app.state as state

            importlib.reload(state)

    def test_python_hosted_workflow_completes_without_matlab(self) -> None:
        sys.modules.pop("matlab.engine", None)
        backend = PythonApyBackend(repo_root())
        config = backend.default_config()
        config.update({"N": 100, "nReps": 5, "seed": 9})
        config = apply_epidemiology_updates(
            config,
            use_default_ltbi_prevalence=False,
            ltbi_prevalence_percent=2.0,
            use_default_active_tb_prevalence=False,
            active_tb_prevalence_percent=(10 / 770) * 0.02 / (47 / 624) * 100,
        )
        report = backend.validate_config(config)
        self.assertTrue(report["isValid"])
        bundle = backend.run_scenario_bundle(config, validation_report=report)
        payload = build_results_workbook(
            config=config,
            bundle=bundle,
            backend_status=backend.status(),
            economics_results=None,
            economics_config=None,
            results_stale=False,
        )
        wb = load_workbook(BytesIO(payload), read_only=True, data_only=True)
        self.assertIn("Scenario_inputs", wb.sheetnames)
        wb.close()
        self.assertNotIn("matlab.engine", sys.modules)


def _median(model_out: dict, metric: str) -> float:
    summary = model_out["results"]["summary"]
    return float(summary.loc[summary["Metric"] == metric, "Median"].iloc[0])


if __name__ == "__main__":
    unittest.main()
