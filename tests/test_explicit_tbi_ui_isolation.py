from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import json
from pathlib import Path
import re
import subprocess
import sys
import unittest
from unittest.mock import patch

from streamlit.testing.v1 import AppTest

from app.explicit_tbi_formatting import display_table, format_aud, format_number, format_percent
from app.explicit_tbi_session import (
    DISPATCH_COUNT_KEY, PAGE_KEY, RESULT_KEY, SCHEMA, SCHEMA_KEY, STATE_KEY,
    initialize_explicit_session,
)
from app.explicit_tbi_workflow import (
    WORKFLOW_PAGES, build_explicit_workflow_config, default_workflow_state,
    discard_incompatible_result, health_signature, results_to_csv,
    run_explicit_workflow_analysis, technical_results_json,
)
from engine.apy.explicit_recent_remote_runner import EXPLICIT_RECENT_REMOTE_PATHWAY_ID

ROOT = Path(__file__).resolve().parents[1]
ENTRY = ROOT / "explicit_tbi_ui/app.py"
START = "136de0556d0ce6ce0ce9cb76b04e22efa5e2cffc"
FROZEN = "03cc16e4a52a55e10019dab16bd4f599571c4b90"


@lru_cache(maxsize=1)
def real_result():
    return run_explicit_workflow_analysis(default_workflow_state())


def app_for(state=None, result=None, page="Set up"):
    app = AppTest.from_file(str(ENTRY), default_timeout=180)
    app.session_state[SCHEMA_KEY] = deepcopy(SCHEMA)
    app.session_state[STATE_KEY] = deepcopy(state or default_workflow_state())
    app.session_state[RESULT_KEY] = deepcopy(result)
    app.session_state[PAGE_KEY] = page
    return app


def workflow(app):
    return next(radio for radio in app.radio if radio.label == "Workflow")


def button(app, label):
    return next(item for item in app.button if item.label == label)


def rendered_text(app):
    return "\n".join(
        str(item.value)
        for kind in ("title", "subheader", "caption", "info", "warning", "error", "markdown")
        for item in getattr(app, kind)
    )


class IsolatedNavigationTests(unittest.TestCase):
    def test_structural_isolation_and_discovery(self):
        self.assertNotEqual(ENTRY.parent, ROOT)
        self.assertFalse((ENTRY.parent / "pages").exists())
        self.assertTrue((ROOT / "pages/4_Economics.py").exists())
        self.assertFalse((ROOT / "explicit_tbi_app.py").exists())
        # Exercise Streamlit's actual discovery implementation in a fresh process.
        code = (
            "import streamlit.source_util as s; "
            "print([p['script_path'] for p in s.get_pages('explicit_tbi_ui/app.py').values()])"
        )
        # 1.39 moved discovery into PagesManager; the live protocol audit covers it.
        import streamlit.source_util as source_util
        if hasattr(source_util, "get_pages"):
            paths = subprocess.check_output(["python", "-c", code], cwd=ROOT, text=True)
            self.assertNotIn("4_Economics", paths)
            self.assertNotIn("1_Scenario", paths)

    def test_exactly_four_rendered_pages(self):
        app = app_for().run()
        self.assertFalse(app.exception)
        self.assertEqual(workflow(app).options, list(WORKFLOW_PAGES))
        self.assertEqual(workflow(app).value, "Set up")

    def test_no_scenario_rendered(self):
        app = app_for().run()
        self.assertNotIn("Scenario", workflow(app).options)

    def test_no_economics_rendered(self):
        app = app_for().run()
        self.assertFalse(any("Economic" in value for value in workflow(app).options))
        self.assertNotIn("Health Economics", rendered_text(app))

    def test_no_dynamic_or_comparison_pages(self):
        app = app_for().run()
        for label in ("Dynamic Model", "Compare", "Dynamic ABM Compare", "Technical Settings", "Decision Analysis"):
            self.assertNotIn(label, workflow(app).options)

    def test_direct_legacy_economics_rejects_explicit_session(self):
        for state_key in (SCHEMA_KEY, "explicit_tbi_results"):
            app = AppTest.from_file(str(ROOT / "pages/4_Economics.py"), default_timeout=180)
            app.session_state[state_key] = deepcopy(SCHEMA)
            app.session_state["economics_results"] = {"stale_cost": 123456789}
            app.run()
            self.assertFalse(app.exception)
            self.assertTrue(app.error)
            self.assertIn("no validated costable event ledger", app.error[0].value)
            self.assertNotIn("economics_results", app.session_state)
            self.assertFalse(app.dataframe)
            self.assertFalse(app.metric)
        # Also reject a direct attempt to put explicit results in a legacy slot.
        app = AppTest.from_file(str(ROOT / "pages/4_Economics.py"), default_timeout=180)
        app.session_state["results_bundle"] = {"analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID}
        app.session_state["economics_results"] = {"stale_cost": 123456789}
        app.run()
        self.assertTrue(app.error)
        self.assertNotIn("results_bundle", app.session_state)
        self.assertNotIn("economics_results", app.session_state)
        self.assertFalse(app.dataframe)

    def test_contaminated_session_opens_clean_setup(self):
        app = app_for(page="Results")
        for key in ("results_bundle", "economics_results", "dynamic_results", "decision_analysis", "explicit_tbi_results"):
            app.session_state[key] = {"old": 123456789}
        app.session_state[STATE_KEY]["recentTBITarget"] = 0.6
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(workflow(app).value, "Set up")
        self.assertIsNone(app.session_state[RESULT_KEY])
        self.assertEqual(app.session_state[STATE_KEY], default_workflow_state())
        self.assertNotIn("economics_results", app.session_state)
        self.assertIn("Incompatible previous session state was discarded", rendered_text(app))

    def test_unmarked_or_incompatible_namespace_rejected(self):
        for marker in ({"version": 0}, {"version": 1, "analysisPathway": "matlab_v9"}):
            session = {SCHEMA_KEY: marker, STATE_KEY: default_workflow_state(), PAGE_KEY: "Results"}
            self.assertTrue(initialize_explicit_session(session))
            self.assertEqual(session[PAGE_KEY], "Set up")
            self.assertIsNone(session[RESULT_KEY])
        self.assertIsNone(discard_incompatible_result({"analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID}))


class PrecisionAndDispatchTests(unittest.TestCase):
    def test_population_input_groups_thousands_and_rejects_fractions(self):
        app = app_for().run()
        population = next(item for item in app.text_input if item.label == "Total modelled population (people)")
        self.assertEqual(population.value, "1,500")
        population.set_value("2,345").run()
        self.assertEqual(app.session_state[STATE_KEY]["population"], 2345)
        next(item for item in app.text_input if item.label == "Total modelled population (people)").set_value("2.5").run()
        self.assertIn("whole number of people", rendered_text(app))
        workflow(app).set_value("Run analysis").run()
        self.assertNotIn("Run explicit TBI analysis", [item.label for item in app.button])
        workflow(app).set_value("Set up").run()
        population = next(item for item in app.text_input if item.label == "Total modelled population (people)")
        self.assertEqual(population.value, "2.5")
        population.set_value("2345").run()
        self.assertEqual(next(item for item in app.text_input if item.label == "Total modelled population (people)").value, "2,345")
        self.assertFalse(app.error)

    def test_percentage_units_and_formats_rendered(self):
        state = default_workflow_state()
        state["baselineActiveTBProportion"] = 0.0001
        app = app_for(state).run()
        self.assertIn("Input known TBI prevalence", rendered_text(app))
        self.assertIn("% of total population", rendered_text(app))
        labels = {item.label: item for item in app.number_input}
        recent = labels["Recently infected within the past 5 years"]
        remote = next(item for item in app.number_input if item.label.startswith("Remote infection only"))
        self.assertEqual(recent.value, 2.0)
        self.assertEqual(remote.value, 8.0)
        self.assertEqual(recent.proto.format, "%.1f")
        self.assertEqual(remote.proto.format, "%.1f")
        baseline = labels["Baseline active TB (% of total population)"]
        self.assertEqual(baseline.proto.format, "%.2f")
        self.assertEqual(baseline.value, 0.01)
        self.assertNotIn("0.000000", rendered_text(app))
        for label in ("Screening coverage (%)", "Treatment initiation (%)", "Treatment completion (%)"):
            self.assertEqual(labels[label].proto.format, "%.1f")

    def test_stored_proportions_preserve_full_precision_on_rerun(self):
        state = default_workflow_state()
        state["recentTBITarget"] = 0.023456789123456
        state["screeningCoverage"] = 0.30123456789123
        app = app_for(state).run()
        for _ in range(2):
            app.run()
            self.assertEqual(app.session_state[STATE_KEY]["recentTBITarget"], state["recentTBITarget"])
            self.assertEqual(app.session_state[STATE_KEY]["screeningCoverage"], state["screeningCoverage"])
        next(item for item in app.number_input if item.label.startswith("Recently infected")).set_value(2.7).run()
        self.assertEqual(app.session_state[STATE_KEY]["recentTBITarget"], 2.7 / 100.0)
        workflow(app).set_value("Results").run()
        workflow(app).set_value("Set up").run()
        self.assertEqual(app.session_state[STATE_KEY]["recentTBITarget"], 2.7 / 100.0)

    def test_formatting_never_alters_calculations_hashes_or_exports(self):
        result = deepcopy(real_result())
        original = deepcopy(result)
        state = default_workflow_state()
        signature = health_signature(state)
        config = build_explicit_workflow_config(state)
        display_table(result["summaryRows"])
        format_percent(state["recentTBITarget"])
        self.assertEqual(result, original)
        self.assertEqual(signature, health_signature(state))
        self.assertEqual(config, build_explicit_workflow_config(state))

    def test_all_ordinary_result_tables_have_at_most_two_decimals(self):
        app = app_for(result=real_result(), page="Results").run()
        self.assertFalse(app.exception)
        self.assertIn("Approximately", rendered_text(app))
        self.assertNotIn("Fitted recent exposure hazard", rendered_text(app))
        self.assertNotIn(real_result()["configHash"], rendered_text(app))
        self.assertIn("Expected value", app.dataframe[0].value.columns)
        self.assertIn("N/A", app.dataframe[0].value.to_string())
        self.assertIn("1,500", app.dataframe[0].value.to_string())
        for frame in app.dataframe:
            for value in frame.value.to_numpy().flat:
                self.assertIsInstance(value, str)
                self.assertFalse(re.search(r"\d+\.\d{3,}", value), value)

    def test_technical_json_csv_preserve_full_precision(self):
        result = deepcopy(real_result())
        value = 0.277932123456789
        result["summaryRows"][0]["Expected value"] = value
        self.assertEqual(display_table(result["summaryRows"]).iloc[0]["Expected value"], "0.28")
        self.assertIn(repr(value), results_to_csv(result))
        self.assertEqual(json.loads(technical_results_json(result))["summaryRows"][0]["Expected value"], value)
        self.assertEqual(result["summaryRows"][0]["Expected value"], value)

    def test_future_currency_contract_is_presentation_only(self):
        for value, expected in ((95052.4, "AUD 95,052"), (-92369.6, "-AUD 92,370"), (0.1, "AUD 0"), (-0.1, "AUD 0"), (1.5, "AUD 2"), (-1.5, "-AUD 2")):
            original = value
            self.assertEqual(format_aud(value), expected)
            self.assertEqual(value, original)
        self.assertEqual(format_aud(1234.56, per_unit="DALY"), "AUD 1,235 per DALY")
        self.assertEqual(format_number(450.0), "450")
        self.assertEqual(format_number(42.619847), "42.62")
        self.assertEqual(format_percent(0.02), "2.0%")

    def test_invalid_prevalence_totals_block_dispatch(self):
        app = app_for().run()
        next(item for item in app.number_input if item.label.startswith("Recently infected")).set_value(95.0).run()
        next(item for item in app.number_input if item.label.startswith("Remote infection only")).set_value(10.0).run()
        self.assertIn("must not exceed 100%", rendered_text(app))
        workflow(app).set_value("Run analysis").run()
        self.assertNotIn("Run explicit TBI analysis", [item.label for item in app.button])
        self.assertEqual(app.session_state[DISPATCH_COUNT_KEY], 0)

    def test_one_click_dispatches_once_and_navigation_survives_reruns(self):
        app = app_for(page="Run analysis").run()
        with patch("app.explicit_tbi_workflow.run_explicit_workflow_analysis", return_value=deepcopy(real_result())) as dispatch:
            button(app, "Run explicit TBI analysis").click().run()
            self.assertFalse(app.exception)
            self.assertEqual(dispatch.call_count, 1)
            self.assertEqual(app.session_state[DISPATCH_COUNT_KEY], 1)
            app.run()
            workflow(app).set_value("Results").run()
            app.run()
            workflow(app).set_value("Run analysis").run()
            self.assertEqual(dispatch.call_count, 1)
            self.assertEqual(app.session_state[DISPATCH_COUNT_KEY], 1)

    def test_non_costable_marker_is_retained(self):
        result = real_result()
        self.assertFalse(result["eventLedger"]["available"])
        self.assertFalse(result["eventLedger"]["validation"]["isValid"])
        self.assertEqual(result["eventLedger"]["annualEvents"], [])
        self.assertNotIn("economics", result)


class FrozenBoundaryTests(unittest.TestCase):
    def test_frozen_navigation_sources_and_artifacts_unchanged(self):
        # The economics page adds only a pre-import guard for explicit sessions.
        paths = ["streamlit_app.py", "abm", "engine", "app/reference_data", "reports/sa_health_apy_health_economic_working_report_v1", "validation/matlab_reference",
                 ":(exclude)engine/apy/explicit_recent_remote_runner.py",
                 ":(exclude)engine/apy/explicit_recent_remote_progression.py",
                 ":(exclude)engine/apy/population_tbi_progression_risk.py"]
        changed = subprocess.check_output(["git", "diff", START, "--name-only", "--", *paths], cwd=ROOT, text=True)
        self.assertEqual(changed.strip(), "")
        for ref in ("release/sa-health-apy-he-v1.0.0", "sa-health-apy-he-v1.0.0^{}"):
            actual = subprocess.check_output(["git", "rev-parse", ref], cwd=ROOT, text=True).strip()
            self.assertEqual(actual, FROZEN)

    def test_frozen_entry_renders_original_navigation_and_results(self):
        import streamlit as st
        if not hasattr(st, "navigation"):
            self.skipTest("Frozen entry requires repository-pinned Streamlit 1.39; covered by the 1.39 validation run")
        # 1.39 AppTest emits an empty tree for this st.navigation router. Use
        # actual live-render messages, including the nine page registry and
        # populated frozen Results tables; do not mistake an empty tree for QA.
        completed = subprocess.run(
            [sys.executable, str(ROOT / "scripts/validate_explicit_tbi_ui.py"),
             "--streamlit-target", str(Path(st.__file__).resolve().parents[1]),
             "--frozen", "--port", "8539", "--output-dir", "abm/output/milestone3g/frozen-integrity"],
            cwd=ROOT, text=True, capture_output=True, timeout=180,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn("PASS:", completed.stdout)


if __name__ == "__main__":
    unittest.main()
