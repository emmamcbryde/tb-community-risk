from __future__ import annotations

from pathlib import Path
import re
import socket
import unittest
from unittest.mock import patch

from streamlit.testing.v1 import AppTest

from app.general.terminology import (
    DETERMINISTIC_LABEL,
    RESTORE_DEFAULTS_LABEL,
    STOCHASTIC_LABEL,
    USER_DEFINED_MARK,
    WORKFLOW_PAGES,
    forbidden_matches,
)
from engine.profiles.demonstration import build_demonstration_profile


ROOT = Path(__file__).resolve().parents[1]
PAGES = [ROOT / path for path, _ in WORKFLOW_PAGES]
GENERAL_SOURCES = [
    ROOT / "general_app.py",
    *PAGES,
    *sorted((ROOT / "app" / "general").glob("*.py")),
    *sorted((ROOT / "engine" / "profiles").glob("*.py")),
    *[p for p in sorted((ROOT / "engine" / "who_incidence").glob("*.py")) if p.name != "adapters.py"],
]
ELEMENT_KINDS = (
    "title", "header", "subheader", "markdown", "caption", "info", "warning", "success", "error", "text",
    "button", "selectbox", "radio", "number_input", "slider", "checkbox", "expander", "code",
)


def _no_network(*args, **kwargs):
    raise AssertionError(f"Network access attempted during page rendering: {args!r}")


def _visible_text(app: AppTest) -> list[str]:
    texts: list[str] = []
    for kind in ELEMENT_KINDS:
        for element in getattr(app, kind, []):
            for attr in ("value", "label", "body", "help"):
                value = getattr(element, attr, None)
                if isinstance(value, str):
                    texts.append(value)
            options = getattr(element, "options", None)
            if options:
                texts.extend(str(option) for option in options)
    for frame in getattr(app, "dataframe", []):
        value = frame.value
        texts.extend(str(column) for column in value.columns)
        texts.extend(str(cell) for cell in value.astype(str).to_numpy().ravel())
    return texts


def _render(path: Path, session: dict | None = None) -> AppTest:
    app = AppTest.from_file(str(path), default_timeout=300)
    for key, value in (session or {}).items():
        app.session_state[key] = value
    with patch.object(socket.socket, "connect", _no_network), patch("socket.create_connection", _no_network):
        app.run(timeout=300)
    return app


class GeneralSourceGuardTests(unittest.TestCase):
    def test_general_sources_have_no_setting_specific_user_text_or_network(self) -> None:
        literal = re.compile(r"(['\"])(?:(?!\1).)*\1")
        for path in GENERAL_SOURCES:
            source = path.read_text(encoding="utf-8")
            with self.subTest(file=path.name):
                self.assertNotRegex(source, r"import (requests|urllib|http\.client|socket)\b")
                self.assertNotIn("gtbreport2025", source)
                if path in PAGES or path.name == "general_app.py":
                    for match in literal.finditer(source):
                        self.assertEqual(forbidden_matches(match.group(0)), [], match.group(0))


class GeneralPageRenderingTests(unittest.TestCase):
    """Render each ordinary page offline and inspect the visible text."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.rendered = {path.name: _render(path) for path in PAGES}

    def test_pages_render_without_exceptions_or_forbidden_terms(self) -> None:
        for name, app in self.rendered.items():
            with self.subTest(page=name):
                self.assertFalse(app.exception, [e.value for e in app.exception])
                found = [m for text in _visible_text(app) for m in forbidden_matches(text)]
                self.assertEqual(found, [])

    def test_population_defaults_to_10000_with_demonstration_warning(self) -> None:
        app = self.rendered["1_Set_up_population.py"]
        size = next(item for item in app.number_input if item.label == "Population size to simulate")
        self.assertEqual(size.value, 10_000)
        warnings = " ".join(item.value for item in app.warning)
        self.assertIn("demonstration working defaults, not evidence for any particular country", warnings)
        self.assertIn(RESTORE_DEFAULTS_LABEL, [button.label for button in app.button])
        self.assertEqual([button.label for button in app.button].count(RESTORE_DEFAULTS_LABEL), 1)

    def test_stochastic_default_is_1000_and_deterministic_preview_available(self) -> None:
        app = self.rendered["3_Run_analysis.py"]
        radio = app.radio[0]
        self.assertEqual(radio.value, "agent_based")
        self.assertIn(STOCHASTIC_LABEL, radio.options)
        self.assertIn(DETERMINISTIC_LABEL, radio.options)
        sims = next(item for item in app.number_input if item.label == "Number of simulated populations")
        self.assertEqual(sims.value, 1_000)
        radio.set_value("expected_value").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["general_analysis"]["analysisMethod"], "expected_value")
        self.assertNotIn("Number of simulated populations", [item.label for item in app.number_input])

    def test_country_selector_offline_and_overrides_marked_then_restored(self) -> None:
        app = _render(ROOT / "general_pages" / "1_Set_up_population.py")
        selector = app.selectbox(key="general_country_select")
        country = next(option for option in selector.options if "(ZAF)" in option)
        with patch.object(socket.socket, "connect", _no_network):
            selector.set_value(country).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["general_profile"]["location"]["iso3"], "ZAF")
        captions = " ".join(item.value for item in app.caption)
        self.assertIn("WHO report year 2025", captions)
        self.assertIn("not the complete WHO dataset", captions)

        app.number_input(key="general_population_input").set_value(5000).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["general_profile"]["populationSize"]["provenance"], "user_defined")
        visible = " ".join(_visible_text(app))
        self.assertIn(USER_DEFINED_MARK, visible)

        next(b for b in app.button if b.label == RESTORE_DEFAULTS_LABEL).click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["general_profile"], build_demonstration_profile().to_dict())
        self.assertEqual(app.number_input(key="general_population_input").value, 10_000)
        self.assertNotIn(USER_DEFINED_MARK, " ".join(_visible_text(app)))


class GeneralEndToEndTests(unittest.TestCase):
    """One small offline run of a profile without risk factors, then results pages."""

    @classmethod
    def setUpClass(cls) -> None:
        from adapters.python_apy_backend import PythonApyBackend
        from engine.profiles.demonstration import demonstration_profile_without_risk_factors
        from engine.profiles.engine_mapping import build_engine_config
        from engine.profiles.population_profile import with_population_size

        profile = with_population_size(demonstration_profile_without_risk_factors(), 200)
        cls.profile = profile
        cls.config = build_engine_config(profile, analysis={"analysisMethod": "agent_based", "nReps": 2, "seed": 1})
        backend = PythonApyBackend(ROOT)
        with patch.object(socket.socket, "connect", _no_network):
            cls.bundle = backend.run_scenario_bundle(cls.config, validation_report=backend.validate_config(cls.config))

    def test_profile_without_risk_factors_runs(self) -> None:
        self.assertTrue(self.bundle["headline"]["keyMetricsRows"])
        link = self.bundle["technical"]["interfaceConfig"]["generalProfileLink"]
        self.assertEqual(link["profileHash"], self.profile.profile_hash())

    def test_results_and_economics_pages_render_without_forbidden_terms(self) -> None:
        session = {
            "general_profile": self.profile.to_dict(),
            "general_results_bundle": self.bundle,
            "general_results_config": self.config,
            "general_results_stale": False,
        }
        results = _render(ROOT / "general_pages" / "4_Results.py", session)
        self.assertFalse(results.exception)
        self.assertEqual([m for t in _visible_text(results) for m in forbidden_matches(t)], [])

        economics = _render(ROOT / "general_pages" / "5_Health_economics.py", session)
        with patch.object(socket.socket, "connect", _no_network):
            next(b for b in economics.button if b.label == "Calculate health economics").click().run(timeout=300)
        self.assertFalse(economics.exception)
        self.assertIsNotNone(economics.session_state["general_economics_results"])
        visible = _visible_text(economics)
        self.assertIn("DALYs averted", visible)
        self.assertEqual([m for t in visible for m in forbidden_matches(t)], [])


if __name__ == "__main__":
    unittest.main()
