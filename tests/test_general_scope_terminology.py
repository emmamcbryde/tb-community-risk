"""Guard the no-transmission model boundary in general-application wording and exports."""

from __future__ import annotations

from pathlib import Path
import re
import unittest

from app.general.provenance_export import LIMITATIONS
from app.general.terminology import (
    DIRECT_ACTIVE_TB_AVERTED,
    DIRECT_EFFECTS_CAPTION,
    GENERAL_METRIC_LABELS,
    METRIC_DEFINITIONS,
    relabel_outcome_rows,
)
from app.results_page_display import FRIENDLY_METRIC_LABELS
from engine import model_scope
from engine.model_scope import prohibited_claims


ROOT = Path(__file__).resolve().parents[1]
ORDINARY_SOURCES = [
    ROOT / "general_app.py",
    *sorted((ROOT / "general_pages").glob("*.py")),
    *sorted((ROOT / "app" / "general").glob("*.py")),
    ROOT / "engine" / "profiles" / "country.py",
    ROOT / "engine" / "profiles" / "local_incidence.py",
    ROOT / "engine" / "who_incidence" / "trend.py",
]
STRING_LITERAL = re.compile(r'"(?:[^"\\\n]|\\.)*"|\'(?:[^\'\\\n]|\\.)*\'')


def _literals(path: Path) -> str:
    source = path.read_text(encoding="utf-8")
    return " ".join(match.group(0)[1:-1] for match in STRING_LITERAL.finditer(source))


class ScopeTerminologyTests(unittest.TestCase):
    def test_prohibited_patterns_detect_known_claims(self) -> None:
        for claim in (
            "Secondary infections prevented",
            "reduces community incidence",
            "herd protection",
            "reduces the force of infection",
            "progress towards elimination",
            "outbreaks prevented",
            "Transmission-mediated benefits are not yet included.",
            "They are not yet used to infer infection pressure or transmission.",
        ):
            with self.subTest(claim=claim):
                self.assertTrue(prohibited_claims(claim))
        for allowed in (DIRECT_EFFECTS_CAPTION, model_scope.DIRECT_EFFECTS_STATEMENT, model_scope.MODEL_IDENTITY):
            with self.subTest(allowed=allowed):
                self.assertEqual(prohibited_claims(allowed), [])

    def test_ordinary_sources_make_no_transmission_effect_claims(self) -> None:
        for path in ORDINARY_SOURCES:
            with self.subTest(file=path.relative_to(ROOT).as_posix()):
                text = _literals(path)
                self.assertEqual(prohibited_claims(text), [])
                self.assertNotIn("screened and treated", text)
                self.assertNotIn("Transmission-mediated benefits", text)

    def test_canonical_statements_are_reused(self) -> None:
        from engine.profiles.country import INCIDENCE_LINK_NOTE
        from engine.who_incidence.trend import DESCRIPTIVE_STATEMENT, INCIDENCE_TO_INFECTION_POLICY

        self.assertEqual(INCIDENCE_LINK_NOTE, model_scope.INCIDENCE_DESCRIPTIVE_STATEMENT)
        self.assertEqual(DESCRIPTIVE_STATEMENT, model_scope.INCIDENCE_DESCRIPTIVE_STATEMENT)
        self.assertEqual(INCIDENCE_TO_INFECTION_POLICY, model_scope.INCIDENCE_TO_INFECTION_POLICY)
        self.assertNotIn("not yet", model_scope.INCIDENCE_DESCRIPTIVE_STATEMENT)

    def test_limitations_export_keeps_technical_limitations_visible(self) -> None:
        text = " ".join(LIMITATIONS.split())
        self.assertEqual(prohibited_claims(text), [])
        for phrase in (
            model_scope.MODEL_IDENTITY,
            "direct outcomes among the modelled population",
            "Secondary infections or cases prevented, changes in force of infection and changes in community incidence are not estimated",
            "`dynamicComparison`",
            "contain no transmission feedback",
            "new infection and reinfection during follow-up are not modelled",
            "Stochastic simulation intervals",
            "progression-hazard multipliers",
        ):
            with self.subTest(phrase=phrase):
                self.assertIn(phrase, text)

    def test_one_metric_has_one_meaning(self) -> None:
        by_shared: dict[str, set[str]] = {}
        for metric, shared in FRIENDLY_METRIC_LABELS.items():
            if metric in GENERAL_METRIC_LABELS:
                by_shared.setdefault(shared, set()).add(GENERAL_METRIC_LABELS[metric])
        self.assertTrue(all(len(labels) == 1 for labels in by_shared.values()), by_shared)
        for metric in ("cumulative_cases_averted", "nPreventedActiveTB", "activeTBCasesPrevented"):
            self.assertEqual(GENERAL_METRIC_LABELS[metric], DIRECT_ACTIVE_TB_AVERTED)
        for label in METRIC_DEFINITIONS:
            self.assertIn(label, GENERAL_METRIC_LABELS.values())
        self.assertTrue(all(len(label) <= 60 for label in GENERAL_METRIC_LABELS.values()))

    def test_relabel_keeps_unmapped_outcomes(self) -> None:
        rows = [{"Outcome": "People screened", "Median": 1.0}, {"Outcome": "Active TB cases averted", "Median": 2.0}]
        self.assertEqual(
            relabel_outcome_rows(rows, FRIENDLY_METRIC_LABELS),
            [{"Outcome": "People screened", "Median": 1.0}, {"Outcome": DIRECT_ACTIVE_TB_AVERTED, "Median": 2.0}],
        )

    def test_shared_release_labels_are_unchanged(self) -> None:
        self.assertEqual(FRIENDLY_METRIC_LABELS["nPreventedActiveTB"], "Active TB cases averted")
        self.assertEqual(FRIENDLY_METRIC_LABELS["relative_reduction_cumulative_active_tb_cases"], "Relative reduction in active TB")


if __name__ == "__main__":
    unittest.main()
