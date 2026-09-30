"""Guard the no-transmission scope documentation and its cross-references."""

from __future__ import annotations

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCOPE_DOC = ROOT / "docs" / "no_transmission_model_scope.md"
SCOPE_LINK = "no_transmission_model_scope.md"


def _text(path: Path) -> str:
    lines = (line.lstrip("> ") for line in path.read_text(encoding="utf-8").splitlines())
    return " ".join(" ".join(lines).split())


class NoTransmissionScopeDocsTests(unittest.TestCase):
    def test_scope_document_states_required_scope(self) -> None:
        text = _text(SCOPE_DOC)
        for phrase in (
            "screening, diagnostic, preventive-treatment and health-economic strategies",
            "does not materially affect the comparison",
            "direct outcomes among the modelled population",
            "herd",
            "secondary infections",
            "force of infection",
            "should not be selected on incidence alone",
            "Starsim",
            "requires scientific review",
            "provisional",
        ):
            with self.subTest(phrase=phrase):
                self.assertIn(phrase.lower(), text.lower())

    def test_scope_document_lists_prohibited_outputs(self) -> None:
        text = _text(SCOPE_DOC).lower()
        self.assertIn("outputs this model must never describe as estimated", text)
        for output in ("secondary infections prevented", "population-level or community-level transmission reduction"):
            with self.subTest(output=output):
                self.assertIn(output, text)

    def test_no_incidence_value_presented_as_validated_cut_off(self) -> None:
        text = _text(SCOPE_DOC).lower()
        self.assertIn("no incidence value is presented here as a validated cut-off", text)
        self.assertIsNone(re.search(r"validated (universal )?(cut-off|threshold) of \d", text))

    def test_readme_and_general_app_docs_cross_reference_scope(self) -> None:
        readme = _text(ROOT / "README.md")
        self.assertIn("Scope and intended use", readme)
        self.assertIn(f"docs/{SCOPE_LINK}", readme)
        for name in ("general_app_milestone1.md", "general_app_milestone2.md", "dynamic_model_readiness_spec.md"):
            with self.subTest(doc=name):
                self.assertIn(SCOPE_LINK, _text(ROOT / "docs" / name))

    def test_starsim_is_not_a_dependency(self) -> None:
        requirements = (ROOT / "requirements.txt").read_text(encoding="utf-8").lower()
        self.assertNotIn("starsim", requirements)


if __name__ == "__main__":
    unittest.main()
