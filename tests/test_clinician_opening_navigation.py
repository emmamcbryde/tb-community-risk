from __future__ import annotations

import ast
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
IMPLEMENTATION_TERMS = ("MATLAB", "Python", "backend", "v9", "port", "ABM")


def _contains_implementation_term(text: str) -> str | None:
    for term in IMPLEMENTATION_TERMS:
        if re.search(rf"\b{re.escape(term)}\b", text, flags=re.IGNORECASE):
            return term
    return None


def _streamlit_call_strings(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    strings: list[str] = []
    visible_functions = {
        "title",
        "subheader",
        "markdown",
        "write",
        "info",
        "warning",
        "success",
        "error",
        "button",
        "download_button",
        "form_submit_button",
        "text_input",
        "number_input",
        "checkbox",
        "radio",
        "selectbox",
        "page_link",
        "caption",
    }
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (
            isinstance(func, ast.Attribute)
            and isinstance(func.value, ast.Name)
            and func.value.id == "st"
            and func.attr in visible_functions
        ):
            continue
        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                strings.append(arg.value)
        for keyword in node.keywords:
            if (
                keyword.arg in {"label", "page_title"}
                and isinstance(keyword.value, ast.Constant)
                and isinstance(keyword.value.value, str)
            ):
                strings.append(keyword.value.value)
    return strings


class ClinicianOpeningNavigationTests(unittest.TestCase):
    def test_opening_page_visible_text_avoids_implementation_terms(self) -> None:
        path = ROOT / "pages" / "0_Start.py"
        visible_text = "\n".join(_streamlit_call_strings(path))
        offending = _contains_implementation_term(visible_text)
        self.assertIsNone(offending, visible_text)
        self.assertIn("LTBI Screening Decision Tool", visible_text)
        self.assertIn("Use default parameters", visible_text)
        self.assertIn("Review or change parameters", visible_text)
        self.assertIn("Run with these defaults", visible_text)
        self.assertIn("Continue to current results", visible_text)

    def test_standard_navigation_labels_avoid_implementation_terms(self) -> None:
        text = (ROOT / "streamlit_app.py").read_text(encoding="utf-8")
        tree = ast.parse(text)
        standard_labels: list[str] = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Dict):
                continue
            for key, value in zip(node.keys, node.values):
                if not (
                    isinstance(key, ast.Constant)
                    and key.value == "LTBI Screening Tool"
                    and isinstance(value, ast.List)
                ):
                    continue
                standard_labels.append(key.value)
                for page_call in value.elts:
                    if not isinstance(page_call, ast.Call):
                        continue
                    for keyword in page_call.keywords:
                        if (
                            keyword.arg == "title"
                            and isinstance(keyword.value, ast.Constant)
                            and isinstance(keyword.value.value, str)
                        ):
                            standard_labels.append(keyword.value.value)
        self.assertEqual(
            standard_labels,
            [
                "LTBI Screening Tool",
                "Start",
                "Define Strategy",
                "Run Analysis",
                "Results",
                "Health Economics",
                "Explore Decisions",
                "Evidence & Assumptions",
            ],
        )
        visible_nav = "\n".join(standard_labels)
        offending = _contains_implementation_term(visible_nav)
        self.assertIsNone(offending, visible_nav)

    def test_standard_workflow_visible_text_avoids_implementation_terms(self) -> None:
        standard_pages = [
            ROOT / "pages" / "0_Start.py",
            ROOT / "pages" / "1_Scenario.py",
            ROOT / "pages" / "2_Run_Model.py",
            ROOT / "pages" / "3_Results.py",
            ROOT / "pages" / "4_Economics.py",
            ROOT / "pages" / "5_Decision_Analysis.py",
            ROOT / "pages" / "6_Evidence_Assumptions.py",
        ]
        visible_text = "\n".join(
            text
            for path in standard_pages
            for text in _streamlit_call_strings(path)
        )
        offending = _contains_implementation_term(visible_text)
        self.assertIsNone(offending, visible_text)

    def test_health_economics_and_evidence_pages_are_separate(self) -> None:
        text = (ROOT / "streamlit_app.py").read_text(encoding="utf-8")

        self.assertIn('st.Page("pages/4_Economics.py", title="Health Economics")', text)
        self.assertIn(
            'st.Page("pages/6_Evidence_Assumptions.py", title="Evidence & Assumptions")',
            text,
        )
        self.assertLess(text.index('title="Results"'), text.index('title="Health Economics"'))
        self.assertLess(text.index('title="Health Economics"'), text.index('title="Explore Decisions"'))

    def test_results_links_to_health_economics(self) -> None:
        text = (ROOT / "pages" / "3_Results.py").read_text(encoding="utf-8")

        self.assertIn("Open Health Economics", text)
        self.assertNotIn("Open Evidence & Assumptions", text)

    def test_health_economics_page_uses_event_ledger_and_runs_existing_engine(self) -> None:
        text = (ROOT / "pages" / "4_Economics.py").read_text(encoding="utf-8")

        self.assertIn('st.title("Health Economics")', text)
        self.assertIn("eventLedger", text)
        self.assertIn("backend.run_economics(results_bundle, econ_config)", text)
        self.assertIn("Run the screening analysis before recalculating health economics", text)

    def test_evidence_assumptions_page_is_readiness_focused(self) -> None:
        text = (ROOT / "pages" / "6_Evidence_Assumptions.py").read_text(encoding="utf-8")

        self.assertIn('st.title("Evidence & Assumptions")', text)
        self.assertIn("assess_apy_reference_readiness", text)
        self.assertIn("load_apy_evidence_registry", text)

    def test_standard_pages_do_not_use_newer_streamlit_width_strings(self) -> None:
        """Keep deployed Streamlit 1.39-compatible dataframe/chart calls."""
        offenders: list[tuple[str, str, str]] = []
        standard_pages = [
            ROOT / "streamlit_app.py",
            ROOT / "pages" / "0_Start.py",
            ROOT / "pages" / "1_Scenario.py",
            ROOT / "pages" / "2_Run_Model.py",
            ROOT / "pages" / "3_Results.py",
            ROOT / "pages" / "4_Economics.py",
            ROOT / "pages" / "5_Decision_Analysis.py",
            ROOT / "pages" / "6_Evidence_Assumptions.py",
        ]
        for path in standard_pages:
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                func = node.func
                if not (
                    isinstance(func, ast.Attribute)
                    and isinstance(func.value, ast.Name)
                    and func.value.id == "st"
                    and func.attr in {"dataframe", "altair_chart"}
                ):
                    continue
                for keyword in node.keywords:
                    if (
                        keyword.arg == "width"
                        and isinstance(keyword.value, ast.Constant)
                        and keyword.value.value in {"stretch", "content"}
                    ):
                        offenders.append(
                            (
                                path.relative_to(ROOT).as_posix(),
                                func.attr,
                                str(keyword.value.value),
                            )
                        )
        self.assertEqual(offenders, [])

    def test_every_streamlit_form_contains_submit_button(self) -> None:
        offenders: list[str] = []
        for path in [*ROOT.glob("pages/**/*.py"), ROOT / "streamlit_app.py"]:
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if not isinstance(node, ast.With):
                    continue
                is_form = any(
                    isinstance(item.context_expr, ast.Call)
                    and isinstance(item.context_expr.func, ast.Attribute)
                    and isinstance(item.context_expr.func.value, ast.Name)
                    and item.context_expr.func.value.id == "st"
                    and item.context_expr.func.attr == "form"
                    for item in node.items
                )
                if not is_form:
                    continue
                has_submit = any(
                    isinstance(child, ast.Call)
                    and isinstance(child.func, ast.Attribute)
                    and isinstance(child.func.value, ast.Name)
                    and child.func.value.id == "st"
                    and child.func.attr == "form_submit_button"
                    for child in ast.walk(node)
                )
                if not has_submit:
                    offenders.append(path.relative_to(ROOT).as_posix())
        self.assertEqual(offenders, [])


if __name__ == "__main__":
    unittest.main()
