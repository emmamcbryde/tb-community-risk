from __future__ import annotations

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import check_environment  # noqa: E402


class EnvironmentCheckTests(unittest.TestCase):
    def test_reports_pins_and_supported_range(self) -> None:
        pins = check_environment.pinned_requirements()
        self.assertEqual(pins["streamlit"], "1.39.0")
        errors, warnings, info = check_environment.check()
        self.assertIn("streamlit", info)
        self.assertEqual(errors, [], errors)
        self.assertIsInstance(warnings, list)

    def test_matlab_not_required(self) -> None:
        self.assertNotIn("matlab", " ".join(check_environment.pinned_requirements()).lower())


if __name__ == "__main__":
    unittest.main()
