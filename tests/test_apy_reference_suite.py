from __future__ import annotations

from pathlib import Path
import unittest

from engine.apy.reference_loader import (
    list_reference_scenario_dirs,
    load_reference_dir,
    load_reference_scenario_suite,
)


REFERENCE_ROOT = Path("validation") / "matlab_reference"
SUITE_FILE = REFERENCE_ROOT / "scenario_suite_v1.json"
DEFAULT_FIXTURE = REFERENCE_ROOT / "default_random_igra_3hp_N1500_seed1"


class ApyReferenceSuiteTests(unittest.TestCase):
    def test_scenario_suite_json_loads(self) -> None:
        suite = load_reference_scenario_suite(SUITE_FILE)

        self.assertEqual(len(suite), 8)
        self.assertEqual(
            suite[0]["scenario_id"],
            "default_random_igra_3hp_N1500_seed1",
        )
        self.assertIn("config_overrides", suite[0])

    def test_existing_default_fixture_still_loads(self) -> None:
        loaded = load_reference_dir(DEFAULT_FIXTURE)

        self.assertIn("summary", loaded)
        self.assertFalse(loaded["summary"].empty)

    def test_reference_scenario_dirs_only_lists_committed_fixtures(self) -> None:
        dirs = list_reference_scenario_dirs(REFERENCE_ROOT)
        names = {path.name for path in dirs}

        self.assertIn("default_random_igra_3hp_N1500_seed1", names)
        self.assertNotIn("scenario_suite_v1.json", names)

    def test_missing_fixture_directory_reports_clear_error(self) -> None:
        with self.assertRaisesRegex(FileNotFoundError, "MATLAB reference directory"):
            load_reference_dir(REFERENCE_ROOT / "missing_fixture")


if __name__ == "__main__":
    unittest.main()
