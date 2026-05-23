from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class PythonApyMigrationCheckTest(unittest.TestCase):
    def test_check_script_reports_success(self) -> None:
        result = subprocess.run(
            [sys.executable, "scripts/check_python_apy_migration.py"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        output = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, output)
        self.assertIn("PASSED 2/2 checks", output)


if __name__ == "__main__":
    unittest.main()
