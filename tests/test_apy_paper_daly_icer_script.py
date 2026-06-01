from __future__ import annotations

import shutil
from pathlib import Path
import subprocess
import sys
import unittest


class ApyPaperDalyIcerScriptTests(unittest.TestCase):
    def test_quick_mode_writes_expected_outputs(self) -> None:
        output_dir = Path("outputs") / "test_apy_paper_daly_icer_quick"
        if output_dir.exists():
            shutil.rmtree(output_dir)
        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "scripts/run_apy_paper_daly_icer.py",
                    "--quick",
                    "--output-dir",
                    str(output_dir),
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=180,
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)
            expected_files = {
                "apy_daly_icer_strategy_summary.csv",
                "apy_daly_icer_health_outcomes.csv",
                "apy_daly_icer_costs.csv",
                "apy_daly_icer_sensitivity.csv",
                "apy_daly_icer_notes.md",
            }
            self.assertEqual(
                {path.name for path in output_dir.iterdir()},
                expected_files,
            )
            strategy_text = (output_dir / "apy_daly_icer_strategy_summary.csv").read_text(
                encoding="utf-8"
            )
            self.assertIn("costPerDALYAverted", strategy_text)
            self.assertIn("costPerQALYGained", strategy_text)
            notes = (output_dir / "apy_daly_icer_notes.md").read_text(encoding="utf-8")
            self.assertIn("direct prevented TB cases", notes)
        finally:
            if output_dir.exists():
                shutil.rmtree(output_dir)


if __name__ == "__main__":
    unittest.main()
