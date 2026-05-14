from pathlib import Path
import unittest


class ApySourceFormattingTests(unittest.TestCase):
    def test_apy_python_sources_are_not_minified(self):
        root = Path(__file__).resolve().parents[1]
        source_dir = root / "engine" / "apy"

        for path in sorted(source_dir.glob("*.py")):
            if path.name == "__init__.py":
                continue

            with self.subTest(file=path.name):
                lines = path.read_text(encoding="utf-8").splitlines()

                self.assertGreater(
                    len(lines),
                    10,
                    f"{path} appears collapsed into too few lines",
                )

                longest_line = max((len(line) for line in lines), default=0)
                self.assertLessEqual(
                    longest_line,
                    240,
                    f"{path} contains an unexpectedly long line",
                )


if __name__ == "__main__":
    unittest.main()
