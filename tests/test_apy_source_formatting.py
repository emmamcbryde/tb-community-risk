from pathlib import Path
import unittest


def _cr_only_count(raw: bytes) -> int:
    return sum(
        1
        for index, byte in enumerate(raw)
        if byte == 13 and (index + 1 >= len(raw) or raw[index + 1] != 10)
    )


class ApySourceFormattingTests(unittest.TestCase):
    def test_apy_python_sources_are_not_minified(self):
        root = Path(__file__).resolve().parents[1]
        source_dir = root / "engine" / "apy"

        for path in sorted(source_dir.glob("*.py")):
            if path.name == "__init__.py":
                continue

            with self.subTest(file=path.name):
                raw = path.read_bytes()
                lines = raw.decode("utf-8").splitlines()

                self.assertGreater(
                    raw.count(b"\n"),
                    0,
                    f"{path} contains no LF line endings",
                )
                self.assertEqual(
                    _cr_only_count(raw),
                    0,
                    f"{path} contains CR-only line endings",
                )

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

    def test_apy_distributional_validation_test_is_not_minified(self):
        root = Path(__file__).resolve().parents[1]
        path = root / "tests" / "test_apy_distributional_validation.py"

        self._assert_multiline_lf_source(path)

    def test_gitattributes_uses_detectable_line_endings(self):
        root = Path(__file__).resolve().parents[1]
        path = root / ".gitattributes"
        raw = path.read_bytes()

        self.assertGreater(raw.count(b"\n"), 0)
        self.assertEqual(_cr_only_count(raw), 0)

    def _assert_multiline_lf_source(self, path: Path):
        raw = path.read_bytes()
        lines = raw.decode("utf-8").splitlines()

        self.assertGreater(
            raw.count(b"\n"),
            0,
            f"{path} contains no LF line endings",
        )
        self.assertEqual(
            _cr_only_count(raw),
            0,
            f"{path} contains CR-only line endings",
        )
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
