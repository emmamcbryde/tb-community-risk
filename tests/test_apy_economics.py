from __future__ import annotations

import json
import unittest

from engine.apy.economics import calculate_economics


class ApyEconomicsTests(unittest.TestCase):
    def test_calculate_economics_returns_stable_placeholder_payload(self) -> None:
        payload = calculate_economics(
            {"metadata": {"backend": "python"}, "results": {"summary": []}},
            {"metadata": {}, "costs": {}},
        )

        self.assertEqual(
            list(payload.keys()),
            ["available", "message", "metadata", "results"],
        )
        self.assertIs(payload["available"], False)
        self.assertIn("not yet ported", payload["message"])
        self.assertEqual(payload["results"], {})
        json.dumps(payload)

    def test_result_bundle_must_be_mapping(self) -> None:
        with self.assertRaisesRegex(TypeError, "result_bundle must be a mapping"):
            calculate_economics([], {})

    def test_economics_config_must_be_mapping(self) -> None:
        with self.assertRaisesRegex(TypeError, "economics_config must be a mapping"):
            calculate_economics({"metadata": {}, "results": {}}, [])

    def test_result_bundle_requires_metadata_and_results(self) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "missing required top-level key\\(s\\): results",
        ):
            calculate_economics({"metadata": {}}, {})

    def test_result_bundle_sections_must_be_mappings(self) -> None:
        with self.assertRaisesRegex(ValueError, "result_bundle.results must be a mapping"):
            calculate_economics({"metadata": {}, "results": []}, {})


if __name__ == "__main__":
    unittest.main()
