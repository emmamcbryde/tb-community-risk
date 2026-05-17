from __future__ import annotations

import unittest

from engine.apy.regimen import (
    apply_regimen_overrides,
    default_regimen_library,
    get_regimen_from_library,
    regimen_partial_efficacy,
    validate_regimen,
)


class ApyRegimenTests(unittest.TestCase):
    def test_default_library_contains_all_regimens(self) -> None:
        library = default_regimen_library()

        self.assertEqual(set(library), {"r3HP", "r4R", "r3HR", "r6H", "r9H"})

    def test_3hp_expected_fields_and_validation(self) -> None:
        reg = get_regimen_from_library("3HP", default_regimen_library())

        self.assertEqual(reg["label"], "3HP")
        self.assertEqual(reg["months"], 3)
        self.assertEqual(reg["targetDoses"], 12)
        self.assertEqual(reg["pComplete"], 0.80)
        self.assertEqual(reg["pADRstop"], 0.05)
        self.assertEqual(reg["effFull"], 0.85)
        self.assertTrue(validate_regimen(reg)["isValid"])

    def test_regimen_overrides_work(self) -> None:
        reg = get_regimen_from_library("3HP", default_regimen_library())

        updated = apply_regimen_overrides(
            reg,
            {
                "regimenPComplete": 0.7,
                "regimenADRstop": 0.1,
                "regimenEffFull": 0.8,
            },
        )

        self.assertEqual(updated["pComplete"], 0.7)
        self.assertEqual(updated["pADRstop"], 0.1)
        self.assertEqual(updated["effFull"], 0.8)

    def test_invalid_completion_and_adr_combination_fails(self) -> None:
        reg = get_regimen_from_library("3HP", default_regimen_library())
        reg["pComplete"] = 0.96
        reg["pADRstop"] = 0.05

        report = validate_regimen(reg)

        self.assertFalse(report["isValid"])

    def test_partial_efficacy_none_linear_and_threshold(self) -> None:
        reg = get_regimen_from_library("3HP", default_regimen_library())

        reg["partialMode"] = "none"
        self.assertEqual(regimen_partial_efficacy(reg, 1.0), 0.0)

        reg["partialMode"] = "linear"
        self.assertAlmostEqual(regimen_partial_efficacy(reg, 0.5), 0.425)

        reg["partialMode"] = "threshold80"
        self.assertEqual(regimen_partial_efficacy(reg, 0.8), 0.0)
        self.assertEqual(regimen_partial_efficacy(reg, 11 / 12), 0.85)

    def test_partial_efficacy_knots(self) -> None:
        reg = get_regimen_from_library("6H", default_regimen_library())

        self.assertAlmostEqual(regimen_partial_efficacy(reg, 0.5), 0.30)
        self.assertAlmostEqual(regimen_partial_efficacy(reg, 1.0), 0.65)

    def test_get_regimen_from_library_is_case_insensitive(self) -> None:
        reg = get_regimen_from_library("3hp", default_regimen_library())

        self.assertEqual(reg["label"], "3HP")


if __name__ == "__main__":
    unittest.main()
