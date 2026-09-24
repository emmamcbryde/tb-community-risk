"""The reference-calibration memoisation must not change any calibration output."""

from __future__ import annotations

import json
import unittest
from unittest.mock import patch

from adapters.serialization import to_json_like
from engine.apy import calibration_policy as policy
from engine.apy.calibration import calibrate_from_config
from engine.apy.working_defaults import build_unified_working_default_preset


def canonical(value) -> str:
    return json.dumps(to_json_like(value), sort_keys=True, default=str)


class CalibrationMemoisationTests(unittest.TestCase):
    def test_full_reference_calibration_runs_once_and_matches_uncached_result(self) -> None:
        config = build_unified_working_default_preset()["config"]
        policy._CALIBRATION_RESULT_CACHE.clear()
        policy._REFERENCE_ARTIFACT_CACHE.clear()
        with patch.object(policy, "calibrate_from_config", wraps=calibrate_from_config) as wrapped:
            resolved = policy.resolve_calibration_for_config(config)
            artifact = policy.build_reference_calibration_artifact(config)
        self.assertEqual(wrapped.call_count, 1)

        uncached = calibrate_from_config(policy.normalise_config(config))
        self.assertEqual(canonical({k: resolved[k] for k in uncached}), canonical(uncached))
        self.assertEqual(artifact["infectionIntercept"], uncached["ageInfLogLambda"])
        self.assertEqual(artifact["earlyProgressionHazard"], uncached["lambdaEarly"])
        self.assertEqual(artifact["lateProgressionHazard"], uncached["lambdaLate"])

    def test_cached_results_are_isolated_copies(self) -> None:
        config = build_unified_working_default_preset()["config"]
        first = policy._calibrate_once(policy.normalise_config(config))
        first["ageInfLogLambda"] = 123.0
        second = policy._calibrate_once(policy.normalise_config(config))
        self.assertNotEqual(second["ageInfLogLambda"], 123.0)


if __name__ == "__main__":
    unittest.main()
