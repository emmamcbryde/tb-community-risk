"""Guard the frozen SA Health release against changes from general-application work."""

from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]
RELEASE_TAG = "sa-health-apy-he-v1.0.0"
RELEASE_BRANCH = "release/sa-health-apy-he-v1.0.0"
FROZEN_COMMIT = "03cc16e4a52a55e10019dab16bd4f599571c4b90"
# SHA-256 of the committed blobs at the release tag. Text files may be checked out
# with CRLF line endings on Windows, so they are compared after CRLF -> LF.
FROZEN_REFERENCE_SHA256 = {
    "app/reference_data/sa_health_frozen_event_ledger_annual.csv.gz": "3c2f10ba93fb17db3d09976c0ee0dadf0f6e8b0bfe491628445ec363be1ad038",
    "app/reference_data/sa_health_frozen_event_ledger_totals.csv.gz": "bbe3ae6e0ab305577c2c026c24641e63cc518b69118144f6f62841aefe707f98",
    "app/reference_data/sa_health_frozen_primary_economic_annual_by_arm.csv.gz": "34b0df45ea3aaaee99df7aa8767bd983cad5960b108961232fe3ee44086130fa",
    "app/reference_data/sa_health_frozen_primary_economic_replicates.csv.gz": "a4b86f96223c09c0fdeb7e8fe23172c39f49fc390790b5310e690db21c74d68c",
    "app/reference_data/sa_health_frozen_reference_config.json": "de3f1f2f6eb3401b025a0cc90e2f540462f20219d511024e66e38066af74e54a",
    "app/reference_data/sa_health_frozen_reference_economics_config.json": "a24029527f177883f7ee8d15c72338a27aa7d499c0254693acf76220cf679502",
    "app/reference_data/sa_health_frozen_reference_manifest.json": "dfda0a26aa4fa32e6a6afe6baf8571c82c48a3c7be61905821aa9766a64e4e33",
    "app/reference_data/sa_health_report_reference_economics.json": "0eead1c56807dca3520b64567a7006aa47c66bc2fedbe1a418bcff787c6f5451",
}
WORKING_DEFAULT_PRESET_HASH = "b987841c073e2d33ccaf72f053a00b0a4d0af2e4fb835c595ac3767655b9b056"
# Files that existed at the release tag and are intentionally edited on the general branch.
# engine/apy/calibration_policy.py: memoises a duplicated, identical calibration call
# (tests/test_calibration_memoisation.py proves outputs are unchanged).
ALLOWED_MODIFIED_RELEASE_FILES = {".gitignore", "README.md", "engine/apy/calibration_policy.py"}


def _git(*args: str) -> str | None:
    if shutil.which("git") is None:
        return None
    try:
        result = subprocess.run(["git", "-C", str(ROOT), *args], check=True, capture_output=True, text=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


class FrozenReleaseIntegrityTests(unittest.TestCase):
    def test_frozen_reference_artifacts_are_byte_identical(self) -> None:
        for relative, expected in FROZEN_REFERENCE_SHA256.items():
            with self.subTest(file=relative):
                raw = (ROOT / relative).read_bytes()
                normalised = raw.replace(bytes([13, 10]), bytes([10]))
                digests = {hashlib.sha256(raw).hexdigest(), hashlib.sha256(normalised).hexdigest()}
                self.assertIn(expected, digests)

    def test_release_tag_and_branch_point_at_frozen_commit(self) -> None:
        tag_commit = _git("rev-parse", "--verify", "--quiet", f"refs/tags/{RELEASE_TAG}^{{commit}}")
        if not tag_commit:
            self.skipTest("Release tag is not available in this checkout.")
        self.assertEqual(tag_commit, FROZEN_COMMIT)
        for ref in (f"refs/heads/{RELEASE_BRANCH}", f"refs/remotes/origin/{RELEASE_BRANCH}"):
            commit = _git("rev-parse", "--verify", "--quiet", ref)
            if commit:
                with self.subTest(ref=ref):
                    self.assertEqual(commit, FROZEN_COMMIT)

    def test_release_files_unchanged_except_documented_exceptions(self) -> None:
        if not _git("rev-parse", "--verify", "--quiet", f"refs/tags/{RELEASE_TAG}^{{commit}}"):
            self.skipTest("Release tag is not available in this checkout.")
        changed = _git("diff", "--name-only", "--diff-filter=MDRT", RELEASE_TAG, "--")
        self.assertIsNotNone(changed)
        modified = {line for line in changed.splitlines() if line}
        self.assertEqual(modified - ALLOWED_MODIFIED_RELEASE_FILES, set())

    def test_general_work_does_not_alter_frozen_numerical_reference(self) -> None:
        from engine.apy.frozen_reference import is_frozen_sa_health_reference_eligible, load_frozen_reference_results
        from engine.apy.working_defaults import build_unified_working_default_preset
        from engine.profiles.demonstration import build_demonstration_profile
        from engine.profiles.engine_mapping import build_engine_config

        before = load_frozen_reference_results()["resultsBundle"]["headline"]["keyMetricsRows"]
        general_config = build_engine_config(build_demonstration_profile())
        preset = build_unified_working_default_preset()
        after = load_frozen_reference_results()["resultsBundle"]["headline"]["keyMetricsRows"]

        self.assertEqual(before, after)
        self.assertEqual(preset["configurationHash"], WORKING_DEFAULT_PRESET_HASH)
        self.assertTrue(is_frozen_sa_health_reference_eligible(preset["config"]))
        self.assertFalse(is_frozen_sa_health_reference_eligible(general_config))


if __name__ == "__main__":
    unittest.main()
