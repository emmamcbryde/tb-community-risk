from __future__ import annotations

from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def abm_dir(root: Path | None = None) -> Path:
    return (root or repo_root()) / "abm"


def output_dir(root: Path | None = None) -> Path:
    return abm_dir(root) / "output"


def scenarios_dir(root: Path | None = None) -> Path:
    return abm_dir(root) / "scenarios"
