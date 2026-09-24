"""Report whether the active Python environment matches the supported configuration.

    python scripts/check_environment.py

Exit status is 0 when supported, 1 when a hard requirement is not met. Warnings
(for example x86-64 Python emulated on an ARM64 machine, which makes calibration
several times slower) do not change the exit status.
"""

from __future__ import annotations

from importlib import metadata
from pathlib import Path
import platform
import re
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SUPPORTED_PYTHON = ((3, 12), (3, 12))
SUPPORTED_STREAMLIT = ((1, 39, 0), (1, 54, 0))


def _version_tuple(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in re.findall(r"\d+", text)[:3])


def pinned_requirements(path: Path = REPO_ROOT / "requirements.txt") -> dict[str, str | None]:
    pins: dict[str, str | None] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name, _, version = line.partition("==")
        pins[name.strip().lower()] = version.strip() or None
    return pins


def check() -> tuple[list[str], list[str], dict[str, str]]:
    errors: list[str] = []
    warnings: list[str] = []
    info = {
        "python": sys.version.split()[0],
        "pythonBuild": platform.architecture()[0] + " " + ("AMD64" if "AMD64" in sys.version else platform.machine()),
        "machine": platform.machine(),
        "platform": platform.platform(),
    }
    if not SUPPORTED_PYTHON[0] <= sys.version_info[:2] <= SUPPORTED_PYTHON[1]:
        errors.append(f"Python {info['python']} is outside the supported range 3.12.x.")
    if platform.machine().upper() in {"ARM64", "AARCH64"} and "AMD64" in sys.version:
        warnings.append(
            "x86-64 Python is running under emulation on an ARM64 CPU; calibration-heavy tests are several times slower."
        )
    for name, pinned in pinned_requirements().items():
        try:
            installed = metadata.version(name)
        except metadata.PackageNotFoundError:
            errors.append(f"{name} is not installed.")
            continue
        info[name] = installed
        if name == "streamlit":
            if not SUPPORTED_STREAMLIT[0] <= _version_tuple(installed) <= SUPPORTED_STREAMLIT[1]:
                errors.append(f"Streamlit {installed} is outside the tested range 1.39.0-1.54.0.")
            elif pinned and installed != pinned:
                warnings.append(f"Streamlit {installed} differs from the pinned {pinned} (tested range 1.39.0-1.54.0).")
        elif pinned and installed != pinned:
            warnings.append(f"{name} {installed} differs from the pinned {pinned}.")
    return errors, warnings, info


def main() -> int:
    errors, warnings, info = check()
    for key, value in info.items():
        print(f"{key}: {value}")
    for message in warnings:
        print(f"WARNING: {message}")
    for message in errors:
        print(f"ERROR: {message}")
    print("Environment supported." if not errors else "Environment NOT supported.")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
