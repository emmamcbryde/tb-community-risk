from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import streamlit as st


ACTIVITY_ICONS = {
    "initialising": "[ ]",
    "expected_value": "[=]",
    "simulating": "[>]",
    "summarising": "[#]",
    "economics": "[$]",
    "finalising": "[✓]",
}


@dataclass(frozen=True)
class ProgressStatus:
    activity: str
    label: str
    progress: float
    replicate: int | None = None
    total_replicates: int | None = None

    @property
    def message(self) -> str:
        if self.replicate is not None and self.total_replicates:
            return f"{self.label}: replicate {self.replicate} of {self.total_replicates}"
        return self.label


def initialising_status() -> ProgressStatus:
    return ProgressStatus("initialising", "Initialising analysis", 0.05)


def expected_value_status() -> ProgressStatus:
    return ProgressStatus("expected_value", "Running expected-outcomes analysis", 0.5)


def replicate_status(replicate: int, total_replicates: int) -> ProgressStatus:
    total = max(int(total_replicates), 1)
    current = max(0, min(int(replicate), total))
    progress = 0.1 + 0.7 * (current / total)
    return ProgressStatus(
        "simulating",
        "Simulating cohort",
        progress,
        replicate=current,
        total_replicates=total,
    )


def summarising_status() -> ProgressStatus:
    return ProgressStatus("summarising", "Summarising outcomes", 0.85)


def economics_status() -> ProgressStatus:
    return ProgressStatus("economics", "Calculating economics", 0.92)


def finalising_status() -> ProgressStatus:
    return ProgressStatus("finalising", "Finalising results", 1.0)


def status_from_event(event: dict[str, Any]) -> ProgressStatus:
    stage = str(event.get("stage") or "")
    if stage == "replicate_completed":
        return replicate_status(
            int(event.get("replicate", 0)),
            int(event.get("totalReplicates", 1)),
        )
    if stage == "summarising":
        return summarising_status()
    if stage == "finalising":
        return finalising_status()
    if stage == "expected_value":
        return expected_value_status()
    return initialising_status()


class StreamlitProgressDisplay:
    def __init__(self) -> None:
        self._status = st.empty()
        self._bar = st.progress(0.0)

    def update(self, status: ProgressStatus) -> None:
        icon = ACTIVITY_ICONS.get(status.activity, "[ ]")
        self._status.write(f"{icon} {status.message}")
        self._bar.progress(max(0.0, min(float(status.progress), 1.0)))

    def callback(self, event: dict[str, Any]) -> None:
        self.update(status_from_event(event))
