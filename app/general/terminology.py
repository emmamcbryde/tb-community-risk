"""User-facing terminology for the general application.

Ordinary general-application pages must not show setting-specific or internal
identifiers. ``FORBIDDEN_TERMS`` is used by tests and by ``display_text`` to keep
engine provenance strings out of the ordinary interface; the underlying
configuration and technical downloads keep full provenance.
"""

from __future__ import annotations

import re


APP_TITLE = "Community TB Screening Decision Support"
WORKFLOW_PAGES = (
    ("general_pages/1_Set_up_population.py", "Set up population"),
    ("general_pages/2_Define_intervention.py", "Define intervention"),
    ("general_pages/3_Run_analysis.py", "Run analysis"),
    ("general_pages/4_Results.py", "Results"),
    ("general_pages/5_Health_economics.py", "Health economics"),
    ("general_pages/6_Evidence_and_technical_information.py", "Evidence and technical information"),
)
RESTORE_DEFAULTS_LABEL = "Restore demonstration defaults"
USER_DEFINED_MARK = "✎ User-defined"
STOCHASTIC_LABEL = "Stochastic analysis - simulated populations"
DETERMINISTIC_LABEL = "Deterministic expected-value preview"
PLANNING_PRINCIPLE = (
    "This tool supports planning and sequencing of community TB screening. "
    "It is not a diagnostic tool and must not be used to deny care to anyone."
)

FORBIDDEN_PATTERNS = (
    re.compile(r"\bAPY\b", re.IGNORECASE),
    re.compile(r"APY\s+Lands", re.IGNORECASE),
    re.compile(r"SA\s+Health", re.IGNORECASE),
    re.compile(r"South\s+Australia", re.IGNORECASE),
    re.compile(r"MATLAB", re.IGNORECASE),
    re.compile(r"\bv9\b", re.IGNORECASE),
    re.compile(r"event[\s_-]*ledger", re.IGNORECASE),
    re.compile(r"contract[\s_-]*version", re.IGNORECASE),
    re.compile(r"\b(Rising|Steady|Falling)\b"),
)


def forbidden_matches(text: str) -> list[str]:
    return [match.group(0) for pattern in FORBIDDEN_PATTERNS for match in pattern.finditer(str(text))]


def display_text(text: str, *, fallback: str) -> str:
    """Return ``text`` unless it contains an internal or setting-specific identifier."""
    return fallback if forbidden_matches(text) else str(text)
