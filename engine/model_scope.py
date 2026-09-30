"""Canonical scope statements for the no-transmission screening model.

Interface pages, exports and profile notes reuse these strings so the model
boundary is described the same way everywhere. See
``docs/no_transmission_model_scope.md`` and
``docs/catalytic_infection_pressure_spec.md``.
"""

from __future__ import annotations

import re


MODEL_IDENTITY = (
    "A screening and preventive-treatment model with no endogenous transmission feedback."
)
DIRECT_EFFECTS_STATEMENT = (
    "Results are direct outcomes among the modelled population, compared with no screening. "
    "The model does not estimate onward transmission, herd or other indirect effects."
)
INCIDENCE_DESCRIPTIVE_STATEMENT = (
    "Country incidence data describe the TB disease burden and trend. They do not set "
    "infection pressure or transmission in this model and do not change its epidemiology."
)
INCIDENCE_TO_INFECTION_POLICY = (
    "Estimated TB disease incidence is never used as infection pressure. Any background exposure "
    "to TB infection must be supplied as its own input, with its own source."
)
NO_FEEDBACK_ANALYSIS_LABEL = "Individual-based no-feedback analysis"

# Claims this model must never make about its own outputs (docs/no_transmission_model_scope.md
# section 6). Tests scan ordinary interface text and exports with these patterns.
PROHIBITED_CLAIM_PATTERNS = (
    re.compile(r"secondary (infections?|cases?) (prevented|averted)", re.IGNORECASE),
    re.compile(r"(prevented|averted) (secondary|onward)", re.IGNORECASE),
    re.compile(r"(reduc\w+|lower\w*) (in )?(community|population)[- ](level )?(tb )?(incidence|transmission)", re.IGNORECASE),
    re.compile(r"(community|population)[- ]level (transmission|incidence) reduction", re.IGNORECASE),
    re.compile(r"herd (effect|protection|immunity)", re.IGNORECASE),
    re.compile(r"reduc\w+ (the )?force of infection", re.IGNORECASE),
    re.compile(r"(progress|contribution) towards? elimination", re.IGNORECASE),
    re.compile(r"(outbreaks?|clusters?) (prevented|averted|reduced)", re.IGNORECASE),
    re.compile(r"transmission[- ](mediated )?benefits? (are|is) not yet", re.IGNORECASE),
    re.compile(r"not yet (used|included|modelled) .{0,40}transmission", re.IGNORECASE),
)


def prohibited_claims(text: str) -> list[str]:
    """Return any prohibited transmission-effect claims found in ``text``."""
    return [match.group(0) for pattern in PROHIBITED_CLAIM_PATTERNS for match in pattern.finditer(str(text))]
