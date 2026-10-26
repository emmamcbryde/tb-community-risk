"""Session boundary for the isolated explicit-TBI app, never legacy state."""
from __future__ import annotations

from collections.abc import MutableMapping
from copy import deepcopy
from math import isfinite
from typing import Any

from app.explicit_tbi_workflow import (
    EXPLICIT_RECENT_REMOTE_PATHWAY_ID, WORKFLOW_PAGES, STOCHASTIC_LABEL,
    default_workflow_state, discard_incompatible_result,
)

PREFIX = "explicit_tbi_ui_v1:"
SCHEMA_KEY = PREFIX + "schema"
STATE_KEY = PREFIX + "state"
RESULT_KEY = PREFIX + "result"
STALE_KEY = PREFIX + "stale"
PAGE_KEY = PREFIX + "page"
WIDGET_VERSION_KEY = PREFIX + "widget_version"
DISPATCH_COUNT_KEY = PREFIX + "analysis_dispatch_count"
RENDER_TIMING_KEY = PREFIX + "render_seconds"
INPUT_ERRORS_KEY = PREFIX + "input_errors"
POPULATION_DRAFT_KEY = PREFIX + "population_draft"
WIDGET_PREFIX = PREFIX + "widget:"
SCHEMA = {"version": 1, "analysisPathway": EXPLICIT_RECENT_REMOTE_PATHWAY_ID}
OWN_KEYS = {
    SCHEMA_KEY, STATE_KEY, RESULT_KEY, STALE_KEY, PAGE_KEY, WIDGET_VERSION_KEY,
    DISPATCH_COUNT_KEY, RENDER_TIMING_KEY, INPUT_ERRORS_KEY, POPULATION_DRAFT_KEY,
}


def initialize_explicit_session(session: MutableMapping[str, Any]) -> bool:
    """Reject contaminated/unmarked state and open a clean Set up workflow.

    Sessions belong to this dedicated application. Unknown legacy keys (including
    economics, dynamic results and old explicit-app keys) are discarded together
    rather than imported or relabelled. Widget state is also invalidated on reset.
    """
    foreign_keys = [key for key in session if key not in OWN_KEYS and not key.startswith(WIDGET_PREFIX)]
    state = session.get(STATE_KEY)
    defaults = default_workflow_state()
    valid_state = isinstance(state, dict) and set(state) == set(defaults)
    if valid_state:
        valid_state = all(
            isinstance(state[key], type(default))
            or (isinstance(default, float) and isinstance(state[key], (int, float)))
            for key, default in defaults.items()
        )
        valid_state = valid_state and all(
            isfinite(float(state[key]))
            for key, default in defaults.items()
            if isinstance(default, (int, float))
        )
        valid_state = valid_state and state.get("analysisMethod") in (
            defaults["analysisMethod"], STOCHASTIC_LABEL,
        )
    reset = bool(foreign_keys) or session.get(SCHEMA_KEY) != SCHEMA or not valid_state
    rejected = bool(foreign_keys) or (STATE_KEY in session and reset)
    if reset:
        for key in list(session):
            del session[key]
        session[SCHEMA_KEY] = deepcopy(SCHEMA)
        session[STATE_KEY] = defaults
        session[RESULT_KEY] = None
        session[STALE_KEY] = False
        session[PAGE_KEY] = "Set up"
        session[WIDGET_VERSION_KEY] = 0
        session[DISPATCH_COUNT_KEY] = 0
        session[INPUT_ERRORS_KEY] = {}
    else:
        session.setdefault(RESULT_KEY, None)
        session.setdefault(STALE_KEY, False)
        session.setdefault(PAGE_KEY, "Set up")
        session.setdefault(WIDGET_VERSION_KEY, 0)
        session.setdefault(DISPATCH_COUNT_KEY, 0)
        session.setdefault(INPUT_ERRORS_KEY, {})
        result = session.get(RESULT_KEY)
        session[RESULT_KEY] = discard_incompatible_result(result)
        if result is not None and session[RESULT_KEY] is None:
            rejected = True
            session[PAGE_KEY] = "Set up"
        if session.get(PAGE_KEY) not in WORKFLOW_PAGES:
            session[PAGE_KEY] = "Set up"
    return rejected
