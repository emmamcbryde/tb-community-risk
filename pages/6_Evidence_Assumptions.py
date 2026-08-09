from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.state import init_session_state
from engine.apy.evidence import assess_apy_reference_readiness, load_apy_evidence_registry


init_session_state()

st.title("Evidence & Assumptions")
st.caption("Review readiness, provenance, conflicts and unresolved evidence for the APY demonstration.")

config = st.session_state.get("config") or {}
economics_config = st.session_state.get("economics_config") or {}

readiness = assess_apy_reference_readiness(config, economics_config)
status_rows = [
    {"Category": "Epidemiology", "Ready": readiness.get("epidemiologyReady")},
    {"Category": "Cost", "Ready": readiness.get("costReady")},
    {"Category": "Health outcomes", "Ready": readiness.get("dalyReady")},
    {"Category": "Decision benchmark", "Ready": readiness.get("thresholdReady")},
    {"Category": "Overall clinician-ready", "Ready": readiness.get("overallClinicianReady")},
]
st.subheader("Readiness")
st.dataframe(arrow_safe_dataframe(status_rows), width="content", hide_index=True)

if not readiness.get("overallClinicianReady"):
    st.warning(
        "Some assumptions remain unresolved or provisional. Results should be "
        "interpreted as modelled consequences, not clinician-ready conclusions."
    )

messages = readiness.get("messages") or []
if messages:
    st.subheader("Messages")
    st.dataframe(
        arrow_safe_dataframe([{"Message": str(message)} for message in messages]),
        width="stretch",
        hide_index=True,
    )

rows = readiness.get("readinessRows") or []
unresolved = [
    row for row in rows
    if not row.get("ready")
]
if unresolved:
    st.subheader("Unresolved Evidence")
    st.dataframe(
        arrow_safe_dataframe(
            [
                {
                    "Category": row.get("category"),
                    "Assumption": row.get("assumptionId"),
                    "Status": row.get("reviewStatus"),
                    "Source": row.get("sourceCitation"),
                    "Reason": row.get("unresolvedReason"),
                }
                for row in unresolved
            ]
        ),
        width="stretch",
        hide_index=True,
    )

conflict_rows = []
for key in ("evidenceConflicts", "bundlingConflicts"):
    for item in readiness.get(key) or []:
        conflict_rows.append({"Type": key, "Message": str(item)})
if conflict_rows:
    st.subheader("Conflicts")
    st.dataframe(arrow_safe_dataframe(conflict_rows), width="stretch", hide_index=True)

with st.expander("Evidence registry", expanded=False):
    registry = load_apy_evidence_registry()
    st.dataframe(arrow_safe_dataframe(registry), width="stretch", hide_index=True)
