from __future__ import annotations

import streamlit as st

from app.display import arrow_safe_dataframe
from app.state import init_session_state
from engine.apy.evidence import assess_apy_reference_readiness, load_apy_evidence_registry


init_session_state()

st.title("Evidence & Assumptions")

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
st.dataframe(arrow_safe_dataframe(status_rows), use_container_width=True, hide_index=True)

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
        use_container_width=True,
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
        use_container_width=True,
        hide_index=True,
    )

conflict_rows = []
for key in ("evidenceConflicts", "bundlingConflicts"):
    for item in readiness.get(key) or []:
        conflict_rows.append({"Type": key, "Message": str(item)})
if conflict_rows:
    st.subheader("Conflicts")
    st.dataframe(arrow_safe_dataframe(conflict_rows), use_container_width=True, hide_index=True)

with st.expander("Evidence registry", expanded=False):
    registry = load_apy_evidence_registry()
    st.dataframe(arrow_safe_dataframe(registry), use_container_width=True, hide_index=True)

with st.expander("Caveats & technical information", expanded=False):
    st.markdown(
        """
        - The primary SA Health working reference uses the frozen software-compatible stochastic anchor; it excludes dynamic transmission effects.
        - Deterministic expected-value runs do not use repetitions or random seeds; stochastic runs use the selected repetition count and seed.
        - The provisional working route uses a compatibility placeholder for unresolved recent-versus-remote LTBI assumptions; this does not promote evidence readiness.
        - The inherited `10/770` active-TB calibration quantity remains unresolved and must not be interpreted as validated future progression from LTBI.
        - The implicit early/late progression structure is retained for compatibility with earlier APY analysis, not as measured recent-LTBI composition.
        - Disease-risk odds ratios are applied as multiplicative hazard multipliers for compatibility; this remains scientifically provisional.
        - Active or near-baseline TB is not fully separated from future incident, preventable TB in the compatibility anchor.
        - Simulation intervals describe finite-population stochastic variation, not full parameter uncertainty.
        - Programme setup, running, travel, outreach and staff-support costs remain not locally costed unless a user supplies local values.
        - DALY and ICER outputs are provisional; no willingness-to-pay threshold has been supplied, so NMB and probability cost-effective remain unavailable.
        - Contract versions, package hashes and release commits are retained in exported manifests and workbooks.
        """
    )
