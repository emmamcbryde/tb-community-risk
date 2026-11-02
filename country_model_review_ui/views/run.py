from time import perf_counter
import streamlit as st
from country_model_review_ui.state import SOURCE_MODES,material_hash
from country_model_review_ui.model import execute_user_run
from country_model_review_ui.gates import execute_diagnostic
from country_model_review_ui.presentation import ensure_context


def render(config):
    st.header('Run diagnostic');st.caption('Deterministic expected values first. No stochastic,intervention or economic outcomes are exposed.')
    if config['example']:st.warning(config['example'])
    user=config['tbi']['mode']==SOURCE_MODES[0]
    if user:
        st.write('Selected source: **Observed or user-specified TBI**')
        st.caption('A run retains source inputs,model-derived comparisons,unresolved populations and scientific hashes. Natural-history sensitivities do not calibrate TBI to incidence.')
        if st.button('Run deterministic review',type='primary',key='run_user'):
            with st.spinner('Constructing source-preserving review outputs…'):
                result,timings=execute_user_run(config,ensure_context(config));st.session_state['review_result']=result
                st.session_state['review_timings'].update(timings)
            st.success('Review run recorded. Open Results and validation.')
            if result['risk_allocation'].get('input_error'):st.warning('Risk allocation unresolved: '+result['risk_allocation']['input_error'])
            if result['projection']['status']=='calculation_blocked':st.warning('Forward calculation blocked: '+result['projection']['reason'])
    else:
        st.warning('Experimental diagnostic — requires reviewed natural-history assumptions')
        receipt=st.session_state.get('review_numerical_clearance')
        ready=config['diagnostic']['enabled'] and receipt is not None and receipt.get('passed') and receipt['configuration_hash']==material_hash(config)
        if not ready:st.error('Calculation blocked — enable and resolve every scientific and numerical gate first. Direct navigation does not bypass these checks.')
        if st.button('Run guarded incidence diagnostic',type='primary',key='run_incidence',disabled=not ready):
            start=perf_counter()
            with st.spinner('Rechecking bound review receipts and native3C numerical gates…'):
                result=execute_diagnostic(config,ensure_context(config),numerical_clearance=receipt)
            st.session_state['review_result']=result;st.session_state['review_timings']['incidence_diagnostic_seconds']=perf_counter()-start
            st.warning('Diagnostic only — not approved for production. '+result['status'])
    with st.expander('Recorded performance timings'):st.json(st.session_state['review_timings'])
