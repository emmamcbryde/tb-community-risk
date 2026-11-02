import pandas as pd
import streamlit as st
from country_model_review_ui.state import GATE_ITEMS,GATE_STATES,review_binding_hash,material_hash
from country_model_review_ui.presentation import text,choice,check,ensure_context
from country_model_review_ui.gates import gate_panel,audit_diagnostic
from country_model_review_ui.forward import RR74_WARNING


def render(config):
    st.header('Natural history and scientific gates');n=config['natural']
    st.warning('No country central parameter configuration has been approved. Sensitivity selections do not resolve incidence-inference gates.')
    n['adult']=choice('Adult progression curve',['External required','Conservative four-point sensitivity','Central-candidate four-point sensitivity','Higher-progression four-point sensitivity'],n['adult'])
    if n['adult']!='External required':st.caption('Menzies focal adult model at infection age20; historical cohort inputs and competing-risk transport unresolved. The2-year5% knot is not adopted.')
    n['children']=choice('Pediatric policy',['External required','Exclude children from reviewed diagnostic target scope'],n['children'])
    st.caption('No adult curve is applied below15. An exclusion requires matching population/observation scope and explicit review; child counts remain visible.')
    n['tail']=choice('Remote-infection tail',['No extrapolation','Fixed remote-hazard sensitivity','Constant-final-segment sensitivity'],n['tail'])
    if n['tail']=='Fixed remote-hazard sensitivity':n['remote_hazard']=text('Explicit remote hazard per year',n['remote_hazard'])
    else:st.caption('No extrapolation means unresolved beyond25years,never zero hazard.' if n['tail']=='No extrapolation' else 'Continuing the final segment is an unapproved tail sensitivity,not an evidence-based default.')
    st.selectbox('Reinfection clock',['Most recent reinfection resets clock'],key='edit_clock')
    st.caption('Repeated events remain possible. Recent reinfection overrides remote-only; prior remote exposure is retained without double counting.')
    n['mortality']=choice('Competing mortality',['No mortality — diagnostic only','ABS2022–2024 validation snapshot — sensitivity'],n['mortality'])
    if n['mortality']!='No mortality — diagnostic only':
        st.warning('ABS validation snapshot is not an adopted national default. It cannot be transferred to another country.100+ remains unresolved; all-cause/TB-death scope requires review.')
        n['mortality_proxy']=check('Explicitly use ABS mortality as a sensitivity proxy for this custom Australian population',n.get('mortality_proxy',False))
    n['active_baseline']=choice('Baseline active-TB risk-set policy',['Unresolved','Assume no baseline active TB — diagnostic only'],n['active_baseline'])
    n['horizon']=choice('Forward diagnostic horizon',['Instantaneous baseline','One-year prospective (bounded)'],n['horizon'])
    n['compare_scenarios']=check('Compare conservative,central-candidate and higher-progression sensitivities',n['compare_scenarios'])
    if n['compare_scenarios']:st.caption('Also compares unresolved versus final-segment remote tails and shows that future mortality does not change instantaneous baseline hazard. No mortality scenario is silently adopted.')
    st.caption('Instantaneous rates do not change with future mortality. Optional one-year first-episode projection uses the unchanged forward operator with depletion and competing death. No long-horizon country result is implied.')
    if config['risk']['legacy_effects']:st.warning(RR74_WARNING)
    st.subheader('Incidence-to-TBI scientific gate')
    d=config['diagnostic'];d['enabled']=check('Explicitly enable scientific incidence-to-TBI diagnostic',d['enabled'])
    panel=gate_panel(config);st.dataframe(pd.DataFrame(panel['rows']),hide_index=True,use_container_width=True)
    st.caption('Each gate is resolved/reviewed,sensitivity-only,unresolved,incompatible or blocked. Current sensitivity controls do not authorize inverse calibration.')
    if not d['enabled']:st.error('Calculation blocked — incidence inversion is disabled.');return
    st.warning('Experimental diagnostic only. Every required scientific gate,then root uniqueness and numerical refinement,must pass. Production availability remains false.')
    uploaded=st.file_uploader('Reviewed diagnostic JSON bundle',type=['json'],key='edit_diagnostic_upload')
    if uploaded and st.button('Use diagnostic JSON upload',key='use_diagnostic_upload'):
        st.session_state['edit_Reviewed diagnostic bundle JSON']=uploaded.getvalue().decode('utf-8-sig')
    d['bundle']=st.text_area('Reviewed diagnostic bundle JSON',value=d['bundle'],height=180,key='edit_Reviewed diagnostic bundle JSON')
    st.caption('The upload explicitly supplies the inverse diagnostic curves,mortality and observation contracts. It does not replace forward sensitivities or user TBI. Population,risk marginals and association must match this review configuration.')
    a,b=st.columns(2)
    with a:d['lower']=text('Scale search lower',d['lower'])
    with b:d['upper']=text('Scale search upper',d['upper'])
    d['grid_intervals']=st.selectbox('Identifiability grid intervals',[64,128,256,512],index=[64,128,256,512].index(d['grid_intervals']),key='edit_grid')
    st.caption('Review binding hash is in technical details. Uploads need native contracts,all15source/reviewer/rationale receipts,a compatible first-episode instantaneous target,approved exclusion threshold and a reviewed numerical protocol. No fictional defaults are supplied.')
    if st.button('Check scientific and numerical gates',type='primary',key='audit_diagnostic_gates'):
        context=ensure_context(config);report=audit_diagnostic(config,context);st.session_state['review_numerical_clearance']=report
        if report['passed']:st.success('Scientific pre-fit checks,root uniqueness and numerical refinement passed for this explicitly supplied diagnostic configuration.')
        else:st.error('Calculation blocked — numerical or identifiability gate failed; no TBI profile was accepted.')
    result=st.session_state.get('review_numerical_clearance')
    if result and result['configuration_hash']==material_hash(config):
        st.write('Numerical gate: '+('passed' if result['passed'] else 'blocked'))
        if result['passed']:
            st.caption('Uploaded inversion configuration: all required scientific receipts and numerical checks passed. Forward controls above remain sensitivity selections; no country or production approval is conferred.')
            st.dataframe(pd.DataFrame([{'Gate':item,'Uploaded diagnostic state':'Resolved and reviewed','Scope':'Source-bound reviewed configuration / checked numerical protocol'} for item in GATE_ITEMS]),hide_index=True)
        with st.expander('Gate audit details'):st.json(result)
    with st.expander('Technical configuration binding'):st.code(review_binding_hash(config))
