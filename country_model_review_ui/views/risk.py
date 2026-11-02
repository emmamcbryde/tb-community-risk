from copy import deepcopy
import pandas as pd
import streamlit as st
from country_model_review_ui.model import catalog
from engine.country_tb.risk_factor_catalog import catalogue_records
from engine.country_tb.risk_policies import feasible_association
from engine.country_tb.data_common import canonical
from country_model_review_ui.state import RISK_COLUMNS,ASSOCIATION_COLUMNS,csv_text,percent_input
from country_model_review_ui.presentation import text,choice,check,csv_input,percent
from country_model_review_ui.forward import RR74_WARNING


def render(config):
    st.header('Risk factors');r=config['risk'];entries=catalogue_records(catalog())
    st.info('Central multimorbidity: highest applicable reviewed relative risk. Factors are not added or multiplied. Adult evidence is not automatically applied below15.')
    r['no_factors']=check('Explicit RR1 reference: exclude all non-recency factors',r['no_factors'],help='A diagnostic reference scenario,not an assertion that a country has no risk factors.')
    r['legacy_effects']=check('Enable legacy adult effects for diagnostic sensitivity only',r['legacy_effects'])
    if r['legacy_effects']:st.warning('Legacy effects are unreviewed as conditional progression hazards. Highest-RR combination does not approve their scientific use.')
    if not r['no_factors']:
        st.subheader('Enter an exposure record')
        labels={e.get('display_label',k.replace('_',' '))+' ('+k+')':k for k,e in entries.items()}
        fid=labels[st.selectbox('Risk factor',list(labels),key='edit_factor')]
        entry=entries[fid];definition=entry['exposure_definition']
        st.caption('Catalogue definition: '+definition)
        if fid=='transplant_immunosuppression':st.warning(RR74_WARNING)
        groups=[g['subgroup_id'] for g in config['population']['subgroups']] or ['main']
        subgroup=st.selectbox('Risk-record subgroup',groups+(['all'] if len(groups)>1 else []),key='edit_risk_subgroup')
        sex=choice('Risk-record sex',['both','female','male'],'both')
        lo=text('Exposure age lower','15');hi=text('Exposure age upper exclusive (blank=open)','')
        status=choice('Prevalence state',['value','missing','not_applicable','excluded','not_reviewed','unavailable_due_to_definition_mismatch'],'missing')
        value=text('Risk-factor prevalence (%)','')
        supplied=text('Supplied exposure definition',definition)
        source=text('Risk prevalence source','Reviewer entry');review=choice('Risk input review status',['not_reviewed','reviewed','diagnostic'],'not_reviewed')
        reason=text('Missing/exclusion/applicability reason','')
        if supplied!=definition:st.warning('Definition mismatch — no automatic mapping. Current smoking is not>one pack/day; BMI is not percent ideal weight.')
        if st.button('Add exposure record',key='add_risk_record'):
            from country_model_review_ui.state import require
            require(lo.isdigit() and (not hi or hi.isdigit() and int(hi)>int(lo)),'Use aligned integer age boundaries.')
            prevalence=percent_input(value,required=status=='value')
            record=dict(record_id='risk_'+str(len(r['records'])+1),subgroup_id=subgroup,factor_id=fid,age_lower=int(lo),age_upper_exclusive=None if not hi else int(hi),
                sex=sex,state=status,prevalence=prevalence,allocation_policy='piecewise_constant_within_supplied_band_v1',
                supplied_exposure_definition=supplied,source=source,review_status=review,reason=reason)
            r['records'].append(record)
        if r['records']:
            st.dataframe(pd.DataFrame(r['records']),hide_index=True,use_container_width=True)
            st.button('Clear entered exposure records',key='clear_risk_records',on_click=lambda:r.update(records=[]))
        r['csv']=csv_input('Risk factors',r['csv'],'risk')
    elif r['records'] or r['csv']:
        st.error('RR1 reference conflicts with existing risk records. Clear the records or unselect RR1; there is no silent override.')
        def clear_reference():r.update(records=[],csv='');st.session_state.pop('edit_Risk factors CSV',None)
        st.button('Clear risk inputs for RR1 reference',key='clear_reference_risks',on_click=clear_reference)
    st.subheader('TBI–any-risk association within age/sex strata')
    a=config['association'];a['global_control']=text('Normalized association c (−1 to +1)',a['global_control'],help='−1 minimum feasible overlap;0 independence;+1 maximum feasible overlap. This is not achieved Pearson correlation.')
    st.caption('The control preserves TBI,each factor,any-risk union,population and recent/remote totals. Residual association is applied after age/sex structure.')
    a['csv']=csv_input('Association overrides',a['csv'],'association')
    last=st.session_state.get('review_result')
    if last and 'risk_allocation' in last:
        prior=deepcopy(last['configuration']);current=deepcopy(config);prior.pop('association');current.pop('association')
        prior.pop('display',None);current.pop('display',None)
        if canonical(prior)==canonical(current):
            values=[]
            for cell in last['risk_allocation']['rows']:
                if cell['tbi_prevalence'] is None or cell['any_risk_prevalence'] is None:continue
                result=feasible_association(cell['tbi_prevalence'],cell['any_risk_prevalence'],a['global_control'])
                values.append(dict(Age=str(cell['age_lower']),Sex=cell['sex'],Subgroup=cell['subgroup_id'],
                    Requested_c=a['global_control'],Achieved_phi=result['achieved_phi'] or 'Undefined',
                    Feasible_lower=percent(result['lower_bound']),Feasible_upper=percent(result['upper_bound']),Joint=percent(result['joint_prevalence'])))
            if values:st.caption('Association preview on the last computed marginals; rerun to record changes. CSV overrides are applied on execution.');st.dataframe(pd.DataFrame(values),hide_index=True)
    with st.expander('Candidate catalogue and exposure definitions'):
        st.warning(RR74_WARNING)
        st.dataframe(pd.DataFrame([dict(Factor=k,Definition=e['exposure_definition'],Candidate_RR=e['candidate_effect']['value'] if isinstance(e.get('candidate_effect'),dict) else e.get('candidate_effect'),
            Status='Diagnostic legacy; contemporary review pending') for k,e in entries.items()]),hide_index=True)
        st.caption('Recent infection and close contact are excluded from ordinary factor inputs. Recency comes from the infection model.')
    with st.expander('Risk and association CSV schemas'):
        st.code(','.join(RISK_COLUMNS),language='text');st.code(','.join(ASSOCIATION_COLUMNS),language='text')
        st.caption('Explicit allocation policy required; exact single-year,banded,both-sex or sex-specific values. Override precedence is exact age/sex,exact age/both,band/sex,band/both,sex,global. Overlaps at equal precedence fail.')
