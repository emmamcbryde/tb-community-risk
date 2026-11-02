"""Display-only rounding and compact tables; scientific values stay unchanged."""
from collections import defaultdict
from copy import deepcopy
from math import fsum
import json
import pandas as pd
import altair as alt
import streamlit as st
from engine.country_tb.data_common import canonical,digest,seal
from engine.country_tb.risk_common import ratio
from country_model_review_ui.state import material_hash,stale,parse_csv,require

LABELS=('Observed','User supplied','Model derived','Diagnostic sensitivity','Unresolved','Calculation blocked')
COLORS=['#64748b','#177e89','#6247aa','#ba6b14','#94a3b8','#b91c1c']


def scalar(value):
    if value is None:return None
    return float(ratio(value))


def percent(value,digits=2):return 'Unavailable' if value is None else f'{scalar(value)*100:,.{digits}f}%'


def count(value):return 'Unavailable' if value is None else f'{scalar(value):,.0f}'


def text(label,default='',*,help=None):return st.text_input(label,value=str(default or ''),key='edit_'+label,help=help)


def choice(label,options,current,*,help=None):
    return st.selectbox(label,options,index=options.index(current) if current in options else 0,key='edit_'+label,help=help)


def check(label,value=False,*,help=None):return st.checkbox(label,value=value,key='edit_'+label,help=help)


def csv_input(label,current,kind):
    uploaded=st.file_uploader(label+' upload',type=['csv'],key='edit_'+label+'_upload')
    if uploaded is not None and st.button('Use '+label+' upload',key='use_'+kind):
        value=uploaded.getvalue();parse_csv(value,kind);st.session_state['edit_'+label+' CSV']=value.decode('utf-8-sig')
    return st.text_area(label+' CSV',value=current,height=180,key='edit_'+label+' CSV',help='UTF-8,comma-separated. Probabilities use0–1; blank means missing. Source/review status required.')


def ensure_context(config):
    from country_model_review_ui.model import select_population
    from time import perf_counter
    key=digest(canonical(config['population']))
    if st.session_state.get('review_population_key')!=key:
        start=perf_counter();st.session_state['review_population']=select_population(config)
        st.session_state['review_population_key']=key;st.session_state['review_timings']['population_selection_seconds']=perf_counter()-start
    return st.session_state['review_population']


def result_ready():
    if 'review_result' not in st.session_state:
        st.info('No diagnostic run yet. Enter inputs,then use Run diagnostic.');return None
    if stale(st.session_state):
        st.warning('Results stale — material inputs changed. Run again before interpreting results or exporting them.');return None
    return st.session_state['review_result']


def age_label(row):
    lo,hi=row['age_lower'],row['age_upper_exclusive']
    return f'{lo}+' if hi is None else str(lo) if hi==lo+1 else f'{lo}–{hi-1}'


def demographic_chart(rows):
    values=[]
    for r in rows:
        band=next((label for label,lo,hi in [('0–4',0,5),('5–14',5,15),('15–44',15,45),('45–64',45,65),('65+',65,None)]
            if r['age_lower']>=lo and (hi is None or r['age_upper_exclusive'] is not None and r['age_upper_exclusive']<=hi)),'Other supplied band')
        values.append(dict(Age=band,Sex=r['sex'],Population=float(r['population'])))
    df=pd.DataFrame(values).groupby(['Age','Sex'],sort=False,as_index=False).sum()
    chart=alt.Chart(df).mark_bar().encode(y=alt.Y('Age:N',sort=['0–4','5–14','15–44','45–64','65+']),
        x=alt.X('Population:Q',title='Population (people)'),color=alt.Color('Sex:N',scale=alt.Scale(range=['#177e89','#6247aa'])),
        tooltip=['Age','Sex',alt.Tooltip('Population:Q',format=',.0f')]).properties(height=230)
    st.altair_chart(chart,use_container_width=True)


def tbi_table(rows):
    return pd.DataFrame([{'Subgroup':r.get('subgroup_id','main'),'Age':age_label(r),'Sex':r['sex'],'Population':count(r['population']),
        'Any TBI':percent(r['any_tbi']),'Recent':percent(r['recent']),'Remote only':percent(r['remote_only']),
        'Prior remote + recent (subset)':percent(r['prior_remote_plus_recent'])} for r in rows])


def baseline_charts(rows):
    values=[dict(Age_cell=age_label(r),Age_lower=r['age_lower'],Sex=r['sex'],Subgroup=r.get('subgroup_id','main'),
        TBI_percent=float(r['any_tbi'])*100) for r in rows if r['any_tbi'] is not None]
    if values:
        st.altair_chart(alt.Chart(pd.DataFrame(values)).mark_point(filled=True,size=35).encode(
            x=alt.X('Age_lower:Q',title='Age cell lower boundary (years)'),y=alt.Y('TBI_percent:Q',title='Any TBI (%)'),
            color='Sex:N',shape='Subgroup:N',tooltip=['Age_cell','Sex','Subgroup',alt.Tooltip('TBI_percent:Q',format='.2f')]),use_container_width=True)
    state_rows=[dict(Age=age_label(r),Sex=r['sex'],State=label,Percent=float(r[key])*100,Population=float(r['population']),Subgroup=r.get('subgroup_id','main'))
        for r in rows for key,label in [('recent','Recent'),('remote_only','Remote only')] if r[key] is not None]
    if state_rows:
        df=pd.DataFrame(state_rows)
        df['Weighted']=df['Percent']*df['Population'];df=df.groupby(['Age','Sex','State'],sort=False,as_index=False)[['Weighted','Population']].sum()
        df['Percent']=df['Weighted']/df['Population']
        st.altair_chart(alt.Chart(df).mark_bar().encode(x=alt.X('Age:N',sort=None),y=alt.Y('Percent:Q',title='Covered population-weighted TBI state (%)'),
            color='State:N',column='Sex:N',tooltip=['Age','Sex','State',alt.Tooltip('Percent:Q',format='.2f')]),use_container_width=True)
        st.caption('Recent and remote-only are exclusive. Prior remote + recent is a subset of recent,never an extra stacked state. Subgroups are combined using population counts; unresolved cells are listed separately.')
    else:st.warning('Unresolved — recent/remote composition was not supplied or explicitly model-allocated.')


def prepare_exports(config,result):
    require(result['configuration_hash']==material_hash(config),'Cannot export a stale run as current.')
    bundle=seal(dict(schema_id='country_tb.review_ui_export.v1',configuration=deepcopy(config),run=deepcopy(result),
        scientific_approval='pending',production_available=False,externally_validated=False,
        gate_outcomes=result.get('gate') or __import__('country_model_review_ui.gates',fromlist=['gate_panel']).gate_panel(config),
        precision='source_strings_and_native_exact_rational_pairs_retained; no_display_rounding_applied'))
    rows=[]
    if 'baseline' in result:
        cells={(r['subgroup_id'],r['sex'],r['age_lower'],r['age_upper_exclusive']):r for r in result['risk_allocation']['rows']}
        for r in result['baseline']['rows']:
            row={k:r.get(k) for k in ('subgroup_id','sex','age_lower','age_upper_exclusive','population','any_tbi','recent','remote_only','never_infected','prior_remote_plus_recent')}
            cell=cells.get((r['subgroup_id'],r['sex'],r['age_lower'],r['age_upper_exclusive']))
            if cell:
                a=cell['association'];row.update(any_non_recency_risk_factor=None if cell['any_risk_prevalence'] is None else cell['any_risk_prevalence']['value'],
                    joint_TBI_any_risk=None if cell['joint_tbi_any_risk'] is None else cell['joint_tbi_any_risk']['value'],
                    requested_c=cell['normalized_association_input'],achieved_phi=None if a is None else a['achieved_phi'],
                    effective_RR_distribution_json=canonical(cell['effective_rr_distribution']).decode('utf-8'),
                    factor_marginals_json=canonical(cell['factor_marginals']).decode('utf-8'))
            rows.append(row)
    elif result.get('diagnostic',{}).get('profile'):
        rows=[dict(r,age_lower=r.get('age'),age_upper_exclusive=None) for r in result['diagnostic']['profile']['rows']]
    from country_model_review_ui.state import csv_text
    columns=('subgroup_id','sex','age_lower','age_upper_exclusive','population','any_tbi','recent','remote_only','never_infected','prior_remote_plus_recent',
        'any_non_recency_risk_factor','joint_TBI_any_risk','requested_c','achieved_phi','effective_RR_distribution_json','factor_marginals_json','run_scientific_hash','classification')
    for row in rows:row.update(run_scientific_hash=result['profile_hash'],classification='internal_review_not_country_estimate')
    return dict(json=canonical(bundle)+b'\n',csv=csv_text(rows,columns).encode('utf-8'),scientific_hash=bundle['profile_hash'])
