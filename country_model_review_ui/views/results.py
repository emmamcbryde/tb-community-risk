from collections import defaultdict
from time import perf_counter
import pandas as pd
import altair as alt
import streamlit as st
from engine.country_tb.composition_common import profile_summary
from country_model_review_ui.state import material_hash
from country_model_review_ui.presentation import result_ready,percent,count,scalar,demographic_chart,tbi_table,baseline_charts,prepare_exports,age_label
from country_model_review_ui.forward import RR74_WARNING


def display_rows(rows,combined=False):
    if not combined:return rows
    groups=defaultdict(list)
    for r in rows:groups[r.get('subgroup_id','main'),r['age_lower'],r['age_upper_exclusive']].append(r)
    output=[]
    for (g,lo,hi),values in groups.items():
        summary=profile_summary(values)
        output.append(dict(subgroup_id=g,age_lower=lo,age_upper_exclusive=hi,sex='derived both sexes',population=summary['population'],
            **{k:v['prevalence_full'] for k,v in summary['components'].items()}))
    return output


def render(config):
    start=perf_counter();st.header('Results and validation');result=result_ready()
    if result is None:return
    if config['example']:st.warning(config['example'])
    st.caption('Observed / User supplied / Model derived / Diagnostic sensitivity / Unresolved / Calculation blocked')
    st.warning('Internal scientific review only. These outputs are not externally validated national TBI estimates.')
    if 'diagnostic' in result:
        value=result['diagnostic'];st.write('Experimental diagnostic status: **'+value['status']+'**')
        st.caption('Public production incidence-inference status remains '+result['public_status'])
        st.json({k:value.get(k) for k in ('fitted_scale','achieved_incidence','target_incidence','failures')})
        if value.get('profile'):
            st.dataframe(pd.DataFrame(value['profile']['rows']),hide_index=True)
            st.caption('Model-derived counterfactual for the explicitly uploaded diagnostic scope; not an empirical country estimate.')
    else:
        baseline=result['baseline'];risk=result['risk_allocation'];projection=result['projection'];s=baseline['summary']
        a,b,c=st.columns(3);a.metric('Population retained',count(s['population']));b.metric('Any TBI — full supplied population',percent(s['components']['any_tbi']['prevalence_full']));c.metric('Unresolved total-TBI population',count(s['unresolved_population']))
        st.caption('Any TBI over covered population: '+percent(s['components']['any_tbi']['prevalence_covered'])+'. Full-population prevalence is unavailable when coverage is incomplete.')
        st.subheader('Age and sex population');demographic_chart(result['context']['demography']['rows'])
        st.caption('Population source: '+('User-supplied CSV' if result['context']['source_profile'] is None else 'Frozen ABS / UN demographic source;2024 baseline')+'. Female/male records remain separate internally. Snapshot IDs are in technical provenance.')
        st.subheader('TBI profile')
        t=config['tbi']
        source_label='User supplied anchor; model-derived age and recency' if t['derive'] and t['input_kind']=='Overall prevalence' else 'User supplied totals; model-derived hybrid recency' if t['hybrid_recency'] else 'Supplied assay positivity; model-derived latent TBI' if t['assay']!='Direct latent TBI' else 'User supplied total or direct states'
        st.caption(source_label+' — supplied source values and model comparison profiles are retained separately. Assay positivity is not latent infection.')
        shown=display_rows(baseline['rows'],config['display']['sexes']=='Combined display')
        st.dataframe(tbi_table(shown),hide_index=True,use_container_width=True);baseline_charts(baseline['rows'])
        st.subheader('Population-weighted subgroup contributions')
        groups=[]
        for g in baseline['subgroups']:
            total=g['tbi_profile']['summary']['components']['any_tbi']
            groups.append({'Subgroup':g['label'],'Population':count(g['population']),'TBI prevalence':percent(total['prevalence_full']),
                'Covered TBI people':count(total['estimated_count']),'Unresolved TBI people':count(total['unresolved_population'])})
        for g in baseline['unresolved_subgroups']:
            groups.append({'Subgroup':g['subgroup_id'],'Population':count(g['population']),'TBI prevalence':'Unresolved',
                'Covered TBI people':'Unavailable','Unresolved TBI people':count(g['population'])})
        st.dataframe(pd.DataFrame(groups),hide_index=True,use_container_width=True)
        if baseline['unresolved_subgroups']:st.error('Incomplete mixture — a complete population TBI result is unavailable. Unresolved subgroup denominators are retained.')
        st.subheader('Risk factors and association')
        if risk.get('input_error'):st.error('Risk allocation unresolved: '+risk['input_error'])
        associations=[];factors=[];rr=defaultdict(float)
        for cell in risk['rows']:
            a=cell['association']
            associations.append(dict(Subgroup=cell['subgroup_id'],Age=age_label(cell),Sex=cell['sex'],
                Any_risk=percent(cell['any_risk_prevalence']),Requested_c=cell['normalized_association_input'],
                Achieved_phi='Undefined / unavailable' if a is None or a['achieved_phi'] is None else f"{float(a['achieved_phi']):.2f}",
                Feasible_lower='Unavailable' if a is None else percent(a['lower_bound']),Feasible_upper='Unavailable' if a is None else percent(a['upper_bound']),
                Joint='Unavailable' if a is None else percent(a['joint_prevalence'])))
            for f in cell['factor_marginals']:
                factors.append(dict(Factor=f['factor_id'],Subgroup=cell['subgroup_id'],Age=age_label(cell),Sex=cell['sex'],
                    Prevalence=percent(f['prevalence']),State=f['state'],Definition=f['compatibility_status']))
            for row in cell['effective_rr_distribution'] or []:rr['Unresolved' if row['effective_rr'] is None else row['effective_rr']['value']]+=scalar(row['expected_population'])
        if associations:
            st.dataframe(pd.DataFrame(associations),hide_index=True,use_container_width=True)
            st.caption('Requested c is a normalized feasible control,not achieved phi. Phi is undefined for degenerate margins. All controls are within exact source age/sex/subgroup cells.')
        if factors:
            with st.expander('Exposure prevalences,states and definition compatibility'):st.dataframe(pd.DataFrame(factors),hide_index=True,use_container_width=True)
            values=[]
            for cell in risk['rows']:
                for f in cell['factor_marginals']:
                    if f['state']=='value' and f['prevalence'] is not None:
                        values.append(dict(Age=cell['age_lower'],Sex=cell['sex'],Subgroup=cell['subgroup_id'],Factor=f['factor_id'],Percent=scalar(f['prevalence'])*100))
            if values:st.altair_chart(alt.Chart(pd.DataFrame(values)).mark_point(filled=True).encode(x=alt.X('Age:Q',title='Age cell lower boundary'),y=alt.Y('Percent:Q',title='Risk-factor prevalence (%)'),color='Sex:N',shape='Factor:N',tooltip=['Age','Sex','Factor','Subgroup',alt.Tooltip('Percent:Q',format='.2f')]),use_container_width=True)
        if rr:
            st.subheader('Highest applicable RR distribution')
            st.bar_chart(pd.DataFrame({'Highest RR':list(rr),'Expected people':list(rr.values())}).set_index('Highest RR'))
            st.caption('Population counts by maximum applicable RR; no additive or multiplicative central combination.')
            if '74' in rr:st.warning(RR74_WARNING)
        st.subheader('Active-TB forward sensitivity and contextual incidence')
        st.caption('First episode after baseline only. Baseline active TB is excluded by the selected diagnostic assumption; prior TB,recurrence,imported disease,migration and historical depletion remain unrepresented.')
        value=projection.get('instantaneous')
        if value is None:st.error('Calculation blocked — '+projection.get('reason','Unresolved forward inputs'))
        else:
            st.metric('Diagnostic instantaneous rate — complete supplied population',
                'Unavailable' if value['full_population_rate_per100000PY'] is None else f"{value['full_population_rate_per100000PY']:,.2f} /100,000 PY")
            st.metric('Unresolved progression population',count(value['unresolved_population']))
            if value['unresolved_population']:
                st.warning('Partial sensitivity only. Missing child,remote-tail,effect or infection-time hazards are unresolved,not zero. Full-population incidence is unavailable.')
            st.caption('Covered case flow: '+f"{value['covered_instantaneous_case_flow_per_year']:.2f}"+' expected first episodes/year at the baseline instant. This is not cumulative one-year cases.')
            target=result['context']['observed_incidence'];observed=None if target is None or target['central'] is None else float(target['central'])
            comparison=[]
            if value['full_population_rate_per100000PY'] is not None:comparison.append({'Source':'Diagnostic sensitivity — supplied population','Rate':value['full_population_rate_per100000PY'],'Type':'Diagnostic sensitivity'})
            if observed is not None:comparison.append({'Source':'Observed WHO national estimate — context only','Rate':observed,'Type':'Observed'})
            if comparison:
                st.altair_chart(alt.Chart(pd.DataFrame(comparison)).mark_bar().encode(y='Source:N',x=alt.X('Rate:Q',title='Rate per100,000 person-years'),color='Type:N',tooltip=['Source',alt.Tooltip('Rate:Q',format='.2f')]),use_container_width=True)
            else:st.warning('Observed incidence unavailable; no contextual value is invented.')
            st.warning('Validation only,not calibration. WHO annual new/recurrent country incidence is not the same population/episode/time estimand as this instantaneous first-episode diagnostic. Disagreement never changes supplied TBI.')
            if projection['scenarios']:
                st.subheader('Progression,remote-tail and baseline-mortality scenarios')
                st.dataframe(pd.DataFrame([{'Scenario':s['parameter_set_id'],'Full-population instantaneous rate':s['full_population_rate_per100000PY'],
                    'Covered case flow/year':s['covered_instantaneous_case_flow_per_year'],'Unresolved population':s['unresolved_population']} for s in projection['scenarios']]),hide_index=True,
                    column_config={'Full-population instantaneous rate':st.column_config.NumberColumn(format='%.2f'),
                        'Covered case flow/year':st.column_config.NumberColumn(format='%.2f'),'Unresolved population':st.column_config.NumberColumn(format='%.0f')})
            st.caption(value['mortality_interpretation'])
            if projection.get('prospective'):
                st.subheader('One-year prospective diagnostic — first episode only')
                year=projection['prospective']['projections'][0]
                a,b,c,d=st.columns(4)
                a.metric('Expected one-year cases — covered scope',f"{float(year['expected_cases_covered']):,.2f}")
                b.metric('One-year incidence /100,000 PY — full scope','Unavailable' if year['incidence_per_100000_person_years_full'] is None else f"{float(year['incidence_per_100000_person_years_full']):,.2f}")
                c.metric('Competing deaths — covered scope',f"{float(year['competing_deaths_covered']):,.2f}")
                d.metric('Surviving disease-free — covered scope',f"{float(year['surviving_disease_free_covered']):,.2f}")
                st.caption('Mortality policy: '+config['natural']['mortality']+'. Cases and deaths leave the first-episode risk set; missing population remains unresolved.')
                with st.expander('Prospective counts,depletion and competing deaths'):st.json(projection['prospective'])
            elif projection.get('prospective_unavailable_reason'):st.warning('Prospective calculation blocked: '+projection['prospective_unavailable_reason'])
        with st.expander('Assumptions,unresolved rows and source comparisons'):
            st.json({'natural_history':config['natural'],'forward_unresolved':None if value is None else value['rows'],
                'baseline_warnings':[g['tbi_profile']['warnings'] for g in baseline['subgroups']]})
    st.subheader('Review exports')
    st.caption('Full-precision source values,native contracts,gate outcomes,unresolved populations and provenance. Preparation is explicit; no expensive export runs during ordinary rendering.')
    if st.button('Prepare full-precision exports',key='prepare_review_exports'):
        st.session_state['review_exports']=prepare_exports(config,result);st.session_state['review_exports_hash']=material_hash(config)
    if st.session_state.get('review_exports_hash')==material_hash(config) and st.session_state.get('review_exports'):
        e=st.session_state['review_exports'];a,b=st.columns(2)
        a.download_button('Download review JSON',e['json'],'country_TB_scientific_review.json','application/json',key='download_review_json')
        b.download_button('Download full-precision profile CSV',e['csv'],'country_TB_review_profile.csv','text/csv',key='download_review_csv')
    with st.expander('Technical provenance and timings'):
        st.json(result.get('provenance',{'run_hash':result['profile_hash'],'configuration_hash':result['configuration_hash']}));st.json(st.session_state['review_timings'])
    st.session_state['review_timings']['results_rendering_seconds']=perf_counter()-start
