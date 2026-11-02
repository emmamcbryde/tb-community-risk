from copy import deepcopy
import pandas as pd
import streamlit as st
from country_model_review_ui.state import working_example,restore,POP_COLUMNS,csv_text
from country_model_review_ui.model import countries,geography
from country_model_review_ui.presentation import text,choice,check,csv_input,ensure_context,count,demographic_chart


def render(config):
    st.header('Set up population');st.caption('Frozen demographic records or an explicitly defined review population. Baseline 2024.')
    c1,c2=st.columns(2)
    c1.button('Restore working example',key='restore_working_example',on_click=lambda:restore(st.session_state,working_example()))
    c2.button('Load synthetic subgroup example',key='restore_subgroup_example',on_click=lambda:restore(st.session_state,working_example(mixture=True)))
    p=config['population'];index={c['country_id']:c['display']['name'] for c in countries()}
    labels={name+' ('+code+')':code for code,name in index.items()}
    selected=st.selectbox('Country or territory',list(labels),index=list(labels.values()).index(p['country']),key='edit_country')
    p['country']=labels[selected]
    p['scope']=choice('Population scope',['National','Custom / subnational'],p['scope'])
    p['label']=text('Population label',p['label'])
    if p['scope']!='National':p['custom_id']=text('Custom geography ID',p['custom_id'],help='Use custom: or subnational:. Current residence is not assumed to be historical exposure geography.')
    sources=['Recommended frozen source']+(['ABS ERP June2024'] if p['country']=='AUS' else [])+['UN WPP2024 medium','User CSV']
    p['source']=choice('Demographic source',sources,p['source'])
    st.selectbox('Baseline year',[2024],key='edit_baseline_year');st.caption('2024 is the latest observation year in the frozen WHO2025 snapshot. Other baseline years require a future contract.')
    if p['source']=='User CSV':
        p['csv']=csv_input('Population',p['csv'],'population');p['size']=text('Declared population (optional check against CSV)',p['size'])
        st.caption('Population CSV must retain separate female and male counts,including open ages. See the documented upload schema.')
    else:
        p['scale_source']=check('Scale source counts to a reviewer-defined population',p['scale_source'])
        if p['scale_source'] or p['scope']!='National':p['size']=text('Reviewer-defined population size',p['size'])
        if p['scope']!='National':p['national_proxy']=check('Use national age/sex proportions as an explicit diagnostic proxy',p['national_proxy'],
            help='This does not establish that national demography or incidence represents the subgroup.')
        def restore_source():
            p.update(scale_source=False,size='',national_proxy=False);st.session_state.pop('edit_Scale source counts to a reviewer-defined population',None)
        st.button('Restore source population',key='restore_source_population',on_click=restore_source)
    config['display']['sexes']=choice('Display sexes',['Separate female and male','Combined display'],config['display']['sexes'])
    first_fast_source=p['scope']=='National' and p['country']=='AUS' and p['source'] in {'Recommended frozen source','ABS ERP June2024'}
    if p['source']=='UN WPP2024 medium' or p['source']=='Recommended frozen source' and p['country']!='AUS':
        st.caption('First UN WPP validation can take about30seconds; it validates the complete frozen snapshot locally and is cached. Loading is explicit,not triggered by every edit.')
    if next(c for c in countries() if c['country_id']==p['country'])['completeness']=='unavailable':
        st.warning('Country incidence unavailable — PRK/missing incidence remains blocked. No value is invented.')
    if st.button('Load / check population',type='primary',key='load_population') or first_fast_source:
        with st.spinner('Checking frozen demographic records and coverage…'):
            context=ensure_context(config);st.session_state['review_population_summary']=context
    context=st.session_state.get('review_population')
    if context is not None:
        from engine.country_tb.data_common import canonical,digest
        if context['population_definition_hash']==digest(canonical(p)):
            st.subheader('Population coverage');a,b=st.columns(2);a.metric('Population',count(context['demography']['population']));b.metric('Open-age population retained',count(context['open_age_population']))
            original=context['source_profile'];st.caption('Source: '+('ABS ERP June2024' if original and original['provenance']['snapshot_id'].startswith('abs-') else 'UN World Population Prospects2024' if original else 'User-supplied population CSV'))
            if original:st.caption('Source status: '+', '.join(sorted({r['estimate_status']+' / '+r['projection_scenario'] for r in original['rows']})))
            target=context['observed_incidence']
            if target and target['central'] is not None:st.info(f"Observed WHO2024 incidence: {float(target['central']):,.2f} per100,000 person-years. Context only for custom/subgroups; not an infection rate.")
            else:st.warning('Country incidence unavailable — PRK/missing incidence remains blocked. No value is invented.')
            demographic_chart(context['demography']['rows'])
            st.caption('Female and male records remain separate internally. Published both-sex rounding differences are preserved in source reconciliation.')
            with st.expander('Technical source details'):st.json(dict(snapshot_identifiers=context['source_identifiers'],demography_hash=context['demography']['profile_hash'],
                published_both_reconciliation=None if original is None else original.get('published_both_checks')))
    st.subheader('Optional population subgroups')
    st.caption('Define a non-overlapping partition. Population counts determine weights. Birth,exposure and arrival fields are metadata,not a migration model.')
    if st.button('Add subgroup',key='add_subgroup'):
        i=1
        while any(g['subgroup_id']=='group_'+str(i) for g in p['subgroups']):i+=1
        p['subgroups'].append(dict(subgroup_id='group_'+str(i),label='New subgroup',population='',tbi_percent='',source_mode='User supplied',
            exposure_geography='',country_of_birth='',previous_residence='',arrival_year='',arrival_date='',residence_duration='',
            prior_tb_prevalence='',imported_active_contribution='',recurrence_contribution=''))
    if p['subgroups']:
        editor=st.data_editor(pd.DataFrame(p['subgroups']),hide_index=True,use_container_width=True,key='edit_subgroup_table',
            column_config={'source_mode':st.column_config.SelectboxColumn(options=['User supplied','Unresolved','Incidence inference — blocked'])})
        p['subgroups']=editor.fillna('').astype(str).to_dict('records')
        selected=st.selectbox('Subgroup to duplicate or remove',[g['subgroup_id'] for g in p['subgroups']],key='edit_selected_group')
        a,b=st.columns(2)
        def duplicate():
            g=deepcopy(next(g for g in p['subgroups'] if g['subgroup_id']==selected));base=g['subgroup_id']+'_copy';i=1;identifier=base
            while any(x['subgroup_id']==identifier for x in p['subgroups']):i+=1;identifier=base+'_'+str(i)
            g['subgroup_id']=identifier;g['population']='';p['subgroups'].append(g)
            st.session_state.pop('edit_subgroup_table',None)
        def remove():
            p['subgroups']=[g for g in p['subgroups'] if g['subgroup_id']!=selected];st.session_state.pop('edit_subgroup_table',None)
        a.button('Duplicate subgroup',key='duplicate_subgroup',on_click=duplicate)
        b.button('Remove subgroup',key='remove_subgroup',on_click=remove)
        p['partition_confirmed']=check('I declare these subgroup memberships non-overlapping',p['partition_confirmed'])
        p['subgroup_proportions']=check('Use parent age/sex proportions separately within every subgroup (diagnostic assumption)',p['subgroup_proportions'])
        if not p['subgroup_proportions']:p['subgroup_csv']=csv_input('Subgroup population',p['subgroup_csv'],'subgroup_population')
        st.caption('tbi_percent is in%; empty values never mean zero. Unresolved or incidence-selected groups prevent a complete mixture result.')
