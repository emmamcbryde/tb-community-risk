import streamlit as st
from country_model_review_ui.state import SOURCE_MODES,TBI_COLUMNS,csv_text,working_example,restore,material_hash
from country_model_review_ui.model import build_baseline,geography
from country_model_review_ui.presentation import text,choice,check,csv_input,ensure_context,tbi_table,percent,count


def render(config):
    st.header('Baseline TBI');t=config['tbi']
    t['mode']=st.radio('How should baseline TBI prevalence be determined?',SOURCE_MODES,index=SOURCE_MODES.index(t['mode']),key='edit_TBI_mode')
    if t['mode']==SOURCE_MODES[1]:
        st.warning('Experimental diagnostic — requires reviewed natural-history assumptions')
        st.error('Calculation blocked by default. The public incidence source remains unavailable_pending_progression_and_risk_forward_operator.')
        st.caption('No fallback to user values,fictional country results,synthetic parameters or legacy calibration. Use Natural history and scientific gates to inspect requirements.')
        return
    st.info('Recommended operational pathway: observed or user-specified TBI. Source values remain separate from model-derived comparisons.')
    t['input_kind']=choice('TBI input format',['Overall prevalence','Age-specific total profile','Age-by-sex total profile','Direct recent and remote'],t['input_kind'])
    t['source']=text('TBI source / citation',t['source']);t['review_status']=choice('TBI scientific review status',['not_reviewed','reviewed','diagnostic'],t['review_status'])
    if t['input_kind']=='Overall prevalence':
        t['overall_percent']=text('Overall prevalence (%)',t['overall_percent'])
        t['derive']=check('Derive age and recent/remote composition using the catalytic model',t['derive'],
            help='Explicit opt-in. One common scale is fitted to covered demography; pressure shape,backcast and reinfection assumptions condition the result.100+ remains unresolved.')
        if not t['derive']:st.caption('Total-only: applies the supplied total explicitly across the selected population. Recent/remote remains unavailable unless explicitly allocated below.')
    elif t['input_kind']=='Direct recent and remote':
        a,b,c=st.columns(3)
        with a:t['recent_percent']=text('Recent (%)',t['recent_percent'])
        with b:t['remote_percent']=text('Remote-only (%)',t['remote_percent'])
        with c:t['prior_percent']=text('Prior remote + recent subset (%) — optional',t['prior_percent'])
        st.caption('Never infected is1−recent−remote-only. Prior remote + recent is within recent and is never added again. Use CSV for age/sex-specific direct states.')
        t['csv']=csv_input('TBI profile',t['csv'],'tbi')
    else:
        t['csv']=csv_input('TBI profile',t['csv'],'tbi')
        st.caption('Direct profiles preserve supplied proportions exactly. Gaps remain unresolved; overlaps fail. No silent interpolation.')
    t['assay']=choice('Input represents',['Direct latent TBI','QFT-Plus','QFT-GIT','TST'],t['assay'])
    if t['assay']!='Direct latent TBI':
        st.warning('Observed test positivity is not latent TBI. Supply fixed,product/threshold-specific assay performance and provenance; it is not fitted with infection pressure.')
        t['sensitivity_percent']=text('Assay sensitivity (%)',t['sensitivity_percent']);t['specificity_percent']=text('Assay specificity (%)',t['specificity_percent'])
        t['assay_source']=text('Assay performance source',t['assay_source']);t['threshold']=text('Assay positivity threshold',t['threshold'])
        st.caption('Supplied performance is applied uniformly across included ages as a diagnostic assumption. Pediatric,BCG and immunosuppression transport remain unapproved.')
        if t['input_kind']!='Overall prevalence':st.error('Age/sex direct CSV values must be latent TBI. Convert with a reviewed observation model before uploading; no automatic mixed-assay conversion.')
    if not (t['input_kind']=='Overall prevalence' and t['derive']):
        t['hybrid_recency']=check('Explicit hybrid: allocate total-only TBI using model-conditional recent/remote composition',t['hybrid_recency'])
        t['conditional_clock']=check('Use model-conditional infection-time distributions for direct states (diagnostic only)',t['conditional_clock'])
    if t['derive'] and t['input_kind']=='Overall prevalence' or t['hybrid_recency'] or t['conditional_clock']:
        t['pressure']=choice('Pressure-history shape',['WHO historical shape','Explicit constant-shape diagnostic'],t['pressure'])
        t['proxy_confirmed']=check('Explicitly use this pressure history as a diagnostic exposure proxy',t['proxy_confirmed'])
        st.caption('WHO shape is a smoothed incidence proxy,not an estimated force of infection. Pre2000 uses the earliest normalized value. Residence does not determine historical exposure; no post2024 extrapolation.')
    if st.button('Build / check TBI profile',type='primary',key='build_TBI'):
        context=ensure_context(config);value=build_baseline(config,context)
        st.session_state['review_TBI_preview']={'configuration_hash':material_hash(config),'baseline':value}
    preview=st.session_state.get('review_TBI_preview')
    if preview and preview['configuration_hash']==material_hash(config):
        baseline=preview['baseline'];s=baseline['summary']
        if s['total_tbi_complete']:st.success('TBI profile constructed for internal review; scientific approval remains pending.')
        else:st.warning('TBI profile has unresolved coverage. The complete population denominator is retained; scientific approval remains pending.')
        a,b=st.columns(2);a.metric('Total TBI — full population',percent(s['components']['any_tbi']['prevalence_full']));b.metric('Unresolved TBI population',count(s['unresolved_population']))
        st.dataframe(tbi_table(baseline['rows']),hide_index=True,use_container_width=True)
    with st.expander('CSV schema and template'):
        st.code(','.join(TBI_COLUMNS),language='text')
        template=dict(record_id='adult_profile',subgroup_id='main',age_lower=0,age_upper_exclusive='',sex='both',any_tbi='',recent='',remote_only='',never_infected='',prior_remote_plus_recent='',
            allocation_policy='piecewise_constant_within_band_v1',source='Enter source',review_status='not_reviewed',geography_id=geography(config),baseline_year=2024)
        st.download_button('Download blank TBI CSV template',csv_text([template],TBI_COLUMNS),'TBI_template.csv','text/csv')
        st.caption('CSV probabilities use proportions0–1. Percent input fields above use0–100. Missing values stay blank. Explicit100+ direct values are permitted without inventing an age.')
    st.button('Restore working example',key='restore_working_example_baseline',on_click=lambda:restore(st.session_state,working_example()))
