"""Explicit baseline-TBI source selection; incidence inference never runs here."""
from copy import deepcopy
from decimal import Decimal

from engine.country_tb.calibration_common import reviewed, scientific_context, probability
from engine.country_tb.composition_common import (
    artifact,verify,metadata,require,identity,validate_demography,as_catalytic_demography,
    profile_summary,uncertainty,counts,STATES,
)
from engine.country_tb.data_common import canonical,decimal,numeric_text,digest
from engine.country_tb.active_tb_targets import TARGETS,BLOCKED,DEPENDENCIES,validate_active_tb_targets
from engine.country_tb.direct_tbi_profile import (
    DIRECT,UNAVAILABLE_RECENCY,CONDITIONAL_RECENCY,allocate_direct_profile,
)

SOURCE='country_tb.baseline_tbi_source.v1'
RESULT='country_tb.baseline_tbi_resolution.v1'
INCIDENCE='active_tb_incidence_inference_v1'
OBSERVED='observed_or_user_tbi_prevalence_v1'
OVERALL='overall_prevalence_anchor_v1'
AGE_ANCHORS='age_specific_prevalence_anchors_v1'
DIRECT_TOTAL='direct_age_specific_profile_v1'
DIRECT_STATES='direct_recent_remote_profile_v1'
MODE_TEXT={
    INCIDENCE:{'label':'Estimate from active TB incidence',
        'description':'Infer the baseline TBI profile using active TB incidence, demographic structure, infection-pressure history, progression natural history and population risk factors.'},
    OBSERVED:{'label':'Use observed or user-specified TBI prevalence',
        'description':'Use a reviewed TBI prevalence estimate or profile for this population or subgroup. Where appropriate, the model can use the historical infection-pressure shape to distribute an overall prevalence across age and recent/remote infection states.'}}


def baseline_tbi_source(*,mode,geography_id,geography_level,subgroup_id,source,provenance,
                        review,profile_kind,configuration,baseline_year=2024):
    require(isinstance(mode,str) and mode in MODE_TEXT,'Exactly one known primary baseline-TBI source mode is required.')
    identity(geography_id,geography_level,subgroup_id,baseline_year)
    metadata(source,provenance,review,required_review=False)
    require(profile_kind in {'empirical','user_defined','model_inferred','diagnostic'},'Explicit profile classification required.')
    require(isinstance(configuration,dict),'Explicit mode-specific configuration required.')
    config=deepcopy(configuration)
    if mode==INCIDENCE:
        require(set(config)=={'targets','pressure_history','forward_operator'},'Incidence mode requires targets, explicit nullable history and pending operator; no prevalence fallback.')
        validate_active_tb_targets(config['targets'])
        require(config['forward_operator'] is None,'No forward operator is approved in this milestone; supplying one cannot enable inference.')
        for row in config['targets']['rows']:
            require(row['geography_id']==geography_id and row['geography_level']==geography_level,
                    'Active-TB targets must match current analysis geography and level.')
        if config['pressure_history'] is not None:
            from engine.country_tb.pressure_shape import validate_history
            validate_history(config['pressure_history'])
        recency='unavailable'
    else:
        submode=config.get('submode')
        require(submode in {OVERALL,AGE_ANCHORS,DIRECT_TOTAL,DIRECT_STATES},'Unknown observed/user-input submode.')
        common={'submode','baseline_time','recent_window_years'}
        required=common|({'pressure_history','calibration_policy','anchor'} if submode==OVERALL else
            {'pressure_history','calibration_policy','observations'} if submode==AGE_ANCHORS else
            {'direct_profile','recency_policy','comparison_profile'})
        optional={'tail_distributions','age_policy','age_review','history_application_policy'} if submode in {OVERALL,AGE_ANCHORS} else set()
        require(required<=config.keys() and not set(config)-required-optional,
                'Incomplete or mixed baseline-TBI submode inputs; no silent fallback.')
        require(2024<=decimal(config['baseline_time'])<2025 and decimal(config['recent_window_years'])>0,
                'Baseline time must lie in 2024 and recent window must be positive.')
        if submode in {OVERALL,AGE_ANCHORS}:
            from engine.country_tb.pressure_shape import validate_history
            from engine.country_tb.calibration_common import verify as verify_m2c
            from engine.country_tb.catalytic_calibration import POLICY
            validate_history(config['pressure_history'])
            verify_m2c(config['calibration_policy'],POLICY)
            if submode==OVERALL:
                require(config['calibration_policy']['method']=='overall_anchor','Overall TBI anchor requires one-scale bounded root policy.')
                anchor=config['anchor']
                require(isinstance(anchor,dict) and {'prevalence','uncertainty','numerator','denominator','year','source','review',
                    'assay_operator','population_denominator','population_or_subgroup_id','coverage_policy','provenance'}<=anchor.keys(),
                    'Incomplete overall anchor contract.')
                metadata(anchor['source'],anchor['provenance'],anchor['review'])
                require(type(anchor['year']) is int and anchor['year']==2024,'Overall helper requires 2024 weights and anchor year; no automatic redating.')
                require(anchor['population_or_subgroup_id']==subgroup_id and
                    anchor['coverage_policy']=='apply_anchor_to_model_covered_population_v1','Overall anchor must explicitly identify subgroup and covered-population policy.')
                with scientific_context():
                    p=probability(anchor['prevalence']);counts(anchor,p);uncertainty(anchor['uncertainty'],p)
            else:
                from engine.country_tb.tbi_observations import OBSERVATIONS
                verify_m2c(config['observations'],OBSERVATIONS)
                require(all(r.get('subgroup_id')==subgroup_id for r in config['observations']['rows']),
                        'Age-specific observations must explicitly identify this subgroup.')
            recency='model_derived'
        else:
            verify(config['direct_profile'],DIRECT)
            require(config['recency_policy'] in {UNAVAILABLE_RECENCY,CONDITIONAL_RECENCY},'Unknown recency allocation policy.')
            leaves=[leaf for r in config['direct_profile']['rows'] for leaf in r.get('within_band_values',[r])]
            if submode==DIRECT_STATES:
                require(all(r.get('recent') is not None for r in leaves),'Direct recent/remote submode requires all exclusive states.')
            if all(r.get('recent') is not None for r in leaves):
                recency='observed' if profile_kind=='empirical' else 'user_defined'
            else:
                recency='model_derived' if config['recency_policy']==CONDITIONAL_RECENCY else 'unavailable'
    return artifact(SOURCE,mode=mode,geography_id=geography_id,geography_level=geography_level,
        subgroup_id=subgroup_id,baseline_year=baseline_year,source=source,provenance=provenance,review=review,
        provenance_hash=digest(canonical({'source':source,'provenance':provenance,'review':review})),
        profile_kind=profile_kind,recent_remote_composition=recency,configuration=config,
        status='selected_pending_resolution',usable_for_downstream_simulation=False,
        fallback_policy='forbidden',display=MODE_TEXT[mode])


def _history_for_population(selection):
    """Reviewed transfer of a fixed proxy; no arrival-time or migration model."""
    history=selection['configuration']['pressure_history']
    if history['geography_id']==selection['geography_id']:
        return history
    policy=selection['configuration'].get('history_application_policy')
    require(isinstance(policy,dict) and policy.get('mode')=='reviewed_subgroup_exposure_proxy_v1',
            'Exposure-history geography differs from residence; an explicit reviewed subgroup proxy application is required.')
    metadata(policy.get('source'),policy.get('provenance'),policy.get('review'))
    require(policy.get('rationale'),'Explain exposure-proxy transport; current residence does not imply lifetime exposure geography.')
    from engine.country_tb.pressure_shape import supplied_shape,apply_backcast
    # Preserve all original curves/history in provenance. This new shape is an
    # explicitly reviewed scenario for the receiving subgroup, not a relabelled
    # source incidence series or a reconstruction of migration histories.
    curves=history['curves']
    q=supplied_shape(curves['common'],geography_id=selection['geography_id'],review=policy['review'],
        provenance={'source_pressure_history':deepcopy(history),'source_pressure_history_hash':history['profile_hash'],
                    'subgroup_id':selection['subgroup_id'],'application_policy':deepcopy(policy)},
        normalization_time='2024',sex_rows={k:v for k,v in curves.items() if k!='common'})
    return apply_backcast(q,start_time=curves['common'][0]['start'])


def _overall_observation(selection,demography):
    from engine.country_tb.tbi_observations import observations_profile
    c=selection['configuration'];a=c['anchor']
    with scientific_context():
        population=sum(decimal(r['population']) for r in demography['rows'])
        require(decimal(a['population_denominator'])==population,'Overall anchor population denominator does not match subgroup demography.')
        strata=[]
        tails=c.get('tail_distributions') or {}
        for r in demography['rows']:
            pop=decimal(r['population'])
            if pop==0:continue
            if r['age_upper_exclusive'] is not None:
                strata.append({'age_lower':str(r['age_lower']),'age_upper_exclusive':str(r['age_upper_exclusive']),
                    'age_policy':c.get('age_policy','uniform_within_single_year'),'age_review':c.get('age_review'),
                    'sex':r['sex'],'weight':numeric_text(pop)})
            elif (tail:=tails.get(r['sex']+':'+str(r['age_lower']))) is not None:
                reviewed(tail.get('review'));require(tail.get('provenance') and tail.get('ages'),'Reviewed tail ages and provenance required.')
                total=sum(decimal(x['weight']) for x in tail['ages'])
                require(total>0,'Tail weights must have positive mass.')
                for x in tail['ages']:
                    require(decimal(x['age'])>=r['age_lower'] and decimal(x['weight'])>0,'Invalid terminal-age weighting.')
                    strata.append({'age_lower':numeric_text(x['age']),'age_upper_exclusive':None,'age_policy':'exact_age',
                                   'sex':r['sex'],'weight':numeric_text(pop*decimal(x['weight'])/total)})
        require(bool(strata),'No age-covered population is available for overall calibration.')
        sexes={r['sex'] for r in strata}
        hi=None if any(r['age_upper_exclusive'] is None for r in strata) else numeric_text(max(decimal(r['age_upper_exclusive']) for r in strata))
        row={**deepcopy(a),'observation_id':a.get('observation_id','overall:'+selection['subgroup_id']),
            'geography_id':selection['geography_id'],'geography_level':selection['geography_level'],
            'subgroup_id':selection['subgroup_id'],'survey_year':a['year'],'survey_time':c['baseline_time'],
            'age_lower':numeric_text(min(decimal(r['age_lower']) for r in strata)),'age_upper_exclusive':hi,
            'sex':next(iter(sexes)) if len(sexes)==1 else 'both','sampling_frame':a.get('sampling_frame','explicit supplied subgroup'),
            'survey_design':a.get('survey_design',{'kind':'reviewed_prevalence_anchor'}),
            'observation_policy':a['coverage_policy'],
            'weight_source':{'demography_hash':demography['profile_hash'],'coverage_policy':a['coverage_policy']},'strata':strata}
        return observations_profile([row],review=a['review'])


def resolve_baseline_tbi(selection,demography):
    verify(selection,SOURCE);validate_demography(demography)
    # Revalidate dispatch and mode-specific contracts even if a caller re-seals
    # a hand-edited selection. Hashes identify content, not scientific approval.
    rebuilt=baseline_tbi_source(**{k:selection[k] for k in ('mode','geography_id','geography_level','subgroup_id',
        'source','provenance','review','profile_kind','configuration','baseline_year')})
    require(rebuilt['profile_hash']==selection['profile_hash'],'Source selection does not reproduce from validated inputs.')
    require(demography['geography_id']==selection['geography_id'] and demography['baseline_year']==2024,
            'Source selection and population must share residence geography and baseline.')
    levels={r.get('geography_level',demography.get('geography_level')) for r in demography['rows']}
    require(levels=={selection['geography_level']},'Source and population geography levels must agree.')
    with scientific_context():
        rows=[{**{k:r[k] for k in ('age_lower','age_upper_exclusive','sex','population')},
            **{k:None for k in STATES},'most_recent_infection_time':None,
            'recent_time_summary':None,'remote_only_time_summary':None,'status':'unresolved'} for r in demography['rows']]
        model=None;calibration=None;gaps=[];warnings=[];dependencies=[]
        if selection['mode']==INCIDENCE:
            status=BLOCKED;dependencies=list(DEPENDENCIES)
            warnings=['active_TB_cannot_identify_TBI_without_an_approved_progression_and_population_risk_forward_operator']
        elif selection['review']['status']!='reviewed':
            status='unavailable_pending_source_review';warnings=['selected_source_is_not_reviewed_no_fallback']
        else:
            config=selection['configuration'];submode=config['submode']
            if submode in {OVERALL,AGE_ANCHORS}:
                from engine.country_tb.catalytic_calibration import calibrate_scale
                from engine.country_tb.tbi_baseline import baseline_tbi_profile
                history=_history_for_population(selection)
                observations=_overall_observation(selection,demography) if submode==OVERALL else config['observations']
                if submode==OVERALL:
                    require(set(config['calibration_policy']['observation_weights'])=={observations['rows'][0]['observation_id']},
                            'Overall policy weights must identify the derived overall observation ID.')
                require(all(r['geography_id']==selection['geography_id'] and r['geography_level']==selection['geography_level'] for r in observations['rows']),
                        'TBI observations must match analysis residence geography and level.')
                calibration=calibrate_scale(history,config['calibration_policy'],observations=observations)
                status=calibration['fit_status']
                if status=='accepted':
                    model=baseline_tbi_profile(as_catalytic_demography(demography),calibration,
                        recent_window=config['recent_window_years'],baseline_time=config['baseline_time'],
                        tail_distributions=config.get('tail_distributions'),age_policy=config.get('age_policy','uniform_within_single_year'),
                        age_review=config.get('age_review'),scenario_id='source:'+selection['subgroup_id'])
                    rows=[deepcopy(r) for r in model['rows']]
                    for r in rows:r['recency_composition']='model_derived' if r['any_tbi'] is not None else 'unavailable'
                    warnings=list(model['warnings'])+['age_and_recency_distributions_are_conditional_model_derivations']
                else:warnings=['one_common_scale_cannot_satisfy_declared_anchor_compatibility_no_fallback']
            else:
                model=config['comparison_profile']
                if model is not None:
                    require(decimal(model['baseline_time'])==decimal(config['baseline_time']) and
                        decimal(model['recent_window_years'])==decimal(config['recent_window_years']),
                        'Comparison model baseline time and recent window must match direct inputs.')
                rows,gaps=allocate_direct_profile(config['direct_profile'],demography,
                    recency_policy=config['recency_policy'],model_profile=model)
                if selection['profile_kind']=='empirical':
                    for row in rows:
                        if row['recency_composition']=='user_defined':row['recency_composition']='observed'
                status='direct_profile_resolved' if not gaps else 'partial_direct_profile_coverage'
                if model is None:warnings.append('model_comparison_unavailable_no_model_profile_supplied')
        summary=profile_summary(rows)
        kind='diagnostic' if selection['profile_kind']=='diagnostic' else (
            'model_inferred' if selection['mode']==OBSERVED and selection['configuration']['submode'] in {OVERALL,AGE_ANCHORS}
            else selection['profile_kind'])
        available={r.get('recency_composition','unavailable') for r in rows if decimal(r['population'])>0}
        composition=next(iter(available)) if len(available)==1 else 'mixed_or_partial'
        usable=(summary['total_tbi_complete'] and summary['recency_complete'] and kind!='diagnostic' and
                selection['review']['status']=='reviewed' and selection['mode']==OBSERVED)
        return artifact(RESULT,mode=selection['mode'],geography_id=selection['geography_id'],geography_level=selection['geography_level'],
            subgroup_id=selection['subgroup_id'],baseline_year=2024,selection=selection,selection_hash=selection['profile_hash'],
            demography=demography,demography_hash=demography['profile_hash'],status=status,profile_kind=kind,
            evidence_kind=selection['profile_kind'],recent_remote_composition=composition,rows=sorted(rows,key=canonical),
            summary=summary,model_derived_profile=model,calibration=calibration,direct_profile_gaps=gaps,
            missing_approved_components=dependencies,tbi_estimate=None if selection['mode']==INCIDENCE else summary['components']['any_tbi']['prevalence_full'],
            usable_for_downstream_simulation=usable,downstream_scope='baseline_exclusive_state_initialization_only_no_progression_or_interventions',
            warnings=sorted(set(warnings)),source=selection['source'],provenance=selection['provenance'],review=selection['review'],
            interpretation='conditional_on_selected_source_history_backcast_reinfection_demography_and_explicit_exposure_proxy_policy',
            migration_policy='metadata_only_no_migration_model',progression_forward_operator=None)


def interface_specification():
    """Pure presentation state; no Streamlit import or operational inference UI."""
    return artifact('country_tb.baseline_tbi_interface_specification.v1',
        question='How should baseline TBI prevalence be determined?',
        options=[{'mode':mode,**MODE_TEXT[mode],'operational':mode==OBSERVED,
            'blocked_status':BLOCKED if mode==INCIDENCE else None} for mode in (INCIDENCE,OBSERVED)],
        user_subchoices=['Overall TBI prevalence','Age-specific TBI prevalence','Upload age/sex-specific profile',
                        'Enter recent and remote infection directly'],
        subgroup_actions=['Add subgroup','Duplicate subgroup','Remove subgroup','Edit subgroup',
                          'Show subgroup results separately','Show population-weighted combined result'],
        duplication_rule='new_stable_ID_and_reviewed_membership_required_before_mixture_inclusion',
        exposure='developer_specification_only_no_production_navigation')
