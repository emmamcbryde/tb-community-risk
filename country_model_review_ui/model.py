"""Thin offline consumers of retained2A–3D contracts; no production dispatch."""
from copy import deepcopy
from decimal import Decimal,localcontext
from functools import lru_cache
from time import perf_counter
import json

from engine.country_tb.data_common import canonical,digest,seal,decimal,numeric_text,DataValidationError
from engine.country_tb.composition_common import population_demography,profile_summary
from engine.country_tb.incidence import load_who_snapshot,WHO_2025,country_index,national_profile
from engine.country_tb.demography import load_demography_snapshot,demographic_profile
from engine.country_tb.pressure_shape import shape_from_incidence,supplied_shape,apply_backcast
from engine.country_tb.catalytic_calibration import calibration_policy
from engine.country_tb.tbi_observations import assay_operator
from engine.country_tb.baseline_tbi_source import baseline_tbi_source,resolve_baseline_tbi,OBSERVED,OVERALL,DIRECT_TOTAL,DIRECT_STATES
from engine.country_tb.direct_tbi_profile import direct_tbi_profile,CONDITIONAL_RECENCY,UNAVAILABLE_RECENCY
from engine.country_tb.population_subgroups import population_subgroup,population_mixture
from engine.country_tb.risk_factor_catalog import progression_risk_factor_catalog,catalogue_records
from engine.country_tb.risk_prevalence import risk_factor_prevalence_profile,CONSTANT,UNIFORM
from engine.country_tb.risk_policies import overlap_policy,association_policy,multimorbidity_policy
from engine.country_tb.risk_allocation import allocate_baseline_risk,combine_risk_allocations
from engine.country_tb.risk_common import ratio
from country_model_review_ui.state import require,parse_csv,percent_input,material_hash,SOURCE_MODES
from country_model_review_ui import APP_VERSION

REVIEW={'status':'reviewed','source_id':'review_ui_internal_contract_processing',
    'rationale':'Internal source transcription and explicitly requested diagnostic processing only; no external scientific approval.'}
SOURCE={'source_id':'reviewer_entry','citation':'Reviewer inputs; original source and scientific review status retained separately'}
PROVENANCE={'application_version':APP_VERSION,'scientific_approval':'pending','production_available':False}
META=dict(source=SOURCE,provenance=PROVENANCE,review=REVIEW)


@lru_cache(maxsize=1)
def who():return load_who_snapshot(WHO_2025)


@lru_cache(maxsize=1)
def countries():return country_index(who())


@lru_cache(maxsize=12)
def incidence(country):return national_profile(who(),country)


@lru_cache(maxsize=2)
def demographic_snapshot(identifier):return load_demography_snapshot(identifier)


@lru_cache(maxsize=12)
def source_demography(country,source):
    require(source in {'Recommended frozen source','ABS ERP June2024','UN WPP2024 medium'},'Unknown demographic source; no fallback.')
    identifier='abs-erp-jun2024-07ca11bcb092' if country=='AUS' and source in {'Recommended frozen source','ABS ERP June2024'} else 'un-wpp2024-medium-8028a7a3e163'
    require(source!='ABS ERP June2024' or country=='AUS','ABS is not a demographic fallback for other countries.')
    return demographic_profile('who:'+country,[demographic_snapshot(identifier)],year=2024,preferred_snapshot_id=identifier)


@lru_cache(maxsize=1)
def catalog():return progression_risk_factor_catalog()


def geography(config):
    p=config['population'];return 'who:'+p['country'] if p['scope']=='National' else p['custom_id']


def scaled_rows(rows,population):
    """Only explicit user scaling, retaining source counts independently."""
    with localcontext() as ctx:
        ctx.prec=80;target=decimal(population);total=sum(decimal(r['population']) for r in rows)
        require(target>0 and total>0,'Positive population required.')
        values=[{k:r[k] for k in ('age_lower','age_upper_exclusive','sex','population')} for r in rows]
        for r in values:r['population']=numeric_text(target*decimal(r['population'])/total)
        correction=target-sum(decimal(r['population']) for r in values)
        values[-1]['population']=numeric_text(decimal(values[-1]['population'])+correction)
    return values,numeric_text(correction)


def select_population(config):
    p=config['population'];geo=geography(config);require(p['baseline_year']==2024,'Only frozen2024 contracts are supported.')
    require(p['scope'] in {'National','Custom / subnational'},'Unknown population scope; no fallback.')
    require(geo and (p['scope']=='National' or geo.startswith(('custom:','subnational:'))),'Custom geography requires custom: or subnational: namespace.')
    original=None;receipt=None;correction='0'
    if p['source']=='User CSV':
        receipt=parse_csv(p['csv'],'population',geography_id=geo);rows=receipt['rows']
        with localcontext() as ctx:
            ctx.prec=80;total=sum(decimal(r['population']) for r in rows)
        require(not p['size'] or total==decimal(p['size']),'CSV population disagrees with the declared denominator.')
    else:
        original=source_demography(p['country'],p['source']);rows=original['rows'];total=decimal(original['combined_population'])
        require(p['scope']=='National' or p['national_proxy'],'Explicitly select national age/sex proportions as a diagnostic proxy,or supply subgroup data.')
        if p['scope']!='National' or p['scale_source']:
            require(p['size'],'Enter the population denominator for explicit source scaling.')
            rows,correction=scaled_rows(rows,p['size']);total=decimal(p['size'])
        else:rows=[{k:r[k] for k in ('age_lower','age_upper_exclusive','sex','population')} for r in rows]
    proof={'source_profile_hash':None if original is None else original['profile_hash'],'original_source':original,
        'CSV_receipt':receipt,'national_proxy_explicitly_selected':p['national_proxy'],
        'scaling_roundoff_adjustment_persons':correction,'scaling_scope':'user_defined_population_only_source_rounding_never_forced_to_zero'}
    level='national' if p['scope']=='National' else 'subnational' if geo.startswith('subnational:') else 'custom'
    demo=population_demography(rows,geography_id=geo,geography_level=level,
        population=numeric_text(total),source=SOURCE,provenance={**PROVENANCE,**proof},review=REVIEW)
    target=next((r for r in incidence(p['country'])['rows'] if r['year']==2024),None)
    return dict(demography=demo,source_profile=original,incidence=incidence(p['country']),observed_incidence=target,
        source_identifiers={'WHO':WHO_2025,'demography':None if original is None else original['provenance']['snapshot_id']},
        open_age_population=numeric_text(sum(decimal(r['population']) for r in demo['rows'] if r['age_upper_exclusive'] is None)),
        national_incidence_applies_to_custom=False,population_definition=deepcopy(p),population_definition_hash=digest(canonical(p)))


def _groups(config,context):
    p=config['population'];groups=deepcopy(p['subgroups'])
    if not groups:groups=[dict(subgroup_id='main',label=p['label'],population=context['demography']['population'],tbi_percent='',source_mode='User supplied',
        exposure_geography='',arrival_year='',prior_tb_prevalence='',imported_active_contribution='',recurrence_contribution='')]
    require(len(groups)<=12 and len({g['subgroup_id'] for g in groups})==len(groups),'At most12 uniquely identified non-overlapping subgroups.')
    if len(groups)>1:require(p['partition_confirmed'],'Explicit non-overlapping membership declaration required.')
    with localcontext() as ctx:
        ctx.prec=80
        require(sum(decimal(g['population']) for g in groups)==decimal(context['demography']['population']),'Subgroup denominators must sum to the full population.')
    return groups


def history(config,context,geo):
    t=config['tbi'];require(t['proxy_confirmed'],'Explicit pressure-shape/exposure-history selection required for model-derived distributions.')
    finite=[r['age_upper_exclusive'] for r in context['demography']['rows'] if r['age_upper_exclusive'] is not None]
    start=numeric_text(min(Decimal('1900'),Decimal('2024')-Decimal(max(finite,default=100))-1))
    if t['pressure']=='Explicit constant-shape diagnostic':
        q=supplied_shape([{'start':'1900','end':'2025','shape':'1'}],geography_id=geo,normalization_time='2024',
            review=REVIEW,provenance={**PROVENANCE,'warning':'explicit_synthetic_shape_not_country_exposure_history'})
        return apply_backcast(q,start_time=start)
    q=shape_from_incidence(context['incidence'],review=REVIEW)
    if q['geography_id']!=geo:
        original=q
        # An explicit diagnostic transfer,never a residence=exposure assumption.
        curves=q['curves'] if 'curves' in q else None
        from engine.country_tb.catalytic_core import history_curve
        h=apply_backcast(q,start_time=start)
        q=supplied_shape(history_curve(h,'female'),geography_id=geo,normalization_time='2024',review=REVIEW,
            provenance={**PROVENANCE,'parent_shape_hash':original['profile_hash'],'explicit_subgroup_pressure_proxy':True,
                'historical_exposure_geography_not_identified':True})
    return apply_backcast(q,start_time=start)


def observation_operator(config):
    t=config['tbi'];direct=t['assay']=='Direct latent TBI'
    require(direct or t['assay_source'] and t['threshold'],'Assay sensitivity/specificity require a source and threshold.')
    se='1' if direct else percent_input(t['sensitivity_percent']);sp='1' if direct else percent_input(t['specificity_percent'])
    require(direct or decimal(se)+decimal(sp)>1,'An informative assay operator requires Se+Sp>1; no clipping.')
    return assay_operator(assay_class='direct_latent_tbi' if direct else 'TST' if t['assay']=='TST' else 'IGRA',name=t['assay'],
        positivity_threshold='not_applicable' if direct else t['threshold'],sensitivity=se,specificity=sp,review=REVIEW,
        provenance={**PROVENANCE,'source':t['assay_source'] or t['source'],'input_review_status':t['review_status'],
            'performance_fixed_not_fitted':True,'population_transport_not_approved':True})


def build_baseline(config,context):
    require(context['population_definition_hash']==digest(canonical(config['population'])),'Population context is stale; reselect population before building TBI.')
    require(config['tbi']['mode']==SOURCE_MODES[0],'Incidence mode cannot fall back to supplied TBI.')
    t=config['tbi'];groups=_groups(config,context);geo=context['demography']['geography_id'];level=context['demography']['geography_level']
    require(t['input_kind']=='Overall prevalence' or t['assay']=='Direct latent TBI',
        'Direct age/sex/state profiles must be latent TBI,not untransformed assay positivity.')
    ids=[g['subgroup_id'] for g in groups];upload=parse_csv(t['csv'],'tbi',geography_id=geo,subgroup_ids=ids) if t['csv'] else None
    pop_upload=parse_csv(config['population']['subgroup_csv'],'subgroup_population',geography_id=geo,subgroup_ids=ids) if config['population']['subgroup_csv'] else None
    subgroup_results=[];missing=[];all_rows=[]
    for g in groups:
        if pop_upload:
            rows=[r for r in pop_upload['rows'] if r['subgroup_id']==g['subgroup_id']]
            require(rows,'Every subgroup needs population rows.')
            demo=population_demography(rows,geography_id=geo,geography_level=level,population=g['population'],
                source=SOURCE,provenance={**PROVENANCE,'CSV_receipt':pop_upload},review=REVIEW)
        elif len(groups)==1:demo=context['demography']
        else:
            require(config['population']['subgroup_proportions'],'Explicit subgroup age/sex distributions or selected parent-proportion policy required.')
            rows,adjustment=scaled_rows(context['demography']['rows'],g['population'])
            demo=population_demography(rows,geography_id=geo,geography_level=level,population=g['population'],
                source=SOURCE,review=REVIEW,provenance={**PROVENANCE,'parent_demography_hash':context['demography']['profile_hash'],
                    'age_sex_policy':'explicit_parent_proportions_diagnostic','roundoff_adjustment':adjustment})
        comparison=None;selection=None
        supplied_p=g.get('tbi_percent') if len(groups)>1 else g.get('tbi_percent') or t['overall_percent']
        missing_input=(t['input_kind']=='Overall prevalence' and not str(supplied_p or '').strip()) or (
            t['input_kind'] in {'Age-specific total profile','Age-by-sex total profile'} and
            (upload is None or not any(r['subgroup_id']==g['subgroup_id'] for r in upload['rows'])))
        if g.get('source_mode') in {'Unresolved','Incidence inference — blocked'} or missing_input:
            missing.append(dict(subgroup_id=g['subgroup_id'],population=demo['population'],source_mode=g['source_mode'],
                status='unavailable_pending_progression_and_risk_forward_operator' if g['source_mode']=='Incidence inference — blocked' else 'unresolved_user_profile',
                original_input=g,demography=demo))
            all_rows.extend({**r,'subgroup_id':g['subgroup_id'],**{k:None for k in ('any_tbi','recent','remote_only','never_infected','prior_remote_plus_recent')}} for r in demo['rows'])
            continue
        else:
            p=percent_input(supplied_p,required=t['input_kind']=='Overall prevalence')
            if t['input_kind']=='Overall prevalence' and t['derive']:
                identifier='review:'+g['subgroup_id'];h=history(config,context,geo)
                anchor=dict(observation_id=identifier,population_or_subgroup_id=g['subgroup_id'],population_denominator=demo['population'],prevalence=p,
                    numerator=None,denominator=None,year=2024,uncertainty={'kind':'unspecified'},assay_operator=observation_operator(config),
                    coverage_policy='apply_anchor_to_model_covered_population_v1',**META)
                selection=baseline_tbi_source(mode=OBSERVED,geography_id=geo,geography_level=level,subgroup_id=g['subgroup_id'],
                    profile_kind='diagnostic',**META,configuration=dict(submode=OVERALL,baseline_time='2024.5',recent_window_years='5',
                    pressure_history=h,calibration_policy=calibration_policy(method='overall_anchor',alpha_upper='2',observation_weights={identifier:'1'},
                    residual_tolerance='.000001',review=REVIEW),anchor=anchor))
                comparison=resolve_baseline_tbi(selection,demo)
                if t['input_kind']=='Overall prevalence':
                    resolved=comparison
            else:
                if upload:
                    records=[{**r,'age_allocation_policy':r['allocation_policy']} for r in upload['rows'] if r['subgroup_id']==g['subgroup_id']]
                elif t['input_kind']=='Direct recent and remote':
                    recent=percent_input(t['recent_percent']);remote=percent_input(t['remote_percent'])
                    with localcontext() as ctx:
                        ctx.prec=80;total=decimal(recent)+decimal(remote);require(total<=1,'Recent + remote-only exceeds100%; no clipping.')
                        never=numeric_text(1-total)
                    records=[dict(record_id='states',age_lower=0,age_upper_exclusive=None,sex='both',
                        age_allocation_policy='piecewise_constant_within_band_v1',recent=recent,remote_only=remote,
                        never_infected=never,any_tbi=numeric_text(total),prior_remote_plus_recent=percent_input(t['prior_percent'],required=False))]
                else:
                    require(t['input_kind']=='Overall prevalence','Upload or paste the age/sex profile CSV before running.')
                    if t['assay']!='Direct latent TBI':
                        op=observation_operator(config);se=decimal(op['sensitivity']);sp=decimal(op['specificity'])
                        with localcontext() as ctx:
                            ctx.prec=80;latent=(decimal(p)-(1-sp))/(se+sp-1)
                        require(0<=latent<=1,'Assay positivity is infeasible under supplied performance; no clipping.');p=numeric_text(latent)
                    records=[dict(record_id='total',age_lower=0,age_upper_exclusive=None,sex='both',
                        age_allocation_policy='piecewise_constant_within_band_v1',any_tbi=p)]
            if selection is not None:
                pass
        if selection is None:
            require(records,'No TBI records for this subgroup; mark it explicitly unresolved instead of inventing values.')
            profile=direct_tbi_profile(records,**META)
            if t['hybrid_recency'] or t['conditional_clock']:
                known=[r for r in records if r.get('any_tbi') is not None]
                require(known,'A total-TBI profile is required for conditional model allocation.')
                # Fit only a population-weighted overall comparison; never alter supplied records.
                provisional=baseline_tbi_source(mode=OBSERVED,geography_id=geo,geography_level=level,subgroup_id=g['subgroup_id'],profile_kind='diagnostic',**META,
                    configuration=dict(submode=DIRECT_TOTAL,baseline_time='2024.5',recent_window_years='5',direct_profile=profile,
                        recency_policy=UNAVAILABLE_RECENCY,comparison_profile=None))
                temp=resolve_baseline_tbi(provisional,demo);mean=temp['summary']['components']['any_tbi']['prevalence_covered']
                require(mean is not None,'No covered total-TBI values for a conditional comparison.')
                cc=deepcopy(config);cc['tbi'].update(input_kind='Overall prevalence',derive=True,overall_percent=numeric_text(decimal(mean)*100),csv='',hybrid_recency=False,conditional_clock=False)
                cc['population'].update(subgroups=[])
                nested=build_baseline(cc,{**context,'demography':demo,'population_definition_hash':digest(canonical(cc['population']))})
                comparison=nested['subgroups'][0]['tbi_profile']
            selection=baseline_tbi_source(mode=OBSERVED,geography_id=geo,geography_level=level,subgroup_id=g['subgroup_id'],profile_kind='diagnostic',**META,
                configuration=dict(submode=DIRECT_STATES if t['input_kind']=='Direct recent and remote' else DIRECT_TOTAL,
                    baseline_time='2024.5',recent_window_years='5',direct_profile=profile,
                    recency_policy=CONDITIONAL_RECENCY if t['hybrid_recency'] else UNAVAILABLE_RECENCY,
                    comparison_profile=None if comparison is None else comparison['model_derived_profile']))
            resolved=resolve_baseline_tbi(selection,demo)
        all_rows.extend({**r,'subgroup_id':g['subgroup_id']} for r in resolved['rows'])
        membership=dict(partition_id='reviewer_declared_partition',category_id=g['subgroup_id'],definition=g['label'],person_ids=None,**META)
        subgroup_results.append(population_subgroup(subgroup_id=g['subgroup_id'],label=g['label'],description='Reviewer-defined population; no migration history inferred',
            current_residence_geography=geo,geography_level=level,population=demo['population'],demography=demo,tbi_source=selection,
            membership=membership,resolved_profile=resolved,exposure_history_geography=None if not g.get('exposure_geography') else
                {'geography_id':g['exposure_geography'],'geography_level':'national'},
            arrival_year=None if not g.get('arrival_year') else int(g['arrival_year']),
            arrival_date=g.get('arrival_date') or None,duration_in_current_setting_years=g.get('residence_duration') or None,
            country_or_region_of_birth=None if not g.get('country_of_birth') else {'geography_id':g['country_of_birth'],'geography_level':'national'},
            previous_residence=None if not g.get('previous_residence') else {'geography_id':g['previous_residence'],'geography_level':'national'},
            source=SOURCE,review=REVIEW,provenance={**PROVENANCE,'original_subgroup_input':g,'original_TBI_input':deepcopy(t),'CSV_receipt':upload}))
    reconciled={};reference={(r['sex'],r['age_lower'],r['age_upper_exclusive']):decimal(r['population']) for r in context['demography']['rows']}
    with localcontext() as ctx:
        ctx.prec=80
        for r in all_rows:
            key=r['sex'],r['age_lower'],r['age_upper_exclusive'];reconciled[key]=reconciled.get(key,Decimal(0))+decimal(r['population'])
        require(set(reconciled)==set(reference),'Subgroup age/sex coverage differs from the parent grid.')
        differences={str(k):numeric_text(reconciled[k]-v) for k,v in reference.items()}
        require(all(abs(reconciled[k]-v)<=max(Decimal('1e-60'),abs(v)*Decimal('1e-70')) for k,v in reference.items()),
            'Subgroup age/sex counts do not reconcile to the parent population. No source counts are adjusted to hide a discrepancy.')
    native=None
    if subgroup_results:
        policy=dict(mode='reviewed_disjoint_partition_v1',partition_id='reviewer_declared_partition',subgroup_ids=[g['subgroup_id'] for g in subgroup_results],**META)
        with localcontext() as ctx:
            ctx.prec=80;covered=numeric_text(sum(decimal(g['population']) for g in subgroup_results))
        native=population_mixture(subgroup_results,mixture_id='review_population_covered_source_scope',geography_id=geo,geography_level=level,
            total_population=covered,membership_policy=policy,**META)
    summary=profile_summary(all_rows)
    return seal(dict(schema_id='country_tb.review_ui_baseline_result.v1',geography_id=geo,geography_level=level,
        subgroups=subgroup_results,rows=all_rows,summary=summary,unresolved_subgroups=missing,
        native_mixture=native,declared_population=context['demography']['population'],complete_mixture_available=not missing and summary['total_tbi_complete'] and summary['recency_complete'],
        source_composition=groups,age_sex_reconciliation_residuals=differences,production_available=False,externally_validated=False))


def allocate_risk(config,baseline):
    risk=config['risk'];cat=catalog();entries=catalogue_records(cat);ids=[g['subgroup_id'] for g in baseline['subgroups']]
    upload=parse_csv(risk['csv'],'risk',subgroup_ids=ids) if risk['csv'] else None
    overrides=parse_csv(config['association']['csv'],'association',subgroup_ids=ids) if config['association']['csv'] else None
    all_allocations=[]
    for g in baseline['subgroups']:
        raw=deepcopy(upload['rows'] if upload else risk['records']);raw=[r for r in raw if r.get('subgroup_id','main') in {g['subgroup_id'],'all'}]
        if risk['no_factors']:
            require(not raw,'Clear risk records before selecting the explicit RR1 reference scenario.')
            raw=[dict(record_id=f,subgroup_id=g['subgroup_id'],factor_id=f,age_lower=15,age_upper_exclusive=None,sex='both',
                state='excluded',prevalence=None,allocation_policy=UNIFORM,supplied_exposure_definition=e['exposure_definition'],
                source='Reviewer explicitly excluded all factors: RR1 reference,not a population fact',review_status='diagnostic',reason='Explicit RR1 reference scope') for f,e in entries.items()]
        scope=sorted({r['factor_id'] for r in raw});require(scope,'Enter risk data or explicitly choose the RR1 reference scenario; missing factors are not zero.')
        records=[]
        for r in raw:
            require(r['factor_id'] in entries,'Unknown factor; recent infection and close contact are not ordinary factors.')
            records.append({**r,'source':{'source_id':r['source']},'review':REVIEW,'provenance':{**PROVENANCE,'original_record':r},
                'definition_mapping':None,'pattern':None,'reason':r.get('reason') or None})
        profile=risk_factor_prevalence_profile(records,catalog=cat,factor_ids=scope,subgroup_id=g['subgroup_id'],geography_id=baseline['geography_id'],**META)
        rows=[]
        for r in overrides['rows'] if overrides else []:
            if r['subgroup_id'] not in {g['subgroup_id'],'all'}:continue
            rows.append({k:v for k,v in {**r,'source':{'source_id':r['source']},'review':REVIEW,'provenance':{**PROVENANCE,'original_record':r}}.items() if k not in {'subgroup_id','review_status'}})
        association=association_policy(global_control=config['association']['global_control'],overrides=rows,**META)
        policy=multimorbidity_policy(effect_use='diagnostic_provisional_adult_effects_v1' if risk['legacy_effects'] else 'reviewed_effects_only_v1',**META)
        all_allocations.append(allocate_baseline_risk(g['tbi_profile'],catalog=cat,prevalence=profile,overlap=overlap_policy(**META),
            association=association,progression_policy=policy,diagnostic_only=True,**META))
    native=combine_risk_allocations(baseline['native_mixture'],all_allocations,**META) if all_allocations else None
    return seal(dict(schema_id='country_tb.review_ui_risk_result.v1',geography_id=baseline['geography_id'],
        subgroups=all_allocations,rows=[r for a in all_allocations for r in a['rows']],native_mixture=native,
        declared_population=baseline['declared_population'],unresolved_subgroups=baseline['unresolved_subgroups'],
        complete_mixture_available=baseline['complete_mixture_available'],production_available=False,externally_validated=False))


def execute_user_run(config,context=None):
    require(config['tbi']['mode']==SOURCE_MODES[0],'No fallback from incidence mode to user TBI.')
    from country_model_review_ui.forward import forward_diagnostic
    timings={};start=perf_counter();provided=context is not None;context=context or select_population(config)
    timings['population_context_reuse_seconds' if provided else 'population_selection_seconds']=perf_counter()-start
    start=perf_counter();baseline=build_baseline(config,context);timings['TBI_construction_seconds']=perf_counter()-start
    start=perf_counter()
    try:allocation=allocate_risk(config,baseline)
    except DataValidationError as exc:
        allocation=seal(dict(schema_id='country_tb.review_ui_risk_result.v1',geography_id=baseline['geography_id'],
            subgroups=[],rows=[],native_mixture=None,declared_population=baseline['declared_population'],
            unresolved_subgroups=[dict(subgroup_id='risk_inputs_unresolved',population=baseline['declared_population'],status=str(exc))],
            complete_mixture_available=False,input_error=str(exc),production_available=False,externally_validated=False))
    timings['risk_allocation_seconds']=perf_counter()-start
    start=perf_counter();projection=forward_diagnostic(config,allocation);timings['deterministic_forward_seconds']=perf_counter()-start
    result=seal(dict(schema_id='country_tb.review_ui_run.v1',application_version=APP_VERSION,configuration_hash=material_hash(config),
        configuration=deepcopy(config),context=context,baseline=baseline,risk_allocation=allocation,projection=projection,
        production_available=False,externally_validated=False,classification='internal_scientific_review_only',
        provenance={'geography':baseline['geography_id'],'population_definition':deepcopy(config['population']),
            'snapshot_identifiers':{**context['source_identifiers'],'selected_mortality':
                'abs-life2022-2024-6785449a5200' if config['natural']['mortality']!='No mortality — diagnostic only' else None},
            'baseline_TBI_source':deepcopy(config['tbi']),
            'risk_profile_hash':allocation['profile_hash'],'association':deepcopy(config['association']),
            'progression_policies':deepcopy(config['natural']),'numerical_solver':projection.get('numerical_policy'),
            'scientific_hashes':{'population':context['demography']['profile_hash'],'baseline':baseline['profile_hash'],
                'risk':allocation['profile_hash'],'projection':projection['profile_hash']},
            'application_version':APP_VERSION,'contracts':{'baseline':'2D/v1','risk':'3A/v1','progression':'3B/v1','inversion':'3C/v1','evidence':'3D/v1'}}))
    return result,timings
