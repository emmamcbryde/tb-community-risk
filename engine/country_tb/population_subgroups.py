"""Reviewed disjoint population partitions and population-weighted TBI mixtures."""
from copy import deepcopy
from datetime import date
from decimal import Decimal

from engine.country_tb.calibration_common import reviewed,scientific_context,probability
from engine.country_tb.composition_common import (
    artifact,verify,require,metadata,identity,validate_demography,profile_summary,STATES,
)
from engine.country_tb.data_common import canonical,decimal,numeric_text
from engine.country_tb.baseline_tbi_source import SOURCE,RESULT,INCIDENCE,resolve_baseline_tbi,baseline_tbi_source

SUBGROUP='country_tb.population_subgroup.v1'
MIXTURE='country_tb.population_mixture.v1'


def _geo(value):
    require(isinstance(value,dict) and value.get('geography_id') and
        value.get('geography_level') in {'national','subnational','custom'},'Explicit optional geography object required, or null.')


def validate_resolution(value):
    verify(value,RESULT)
    selection=value['selection'];verify(selection,SOURCE)
    require(value['selection_hash']==selection['profile_hash'] and value['demography_hash']==value['demography']['profile_hash'],
            'Resolution source/demography references disagree.')
    require(value['mode']==selection['mode'] and value['subgroup_id']==selection['subgroup_id'] and
        value['geography_id']==selection['geography_id'],'Resolution identity differs from its selected source.')
    validate_demography(value['demography'])
    expected={(r['sex'],r['age_lower'],r['age_upper_exclusive']):decimal(r['population']) for r in value['demography']['rows']}
    actual={(r['sex'],r['age_lower'],r['age_upper_exclusive']):decimal(r['population']) for r in value['rows']}
    require(expected==actual and len(actual)==len(value['rows']),'Resolved rows duplicate or lose demographic cells.')
    with scientific_context():
        for row in value['rows']:
            for field in STATES:
                if row[field] is not None:probability(row[field])
            if row['any_tbi'] is not None:
                require(row['never_infected'] is not None and decimal(row['any_tbi'])+decimal(row['never_infected'])==1,
                        'Never infected and any TBI do not reconcile.')
            require((row['recent'] is None)==(row['remote_only'] is None),'Partial direct recency states are invalid.')
            if row['recent'] is not None:
                require(decimal(row['recent'])+decimal(row['remote_only'])==decimal(row['any_tbi']),
                        'Recent plus remote-only must equal any TBI without adding the prior-remote subset.')
            if row['prior_remote_plus_recent'] is not None:
                require(row['recent'] is not None and decimal(row['prior_remote_plus_recent'])<=decimal(row['recent']),
                        'Prior remote plus recent must remain a subset of recent.')
        require(value['summary']==profile_summary(value['rows']),'Resolved summary does not reconcile with rows.')
        expected_usable=(value['summary']['total_tbi_complete'] and value['summary']['recency_complete'] and
            value['profile_kind']!='diagnostic' and selection['profile_kind']!='diagnostic' and
            selection['review']['status']=='reviewed' and selection['mode']!=INCIDENCE)
        require(value['usable_for_downstream_simulation']==expected_usable,'Resolution usability is inconsistent with component coverage, review or diagnostic status.')
        if value['mode']==INCIDENCE:
            require(all(r[k] is None for r in value['rows'] for k in STATES) and value['tbi_estimate'] is None and
                not value['usable_for_downstream_simulation'],'Incidence mode cannot contain a TBI estimate in this milestone.')
    return value


def population_subgroup(*,subgroup_id,label,description,current_residence_geography,geography_level,
    population,demography,tbi_source,membership,source,provenance,review,baseline_year=2024,
    resolved_profile=None,exposure_history_source=None,country_or_region_of_birth=None,previous_residence=None,
    exposure_history_geography=None,arrival_date=None,arrival_year=None,duration_in_current_setting_years=None,
    risk_factor_profile=None):
    identity(current_residence_geography,geography_level,subgroup_id,baseline_year)
    require(isinstance(label,str) and label and isinstance(description,str) and description,'Subgroup label and description required.')
    metadata(source,provenance,review)
    validate_demography(demography);verify(tbi_source,SOURCE)
    expected=baseline_tbi_source(**{k:tbi_source[k] for k in ('mode','geography_id','geography_level','subgroup_id',
        'source','provenance','review','profile_kind','configuration','baseline_year')})
    require(expected['profile_hash']==tbi_source['profile_hash'],'Invalid subgroup source selection.')
    require((tbi_source['subgroup_id'],tbi_source['geography_id'],tbi_source['geography_level'],tbi_source['baseline_year'])==
        (subgroup_id,current_residence_geography,geography_level,baseline_year),'Subgroup and TBI-source identity mismatch.')
    require(demography['geography_id']==current_residence_geography and demography['baseline_year']==baseline_year,
            'Subgroup current-residence demography does not match.')
    with scientific_context():
        total=sum(decimal(r['population']) for r in demography['rows'])
        require(decimal(population)>0 and decimal(population)==total,'Subgroup denominator must match all age/sex cells, including the open tail.')
    require(isinstance(membership,dict) and {'partition_id','category_id','definition','source','provenance','review','person_ids'}<=membership.keys(),
            'Explicit disjoint-partition membership metadata required, with nullable person IDs.')
    metadata(membership['source'],membership['provenance'],membership['review'])
    require(all(isinstance(membership[k],str) and membership[k] for k in ('partition_id','category_id','definition')),
            'Membership partition/category IDs and definition required.')
    members=deepcopy(membership)
    if members['person_ids'] is not None:
        ids=members['person_ids']
        require(isinstance(ids,list) and all(isinstance(x,str) and x for x in ids) and len(set(ids))==len(ids),
                'Enumerated person IDs must be unique and nonempty.')
        require(Decimal(len(ids))==total,'Enumerated membership must account for the complete subgroup denominator.')
        members['person_ids']=sorted(ids)
    require(risk_factor_profile is None,'Risk-factor population allocation is deferred; placeholder must remain null.')
    for geo in (country_or_region_of_birth,exposure_history_geography):
        if geo is not None:_geo(geo)
    if previous_residence is not None:
        require(isinstance(previous_residence,list),'Previous residence must be an explicit list or null.')
        for geo in previous_residence:_geo(geo)
    if exposure_history_source is not None:
        require(isinstance(exposure_history_source,dict) and exposure_history_source.get('source_id'),
                'Exposure-history source requires source_id or explicit null.')
    if arrival_year is not None:require(type(arrival_year) is int and 0<arrival_year<=2024,'Invalid arrival year.')
    if arrival_date is not None:
        try:parsed=date.fromisoformat(arrival_date)
        except (ValueError,TypeError):require(False,'Arrival date must be an ISO date.')
        require(parsed.year<=2024 and (arrival_year is None or parsed.year==arrival_year),'Arrival date and year disagree or postdate baseline.')
    if duration_in_current_setting_years is not None:
        duration=decimal(duration_in_current_setting_years)
        if arrival_year is not None:require(duration<=2025-arrival_year,'Residence duration contradicts the supplied arrival year.')
    resolution=resolve_baseline_tbi(tbi_source,demography) if resolved_profile is None else validate_resolution(resolved_profile)
    require(resolution['selection_hash']==tbi_source['profile_hash'] and resolution['demography_hash']==demography['profile_hash'],
            'Subgroup resolved profile belongs to a different source or population.')
    with scientific_context():
        age={};sex={}
        for row in demography['rows']:
            key=(row['age_lower'],row['age_upper_exclusive'])
            age[key]=age.get(key,Decimal(0))+decimal(row['population'])
            sex[row['sex']]=sex.get(row['sex'],Decimal(0))+decimal(row['population'])
        return artifact(SUBGROUP,subgroup_id=subgroup_id,label=label,description=description,
            current_residence_geography=current_residence_geography,geography_level=geography_level,baseline_year=baseline_year,
            population=numeric_text(total),demography=demography,demography_hash=demography['profile_hash'],
            age_distribution=[{'age_lower':lo,'age_upper_exclusive':hi,'population':numeric_text(n),'weight':numeric_text(n/total)}
                for (lo,hi),n in sorted(age.items(),key=lambda x:x[0][0])],
            sex_distribution=[{'sex':s,'population':numeric_text(n),'weight':numeric_text(n/total)} for s,n in sorted(sex.items())],
            tbi_source_mode=tbi_source['mode'],tbi_source=tbi_source,tbi_profile=resolution,tbi_profile_hash=resolution['profile_hash'],
            membership=members,exposure_history_source=exposure_history_source,
            residence_metadata={'country_or_region_of_birth':country_or_region_of_birth,'previous_residence':previous_residence,
                'exposure_history_geography':exposure_history_geography,'arrival_date':arrival_date,'arrival_year':arrival_year,
                'duration_in_current_setting_years':numeric_text(duration_in_current_setting_years) if duration_in_current_setting_years is not None else None},
            risk_factor_profile={'status':'deferred','profile':None},source=source,provenance=provenance,review=review,
            migration_policy='metadata_only_no_inferred_arrival_exposure_or_migration_model')


def validate_subgroup(value):
    verify(value,SUBGROUP);validate_resolution(value['tbi_profile'])
    require(value['tbi_profile_hash']==value['tbi_profile']['profile_hash'],'Subgroup profile reference mismatch.')
    kwargs={k:value[k] for k in ('subgroup_id','label','description','current_residence_geography','geography_level',
        'population','demography','tbi_source','membership','source','provenance','review','baseline_year','exposure_history_source')}
    kwargs.update(value['residence_metadata']);kwargs['resolved_profile']=value['tbi_profile']
    rebuilt=population_subgroup(**kwargs)
    require(rebuilt['profile_hash']==value['profile_hash'],'Subgroup does not reproduce from its validated source, population and membership.')
    return value


def population_mixture(subgroups, *,mixture_id,geography_id,geography_level,total_population,
    membership_policy,source,provenance,review,baseline_year=2024,subgroup_weights=None,
    reference_age_sex_distribution=None,geography_policy=None):
    identity(geography_id,geography_level,mixture_id,baseline_year)
    metadata(source,provenance,review)
    require(isinstance(subgroups,list) and bool(subgroups),'At least one population subgroup required.')
    for subgroup in subgroups:validate_subgroup(subgroup)
    ids=[g['subgroup_id'] for g in subgroups]
    require(len(set(ids))==len(ids),'Subgroup IDs must be unique; duplicated subgroups double-count membership.')
    require(isinstance(membership_policy,dict) and membership_policy.get('mode')=='reviewed_disjoint_partition_v1' and
        {'partition_id','subgroup_ids','source','provenance','review'}<=membership_policy.keys(),
        'A reviewed disjoint membership partition is required; group labels alone cannot establish non-overlap.')
    metadata(membership_policy['source'],membership_policy['provenance'],membership_policy['review'])
    require(isinstance(membership_policy['subgroup_ids'],list) and all(isinstance(x,str) for x in membership_policy['subgroup_ids']),
            'Partition subgroup IDs must be an explicit list.')
    require(sorted(membership_policy['subgroup_ids'])==sorted(ids),'Partition declaration must cover exactly the selected subgroup IDs.')
    require(all(g['membership']['partition_id']==membership_policy['partition_id'] for g in subgroups),'Membership partition mismatch.')
    require(len({g['membership']['category_id'] for g in subgroups})==len(subgroups),'Duplicate subgroup membership category.')
    if geography_policy is not None:
        require(isinstance(geography_policy,dict) and geography_policy.get('mode')=='reviewed_containment_v1' and
            isinstance(geography_policy.get('contained_geography_ids'),list) and
            all(isinstance(x,str) and x for x in geography_policy['contained_geography_ids']) and
            len(set(geography_policy['contained_geography_ids']))==len(geography_policy['contained_geography_ids']),
            'Reviewed geography containment requires an explicit unique geography-ID list.')
        metadata(geography_policy.get('source'),geography_policy.get('provenance'),geography_policy.get('review'))
    seen=set()
    for g in subgroups:
        members=set(g['membership']['person_ids'] or [])
        require(not members&seen,'Duplicated person membership between subgroups.')
        seen|=members
        require(g['baseline_year']==baseline_year,'Subgroup baseline years differ.')
        if (g['current_residence_geography'],g['geography_level'])!=(geography_id,geography_level):
            require(isinstance(geography_policy,dict) and geography_policy.get('mode')=='reviewed_containment_v1' and
                geography_policy.get('parent_geography_id')==geography_id and
                g['current_residence_geography'] in geography_policy.get('contained_geography_ids',[]),
                'Subgroup geography is incompatible without explicit reviewed containment evidence.')
            metadata(geography_policy.get('source'),geography_policy.get('provenance'),geography_policy.get('review'))
    grids=[{(r['sex'],r['age_lower'],r['age_upper_exclusive']) for r in g['demography']['rows']} for g in subgroups]
    require(all(grid==grids[0] for grid in grids),'Mixture age/sex grids must align; no implicit population rebinning or sex split.')
    with scientific_context():
        total=decimal(total_population)
        require(total>0 and sum(decimal(g['population']) for g in subgroups)==total,'Subgroup populations do not reconcile to the mixture denominator.')
        derived={g['subgroup_id']:decimal(g['population'])/total for g in subgroups}
        if subgroup_weights is not None:
            require(isinstance(subgroup_weights,dict) and set(subgroup_weights)==set(ids),'Explicit weights must cover every subgroup exactly.')
            require(abs(sum(probability(x) for x in subgroup_weights.values())-1)<=Decimal('1e-12'), 'Subgroup weights must sum to one.')
            require(all(abs(probability(subgroup_weights[k])-v)<=Decimal('1e-12') for k,v in derived.items()),
                    'Weights disagree with population counts; counts determine mixture arithmetic.')
        cells={key:[] for key in grids[0]}
        for g in sorted(subgroups,key=lambda x:x['subgroup_id']):
            for row in g['tbi_profile']['rows']:
                cells[row['sex'],row['age_lower'],row['age_upper_exclusive']].append(row)
        rows=[]
        for (sex,lo,hi),values in sorted(cells.items(),key=lambda x:(x[0][0],x[0][1])):
            s=profile_summary(values)
            rows.append({'sex':sex,'age_lower':lo,'age_upper_exclusive':hi,'population':s['population'],
                **{field:s['components'][field]['prevalence_full'] for field in STATES},'partial_component_summary':s})
        if reference_age_sex_distribution is not None:
            reference={(r['sex'],r['age_lower'],r['age_upper_exclusive']):decimal(r['population']) for r in reference_age_sex_distribution}
            require(len(reference)==len(reference_age_sex_distribution) and reference=={
                (r['sex'],r['age_lower'],r['age_upper_exclusive']):decimal(r['population']) for r in rows},
                'Combined age/sex population does not reconcile to the supplied reference distribution.')
        # Aggregate original subgroup cells: partially covered mixture cells
        # must not lose the known subgroup contributions when becoming null.
        summary=profile_summary([r for g in subgroups for r in g['tbi_profile']['rows']])
        require(sum(decimal(r['population']) for r in rows)==total,'Combined age/sex denominator does not reconcile.')
        complete=summary['total_tbi_complete'] and summary['recency_complete']
        unresolved=[{'subgroup_id':g['subgroup_id'],'population':g['population'],
            'status':g['tbi_profile']['status'],'unresolved_population':g['tbi_profile']['summary']['unresolved_population'],
            'recency_unresolved_population':g['tbi_profile']['summary']['components']['recent']['unresolved_population']}
            for g in subgroups if not(g['tbi_profile']['summary']['total_tbi_complete'] and g['tbi_profile']['summary']['recency_complete'])]
        policy=deepcopy(membership_policy);policy['subgroup_ids']=sorted(policy['subgroup_ids'])
        geo_policy=deepcopy(geography_policy)
        if geo_policy is not None:geo_policy['contained_geography_ids']=sorted(geo_policy['contained_geography_ids'])
        return artifact(MIXTURE,mixture_id=mixture_id,geography_id=geography_id,geography_level=geography_level,baseline_year=baseline_year,
            subgroups=sorted(subgroups,key=lambda g:g['subgroup_id']),rows=rows,summary=summary,
            total_population=numeric_text(total),subgroup_weights={k:numeric_text(v) for k,v in sorted(derived.items())},
            supplied_subgroup_weights={k:numeric_text(v) for k,v in sorted(subgroup_weights.items())} if subgroup_weights is not None else None,
            weight_rounding_tolerance='0.000000000001',membership_policy=policy,geography_policy=geo_policy,
            membership_validation='reviewed_disjoint_partition_plus_available_enumerated_person_ID_checks',
            status='complete_baseline_state_coverage' if complete else 'partial_component_coverage',
            complete_mixture_result=deepcopy(summary) if complete else None,partial_summaries=summary,
            unresolved_subgroups=sorted(unresolved,key=lambda r:r['subgroup_id']),
            source_composition=[{'subgroup_id':g['subgroup_id'],'mode':g['tbi_source_mode'],'population':g['population'],
                'weight':numeric_text(derived[g['subgroup_id']]),'profile_kind':g['tbi_profile']['profile_kind'],
                'profile_hash':g['tbi_profile_hash'],'recent_remote_composition':g['tbi_profile']['recent_remote_composition']}
                for g in sorted(subgroups,key=lambda g:g['subgroup_id'])],
            contains_diagnostic_profiles=any(g['tbi_profile']['profile_kind']=='diagnostic' for g in subgroups),
            usable_for_downstream_simulation=complete and all(g['tbi_profile']['usable_for_downstream_simulation'] for g in subgroups),
            downstream_scope='baseline_exclusive_state_initialization_only_no_progression_or_interventions',
            source=source,provenance=provenance,review=review,
            migration_policy='current_residence_is_not_assumed_to_be_historical_exposure')
