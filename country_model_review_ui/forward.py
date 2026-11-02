"""Bounded review adapters to3B natural-history and forward contracts.

Instantaneous integration uses the retained2C conditional last-event CDF and3B
piecewise hazards. It never substitutes mean infection time or an age midpoint.
Optional prospective calculations call the unchanged3B operator,with a node
budget checked before constructing its expanded cohort.
"""
from collections import defaultdict
from copy import deepcopy
from math import fsum,ceil,floor
from engine.country_tb.data_common import seal,decimal,numeric_text,DataValidationError
from engine.country_tb.risk_common import ratio
from engine.country_tb.catalytic_core import history_curve
from engine.country_tb.most_recent_event import event_distribution,conditional_cdf
from engine.country_tb.progression_natural_history import audited_adult_candidate,curve_hazard
from engine.country_tb.progression_baseline import cohort_from_allocation
from engine.country_tb.progression_forward import forward_operator,project_active_tb
from engine.country_tb.progression_policies import effect_application,reinfection_clock,competing_mortality,future_pressure
from engine.country_tb.review_parameters_3d import review_parameter_set
from engine.country_tb.validation_data import load_mortality_validation_snapshot
from country_model_review_ui.state import require
from country_model_review_ui.model import META

RR74_WARNING='RR74 sensitivity gives94.31% one-year progression under the3.8% adult surrogate. This is a scientific warning,not evidence of plausibility.'


def selected_history(config,*,set_id=None):
    n=config['natural'];require(n['adult']!='External required','Adult progression curve is unresolved.')
    if set_id is None:
        set_id={'Conservative four-point sensitivity':'conservative','Central-candidate four-point sensitivity':'central_candidate',
            'Higher-progression four-point sensitivity':'higher_progression'}[n['adult']]
    if n['tail']=='No extrapolation':return audited_adult_candidate(tail='no_extrapolation_v1') if set_id=='central_candidate' else _with_tail(set_id,'no_extrapolation_v1',{})
    if n['tail']=='Constant-final-segment sensitivity':return _with_tail(set_id,'constant_final_segment_hazard_v1',{})
    h=decimal(n['remote_hazard']);return _with_tail(set_id,'fixed_remote_hazard_v1',{'hazard_per_year':numeric_text(h)})


def _with_tail(identifier,tail,parameters):
    from engine.country_tb.progression_natural_history import curve_specification,natural_history
    original=review_parameter_set(identifier)['natural_history']['group_curves']['15+']
    curve=curve_specification([{'time_years':r['time_years'],'cumulative_risk':r['cumulative_risk']} for r in original['knots']],
        tail=tail,tail_parameters=parameters,status='diagnostic_adult_surrogate',
        risk_interpretation='diagnostic_cumulative_risk_surrogate_pending_cause_specific_review_v1',**META)
    return natural_history({'0-4':None,'5-14':None,'15+':curve},**META)


def _conditional_hazard(model,cell,state,curve):
    distribution=event_distribution(history_curve(model['calibration']['pressure_history'],cell['sex']),
        alpha=model['calibration']['alpha'],baseline_time=model['baseline_time'],age_lower=cell['age_lower'],
        age_upper=cell['age_upper_exclusive'],recent_window=model['recent_window_years'],group=state)
    time=decimal(model['baseline_time']);end=decimal(cell['age_upper_exclusive'])
    def mass(lo,hi):
        a=conditional_cdf(distribution,time-decimal(str(lo)));b=conditional_cdf(distribution,time-decimal(str(hi)))
        value=decimal(a)-decimal(b);require(value>=0,'Conditional CDF mass became negative; no clipping.')
        return float(value)
    value=0.;unresolved=0.
    for segment in curve['segments']:
        lo=float(segment['start']);hi=min(float(segment['end']),float(end))
        if hi>lo:value+=mass(lo,hi)*float(segment['hazard_per_year'])
    last=float(curve['knots'][-1]['time_years'])
    if float(end)>last:
        weight=mass(last,float(end))
        if curve['tail']=='no_extrapolation_v1':unresolved=weight
        else:value+=weight*curve_hazard(curve,last+1e-8)
    return value,unresolved


def instantaneous(config,allocation,nh):
    outputs=[];known_flow=[];unresolved=[];population=float(allocation['declared_population'])
    for group in allocation['subgroups']:
        baseline=group['resolved_baseline'];model=baseline.get('model_derived_profile')
        for cell in group['rows']:
            N=float(ratio(cell['population']));key=dict(subgroup_id=group['subgroup_id'],sex=cell['sex'],age_lower=cell['age_lower'],age_upper_exclusive=cell['age_upper_exclusive'])
            if not N:continue
            p=None if cell['tbi_prevalence'] is None else float(ratio(cell['tbi_prevalence']))
            reasons=[];missing=flow=0.
            if p is None:missing=N;reasons=['TBI_unresolved']
            elif cell['age_upper_exclusive'] is None:missing=N*p;reasons=['open_age_distribution_unresolved']
            elif cell['age_lower']<15:missing=N*p;reasons=['pediatric_progression_external_required']
            elif model is None or cell['tbi_state_prevalences']['recent'] is None:
                missing=N*p;reasons=['conditional_most_recent_infection_distribution_unavailable']
            elif not cell['compressed_cells']:
                missing=N*p;reasons=['risk_prevalence_or_definition_unresolved']
            else:
                curve=nh['group_curves']['15+'];hazards={}
                for state in ('recent','remote_only'):
                    value=cell['tbi_state_prevalences'][state]
                    if value is None or float(ratio(value))==0:hazards[state]=(0.,0.);continue
                    hazards[state]=_conditional_hazard(model,cell,state,curve)
                for part in cell['compressed_cells']:
                    rr=None if part['effective_rr'] is None else float(ratio(part['effective_rr']))
                    for state in ('recent','remote_only'):
                        count=part['tbi_state_counts'][state]
                        if count is None:continue
                        mass=float(ratio(count))
                        if not mass:continue
                        if rr is None:missing+=mass;reasons.append('unreviewed_effect_multiplier');continue
                        h,gap=hazards[state];flow+=mass*rr*h;missing+=mass*gap
                        if gap:reasons.append('remote_tail_external_required_not_zero')
            require(missing<=N+max(1e-9,N*1e-12),'Unresolved population exceeded stratum denominator.')
            known_flow.append(flow);unresolved.append(missing)
            outputs.append({**key,'population':N,'covered_instantaneous_case_flow_per_year':flow,
                'unresolved_population':missing,'complete_rate_per100000PY':None if missing>1e-10 else flow/N*100000,
                'unresolved_reasons':sorted(set(reasons)),'label':'Diagnostic sensitivity'})
    for g in allocation['unresolved_subgroups']:
        missing=float(g['population']);unresolved.append(missing)
        outputs.append(dict(subgroup_id=g['subgroup_id'],population=missing,sex='unresolved',age_lower=None,
            age_upper_exclusive=None,covered_instantaneous_case_flow_per_year=None,unresolved_population=missing,
            complete_rate_per100000PY=None,unresolved_reasons=[g['status']],label='Unresolved'))
    gap=fsum(unresolved);flow=fsum(known_flow)
    return dict(rows=outputs,population=population,unresolved_population=gap,covered_population=population-gap,
        covered_instantaneous_case_flow_per_year=flow,full_population_rate_per100000PY=None if gap>0 else flow/population*100000,
        covered_population_rate_per100000PY=None if population<=gap else flow/(population-gap)*100000,
        incidence_units='instantaneous_first_episode_case_flow_per100000_disease_free_person_years_not_cumulative_one_year_cases',
        mortality_interpretation='Future survival has no effect on an instantaneous baseline hazard; no long-horizon country estimate.',
        classification='diagnostic_sensitivity_not_country_estimate')


def prospective(config,allocation,nh):
    native=allocation['native_mixture'];require(native is not None,'No resolved subgroup risk allocation for prospective projection.')
    estimate=0
    for g in allocation['subgroups']:
        for r in g['rows']:
            if r['age_upper_exclusive'] is not None and float(ratio(r['population']))>0:
                estimate+=8*(8*ceil(r['age_upper_exclusive'])+1)*max(1,len(r['effective_rr_distribution'] or []))
    require(estimate<=50000,'Prospective support exceeds the declared50000-node review budget. Use instantaneous results or an explicitly defined smaller review population; no population is dropped.')
    cohort=cohort_from_allocation(native,active_TB_exclusion_policy='assumed_zero_baseline_active_TB_diagnostic_v1',max_time_interval_years='1',**META)
    n=config['natural'];geo=allocation['geography_id']
    if n['mortality']=='No mortality — diagnostic only':mortality=competing_mortality(geography_id=geo,mode='no_mortality_diagnostic_v1',**META)
    else:
        require(config['population']['country']=='AUS','ABS mortality is unavailable for other countries.')
        require(geo=='who:AUS' or n.get('mortality_proxy',False),'Explicit national mortality proxy selection required for a custom population.')
        frozen=load_mortality_validation_snapshot('abs-life2022-2024-6785449a5200',geography_id='who:AUS')
        metadata=frozen['metadata']
        needed={(age,r['sex']) for r in cohort['rows'] if r['age'] is not None and not r['unresolved_reasons']
            for age in range(floor(float(r['age'])),floor(float(r['age'])+1)+1)}
        rows=[r for r in frozen['rows'] if (r['age_lower'],r['sex']) in needed]
        require(rows,'No mortality intervals cover the represented cohort; no fallback.')
        mortality=competing_mortality(geography_id=geo,mode='external_age_sex_survival_v1',rows=rows,
            snapshot_hash=frozen['profile_hash'],licensing={'status':'reviewed','licence':metadata['licence'],
                'scope':'permission_review_only_not_scientific_adoption'},
            **{**META,'source':{'source_id':'abs-life2022-2024-6785449a5200','citation':metadata['citation'],'url':metadata['url']},
                'provenance':{**META['provenance'],'immutable_snapshot_hash':frozen['profile_hash'],
                    'use':'explicit_validation_sensitivity_not_adopted_competing_hazard','all_cause_TB_death_scope':'unresolved',
                    'period_table_held_fixed_one_year_sensitivity':True,'source_interval_count':len(frozen['rows']),
                    'represented_age_sex_intervals_only_no_interpolation_or_aggregation':True,'selected_interval_count':len(rows)}})
    operator=forward_operator(natural_history=nh,effect=effect_application(**META),clock=reinfection_clock(**META),mortality=mortality,
        future=future_pressure(geography_id=geo,scenario='zero',horizon_years='1',**META),mode='baseline_cohort_only_v1',
        step_years='.25',horizons=['1'],diagnostic_only=True,**META)
    return project_active_tb(cohort,operator)


def forward_diagnostic(config,allocation):
    n=config['natural'];base=dict(schema_id='country_tb.review_ui_forward_result.v1',production_available=False,externally_validated=False,
        assumptions=n,numerical_policy={'instantaneous':'exact_conditional_last_event_CDF_mass_times_piecewise_hazard',
            'age_integration':'retained_uniform_continuous_within_single_year_no_midpoints',
            'prospective':'unchanged3B_GL8_support_max1year,step.25year,50000node_budget'},
        warning='Sensitivity assumptions only; active baseline disease,prior TB,recurrence,migration and country transport unresolved.')
    require(n['reset']=='Most recent reinfection resets clock' and n['future']=='No future infection — baseline cohort only',
        'This UI requires explicit central clock reset and zero future pressure; unknown policies do not fall back.')
    if n['adult']=='External required' or n['active_baseline']=='Unresolved':
        return seal(dict(**base,status='calculation_blocked',reason='Select an explicit adult sensitivity curve and baseline-active-TB exclusion assumption.',
            unresolved_population=allocation['declared_population'],instantaneous=None,scenarios=[],prospective=None))
    if n['mortality']!='No mortality — diagnostic only' and config['population']['country']!='AUS':
        return seal(dict(**base,status='calculation_blocked',reason='ABS mortality is incompatible with this country; choose an appropriate explicit mortality policy.',
            unresolved_population=allocation['declared_population'],instantaneous=None,scenarios=[],prospective=None))
    nh=selected_history(config);value=instantaneous(config,allocation,nh);scenarios=[];projection=None;projection_error=None
    if n['compare_scenarios']:
        for identifier in ('conservative','central_candidate','higher_progression'):
            h=selected_history(config,set_id=identifier);v=instantaneous(config,allocation,h)
            scenarios.append(dict(parameter_set_id=identifier,**v,scientific_hash=h['profile_hash']))
        for tail in ('No extrapolation','Constant-final-segment sensitivity'):
            alternative=deepcopy(config);alternative['natural']['tail']=tail;h=selected_history(alternative)
            scenarios.append(dict(parameter_set_id='Remote tail: '+tail,**instantaneous(alternative,allocation,h),scientific_hash=h['profile_hash']))
        scenarios.append(dict(parameter_set_id='Mortality: baseline hazard invariant to future survival',**value,
            scientific_hash=nh['profile_hash'],interpretation='Instantaneous comparison only; not a long-horizon or mortality-adjusted cumulative risk.'))
    if n['horizon']=='One-year prospective (bounded)':
        try:projection=prospective(config,allocation,nh)
        except DataValidationError as exc:projection_error=str(exc)
    return seal(dict(**base,status='diagnostic_sensitivity' if value['unresolved_population']==0 else 'partial_diagnostic_sensitivity',
        instantaneous=value,scenarios=scenarios,prospective=projection,prospective_unavailable_reason=projection_error,
        unresolved_population=value['unresolved_population'],natural_history=nh))
