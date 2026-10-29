"""Reviewed/external age curves and an explicitly diagnostic adult surrogate."""
from bisect import bisect_right
from copy import deepcopy
from decimal import Decimal
from math import exp,expm1,log1p

from engine.country_tb.progression_common import (
    artifact,verify,require,metadata,number,num,scientific_context,decimal,DataValidationError,
)

NATURAL_HISTORY='country_tb.progression_natural_history.v1'
GROUPS=((0,5,'0-4'),(5,15,'5-14'),(15,None,'15+'))
TAILS={'no_extrapolation_v1','constant_final_segment_hazard_v1','declining_to_remote_plateau_v1','fixed_remote_hazard_v1'}
ADULT_SOURCE={'source_id':'Menzies_NA_2021_Epidemiology','doi':'10.1097/EDE.0000000000001271',
    'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7707158/','figure':'Figure 1B',
    'population':'Modelled individual infected at age 20 with no additional progression risk factors'}


def risk_to_hazard(risk):
    d=decimal(risk);require(d<1,'Finite cumulative hazard needs risk below one.')
    with scientific_context():return num(-(1-d).ln())


def curve_specification(knots,*,tail,source,provenance,review,status='reviewed_external',tail_parameters=None,risk_interpretation=None):
    metadata(source,provenance,review)
    require(status in {'reviewed_external','diagnostic_synthetic','diagnostic_adult_surrogate'},'Unknown curve evidence status.')
    if risk_interpretation is None and status=='diagnostic_synthetic':risk_interpretation='net_risk_without_competing_mortality_v1'
    require(risk_interpretation in {'net_risk_without_competing_mortality_v1','diagnostic_cumulative_risk_surrogate_pending_cause_specific_review_v1'},
        'Cumulative incidence with competing death is not automatically a cause-specific hazard; declare the risk interpretation.')
    require(status!='reviewed_external' or risk_interpretation=='net_risk_without_competing_mortality_v1',
        'Reviewed cause-specific progression requires net-risk inputs, not an unreviewed cumulative-incidence surrogate.')
    require(tail in TAILS,'Explicit supported long-tail policy required.')
    rows=deepcopy(knots);require(isinstance(rows,list) and len(rows)>=2,'At least two cumulative-risk knots required.')
    rows=sorted(rows,key=lambda x:decimal(x['time_years']))
    times=[decimal(x['time_years']) for x in rows];risks=[decimal(x['cumulative_risk']) for x in rows]
    require(times[0]==0 and risks[0]==0 and all(b>a for a,b in zip(times,times[1:])),'Unique increasing times starting at zero required.')
    require(all(0<=r<1 for r in risks) and all(b>=a for a,b in zip(risks,risks[1:])),'Cumulative risks must be monotonic in [0,1).')
    hazards=[risk_to_hazard(r) for r in risks]
    with scientific_context():
        segments=[{'start':num(a),'end':num(b),'hazard_per_year':num((decimal(y)-decimal(x))/(b-a))}
            for a,b,x,y in zip(times,times[1:],hazards,hazards[1:])]
    pars=deepcopy(tail_parameters or {})
    if tail=='fixed_remote_hazard_v1':require(set(pars)=={'hazard_per_year'} and number(pars['hazard_per_year'])>=0,'Explicit fixed remote hazard required.')
    elif tail=='declining_to_remote_plateau_v1':
        require(set(pars)=={'plateau_hazard_per_year','decay_rate_per_year'} and 0<=number(pars['plateau_hazard_per_year'])<=number(segments[-1]['hazard_per_year'])
            and number(pars['decay_rate_per_year'])>0,'Explicit nonnegative plateau and positive decline parameter required.')
    else:require(not pars,'Unused tail parameters are not accepted.')
    return {'status':status,'risk_interpretation':risk_interpretation,'knots':[{'time_years':num(t),'cumulative_risk':num(r),'cumulative_hazard':h}
        for t,r,h in zip(times,risks,hazards)],'segments':segments,'tail':tail,'tail_parameters':{k:num(v) for k,v in pars.items()},
        'source':source,'provenance':provenance,'review':review,
        'interpolation':'piecewise_linear_cumulative_hazard_not_linear_probability',
        'tail_scientific_status':'not_extrapolated' if tail=='no_extrapolation_v1' else 'explicit_sensitivity_assumption_not_evidence_extension'}


def natural_history(group_curves,*,source,provenance,review,sex_curves=None):
    metadata(source,provenance,review)
    require(isinstance(group_curves,dict) and set(group_curves)=={g[2] for g in GROUPS},'All three age policies must be explicit, including external-required children.')
    def validate(curves):
        require(set(curves)==set(group_curves),'Sex curves must declare all three age policies.')
        for c in curves.values():
            if c is None:continue
            rebuilt=curve_specification(c['knots'],tail=c['tail'],tail_parameters=c['tail_parameters'],source=c['source'],
                provenance=c['provenance'],review=c['review'],status=c['status'],risk_interpretation=c['risk_interpretation'])
            require(rebuilt==c,'Curve does not reproduce from its source knots.')
    validate(group_curves)
    variants=deepcopy(sex_curves or {})
    require(set(variants)<= {'female','male'},'Unknown sex-specific natural-history policy.')
    for c in variants.values():validate(c)
    return artifact(NATURAL_HISTORY,group_curves=group_curves,sex_curves=variants,
        age_policy='attained_age_selects_curve_time_since_most_recent_event_age_at_event_retained_v1',
        age_groups=[{'age_lower':a,'age_upper_exclusive':b,'label':g,'status':'external_required' if group_curves[g] is None else group_curves[g]['status']}
                    for a,b,g in GROUPS],sex_effect='reviewed_sex_specific_only' if variants else 'common_baseline_curve_female_male_demography_retained',
        source=source,provenance=provenance,review=review)


def audited_adult_candidate(*,tail='no_extrapolation_v1',tail_parameters=None):
    review={'status':'reviewed','source_id':'M3B_source_audit','rationale':'Four published cumulative-risk knots verified; cause-specific-hazard transport not approved.'}
    provenance={'primary_xml_sha256':'167a755a17de17b641222adb9513814a1af72bbf9e0cbeefa30279bd9c558f3d',
        'diagnostic_only':True,'two_year_005_knot':'unverified_not_adopted','five_knot_candidate':'rejected_for_unqualified_adoption',
        'reported_95_uncertainty_intervals':{'1':['.032','.045'],'5':['.058','.077'],'10':['.064','.083'],'25':['.070','.089']},
        'mortality':'Paper specifies US mortality for remaining-lifetime results; four finite-horizon risks are not established here as portable cause-specific hazards.',
        'baseline_case_exclusion':'Infection-at-time-zero model output; not a raw baseline disease-free contact cohort; source cohorts require separate ascertainment review.'}
    adult=curve_specification([{'time_years':t,'cumulative_risk':r} for t,r in [('0','0'),('1','.038'),('5','.066'),('10','.072'),('25','.079')]],
        tail=tail,tail_parameters=tail_parameters,source=ADULT_SOURCE,provenance=provenance,review=review,status='diagnostic_adult_surrogate',
        risk_interpretation='diagnostic_cumulative_risk_surrogate_pending_cause_specific_review_v1')
    return natural_history({'0-4':None,'5-14':None,'15+':adult},source=ADULT_SOURCE,provenance=provenance,review=review)


def select_curve(history,age,sex):
    groups=history['sex_curves'].get(sex,history['group_curves'])
    label='0-4' if age<5 else '5-14' if age<15 else '15+'
    curve=groups[label]
    if curve is None:raise DataValidationError('external_required_natural_history_'+label)
    return curve


def curve_hazard(curve,tau):
    require(tau>=0,'Negative time since infection.')
    knots=[number(x['time_years']) for x in curve['knots']];last=knots[-1]
    if tau<last:return number(curve['segments'][max(0,bisect_right(knots,tau)-1)]['hazard_per_year'])
    tail=curve['tail'];r=number(curve['segments'][-1]['hazard_per_year'])
    if tail=='no_extrapolation_v1':raise DataValidationError('natural_history_tail_external_required')
    if tail=='fixed_remote_hazard_v1':return number(curve['tail_parameters']['hazard_per_year'])
    if tail=='declining_to_remote_plateau_v1':
        plateau=number(curve['tail_parameters']['plateau_hazard_per_year']);k=number(curve['tail_parameters']['decay_rate_per_year'])
        return plateau+(r-plateau)*exp(-k*(tau-last))
    return r


def curve_integral(curve,start,end):
    require(0<=start<=end,'Inverted progression interval.')
    last=number(curve['knots'][-1]['time_years'])
    if end>last and curve['tail']=='no_extrapolation_v1':raise DataValidationError('natural_history_tail_external_required')
    total=0.
    for s in curve['segments']:
        lo=max(start,number(s['start']));hi=min(end,number(s['end']))
        if hi>lo:total+=(hi-lo)*number(s['hazard_per_year'])
    if end>last:
        lo=max(start,last);duration=end-lo;tail=curve['tail'];r=number(curve['segments'][-1]['hazard_per_year'])
        if tail=='fixed_remote_hazard_v1':total+=duration*number(curve['tail_parameters']['hazard_per_year'])
        elif tail=='declining_to_remote_plateau_v1':
            p=number(curve['tail_parameters']['plateau_hazard_per_year']);k=number(curve['tail_parameters']['decay_rate_per_year'])
            total+=p*duration+(r-p)*exp(-k*(lo-last))*(-expm1(-k*duration))/k
        else:total+=r*duration
    return total


def age_breaks(age,duration):return [x-age for x in (5,15) if 0<x-age<duration]


def interval_hazard(history,*,tau,age_current,age_at_infection,sex,duration):
    require(tau>=0 and age_current>=0 and 0<=age_at_infection<=age_current and duration>=0,'Invalid age or clock interval.')
    points=[0.,*age_breaks(age_current,duration),duration];total=0.
    for lo,hi in zip(points,points[1:]):
        if hi>lo:total+=curve_integral(select_curve(history,age_current+lo,sex),tau+lo,tau+hi)
    return total


def cumulative_hazard(history,*,tau,age_at_infection,age_current,sex):
    require(abs(age_current-age_at_infection-tau)<1e-9,'Cumulative baseline hazard needs consistent actual infection age and elapsed time.')
    return interval_hazard(history,tau=0,age_current=age_at_infection,age_at_infection=age_at_infection,sex=sex,duration=tau)


def instantaneous_hazard(history,*,tau,age_current,age_at_infection,sex):
    require(0<=age_at_infection<=age_current and 0<=tau<=age_current,'Invalid infection clock.')
    return curve_hazard(select_curve(history,age_current,sex),tau)
