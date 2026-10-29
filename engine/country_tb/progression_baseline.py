"""Source-preserving cohort boundary between 3A allocation and forward progression."""
from collections import defaultdict
from copy import deepcopy
from math import exp,fsum,log

from engine.country_tb.progression_common import artifact,verify,require,metadata,number,num,quadrature,canonical
from engine.country_tb.risk_common import ratio,text_number
from engine.country_tb.data_common import profile_hash

COHORT='country_tb.progression_baseline_cohort.v1'


def baseline_cohort(rows,*,population,geography_id,baseline_time,source,provenance,review,diagnostic_only,
    active_TB_exclusion_policy,source_allocation=None,numerical_policy=None,recent_window_years='5'):
    metadata(source,provenance,review)
    require(type(diagnostic_only) is bool and isinstance(geography_id,str) and bool(geography_id),'Explicit diagnostic flag and residence geography required.')
    require(active_TB_exclusion_policy in {'reviewed_active_TB_exclusion_v1','assumed_zero_baseline_active_TB_diagnostic_v1'},'Explicit baseline active-TB risk-set policy required.')
    require(diagnostic_only or active_TB_exclusion_policy!='assumed_zero_baseline_active_TB_diagnostic_v1','Assumed absence of active TB is diagnostic only.')
    require(isinstance(rows,list) and rows and number(population)>0,'Positive, explicit cohort population required.')
    require(2024<=number(baseline_time)<2025,'Baseline time must be in 2024.')
    require(number(recent_window_years)>0,'Positive recent window required.')
    records=deepcopy(rows);ids=set()
    for r in records:
        required={'cohort_id','subgroup_id','sex','population','age','baseline_tbi_state','tau','effective_rr',
            'responsibility_weights','prior_remote_plus_recent_fraction','baseline_active_tb_population','effect_status',
            'pediatric_effect_approved','unresolved_reasons'}
        require(required<=r.keys(),'Incomplete forward cohort; exact ages and infection times cannot be invented.')
        require(isinstance(r['cohort_id'],str) and r['cohort_id'] and r['cohort_id'] not in ids,'Unique cohort IDs required.');ids.add(r['cohort_id'])
        require(r['sex'] in {'female','male'} and r['subgroup_id'] and r['baseline_tbi_state'] in {'never_infected','recent','remote_only','unresolved'},'Invalid cohort identities/state.')
        n=number(r['population']);active=number(r['baseline_active_tb_population']);require(0<=active<=n,'Baseline active TB exceeds cohort population.')
        require(isinstance(r['unresolved_reasons'],list) and all(isinstance(x,str) for x in r['unresolved_reasons']),'Explicit unresolved reasons required.')
        age=None if r['age'] is None else number(r['age']);tau=None if r['tau'] is None else number(r['tau'])
        rr=None if r['effective_rr'] is None else number(r['effective_rr']);prior=number(r['prior_remote_plus_recent_fraction'],probability=True)
        require(age is None or age>=0,'Negative age.');require(rr is None or rr>=1,'Progression risk multiplier must be at least one.')
        if tau is not None:require(age is not None and 0<=tau<=age,'Infection time cannot precede birth.')
        if r['baseline_tbi_state']=='never_infected':require(tau is None and prior==0,'Never-infected has no infection clock or prior remote subset.')
        if tau is not None and r['baseline_tbi_state']=='recent':require(tau<number(recent_window_years),'Recent infection lies outside the declared recent window.')
        if tau is not None and r['baseline_tbi_state']=='remote_only':require(tau>=number(recent_window_years),'Remote-only infection lies inside the recent window.')
        if r['baseline_tbi_state']!='recent':require(prior==0,'Prior-remote-plus-recent is only a subset of recent.')
        if age is not None and age<15 and rr is not None and rr>1:
            require(r['pediatric_effect_approved'] and r['effect_status']=='reviewed_pediatric','Adult risk effects cannot be applied to children.')
        scope=r.get('effect_age_applicability')
        if scope is None:
            require(rr in (None,1) or (diagnostic_only and age is not None and age>=15),
                'Reviewed nonbaseline or pediatric RR needs explicit follow-up age applicability.')
            scope={'age_lower':'15' if rr is not None and rr>1 else '0','age_upper_exclusive':None,
                'policy':'diagnostic_broad_adult_effect_scope_v1' if rr is not None and rr>1 else 'no_nonbaseline_effect_v1'}
        require(isinstance(scope,dict) and {'age_lower','age_upper_exclusive','policy'}<=scope.keys(),'Explicit effect age scope required.')
        lower=number(scope['age_lower']);upper=None if scope['age_upper_exclusive'] is None else number(scope['age_upper_exclusive'])
        require(lower>=0 and (upper is None or upper>lower) and (age is None or lower<=age and (upper is None or age<upper)),
            'Baseline age lies outside the reviewed effect applicability.')
        r['effect_age_applicability']={**scope,'age_lower':num(lower),'age_upper_exclusive':None if upper is None else num(upper)}
        require(isinstance(r['responsibility_weights'],dict) and all(0<=number(v,probability=True)<=1 for v in r['responsibility_weights'].values()),'Invalid highest-RR responsibility weights.')
        if not r['unresolved_reasons']:
            require(age is not None and rr is not None and (tau is not None or r['baseline_tbi_state']=='never_infected'),'Missing age, clock or RR must be visibly unresolved.')
        if not diagnostic_only:require(r['effect_status'] in {'reviewed','reviewed_pediatric','no_factors'},'Unreviewed effects require diagnostic-only allocation.')
        for k in ('population','baseline_active_tb_population','prior_remote_plus_recent_fraction'):r[k]=num(r[k])
        for k in ('age','tau','effective_rr'):
            if r[k] is not None:r[k]=num(r[k])
        r['age_at_most_recent_infection']=num(age-tau) if age is not None and tau is not None else None
        r['responsibility_weights']={k:num(v) for k,v in sorted(r['responsibility_weights'].items())}
    achieved=fsum(number(r['population']) for r in records)
    require(abs(achieved-number(population))<=max(1e-9,number(population)*1e-12),'Cohort population reconciliation failed.')
    return artifact(COHORT,geography_id=geography_id,baseline_time=num(baseline_time),baseline_year=2024,population=num(population),recent_window_years=num(recent_window_years),
        rows=sorted(records,key=lambda r:r['cohort_id']),source_allocation=source_allocation,
        active_TB_exclusion_policy=active_TB_exclusion_policy,baseline_active_TB_excluded=num(fsum(number(r['baseline_active_tb_population']) for r in records)),
        prior_TB_history='not_represented_treated_prior_TB_may_remain_TBI_positive_historical_depletion_unresolved',
        risk_factor_followup='baseline_memberships_and_effects_fixed_only_within_declared_age_scope_no_future_effect_fallback_or_incident_comorbidity_model',
        numerical_policy=numerical_policy or {'age_time_support':'explicit_reviewed_discrete_support_no_interpolation'},
        diagnostic_only=diagnostic_only,source=source,provenance=provenance,review=review)


def _integral(curve,lo,hi):
    return fsum(number(r['shape'])*(min(hi,number(r['end']))-max(lo,number(r['start'])))
        for r in curve if min(hi,number(r['end']))>max(lo,number(r['start'])))


def _support(model,cell,state,*,max_width):
    """Joint attained-age / last-event density, not a mean infection time."""
    from engine.country_tb.catalytic_core import history_curve
    history=model['calibration']['pressure_history'];alpha=number(model['calibration']['alpha'])
    curve=history_curve(history,cell['sex']);time=number(model['baseline_time']);window=number(model['recent_window_years'])
    lo=cell['age_lower'];hi=cell['age_upper_exclusive'];source=next(r for r in model['rows'] if
        r['sex']==cell['sex'] and r['age_lower']==lo and r['age_upper_exclusive']==hi)
    if hi is None:
        tails=model['configuration'].get('tail_distributions',{});tail=tails.get(cell['sex']+':'+str(lo))
        if tail is None:return None
        total=fsum(number(x['weight']) for x in tail['ages']);ages=[(number(x['age']),number(x['weight'])/total) for x in tail['ages']]
    elif source['age_policy']=='exact_age':ages=[(float(lo),1.)]
    elif source['age_policy'] in {'uniform_within_single_year','reviewed_uniform_within_interval'}:ages=[(a,w/(hi-lo)) for a,w in quadrature(float(lo),float(hi))]
    else:return None
    raw=[]
    for age,weight in ages:
        birth=time-age
        if state=='never_infected':raw.append((age,None,log(weight)-alpha*_integral(curve,birth,time)));continue
        start=max(birth,time-window) if state=='recent' else birth
        end=time if state=='recent' else min(time,time-window)
        if end<=start:continue
        for r in curve:
            a=max(start,number(r['start']));b=min(end,number(r['end']));rate=alpha*number(r['shape'])
            if b<=a or rate==0:continue
            boundaries={a,b,*[time-t for t in (1,5,10,25) if a<time-t<b]}
            points=sorted(boundaries)
            for x,y in zip(points,points[1:]):
                count=max(1,int((y-x)/max_width+.999999999));length=(y-x)/count
                for j in range(count):
                    for event,w in quadrature(x+j*length,x+(j+1)*length):
                        raw.append((age,time-event,log(weight*w*rate)-alpha*_integral(curve,event,time)))
    if not raw:return None
    shift=max(x[2] for x in raw);weights=[exp(x[2]-shift) for x in raw];total=fsum(weights)
    return [(a,t,w/total) for (a,t,_),w in zip(raw,weights)]


def cohort_from_allocation(allocation,*,source,provenance,review,active_TB_exclusion_policy,
    conditional_support=None,max_time_interval_years='1'):
    """Exact 3A state/RR totals, with explicit joint age/time integration support.

    Supplied direct profiles require reviewed conditional age/time nodes. A TBI
    prevalence or a mean/quantile is insufficient to reconstruct a distribution.
    """
    from engine.country_tb.risk_common import verify as verify_risk
    from engine.country_tb.risk_allocation import ALLOCATION,MIXTURE
    verify_risk(allocation,allocation['schema_id']);require(allocation['schema_id'] in {ALLOCATION,MIXTURE},'Expected 3A risk allocation.')
    width=number(max_time_interval_years);require(0<width<=1,'Explicit bounded infection-time quadrature interval required.')
    groups=allocation.get('subgroups',[allocation]);rows=[];support=conditional_support or {}
    for group in groups:
        base=group['resolved_baseline'];model=base.get('model_derived_profile')
        for index,cell in enumerate(group['rows']):
            n=float(ratio(cell['population']));prefix=group['subgroup_id']+':'+cell['sex']+':'+str(cell['age_lower'])
            if not n:continue
            common={'subgroup_id':group['subgroup_id'],'sex':cell['sex'],'baseline_active_tb_population':'0','pediatric_effect_approved':False}
            def unresolved(reasons):
                rows.append({**common,'cohort_id':prefix+':unresolved','population':num(n),'age':None,'tau':None,'effective_rr':None,
                    'baseline_tbi_state':'unresolved','responsibility_weights':{},'effect_status':'unresolved',
                    'prior_remote_plus_recent_fraction':'0','unresolved_reasons':reasons,
                    'source_age_lower':cell['age_lower'],'source_age_upper_exclusive':cell['age_upper_exclusive']})
            if not cell['complete_expected_stratification']:
                unresolved(['baseline_TBI_or_risk_allocation_unresolved']);continue
            # Aggregate over count categories: m and any-risk determine hazard and TBI coupling.
            masses=defaultdict(float)
            for c in cell['compressed_cells']:
                for state in ('never_infected','recent','remote_only'):
                    mass=float(ratio(c['tbi_state_counts'][state]))
                    if mass:masses[number(c['effective_rr']['value']),c['any_non_recency_risk_factor'],state]+=mass
            compiled=[];failed=False
            for state in ('never_infected','recent','remote_only'):
                if not any(key[2]==state for key in masses):continue
                key=prefix+':'+state
                if key in support:
                    supplied=support[key];metadata(supplied['source'],supplied['provenance'],supplied['review'])
                    nodes=[(number(x['age']),None if x['tau'] is None else number(x['tau']),number(x['weight'])) for x in supplied['nodes']]
                    require(nodes and all(w>0 for _,_,w in nodes) and abs(fsum(w for _,_,w in nodes)-1)<1e-12,'Reviewed conditional support must sum to one.')
                    require(all(cell['age_lower']<=a and (cell['age_upper_exclusive'] is None or a<cell['age_upper_exclusive'])
                        and ((t is None and state=='never_infected') or (t is not None and 0<=t<=a and state!='never_infected')) for a,t,_ in nodes),'Conditional support lies outside age/state bounds.')
                elif model is not None and base['calibration'] is not None:nodes=_support(model,cell,state,max_width=width)
                else:nodes=None
                if nodes is None:failed=True;break
                for (rr,risk,s),mass in sorted(masses.items()):
                    if s!=state:continue
                    denominator=sum((ratio(c['probability']) for c in cell['compressed_cells']
                        if number(c['effective_rr']['value'])==rr and c['any_non_recency_risk_factor']==risk),ratio('0'))
                    attribution={x['factor_id']:text_number(ratio(x['probability'])/denominator) for x in cell['highest_rr_contributor_probabilities']
                        if risk and x['probability'] is not None and float(ratio(x['probability']))>0 and
                        next(f for f in cell['factor_marginals'] if f['factor_id']==x['factor_id'])['effect']['rr']==num(rr)}
                    recent=cell['tbi_state_prevalences']['recent'];prior=cell['tbi_state_prevalences']['prior_remote_plus_recent']
                    prior_fraction=float(ratio(prior)/ratio(recent)) if state=='recent' and ratio(recent)>0 else 0.
                    ranges=[]
                    if risk:
                        for f in cell['factor_marginals']:
                            if f['prevalence'] is None or ratio(f['prevalence'])==0:continue
                            catalogue=next(x for x in group['catalogue']['rows'] if x['factor_id']==f['factor_id'])
                            matching=[x for x in catalogue['reviewed_effects'] if x['sex'] in {'both',cell['sex']} and
                                x['age_lower']<=cell['age_lower'] and (x['age_upper_exclusive'] is None or
                                cell['age_upper_exclusive'] is not None and cell['age_upper_exclusive']<=x['age_upper_exclusive'])]
                            ranges.extend(matching or [catalogue['age_applicability']])
                    lower=max([0,*[x['age_lower'] for x in ranges]])
                    upper=min([x['age_upper_exclusive'] for x in ranges if x['age_upper_exclusive'] is not None],default=None)
                    scope={'age_lower':num(lower),'age_upper_exclusive':None if upper is None else num(upper),
                        'policy':'intersection_of_included_positive_factor_baseline_effect_scopes_no_future_effect_fallback_v1'}
                    for j,(age,tau,w) in enumerate(nodes):
                        pediatric=age<15 and bool(ranges) and all(x.get('review',{}).get('status')=='reviewed' for x in ranges)
                        compiled.append({**common,'cohort_id':prefix+':'+state+':'+num(rr)+':'+str(risk)+':'+str(j),
                            'population':num(mass*w),'age':num(age),'tau':None if tau is None else num(tau),'baseline_tbi_state':state,
                            'effective_rr':num(rr),'responsibility_weights':attribution,
                            'effect_status':'reviewed_pediatric' if pediatric else 'diagnostic_provisional' if group['diagnostic_only'] else 'reviewed',
                            'pediatric_effect_approved':pediatric,'prior_remote_plus_recent_fraction':num(prior_fraction),
                            'unresolved_reasons':[],'effect_age_applicability':scope})
            if failed:unresolved(['conditional_age_and_most_recent_infection_distribution_external_required'])
            else:rows.extend(compiled)
    times={g['resolved_baseline']['selection']['configuration'].get('baseline_time','2024.5') for g in groups}
    require(len(times)==1,'Subgroups need a common baseline time for projection.')
    windows={g['resolved_baseline']['selection']['configuration'].get('recent_window_years','5') for g in groups}
    require(len(windows)==1,'Subgroups must use the same explicit recent window in this forward contract.')
    return baseline_cohort(rows,population=allocation['summary']['population']['value'],geography_id=allocation['geography_id'],
        baseline_time=next(iter(times)),source=source,provenance=provenance,review=review,
        diagnostic_only=allocation['diagnostic_only'] or active_TB_exclusion_policy=='assumed_zero_baseline_active_TB_diagnostic_v1',
        active_TB_exclusion_policy=active_TB_exclusion_policy,source_allocation=allocation,
        recent_window_years=next(iter(windows)),
        numerical_policy={'age_time_support':'positive_Gauss_Legendre_8_joint_age_and_last_event_density_within_state_normalized_to_unchanged_3A_state_RR_counts',
            'maximum_time_interval_years':num(width),'support_source_hashes':[g['resolved_baseline_hash'] for g in sorted(groups,key=lambda g:g['subgroup_id'])],
            'reviewed_direct_support':support,'prior_remote_audit':'source_subset_fraction_retained_within_recent_not_an_extra_state_or_hazard',
            'precision':'binary64_full_precision_no_display_rounding_grid_error_requires_refinement'})
