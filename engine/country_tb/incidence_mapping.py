"""One historical Poisson scale mapped through the frozen 3B forward kernels.

Only exact ages or externally expanded reviewed age supports enter this layer.
Piecewise-linear response kernels are integrated against exact exponential
last-event interval masses. Refinement is independent of the alpha search grid.
"""
from bisect import bisect_right
from copy import deepcopy
from math import exp, expm1, fsum, log1p

from engine.country_tb.data_common import canonical, profile_hash, seal, DataValidationError
from engine.country_tb.composition_common import require
from engine.country_tb.progression_common import number, num, verify as verify_forward
from engine.country_tb.progression_forward import OPERATOR, _trajectory
from engine.country_tb.progression_natural_history import select_curve, instantaneous_hazard
from engine.country_tb.pressure_shape import validate_history
from engine.country_tb.catalytic_core import history_curve
from engine.country_tb.incidence_inversion import provenance
from engine.country_tb.risk_factor_catalog import FACTOR_IDS


def _event_mass(alpha, q0, q1):
    return exp(-alpha*q0)*(-expm1(-alpha*(q1-q0)))


def _mean_offset(rate, width):
    z=rate*width
    if z < 1e-4:return width*(.5-z/12+z**3/720)
    if z > 700:return 1/rate
    return 1/rate-width/expm1(z)


def _association_cells(p,q,c):
    """Convex mixture of nonnegative Frechet tables, avoiding subtraction of J.

    Computing p-J then dividing by1-q can produce negative masses near q=1.
    Explicit boundary cells retain their small complements without clipping.
    Order: TBI+risk,TBI only,risk only,neither.
    """
    independent=(p*q,p*(1-q),(1-p)*q,(1-p)*(1-q))
    if c<0:
        boundary=(0.,p,q,1-(p+q)) if p+q<=1 else (p+q-1,1-q,1-p,0.)
    else:
        boundary=(p,0.,q-p,1-q) if p<=q else (q,p-q,0.,1-p)
    amount=abs(c)
    cells=tuple((1-amount)*x+amount*y for x,y in zip(independent,boundary))
    require(all(v>=0 for v in cells),'Association produced a negative mass; no clipping.')
    return cells


class IncidenceMapping:
    """Immutable scenario with deterministic internal caches, no fitted nuisance parameters.

    Future pressure is a separately supplied fixed scenario. Historical alpha
    never changes it. Both baseline instantaneous and prospective person-time
    incidence are supported; these observation measures are not interchangeable.
    """
    def __init__(self, *, history, strata, operator, measure, horizon_years, population_kind,
        geography_id, subgroup_id, sex_scope, age_lower, age_upper_exclusive,
        population_history, source, review, provenance_record, kernel_width=.25):
        provenance(source,review,provenance_record);validate_history(history);verify_forward(operator,OPERATOR)
        require(measure in {'instantaneous_baseline_per_100000_person_years','prospective_per_100000_person_years'}, 'Explicit incidence measure required.')
        require(population_kind in {'synthetic_closed_population','diagnostic_country_population','reviewed_population'}, 'Explicit population evidence kind required.')
        require(geography_id==history['geography_id']==operator['competing_mortality']['geography_id'], 'History, mortality and population geography mismatch.')
        require(sex_scope in {'both','female','male'} and type(age_lower) is int and age_lower>=0 and
            (age_upper_exclusive is None or type(age_upper_exclusive) is int and age_upper_exclusive>age_lower), 'Explicit population age/sex scope required.')
        require(0 < number(horizon_years) <= 20 and 0 < kernel_width <= 1, 'Bounded positive horizon and response grid required.')
        require(number(horizon_years) in {number(h) for h in operator['horizons']},'Mapping horizon must be declared in the bounded forward operator.')
        require(number(operator['future_pressure']['horizon_years'])>=number(horizon_years),'Future pressure does not cover mapping horizon.')
        require(isinstance(strata,list) and strata and len(strata)<=1000,'Explicit population supports, at most 1000, required.')
        rows=deepcopy(strata);ids=set();self.unreviewed_risk_factors=set()
        for r in rows:
            require(set(r)>={'stratum_id','age','sex','subgroup_id','population','factors','association_control','age_support_policy'}, 'Incomplete exact-age/sex stratum.')
            require(r['stratum_id'] not in ids and r['stratum_id'],'Unique stable stratum IDs required.');ids.add(r['stratum_id'])
            require(r['sex'] in {'female','male'} and (sex_scope=='both' or r['sex']==sex_scope),'Population sex mismatch.')
            require(number(r['population'])>0 and -1<=number(r['association_control'])<=1,'Population/control outside bounds.')
            require(r['age_support_policy'] in {'exact_age_supplied','reviewed_within_cell_age_distribution','unresolved_open_age'}, 'No implicit demographic midpoint.')
            if r['age'] is not None:
                a=number(r['age']);require(a>=age_lower and (age_upper_exclusive is None or a<age_upper_exclusive),'Age outside target population.')
                require(r['age_support_policy']!='unresolved_open_age','Unresolved age cannot supply an invented age.')
                if r['age_support_policy']=='reviewed_within_cell_age_distribution':require(r.get('age_distribution_review'),'Reviewed within-cell support provenance required.')
            else:require(r['age_support_policy']=='unresolved_open_age','Unknown age must be visibly unresolved.')
            require(isinstance(r['factors'],list) and len(r['factors'])<=11,'Bounded explicit factor inputs required.')
            fids=set()
            for f in r['factors']:
                require(set(f)>={'factor_id','prevalence','rr','effect_status','definition_compatibility','source'},'Incomplete factor prevalence/effect provenance.')
                require(f['factor_id'] not in fids and f['factor_id'] in FACTOR_IDS,'Duplicate/unknown/recency factor.');fids.add(f['factor_id'])
                require(0<=number(f['prevalence'])<=1 and number(f['rr'])>=1,'Invalid prevalence/RR; no clipping.')
                require(f['definition_compatibility'] in {'exact','reviewed_mapping'},'Definition mismatch cannot enter a mapping.')
                if f['definition_compatibility']=='reviewed_mapping':require(f.get('mapping_review') and f.get('mapping_provenance'),'Reviewed mapping provenance required.')
                require(f['effect_status'] in {'reviewed','synthetic_mathematical','legacy_sensitivity_only'},'Unknown effect status.')
                if f['effect_status']=='reviewed':
                    er=f.get('effect_review',{})
                    require(er.get('status')=='reviewed' and all(er.get(k) for k in
                        ('source_id','provenance','effect_measure','reference_group','latent_infection_definition','treatment_status','age_applicability','horizon')),
                        'Reviewed effect needs its conditional-TBI estimand, applicability and provenance, not a status label alone.')
                require(f['effect_status']!='synthetic_mathematical' or population_kind=='synthetic_closed_population','Synthetic effect cannot authorize real population.')
                require(r['age'] is None or number(r['age'])>=15 or f.get('pediatric_review'),'Adult factors cannot silently apply below 15.')
                if f['effect_status']=='legacy_sensitivity_only' and number(f['prevalence'])>0:self.unreviewed_risk_factors.add(f['factor_id'])
            r['factors']=sorted(r['factors'],key=lambda f:f['factor_id'])
        rows=sorted(rows,key=canonical)
        self.operator=deepcopy(operator);self.history=deepcopy(history);self.rows=rows;self.width=kernel_width
        self.measure=measure;self.horizon=number(horizon_years);self.compiled=[];self.material_curve_gaps=[];self.tail_gaps=[]
        self.specification=seal(dict(schema_id='country_tb.incidence_mapping_scenario.internal.v1',
            geography_id=geography_id,subgroup_id=subgroup_id,sex_scope=sex_scope,age_lower=age_lower,
            age_upper_exclusive=age_upper_exclusive,baseline_year=2024,population_kind=population_kind,
            history_hash=history['profile_hash'],operator_hash=operator['profile_hash'],rows=rows,
            population_history=deepcopy(population_history),measure=measure,horizon_years=self.horizon,
            future_pressure_policy='fixed_scenario_independent_of_historical_scale_no_automatic_extrapolation',
            kernel_width=kernel_width,source=source,review=review,provenance=provenance_record))
        self._constructor=deepcopy(dict(history=history,strata=strata,operator=operator,measure=measure,horizon_years=horizon_years,
            population_kind=population_kind,geography_id=geography_id,subgroup_id=subgroup_id,sex_scope=sex_scope,
            age_lower=age_lower,age_upper_exclusive=age_upper_exclusive,population_history=population_history,
            source=source,review=review,provenance_record=provenance_record))
        for row in rows:self.compiled.append(self._compile(row))
        self.population=fsum(number(r['population']) for r in rows)
        self.maximum_unresolved_share=fsum(number(r['population']) for r,c in zip(rows,self.compiled) if c['gaps'])/self.population
        self._values={}

    def refined(self):
        arguments=deepcopy(self._constructor)
        if self.operator['mode']=='exogenous_future_infection_v1':
            step=number(self.operator['step_years'])/2
            require(self.horizon/step<=10000,'numerical_refinement_budget_exceeded')
            arguments['operator']['step_years']=num(step)
            arguments['operator']=seal(arguments['operator'])
        return IncidenceMapping(**arguments,kernel_width=self.width/2)

    def total_tbi(self,alpha):
        if any('total_exposure' not in c for c in self.compiled):return None
        return fsum(number(r['population'])*(-expm1(-number(alpha)*c['total_exposure']))
            for r,c in zip(self.rows,self.compiled))/self.population

    def _compile(self,r):
        if r['age'] is None:return {'gaps':['unresolved_open_age'],'patterns':[]}
        age=number(r['age']);sex=r['sex'];gaps=[]
        try:curve=select_curve(self.operator['natural_history'],age,sex)
        except DataValidationError as exc:
            self.material_curve_gaps.append(r['stratum_id']);return {'gaps':[str(exc)],'patterns':[]}
        if not self.measure.startswith('instantaneous'):
            for attained in (age+self.horizon, *[b for b in (5,15) if age<b<age+self.horizon]):
                try:select_curve(self.operator['natural_history'],attained,sex)
                except DataValidationError as exc:
                    self.material_curve_gaps.append(r['stratum_id']);return {'gaps':[str(exc)],'patterns':[]}
        if curve['tail']=='no_extrapolation_v1' and age+(0 if self.measure.startswith('instantaneous') else self.horizon)>number(curve['knots'][-1]['time_years']):
            self.tail_gaps.append(r['stratum_id'])
        raw=history_curve(self.history,sex);intervals=[(number(x['start']),number(x['end']),number(x['shape'])) for x in raw]
        require(intervals[0][0]<=2024-age and intervals[-1][1]>=2024,'Complete lifetime history required.')
        points={0.,age,min(age,number(self.operator['recent_window_years']))}
        points.update(2024-x for a,b,q in intervals for x in (a,b) if 0<2024-x<age)
        for c in [v for cs in [self.operator['natural_history']['group_curves'],*self.operator['natural_history']['sex_curves'].values()] for v in cs.values() if v]:
            for knot in c['knots']:
                t=number(knot['time_years'])
                points.update(x for x in (t,t-self.horizon) if 0<x<age)
        cuts=sorted(points);points=set(cuts)
        from math import ceil
        for a,b in zip(cuts,cuts[1:]):
            n=max(1,ceil((b-a)/self.width));points.update(a+(b-a)*j/n for j in range(1,n))
        points=sorted(points);bins=[];cumulative=0.
        for a,b in zip(points,points[1:]):
            q=next(q for lo,hi,q in intervals if lo<=2024-(a+b)/2<hi)
            bins.append((a,b,q,cumulative,cumulative+q*(b-a)));cumulative+=q*(b-a)
        patterns=[{'factor_ids':[],'probability':1.,'rr':1.}]
        for factor in r['factors']:
            p=number(factor['prevalence']);rr=number(factor['rr']);new=[]
            for x in patterns:
                if p<1:new.append(dict(x,probability=x['probability']*(1-p)))
                if p>0:new.append(dict(factor_ids=x['factor_ids']+[factor['factor_id']],probability=x['probability']*p,rr=max(x['rr'],rr)))
            patterns=new
        ps=[number(f['prevalence']) for f in r['factors']]
        q=1. if any(p==1 for p in ps) else -expm1(fsum(log1p(-p) for p in ps))
        # Scalar scans need only RR and union membership. This exact compression
        # is not used for individual-factor audit outputs.
        compressed={}
        for x in patterns:compressed.setdefault((x['rr'],bool(x['factor_ids'])),[]).append(x['probability'])
        fastpatterns=[dict(rr=rr,factor_ids=['__union__'] if positive else [],probability=fsum(masses))
            for (rr,positive),masses in sorted(compressed.items())]
        kernels={};cache={}
        for rr in sorted({x['rr'] for x in patterns}):
            def response(tau):
                row=dict(age=age,sex=sex,subgroup_id=r['subgroup_id'],effective_rr=rr,tau=tau)
                try:
                    if self.measure.startswith('instantaneous'):
                        return (0. if tau is None else rr*instantaneous_hazard(self.operator['natural_history'],tau=tau,age_current=age,age_at_infection=age-tau,sex=sex),1.)
                    value=_trajectory(self.operator,row,self.horizon,cache)
                    return value['cases'],value['person_years']
                except DataValidationError as exc:
                    if 'external_required' not in str(exc):raise
                    return None
            values=[]
            for a,b,shape,qa,qb in bins:
                # Endpoints are one-sided so a hazard knot never contaminates its neighbour.
                eps=min(1e-9,(b-a)*1e-7)
                values.append((response(a+eps),response(b-eps)))
            never=response(None)
            if never is None:gaps.append('mortality_or_future_coverage_external_required')
            kernels[rr]=(values,never)
        return dict(gaps=sorted(set(gaps)),patterns=patterns,fastpatterns=fastpatterns,union=q,bins=bins,kernels=kernels,total_exposure=cumulative)

    def incidence(self,alpha):return self.evaluate(alpha,detailed=False)['incidence']

    def evaluate(self,alpha,*,detailed=True):
        alpha=number(alpha);require(alpha>=0,'Scale must be non-negative; no clipping.')
        if detailed and alpha in self._values:return deepcopy(self._values[alpha])
        results=[];parts=[];total_num=total_den=unresolved=0.
        for r,c in zip(self.rows,self.compiled):
            pop=number(r['population'])
            if c['gaps']:
                unresolved+=pop
                if detailed:results.append(dict(stratum_id=r['stratum_id'],population=pop,unresolved_population=pop,reasons=c['gaps'],any_tbi=None))
                continue
            p=-expm1(-alpha*c['total_exposure']);q=c['union'];control=number(r['association_control'])
            j0=p*q;cells=_association_cells(p,q,control);j=cells[0]
            recent=remote=prior=time_total=unknown=0.;numerator=denominator=0.;factor_tbi={}
            mass_bins=[_event_mass(alpha,qa,qb) for a,b,shape,qa,qb in c['bins']]
            w=number(self.operator['recent_window_years']);recent_exposure=fsum(qb-qa for a,b,shape,qa,qb in c['bins'] if b<=w+1e-10)
            recent=-expm1(-alpha*recent_exposure);remote=p-recent
            prior=(-expm1(-alpha*(c['total_exposure']-recent_exposure)))*recent
            for pattern in c['patterns'] if detailed else c['fastpatterns']:
                mass=pattern['probability'];positive=bool(pattern['factor_ids']);rr=pattern['rr']
                tbi_cond=(cells[0]/q if positive and q else cells[1]/(1-q) if not positive and q<1 else 0.)
                never_cond=(cells[2]/q if positive and q else cells[3]/(1-q) if not positive and q<1 else 0.)
                ti_mass=mass*tbi_cond;never_mass=mass*never_cond;v,never=c['kernels'][rr]
                responsible=[f['factor_id'] for f in r['factors'] if f['factor_id'] in pattern['factor_ids'] and number(f['rr'])==rr]
                def part(weight,kernel,state,tau):
                    if detailed and weight>0:
                        parts.append(dict(stratum_id=r['stratum_id'],subgroup_id=r['subgroup_id'],effective_rr=rr,
                            factor_ids=pattern['factor_ids'],responsible_factors=responsible,tbi_state=state,tau=tau,
                            population=pop*weight,numerator=pop*weight*kernel[0],person_years=pop*weight*kernel[1],expected_risk=kernel[0]))
                n=d=u=0.
                for idx,((a,b,shape,qa,qb),event_mass) in enumerate(zip(c['bins'],mass_bins)):
                    effective_mass=ti_mass*event_mass/p if p else 0.
                    if effective_mass==0:continue
                    pair=v[idx]
                    if pair[0] is None or pair[1] is None:u+=effective_mass;continue
                    offset=_mean_offset(alpha*shape,b-a) if alpha*shape else (b-a)/2
                    f=offset/(b-a)
                    response=(pair[0][0]+f*(pair[1][0]-pair[0][0]),pair[0][1]+f*(pair[1][1]-pair[0][1]))
                    n+=effective_mass*response[0];d+=effective_mass*response[1]
                    part(effective_mass,response,'recent' if b<=w+1e-10 else 'remote_only',a+offset)
                if never is None:u+=never_mass
                else:
                    n+=never_mass*never[0];d+=never_mass*never[1]
                    part(never_mass,never,'never_infected',None)
                numerator+=n;denominator+=d;unknown+=u
                for fid in pattern['factor_ids']:factor_tbi[fid]=factor_tbi.get(fid,0.)+ti_mass
            if detailed:
                for (a,b,shape,qa,qb),event_mass in zip(c['bins'],mass_bins):
                    offset=_mean_offset(alpha*shape,b-a) if alpha*shape else (b-a)/2
                    time_total+=event_mass*(a+offset)
            phi_den=(p*(1-p)*q*(1-q))**.5
            results.append(dict(stratum_id=r['stratum_id'],age=r['age'],sex=r['sex'],subgroup_id=r['subgroup_id'],population=pop,
                any_tbi=p,recent=recent,remote_only=remote,never_infected=1-p,prior_remote_plus_recent=prior,
                expected_time_since_most_recent_infection=None if p==0 else time_total/p,
                any_risk=q,tbi_and_any_risk=j,association_control=control,achieved_phi=None if phi_den==0 else (j-j0)/phi_den,
                factor_marginals={f['factor_id']:number(f['prevalence']) for f in r['factors']},factor_tbi_joint=factor_tbi,
                unresolved_population=pop*unknown,covered_numerator=pop*numerator,covered_person_years=pop*denominator))
            total_num+=pop*numerator;total_den+=pop*denominator;unresolved+=pop*unknown
        if not detailed:return {'incidence':None if unresolved>0 or total_den==0 else 100000*total_num/total_den}
        # Tiny subnormal unresolved mass still prevents a complete estimate; no clipping threshold.
        value=seal(dict(schema_id='country_tb.incidence_mapping_evaluation.internal.v1',scenario_hash=self.specification['profile_hash'],
            alpha=alpha,measure=self.measure,population=self.population,unresolved_population=unresolved,
            incidence=None if unresolved>0 or total_den==0 else 100000*total_num/total_den,
            covered_numerator=total_num,covered_person_years=total_den,
            any_tbi=None if any(x.get('any_tbi') is None for x in results) else fsum(x['population']*x['any_tbi'] for x in results)/self.population,
            recent=None if any(x.get('any_tbi') is None for x in results) else fsum(x['population']*x['recent'] for x in results)/self.population,
            remote_only=None if any(x.get('any_tbi') is None for x in results) else fsum(x['population']*x['remote_only'] for x in results)/self.population,
            rows=results,risk_parts=sorted(parts,key=canonical),production_available=False))
        # Dense scans retain scalars, not thousands of full population decompositions.
        if len(self._values)<8:self._values[alpha]=value
        return deepcopy(value)

    def concentration(self,alpha,top_share):
        require(0<top_share<=1,'Top population share outside (0,1].')
        value=self.evaluate(alpha);rows=sorted(value['risk_parts'],key=lambda r:(-(r['expected_risk'] or 0),canonical(r)))
        budget=self.population*top_share;cases=0.
        for r in rows:
            take=min(budget,r['population']);cases+=take*(r['expected_risk'] or 0);budget-=take
            if budget<=0:break
        total=value['covered_numerator']
        return {'top_population_share':top_share,'case_share':None if total==0 else cases/total,
            'ranking':'expected_first_episode_risk_or_instantaneous_rate_with_fractional_tie_allocation',
            'unresolved_population':value['unresolved_population']}
