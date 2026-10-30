"""Independent one-year numerical audit, not a replacement forward operator.

Within one year after a future event the3B adult surrogate's first segment is
constant. The resulting two transient-state system has a closed-form solution.
This deliberately limited reference diagnoses interpolation/time-grid error.
"""
from math import exp,expm1,fsum,ceil

from engine.country_tb.progression_common import number,quadrature
from engine.country_tb.progression_natural_history import curve_hazard
from engine.country_tb.catalytic_core import history_curve
from engine.country_tb.validation_data import require


def _phi(rate,duration):return duration if rate==0 else -expm1(-rate*duration)/rate


def one_year_kernel(curve,tau,*,rr=1.,mortality_rate=0.,future_pressure=0.):
    """B:original clock/uninfected; F:at least one prospective event.

    B'=-[h_B+mu+lambda]B; F'=lambda*B-[h_F+mu]F.
    Repeated events within F do not change its constant first-year hazard.
    First TB/death are absorbing. All rates are cause-specific mathematical
    assumptions here, not an approval of surrogate clinical interpretation.
    """
    require(rr>=1 and mortality_rate>=0 and future_pressure>=0,'Nonnegative rates and RR>=1 required.')
    require(number(curve['segments'][0]['start'])==0 and number(curve['segments'][0]['end'])>=1,
        'One-year reference requires a constant first-year segment.')
    points={0.,1.}
    if tau is not None:
        require(tau>=0,'Negative infection clock.')
        points.update(number(k['time_years'])-tau for k in curve['knots'] if 0<number(k['time_years'])-tau<1)
    B=1.;F=cases=deaths=years=0.;fresh_h=rr*curve_hazard(curve,0.)
    for lo,hi in zip(sorted(points),sorted(points)[1:]):
        d=hi-lo;h=0. if tau is None else rr*curve_hazard(curve,tau+(lo+hi)/2)
        b=h+mortality_rate+future_pressure;r=fresh_h+mortality_rate;gap=r-b
        pB=B*_phi(b,d)
        if abs(gap*d)<1e-5:
            generated_years=fsum(w*exp(-r*t)*(t if gap==0 else expm1(gap*t)/gap) for t,w in quadrature(0.,d))
        else:generated_years=(_phi(b,d)-_phi(r,d))/gap
        pF=F*_phi(r,d)+future_pressure*B*generated_years
        nextF=F*exp(-r*d)+future_pressure*B*exp(-r*d)*(d if gap==0 else expm1(gap*d)/gap)
        cases+=h*pB+fresh_h*pF;deaths+=mortality_rate*(pB+pF);years+=pB+pF
        B*=exp(-b*d);F=nextF
    require(min(B,F,cases,deaths,years)>=-1e-13 and abs(B+F+cases+deaths-1)<1e-11,
        'Independent competing-event audit did not conserve mass.')
    return {'cases':cases,'deaths':deaths,'person_years':years,'surviving_disease_free':B+F}


def _history_bins(history,sex,age,time,curve,width):
    raw=[(number(r['start']),number(r['end']),number(r['shape'])) for r in history_curve(history,sex)]
    require(raw[0][0]<=time-age and raw[-1][1]>=time,'Complete historical exposure required.')
    points={0.,age}
    points.update(time-x for a,b,q in raw for x in (a,b) if 0<time-x<age)
    for k in curve['knots']:
        t=number(k['time_years']);points.update(x for x in (t,t-1) if 0<x<age)
    cuts=sorted(points)
    for a,b in zip(cuts,cuts[1:]):
        n=max(1,ceil((b-a)/width));points.update(a+(b-a)*i/n for i in range(1,n))
    bins=[];Q=0.
    for a,b in zip(sorted(points),sorted(points)[1:]):
        q=next(q for lo,hi,q in raw if lo<=time-(a+b)/2<hi)
        bins.append((a,b,q,Q));Q+=q*(b-a)
    return bins,Q


def cohort_reference(history,strata,curve,alpha,*,mortality_rate=0.,future_pressure=0.,width=.5,time=2024.,instantaneous=False):
    """Positive GL8 integration independent of3C's linear response interpolation.

    Exact-age mathematical fixtures, independent factors and c=0 only. More
    general clinical/age/state policies remain the responsibility of3B.
    """
    require(0<width<=1 and alpha>=0,'Bounded quadrature width and nonnegative scale required.')
    numerator=denominator=population=0.
    for row in strata:
        age=number(row['age']);N=number(row['population']);require(age>=15 and number(row['association_control'])==0,'Adult independent-association audit domain only.')
        bins,Q=_history_bins(history,row['sex'],age,time,curve,width);never=exp(-alpha*Q)
        patterns=[(1.,1.)]
        for f in row['factors']:
            p=number(f['prevalence']);rr=number(f['rr'])
            patterns=[(mass*weight,max(m,rr) if present else m) for mass,m in patterns for weight,present in [(1-p,False),(p,True)] if weight]
        for mass,rr in patterns:
            n=d=0.
            if not instantaneous:
                y=one_year_kernel(curve,None,rr=rr,mortality_rate=mortality_rate,future_pressure=future_pressure)
                n=never*y['cases'];d=never*y['person_years']
            for a,b,q,Qa in bins:
                for tau,w in quadrature(a,b):
                    density=alpha*q*exp(-alpha*(Qa+q*(tau-a)))
                    if instantaneous:n+=w*density*rr*curve_hazard(curve,tau)
                    else:
                        y=one_year_kernel(curve,tau,rr=rr,mortality_rate=mortality_rate,future_pressure=future_pressure)
                        n+=w*density*y['cases'];d+=w*density*y['person_years']
            numerator+=N*mass*n;denominator+=N*mass*(1. if instantaneous else d)
        population+=N
    return {'incidence':100000*numerator/denominator,'numerator':numerator,'person_years':denominator,'population':population,
        'scope':'independent_one_year_mathematical_audit_not_a_country_estimate'}


def band_instantaneous_reference(history,curve,alpha,*,age_lower,age_upper_exclusive,time,width=1.):
    """Uniform age integration is an explicit diagnostic assumption, not demography.

    This uses a common sex curve and RR1 without mortality/future infection;
    no female/male population split is invented from a both-sex observation.
    """
    require(15<=age_lower<age_upper_exclusive and 0<width<=1,'Finite adult band required.')
    values=[]
    for age in range(age_lower,age_upper_exclusive):
        for exact,w in quadrature(float(age),float(age+1)):
            r=cohort_reference(history,[{'age':exact,'sex':'both','population':1,'association_control':0,'factors':[]}],
                curve,alpha,width=width,time=time,instantaneous=True)
            values.append(w*r['incidence'])
    return fsum(values)/(age_upper_exclusive-age_lower)
