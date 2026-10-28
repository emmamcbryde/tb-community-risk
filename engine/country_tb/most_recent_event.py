"""Conditional last-event distributions for piecewise-constant Poisson pressure.

Calendar event times are deterministic; sampling always needs a caller RNG.
Uniform age cells are integrated over their birth interval, never a midpoint.
"""
from decimal import Decimal
from copy import deepcopy
from functools import lru_cache

from engine.country_tb.calibration_common import scientific_context, probability
from engine.country_tb.data_common import DataValidationError, decimal, numeric_text


def survival(hazard):
    # Below Decimal's useful range this is zero to far better than the declared
    # precision. Avoid Decimal underflow for arbitrarily large finite hazards.
    return Decimal(0) if hazard > 100000 else (-hazard).exp()


def event_probability(hazard):
    if hazard >= Decimal('.01'):
        return 1 - survival(hazard)
    term = hazard
    result = term
    for n in range(2, 200):
        term *= -hazard / n
        previous = result
        result += term
        if previous == result:
            break
    return result


def _log1p(value):
    if abs(value) >= Decimal('.01'):
        return (1+value).ln()
    term = value
    result = term
    for n in range(2,200):
        term *= -value
        previous = result
        result += term/n
        if result == previous:
            break
    return result


def _moments(rate, duration):
    """rate * integral x**k exp(-rate*x) dx, k=0,1,2; stable at zero."""
    h = rate * duration
    if h == 0:
        return (Decimal(0),) * 3
    if h < 1:
        out = []
        for k in range(3):
            term = Decimal(1)
            value = term / (k + 1)
            for n in range(1, 250):
                term *= -h / n
                previous = value
                value += term / (n + k + 1)
                if value == previous:
                    break
            out.append(h * duration**k * value)
        return out
    e = survival(h)
    m0 = event_probability(h)
    return m0, (m0 - h * e) / rate, (2 * m0 - (h*h + 2*h)*e) / rate**2


def _mass(segment, start=None, end=None, *, moment=True):
    lo, hi, rate, after, eligible, slope = segment
    start = lo if start is None else start
    end = hi if end is None else end
    a = eligible - slope * (hi - end)
    factor = survival(after + rate * (hi - end))
    if not moment and slope == 0:
        return factor*a*event_probability(rate*(end-start)), Decimal(0)
    m0, m1, m2 = _moments(rate, end - start)
    return factor * (a*m0 - slope*m1), factor * (a*m1 - slope*m2)


def event_distribution(curve, *, alpha, baseline_time, age_lower, age_upper=None,
                       recent_window=5, group='any_tbi'):
    """Internal density segments: lambda exp(-future hazard) * birth eligibility.

    For remote_only there must also be no event in the recent window; the
    future-hazard factor extends to baseline, including that window.
    """
    with scientific_context():
        alpha, time, age, window = map(decimal, (alpha, baseline_time, age_lower, recent_window))
        upper = age if age_upper is None else decimal(age_upper)
        if upper < age or window <= 0 or group not in {'any_tbi', 'recent', 'remote_only'}:
            raise DataValidationError('Invalid last-event age interval, window or group.')
        birth_lo, birth_hi = time-upper, time-age
        if birth_lo < decimal(curve[0]['start']) or time > decimal(curve[-1]['end']):
            raise DataValidationError('Missing lifetime pressure history for last-event distribution.')
        lower = max(birth_lo, time-window) if group == 'recent' else birth_lo
        end = min(time, time-window) if group == 'remote_only' else time
        knots = sorted({lower, end, *[decimal(r['start']) for r in curve if lower < decimal(r['start']) < end],
                        *([birth_hi] if lower < birth_hi < end else [])})
        segments = []
        for lo, hi in zip(knots, knots[1:]):
            if hi <= lo or lo < lower or hi > end:
                continue
            rate = alpha * next(decimal(r['shape']) for r in curve
                                if decimal(r['start']) <= lo < decimal(r['end']))
            after = alpha * sum((decimal(r['shape']) * (min(time, decimal(r['end'])) - max(hi, decimal(r['start'])))
                                for r in curve if min(time, decimal(r['end'])) > max(hi, decimal(r['start']))), Decimal(0))
            slope = Decimal(1)/(upper-age) if upper > age and hi <= birth_hi else Decimal(0)
            eligible = (hi-birth_lo)/(upper-age) if slope else Decimal(1)
            if rate > 0:
                segments.append((lo, hi, rate, after, eligible, slope))
        log_scale = min((s[3] for s in segments), default=Decimal(0))
        segments = [(*s[:3], s[3]-log_scale, *s[4:]) for s in segments]
        return {'segments': segments, 'baseline_time': time, 'group': group, 'log_mass_scale': -log_scale,
                'weights': [_mass(s)[0] for s in segments]}


def _weights(distribution):
    if 'weights' in distribution:
        return distribution['weights']
    return [_mass(s)[0] for s in distribution['segments']]


def mixture_distribution(weighted_distributions):
    """Mix unconditional group densities before conditioning, including rare groups."""
    with scientific_context():
        total = sum(decimal(w) for w, _ in weighted_distributions)
        scale = max(d['log_mass_scale'] for _, d in weighted_distributions)
        segments = []
        for w, d in weighted_distributions:
            factor = decimal(w)/total*survival(scale-d['log_mass_scale'])
            segments.extend([(*s[:4], s[4]*factor, s[5]*factor) for s in d['segments']])
        segments.sort(key=lambda s: (s[0],s[1]))
        # This representation can overlap, so inverse CDF uses global bisection.
        return {'segments': segments, 'weights': [_mass(s)[0] for s in segments],
                'baseline_time': weighted_distributions[0][1]['baseline_time'],
                'group': weighted_distributions[0][1]['group'], 'log_mass_scale': scale, 'mixture': True}


def conditional_cdf(distribution, event_time):
    """P(most_recent_infection_time <= event_time | group); null at zero mass."""
    with scientific_context():
        when = decimal(event_time)
        segments = distribution['segments']
        weights = _weights(distribution)
        total = sum(weights)
        if total == 0:
            return None
        if when <= min(s[0] for s in segments):
            return Decimal(0)
        if when >= max(s[1] for s in segments):
            return Decimal(1)
        mass = Decimal(0)
        for s, weight in zip(segments, weights):
            if when >= s[1]:
                mass += weight
            elif when > s[0]:
                mass += _mass(s, end=when, moment=False)[0]
        return mass / total


def inverse_cdf(distribution, quantile):
    """Stable bounded inversion within the selected density segment."""
    with scientific_context():
        p = probability(quantile)
        weights = _weights(distribution)
        total = sum(weights)
        if total == 0:
            return None
        active = [(s, w) for s, w in zip(distribution['segments'], weights) if w > 0]
        if p == 0:
            return active[0][0][0]
        if p == 1:
            return active[-1][0][1]
        if distribution.get('mixture'):
            lo, hi = min(s[0] for s,_ in active), max(s[1] for s,_ in active)
            for _ in range(180):
                mid = (lo+hi)/2
                if conditional_cdf(distribution, mid) < p:
                    lo = mid
                else:
                    hi = mid
                if hi-lo <= Decimal('1e-24'):
                    break
            return (lo+hi)/2
        target, accumulated = p * total, Decimal(0)
        for s, weight in active:
            if accumulated + weight >= target:
                local = target - accumulated
                if s[5] == 0:
                    fraction = local/weight
                    if fraction == 0:
                        return s[0]
                    if fraction == 1:
                        return s[1]
                    h = s[2]*(s[1]-s[0])
                    # log(exp(-h)+fraction*(1-exp(-h))) without cancellation.
                    log_tail = _log1p(-(1-fraction)*event_probability(h)) if h < 1 else (
                        survival(h)+fraction*event_probability(h)).ln()
                    return s[1]+log_tail/s[2]
                lo, hi = s[0], s[1]
                # Fixed iteration budget with a time tolerance below published
                # precision; no stochastic optimizer or hidden global RNG.
                for _ in range(180):
                    mid = (lo+hi)/2
                    if _mass(s, end=mid, moment=False)[0] < local:
                        lo = mid
                    else:
                        hi = mid
                    if hi-lo <= Decimal('1e-24'):
                        break
                return (lo+hi)/2
            accumulated += weight
        raise DataValidationError('Last-event inverse CDF failed to bracket.')


def summarize(distribution, quantiles=('.05', '.5', '.95')):
    with scientific_context():
        weights = _weights(distribution)
        mass = sum(weights)
        time = distribution['baseline_time']
        mean = None if not mass else sum((time-s[1])*w + _mass(s)[1]
                    for s, w in zip(distribution['segments'], weights))/mass
        event_quantiles = {numeric_text(p): inverse_cdf(distribution, p) for p in quantiles}
        return {
            'conditional_on': distribution['group'], 'status': 'defined' if mass else 'undefined_zero_probability',
            'group_probability': numeric_text(mass*survival(-distribution.get('log_mass_scale', Decimal(0)))),
            'log_group_probability': numeric_text(mass.ln()+distribution.get('log_mass_scale', Decimal(0))) if mass else None,
            'expected_time_since_most_recent_infection': numeric_text(mean) if mean is not None else None,
            'most_recent_infection_time': {
                'expected': numeric_text(time-mean) if mean is not None else None,
                'quantiles': {p: numeric_text(v) if v is not None else None for p, v in event_quantiles.items()}},
            'time_since_most_recent_infection_quantiles': {
                numeric_text(p): numeric_text(time-v) if (v := event_quantiles.get(numeric_text(1-decimal(p)))) is not None else
                    numeric_text(time-v) if (v := inverse_cdf(distribution, 1-decimal(p))) is not None else None for p in quantiles},
            'method': 'analytic_piecewise_exponential_density_uniform_birth_eligibility_bounded_inverse_cdf',
            'inverse_time_tolerance_years': '0.000000000000000000000001',
        }


def sample_most_recent_infection_time(distribution, *, rng):
    if rng is None or not callable(getattr(rng, 'random', None)):
        raise DataValidationError('An explicit RNG is required.')
    return inverse_cdf(distribution, str(rng.random()))


@lru_cache(maxsize=2048)
def _cached_age_summaries(curve_key, alpha, time, age, upper, window):
    curve = [dict(zip(('start','end','shape'), row)) for row in curve_key]
    return {group: summarize(event_distribution(curve, alpha=alpha, baseline_time=time,
        age_lower=age, age_upper=upper, recent_window=window, group=group))
        for group in ('any_tbi', 'recent', 'remote_only')}


def age_event_summaries(curve, alpha, time, age, upper, window):
    key = tuple(tuple(row[k] for k in ('start','end','shape')) for row in curve)
    return deepcopy(_cached_age_summaries(key, numeric_text(alpha), numeric_text(time), numeric_text(age),
        numeric_text(upper) if upper is not None else None, numeric_text(window)))


def interval_distribution(history, *, sex, alpha, start, end, baseline_time=None):
    """Conditional last-event law on [start,end], with optional later clock date."""
    from engine.country_tb.catalytic_core import history_curve
    with scientific_context():
        c, b = decimal(start), decimal(end)
        if c > b:
            raise DataValidationError('Inverted last-event interval.')
        result = event_distribution(history_curve(history, sex), alpha=alpha,
            baseline_time=b, age_lower=b-c)
        result['baseline_time'] = b if baseline_time is None else decimal(baseline_time)
        if result['baseline_time'] < b:
            raise DataValidationError('Clock baseline precedes the infection interval.')
        return result
