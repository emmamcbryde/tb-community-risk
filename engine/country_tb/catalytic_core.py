"""Repeated-event Poisson catalytic integrals and exact uniform age averaging."""
from __future__ import annotations

from decimal import Decimal

from engine.country_tb.calibration_common import scientific_context, artifact, reviewed, states_from_survivors, verify
from engine.country_tb.data_common import DataValidationError, decimal, numeric_text
from engine.country_tb.pressure_shape import HISTORY, _normalize_intervals, validate_history
from engine.country_tb.most_recent_event import age_event_summaries, survival, event_probability


def history_curve(history, sex):
    validate_history(history)
    if sex not in {"female", "male", "both", "unknown"}:
        raise DataValidationError("Sex must be explicitly identified.")
    if sex in {"both", "unknown"} and sex not in history["curves"] and ({"female", "male"} & history["curves"].keys()):
        raise DataValidationError("Sex-specific pressure requires explicit aggregate evidence or sex-specific population weights.")
    return _normalize_intervals(history["curves"].get(sex, history["curves"]["common"]))


def _integral(curve, start, end):
    if start > end or start < decimal(curve[0]["start"]) or end > decimal(curve[-1]["end"]):
        raise DataValidationError("Missing lifetime pressure history; no implicit extrapolation.")
    return sum((decimal(row["shape"]) * (min(end, decimal(row["end"])) - max(start, decimal(row["start"])))
                for row in curve if min(end, decimal(row["end"])) > max(start, decimal(row["start"]))), Decimal(0))


def _exposure(curve, age, time, window):
    birth = time - age
    split = max(birth, time - window)
    remote, recent = _integral(curve, birth, split), _integral(curve, split, time)
    return remote, recent


def exposure_components(history, *, sex, age_lower, age_upper_exclusive=None,
                        time="2024", recent_window="5", age_policy="exact_age", review=None):
    """Reusable alpha-free exposure pieces, including analytic uniform age cells.

    Open tails require explicit reviewed discrete ages outside this function.
    No age interval is silently represented by its lower bound or midpoint.
    """
    with scientific_context():
        curve = history_curve(history, sex)
        time, lower, window = decimal(time), decimal(age_lower), decimal(recent_window)
        if window <= 0 or time >= 2025:
            raise DataValidationError("Invalid recent window or post-2024 evaluation.")
        if age_policy == "exact_age":
            if age_upper_exclusive is not None:
                raise DataValidationError("Exact age must not also specify an age interval.")
            remote, recent = _exposure(curve, lower, time, window)
            return [{"weight": "1", "remote_lo": numeric_text(remote), "remote_hi": numeric_text(remote),
                     "recent_lo": numeric_text(recent), "recent_hi": numeric_text(recent),
                     "total_lo": numeric_text(remote + recent), "total_hi": numeric_text(remote + recent)}]
        if age_upper_exclusive is None:
            raise DataValidationError("An open age interval needs reviewed external age weights.")
        upper = decimal(age_upper_exclusive)
        if upper <= lower:
            raise DataValidationError("Invalid age interval.")
        if age_policy == "uniform_within_single_year":
            if upper - lower != 1:
                raise DataValidationError("The single-year policy cannot split a broader source interval.")
        elif age_policy == "reviewed_uniform_within_interval":
            reviewed(review)
        else:
            raise DataValidationError("Explicit supported age-averaging policy required.")
        # Changes of slope occur at birth crossing a pressure boundary, or a=w.
        knots = {lower, upper}
        knots.update(time - decimal(row["start"]) for row in curve)
        knots.add(window)
        points = sorted(value for value in knots if lower <= value <= upper)
        pieces = []
        for lo, hi in zip(points, points[1:]):
            r0, h0 = _exposure(curve, lo, time, window)
            r1, h1 = _exposure(curve, hi, time, window)
            pieces.append({"weight": numeric_text((hi - lo) / (upper - lower)),
                "remote_lo": numeric_text(r0), "remote_hi": numeric_text(r1),
                "recent_lo": numeric_text(h0), "recent_hi": numeric_text(h1),
                "total_lo": numeric_text(r0 + h0), "total_hi": numeric_text(r1 + h1)})
        return pieces


def _average_exponential(alpha, lo, hi):
    lo, hi = decimal(lo), decimal(hi)
    if hi < lo:
        raise DataValidationError("Nonmonotone exposure.")
    difference = alpha * (hi - lo)
    if difference == 0:
        factor = Decimal(1)
    elif difference < Decimal("0.00000001"):
        # (1-exp(-d))/d, avoiding cancellation for extremely small exposures.
        term, factor = Decimal(1), Decimal(1)
        for n in range(1, 100):
            term *= -difference / (n + 1)
            previous = factor
            factor += term
            if factor == previous:
                break
    else:
        factor = (1 - (-difference).exp()) / difference
    return survival(alpha * lo) * factor


def survivor_average(components, alpha, *, kind="total"):
    with scientific_context():
        alpha = decimal(alpha)
        if kind not in {"total", "remote", "recent"} or not components:
            raise DataValidationError("Invalid exposure components.")
        weights = [decimal(p["weight"]) for p in components]
        # Decimal division of arbitrary widths can introduce <1e-75 weight error.
        if abs(sum(weights) - 1) > Decimal("1e-70"):
            raise DataValidationError("Age weights do not sum to one.")
        def endpoints(p):
            if kind == 'recent':
                return p['recent_lo'], p['recent_hi']
            return p[kind + '_lo'], p[kind + '_hi']
        return sum(weight * _average_exponential(alpha, *endpoints(p))
                   for weight, p in zip(weights, components)) / sum(weights)


def component_states(components, alpha):
    with scientific_context():
        never = survivor_average(components, alpha)
        no_recent = survivor_average(components, alpha, kind='recent')
        scale = decimal(alpha)
        prior = Decimal(0)
        for p in components:
            if decimal(p['remote_hi']) == 0:
                continue
            recent = decimal(p['recent_hi'])
            remote_any = 1-_average_exponential(scale, p['remote_lo'], p['remote_hi'])
            prior += decimal(p['weight'])*event_probability(scale*recent)*remote_any
        return states_from_survivors(never, no_recent, prior/sum(decimal(p['weight']) for p in components))


def evaluate_age(history, *, alpha, age, sex, time="2024", recent_window="5",
                 age_upper_exclusive=None, age_policy="exact_age", review=None, include_event_times=True):
    with scientific_context():
        pieces = exposure_components(history, sex=sex, age_lower=age, age_upper_exclusive=age_upper_exclusive,
            time=time, recent_window=recent_window, age_policy=age_policy, review=review)
        scale = decimal(alpha)
        total = sum(decimal(p["weight"]) * (decimal(p["total_lo"]) + decimal(p["total_hi"])) / 2 for p in pieces)
        remote = sum(decimal(p["weight"]) * (decimal(p["remote_lo"]) + decimal(p["remote_hi"])) / 2 for p in pieces)
        summaries = {}
        if include_event_times:
            curve = history_curve(history, sex)
            summaries = age_event_summaries(curve, scale, time, age, age_upper_exclusive, recent_window)
        states = component_states(pieces, scale)
        return artifact("country_tb.catalytic_age_result.v1",
            geography_id=history["geography_id"], sex=sex, time=numeric_text(time),
            age_lower=numeric_text(age), age_upper_exclusive=numeric_text(age_upper_exclusive) if age_upper_exclusive is not None else None,
            age_policy=age_policy, alpha=numeric_text(scale), alpha_unit="per_person_year",
            recent_window_years=numeric_text(recent_window), history_hash=history["profile_hash"],
            cumulative_hazard=numeric_text(scale * total), remote_hazard=numeric_text(scale * remote),
            recent_hazard=numeric_text(scale * (total - remote)),
            expected_lifetime_infection_events=numeric_text(scale * total),
            most_recent_infection_time=summaries.get('any_tbi'),
            recent_time_summary=summaries.get('recent'), remote_only_time_summary=summaries.get('remote_only'),
            expected_time_since_most_recent_infection=summaries.get('any_tbi', {}).get('expected_time_since_most_recent_infection'),
            hazard_interpretation="exact_age" if age_policy == "exact_age" else "mean_hazard_not_hazard_of_mean_prevalence",
            probability_rounding_absolute_bound="0.0000000000000000000000000000000000000001",
            **states)
