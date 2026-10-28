"""Isolated synthetic first-infection mathematics, NOT a country calibration.

No reinfection, clearance, disease exits, selective mortality or migration model.
Conditional-survivor interpretation therefore requires the assumptions in the
Milestone 2B specification. No WHO-to-infection mapping is executed here.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
import math

from engine.country_tb.data_common import DataValidationError, decimal, numeric_text, seal

VERSION = "country_tb_catalytic_diagnostic_1.0.0"
QUANTUM = Decimal("1e-40")


def normalized_history(history):
    result = []
    for row in history:
        start = decimal(row["start"], nonnegative=False)
        end = decimal(row["end"], nonnegative=False)
        hazard = decimal(row["hazard_per_person_year"])
        if end <= start:
            raise DataValidationError("History intervals must have positive duration.")
        result.append((start, end, hazard))
    result.sort()
    if not result:
        raise DataValidationError("An explicit infection-pressure history is required.")
    if any(left[1] != right[0] for left, right in zip(result, result[1:])):
        raise DataValidationError("History has an overlap or gap; no implicit interpolation/backcast.")
    return result


def _integral(rows, start, end):
    if start < rows[0][0] or end > rows[-1][1]:
        raise DataValidationError("Missing lifetime history; choose a reviewed backcast explicitly.")
    return sum((rate * (min(end, hi) - max(start, lo)) for lo, hi, rate in rows
                if min(end, hi) > max(start, lo)), Decimal(0))


def first_infection_profile(history, *, age, baseline_time, recent_window, sex,
                            geography_id, scenario_id, provenance):
    """Piecewise-constant calendar hazards; exact ages, closed history coverage.

    Recent = FIRST infection within the recent window. These states are not a
    most-recent-event/reinfection calculation or a validated viable-TBI profile.
    Probabilities round survivor functions to 40 decimal places; residual states
    are algebraic differences, with no clipping or post-hoc normalization.
    """
    if sex not in {"female", "male", "both", "unknown"} or not geography_id or not scenario_id or not provenance:
        raise DataValidationError("Explicit sex, geography, scenario and provenance are required.")
    with localcontext() as context:
        context.prec = 80
        rows = normalized_history(history)
        age = decimal(age)
        baseline = decimal(baseline_time, nonnegative=False)
        window = decimal(recent_window)
        if window <= 0:
            raise DataValidationError("Recent window must be positive.")
        birth = baseline - age
        split = max(birth, baseline - window)
        early = _integral(rows, birth, split)
        recent_hazard = _integral(rows, split, baseline)
        total = early + recent_hazard
        try:
            uninfected = (-total).exp().quantize(QUANTUM)
            before = (-early).exp().quantize(QUANTUM)
        except ArithmeticError as exc:
            raise DataValidationError("Numerical range exceeded; no clipped output.") from exc
        recent, remote = before - uninfected, Decimal(1) - before
        if min(uninfected, recent, remote) < 0 or max(uninfected, recent, remote) > 1:
            raise DataValidationError("Invalid probability distribution.")
        assert uninfected + recent + remote == 1
        return seal({
            "schema_id": "country_tb.catalytic_diagnostic.v1", "contract_version": "1.0.0",
            "implementation_version": VERSION, "status": "synthetic_diagnostic_only",
            "age": numeric_text(age), "baseline_time": numeric_text(baseline),
            "recent_window_years": numeric_text(window), "sex": sex, "geography_id": geography_id,
            "scenario_id": scenario_id, "provenance": provenance,
            "history": [{"start": numeric_text(lo), "end": numeric_text(hi),
                         "hazard_per_person_year": numeric_text(rate)} for lo, hi, rate in rows],
            "cumulative_hazard": numeric_text(total), "recent_interval_hazard": numeric_text(recent_hazard),
            "uninfected": numeric_text(uninfected), "recently_first_infected": numeric_text(recent),
            "remote_only_infected": numeric_text(remote), "total_ever_infected": numeric_text(1 - uninfected),
            "probability_rounding_absolute_bound": "0.0000000000000000000000000000000000000001",
            "exclusions": ["reinfection", "clearance", "disease_exits", "selective_mortality", "migration"],
        })


def synthetic_scale_solution(exposures, weights, target, *, upper_scale):
    """Bisection proof for a PERFECT-ASSAY synthetic anchor, not survey fitting.

    Exposures are known lifetime integrals of a dimensionless shape, in years.
    Only one nonnegative scale is solved. No incidence or empirical data access.
    Strict feasibility uses the supplied finite upper bound; 1 is not attainable
    with finite hazard. Returns the declared bracket accuracy.
    """
    if not exposures or len(exposures) != len(weights):
        raise DataValidationError("Nonempty matched exposures and weights required.")
    with localcontext() as context:
        context.prec = 80
        h = [decimal(v) for v in exposures]
        w = [decimal(v) for v in weights]
        p = decimal(target)
        upper = decimal(upper_scale)
        if p >= 1 or upper <= 0 or sum(w) == 0:
            raise DataValidationError("Invalid target, bound or weights; no clipping.")
        denominator = sum(w)
        def predicted(k):
            return sum(weight * (1 - (-k * exposure).exp()) for weight, exposure in zip(w, h)) / denominator
        if not any(exposure > 0 and weight > 0 for exposure, weight in zip(h, w)):
            raise DataValidationError("Scale unidentifiable under zero weighted exposure.")
        if p > predicted(upper):
            raise DataValidationError("Target incompatible with permitted scale/history.")
        low, high = Decimal(0), upper
        for _ in range(180):
            midpoint = (low + high) / 2
            if predicted(midpoint) < p:
                low = midpoint
            else:
                high = midpoint
        scale = Decimal(0) if p == 0 else (low + high) / 2
        return {"scale": numeric_text(scale), "predicted_prevalence": numeric_text(predicted(scale)),
                "target": numeric_text(p), "scale_bracket_width": numeric_text(high - low),
                "method": "synthetic_perfect_assay_bisection_180", "empirical_fit": False}


def sample_states(profile, n, *, rng):
    """Categorical diagnostic only; caller supplies an RNG, never global state."""
    if type(n) is not int or n <= 0:
        raise DataValidationError("Positive integer sample size required.")
    if profile.get("status") != "synthetic_diagnostic_only":
        raise DataValidationError("Only an explicit synthetic diagnostic is accepted.")
    probabilities = [decimal(profile[k]) for k in ("uninfected", "recently_first_infected", "remote_only_infected")]
    with localcontext() as context:
        context.prec = 80
        if any(p > 1 for p in probabilities) or sum(probabilities) != 1:
            raise DataValidationError("Invalid probabilities; no normalization.")
    thresholds = [float(probabilities[0]), float(probabilities[0] + probabilities[1])]
    counts = [0, 0, 0]
    for _ in range(n):
        value = rng.random()
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or not 0 <= value < 1:
            raise DataValidationError("RNG must return finite uniform values in [0,1).")
        counts[0 if value < thresholds[0] else 1 if value < thresholds[1] else 2] += 1
    return dict(zip(("uninfected", "recently_first_infected", "remote_only_infected"), counts))
