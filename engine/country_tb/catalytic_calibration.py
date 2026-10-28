"""Deterministic bounded one-scale calibration with fixed shapes and assays."""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal

from engine.country_tb.calibration_common import scientific_context, artifact, finite_integer, reviewed, verify
from engine.country_tb.catalytic_core import exposure_components, survivor_average
from engine.country_tb.data_common import DataValidationError, canonical, decimal, digest, numeric_text
from engine.country_tb.pressure_shape import HISTORY, validate_history
from engine.country_tb.tbi_observations import OBSERVATIONS, assay_probability

POLICY = "country_tb.tbi_calibration_policy.v1"
CALIBRATION = "country_tb.catalytic_scale_calibration.v1"


def calibration_policy(*, alpha_upper, observation_weights, residual_tolerance,
                       review, alpha_tolerance="0.000000000000000001",
                       prevalence_tolerance="0.000000000001", max_iterations=200,
                       uncertainty="none", method="overall_anchor"):
    reviewed(review)
    upper, a_tol, p_tol = decimal(alpha_upper), decimal(alpha_tolerance), decimal(prevalence_tolerance)
    residual = decimal(residual_tolerance)
    if upper <= 0 or a_tol <= 0 or a_tol >= upper or not 0 < p_tol <= 1 or residual > 1:
        raise DataValidationError("Finite positive scale bound and numerical tolerances required.")
    finite_integer(max_iterations, minimum=1)
    if max_iterations > 1000 or uncertainty not in {"none", "coherent_observation_intervals"}:
        raise DataValidationError("Unsupported solver/uncertainty policy.")
    if not isinstance(observation_weights, dict):
        raise DataValidationError("Explicit observation weight mapping required.")
    weights = {key: numeric_text(decimal(value)) for key, value in observation_weights.items()}
    if any(not key or decimal(value) <= 0 for key, value in weights.items()):
        raise DataValidationError("Explicit observation weights must be positive.")
    if method not in {"overall_anchor", "binomial", "weighted_least_squares"}:
        raise DataValidationError("Select overall_anchor, binomial or weighted_least_squares explicitly.")
    if uncertainty != 'none' and method != 'overall_anchor':
        raise DataValidationError('Coherent interval inversion requires one overall anchor; provide separate reviewed scenarios for multiple-anchor objectives.')
    if method == 'binomial' and any(decimal(v) != 1 for v in weights.values()):
        raise DataValidationError('Independent binomial likelihood uses actual sample counts, with unit observation weights.')
    return artifact(POLICY, method=method, objective=method,
        free_parameters=["alpha"], alpha_lower="0", alpha_upper=numeric_text(upper),
        alpha_unit="per_person_year", alpha_tolerance=numeric_text(a_tol),
        prevalence_tolerance=numeric_text(p_tol), residual_tolerance=numeric_text(residual),
        max_iterations=max_iterations, observation_weights=weights, uncertainty=uncertainty,
        review=deepcopy(review), solver="bounded_bisection_or_golden_section", likelihood_claim=method == "binomial",
        incompatibility_action="diagnostic_only_no_accepted_profile")


def prepare_observations(history, observations):
    verify(observations, OBSERVATIONS)
    reviewed(observations["review"])
    result = []
    for observation in observations["rows"]:
        receipt = observation["receipt"]
        if digest(canonical(receipt["source_record"])) != receipt["source_record_sha256"]:
            raise DataValidationError("Original observation receipt checksum mismatch.")
        if observation["geography_id"] != history["geography_id"]:
            raise DataValidationError("No automatic geographic transport of a TBI anchor.")
        entries = []
        for row in observation["strata"]:
            components = exposure_components(history, sex=row["sex"], age_lower=row["age_lower"],
                age_upper_exclusive=row["age_upper_exclusive"], time=observation["survey_time"],
                recent_window="5", age_policy=row["age_policy"], review=row.get("age_review"))
            entries.append({"weight": decimal(row["weight"]), "components": components})
        result.append({"observation": observation, "entries": entries})
    return result


def predicted_positivity(prepared, alpha):
    with scientific_context():
        denominator = sum(entry["weight"] for entry in prepared["entries"])
        tbi = sum(entry["weight"] * (1 - survivor_average(entry["components"], alpha))
                  for entry in prepared["entries"]) / denominator
        return assay_probability(tbi, prepared["observation"]["assay_operator"])


def _asymptote(prepared):
    denominator = sum(entry["weight"] for entry in prepared["entries"])
    never = Decimal(0)
    for entry in prepared["entries"]:
        pieces = entry["components"]
        total_weight = sum(decimal(p["weight"]) for p in pieces)
        zero_weight = sum((decimal(p["weight"]) for p in pieces
                           if decimal(p["total_lo"]) == decimal(p["total_hi"]) == 0), Decimal(0))
        never += entry["weight"] * zero_weight / total_weight
    return assay_probability(1 - never / denominator, prepared["observation"]["assay_operator"])


def calibrate_scale(history, policy, *, observations=None, supplied_alpha=None,
                    scale_provenance=None, supplied_alpha_interval=None):
    validate_history(history)
    verify(policy, POLICY)
    reviewed(policy["review"])
    if policy["method"] not in {"overall_anchor", "binomial", "weighted_least_squares"} or policy["free_parameters"] != ["alpha"]:
        raise DataValidationError("Only a reviewed one-scale policy is implemented.")
    with scientific_context():
        upper = decimal(policy["alpha_upper"])
        if observations is None:
            if supplied_alpha is None or not scale_provenance or policy["observation_weights"]:
                raise DataValidationError("Without an anchor, supply alpha and source/scenario provenance; no inferred default.")
            alpha = decimal(supplied_alpha)
            if alpha > upper:
                raise DataValidationError("Supplied scale exceeds the explicit policy bound.")
            interval = {"kind": "not_supplied"}
            if supplied_alpha_interval is not None:
                lo, hi = map(decimal, (supplied_alpha_interval["lower"], supplied_alpha_interval["upper"]))
                if not 0 <= lo <= alpha <= hi <= upper or not supplied_alpha_interval.get("provenance"):
                    raise DataValidationError("Invalid supplied scale uncertainty/provenance.")
                interval = {"kind": "supplied_scenario_interval", "lower": numeric_text(lo), "upper": numeric_text(hi),
                            "provenance": deepcopy(supplied_alpha_interval["provenance"])}
            return artifact(CALIBRATION, fit_status="supplied_scale", alpha=numeric_text(alpha),
                alpha_unit="per_person_year", scale_provenance=deepcopy(scale_provenance),
                pressure_history=deepcopy(history), pressure_history_hash=history["profile_hash"],
                policy=deepcopy(policy), policy_hash=policy["profile_hash"],
                observations=None, observations_hash=None, alpha_uncertainty=interval,
                diagnostics={"identified_from_tbi_data": False, "iterations": 0,
                             "warning": "model_estimated_conditional_on_supplied_scale", "converged": True,
                             "boundary_hit": alpha == 0 or alpha == upper, "objective": None,
                             "feasible_alpha_range": ["0", numeric_text(upper)]}, rows=[])
        if supplied_alpha is not None or supplied_alpha_interval is not None:
            raise DataValidationError("Choose either supplied scale or TBI calibration explicitly.")
        prepared = prepare_observations(history, observations)
        weights = policy["observation_weights"]
        if set(weights) != {p["observation"]["observation_id"] for p in prepared}:
            raise DataValidationError("Policy weights must cover exactly the selected observations.")
        fit_weights = [decimal(weights[p["observation"]["observation_id"]]) for p in prepared]
        denominator = sum(fit_weights)
        target = sum(weight * decimal(p["observation"]["prevalence"]) for weight, p in zip(fit_weights, prepared)) / denominator
        def predict(alpha):
            return sum(weight * predicted_positivity(p, alpha) for weight, p in zip(fit_weights, prepared)) / denominator
        floor = predict(Decimal(0))
        ceiling = sum(weight * _asymptote(p) for weight, p in zip(fit_weights, prepared)) / denominator
        finite_max = predict(upper)
        if ceiling <= floor:
            raise DataValidationError("Scale unidentifiable: all reviewed observation exposure is zero.")
        def solve(value):
            if value < floor or value >= ceiling or value > finite_max:
                raise DataValidationError("Target incompatible with assay/history/finite scale bounds; no clipping.")
            if value == floor:
                return Decimal(0), 0, Decimal(0)
            if value == finite_max:
                return upper, 0, Decimal(0)
            lo, hi = Decimal(0), upper
            for iteration in range(1, policy["max_iterations"] + 1):
                mid = (lo + hi) / 2
                prediction = predict(mid)
                if prediction < value:
                    lo = mid
                else:
                    hi = mid
                if hi - lo <= decimal(policy["alpha_tolerance"]):
                    estimate = (lo + hi) / 2
                    if abs(predict(estimate) - value) <= decimal(policy["prevalence_tolerance"]):
                        return estimate, iteration, hi - lo
            raise DataValidationError("One-scale solver failed the declared convergence tolerances.")
        method = policy['method']
        if method == 'overall_anchor' and len(prepared) != 1:
            raise DataValidationError('Overall root finding requires exactly one reviewed anchor.')
        for p in prepared:
            o = p['observation']
            if method == 'binomial' and (o['numerator'] is None or o['survey_design'].get('kind') not in
                    {'simple_random_sample', 'independent_binomial', 'synthetic_only'}):
                raise DataValidationError('Binomial likelihood requires actual counts and a compatible independent sampling design.')
            if method == 'weighted_least_squares' and o['uncertainty']['kind'] == 'unspecified':
                raise DataValidationError('Weighted least squares requires supplied uncertainty and explicit weights.')
            low, high = predicted_positivity(p, Decimal(0)), predicted_positivity(p, upper)
            asymptote = _asymptote(p)
            if not low <= decimal(o['prevalence']) <= high or (asymptote > low and decimal(o['prevalence']) >= asymptote):
                raise DataValidationError('Target incompatible with assay/history/finite scale bounds; no clipping.')
        def contribution(p, a):
            o = p['observation']
            prediction = predicted_positivity(p, a)
            if method != 'binomial':
                return decimal(weights[o['observation_id']])*(prediction-decimal(o['prevalence']))**2
            x, n = o['numerator'], o['denominator']
            if (prediction == 0 and x) or (prediction == 1 and n-x):
                return Decimal('Infinity')
            return -(Decimal(x)*prediction.ln() if x else Decimal(0)) - (Decimal(n-x)*(1-prediction).ln() if n-x else Decimal(0))
        def objective(a):
            return sum(contribution(p, a) for p in prepared)
        if method == 'overall_anchor':
            alpha, iterations, bracket_width = solve(target)
        else:
            # The one-dimensional objective is evaluated deterministically on a
            # bounded interval. In log-scale coordinates prevalence-only WLS can
            # be nonconvex; scan logarithmic AND linear knots, refine each local
            # basin and compare boundaries. This is not an identifiability claim.
            grid = sorted({Decimal(0), upper, *[upper*Decimal(i)/32 for i in range(1,32)],
                           *[upper*(Decimal(10)**(-i)) for i in range(1,25)]})
            scores = [objective(a) for a in grid]
            candidates = [(scores[0], grid[0]), (scores[-1], grid[-1])]
            ratio = (Decimal(5).sqrt()-1)/2
            iterations = 0
            bracket_width = upper
            for i in range(1,len(grid)-1):
                if scores[i] <= scores[i-1] and scores[i] <= scores[i+1] and (scores[i] < scores[i-1] or scores[i] < scores[i+1]):
                    lo, hi = grid[i-1], grid[i+1]
                    c, d = hi-ratio*(hi-lo), lo+ratio*(hi-lo)
                    fc, fd = objective(c), objective(d)
                    converged = False
                    for iteration in range(1, policy['max_iterations']+1):
                        if hi-lo <= decimal(policy['alpha_tolerance']):
                            converged = True
                            break
                        if fc < fd:
                            hi,d,fd = d,c,fc
                            c = hi-ratio*(hi-lo); fc = objective(c)
                        else:
                            lo,c,fc = c,d,fd
                            d = lo+ratio*(hi-lo); fd = objective(d)
                    if not converged:
                        raise DataValidationError('One-scale objective solver failed convergence tolerances.')
                    iterations += iteration
                    bracket_width = hi-lo
                    a = (lo+hi)/2
                    candidates.append((objective(a),a))
            _,alpha = min(candidates)

        diagnostics = []
        for p in prepared:
            o = p["observation"]
            predicted = predicted_positivity(p, alpha)
            observed = decimal(o["prevalence"])
            uncertainty = o["uncertainty"]
            tol = decimal(policy["prevalence_tolerance"])
            if uncertainty["kind"] == "interval":
                compatible = decimal(uncertainty["lower"]) - tol <= predicted <= decimal(uncertainty["upper"]) + tol
            else:
                allowed = tol if uncertainty["kind"] == "exact" else decimal(policy["residual_tolerance"])
                compatible = abs(predicted - observed) <= allowed
            diagnostics.append({"observation_id": o["observation_id"], "observed": numeric_text(observed),
                "predicted_assay_positivity": numeric_text(predicted), "residual": numeric_text(predicted - observed),
                "compatible": compatible, "objective_contribution": numeric_text(contribution(p, alpha)),
                "numerator": o["numerator"], "denominator": o["denominator"], "uncertainty": deepcopy(uncertainty), "sex": o["sex"],
                'feasible_range': [numeric_text(predicted_positivity(p, Decimal(0))), numeric_text(predicted_positivity(p, upper))],
                "survey_time": o["survey_time"], "assay_operator_hash": o["assay_operator"]["profile_hash"]})
        delta = max(Decimal("0.00000001"), alpha * Decimal("0.0001"))
        left, right = max(Decimal(0), alpha - delta), min(upper, alpha + delta)
        derivative = (predict(right) - predict(left)) / (right - left)
        uncertainty_result = {"kind": "not_requested"}
        if policy["uncertainty"] == "coherent_observation_intervals":
            if not all(p["observation"]["uncertainty"]["kind"] == "interval" for p in prepared):
                raise DataValidationError("Coherent interval sensitivity requires bounds for every observation.")
            uncertainty_result = {"kind": "coherent_observation_interval_sensitivity",
                "interpretation": "endpoint_scenarios_not_a_confidence_or_posterior_interval"}
            for endpoint in ("lower", "upper"):
                value = sum(weight * decimal(p["observation"]["uncertainty"][endpoint])
                            for weight, p in zip(fit_weights, prepared)) / denominator
                try:
                    estimate, _, _ = solve(value)
                    uncertainty_result[endpoint] = numeric_text(estimate)
                    uncertainty_result[endpoint + "_status"] = "solved"
                except DataValidationError:
                    uncertainty_result[endpoint] = None
                    uncertainty_result[endpoint + "_status"] = "outside_declared_feasible_range_not_clipped"
        return artifact(CALIBRATION,
            fit_status="accepted" if all(d["compatible"] for d in diagnostics) and derivative > 0 else "incompatible_anchors",
            alpha=numeric_text(alpha), alpha_unit="per_person_year",
            pressure_history=deepcopy(history), pressure_history_hash=history["profile_hash"],
            policy=deepcopy(policy), policy_hash=policy["profile_hash"],
            observations=deepcopy(observations), observations_hash=observations["profile_hash"],
            alpha_uncertainty=uncertainty_result, rows=diagnostics,
            diagnostics={"identified_from_tbi_data": derivative > 0, "free_parameter_count": 1,
                "conditional_rank": 1 if derivative > 0 else 0, "moment_derivative_per_alpha": numeric_text(derivative),
                "iterations": iterations, "alpha_bracket_width": numeric_text(bracket_width),
                "weighted_observed": numeric_text(target), "weighted_predicted": numeric_text(predict(alpha)),
                "moment_residual": numeric_text(predict(alpha) - target),
                "feasible_positivity_at_zero": numeric_text(floor),
                "feasible_positivity_at_upper_bound": numeric_text(finite_max),
                "asymptotic_positivity_unattainable_at_finite_alpha": numeric_text(ceiling),
                "boundary_hit": alpha == 0 or alpha == upper, "boundary_status": "lower" if alpha == 0 else "upper" if alpha == upper else "interior",
                "converged": True, "objective": numeric_text(objective(alpha)),
                "log_likelihood_without_binomial_constant": numeric_text(-objective(alpha)) if method == "binomial" else None,
                "likelihood_claim": method == "binomial",
                "identifiability": "one_scale_only_conditional_on_pressure_shape_backcast_assay_no_clearance_reinfection_process_demography_and_migration_assumptions; total_TBI_does_not_independently_identify_protection_clearance_progression_clock_or_event_time_detail",
                "warnings": ["conditional_diagnostic_not_country_estimate", "event_time_detail_is_shape_conditional"]})


def restore_scale_calibration(value):
    """Deterministically re-solve against embedded original inputs, offline."""
    verify(value, CALIBRATION)
    if value["pressure_history_hash"] != value["pressure_history"]["profile_hash"] or value["policy_hash"] != value["policy"]["profile_hash"]:
        raise DataValidationError("Calibration source/policy reference mismatch.")
    if value["observations"] is None:
        uncertainty = value["alpha_uncertainty"]
        interval = {key: uncertainty[key] for key in ("lower", "upper", "provenance")} if uncertainty["kind"] == "supplied_scenario_interval" else None
        rebuilt = calibrate_scale(value["pressure_history"], value["policy"], supplied_alpha=value["alpha"],
            scale_provenance=value["scale_provenance"], supplied_alpha_interval=interval)
    else:
        if value["observations_hash"] != value["observations"]["profile_hash"]:
            raise DataValidationError("Calibration observation reference mismatch.")
        rebuilt = calibrate_scale(value["pressure_history"], value["policy"], observations=value["observations"])
    if rebuilt["profile_hash"] != value["profile_hash"]:
        raise DataValidationError("Calibration does not reproduce from its original inputs.")
    return value
