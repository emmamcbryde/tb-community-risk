"""One shared expected age/sex TBI profile; no individual population or runner."""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal
import json

from engine.country_tb.calibration_common import (
    EXCLUSIONS, STATE_FIELDS, artifact, probability, reviewed, states_from_survivors,
    validate_states, verify, scientific_context,
)
from engine.country_tb.catalytic_calibration import CALIBRATION, restore_scale_calibration
from engine.country_tb.catalytic_core import evaluate_age, history_curve
from engine.country_tb.most_recent_event import event_distribution, mixture_distribution, summarize
from engine.country_tb.data_common import DataValidationError, canonical, csv_bytes, decimal, numeric_text, digest
from engine.country_tb.demography import DISPLAY_BANDS, canonical_sexes, normalize_records

PROFILE = "country_tb.calibrated_tbi_age_profile.v1"


def _cell(history, row, alpha, time, window, tail, age_policy, age_review):
    lo, hi, sex = row["age_lower"], row["age_upper_exclusive"], row["sex"]
    if hi is not None:
        return evaluate_age(history, alpha=alpha, age=lo, age_upper_exclusive=hi,
            sex=sex, time=time, recent_window=window, age_policy=age_policy, review=age_review)
    if tail is None:
        return None
    reviewed(tail.get("review"))
    if not tail.get("provenance") or not tail.get("ages"):
        raise DataValidationError("Open-tail ages need explicit source provenance and weights.")
    seen, values = set(), []
    for item in tail["ages"]:
        age, weight = decimal(item["age"]), decimal(item["weight"])
        if age < lo or age in seen or weight <= 0:
            raise DataValidationError("Invalid/duplicate open-tail age or weight.")
        seen.add(age)
        result = evaluate_age(history, alpha=alpha, age=age, sex=sex, time=time, recent_window=window, include_event_times=False)
        values.append((weight, result))
    weight_sum = sum(weight for weight, _ in values)
    never = sum(weight * decimal(result["never_infected"]) for weight, result in values) / weight_sum
    before = sum(weight * (decimal(result["never_infected"]) + decimal(result["remote_only"]))
                 for weight, result in values) / weight_sum
    hazards = {key: numeric_text(sum(weight * decimal(result[key]) for weight, result in values) / weight_sum)
               for key in ("cumulative_hazard", "remote_hazard", "recent_hazard", "expected_lifetime_infection_events")}
    prior = sum(weight*decimal(result['prior_remote_plus_recent']) for weight,result in values)/weight_sum
    curve = history_curve(history, sex)
    summaries = {group: summarize(mixture_distribution([(weight, event_distribution(curve, alpha=alpha,
        baseline_time=time, age_lower=item['age'], recent_window=window, group=group))
        for item,(weight,_) in zip(tail['ages'],values)])) for group in ('any_tbi','recent','remote_only')}
    return {**states_from_survivors(never, before, prior), **hazards,
            'most_recent_infection_time': summaries['any_tbi'], 'recent_time_summary': summaries['recent'],
            'remote_only_time_summary': summaries['remote_only'],
            'expected_time_since_most_recent_infection': summaries['any_tbi']['expected_time_since_most_recent_infection'],
            "hazard_interpretation": "mean_over_reviewed_open_tail_ages"}


def _summary(rows):
    total = sum(decimal(r["population"]) for r in rows)
    covered = sum(decimal(r["population"]) for r in rows if r["any_tbi"] is not None)
    counts = {key + "_count": sum((decimal(r["population"]) * decimal(r[key])
                                  for r in rows if r[key] is not None), Decimal(0))
              for key in (*STATE_FIELDS, "any_tbi", "prior_remote_plus_recent")}
    missing = total - covered
    if sum(counts[key + "_count"] for key in STATE_FIELDS) != covered:
        raise DataValidationError("Expected infection-state counts do not reconcile.")
    return {"population": numeric_text(total), "covered_population": numeric_text(covered),
            "unestimated_population": numeric_text(missing),
            "unestimated_population_proportion": numeric_text(missing/total) if total else None,
            **{key: numeric_text(value) for key, value in counts.items()},
            "any_tbi_prevalence_covered": numeric_text(counts["any_tbi_count"] / covered) if covered else None,
            "any_tbi_prevalence_full": numeric_text(counts["any_tbi_count"] / total) if total and missing == 0 else None,
            "missing_age_coverage_bounds_not_confidence_interval":
                {"lower": numeric_text(counts["any_tbi_count"] / total),
                 "upper": numeric_text((counts["any_tbi_count"] + missing) / total)} if total else None}


def baseline_tbi_profile(demography, calibration, *, recent_window="5", baseline_time="2024",
                         tail_distributions=None, age_policy="uniform_within_single_year", age_review=None, scenario_id="diagnostic"):
    verify(demography, "country_tb.demography.v1")
    verify(calibration, CALIBRATION)
    restore_scale_calibration(calibration)
    if calibration["fit_status"] not in {"accepted", "supplied_scale"}:
        raise DataValidationError("Incompatible anchors cannot produce an accepted baseline profile.")
    history = calibration["pressure_history"]
    if demography["geography_id"] != history["geography_id"] or demography["baseline_year"] != 2024:
        raise DataValidationError("Baseline demography must match the calibrated geography and 2024.")
    with scientific_context():
        time, window = decimal(baseline_time), decimal(recent_window)
        if not 2024 <= time < 2025 or window <= 0:
            raise DataValidationError("Baseline time must be in 2024; recent window must be positive.")
        normalized, duplicate_count = normalize_records(demography["rows"])
        sex_check = canonical_sexes(normalized)
        calculation_keys = {(r["sex"], r["age_lower"], r["age_upper_exclusive"]) for r in normalized}
        if duplicate_count or len(sex_check["rows"]) != len(normalized) or len(calculation_keys) != len(normalized):
            raise DataValidationError("Calculable demography contains duplicates or redundant both-sex rows.")
        tails = deepcopy(tail_distributions or {})
        for tail in tails.values():
            tail["ages"] = sorted([{**r, "age": numeric_text(r["age"]), "weight": numeric_text(r["weight"])}
                                   for r in tail["ages"]], key=lambda r: decimal(r["age"]))
        tail_keys = {r["sex"] + ":" + str(r["age_lower"]) for r in normalized if r["age_upper_exclusive"] is None}
        if set(tails) - tail_keys:
            raise DataValidationError("Unused/mismatched open-tail distributions.")
        alpha = decimal(calibration["alpha"])
        output = []
        uncertainty = calibration["alpha_uncertainty"]
        for row in normalized:
            tail = tails.get(row["sex"] + ":" + str(row["age_lower"]))
            result = _cell(history, row, alpha, time, window, tail, age_policy, age_review)
            population = decimal(row["population"])
            out = {"sex": row["sex"], "age_lower": row["age_lower"], "age_upper_exclusive": row["age_upper_exclusive"],
                "population": numeric_text(population), "source_year": row["year"], "baseline_year": 2024,
                "source_estimate_status": row["estimate_status"], "source_projection_scenario": row["projection_scenario"],
                "geography_id": demography["geography_id"],
                "status": "estimated" if result is not None else "unestimated_open_age_tail",
                "age_policy": age_policy if row["age_upper_exclusive"] is not None else
                              "reviewed_discrete_tail" if tail else "unestimated_no_age_split",
                **{key: result[key] if result is not None else None
                   for key in (*STATE_FIELDS, "any_tbi", "prior_remote_plus_recent", "cumulative_hazard", "remote_hazard", "recent_hazard", "expected_lifetime_infection_events", "expected_time_since_most_recent_infection", "most_recent_infection_time", "recent_time_summary", "remote_only_time_summary")},
                "hazard_interpretation": result["hazard_interpretation"] if result is not None else None,
                "parameter_endpoint_scenarios": {}}
            out.update({
                'pressure_shape_identifier': history['shape_hash'], 'alpha': calibration['alpha'],
                'backcast_policy': history['backcast']['method'], 'recent_window_years': numeric_text(window),
                'calibration_source': calibration['fit_status'], 'scenario_id': scenario_id,
                'reinfection_policy': calibration['reinfection_policy'],
                'progression_clock_policy': calibration['progression_clock_policy'],
                'explicit_exclusions': list(EXCLUSIONS), 'provenance_hash': demography['profile_hash'],
                'scientific_input_hashes': {'history': history['profile_hash'], 'calibration': calibration['profile_hash']},
                'expected_state_counts': {k: numeric_text(population*decimal(out[k])) if result else None for k in (*STATE_FIELDS,'any_tbi','prior_remote_plus_recent')},
                **{alias: out[field] for alias,field in {
                    'neverInfectedProbability': 'never_infected', 'remoteOnlyProbability': 'remote_only',
                    'recentProbability': 'recent', 'anyTbiProbability': 'any_tbi',
                    'priorRemotePlusRecentProbability': 'prior_remote_plus_recent'}.items()}})
            if result is not None:
                validate_states(out)
                for endpoint in ("lower", "upper"):
                    if uncertainty.get(endpoint) is not None:
                        endpoint_result = _cell(history, row, uncertainty[endpoint], time, window, tail, age_policy, age_review)
                        out["parameter_endpoint_scenarios"][endpoint] = {
                            "alpha": uncertainty[endpoint], **{key: endpoint_result[key] for key in (*STATE_FIELDS, "any_tbi", "prior_remote_plus_recent")}}
            out['scientific_hash'] = digest(canonical(out))
            output.append(out)
        output.sort(key=lambda r: (r["sex"], r["age_lower"]))
        total = _summary(output)
        if decimal(total["population"]) <= 0 or total["population"] != numeric_text(demography["combined_population"]):
            raise DataValidationError("Baseline population does not reconcile with source demography.")
        display = []
        for sex in sorted({r["sex"] for r in output}):
            for lo, hi in DISPLAY_BANDS:
                selected = []
                for row in output:
                    upper = row["age_upper_exclusive"]
                    if row["sex"] != sex:
                        continue
                    fits = row["age_lower"] >= lo and (hi is None or (upper is not None and upper <= hi))
                    overlaps = (hi is None or row["age_lower"] < hi) and (upper is None or upper > lo)
                    if overlaps and not fits:
                        raise DataValidationError("A source age cell crosses a display boundary; no population split inferred.")
                    if fits:
                        selected.append(row)
                display.append({"sex": sex, "age_lower": lo, "age_upper_exclusive": hi, **_summary(selected)})
        sex_totals = [{"sex": sex, **_summary([r for r in output if r["sex"] == sex])}
                      for sex in sorted({r["sex"] for r in output})]
        if sum(decimal(r["population"]) for r in display) != decimal(total["population"]):
            raise DataValidationError("Display-band population does not reconcile.")
        warnings = list(demography.get("warnings", []))
        if decimal(total["unestimated_population"]) > 0:
            warnings.append("open_age_tail_unestimated_no_full_population_point_prevalence")
        if "incomplete_age_coverage" in warnings:
            total["any_tbi_prevalence_full"] = None
            warnings.append("source_demography_has_age_gaps_population_scope_is_supplied_records_only")
        if time == 2024:
            warnings.append("integer_2024_time_convention_source_midyear_counts_retained_without_redating")
        return artifact(PROFILE, geography_id=demography["geography_id"], baseline_year=2024,
            baseline_time=numeric_text(time), recent_window_years=numeric_text(window),
            recency_clock="most_recent_infection_time", target_state="any_tbi",
            remote_only_rule="prior_infection_and_no_event_in_recent_window", scenario_id=scenario_id,
            status="complete_age_coverage" if decimal(total["unestimated_population"]) == 0 and "incomplete_age_coverage" not in warnings else "partial_age_coverage",
            population_scope="supplied_demographic_records_only", prevalence_unit="proportion_0_1",
            alpha=calibration["alpha"], alpha_unit="per_person_year",
            calibration=deepcopy(calibration), calibration_hash=calibration["profile_hash"],
            source_demography=deepcopy(demography), demography_hash=demography["profile_hash"],
            rows=output, display_bands=display, sex_totals=sex_totals, combined_totals=total,
            derived_both_display_bands=[{'sex': 'derived_both', 'age_lower': lo, 'age_upper_exclusive': hi,
                **_summary([r for r in output if r['age_lower'] >= lo and (hi is None or
                    (r['age_upper_exclusive'] is not None and r['age_upper_exclusive'] <= hi))])} for lo,hi in DISPLAY_BANDS],
            published_both_reconciliation=deepcopy(demography.get('published_both_checks', [])),
            cross_section_assumptions=['alive_in_2024', 'infection_independent_survival', 'no_historical_deaths_simulated',
                'migration_histories_not_modelled_infection_may_have_been_acquired_elsewhere'],
            configuration={"age_policy": age_policy, "age_review": deepcopy(age_review), "tail_distributions": tails, "scenario_id": scenario_id},
            uncertainty_interpretation="parameter_endpoints_are_scenarios_not_complete_probability_intervals_remote_only_is_not_monotone_in_alpha",
            warnings=sorted(set(warnings)))


def restore_tbi_profile(payload):
    """Rebuild from embedded hashed original inputs, never from a latest source."""
    value = json.loads(payload) if isinstance(payload, (str, bytes)) else deepcopy(payload)
    verify(value, PROFILE)
    rebuilt = baseline_tbi_profile(value["source_demography"], value["calibration"],
        recent_window=value["recent_window_years"], baseline_time=value["baseline_time"], **value["configuration"])
    if rebuilt["profile_hash"] != value["profile_hash"]:
        raise DataValidationError("Saved TBI profile does not reproduce from its original scientific inputs.")
    return value


def uncertainty_scenario_set(profiles, *, review, weights=None):
    reviewed(review)
    if not profiles:
        raise DataValidationError("Explicit named scenario profiles required.")
    for name, profile in profiles.items():
        if not name:
            raise DataValidationError("Scenario name required.")
        verify(profile, PROFILE)
    if len({(p["geography_id"], p["baseline_time"], p["demography_hash"]) for p in profiles.values()}) != 1:
        raise DataValidationError("Scenario comparisons require the same geography, time and demographic weights.")
    if weights is not None:
        if set(weights) != set(profiles):
            raise DataValidationError("Scenario weights must cover all named profiles.")
        with scientific_context():
            if sum(probability(v) for v in weights.values()) != 1:
                raise DataValidationError("Scenario weights must sum to one without clipping.")
    return artifact("country_tb.uncertainty_scenario_set.v1",
        uncertainty_kind="explicit_scenarios_not_population_stochastic_variation_or_posterior",
        uncertainty_components={
            'parameter': 'scale_assay_performance_and_anchor_endpoints_when_explicitly_supplied',
            'source_data': 'frozen_source_identity_and_coherent_WHO_bound_trajectories',
            'scenario': 'pressure_shape_and_historical_backcast_choices',
            'stochastic': 'none_in_deterministic_profile_separate_explicit_RNG_validation_only'},
        review=deepcopy(review), rows=[{"scenario_id": name, "profile_hash": profiles[name]["profile_hash"],
            "profile": deepcopy(profiles[name]), "weight": numeric_text(weights[name]) if weights is not None else None}
            for name in sorted(profiles)], rng=None)


def preview_files(profile):
    verify(profile, PROFILE)
    metadata = {"model_exclusions": ";".join(EXCLUSIONS), "infection_state_semantics": profile["infection_state_semantics"],
                "interpretation": profile["interpretation"]}
    rows = [{**r, "parameter_endpoint_scenarios": json.dumps(r["parameter_endpoint_scenarios"], sort_keys=True),
             **metadata} for r in profile["rows"]]
    bands = [{**r, "missing_age_coverage_bounds_not_confidence_interval": json.dumps(
        r["missing_age_coverage_bounds_not_confidence_interval"], sort_keys=True), **metadata} for r in profile["display_bands"]]
    report = artifact("country_tb.catalytic_preview_diagnostics.v1", profile_hash_reference=profile["profile_hash"],
        combined_totals=profile["combined_totals"], sex_totals=profile["sex_totals"],
        calibration_diagnostics=profile["calibration"]["diagnostics"], rows=profile["calibration"]["rows"],
        warnings=profile["warnings"])
    return {"profile.json": canonical(profile) + b"\n",
            "diagnostics.json": canonical(report) + b"\n",
            "age_sex_tbi.csv": csv_bytes(rows, tuple(rows[0])),
            "display_bands.csv": csv_bytes(bands, tuple(bands[0])),
            "README.md": (
                "# Isolated reinfection-capable TBI diagnostic\n\n"
                "Expected age/sex cross-section; no population allocation or projection.\n\n"
                "Scientific profile SHA-256: " + profile["profile_hash"] + "\n\n"
                "Recent reinfection overrides remote-only classification; prior remote exposure remains auditable.\n"
                "Exclusions: " + "; ".join(EXCLUSIONS) + ".\n\n"
                "Alpha: " + profile["alpha"] + "; origin: " + profile["calibration"]["fit_status"] +
                "; age coverage: " + profile["status"] + ".\n\n"
                "Pressure is a reviewed proxy, not disease-calibrated force of infection.\n"
                "Parameter endpoints are scenarios, not confidence intervals for recency.\n"
            ).encode("utf-8")}
