"""Reviewed assay operators and immutable survey/override observations."""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal

from engine.country_tb.calibration_common import scientific_context, artifact, finite_integer, probability, reviewed, verify
from engine.country_tb.data_common import DataValidationError, canonical, decimal, numeric_text, digest

ASSAY = "country_tb.tbi_assay_operator.v1"
OBSERVATIONS = "country_tb.tbi_observations.v1"


def assay_operator(*, assay_class, name, positivity_threshold, sensitivity,
                   specificity, review, provenance):
    """No evidence defaults. Sensitivity/specificity are supplied and reviewed.

    Identity observation is restricted to explicitly synthetic assay records.
    This affine operator assumes no clearance/reversion under this milestone.
    """
    reviewed(review)
    if assay_class not in {"TST", "IGRA", "TBST", "synthetic_identity", "direct_latent_tbi"}:
        raise DataValidationError("Explicit supported assay class required.")
    if not name or not positivity_threshold or not provenance:
        raise DataValidationError("Assay/product, threshold and source provenance are required.")
    with scientific_context():
        se, sp = probability(sensitivity), probability(specificity)
        if se + sp <= 1:
            raise DataValidationError("Uninformative/inverted assay cannot identify a monotone infection scale.")
        if assay_class in {"synthetic_identity", "direct_latent_tbi"} and (se != 1 or sp != 1):
            raise DataValidationError("Synthetic identity must have sensitivity/specificity exactly one.")
        if assay_class not in {"synthetic_identity", "direct_latent_tbi"} and se == sp == 1:
            raise DataValidationError("Perfect identity is not a default observation model for a real assay.")
        return artifact(ASSAY, assay_class=assay_class, name=name, positivity_threshold=deepcopy(positivity_threshold),
            sensitivity=numeric_text(se), specificity=numeric_text(sp),
            review=deepcopy(review), provenance=deepcopy(provenance),
            performance_provenance={'sensitivity': deepcopy(provenance.get('sensitivity', provenance)),
                                    'specificity': deepcopy(provenance.get('specificity', provenance))},
            observation_target="latent_any_tbi" if assay_class == 'direct_latent_tbi' else "assay_positivity_given_any_tbi")


def assay_probability(tbi_prevalence, operator):
    verify(operator, ASSAY)
    reviewed(operator["review"])
    with scientific_context():
        value = probability(tbi_prevalence)
        se, sp = probability(operator["sensitivity"]), probability(operator["specificity"])
        if se + sp <= 1:
            raise DataValidationError("Assay operator is not positively informative.")
        return se * value + (1 - sp) * (1 - value)


def _strata(record):
    lower = decimal(record["age_lower"])
    upper = None if record["age_upper_exclusive"] is None else decimal(record["age_upper_exclusive"])
    if upper is not None and upper <= lower:
        raise DataValidationError("Invalid observation age range.")
    sex = record["sex"]
    if sex not in {"female", "male", "both", "unknown"} or not record.get("weight_source"):
        raise DataValidationError("Explicit observation sex and age/sex weight provenance required.")
    result, keys = [], set()
    if not isinstance(record['strata'], list):
        raise DataValidationError('Observation strata must be an explicit list.')
    for source in record["strata"]:
        if not isinstance(source, dict) or not {'age_lower','age_upper_exclusive','sex','age_policy','weight'} <= source.keys():
            raise DataValidationError('Incomplete observation weight stratum.')
        row = deepcopy(source)
        lo = decimal(row["age_lower"])
        hi = None if row["age_upper_exclusive"] is None else decimal(row["age_upper_exclusive"])
        if lo < lower or (upper is not None and (lo >= upper or (hi is not None and hi > upper))):
            raise DataValidationError("Observation weights fall outside the reported age range.")
        if hi is not None and hi <= lo:
            raise DataValidationError("Invalid weighted age interval.")
        if row["sex"] not in {"female", "male", "both", "unknown"}:
            raise DataValidationError("Invalid weighted sex.")
        if sex != "both" and row["sex"] != sex:
            raise DataValidationError("Sex-specific observation has mismatched weights.")
        if sex == "both" and row["sex"] == "unknown":
            raise DataValidationError("Unknown sex cannot silently become a both-sex observation.")
        if row["age_policy"] == "exact_age" and hi is not None:
            raise DataValidationError("Exact-age observation cannot also specify an interval.")
        if row["age_policy"] == "uniform_within_single_year" and (hi is None or hi - lo != 1):
            raise DataValidationError("Single-year observation weights need a one-year interval.")
        if row["age_policy"] not in {"exact_age", "uniform_within_single_year", "reviewed_uniform_within_interval"}:
            raise DataValidationError("Unsupported observation age weighting.")
        if row["age_policy"] == "reviewed_uniform_within_interval":
            reviewed(row.get("age_review"))
        weight = decimal(row["weight"])
        if weight <= 0:
            raise DataValidationError("Observation weights must be positive; omit zero-weight cells explicitly.")
        key = (numeric_text(lo), numeric_text(hi) if hi is not None else None, row["sex"], row["age_policy"])
        if key in keys:
            raise DataValidationError("Duplicate observation weight stratum.")
        keys.add(key)
        result.append({**row, "age_lower": numeric_text(lo), "age_upper_exclusive": numeric_text(hi) if hi is not None else None,
                       "weight": numeric_text(weight)})
    if not result:
        raise DataValidationError("Explicit age/sex weights required; none inferred from population size.")
    if "both" in {r["sex"] for r in result} and len({r["sex"] for r in result}) != 1:
        raise DataValidationError("Both-sex and sex-specific weights would double-count.")
    # Age intervals may not overlap within a sex; exact ages cannot duplicate/lie
    # within an interval. No assumption of independent survey rows is introduced.
    for s in {r["sex"] for r in result}:
        values = sorted((r for r in result if r["sex"] == s), key=lambda r: decimal(r["age_lower"]))
        for left, right in zip(values, values[1:]):
            lo = decimal(left["age_lower"])
            hi = decimal(left["age_upper_exclusive"]) if left["age_upper_exclusive"] is not None else lo
            if decimal(right["age_lower"]) < hi or decimal(right["age_lower"]) == lo:
                raise DataValidationError("Overlapping observation weight strata.")
    return sorted(result, key=canonical)


def observations_profile(records, *, review, mode="survey", parent_profile=None):
    reviewed(review)
    if not isinstance(records, (list, tuple)):
        raise DataValidationError('Explicit observation records required.')
    if mode not in {"survey", "reviewed_profile_override"}:
        raise DataValidationError("Unsupported observation mode.")
    if mode == "reviewed_profile_override":
        verify(parent_profile, "country_tb.calibrated_tbi_age_profile.v1")
    elif parent_profile is not None:
        raise DataValidationError("Parent profile belongs only to an explicit override.")
    result, ids = [], set()
    with scientific_context():
        for original in records:
            required = {"observation_id", "geography_id", "age_lower", "age_upper_exclusive", "sex",
                        "survey_year", "survey_time", "numerator", "denominator", "prevalence",
                        "review", "source", "sampling_frame", "survey_design", "uncertainty",
                        "assay_operator", "weight_source", "strata"}
            if not isinstance(original, dict) or not required <= original.keys():
                raise DataValidationError("Incomplete TBI observation contract; missing is explicit, not omitted.")
            if {"recent_prevalence", "remote_prevalence", "recent", "remote_only"} & original.keys():
                raise DataValidationError("Recency is derived from infection history, not supplied as an anchor prevalence.")
            row = deepcopy(original)
            if row.get('geography_level') not in {'national', 'subnational', 'custom'} or not row.get('observation_policy'):
                raise DataValidationError('Explicit geography level and observation policy required.')
            identifier = row.get("observation_id")
            if not identifier or identifier in ids:
                raise DataValidationError("Observation IDs must be nonempty and unique.")
            ids.add(identifier)
            reviewed(row.get("review"))
            if not row.get("geography_id") or not row.get("source") or not row.get("sampling_frame") or not row.get("survey_design"):
                raise DataValidationError("Observation geography/source/frame/design required.")
            year = finite_integer(row["survey_year"], minimum=1)
            time = decimal(row["survey_time"])
            if year > 2024 or not year <= time < year + 1:
                raise DataValidationError("Survey time must match its source year and not project beyond 2024.")
            value = probability(row["prevalence"])
            numerator, denominator = row["numerator"], row["denominator"]
            if (numerator is None) != (denominator is None):
                raise DataValidationError("Numerator and denominator are jointly known or jointly missing.")
            if numerator is not None:
                finite_integer(numerator)
                finite_integer(denominator, minimum=1)
                if numerator > denominator:
                    raise DataValidationError("Numerator exceeds denominator.")
                tolerance = decimal(row.get("count_prevalence_tolerance", "0"))
                discrepancy = Decimal(numerator) / denominator - value
                if abs(discrepancy) > tolerance and not row.get("weighted_count_discrepancy_rationale"):
                    raise DataValidationError("Counts/prevalence disagree; preserve and explicitly review weighting/rounding.")
                row["count_prevalence_difference"] = numeric_text(discrepancy)
            verify(row["assay_operator"], ASSAY)
            reviewed(row["assay_operator"]["review"])
            uncertainty = deepcopy(row["uncertainty"])
            if not isinstance(uncertainty, dict):
                raise DataValidationError('Explicit uncertainty object required.')
            if uncertainty.get("kind") == "interval":
                if not {'lower','upper','method'} <= uncertainty.keys():
                    raise DataValidationError('Uncertainty interval requires lower, upper and method.')
                lo, hi = probability(uncertainty["lower"]), probability(uncertainty["upper"])
                if not lo <= value <= hi or not uncertainty.get("method"):
                    raise DataValidationError("Invalid observation interval/method.")
                uncertainty.update(lower=numeric_text(lo), upper=numeric_text(hi))
                if "level" in uncertainty and not 0 < probability(uncertainty["level"]) < 1:
                    raise DataValidationError("Invalid uncertainty level.")
            elif uncertainty.get("kind") not in {"exact", "unspecified"}:
                raise DataValidationError("Uncertainty must be exact, supplied interval, or explicitly unspecified.")
            if mode == "reviewed_profile_override" and (row["geography_id"] != parent_profile["geography_id"] or year != 2024):
                raise DataValidationError("Override must retain the parent geography and baseline year.")
            result.append({**row, "survey_time": numeric_text(time), "prevalence": numeric_text(value),
                'provenance_hash': digest(canonical({'source': row['source'], 'review': row['review'],
                    'observation_policy': row['observation_policy'], 'assay_hash': row['assay_operator']['profile_hash']})),
                "age_lower": numeric_text(row["age_lower"]),
                "age_upper_exclusive": numeric_text(row["age_upper_exclusive"]) if row["age_upper_exclusive"] is not None else None,
                "strata": _strata(row), "uncertainty": uncertainty,
                "receipt": {"source_record": deepcopy(original), "source_record_sha256": digest(canonical(original))}})
    if not result:
        raise DataValidationError("An observation artifact must contain observations.")
    return artifact(OBSERVATIONS, mode=mode, prevalence_unit="proportion_0_1", review=deepcopy(review), rows=sorted(result, key=lambda r: r["observation_id"]),
        parent_profile_hash=parent_profile["profile_hash"] if parent_profile is not None else None,
        source_preservation="immutable_original_records_retained; override_is_a_new_anchor_not_display_replacement")


def reviewed_tbi_override(parent_profile, records, *, review):
    """Separate reviewed age-specific total-TBI artifact; no invented recency."""
    observations = observations_profile(records, review=review, mode='reviewed_profile_override', parent_profile=parent_profile)
    intervals = {}
    for row in observations['rows']:
        if row['assay_operator']['assay_class'] not in {'direct_latent_tbi', 'synthetic_identity'}:
            raise DataValidationError('A direct TBI override needs reviewed latent-TBI values; assay positivity is a calibration observation.')
        key = row['sex']
        lo = decimal(row['age_lower'])
        hi = decimal(row['age_upper_exclusive']) if row['age_upper_exclusive'] is not None else Decimal('Infinity')
        if any(lo < b and a < hi for a,b in intervals.get(key, [])):
            raise DataValidationError('Overlapping age-specific override rows.')
        intervals.setdefault(key, []).append((lo,hi))
    if 'both' in intervals and len(intervals) > 1:
        raise DataValidationError('Override both-sex rows cannot overlap separate sex strata.')
    return artifact('country_tb.reviewed_tbi_override.v1', parent_profile_hash=parent_profile['profile_hash'],
        geography_id=parent_profile['geography_id'], baseline_year=2024, observations=observations,
        observations_hash=observations['profile_hash'], review=deepcopy(review),
        rows=[{'observation_id': r['observation_id'], 'age_lower': r['age_lower'],
            'age_upper_exclusive': r['age_upper_exclusive'], 'sex': r['sex'], 'any_tbi': r['prevalence'],
            'recent': None, 'remote_only': None, 'prior_remote_plus_recent': None,
            'most_recent_infection_time': None, 'composition_status': 'unestimated_total_only_override'}
            for r in observations['rows']], source_preservation='separate_artifact_parent_and_sources_unchanged')
