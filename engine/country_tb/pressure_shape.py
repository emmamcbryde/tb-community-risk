"""Offline incidence descriptors, explicit normalized shapes and backcasts.

No disease-to-infection conversion or disease calibration is performed.
"""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal

from engine.country_tb.calibration_common import scientific_context, artifact, finite_integer, reviewed, verify
from engine.country_tb.data_common import DataValidationError, decimal, numeric_text

SHAPE = "country_tb.infection_pressure_shape.v1"
HISTORY = "country_tb.infection_pressure_history.v1"
SENSITIVITY_SPECIFICATIONS = (
    "proportional_unsmoothed", "covid_aware_piecewise_reviewed_external_shape",
    "user_supplied_pressure_history",
)


def _positive(value):
    result = decimal(value)
    if result <= 0:
        raise DataValidationError("Log-incidence requires strictly positive values; no epsilon or imputation.")
    return result


def shape_from_incidence(incidence, *, review, normalization_year=2024,
                         method="smoothed_log_local_linear", window_years=5, source_bound='central'):
    """Centred local-linear OLS of log incidence, fixed five-year default window.

    Edge windows truncate explicitly (3 points at a five-year-window endpoint).
    All internal years must be present and strictly positive. Raw records/bounds
    are retained verbatim. No annual uncertainty independence is assumed.
    """
    verify(incidence, "country_tb.incidence_history.v1")
    review = reviewed(review, required=False)
    if source_bound not in {'lower', 'central', 'upper'}:
        raise DataValidationError('Explicit WHO lower/central/upper trajectory required.')
    finite_integer(normalization_year, minimum=1)
    if method not in {"smoothed_log_local_linear", "proportional_unsmoothed"}:
        raise DataValidationError("Unsupported shape method; COVID-aware alternatives require a supplied reviewed shape.")
    finite_integer(window_years, minimum=5)
    if window_years % 2 == 0:
        raise DataValidationError("Local-linear window must be odd.")
    source = sorted(deepcopy(incidence["rows"]), key=lambda r: r["year"])
    years = [finite_integer(row["year"], minimum=1) for row in source]
    if not years or len(set(years)) != len(years) or max(years) > 2024:
        raise DataValidationError("Empty, duplicated or post-2024 incidence years.")
    invalid = []
    for row in source:
        try:
            _positive(row[source_bound])
        except (DataValidationError, KeyError):
            invalid.append(row['year'])
    if invalid:
        label = 'No usable incidence series' if len(invalid) == len(source) else 'Invalid log-incidence series'
        raise DataValidationError(f'{label}; offending years: {invalid}. Missing/zero/nonfinite is not imputed; supply a reviewed alternative history.')
    if years != list(range(min(years), max(years) + 1)):
        raise DataValidationError(f"Missing incidence years: {sorted(set(range(min(years), max(years)+1))-set(years))}; no silent interpolation.")
    if normalization_year not in years:
        raise DataValidationError("Normalization year must be an observed year.")
    if len(source) < (window_years if method == "smoothed_log_local_linear" else 2):
        raise DataValidationError("Insufficient incidence series length for the fixed method.")
    with scientific_context():
        logs = []
        for row in source:
            if row["central"] is None:
                raise DataValidationError("Missing incidence observation; provide a separate reviewed shape instead.")
            central = _positive(row["central"])
            lower, upper = decimal(row["lower"]), decimal(row["upper"])
            if not lower <= central <= upper:
                raise DataValidationError("Invalid source uncertainty ordering.")
            if row.get("quantity") != "estimated_tb_incidence":
                raise DataValidationError("Notifications cannot silently become incidence.")
            logs.append(_positive(row[source_bound]).ln())
        fits, diagnostics = [], []
        radius = window_years // 2
        for i, year in enumerate(years):
            indices = [j for j, y in enumerate(years) if abs(y - year) <= radius]
            if method == "proportional_unsmoothed":
                fitted, slope = logs[i], None
                indices = [i]
            else:
                if len(indices) < 3:
                    raise DataValidationError("Insufficient edge observations.")
                x = [Decimal(years[j] - year) for j in indices]
                y = [logs[j] for j in indices]
                xm, ym = sum(x) / len(x), sum(y) / len(y)
                slope = sum((xx - xm) * (yy - ym) for xx, yy in zip(x, y)) / sum((xx - xm) ** 2 for xx in x)
                fitted = ym - slope * xm
            fits.append(fitted)
            diagnostics.append({
                "year": year, "observed_incidence": row_value(source[i][source_bound]),
                "fitted_incidence": numeric_text(fitted.exp()),
                "log_residual": numeric_text(logs[i] - fitted),
                "incidence_residual": numeric_text(decimal(source[i][source_bound]) - fitted.exp()),
                "local_log_slope_per_year": numeric_text(slope) if slope is not None else None,
                "window_years_used": [years[j] for j in indices],
                "truncated_edge": method != "proportional_unsmoothed" and len(indices) < window_years,
                "covid_period": 2020 <= year <= 2022,
            })
        reference = fits[years.index(normalization_year)]
        rows = [{"year": year, "start": str(year), "end": str(year + 1),
                 "shape": "1" if year == normalization_year else numeric_text((fit - reference).exp()),
                 "fitted_incidence": numeric_text(fit.exp())}
                for year, fit in zip(years, fits)]
        for i, row in enumerate(diagnostics):
            row["fitted_log_change_from_prior_year"] = numeric_text(fits[i] - fits[i-1]) if i else None
        return artifact(SHAPE, kind="infection_pressure_shape_proxy", unit="dimensionless_shape", geography_id=incidence["geography_id"],
            sex_policy="common_curve_no_acquisition_difference_inferred_from_population_counts",
            review=review, source_incidence=deepcopy(incidence), source_incidence_hash=incidence["profile_hash"],
            normalization={"year": normalization_year, "fitted_incidence_value": numeric_text(reference.exp()),
                           "shape_value": "1", "rule": "divide_fitted_incidence_by_fitted_reference_year"},
            method=method, source_bound=source_bound, parameters={"window_years": window_years if method != "proportional_unsmoothed" else None,
                "fit": "unweighted_OLS_log_incidence", "edge_rule": "truncate_window_minimum_three",
                "missing_policy": "reject", "zero_policy": "reject", "post2024_extrapolation": "forbidden"},
            uncertainty_policy="WHO_bounds_retained_not_pressure_confidence_intervals",
            rows=rows, diagnostics=diagnostics, missing_years=[],
            sensitivity_specifications=list(SENSITIVITY_SPECIFICATIONS))


def row_value(value):
    return numeric_text(decimal(value))


def _normalize_intervals(rows):
    normalized = []
    for row in rows:
        if set(row) != {"start", "end", "shape"}:
            raise DataValidationError("Shape rows require exactly start/end/shape; supply uncertainty as separate reviewed scenarios.")
        lo, hi = decimal(row["start"]), decimal(row["end"])
        value = decimal(row["shape"])
        if hi <= lo or hi > 2025:
            raise DataValidationError("Invalid interval or pressure after the 2024 observation year.")
        normalized.append({"start": numeric_text(lo), "end": numeric_text(hi), "shape": numeric_text(value)})
    normalized.sort(key=lambda r: decimal(r["start"]))
    if not normalized or any(decimal(a["end"]) != decimal(b["start"]) for a, b in zip(normalized, normalized[1:])):
        raise DataValidationError("Shape intervals overlap or have missing history.")
    return normalized


def supplied_shape(rows, *, geography_id, review, provenance, normalization_time,
                   quantity="relative_shape", sex_rows=None):
    """Explicit reviewed alternative; no invented age/sex or annual extrapolation.

    All sex curves use the SAME reference value. This preserves supplied sex
    contrasts rather than normalizing them away. Common curve is mandatory;
    a sex-specific override must cover identical time intervals.
    """
    review = reviewed(review)
    if not geography_id or not provenance or quantity not in {"relative_shape", "infection_hazard_per_person_year"}:
        raise DataValidationError("Explicit geography, source and input quantity required.")
    with scientific_context():
        common = _normalize_intervals(rows)
        when = decimal(normalization_time)
        reference = next((decimal(r["shape"]) for r in common if decimal(r["start"]) <= when < decimal(r["end"])), None)
        if reference is None or reference <= 0:
            raise DataValidationError("Normalization time needs a strictly positive supplied value.")
        curves = {"common": common}
        for sex, values in (sex_rows or {}).items():
            if sex not in {"female", "male", "both", "unknown"}:
                raise DataValidationError("Invalid explicit sex curve.")
            curve = _normalize_intervals(values)
            if [(r["start"], r["end"]) for r in curve] != [(r["start"], r["end"]) for r in common]:
                raise DataValidationError("Sex curves must have identical explicit coverage.")
            curves[sex] = curve
        normalized = {sex: [{**r, "shape": numeric_text(decimal(r["shape"]) / reference)} for r in curve]
                      for sex, curve in curves.items()}
        return artifact(SHAPE, kind="infection_pressure_shape_proxy", unit="dimensionless_shape", geography_id=geography_id,
            method="user_supplied_pressure_history", review=review, provenance=deepcopy(provenance),
            original_curves=curves, original_quantity=quantity, curves=normalized,
            normalization={"time": numeric_text(when), "supplied_value": numeric_text(reference), "shape_value": "1",
                           "rule": "common_reference_value_for_all_sexes"},
            sex_policy="explicit_sex_curves" if sex_rows else "common_curve",
            parameters={"post2024_extrapolation": "forbidden"},
            uncertainty_policy="supplied_scenarios_only_no_automatic_distribution")


def apply_backcast(shape, *, start_time, method="constant_earliest", parameters=None):
    """Approved constant central policy; bounded/supplied histories are explicit.

    Anchor-informed history is supplied as a fixed reviewed sensitivity, not
    fitted alongside alpha. No second calibration parameter is introduced.
    """
    verify(shape, SHAPE)
    reviewed(shape["review"])
    parameters = deepcopy(parameters or {})
    with scientific_context():
        start = decimal(start_time)
        curves = shape.get("curves") or {"common": [{k: row[k] for k in ("start", "end", "shape")} for row in shape["rows"]]}
        first = decimal(curves["common"][0]["start"])
        if start > first:
            raise DataValidationError("Backcast start must not remove observed history.")
        if method not in {"constant_earliest", "bounded_backward_trend", "supplied_historical_shape", "anchor_informed_fixed_shape"}:
            raise DataValidationError("Unsupported backcast; scenario uncertainty requires separate named histories.")
        if method == "constant_earliest" and parameters:
            raise DataValidationError("Constant backcast has no slope/bound parameters.")
        if method == "bounded_backward_trend":
            reviewed(parameters.get("review"))
            if not {"slope_per_year", "minimum_shape", "maximum_shape"} <= parameters.keys():
                raise DataValidationError("Bounded backcast requires slope and finite limits.")
        result = {}
        for sex, source in sorted(curves.items()):
            early = decimal(source[0]["shape"])
            before = []
            if start < first and method in {"constant_earliest", "bounded_backward_trend"}:
                slope = Decimal(0) if method == "constant_earliest" else decimal(parameters["slope_per_year"], nonnegative=False)
                lower = decimal(parameters["minimum_shape"]) if method != "constant_earliest" else Decimal(0)
                upper = decimal(parameters["maximum_shape"]) if method != "constant_earliest" else None
                if upper is not None and upper < lower:
                    raise DataValidationError("Inverted explicit backcast bounds.")
                if upper is not None and not lower <= early <= upper:
                    raise DataValidationError("Earliest shape violates the backcast join bounds.")
                cursor = start
                while cursor < first:
                    end = min(first, cursor.to_integral_value(rounding="ROUND_FLOOR") + 1)
                    # Anchor annual bins to the calendar grid, not the requested
                    # oldest birth. A fractional start must not alter later q.
                    grid_year = cursor.to_integral_value(rounding="ROUND_FLOOR")
                    level = early * (slope * (grid_year - first)).exp()
                    if level < lower or (upper is not None and level > upper):
                        raise DataValidationError("Backcast violates explicit bounds; no cap applied.")
                    before.append({"start": numeric_text(cursor), "end": numeric_text(end), "shape": numeric_text(level)})
                    cursor = end
            elif start < first:
                reviewed(parameters.get("review"))
                historical = _normalize_intervals(parameters["curves"].get(sex, parameters["curves"].get("common", [])))
                if decimal(historical[0]["start"]) != start or decimal(historical[-1]["end"]) != first:
                    raise DataValidationError("Supplied history must cover exactly the required backcast.")
                if parameters.get("join_rule") not in {"continuous", "reviewed_discontinuity"}:
                    raise DataValidationError("Explicit reviewed historical join rule required.")
                if parameters["join_rule"] == "continuous" and decimal(historical[-1]["shape"]) != early:
                    raise DataValidationError("Historical shape has an incompatible join.")
                before = historical
            result[sex] = _normalize_intervals(before + deepcopy(source))
        policy = artifact('country_tb.historical_backcast_policy.v1', method=method,
            start_time=numeric_text(start), first_observed_year=numeric_text(first), parameters=parameters,
            reconstructed_period=[numeric_text(start), numeric_text(first)] if start < first else None,
            bounds={'minimum_shape': parameters.get('minimum_shape'), 'maximum_shape': parameters.get('maximum_shape')},
            warnings=['pre_observation_exposure_is_assumed_not_identified'] if start < first else [],
            provenance_hash=shape['profile_hash'])
        return artifact(HISTORY, kind="reviewed_relative_pressure_history", unit="dimensionless_shape",
            geography_id=shape["geography_id"], review=deepcopy(shape["review"]),
            shape=deepcopy(shape), shape_hash=shape["profile_hash"], curves=result,
            backcast_policy=policy, backcast_policy_hash=policy['profile_hash'],
            backcast={"method": method, "start_time": numeric_text(start), "first_source_time": numeric_text(first),
                      "parameters": parameters, "central_policy": "constant_earliest", "applied": start < first},
            future_policy="none_no_post2024_extrapolation", age_exposure_policy="calendar_pressure_only",
            uncertainty_policy="conditional_on_fixed_reviewed_history")


def validate_history(history):
    verify(history, HISTORY)
    reviewed(history["review"])
    verify(history["shape"], SHAPE)
    verify(history['backcast_policy'], 'country_tb.historical_backcast_policy.v1')
    if history['backcast_policy_hash'] != history['backcast_policy']['profile_hash']:
        raise DataValidationError('Backcast policy hash mismatch.')
    if (history["shape_hash"] != history["shape"]["profile_hash"] or history["unit"] != "dimensionless_shape"
            or history["geography_id"] != history["shape"]["geography_id"]):
        raise DataValidationError("Pressure-history source reference/unit mismatch.")
    if "source_incidence" in history["shape"]:
        source = history["shape"]["source_incidence"]
        verify(source, "country_tb.incidence_history.v1")
        if source["profile_hash"] != history["shape"]["source_incidence_hash"]:
            raise DataValidationError("Pressure proxy incidence reference mismatch.")
    curves = history.get("curves")
    if not isinstance(curves, dict) or "common" not in curves or set(curves) - {"common", "female", "male", "both", "unknown"}:
        raise DataValidationError("Invalid pressure-history sex curves.")
    for rows in curves.values():
        _normalize_intervals(rows)
    return history
