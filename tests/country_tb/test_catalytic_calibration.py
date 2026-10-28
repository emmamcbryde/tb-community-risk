"""Focused Milestone 2C checks using invented shapes/anchors, not country evidence."""
from copy import deepcopy
from decimal import Decimal, localcontext
import ast
import csv
import io
import json
import math
import random
import socket
import subprocess

import pytest

from engine.country_tb.calibration_common import EXCLUSIONS, STATE_FIELDS, validate_states
from engine.country_tb.catalytic_calibration import calibration_policy, calibrate_scale
from engine.country_tb.catalytic_core import evaluate_age as _evaluate_age

def evaluate_age(*args, **kwargs):
    return _evaluate_age(*args, include_event_times=False, **kwargs)

from engine.country_tb.data_common import DataValidationError, REPO_ROOT, canonical, digest, profile_hash, seal
from engine.country_tb.demography import demographic_profile
from engine.country_tb.incidence import WHO_2025, load_who_snapshot, national_profile, user_series_profile
from engine.country_tb.pressure_shape import shape_from_incidence, supplied_shape, apply_backcast
from engine.country_tb.tbi_baseline import baseline_tbi_profile, preview_files, restore_tbi_profile, uncertainty_scenario_set
from engine.country_tb.tbi_observations import assay_operator, assay_probability, observations_profile

REVIEW = {"status": "reviewed", "source_id": "synthetic:M2C", "rationale": "Invented software-validation fixture only"}
GEO = "custom:synthetic"


def shape(rows=None, sex_rows=None):
    return supplied_shape(rows or [{"start": "2000", "end": "2025", "shape": "1"}],
        geography_id=GEO, review=REVIEW, provenance={"source_id": "synthetic:M2C"},
        normalization_time="2024", sex_rows=sex_rows)


def history(**kwargs):
    return apply_backcast(shape(**kwargs), start_time="1900")


def operator(kind="synthetic_identity", se="1", sp="1"):
    return assay_operator(assay_class=kind, name=kind, positivity_threshold={"value": "synthetic", "unit": "fixture"},
        sensitivity=se, specificity=sp, review=REVIEW, provenance={"source_id": "synthetic:M2C"})


def observation(h, age=20, sex="female", alpha=".01", identifier=None, assay=None, prevalence=None, uncertainty=None):
    assay = assay or operator()
    value = evaluate_age(h, alpha=alpha, age=age, sex=sex)["any_tbi"]
    return {"observation_id": identifier or f"synthetic:{sex}:{age}", "geography_id": GEO, "geography_level": "custom",
        "observation_policy": "reviewed_direct_or_fixed_assay", "age_lower": str(age), "age_upper_exclusive": str(age + 1), "sex": sex,
        "survey_year": 2024, "survey_time": "2024", "numerator": None, "denominator": None,
        "prevalence": str(assay_probability(value, assay)) if prevalence is None else str(prevalence),
        "review": deepcopy(REVIEW), "source": {"source_id": "synthetic:M2C"}, "sampling_frame": "invented exact ages",
        "survey_design": {"kind": "synthetic_only"}, "uncertainty": uncertainty or {"kind": "exact"},
        "assay_operator": assay, "weight_source": {"source_id": "synthetic:M2C"},
        "strata": [{"age_lower": str(age), "age_upper_exclusive": None, "age_policy": "exact_age",
                    "sex": sex, "weight": "1"}]}


def policy(records=(), **updates):
    arguments = {"alpha_upper": "1", "observation_weights": {r["observation_id"]: "1" for r in records},
                 "residual_tolerance": ".001", "review": REVIEW}
    arguments.update(updates)
    return calibration_policy(**arguments)


def fit(h, records):
    return calibrate_scale(h, policy(records, method="weighted_least_squares" if len(records)>1 else "overall_anchor"), observations=observations_profile(records, review=REVIEW))


def supplied(h=None, alpha=".01", interval=None):
    return calibrate_scale(h or history(), policy(), supplied_alpha=alpha,
                           scale_provenance={"source_id": "synthetic:M2C"}, supplied_alpha_interval=interval)


@pytest.fixture(scope="module")
def demography():
    rows = []
    for sex, count in (("female", 10), ("male", 20)):
        for age in range(101):
            rows.append({"source_id": "synthetic", "source_vintage": "v1", "geography_id": GEO,
                "geography_level": "custom", "year": 2024, "estimate_status": "estimate",
                "projection_scenario": "historical", "sex": sex, "age_lower": age,
                "age_upper_exclusive": age + 1 if age < 100 else None,
                "unit": "persons", "population_variant": "synthetic_residents", "population": str(count)})
    return demo_for_rows(rows)


def demo_for_rows(rows):
    candidate = {"snapshot_id": "synthetic-demo", "rows": rows, "input_source_hash": "0"*64,
        "transformed_file_hash": "1"*64, "manifest": {"source_id": "synthetic", "source_category": "user",
        "source_vintage": "v1", "review_status": "reviewed", "manifest_hash": "2"*64,
        "schema_version": "country_tb.demography_snapshot.v1", "transformation_version": "synthetic_fixture",
        "transformations": ["invented uniform age counts"], "source_url": "urn:synthetic", "import_date": "2026-10-28",
        "terms": {"status": "synthetic"}, "citation": "Invented software data", "projection_boundary": None}}
    return demographic_profile(GEO, [candidate])


def incidence():
    with localcontext() as context:
        context.prec = 80
        rows = []
        for year in range(2000, 2025):
            rate = Decimal(100) * (Decimal(".02") * (year - 2000)).exp()
            rows.append({"year": year, "quantity": "estimated_tb_incidence", "central": str(rate),
                "lower": str(rate / 2), "upper": str(rate * 2), "numerator": None, "denominator": None,
                "original_rate": str(rate), "unit": "cases_per_100000_person_year", "status": "estimated"})
    return user_series_profile("custom_population_series", GEO, rows,
                               provenance={"source_id": "synthetic"}, review_status="reviewed")


def test_log_linear_smoother_normalization_raw_preservation_and_edges():
    source = incidence()
    before = deepcopy(source)
    result = shape_from_incidence(source, review=REVIEW)
    assert source == before and result["source_incidence"] == source
    assert result["kind"] == "infection_pressure_shape_proxy" and result["disease_forward_operator"] is None
    assert result["normalization"]["year"] == 2024 and result["rows"][-1]["shape"] == "1"
    for row, diagnostic in zip(result["rows"], result["diagnostics"]):
        assert float(row["shape"]) == pytest.approx(math.exp(.02 * (row["year"] - 2024)), rel=1e-14)
        assert float(diagnostic["local_log_slope_per_year"]) == pytest.approx(.02, abs=1e-14)
        assert abs(float(diagnostic["log_residual"])) < 1e-14
    assert result["diagnostics"][0]["window_years_used"] == [2000, 2001, 2002]
    assert result["diagnostics"][-1]["window_years_used"] == [2022, 2023, 2024]
    assert sum(r["covid_period"] for r in result["diagnostics"]) == 3


def test_shape_order_hash_alternatives_and_normalization_scale_invariance():
    source = incidence()
    other = deepcopy(source); other["rows"].reverse()
    a = shape_from_incidence(source, review=REVIEW)
    b = shape_from_incidence(other, review=REVIEW)
    assert a["profile_hash"] == b["profile_hash"]
    c = shape_from_incidence(source, review=REVIEW, normalization_year=2000)
    assert a["profile_hash"] != c["profile_hash"]
    h1, h2 = apply_backcast(a, start_time=1900), apply_backcast(c, start_time=1900)
    alpha2 = Decimal(".01") / Decimal(c["rows"][-1]["shape"])
    p1 = evaluate_age(h1, alpha=".01", age=30, sex="female")
    p2 = evaluate_age(h2, alpha=alpha2, age=30, sex="female")
    assert float(p1["any_tbi"]) == pytest.approx(float(p2["any_tbi"]), abs=1e-14)
    unsmoothed = shape_from_incidence(source, review=REVIEW, method="proportional_unsmoothed")
    assert unsmoothed["method"] == "proportional_unsmoothed"
    assert "covid_aware_piecewise_reviewed_external_shape" in a["sensitivity_specifications"]


@pytest.mark.parametrize("kind", ["zero", "negative", "nonfinite", "missing", "gap", "short", "future"])
def test_bad_incidence_log_inputs_fail_without_mutation(kind):
    source = incidence()
    if kind in {"zero", "negative", "nonfinite", "missing"}:
        source["rows"][5]["central"] = {"zero": "0", "negative": "-1", "nonfinite": "NaN", "missing": None}[kind]
    elif kind == "gap":
        source["rows"].pop(5)
    elif kind == "short":
        source["rows"] = source["rows"][-4:]
    else:
        source["rows"][-1]["year"] = 2025
    source = seal(source)
    original = deepcopy(source)
    with pytest.raises(DataValidationError):
        shape_from_incidence(source, review=REVIEW)
    assert source == original


@pytest.mark.parametrize("window", [0, 3, 4, 6, True])
def test_invalid_fixed_smoothing_window_fails(window):
    with pytest.raises(DataValidationError):
        shape_from_incidence(incidence(), review=REVIEW, window_years=window)


def test_real_frozen_who_vnm_and_prk_missingness():
    who = load_who_snapshot(WHO_2025)
    vnm = national_profile(who, "VNM")
    result = shape_from_incidence(vnm, review=REVIEW)
    assert len(result["rows"]) == 25 and result["source_incidence_hash"] == vnm["profile_hash"]
    assert result["source_incidence"]["rows"][-1]["upper"] == vnm["rows"][-1]["upper"]
    with pytest.raises(DataValidationError, match="No usable"):
        shape_from_incidence(national_profile(who, "PRK"), review=REVIEW)


def test_central_backcast_is_constant_source_preserving_and_deterministic():
    source = shape()
    before = deepcopy(source)
    h = apply_backcast(source, start_time=1900)
    assert source == before and h == apply_backcast(source, start_time=1900)
    assert h["backcast"]["method"] == "constant_earliest"
    assert all(r["shape"] == "1" for r in h["curves"]["common"])
    assert h["shape_hash"] == source["profile_hash"]


def test_bounded_and_supplied_backcasts_fail_instead_of_capping():
    source = shape()
    parameters = {"slope_per_year": "-.02", "minimum_shape": "0", "maximum_shape": "10", "review": REVIEW}
    h = apply_backcast(source, start_time=1900, method="bounded_backward_trend", parameters=parameters)
    assert float(h["curves"]["common"][0]["shape"]) == pytest.approx(math.exp(2))
    with pytest.raises(DataValidationError, match="bounds"):
        apply_backcast(source, start_time=1900, method="bounded_backward_trend",
                       parameters={**parameters, "maximum_shape": "2"})
    past = {"curves": {"common": [{"start": "1900", "end": "2000", "shape": "2"}]},
            "review": REVIEW, "join_rule": "reviewed_discontinuity"}
    fixed = apply_backcast(source, start_time=1900, method="anchor_informed_fixed_shape", parameters=past)
    assert fixed["curves"]["common"][0]["shape"] == "2"
    assert policy()["free_parameters"] == ["alpha"]
    with pytest.raises(DataValidationError, match="join"):
        apply_backcast(source, start_time=1900, method="supplied_historical_shape",
                       parameters={**past, "join_rule": "continuous"})


def test_bounded_history_is_coherent_when_oldest_birth_changes_fractionally():
    parameters = {"slope_per_year": "-.02", "minimum_shape": "0", "maximum_shape": "10", "review": REVIEW}
    a = apply_backcast(shape(), start_time="1900", method="bounded_backward_trend", parameters=parameters)
    b = apply_backcast(shape(), start_time="1900.5", method="bounded_backward_trend", parameters=parameters)
    pa = evaluate_age(a, alpha=".01", age="123.5", sex="female")
    pb = evaluate_age(b, alpha=".01", age="123.5", sex="female")
    assert {k: pa[k] for k in (*STATE_FIELDS, "any_tbi", "cumulative_hazard")} == {
        k: pb[k] for k in (*STATE_FIELDS, "any_tbi", "cumulative_hazard")}


@pytest.mark.parametrize("rows", [
    [{"start": 2000, "end": 2026, "shape": 1}],
    [{"start": 2000, "end": 2025, "shape": -1}],
    [{"start": 2000, "end": 2025, "shape": "NaN"}],
    [{"start": 2000, "end": 2010, "shape": 1}, {"start": 2011, "end": 2025, "shape": 1}],
    [{"start": 2000, "end": 2025, "shape": 1, "ignored": 5}],
])
def test_bad_supplied_shapes_rejected(rows):
    with pytest.raises(DataValidationError):
        shape(rows)


@pytest.mark.parametrize("age", [0, ".001", 1, 4, 5, 14, 40, 85, 124])
def test_constant_pressure_first_state_closed_forms_and_zero(age):
    h = history()
    p = evaluate_age(h, alpha=".01", age=age, sex="female")
    a = float(age)
    assert float(p["cumulative_hazard"]) == pytest.approx(.01 * a)
    assert float(p["never_infected"]) == pytest.approx(math.exp(-.01 * a), abs=1e-15)
    assert float(p["remote_only"]) == pytest.approx(math.exp(-.01 * min(a,5)) * -math.expm1(-.01 * max(0,a-5)), abs=1e-15)
    assert float(p["recent"]) == pytest.approx(-math.expm1(-.01 * min(a,5)), abs=1e-15)
    validate_states(p)
    zero = evaluate_age(h, alpha=0, age=age, sex="male")
    assert [zero[k] for k in STATE_FIELDS] == ["1", "0", "0"]


def test_step_pressure_recency_window_and_missing_history():
    h = history(rows=[{"start": 1900, "end": 2000, "shape": 5}, {"start": 2000, "end": 2025, "shape": 1}])
    p = evaluate_age(h, alpha=".01", age=40, sex="male")
    assert float(p["cumulative_hazard"]) == pytest.approx(1.04)
    assert float(p["recent"]) == pytest.approx(-math.expm1(-.05))
    p2 = evaluate_age(h, alpha=".01", age=40, sex="male", recent_window=2)
    assert p["any_tbi"] == p2["any_tbi"] and Decimal(p["recent"]) > Decimal(p2["recent"])
    with pytest.raises(DataValidationError, match="Missing"):
        evaluate_age(h, alpha=".01", age=125, sex="male")
    with pytest.raises(DataValidationError):
        evaluate_age(h, alpha=".01", age=40, sex="male", time=2025)


@pytest.mark.parametrize("alpha", ["0", ".01", "1", "1000000", "0.000000000000000000000000000000000000000000000000001"])
def test_exact_uniform_single_year_age_averaging(alpha):
    p = evaluate_age(history(), alpha=alpha, age=0, age_upper_exclusive=1,
                     age_policy="uniform_within_single_year", sex="female")
    value = float(alpha)
    expected = 0 if value == 0 else 1 + math.expm1(-value) / value
    assert float(p["any_tbi"]) == pytest.approx(expected, abs=2e-15)
    validate_states(p)


def test_age_cell_averaging_splits_birth_calendar_boundary_exactly():
    h = history(rows=[{"start": 1900, "end": 2024, "shape": 2},
                      {"start": 2024, "end": 2025, "shape": 1}])
    p = evaluate_age(h, alpha=".01", age=0, age_upper_exclusive=1,
                     age_policy="uniform_within_single_year", time="2024.5", sex="female")
    survivor = -math.expm1(-.005) / .01 + math.exp(-.005) * -math.expm1(-.01) / .02
    assert float(p["never_infected"]) == pytest.approx(survivor, abs=2e-15)
    assert p["remote_only"] == "0"
    validate_states(p)


def test_age_cell_averaging_splits_fractional_recent_window_exactly():
    p = evaluate_age(history(), alpha=".01", age=2, age_upper_exclusive=3,
                     age_policy="uniform_within_single_year", recent_window="2.5", sex="male")
    survivor_before_recent = .5 + -math.expm1(-.005) / .01
    assert float(p["remote_only"]) == pytest.approx(math.exp(-.025)*(.5 - (-math.expm1(-.005)/.01)), abs=2e-15)
    assert float(p["never_infected"]) == pytest.approx(
        math.exp(-.02) * -math.expm1(-.01) / .01, abs=2e-15)
    validate_states(p)


def test_assay_operator_has_no_hidden_defaults_and_is_assay_specific():
    iga, tst = operator("IGRA", ".8", ".95"), operator("TST", ".8", ".95")
    assert assay_probability(".2", iga) == Decimal(".20")
    assert iga["profile_hash"] != tst["profile_hash"]
    for se, sp in (("1", "1"), (".4", ".6"), ("NaN", ".9"), ("1.1", ".9")):
        with pytest.raises(DataValidationError):
            operator("IGRA", se, sp)


@pytest.mark.parametrize("field,value", [("prevalence", "-.1"), ("prevalence", "1.1"),
    ("prevalence", "NaN"), ("numerator", -1), ("denominator", 0), ("sex", "unknown"),
    ("survey_year", 2025), ("strata", []), ("weight_source", None)])
def test_bad_observations_fail_without_rewriting_source(field, value):
    record = observation(history())
    record[field] = value
    original = deepcopy(record)
    with pytest.raises(DataValidationError):
        observations_profile([record], review=REVIEW)
    assert record == original


def test_observation_source_roundtrip_and_missing_versus_zero():
    h = history()
    original = observation(h, prevalence="0")
    original.update(numerator=0, denominator=100)
    artifact = observations_profile([original], review=REVIEW)
    assert artifact["rows"][0]["numerator"] == 0
    assert artifact["rows"][0]["receipt"]["source_record"] == original
    assert json.loads(canonical(artifact)) == artifact
    missing = observation(h)
    assert observations_profile([missing], review=REVIEW)["rows"][0]["numerator"] is None


@pytest.mark.parametrize("assay", [None, ("IGRA", ".85", ".96"), ("TST", ".7", ".9")])
def test_exact_one_scale_recovery_and_multiple_age_anchors(assay):
    h = history()
    assay = operator(*assay) if assay else None
    records = [observation(h, age=a, sex=s, alpha=".013", assay=assay) for a, s in ((5, "female"), (20, "male"), (65, "female"))]
    result = fit(h, records)
    assert result["fit_status"] == "accepted"
    assert abs(Decimal(result["alpha"]) - Decimal(".013")) < Decimal("1e-18")
    assert result["diagnostics"]["conditional_rank"] == 1 and not result["diagnostics"]["likelihood_claim"]
    assert all(abs(Decimal(r["residual"])) < Decimal("1e-15") for r in result["rows"])
    reverse = fit(h, list(reversed(records)))
    assert reverse["profile_hash"] == result["profile_hash"]


def test_incompatible_age_profile_is_reported_and_cannot_be_used(demography):
    h = history()
    result = fit(h, [observation(h, age=20, prevalence=".5"), observation(h, age=65, prevalence=".01")])
    assert result["fit_status"] == "incompatible_anchors"
    assert any(not r["compatible"] for r in result["rows"])
    with pytest.raises(DataValidationError, match="Incompatible"):
        baseline_tbi_profile(demography, result)


@pytest.mark.parametrize("target", ["1", "-.1", "1.1", "NaN"])
def test_impossible_targets_are_not_clipped(target):
    h = history()
    with pytest.raises(DataValidationError):
        fit(h, [observation(h, prevalence=target)])


def test_zero_exposure_scale_unidentifiability_and_assay_false_positive_floor():
    h = history()
    with pytest.raises(DataValidationError, match="unidentifiable"):
        fit(h, [observation(h, age=0, prevalence="0")])
    with pytest.raises(DataValidationError, match="incompatible"):
        fit(h, [observation(h, prevalence="0", assay=operator("IGRA", ".8", ".95"))])


def test_zero_scale_and_missing_anchor_require_explicit_provenance():
    h = history()
    assert fit(h, [observation(h, prevalence="0")])["alpha"] == "0"
    with pytest.raises(DataValidationError, match="Without an anchor"):
        calibrate_scale(h, policy())
    assert supplied(h, alpha=0)["diagnostics"]["identified_from_tbi_data"] is False
    with pytest.raises(DataValidationError):
        calibrate_scale(h, policy(), supplied_alpha=".01")


def test_supplied_scale_bounds_and_coherent_interval_sensitivity(demography):
    h = history()
    rec = observation(h, prevalence=".2", uncertainty={"kind": "interval", "lower": ".1", "upper": ".3", "method": "synthetic range"})
    obs = observations_profile([rec], review=REVIEW)
    result = calibrate_scale(h, policy([rec], uncertainty="coherent_observation_intervals"), observations=obs)
    u = result["alpha_uncertainty"]
    assert Decimal(u["lower"]) < Decimal(result["alpha"]) < Decimal(u["upper"])
    p = baseline_tbi_profile(demography, result)
    assert p["rows"][0]["parameter_endpoint_scenarios"]
    assert "not_complete_probability_intervals" in p["uncertainty_interpretation"]
    with pytest.raises(DataValidationError):
        supplied(h, interval={"lower": ".02", "upper": ".03", "provenance": {"source_id": "test"}})


def test_below_assay_floor_uncertainty_endpoint_is_unresolved_not_clipped():
    h = history()
    rec = observation(h, prevalence=".2", assay=operator("IGRA", ".8", ".9"),
        uncertainty={"kind": "interval", "lower": ".01", "upper": ".3", "method": "synthetic interval"})
    result = calibrate_scale(h, policy([rec], uncertainty="coherent_observation_intervals"),
                            observations=observations_profile([rec], review=REVIEW))
    assert result["alpha_uncertainty"]["lower"] is None
    assert result["alpha_uncertainty"]["lower_status"] == "outside_declared_feasible_range_not_clipped"


def test_population_sexes_do_not_infer_different_acquisition_and_tail_is_honest(demography):
    p = baseline_tbi_profile(demography, supplied())
    assert p["status"] == "partial_age_coverage"
    assert p["combined_totals"]["population"] == "3030"
    assert p["combined_totals"]["covered_population"] == "3000"
    assert p["combined_totals"]["unestimated_population"] == "30"
    assert p["combined_totals"]["any_tbi_prevalence_full"] is None
    assert {r["sex"] for r in p["rows"]} == {"female", "male"}
    females = [r for r in p["rows"] if r["sex"] == "female"]
    males = [r for r in p["rows"] if r["sex"] == "male"]
    assert all(a["any_tbi"] == b["any_tbi"] for a, b in zip(females, males))
    # Uniform ages [0,100), alpha=.01: mean survivor = (1-exp(-1))/1.
    expected = math.exp(-1)
    assert float(p["combined_totals"]["any_tbi_prevalence_covered"]) == pytest.approx(expected, abs=1e-14)
    with localcontext() as context:
        context.prec = 80
        assert sum(Decimal(r["population"]) for r in p["display_bands"]) == 3030
        total = p["combined_totals"]
        assert sum(Decimal(total[k + "_count"]) for k in STATE_FIELDS) == 3000


def test_explicit_sex_pressure_difference_is_preserved():
    female = [{"start": "2000", "end": "2025", "shape": "2"}]
    h = history(sex_rows={"female": female})
    pf = evaluate_age(h, alpha=".01", age=20, sex="female")
    pm = evaluate_age(h, alpha=".01", age=20, sex="male")
    assert float(pf["any_tbi"]) == pytest.approx(1 - math.exp(-.4))
    assert float(pm["any_tbi"]) == pytest.approx(1 - math.exp(-.2))


def test_reviewed_tail_restoration_exports_and_override_recalibration(demography):
    h = history()
    tails = {sex + ":100": {"review": REVIEW, "provenance": {"source_id": "synthetic"},
             "ages": [{"age": "101", "weight": "1"}]} for sex in ("female", "male")}
    p = baseline_tbi_profile(demography, supplied(h), tail_distributions=tails)
    assert p["status"] == "complete_age_coverage" and p["combined_totals"]["unestimated_population"] == "0"
    assert restore_tbi_profile(canonical(p)) == p
    original = deepcopy(p)
    rows = [observation(h, age=20, alpha=".02")]
    override = observations_profile(rows, review=REVIEW, mode="reviewed_profile_override", parent_profile=p)
    assert p == original and override["parent_profile_hash"] == p["profile_hash"]
    updated = calibrate_scale(h, policy(rows), observations=override)
    updated_profile = baseline_tbi_profile(demography, updated, tail_distributions=tails)
    assert abs(Decimal(updated["alpha"]) - Decimal(".02")) < Decimal("1e-18")
    assert p["profile_hash"] != updated_profile["profile_hash"]
    files = preview_files(updated_profile)
    for name in ("profile.json", "diagnostics.json"):
        assert json.loads(files[name])["exclusions"] == list(EXCLUSIONS)
    for name in ("age_sex_tbi.csv", "display_bands.csv"):
        records = list(csv.DictReader(io.StringIO(files[name].decode())))
        assert records and all(r["model_exclusions"] == ";".join(EXCLUSIONS) for r in records)
    bad = deepcopy(p); bad["alpha"] = ".9"
    with pytest.raises(DataValidationError):
        restore_tbi_profile(canonical(bad))


def test_scenarios_and_material_hashes(demography):
    first = baseline_tbi_profile(demography, supplied(alpha=".01"))
    second = baseline_tbi_profile(demography, supplied(alpha=".02"))
    assert first["profile_hash"] != second["profile_hash"]
    assert profile_hash({**first, "display": {"name": "label"}}) == first["profile_hash"]
    scenarios = uncertainty_scenario_set({"a": first, "b": second}, review=REVIEW, weights={"a": ".5", "b": ".5"})
    assert scenarios["rng"] is None and len(scenarios["rows"]) == 2
    with pytest.raises(DataValidationError):
        uncertainty_scenario_set({"a": first, "b": second}, review=REVIEW, weights={"a": ".6", "b": ".5"})


def test_new_pathway_is_offline_and_unwired(monkeypatch, demography):
    def forbidden(*args, **kwargs):
        raise AssertionError("Network prohibited")
    monkeypatch.setattr(socket, "socket", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)
    p = baseline_tbi_profile(demography, supplied())
    assert p["profile_hash"] == profile_hash(p)
    for folder in ("app", "ui", "general_ui", "explicit_tbi_ui", "engine/apy"):
        for path in (REPO_ROOT / folder).rglob("*.py"):
            source = path.read_text(encoding="utf-8-sig")
            assert all(name not in source for name in ("catalytic_calibration", "tbi_baseline", "pressure_shape"))


def test_both_only_remains_both_and_never_splits_or_double_counts(demography):
    rows = [{**r, "sex": "both", "population": "30"} for r in demography["rows"] if r["sex"] == "female"]
    both = demo_for_rows(rows)
    p = baseline_tbi_profile(both, supplied())
    assert {r["sex"] for r in p["rows"]} == {"both"}
    assert p["combined_totals"]["population"] == "3030"
    bad = deepcopy(demography)
    bad["rows"].extend(rows)
    with pytest.raises(DataValidationError):
        baseline_tbi_profile(seal(bad), supplied())


def test_sex_specific_pressure_cannot_invent_aggregate_sex_weights():
    h = history(sex_rows={"female": [{"start": 2000, "end": 2025, "shape": 2}]})
    with pytest.raises(DataValidationError, match="aggregate"):
        evaluate_age(h, alpha=".01", age=20, sex="both")


def test_observation_weight_order_and_source_receipt_integrity():
    h = history()
    row = observation(h, sex="both")
    row["strata"] = [{**row["strata"][0], "sex": sex, "weight": str(weight)}
                     for sex, weight in (("female", 10), ("male", 20))]
    a = observations_profile([row], review=REVIEW)
    row["strata"].reverse()
    b = observations_profile([row], review=REVIEW)
    assert a["profile_hash"] == b["profile_hash"]
    bad = deepcopy(a)
    bad["rows"][0]["receipt"]["source_record"]["prevalence"] = ".9"
    with pytest.raises(DataValidationError, match="receipt"):
        calibrate_scale(h, policy([row]), observations=bad)


def test_exact_anchors_resealed_calibration_and_version_failures(demography):
    h = history()
    rows = [observation(h)]
    result = fit(h, rows)
    result["alpha"] = ".5"
    with pytest.raises(DataValidationError, match="reproduce"):
        baseline_tbi_profile(demography, seal(result))
    old = deepcopy(h); old["contract_version"] = "1.0.0-design.1"
    with pytest.raises(DataValidationError, match="migration"):
        evaluate_age(seal(old), alpha=".01", age=20, sex="female")


@pytest.mark.parametrize("case", ["negative_alpha", "finite_bound", "iterations", "weights", "geography", "unreviewed", "recency_input"])
def test_calibration_failures_are_explicit(case):
    h = history()
    rows = [observation(h)]
    if case == "negative_alpha":
        with pytest.raises(DataValidationError):
            supplied(h, alpha="-1")
    elif case == "finite_bound":
        with pytest.raises(DataValidationError, match="finite scale"):
            calibrate_scale(h, policy(rows, alpha_upper=".0001", alpha_tolerance=".00000001"),
                            observations=observations_profile(rows, review=REVIEW))
    elif case == "iterations":
        with pytest.raises(DataValidationError, match="convergence"):
            calibrate_scale(h, policy(rows, max_iterations=1), observations=observations_profile(rows, review=REVIEW))
    elif case == "weights":
        with pytest.raises(DataValidationError, match="cover exactly"):
            calibrate_scale(h, policy(), observations=observations_profile(rows, review=REVIEW))
    elif case == "geography":
        rows[0]["geography_id"] = "custom:another"
        with pytest.raises(DataValidationError, match="geographic"):
            fit(h, rows)
    elif case == "unreviewed":
        with pytest.raises(DataValidationError, match="reviewed"):
            observations_profile(rows, review={**REVIEW, "status": "not_reviewed"})
    else:
        rows[0]["recent_prevalence"] = ".1"
        with pytest.raises(DataValidationError, match="Recency"):
            observations_profile(rows, review=REVIEW)


def test_cell_averaging_cannot_silently_split_coarse_or_open_intervals():
    h = history()
    with pytest.raises(DataValidationError):
        evaluate_age(h, alpha=".01", age=0, age_upper_exclusive=5,
                     age_policy="uniform_within_single_year", sex="female")
    with pytest.raises(DataValidationError, match="open"):
        evaluate_age(h, alpha=".01", age=100, age_upper_exclusive=None,
                     age_policy="uniform_within_single_year", sex="female")
    with pytest.raises(DataValidationError, match="parameters"):
        apply_backcast(shape(), start_time=1900, parameters={"slope_per_year": "-.02"})


def test_executable_schema_and_cross_field_validation(demography):
    from jsonschema import Draft202012Validator, ValidationError
    schema = json.loads((REPO_ROOT / "schemas/country_tb/v1/catalytic_calibration.schema.json").read_text())
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    q = shape()
    h = apply_backcast(q, start_time=1900)
    records = [observation(h)]
    obs = observations_profile(records, review=REVIEW)
    config = policy(records)
    calibration = calibrate_scale(h, config, observations=obs)
    p = baseline_tbi_profile(demography, calibration)
    for value in (q, h, operator(), obs, config, calibration, p,
                  evaluate_age(h, alpha=".01", age=20, sex="female")):
        validator.validate(value)
    files = preview_files(p)
    validator.validate(json.loads(files["diagnostics.json"]))
    invalid = deepcopy(p); invalid["rows"][0]["any_tbi"] = "1.1"
    with pytest.raises(ValidationError):
        validator.validate(invalid)
    invalid = deepcopy(p); invalid["rows"][0]["never_infected"] = "0.9"
    with pytest.raises(DataValidationError):
        validate_states(invalid["rows"][0])
    invalid = deepcopy(p); invalid["exclusions"] = []
    with pytest.raises(ValidationError):
        validator.validate(invalid)


def test_protected_worktrees_refs_and_all_frozen_inputs_are_preserved():
    receipt = REPO_ROOT / "abm/output/country_tb_m2c_audit/starting_state.json"
    if not receipt.is_file():
        pytest.skip("Live session protection receipt unavailable on this checkout.")
    state = json.loads(receipt.read_text())
    assert state["starting_head"] == "a1690c65e579d7ee375a3852436ac56aec556705"
    for filename, sha in state["frozen_files"].items():
        assert digest((REPO_ROOT / filename).read_bytes()) == sha, filename
    for name, expected in state["protected"].items():
        folder = REPO_ROOT.parent / name
        def git(*args):
            return subprocess.check_output(["git", "-C", str(folder), *args], stderr=subprocess.DEVNULL)
        assert git("rev-parse", "HEAD").decode().strip() == expected["head"]
        assert git("branch", "--show-current").decode().strip() == expected["branch"]
        assert git("status", "--porcelain=v1").decode().splitlines() == expected["status"]
        assert digest(git("diff", "--binary", "HEAD")) == expected["tracked_diff_sha256"]
        assert digest(git("diff", "--binary", "--cached")) == expected["index_diff_sha256"]
        names = set(git("diff", "--name-only", "HEAD", "-z").decode().strip("\0").split("\0")) - {""}
        names |= set(git("ls-files", "--others", "--exclude-standard", "-z").decode().strip("\0").split("\0")) - {""}
        assert {n: digest((folder/n).read_bytes()) if (folder/n).is_file() else "absent"
                for n in sorted(names)} == expected["changed_files"]
    refs = dict(line.split(" ",1) for line in subprocess.check_output(
        ["git","for-each-ref","--format=%(refname) %(objectname)"]).decode().splitlines())
    allowed = {"refs/heads/feature/country-subnational-tb-model-v1", "refs/remotes/origin/feature/country-subnational-tb-model-v1"}
    assert {k:v for k,v in refs.items() if k not in allowed} == {k:v for k,v in state["refs"].items() if k not in allowed}


def test_developer_preview_recovery_population_and_all_export_schemas():
    from jsonschema import Draft202012Validator
    from scripts.preview_country_tb_calibration import synthetic_example, shape_files
    q, h, obs, config, calibration, profile = synthetic_example()
    assert abs(Decimal(calibration["alpha"]) - Decimal(".013")) < Decimal("1e-18")
    assert profile["combined_totals"]["population"] == "190018"
    assert profile["combined_totals"]["covered_population"] == "190000"
    assert profile["combined_totals"]["unestimated_population"] == "18"
    schema = json.loads((REPO_ROOT / "schemas/country_tb/v1/catalytic_calibration.schema.json").read_text())
    validator = Draft202012Validator(schema)
    for value in (q, h, obs, config, calibration, profile):
        validator.validate(value)
    assert "shape_diagnostics.csv" in shape_files(q)


@pytest.mark.parametrize("age", [0, 2, 5, 40])
def test_reinfection_states_and_event_count(age):
    p = _evaluate_age(history(), alpha='.02', age=age, sex='female')
    hr, hn = .02*max(age-5,0), .02*min(age,5)
    assert float(p['recent']) == pytest.approx(-math.expm1(-hn))
    assert float(p['remote_only']) == pytest.approx(math.exp(-hn)*-math.expm1(-hr))
    assert float(p['prior_remote_plus_recent']) == pytest.approx(-math.expm1(-hr)*-math.expm1(-hn))
    assert p['expected_lifetime_infection_events'] == p['cumulative_hazard']
    assert Decimal(p['prior_remote_plus_recent']) <= Decimal(p['recent'])
    validate_states(p)
    if age <= 5:
        assert p['remote_only'] == p['prior_remote_plus_recent'] == '0'


def test_last_event_analytic_cdf_inverse_mean_and_remote_condition():
    from engine.country_tb.most_recent_event import interval_distribution, conditional_cdf, inverse_cdf, summarize
    d = interval_distribution(history(), sex='female', alpha='.02', start=1984, end=2024)
    assert float(conditional_cdf(d, 2000)) == pytest.approx(math.exp(-.48)*-math.expm1(-.32)/-math.expm1(-.8))
    result = summarize(d)
    assert float(result['expected_time_since_most_recent_infection']) == pytest.approx(50-40/math.expm1(.8))
    for q in ('0','.000001','.05','.5','.95','.999999','1'):
        t = inverse_cdf(d,q)
        assert float(conditional_cdf(d,t)) == pytest.approx(float(q),abs=1e-14)
    p = _evaluate_age(history(),alpha='.02',age=40,sex='female')
    assert float(p['remote_only_time_summary']['expected_time_since_most_recent_infection']) == pytest.approx(5+50-35/math.expm1(.7))
    assert float(p['recent_time_summary']['expected_time_since_most_recent_infection']) < 5


@pytest.mark.parametrize('alpha', ['0', '0.00000000000000000000000000000000000000000000000001', '1000000'])
def test_last_event_extreme_hazards_and_rare_remote_group(alpha):
    from engine.country_tb.most_recent_event import interval_distribution, conditional_cdf, inverse_cdf, summarize
    d = interval_distribution(history(),sex='female',alpha=alpha,start=1984,end=2024)
    s = summarize(d)
    if Decimal(alpha) == 0:
        assert s['status'] == 'undefined_zero_probability'
        assert inverse_cdf(d,'.5') is None and conditional_cdf(d,2000) is None
    else:
        assert float(conditional_cdf(d,inverse_cdf(d,'.5'))) == pytest.approx(.5,abs=1e-14)
        expected = 20 if Decimal(alpha)<Decimal('1e-40') else .000001
        assert float(s['expected_time_since_most_recent_infection']) == pytest.approx(expected,rel=1e-12)
    if Decimal(alpha)>1:
        p = _evaluate_age(history(),alpha=alpha,age=40,sex='female')
        assert p['remote_only_time_summary']['status']=='defined'
        assert float(p['remote_only_time_summary']['expected_time_since_most_recent_infection']) == pytest.approx(5.000001)


@pytest.mark.parametrize('shapes', [(1,3),(3,1),(1,1)])
def test_rising_declining_step_last_event_distributions(shapes):
    from engine.country_tb.most_recent_event import interval_distribution,conditional_cdf,inverse_cdf
    h=history(rows=[{'start':1900,'end':2010,'shape':shapes[0]}, {'start':2010,'end':2025,'shape':shapes[1]}])
    d=interval_distribution(h,sex='male',alpha='.03',start=1984,end=2024)
    # Independent fine-grid midpoint integration of the analytic density.
    curve=h['curves']['common']; total=float(evaluate_age(h,alpha='.03',age=40,sex='male')['cumulative_hazard'])
    numerator=0
    for i in range(16000):
        u=1984+(i+.5)*.001
        future=.03*((2010-u)*shapes[0]/shapes[1]+14) if u<2010 else .03*(2024-u)
        rate=.03*shapes[0]/shapes[1]
        numerator+=rate*math.exp(-future)*.001
    assert float(conditional_cdf(d,2000))==pytest.approx(numerator/-math.expm1(-total),rel=1e-8)
    for q in ('.05','.5','.95'):
        assert float(conditional_cdf(d,inverse_cdf(d,q)))==pytest.approx(float(q),abs=1e-13)


def test_uniform_age_event_time_density_agrees_with_probability():
    p=_evaluate_age(history(),alpha='.03',age=0,age_upper_exclusive=1,
        age_policy='uniform_within_single_year',sex='female')
    summary=p['most_recent_infection_time']
    assert float(summary['group_probability'])==pytest.approx(float(p['any_tbi']),abs=1e-15)
    numerator=sum(.03*math.exp(-.03*((i+.5)/10000))*(1-(i+.5)/10000)*(i+.5)/10000/10000 for i in range(10000))
    assert float(summary['expected_time_since_most_recent_infection'])==pytest.approx(numerator/float(p['any_tbi']),abs=2e-9)


def test_explicit_rng_last_event_and_poisson_reconciliation():
    from engine.country_tb.most_recent_event import interval_distribution,sample_most_recent_infection_time
    h=history();d=interval_distribution(h,sex='female',alpha='.03',start=1984,end=2024)
    a,b=random.Random(712),random.Random(712)
    assert [sample_most_recent_infection_time(d,rng=a) for _ in range(10)] == [sample_most_recent_infection_time(d,rng=b) for _ in range(10)]
    with pytest.raises(DataValidationError):sample_most_recent_infection_time(d,rng=None)
    rng=random.Random(1891);counts={'never_infected':0,'recent':0,'remote_only':0,'prior_remote_plus_recent':0};events=0;since=[]
    n=25000
    for _ in range(n):
        u=0;times=[]
        while True:
            u+=rng.expovariate(.03)
            if u>40:break
            times.append(u)
        events+=len(times)
        remote=any(t<35 for t in times);recent=any(t>=35 for t in times)
        counts['recent' if recent else 'remote_only' if remote else 'never_infected']+=1
        counts['prior_remote_plus_recent']+=remote and recent
        if times:since.append(40-times[-1])
    p=_evaluate_age(h,alpha='.03',age=40,sex='female')
    for key,value in counts.items():assert value/n==pytest.approx(float(p[key]),abs=.012)
    assert events/n==pytest.approx(float(p['cumulative_hazard']),abs=.025)
    assert sum(since)/len(since)==pytest.approx(float(p['expected_time_since_most_recent_infection']),abs=.35)


def test_binomial_likelihood_counts_scale_recovery_and_no_pseudocounts():
    h=history();alpha=str(-math.log(.8)/20)
    rows=[observation(h,age=20,sex=s,prevalence='.2',uncertainty={'kind':'unspecified'}) for s in ('female','male')]
    for row in rows:row.update(numerator=200,denominator=1000)
    p=calibrate_scale(h,policy(rows,method='binomial'),observations=observations_profile(rows,review=REVIEW))
    assert float(p['alpha'])==pytest.approx(float(alpha),abs=1e-14)
    assert p['diagnostics']['likelihood_claim'] and p['diagnostics']['converged']
    assert sum(float(r['objective_contribution']) for r in p['rows'])==pytest.approx(float(p['diagnostics']['objective']))
    rows[0]['numerator']=rows[0]['denominator']=None
    with pytest.raises(DataValidationError,match='actual counts'):
        calibrate_scale(h,policy(rows,method='binomial'),observations=observations_profile(rows,review=REVIEW))


def test_overall_aggregate_anchor_and_near_one_recovery():
    h=history();row=observation(h,sex='both')
    row['age_upper_exclusive']='81';row['strata']=[{'age_lower':str(a),'age_upper_exclusive':None,
        'age_policy':'exact_age','sex':s,'weight':str(w)} for a,s,w in ((5,'female',3),(40,'male',7),(80,'female',2))]
    row['age_lower']='5'
    with localcontext() as context:
        context.prec=80
        row['prevalence']=str(sum(Decimal(r['weight'])*Decimal(evaluate_age(h,alpha='.017',age=r['age_lower'],sex=r['sex'])['any_tbi']) for r in row['strata'])/12)
    assert float(fit(h,[row])['alpha'])==pytest.approx(.017,abs=1e-15)
    near=observation(h,age=20,prevalence='.999999999')
    result=calibrate_scale(h,policy([near],alpha_upper='2'),observations=observations_profile([near],review=REVIEW))
    assert float(result['alpha'])==pytest.approx(-math.log(1e-9)/20,abs=1e-9)


def test_total_only_override_leaves_recency_unestimated(demography):
    from engine.country_tb.tbi_observations import reviewed_tbi_override
    parent=baseline_tbi_profile(demography,supplied());before=deepcopy(parent)
    records=[observation(history(),prevalence='.3')]
    result=reviewed_tbi_override(parent,records,review=REVIEW)
    assert parent==before and result['profile_hash']!=parent['profile_hash']
    assert result['rows'][0]['any_tbi']=='0.3'
    assert all(result['rows'][0][key] is None for key in ('recent','remote_only','prior_remote_plus_recent','most_recent_infection_time'))
    with pytest.raises(DataValidationError,match='Overlapping'):
        reviewed_tbi_override(parent,[*records,{**records[0],'observation_id':'another'}],review=REVIEW)


def test_missing_year_diagnostic_and_observation_contract_fields():
    source=incidence();source['rows'][5]['central']='0';source['rows'][6]['central']=None
    with pytest.raises(DataValidationError,match=r'2005.*2006'):
        shape_from_incidence(seal(source),review=REVIEW)
    row=observation(history())
    for key in ('geography_level','observation_policy'):
        bad=deepcopy(row);bad.pop(key)
        with pytest.raises(DataValidationError):observations_profile([bad],review=REVIEW)
    direct=operator('direct_latent_tbi')
    assert assay_probability('.3',direct)==Decimal('.3')


def test_inverse_cdf_sampling_reconciles_with_conditional_law():
    from engine.country_tb.most_recent_event import interval_distribution,sample_most_recent_infection_time,summarize,conditional_cdf
    d=interval_distribution(history(),sex='female',alpha='.03',start=1984,end=2024)
    rng=random.Random(4127)
    times=[sample_most_recent_infection_time(d,rng=rng) for _ in range(5000)]
    expected=float(summarize(d)['expected_time_since_most_recent_infection'])
    assert sum(2024-float(t) for t in times)/len(times)==pytest.approx(expected,abs=.5)
    assert sum(t<=2000 for t in times)/len(times)==pytest.approx(float(conditional_cdf(d,2000)),abs=.025)


@pytest.mark.parametrize('change',[{'uncertainty':None},{'uncertainty':{'kind':'interval'}},{'strata':[{}]},{'strata':None}])
def test_malformed_observation_structures_fail_explicitly(change):
    row=observation(history());row.update(change)
    with pytest.raises(DataValidationError):observations_profile([row],review=REVIEW)
