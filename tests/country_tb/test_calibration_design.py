"""Milestone 2B independent evidence checks and synthetic mathematics only."""
from copy import deepcopy
from decimal import Decimal, localcontext
import ast
import csv
import io
import json
import math
from pathlib import Path
import random
import socket
import subprocess

from jsonschema import Draft202012Validator, ValidationError
import pytest

from engine.country_tb.catalytic_diagnostic import (
    first_infection_profile, normalized_history, sample_states, synthetic_scale_solution,
)
from engine.country_tb.data_common import DataValidationError, REPO_ROOT, canonical, digest, profile_hash
from engine.country_tb.data_profile import baseline_data_profile, load_candidates, restore_baseline_data
from engine.country_tb.incidence import WHO_2025, country_index, load_who_snapshot
from scripts.recheck_country_tb_frozen import independent_scientific_hash, legacy_check, recheck
from scripts.preview_country_tb_catalytic import cumulative, examples, exponential_integral

ROOT = REPO_ROOT
BASE = "a4fceee5af7c2871943061ad9635084fd18cec89"
FIELDS = ("uninfected", "recently_first_infected", "remote_only_infected")


def history(rate="0.01", start="1900", end="2025"):
    return [{"start": start, "end": end, "hazard_per_person_year": rate}]


def profile(curve=None, **overrides):
    options = dict(age=40, baseline_time=2024, recent_window=2, sex="female",
                   geography_id="synthetic:G", scenario_id="test",
                   provenance={"source_id": "synthetic", "sha256": "0" * 64})
    options.update(overrides)
    return first_infection_profile(curve if curve is not None else history(), **options)


@pytest.fixture(scope="module")
def independent():
    return recheck()


@pytest.fixture(scope="module")
def sources():
    return load_who_snapshot(WHO_2025), load_candidates()


def test_independent_who_and_native_artifact_reconstruction(independent):
    who = independent["who"]
    assert (who["rows"], who["countries_areas"], who["with_estimates"]) == (5347, 217, 216)
    assert who["years"] == list(range(2000, 2025)) and who["report_year"] == 2025
    assert who["transformed_file_sha256"] == "6997e3b087d26788b21a47d09968eca87b5b03c3cc50644ae4dcdcd860acfd07"
    assert who["manifest_analytical_hash"] == "32f228fc4338a9a61d56fa4770799c2274a3b9696f814d2304b3e2a3ebdee765"
    assert who["missing_triples"] == {"incidence_per_100k": 25, "incident_cases": 25,
                                     "hiv_positive_incidence_per_100k": 895}
    assert who["PRK_all_missing"]
    assert independent["demography"][0]["independently_reconstructed_rows"] == 303
    assert independent["demography"][1]["independently_reconstructed_rows"] == 502677


def test_independent_complete_duplicate_decomposition(independent):
    result = independent["legacy"]
    assert (result["source_rows"], result["unique_keys"], result["duplicate_groups"],
            result["excess_rows"]) == (21917, 20705, 1111, 1212)
    assert result["groups_by_multiplicity"] == {2: 1010, 3: 101}
    assert result["maximum_difference"] == "0" and result["all_repeated_stored_fields_identical"]
    assert result["per_country_excess"]["VCT"] == 202
    assert sum(v for k, v in result["per_country_excess"].items() if k != "VCT") == 1010


def test_independent_dedup_order_idempotence_and_conflicts():
    raw = (ROOT / "data/population_age_latest.csv").read_bytes()
    _, unique = legacy_check(raw)
    assert legacy_check(unique)[1] == unique
    lines = raw.splitlines()
    rng = random.Random(821)
    shuffled = lines[1:]
    rng.shuffle(shuffled)
    assert legacy_check(b"\n".join([lines[0], *shuffled]) + b"\n")[1] == unique
    with pytest.raises(ValueError, match="Conflicting"):
        legacy_check(b"location,iso_code,year,age,population\nX,TST,2021,0,1\nX,TST,2021,0,2\n")


@pytest.mark.parametrize("country,total,status", [
    ("AUS", "27204809", "estimate"), ("VNM", "100987682", "projection"),
    ("PRK", "26498835", "projection"),
])
def test_native_source_vs_loader_sex_age_year_and_restore(independent, sources, monkeypatch, country, total, status):
    who, candidates = sources
    p = baseline_data_profile(country, who, candidates)
    d = p["demography"]
    selected = independent["demography"][0 if country == "AUS" else 1]
    native = selected["examples"]["who:" + country]
    assert d["combined_population"] == native["derived_combined"] == total
    assert {r["estimate_status"] for r in d["rows"]} == {status}
    assert {r["year"] for r in d["rows"]} == {2024}
    assert {r["sex"] for r in d["rows"]} == {"female", "male"}
    assert len(d["rows"]) == 202 and len(d["combined_rows"]) == 101
    assert sum(Decimal(r["population"]) for r in d["combined_rows"]) == Decimal(total)
    for sex in ("female", "male"):
        assert d["sex_totals"][sex] == native["sex_totals"][sex]
        assert sum(Decimal(r["population"]) for r in d["display_bands"] if r["sex"] == sex) == Decimal(d["sex_totals"][sex])
    # Restoration binds original IDs, never consults/reranks a "latest" catalogue.
    def forbidden_catalogue(*args, **kwargs):
        raise AssertionError("Restoration must not query current source candidates")
    monkeypatch.setattr("engine.country_tb.data_profile.load_candidates", forbidden_catalogue)
    restored = restore_baseline_data(canonical(p))
    assert restored == p
    assert d["provenance"]["snapshot_id"] == selected["snapshot_id"]
    assert d["provenance"]["manifest_hash"] == selected["manifest_hash"]
    for component in (p, p["incidence"], p["demography"]):
        assert component["profile_hash"] == independent_scientific_hash(component) == profile_hash(component)


def test_frozen_scientific_profile_hashes_still_match_m2a(sources):
    who, candidates = sources
    expected = {
        "AUS": "11c7e058cbe54db4af1d302c0a95044807bacdabf36a2ea5003449b7c8fc8a29",
        "VNM": "f08621c9e3574eeddc01e83629578962bcd7fa4a60bec4eaf7b5d2127ba2f08a",
        "PRK": "78fd3c8c8ce0c98faf2842fefcd0d6e2d4c5602312d197def6653d0f0d005954",
    }
    for code, wanted in expected.items():
        assert baseline_data_profile(code, who, candidates)["profile_hash"] == wanted


def test_territories_missing_values_and_offline_loading(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Network forbidden in ordinary load")
    monkeypatch.setattr(socket, "socket", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)
    who = load_who_snapshot(WHO_2025)
    index = {r["country_id"]: r for r in country_index(who)}
    assert len(index) == 217
    assert all(code in index for code in ("ANT", "SCG", "CUW", "SXM", "PRK"))
    assert index["PRK"]["completeness"] == "unavailable"
    assert index["ANT"]["record_years"] == list(range(2000, 2010))
    assert index["SCG"]["record_years"] == list(range(2000, 2005))
    baseline_data_profile("PRK", who, load_candidates())


@pytest.mark.parametrize("age", [0, "0.001", 1, 4, 14, 40, 85, 124])
def test_zero_and_constant_pressure_closed_form(age):
    zero = profile(history("0"), age=age)
    assert [zero[k] for k in FIELDS] == ["1", "0", "0"]
    p = profile(age=age)
    a = float(age)
    assert float(p["cumulative_hazard"]) == pytest.approx(.01 * a)
    assert float(p["uninfected"]) == pytest.approx(math.exp(-.01 * a), abs=1e-15)
    assert float(p["recently_first_infected"]) == pytest.approx(
        math.exp(-.01 * max(0, a - 2)) * -math.expm1(-.01 * min(a, 2)), abs=1e-15)
    with localcontext() as ctx:
        ctx.prec = 80
        assert sum(Decimal(p[k]) for k in FIELDS) == 1
        assert Decimal(p["total_ever_infected"]) == Decimal(p[FIELDS[1]]) + Decimal(p[FIELDS[2]])


@pytest.mark.parametrize("rates", [("0.02", "0.01"), ("0.01", "0.02"), ("0", "0.03")])
def test_declining_rising_step_curves(rates):
    curve = [{"start": "1900", "end": "2000", "hazard_per_person_year": rates[0]},
             {"start": "2000", "end": "2025", "hazard_per_person_year": rates[1]}]
    p = profile(curve, age=40)
    expected = 16 * float(rates[0]) + 24 * float(rates[1])
    assert float(p["cumulative_hazard"]) == pytest.approx(expected)
    assert float(p["total_ever_infected"]) == pytest.approx(-math.expm1(-expected))
    assert profile(list(reversed(curve)), age=40) == p


def test_first_not_last_recency_is_explicit():
    p = profile(age=40)
    last_event_recent = -math.expm1(-.02)
    assert float(p["recently_first_infected"]) < last_event_recent
    assert "reinfection" in p["exclusions"]
    assert p["status"] == "synthetic_diagnostic_only"
    assert p["profile_hash"] == profile_hash(p)


@pytest.mark.parametrize("change", [
    {"age": -1}, {"age": True}, {"age": "NaN"}, {"recent_window": 0},
    {"sex": ""}, {"provenance": {}}, {"age": 125},
])
def test_invalid_or_missing_history_fails(change):
    with pytest.raises(DataValidationError):
        profile(**change)


@pytest.mark.parametrize("curve", [
    [], history("-1"), history("NaN"), history("Infinity"),
    [{"start": 2000, "end": 2000, "hazard_per_person_year": 1}],
    [{"start": 1900, "end": 2000, "hazard_per_person_year": 1},
     {"start": 2001, "end": 2025, "hazard_per_person_year": 1}],
    [{"start": 1900, "end": 2001, "hazard_per_person_year": 1},
     {"start": 2000, "end": 2025, "hazard_per_person_year": 1}],
])
def test_invalid_history_fail_without_clipping(curve):
    with pytest.raises(DataValidationError):
        profile(curve)


def test_extreme_valid_pressure_and_hash_materiality():
    p = profile(history("100"), age=100)
    assert p["uninfected"] == "0"  # Declared 40-place numerical rounding, not input clipping.
    assert p["total_ever_infected"] == "1"
    ordinary = profile()
    for change in ({"sex": "male"}, {"age": 41}, {"scenario_id": "alternative"},
                   {"recent_window": 5}, {"geography_id": "synthetic:H"},
                   {"provenance": {"source_id": "new", "sha256": "1" * 64}}):
        assert profile(**change)["profile_hash"] != ordinary["profile_hash"]
    assert profile(history(".02"))["profile_hash"] != ordinary["profile_hash"]
    labelled = {**ordinary, "display": {"name": "presentation only"}}
    assert profile_hash(labelled) == ordinary["profile_hash"]
    assert profile(history(".010")) == ordinary


def test_exact_synthetic_anchor_scale_and_age_recovery():
    true_scale = Decimal(".013")
    with localcontext() as ctx:
        ctx.prec = 80
        ages, weights = [5, 20, 65], [10, 30, 60]
        target = sum(Decimal(w) * (1 - (-true_scale * a).exp()) for a, w in zip(ages, weights)) / 100
    fit = synthetic_scale_solution(ages, weights, target, upper_scale="1")
    scale = Decimal(fit["scale"])
    assert abs(scale - true_scale) < Decimal("1e-45")
    assert not fit["empirical_fit"]
    for age in ages:
        actual = profile(history(str(scale)), age=age)
        expected = profile(history(str(true_scale)), age=age)
        assert abs(Decimal(actual["total_ever_infected"]) - Decimal(expected["total_ever_infected"])) <= Decimal("1e-39")
    assert synthetic_scale_solution([10], [1], 0, upper_scale=1)["scale"] == "0"


@pytest.mark.parametrize("exposures,weights,target,bound", [
    ([10], [1], 1, 1), ([10], [1], -1, 1), ([0], [1], 0, 1),
    ([10], [0], ".2", 1), ([10], [1], ".99", ".001"),
    ([10], [1], ".2", 0), ([1], [1, 2], ".2", 1), ([], [], ".2", 1),
    ([10], [1], "NaN", 1),
])
def test_infeasible_or_unidentifiable_synthetic_targets_fail(exposures, weights, target, bound):
    with pytest.raises(DataValidationError):
        synthetic_scale_solution(exposures, weights, target, upper_scale=bound)


def test_explicit_rng_sampling_and_deterministic_reconciliation():
    p, n = profile(age=65), 100000
    counts = sample_states(p, n, rng=random.Random(20241028))
    assert counts == sample_states(p, n, rng=random.Random(20241028))
    assert sum(counts.values()) == n
    for field in FIELDS:
        probability = float(p[field])
        tolerance = 6 * math.sqrt(probability * (1 - probability) / n) + 1 / n
        assert abs(counts[field] / n - probability) <= tolerance
    class InvalidRng:
        def random(self):
            return 1
    with pytest.raises(DataValidationError):
        sample_states(p, 1, rng=InvalidRng())
    with pytest.raises(TypeError):
        sample_states(p, 1)


def test_male_female_weighted_counts_reconcile():
    female, male = profile(sex="female"), profile(sex="male", curve=history(".02"))
    with localcontext() as ctx:
        ctx.prec = 80
        combined = [Decimal(female[k]) * 60 + Decimal(male[k]) * 40 for k in FIELDS]
        assert sum(combined) == 100
    assert female["sex"] == "female" and male["sex"] == "male"


def test_backcast_synthetic_example_independent_integrals():
    rows = examples()
    assert len(rows) == 40
    # An independent fine midpoint quadrature checks the analytic rising/declining integral.
    for slope in (-.02, 0, .02):
        dt = .001
        numeric = sum(.01 * math.exp(-slope * (1960 + (i + .5) * dt - 2000)) * dt for i in range(40000))
        assert exponential_integral(1960, 2000, slope) == pytest.approx(numeric, rel=1e-9)
    for age in (4, 24):
        assert len({r["total_ever_infected"] for r in rows if r["age"] == age}) == 1
    older = {r["policy"]: r["total_ever_infected"] for r in rows if r["age"] == 85}
    assert older["bounded_early_trend"] > older["constant_earliest"]
    assert cumulative(65, "anchor_conditional_slope") == pytest.approx(.6974260892097832, rel=1e-12)


def test_proposed_schema_examples_hashes_and_failures():
    schema = json.loads((ROOT / "schemas/country_tb/v1/calibration_design.schema.json").read_text())
    examples = json.loads((ROOT / "schemas/country_tb/v1/calibration_design_examples.json").read_text())
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    assert len(examples) == 6
    for example in examples:
        validator.validate(example)
        assert example["profile_hash"] == profile_hash(example)
        assert json.loads(canonical(example)) == example
        changed = deepcopy(example)
        changed["provenance"][0]["sha256"] = "1" * 64
        assert profile_hash(changed) != example["profile_hash"]
        changed["contract_version"] = "1.0.0"
        with pytest.raises(ValidationError):
            validator.validate(changed)
    obs = next(e for e in examples if e["schema_id"] == "country_tb.tbi_observations.v1")
    bad = deepcopy(obs)
    bad["data"]["rows"][0]["prevalence"] = "1.1"
    with pytest.raises(ValidationError):
        validator.validate(bad)
    bad["data"]["rows"][0]["prevalence"] = "NaN"
    with pytest.raises(ValidationError):
        validator.validate(bad)
    pressure = next(e for e in examples if e["schema_id"] == "country_tb.infection_pressure_history.v1")
    bad = deepcopy(pressure)
    bad["data"]["rows"][0]["hazard_per_person_year"] = "-1"
    with pytest.raises(ValidationError):
        validator.validate(bad)


def test_protected_starting_receipt_and_frozen_files():
    receipt = ROOT / "abm/output/country_tb_m2b_audit/starting_state.json"
    if not receipt.is_file():
        pytest.skip("Live session fingerprint receipt unavailable on this checkout.")
    state = json.loads(receipt.read_text())
    assert state["starting_head"] == BASE
    for filename, expected in state["frozen_files"].items():
        assert digest((ROOT / filename).read_bytes()) == expected, filename
    for name, expected in state["protected"].items():
        directory = ROOT.parent / name
        def git(*args):
            return subprocess.check_output(["git", "-C", str(directory), *args], stderr=subprocess.DEVNULL)
        assert git("rev-parse", "HEAD").decode().strip() == expected["head"]
        assert git("branch", "--show-current").decode().strip() == expected["branch"]
        assert git("status", "--porcelain=v1").decode().splitlines() == expected["status"]
        assert digest(git("diff", "--binary", "HEAD")) == expected["tracked_diff_sha256"]
        assert digest(git("diff", "--binary", "--cached")) == expected["index_diff_sha256"]
        changed_names = set(git("diff", "--name-only", "HEAD", "-z").decode().strip("\0").split("\0")) - {""}
        untracked_names = set(git("ls-files", "--others", "--exclude-standard", "-z").decode().strip("\0").split("\0")) - {""}
        actual = {name: digest((directory / name).read_bytes()) if (directory / name).is_file() else "absent"
                  for name in sorted(changed_names | untracked_names)}
        assert actual == expected["changed_files"]
    refs = dict(line.split(" ", 1) for line in subprocess.check_output(
        ["git", "for-each-ref", "--format=%(refname) %(objectname)"]).decode().splitlines())
    for name, value in state["refs"].items():
        if name not in {"refs/heads/feature/country-subnational-tb-model-v1",
                        "refs/remotes/origin/feature/country-subnational-tb-model-v1"}:
            assert refs[name] == value


def test_new_diagnostic_is_unwired_and_no_scientific_pipeline():
    module = ROOT / "engine/country_tb/catalytic_diagnostic.py"
    tree = ast.parse(module.read_text())
    modules = [node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
    assert not any(name and any(word in name for word in ("incidence", "demography", "apy", "static", "dynamic", "streamlit")) for name in modules)
    for folder in ("app", "ui", "general_ui", "explicit_tbi_ui", "engine/apy"):
        for path in (ROOT / folder).rglob("*.py"):
            assert "catalytic_diagnostic" not in path.read_text(encoding="utf-8-sig")
    changed = subprocess.check_output(["git", "diff", BASE, "--name-only"]).decode().splitlines()
    assert not any(p.startswith(("abm/", "app/", "ui/", "data/", "engine/apy/", "validation/")) for p in changed)
