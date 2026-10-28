"""Focused software checks for the isolated architecture proposals."""

from copy import deepcopy
import ast
import hashlib
import json
import math
from pathlib import Path
import random

import pytest
from jsonschema import Draft202012Validator, ValidationError

from engine.country_tb.contracts import (
    FACTOR_IDS, association_summary, canonical_bytes, content_hash,
    feasible_association, reconcile_age_counts, validate_age_bands,
    validate_overlap, validate_prevalence,
)


ROOT = Path(__file__).resolve().parents[2]
SCHEMA_DIR = ROOT / "schemas" / "country_tb" / "v1"
SCHEMA = json.loads((SCHEMA_DIR / "contracts.schema.json").read_text(encoding="utf-8"))
EXAMPLES = json.loads((SCHEMA_DIR / "review_examples.json").read_text(encoding="utf-8"))
BANDS = [
    {"age_lower": lo, "age_upper_exclusive": hi}
    for lo, hi in ((0, 5), (5, 15), (15, 45), (45, 65), (65, None))
]


def prevalence(**updates):
    return {
        "factor_id": "diabetes", "state": "value", "value_percent": 0,
        "lower_percent": None, "upper_percent": None, "reason": None, **updates,
    }


def test_schema_and_all_seven_provenance_round_trips():
    Draft202012Validator.check_schema(SCHEMA)
    validator = Draft202012Validator(SCHEMA)
    assert len(EXAMPLES) == 7
    for payload in EXAMPLES.values():
        validator.validate(payload)
        restored = json.loads(canonical_bytes(payload))
        assert restored == payload
        assert restored["provenance"] == payload["provenance"]
        assert content_hash(restored) == content_hash(payload)


def test_canonical_serialization_and_hash_have_fixed_bytes():
    payload = {"z": 0, "a": ["é", 0.12345678901234568, None]}
    expected = b'{"a":["\\u00e9",0.12345678901234568,null],"z":0}'
    assert canonical_bytes(payload) == expected
    assert content_hash(payload) == hashlib.sha256(expected).hexdigest()
    assert content_hash(dict(reversed(list(payload.items())))) == content_hash(payload)
    assert content_hash({"a": [1, 2]}) != content_hash({"a": [2, 1]})
    assert content_hash({"x": 1}) != content_hash({"x": 1.0})


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), (1, 2), {1: "bad"}])
def test_non_json_hash_inputs_are_rejected(bad):
    with pytest.raises(ValueError):
        canonical_bytes({"bad": bad})


def test_provenance_changes_are_analytical():
    original = EXAMPLES["geography"]
    changed = deepcopy(original)
    changed["provenance"][0]["review_status"] = "reviewed"
    assert content_hash(changed) != content_hash(original)
    changed = deepcopy(original)
    changed["provenance"][0]["sha256"] = "1" * 64
    assert content_hash(changed) != content_hash(original)


@pytest.mark.parametrize("value,expected", [(0, 0), (100, 1), (3.141592653589793, .03141592653589793)])
def test_percentage_limits_and_precision(value, expected):
    row = prevalence(value_percent=value)
    assert validate_prevalence(row) == expected
    assert json.loads(canonical_bytes(row))["value_percent"] == value


@pytest.mark.parametrize("value", [-0.01, 100.01, None, True, "10", float("nan"), float("inf")])
def test_invalid_prevalence_rejected_without_clipping(value):
    with pytest.raises(ValueError):
        validate_prevalence(prevalence(value_percent=value))


def test_missing_zero_excluded_and_not_applicable_are_distinct():
    rows = [prevalence()]
    for state in ("missing", "excluded", "not_applicable"):
        rows.append(prevalence(state=state, value_percent=None, reason="Declared test reason"))
    assert validate_prevalence(rows[0]) == 0
    assert [validate_prevalence(row) for row in rows[1:]] == [None, None, None]
    assert len({content_hash(row) for row in rows}) == 4
    with pytest.raises(ValueError):
        validate_prevalence(prevalence(state="missing"))
    with pytest.raises(ValueError):
        validate_prevalence(prevalence(state="excluded", value_percent=None, reason=" "))


@pytest.mark.parametrize("factor", ["recent_infection", "close_contact", "aids"])
def test_recency_contact_and_duplicate_hiv_factor_are_not_inputs(factor):
    with pytest.raises(ValueError):
        validate_prevalence(prevalence(factor_id=factor))
    payload = deepcopy(EXAMPLES["risk_factor_prevalence"])
    payload["rows"][0]["factor_id"] = factor
    with pytest.raises(ValidationError):
        Draft202012Validator(SCHEMA).validate(payload)
    assert len(FACTOR_IDS) == 11


@pytest.mark.parametrize("lower,upper", [(None, 5), (6, 10), (-1, 10), (0, 101), (0, 4)])
def test_uncertainty_bounds_need_both_and_must_enclose(lower, upper):
    with pytest.raises(ValueError):
        validate_prevalence(prevalence(value_percent=5, lower_percent=lower, upper_percent=upper))


def test_valid_prevalence_bounds():
    assert validate_prevalence(prevalence(value_percent=5, lower_percent=2, upper_percent=8)) == .05


def test_age_counts_reconcile_without_splitting_or_losing_tail():
    source = [{"age_lower": age, "age_upper_exclusive": age + 1, "count": 1} for age in range(65)]
    source.append({"age_lower": 65, "age_upper_exclusive": None, "count": 35})
    result = reconcile_age_counts(source, BANDS)
    assert [row["count"] for row in result] == [5, 10, 30, 20, 35]
    assert sum(row["count"] for row in result) == 100
    assert sum(row["proportion"] for row in result) == 1
    assert result[-1]["age_upper_exclusive"] is None


@pytest.mark.parametrize("bands", [
    [], [{"age_lower": 1, "age_upper_exclusive": None}],
    [{"age_lower": 0, "age_upper_exclusive": 5}, {"age_lower": 6, "age_upper_exclusive": None}],
    [{"age_lower": 0, "age_upper_exclusive": 6}, {"age_lower": 5, "age_upper_exclusive": None}],
    [{"age_lower": 0, "age_upper_exclusive": 5}],
    [{"age_lower": False, "age_upper_exclusive": None}],
])
def test_invalid_age_partitions_fail(bands):
    with pytest.raises(ValueError):
        validate_age_bands(bands)


def test_coarse_age_band_requires_explicit_disaggregation():
    source = [
        {"age_lower": 0, "age_upper_exclusive": 15, "count": 20},
        {"age_lower": 15, "age_upper_exclusive": None, "count": 80},
    ]
    with pytest.raises(ValueError, match="allocation rule"):
        reconcile_age_counts(source, BANDS)
    with pytest.raises(ValueError):
        reconcile_age_counts([{**row, "count": 0} for row in BANDS], BANDS)
    with pytest.raises(ValueError):
        reconcile_age_counts([{**row, "count": -1} for row in BANDS], BANDS)


@pytest.mark.parametrize("p,q,low,high", [
    (0, .2, 0, 0), (1, .2, .2, .2), (.3, 0, 0, 0), (.3, 1, .3, .3),
    (0, 0, 0, 0), (1, 1, 1, 1), (.3, .2, 0, .2),
    (.9, .8, .7, .8), (1e-300, 1, 1e-300, 1e-300),
])
def test_frechet_endpoints_independence_and_phi(p, q, low, high):
    neg, zero, pos = [feasible_association(p, q, a) for a in (-1, 0, 1)]
    assert neg["overlap"] == low
    assert zero["overlap"] == pytest.approx(p * q)
    assert pos["overlap"] == high
    for result in (neg, zero, pos):
        assert min(result["cells"].values()) >= 0
        assert sum(result["cells"].values()) == pytest.approx(1)
    if p in (0, 1) or q in (0, 1):
        assert zero["phi"] is None
        assert zero["phi_status"] == "undefined_degenerate_marginal"
    else:
        assert zero["phi"] == 0
        expected = (pos["overlap"] - p * q) / math.sqrt(p * (1-p) * q * (1-q))
        assert pos["phi"] == pytest.approx(expected)


def test_scale_is_not_pearson_and_rejects_invalid_joint():
    result = feasible_association(.3, .2, 1)
    assert result["requested_alpha"] == 1
    assert result["overlap"] == .2
    assert result["phi"] == pytest.approx(.7637626158259733)
    assert feasible_association(.3, .2, -1)["phi"] == pytest.approx(-.32732683535398854)
    with pytest.raises(ValueError):
        validate_overlap(.3, .2, .21)
    with pytest.raises(ValueError):
        validate_overlap(.9, .8, .69)


@pytest.mark.parametrize("p,q,alpha", [
    (-.1,.2,0), (1.1,.2,0), (.2,float("nan"),0), (.2,.3,1.01),
    (.2,.3,-1.01), (.2,.3,True), (.2,.3,"0"), (True,.3,0),
])
def test_invalid_association_inputs_fail(p, q, alpha):
    with pytest.raises(ValueError):
        feasible_association(p, q, alpha)


def test_age_mediated_association_survives_zero_within_age_control():
    result = association_summary([
        {"count": 50, "p": .1, "q": .1, "alpha": 0},
        {"count": 50, "p": .5, "q": .5, "alpha": 0},
    ])
    assert [row["phi"] for row in result["strata"]] == [0, 0]
    assert result["overall"]["overlap"] == .13
    assert result["overall"]["phi"] == pytest.approx(.04 / .21)
    assert result["within_age_covariance"] == 0
    assert result["age_driven_covariance"] == pytest.approx(.04)


def test_generated_joint_grid_preserves_marginals():
    rng = random.Random(723)
    for _ in range(100):
        p, q, alpha = rng.random(), rng.random(), rng.uniform(-1, 1)
        result = feasible_association(p, q, alpha)
        c = result["cells"]
        assert sum(c.values()) == pytest.approx(1)
        assert c["tbi_and_risk"] + c["tbi_only"] == pytest.approx(p)
        assert c["tbi_and_risk"] + c["risk_only"] == pytest.approx(q)
        assert -1 <= result["phi"] <= 1


def test_structural_schema_rejects_percentage_and_source_corruption():
    validator = Draft202012Validator(SCHEMA)
    for mutation in ("percent", "missing_value", "source_hash", "recency"):
        payload = deepcopy(EXAMPLES["risk_factor_prevalence"])
        if mutation == "percent":
            payload["rows"][0]["value_percent"] = 101
        elif mutation == "missing_value":
            payload["rows"][0]["state"] = "missing"
        elif mutation == "source_hash":
            payload["provenance"][0]["sha256"] = "bad"
        else:
            payload["rows"][0]["factor_id"] = "recent_infection"
        with pytest.raises(ValidationError):
            validator.validate(payload)


def test_no_protected_app_runner_or_engine_imports_country_prototype():
    for directory in ("app", "ui", "scripts", "engine"):
        for path in (ROOT / directory).rglob("*.py"):
            if (ROOT / "engine" / "country_tb") in path.parents:
                continue
            if path in {ROOT / "scripts" / "preview_country_tb_data.py",
                        ROOT / "scripts" / "preview_country_tb_calibration.py"}:
                continue  # Explicitly authorized isolated developer CLIs.
            source = path.read_text(encoding="utf-8-sig")
            assert "engine.country_tb" not in source, path
            for node in ast.walk(ast.parse(source, filename=str(path))):
                if isinstance(node, ast.Import):
                    assert all("country_tb" not in alias.name for alias in node.names), path
                elif isinstance(node, ast.ImportFrom):
                    assert "country_tb" not in (node.module or ""), path
