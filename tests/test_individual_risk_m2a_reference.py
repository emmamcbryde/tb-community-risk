"""Independent publication arithmetic, package isolation and explicit limitations."""
from copy import deepcopy
from dataclasses import FrozenInstanceError,replace
from decimal import Decimal,localcontext
import csv
import json
import math
from pathlib import Path
import socket
import pytest
from tb_individual_risk import IndividualTBRiskModel,load_parameter_package,ParameterPackage,ValidationError,export_results
from tb_individual_risk.provenance import canonical,scientific_hash
from tb_individual_risk.menzies_2008 import prevalence_from_ari,ppv_from_components,linear_rate_projection
from tb_individual_risk.reference_comparison import reproduce_published_examples,formulation_comparisons,published_examples
from examples.individual_risk.reference_fixtures import fictional_reference_people

@pytest.fixture(scope="module")
def package():return load_parameter_package("menzies_2008_published_reference_v1")
@pytest.fixture(scope="module")
def historical(package):return IndividualTBRiskModel(package,formulation="menzies_2008_compatibility_v1")
@pytest.fixture
def person():return fictional_reference_people()[0]


def test_exact_published_equations_independent_decimal(historical,person):
    r=historical.predict(person)
    with localcontext() as context:
        context.prec=50
        tp=1-(1-Decimal('.001'))**25
        expected=tp/(tp+Decimal('.027')*Decimal('.1'))
    low=historical.predict(dict(person,tst_mm=12))
    assert low["tbi_probability"]==pytest.approx(float(expected),abs=1e-15)
    assert r["tbi_probability"]==1 and r["annual_risk"]==.001
    assert r["central_absolute_risk"]==.005 and r["historical_cumulative_risk_to_age80"]==.055


def test_linear_is_not_exponential(historical,person,monkeypatch):
    import tb_individual_risk.progression as progression
    monkeypatch.setattr(progression,"hazard_to_risk",lambda *a,**k:pytest.fail("Historical progression used exponential hazard helper"))
    p=dict(person,hiv=True,horizon_years=5)
    r=historical.predict(p)
    assert r["central_absolute_risk"]==.25 and r["central_absolute_risk"]!=pytest.approx(-math.expm1(-.25))
    assert r["historical_cumulative_risk_to_age80"] is None
    assert r["historical_arithmetic"]["cumulative_to_age80_literal_raw"]==2.75


def test_default_formulation_and_parameters_not_changed():
    m=IndividualTBRiskModel();r=m.predict(fictional_reference_people()[0])
    assert m.formulation=="hazard_based_progression_v1" and r["central_absolute_risk"] is None
    assert r["classification"]=="external_required"

@pytest.mark.parametrize("name,form",[("menzies_2008_published_reference_v1","hazard_based_progression_v1"),("synthetic_demonstration_v1","menzies_2008_compatibility_v1"),("contemporary_parameters_external_required_v1","menzies_2008_compatibility_v1")])
def test_packages_cannot_cross_formulations(name,form):
    with pytest.raises(ValidationError):IndividualTBRiskModel(load_parameter_package(name),formulation=form)


def test_immutable_packages_and_full_metadata(package):
    record=package.record();record["records"][0]["value"]=123
    assert package.value("baseline_annual_rate")==.001
    with pytest.raises(FrozenInstanceError):package._json="{}"
    with pytest.raises(ValidationError):ParameterPackage(canonical(record))
    required={"name","value","unit","definition","equation","domain","source_citation","source_page_table","source_population","publication_year","scope","review_status","limitations","adoption"}
    assert all(required<=set(r) for r in package.record()["records"])
    assert all(r["value"] is None for r in load_parameter_package("contemporary_parameters_external_required_v1").record()["records"])


def test_historical_marker_cannot_enter_hazard_config():
    p=IndividualTBRiskModel(load_parameter_package("synthetic_demonstration_v1")).parameters
    c=deepcopy(p.components);c["medical_risk_effects"].config["factors"]["diabetes"]["adoption"]="historical_only"
    with pytest.raises(ValidationError):replace(p,components=c)

@pytest.mark.parametrize("age,lifetime",[(25,.055),(79,.001),(80,0)])
def test_age80_equation(historical,person,age,lifetime):
    r=historical.predict(dict(person,current_age=age,horizon_years=0))
    assert r["historical_cumulative_risk_to_age80"]==pytest.approx(lifetime)


def test_above80_never_clipped(historical,person):
    r=historical.predict(dict(person,current_age=85,horizon_years=1))
    assert r["central_absolute_risk"] is None and r["historical_cumulative_risk_to_age80"] is None
    assert r["historical_arithmetic"]["cumulative_to_age80_literal_raw"]==-.005
    assert any("above80" in w for w in r["warnings"])


def test_recent_contact_and_conversion_window_and_conflict(historical,person):
    remote=historical.predict(person)
    for p in [dict(person,historical_contact_type="close",close_contact=True,time_since_infection=0),dict(person,recent_conversion=True,time_since_infection=0)]:
        r=historical.predict(dict(p,horizon_years=2))
        assert r["annual_risk"]==.015 and r["two_year_risk"]==.03 and r["central_absolute_risk"]==.03
        assert r["five_year_risk"] is None and r["historical_cumulative_risk_to_age80"] is None
        assert r["historical_arithmetic"]["footnote_schedule_raw_5.0"]==pytest.approx(.033)
    assert remote["annual_risk"]==.001


def test_unverified_general_combination_not_modern_max(historical,person):
    r=historical.predict(dict(person,diabetes=True,smoking=True,smoking_packs_per_day=2))
    assert r["central_absolute_risk"] is None and r["effective_rr"] is None
    assert "unresolved" in r["combination_policy"]
    with pytest.raises(ValidationError):historical.compare_policies([person])

@pytest.mark.parametrize("bcg_age",[.5,6])
def test_bcg_mapping_unresolved_not_progression_effect(historical,person,bcg_age):
    r=historical.predict(dict(person,bcg=True,bcg_age=bcg_age))
    assert r["tbi_probability"] is None and r["central_absolute_risk"] is None
    assert r["historical_arithmetic"]["annual_conditional_rate"]==.001
    assert any("BCG aggregate" in w for w in r["warnings"])

@pytest.mark.parametrize("state",["unknown","excluded","not_applicable"])
def test_no_compatibility_blank_defaults(historical,person,state):
    r=historical.predict(dict(person,diabetes={"state":state}))
    assert r["central_absolute_risk"] is None and r["input_states"]["diabetes"]["state"]==state
    assert r["lower_missing_scenario"] is None


def test_invalid_state_stays_structured(historical,person):
    r=historical.predict(dict(person,diabetes={"state":"invalid"}))
    assert r["status"]=="invalid" and r["validation_errors"]

@pytest.mark.parametrize("assay",["positive","negative"])
def test_igra_not_imported_from_current_site(historical,person,assay):
    r=historical.predict(dict(person,tst_result="not_done",tst_mm=None,igra_result=assay))
    assert r["tbi_probability"] is None and r["central_absolute_risk"] is None
    assert any("IGRA" in w for w in r["warnings"])


def test_historical_incidence_unit_conversion_not_country_download(historical,person):
    p=dict(person,historical_ari_probability=None,historical_smear_positive_incidence_per_100k=49)
    r=historical.predict(p)
    assert r["infection_calculation"]["ari"]==.01
    assert r["infection_calculation"]["p_ltbi"]==pytest.approx(1-.99**25)


def test_no_incidence_or_ntm_invented(historical,person):
    assert historical.predict(dict(person,historical_ari_probability=None))["tbi_probability"] is None
    assert historical.predict(dict(person,tst_mm=12,historical_ntm_sensitization_prevalence=None))["tbi_probability"] is None


def test_all_published_examples_and_partial_only(package):
    r=reproduce_published_examples(package)
    assert len(r["rows"])==42 and r["full_algorithm_matches"]==0
    assert all(not x["full_algorithm_reproduced"] for x in r["rows"])
    figure=next(x for x in r["rows"] if x["example_id"].endswith("_annual"))
    assert figure["reproduced_value"]==pytest.approx(.001811916)
    assert figure["absolute_difference"]==pytest.approx(.000001916)
    broken=next(x for x in r["rows"] if x["example_id"]=="T3_16")
    assert broken["within_declared_tolerance"] is False
    combined=[x for x in r["rows"] if x.get("group")=="granuloma_underweight"]
    assert all(x["reproduced_value"] is None for x in combined)


def test_linear_exponential_differences_not_hidden():
    r=formulation_comparisons();large=next(x for x in r["rows"] if x["scenario"]=="high_rr_long")
    assert large["historical_linear_raw"]==pytest.approx(3.33)
    assert large["linear_probability"] is None and 0<large["exponential_math_sensitivity"]<1
    assert not r["contemporary_adoption"]


def test_fixture_coverage_and_export(package,historical,tmp_path):
    people=fictional_reference_people();assert len(people)==26
    rows=[historical.predict(p) for p in people]
    c=tmp_path/"references.csv";j=tmp_path/"references.json";export_results(rows,csv_path=c,json_path=j)
    assert json.loads(j.read_text())==rows
    loaded=list(csv.DictReader(c.open()))
    assert [json.loads(x["full_result_json"]) for x in loaded]==rows
    assert all(x["formulation"]=="menzies_2008_compatibility_v1" for x in loaded)
    for row in rows:
        data=dict(row);data.pop("scientific_hash");assert scientific_hash(data)==row["scientific_hash"]


def test_no_network_either_formulation(historical,person,monkeypatch):
    def denied(*a,**k):pytest.fail("Network during ordinary prediction")
    monkeypatch.setattr(socket.socket,"connect",denied)
    assert historical.predict(person)["annual_risk"]==.001
    assert IndividualTBRiskModel(load_parameter_package("synthetic_demonstration_v1")).predict(person)["central_absolute_risk"] is not None


def test_source_version_and_package_in_all_outputs(historical,person):
    r=historical.predict(person)
    assert r["parameter_package"]=="menzies_2008_published_reference_v1"
    assert r["source_version"]=="Menzies_2008_publication_12_5_498_505"
    assert r["formulation"]=="menzies_2008_compatibility_v1"



def test_historical_input_units_and_known_risk_factors(historical):
    people={p["record_id"]:p for p in fictional_reference_people()}
    for name,rr in [("underweight",3),("heavy_smoking",3),("haemodialysis",25.3),("transplant",74),("infliximab",4)]:
        r=historical.predict(people[name])
        assert r["status"]!="invalid",r
        assert r["annual_risk"]==pytest.approx(.001*rr)
    assert prevalence_from_ari(0,1,.4)==.4



def test_explicit_modern_bcg_assay_changes_only_infection(person):
    p=IndividualTBRiskModel(load_parameter_package("synthetic_demonstration_v1")).parameters
    c=deepcopy(p.components)
    c["infection_probability"].config.update(mode="bayes_assays",assays=["tst"])
    c["tst_interpretation"].config.update(by_bcg={"false":{"sensitivity":.8,"specificity":.95},"true":{"sensitivity":.8,"specificity":.7}})
    m=IndividualTBRiskModel(replace(p,components=c))
    a=m.predict(dict(person,tst_mm=12,bcg=False));b=m.predict(dict(person,tst_mm=12,bcg=True))
    assert a["tbi_probability"]!=b["tbi_probability"]
    assert a["conditional_progression_risk"]==pytest.approx(b["conditional_progression_risk"])


def test_hazard_policy_comparison_keeps_package_provenance(person):
    m=IndividualTBRiskModel(load_parameter_package("synthetic_demonstration_v1"))
    comparisons=m.compare_policies([dict(person,diabetes=True,smoking=True)])
    for policy,rows in comparisons.items():
        assert rows[0]["formulation"]=="hazard_based_progression_v1"
        assert rows[0]["parameter_package"]=="synthetic_demonstration_v1"
        assert rows[0]["combination_policy"]==policy
        assert rows[0]["parameter_record"]["combination_policy"]==policy



def test_printed_smoking_and_exact_recency_boundaries_unresolved(historical,person,package):
    r=historical.predict(dict(person,smoking=True,smoking_packs_per_day=1))
    assert r["central_absolute_risk"] is None
    r=historical.predict(dict(person,time_since_infection=2))
    assert r["central_absolute_risk"] is None and any("boundary" in w for w in r["warnings"])
    entry=next(r for r in package.record()["records"] if r["name"]=="rr_smoking")
    assert "Figure1" in entry["definition"] and ">=1" in entry["definition"]
