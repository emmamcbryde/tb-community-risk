"""Bounded independent arithmetic, I/O and protection tests for isolated M1."""
from copy import deepcopy
from dataclasses import replace
from decimal import Decimal,localcontext
import csv
import hashlib
import json
import math
from pathlib import Path
import socket
import subprocess
import pytest

from tb_individual_risk import (IndividualTBRiskModel,ModelParameters,Parameter,PopulationPrior,
    PersonRecord,ValidationError,Observation,read_csv,read_json,write_records,export_results,SyntheticPopulationContract)
from tb_individual_risk.provenance import canonical,seal,scientific_hash
from tb_individual_risk.progression import hazard_to_risk,cumulative_hazard
from tb_individual_risk.risk_factors import combine_rr
from tb_individual_risk.infection_probability import Unavailable
from examples.individual_risk.synthetic_parameters import synthetic_parameters,fictional_people

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture
def parameters():return synthetic_parameters()

@pytest.fixture
def person():return fictional_people()[0]


def with_component(parameters,name,config):
    c=deepcopy(parameters.components);c[name]=replace(c[name],config=config)
    return replace(parameters,components=c)


def test_absent_unknown_and_explicit_states(person):
    absent=PersonRecord.from_dict(person).fields["diabetes"]
    missing=PersonRecord.from_dict(dict(person,diabetes=None)).fields["diabetes"]
    assert absent.state=="known_absent" and absent.value is False
    assert missing.state=="unknown" and missing.value is None
    for state in ("unknown","excluded","not_applicable"):
        assert PersonRecord.from_dict(dict(person,diabetes={"state":state})).fields["diabetes"].state==state

@pytest.mark.parametrize("value",[None,"unknown",{"state":"unknown"}])
def test_missing_boolean_not_false(person,value):
    assert PersonRecord.from_dict(dict(person,diabetes=value)).values()["diabetes"] is None


def test_complete_deterministic_and_parameter_roundtrip(parameters,person):
    model=IndividualTBRiskModel(parameters)
    assert canonical(model.predict(person))==canonical(model.predict(person))
    assert ModelParameters.from_dict(parameters.record()).record()==parameters.record()
    snapshot=model.parameters;snapshot.components["infection_probability"].config["mode"]="bad"
    assert model.predict(person)["central_absolute_risk"] is not None


def test_independent_worked_bayes_and_hazard(parameters,person):
    r=IndividualTBRiskModel(parameters).predict(person)
    with localcontext() as ctx:
        ctx.prec=50
        p=Decimal(".16")/(Decimal(".16")+Decimal(".08"))
        conditional=1-(-Decimal(".04")).exp()
        expected=p*conditional
    assert r["tbi_probability"]==pytest.approx(float(p),abs=1e-15)
    assert r["conditional_progression_risk"]==pytest.approx(float(conditional),abs=1e-15)
    assert r["central_absolute_risk"]==pytest.approx(float(expected),abs=1e-15)
    assert r["tbi_probability"]!=1 and r["tbi_probability"]!=r["conditional_progression_risk"]


def test_missing_only_through_prior_and_bounds(parameters,person):
    person["diabetes"]=None
    result=IndividualTBRiskModel(parameters).predict(person)
    expected=(.8*hazard_to_risk(.04,1)+.2*hazard_to_risk(.04,2))*(2/3)
    assert result["central_absolute_risk"]==pytest.approx(expected,abs=1e-15)
    assert result["lower_missing_scenario"]==pytest.approx((2/3)*hazard_to_risk(.04,1))
    assert result["upper_missing_scenario"]==pytest.approx((2/3)*hazard_to_risk(.04,2))
    assert result["input_contributions"]["diabetes"]["resolution"]=="marginalised_explicit_prior"
    no_prior=replace(parameters,population_prior=PopulationPrior())
    r=IndividualTBRiskModel(no_prior).predict(person)
    assert r["central_absolute_risk"] is None and r["lower_missing_scenario"] is None
    bounds=PopulationPrior(source="invented",status="synthetic_not_clinical",bounds={"diabetes":[False,True]})
    r=IndividualTBRiskModel(replace(parameters,population_prior=bounds)).predict(person)
    assert r["status"]=="bounded_only" and r["central_absolute_risk"] is None
    assert r["lower_missing_scenario"]>0


def test_two_factor_marginalisation_not_average_rr(parameters,person):
    r=IndividualTBRiskModel(parameters).predict(dict(person,diabetes=None,smoking=None))
    expected=(.56*hazard_to_risk(.04,1)+.14*hazard_to_risk(.04,2)+.3*hazard_to_risk(.04,3))*2/3
    assert r["central_absolute_risk"]==pytest.approx(expected,abs=1e-15)
    assert r["effective_rr"] is None and r["expected_effective_rr"]==pytest.approx(1.74)
    assert r["central_absolute_risk"]!=pytest.approx(hazard_to_risk(.04,1.74)*2/3,abs=1e-6)


def test_joint_marginalisation_preserves_conditional_infection_weighting(parameters,person):
    p=with_component(parameters,"infection_probability",{"mode":"supplied_probability"})
    states=({"values":{"tbi_probability":.1,"diabetes":False},"probability":.5},
            {"values":{"tbi_probability":.9,"diabetes":True},"probability":.5})
    prior=PopulationPrior(source="invented joint",status="synthetic_not_clinical",dependence_policy="explicit_joint",joint_states=states,joint_bound_states=tuple({"values":s["values"]} for s in states))
    model=IndividualTBRiskModel(replace(p,population_prior=prior));r=model.predict(dict(person,diabetes=None))
    expected=.5*.1*hazard_to_risk(.04,1)+.5*.9*hazard_to_risk(.04,2)
    assert r["central_absolute_risk"]==pytest.approx(expected)
    assert r["conditional_progression_risk"]==pytest.approx(expected/.5)
    conditioned=model.predict(dict(person,diabetes=True))
    assert conditioned["tbi_probability"]==.9


def test_explicit_dependence_and_scenario_budget(parameters,person):
    with pytest.raises(ValidationError):PopulationPrior(source="x",status="reviewed",distributions={"diabetes":[{"value":False,"probability":1}]})
    r=IndividualTBRiskModel(replace(parameters,max_scenarios=1)).predict(dict(person,diabetes=None))
    assert r["central_absolute_risk"] is None and any("budget" in w for w in r["warnings"])


def test_uninfected_no_reactivation_even_unresolved_progression(parameters,person):
    p=with_component(parameters,"infection_probability",{"mode":"supplied_probability"})
    r=IndividualTBRiskModel(p).predict(dict(person,tbi_probability=0,current_age=None,diabetes=None))
    assert r["central_absolute_risk"]==0 and r["conditional_progression_risk"] is None
    assert any("No future infection" in w for w in r["warnings"])

@pytest.mark.parametrize("policy,expected",[("highest_applicable_rr",3),("additive_excess_rr",4),("multiplicative_rr",6)])
def test_policies_and_hazard_application(parameters,person,policy,expected):
    r=IndividualTBRiskModel(replace(parameters,combination_policy=policy)).predict(dict(person,diabetes=True,smoking=True))
    assert r["effective_rr"]==expected
    assert r["central_absolute_risk"]==pytest.approx((2/3)*hazard_to_risk(.04,expected),abs=1e-15)
    assert r["central_absolute_risk"]!=pytest.approx((2/3)*hazard_to_risk(.04)*expected,abs=1e-6)


def test_highest_ties_and_no_silent_policy_switch():
    r=combine_rr({"a":3,"b":3,"c":2})
    assert r["effective_rr"]==3 and r["highest_contributors"]==["a","b"]
    assert combine_rr({})["effective_rr"]==1
    with pytest.raises(ValidationError):combine_rr({"a":3},"unknown")
    with pytest.raises(Unavailable):combine_rr({"a":1e300,"b":1e300},"multiplicative_rr")

@pytest.mark.parametrize("h,m",[(0,1),(1e-15,74),(.04,74),(100,1e100)])
def test_probability_stability(h,m):
    r=hazard_to_risk(h,m);assert 0<=r<=1
    assert r==pytest.approx(-math.expm1(-h*m))


def test_multiple_extreme_effects_warning(parameters,person):
    p=deepcopy(parameters);p.components["medical_risk_effects"].config["factors"]["hiv"]["rr"]=50
    p.components["medical_risk_effects"].config["factors"]["transplantation"]["rr"]=74
    r=IndividualTBRiskModel(p).predict(dict(person,hiv=True,transplantation=True))
    assert r["effective_rr"]==74 and any("Multiple extreme" in w for w in r["warnings"])
    assert r["central_absolute_risk"]<=1 and r["scenarios"][0]["factors"]["cap"] is None


def test_treated_hiv_and_pediatric_not_borrowed(parameters,person):
    model=IndividualTBRiskModel(parameters)
    assert model.predict(dict(person,hiv=True,hiv_treatment=True))["central_absolute_risk"] is None
    assert model.predict(dict(person,current_age=8,diabetes=True))["central_absolute_risk"] is None
    assert any("adult RR" in w for w in model.predict(dict(person,current_age=8,diabetes=True))["warnings"])
    assert model.predict(dict(person,current_age=8))["central_absolute_risk"] is None


def test_tail_never_silently_zero(parameters,person):
    p=with_component(parameters,"remote_infection_progression",{"tail_policy":"no_extrapolation"})
    r=IndividualTBRiskModel(p).predict(dict(person,time_since_infection=29))
    assert r["central_absolute_risk"] is None and any("not zero" in w for w in r["warnings"])
    assert cumulative_hazard(2,[[0,0],[2,.04],[10,.08]])==.04
    with pytest.raises(ValidationError):cumulative_hazard(1,[[0,0],[1,.1],[2,.05]])


def test_competing_death_analytic_and_depletion(parameters,person):
    p=with_component(parameters,"competing_mortality",{"mode":"constant_external_hazard","annual_hazard":.1})
    r=IndividualTBRiskModel(p).predict(person);s=r["scenarios"][0]["progression"]
    assert sum(s[k] for k in ("probability","competing_death_probability","surviving_disease_free_probability"))==pytest.approx(1)
    assert r["central_absolute_risk"]<IndividualTBRiskModel(parameters).predict(person)["central_absolute_risk"]

@pytest.mark.parametrize("field",["previous_active_tb","previous_preventive_treatment","active_tb_at_baseline"])
def test_history_requires_external_policy(parameters,person,field):
    r=IndividualTBRiskModel(parameters).predict(dict(person,**{field:True}));assert r["central_absolute_risk"] is None


def test_ranking_absolute_and_deterministic_ties(parameters,person):
    rows=[dict(person,record_id="z"),dict(person,record_id="a"),dict(person,record_id="higher",diabetes=True)]
    model=IndividualTBRiskModel(parameters);ranked=model.rank(rows)
    assert [r["record_id"] for r in ranked]==["higher","a","z"]
    assert [r["rank"] for r in ranked]==[1,2,2]
    assert [r["record_id"] for r in model.rank(list(reversed(rows)))]==["higher","a","z"]
    assert ranked[1]["overlapping_missing_information_ranges"] and ranked[1]["ranking_uncertain"]
    cases=model.rank(rows,by="expected_cases_per_1000");assert cases[0]["ranking_score"]==1000*ranked[0]["central_absolute_risk"]
    with pytest.raises(ValidationError):model.rank(rows,by="conditional_progression_risk")
    assert model.rank(rows,by="conditional_progression_risk",likely_infected_threshold=.5)[0]["rank"]==1


def test_absolute_default_differs_from_conditional(parameters,person):
    rows=[dict(person,record_id="high_infection",pretest_tbi_probability=.8),dict(person,record_id="low_infection_high_progression",pretest_tbi_probability=.001,haemodialysis=True)]
    m=IndividualTBRiskModel(parameters)
    assert m.rank(rows)[0]["record_id"]=="high_infection"
    assert m.rank(rows,by="conditional_progression_risk",likely_infected_threshold=0)[0]["record_id"]=="low_infection_high_progression"


def test_overlapping_bounds_not_meaningful_order(parameters,person):
    rows=[dict(person,record_id="a",diabetes=None),dict(person,record_id="b",smoking=None)]
    result=IndividualTBRiskModel(parameters).rank(rows)
    assert all(r["ranking_uncertain"] and r["overlapping_missing_information_ranges"] for r in result)


def test_mixed_horizons_not_silently_ranked(parameters,person):
    rows=[dict(person,record_id="a",horizon_years=1),dict(person,record_id="b",horizon_years=5)]
    model=IndividualTBRiskModel(parameters)
    assert all(r["rank"] is None for r in model.rank(rows))
    assert all(r["rank"] is not None for r in model.rank(rows,horizon_years=5))

@pytest.mark.parametrize("changes",[{"diabetes":"invalid"},{"current_age":-1},{"time_since_infection":100},{"pretest_tbi_probability":1.1},{"diabetes":{"state":"known_absent","value":True}},{"birth_year":2020},{"typo_factor":True},{"current_age":float("inf")}])
def test_invalid_not_discarded(parameters,person,changes):
    result=IndividualTBRiskModel(parameters).rank([dict(person,**changes)])
    assert len(result)==1 and result[0]["status"]=="invalid" and result[0]["validation_errors"]
    assert result[0]["rank"] is None


def test_duplicate_ids_preserved_as_invalid(parameters,person):
    rows=IndividualTBRiskModel(parameters).rank([person,person])
    assert len(rows)==2 and all(r["status"]=="invalid" and r["rank"] is None for r in rows)


def test_csv_json_and_full_precision_roundtrip(tmp_path,parameters,person):
    people=[dict(person,record_id="001",diabetes=None),dict(person,record_id="comma,quoted",current_age="40.00000000000000001")]
    c=tmp_path/"people.csv";j=tmp_path/"people.json";write_records(people,csv_path=c,json_path=j)
    assert read_csv(c)==people and read_json(j)==people
    model=IndividualTBRiskModel(parameters);rows=model.rank(read_csv(c))
    export_results(rows,csv_path=c,json_path=j)
    assert json.loads(j.read_text())==rows
    exported=list(csv.DictReader(c.open()))
    assert [json.loads(r["full_result_json"]) for r in exported]==rows
    assert all(float(r["central_absolute_risk"])==original["central_absolute_risk"] for r,original in zip(exported,rows))


def test_dataframe_missing_preserved(parameters,person):
    pd=pytest.importorskip("pandas")
    frame=pd.DataFrame([dict(person,record_id="a"),dict(person,record_id="b",diabetes=None)])
    rows=IndividualTBRiskModel(parameters).rank(frame)
    assert len(rows)==2 and next(r for r in rows if r["record_id"]=="b")["missing_critical_variables"]==["diabetes"]


def test_external_parameters_and_synthetic_cannot_be_reviewed(parameters,person):
    r=IndividualTBRiskModel().predict(person)
    assert r["status"]=="unavailable" and r["central_absolute_risk"] is None
    with pytest.raises(ValidationError):replace(parameters,classification="reviewed")
    r=IndividualTBRiskModel(parameters).predict(person)
    assert r["classification"]=="synthetic_not_clinical" and not r["clinical_use_approved"]
    assert any("INVENTED" in w for w in r["warnings"])
    tampered=parameters.record();tampered["max_scenarios"]=2
    with pytest.raises(ValidationError):ModelParameters.from_dict(tampered)


def test_population_design_requires_joint_policy_and_provenance():
    p=SyntheticPopulationContract("example",{"diabetes":.1},{"diabetes":"invented"})
    assert p.record()["status"]=="design_only_not_implemented" and p.joint_policy=="external_required"
    with pytest.raises(ValidationError):SyntheticPopulationContract("x",{"diabetes":.1},{})


def test_no_prediction_network(parameters,person,monkeypatch):
    def denied(*a,**k):pytest.fail("Unexpected network access")
    monkeypatch.setattr(socket,"create_connection",denied);monkeypatch.setattr(socket.socket,"connect",denied)
    assert IndividualTBRiskModel(parameters).predict(person)["central_absolute_risk"]>0


def test_frozen_files_and_no_existing_ui_import():
    baseline=json.loads((ROOT/"docs/individual_risk_model/protected_start.json").read_text())
    for name,expected in baseline["tracked_files"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
        if name.endswith(".py"):
            assert "tb_individual_risk" not in (ROOT/name).read_text(encoding="utf-8"),name
    assert subprocess.check_output(["git","rev-parse","sa-health-apy-he-v1.0.0^{}"],cwd=ROOT).decode().strip()==baseline["starting_head"]


def test_known_unmodelled_finding_not_silently_ignored(parameters,person):
    r=IndividualTBRiskModel(parameters).predict(dict(person,factors={"new_condition":True}))
    assert r["central_absolute_risk"] is None and any("lack an effect" in w for w in r["warnings"])


def test_parameter_support_validation(parameters):
    with pytest.raises(ValidationError):
        IndividualTBRiskModel(replace(parameters,population_prior=PopulationPrior(source="x",status="synthetic_not_clinical",dependence_policy="explicit_independence",distributions={"diabetes":[{"value":0,"probability":1}]})))


def test_example_policy_comparison_and_ids(parameters):
    rows=fictional_people();result=IndividualTBRiskModel(parameters).compare_policies(rows)
    assert len(rows)>=10 and len(result)==3
    assert all(len(v)==len(rows) for v in result.values())
    assert {r["record_id"] for r in result["highest_applicable_rr"]}=={r["record_id"] for r in rows}


@pytest.mark.parametrize("state",["excluded","not_applicable"])
def test_excluded_not_overwritten_by_prior(parameters,person,state):
    r=IndividualTBRiskModel(parameters).predict(dict(person,diabetes={"state":state}))
    assert r["central_absolute_risk"] is None and r["lower_missing_scenario"] is None
    assert r["input_contributions"]["diabetes"]["input_state"]==state


def test_tst_and_igra_are_distinct_observation_contracts(parameters,person):
    p=with_component(parameters,"infection_probability",{"mode":"bayes_assays","assays":["tst"]})
    r=IndividualTBRiskModel(p).predict(dict(person,tst_mm=12,tst_result=None))
    assert r["tbi_probability"]==pytest.approx(.16/(.16+.12))
    mismatch=IndividualTBRiskModel(p).predict(dict(person,tst_mm=12,tst_result="negative"))
    assert mismatch["central_absolute_risk"] is None
    assert IndividualTBRiskModel(parameters).predict(dict(person,igra_product="OTHER"))["central_absolute_risk"] is None
    assert IndividualTBRiskModel(parameters).predict(dict(person,igra_result="indeterminate"))["central_absolute_risk"] is None
    p=with_component(parameters,"infection_probability",{"mode":"bayes_assays","assays":["tst","igra"]})
    assert IndividualTBRiskModel(p).predict(dict(person,tst_result="positive"))["central_absolute_risk"] is None
    p=with_component(p,"infection_probability",{"mode":"bayes_assays","assays":["tst","igra"],"joint_policy":"conditionally_independent_given_infection"})
    assert IndividualTBRiskModel(p).predict(dict(person,tst_result="positive"))["central_absolute_risk"] is not None


def test_unresolved_effect_component_never_defaults_to_rr_one(parameters,person):
    c=deepcopy(parameters.components);c["medical_risk_effects"]=Parameter("medical_risk_effects")
    r=IndividualTBRiskModel(replace(parameters,components=c)).predict(person)
    assert r["central_absolute_risk"] is None and any("external_parameter_required" in w for w in r["warnings"])


def test_original_unknown_and_sources_retained(parameters,person):
    person["diabetes"]={"state":"unknown","source":{"table":"fictional_table","row":"017"}}
    r=IndividualTBRiskModel(parameters).predict(person)
    assert r["original_record"]==person
    assert r["input_contributions"]["diabetes"]["source"]==person["diabetes"]["source"]
    source=deepcopy(person);source["current_age"]=40.5;source["birth_year"]=1986
    assert PersonRecord.from_dict(source).fields["birth_year"].known


def test_csv_invalid_rows_and_duplicate_headers(tmp_path,parameters):
    p=tmp_path/"bad.csv";p.write_text("record_id,current_age,diabetes\n001,20,not-a-boolean\n002,40,False,extra\n")
    rows=IndividualTBRiskModel(parameters).rank(read_csv(p))
    assert len(rows)==2 and all(r["status"]=="invalid" for r in rows)
    p.write_text("record_id,record_id\na,b\n")
    with pytest.raises(ValidationError):read_csv(p)


def test_bounds_not_implied_from_prior_support(parameters,person):
    p=replace(parameters,population_prior=replace(parameters.population_prior,bounds={}))
    r=IndividualTBRiskModel(p).predict(dict(person,diabetes=None))
    assert r["central_absolute_risk"] is not None and r["lower_missing_scenario"] is None
    assert r["upper_missing_scenario"] is None


def test_evidence_table_and_all_parameter_families():
    rows=list(csv.DictReader((ROOT/"docs/individual_risk_model/evidence.csv").open(encoding="utf-8")))
    assert len(rows)>=20 and all(r["source"] and r["review_status"] and r["pediatric_applicability"] for r in rows)
    assert len(ModelParameters().components)==11
    assert all(c.status=="external_parameter_required" for c in ModelParameters().components.values())



def test_additive_overflow_is_structured_unavailable():
    with pytest.raises(Unavailable):combine_rr({"a":1e308,"b":1e308},"additive_excess_rr")


def test_standard_csv_boolean_spelling(tmp_path):
    p=tmp_path/"bool.csv";p.write_text("record_id,diabetes,smoking\na,True,False\nb,true,false\n")
    rows=read_csv(p)
    assert all(r["diabetes"] is True and r["smoking"] is False for r in rows)


def test_saturated_infection_marginal_stays_within_bounds(parameters,person):
    p=with_component(parameters,"infection_probability",{"mode":"supplied_probability"})
    r=IndividualTBRiskModel(p).predict(dict(person,tbi_probability=1,diabetes=None,smoking=None))
    assert r["tbi_probability"]==1 and 0<=r["central_absolute_risk"]<=1



def test_heterogeneous_csv_preserves_omitted_versus_explicit_unknown(tmp_path):
    rows=[{"record_id":"001","diabetes":None},{"record_id":"002","smoking":False},{"record_id":"003","factors":{"extra":True}}]
    c=tmp_path/"heterogeneous.csv";j=tmp_path/"heterogeneous.json"
    write_records(rows,csv_path=c,json_path=j)
    assert read_csv(c)==read_json(j)==rows



def test_effects_require_explicit_age_applicability(parameters,person):
    p=deepcopy(parameters);p.components["medical_risk_effects"].config["factors"]["diabetes"].pop("age_min")
    r=IndividualTBRiskModel(p).predict(dict(person,diabetes=True))
    assert r["central_absolute_risk"] is None and any("age applicability" in w for w in r["warnings"])


def test_unrecognized_effect_and_curve_status_not_accepted(parameters,person):
    p=deepcopy(parameters);p.components["medical_risk_effects"].config["factors"]["diabetes"]["status"]="not_reviewed"
    assert IndividualTBRiskModel(p).predict(dict(person,diabetes=True))["central_absolute_risk"] is None
    p=deepcopy(parameters);p.components["baseline_progression_hazard"].config["age_curves"][0].pop("status")
    assert IndividualTBRiskModel(p).predict(person)["central_absolute_risk"] is None


@pytest.mark.parametrize("identifier",[None,{"state":"unknown"},{"state":"excluded"},17])
def test_identifier_unknown_supported_without_invention(parameters,person,identifier):
    r=IndividualTBRiskModel(parameters).predict(dict(person,record_id=identifier))
    assert r["central_absolute_risk"] is not None
    assert r["record_id"]==(17 if identifier==17 else None)
    assert r["original_record"]["record_id"]==identifier


def test_identifier_types_and_missing_identifier_csv_roundtrip(tmp_path,parameters,person):
    rows=[dict(person,record_id=17),dict(person,record_id=None),dict(person,record_id={"state":"unknown"})]
    no_id=dict(person);no_id.pop("record_id");rows.append(no_id)
    c=tmp_path/"id.csv";j=tmp_path/"id.json";write_records(rows,csv_path=c,json_path=j)
    assert read_csv(c)==read_json(j)==rows
    assert len(IndividualTBRiskModel(parameters).rank(rows))==len(rows)



def test_inconsistent_identifier_state_is_structured_error(parameters,person):
    r=IndividualTBRiskModel(parameters).rank([dict(person,record_id={"state":"known_present"})])
    assert len(r)==1 and r[0]["status"]=="invalid" and r[0]["validation_errors"]


def test_csv_identifier_metadata_cannot_hide_supplied_id(tmp_path,parameters):
    p=tmp_path/"bad_metadata.csv";p.write_text("record_id,__identifier_kind\noriginal,unknown\n")
    r=IndividualTBRiskModel(parameters).rank(read_csv(p))
    assert len(r)==1 and r[0]["status"]=="invalid"
