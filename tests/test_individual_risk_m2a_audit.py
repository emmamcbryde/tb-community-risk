"""Independent M1 audit regressions, written before reference parameters."""
from copy import deepcopy
from dataclasses import replace
import json
import socket
import pytest
from tb_individual_risk import IndividualTBRiskModel,ValidationError,ModelParameters
from tb_individual_risk.provenance import canonical
from examples.individual_risk.synthetic_parameters import synthetic_parameters,fictional_people


def test_unclassified_radiograph_does_not_become_reference_risk():
    result=IndividualTBRiskModel(synthetic_parameters()).predict(dict(fictional_people()[0],cxr="other"))
    assert result["central_absolute_risk"] is None
    assert any("radiograph" in w.lower() for w in result["warnings"])


def test_nested_synthetic_evidence_cannot_be_relabelled_reviewed():
    p=synthetic_parameters()
    components={k:replace(v,status="reviewed") for k,v in p.components.items()}
    with pytest.raises(ValidationError):
        replace(p,components=components,population_prior=replace(p.population_prior,status="reviewed"),classification="reviewed")


def test_conditional_ranking_uses_conditional_ranges():
    base=fictional_people()[0]
    rows=[dict(base,record_id="high",pretest_tbi_probability=.8,smoking=None),dict(base,record_id="low",pretest_tbi_probability=.001,smoking=None)]
    result=IndividualTBRiskModel(synthetic_parameters()).rank(rows,by="conditional_progression_risk",likely_infected_threshold=0)
    assert all(r["overlapping_missing_information_ranges"] for r in result)
    assert all(r["lower_conditional_missing_scenario"] is not None for r in result)


def test_parameter_provenance_can_reconstruct_prediction():
    p=fictional_people()[0];m=IndividualTBRiskModel(synthetic_parameters());r=m.predict(p)
    assert canonical(IndividualTBRiskModel(ModelParameters.from_dict(r["parameter_record"])).predict(p))==canonical(r)


def test_tie_order_preserves_identifier_types():
    base=fictional_people()[0];rows=[dict(base,record_id=1),dict(base,record_id="1")]
    model=IndividualTBRiskModel(synthetic_parameters())
    assert [canonical(r["record_id"]) for r in model.rank(rows)]==[canonical(r["record_id"]) for r in model.rank(list(reversed(rows)))]


def test_offline_prediction_and_joint_product_identity(monkeypatch):
    def denied(*a,**k):pytest.fail("Network during prediction")
    monkeypatch.setattr(socket.socket,"connect",denied)
    r=IndividualTBRiskModel(synthetic_parameters()).predict(dict(fictional_people()[0],diabetes=None))
    assert r["central_absolute_risk"]==pytest.approx(r["tbi_probability"]*r["conditional_progression_risk"])



def test_adoption_metadata_is_not_mistaken_for_evidence_status():
    from tb_individual_risk import Parameter,PopulationPrior
    base=ModelParameters()
    components={name:Parameter(name,status="reviewed",source="test_only_contract_evidence",config={"adoption":"adopted"}) for name in base.components}
    result=replace(base,components=components,population_prior=PopulationPrior(source="no_imputation_test_only",status="reviewed"),classification="reviewed")
    assert result.classification=="reviewed"
