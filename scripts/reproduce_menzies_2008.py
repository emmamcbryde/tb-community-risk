"""Offline fictional reproduction and formulation comparison. No clinical submissions."""
import argparse
import csv
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tb_individual_risk import IndividualTBRiskModel,load_parameter_package,export_results
from tb_individual_risk.provenance import canonical,seal
from tb_individual_risk.reference_comparison import reproduce_published_examples,formulation_comparisons
from tb_individual_risk.risk_factors import combine_rr
from examples.individual_risk.reference_fixtures import fictional_reference_people


def run(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    package=load_parameter_package("menzies_2008_published_reference_v1")
    model=IndividualTBRiskModel(package,formulation="menzies_2008_compatibility_v1")
    hazard=IndividualTBRiskModel(load_parameter_package("synthetic_demonstration_v1"),formulation="hazard_based_progression_v1")
    people=fictional_reference_people();results=[]
    for person in people:
        reference=model.predict(person)
        modern=hazard.predict(person)
        rates=reference.get("factor_calculation",{}).get("individual_rr",{})
        policies={p:combine_rr(rates,p) for p in ("highest_applicable_rr","additive_excess_rr","multiplicative_rr")}
        horizons={}
        for t in (1,2,5):
            horizons[str(t)]={"historical":model.predict(dict(person,horizon_years=t)),"synthetic_hazard":hazard.predict(dict(person,horizon_years=t))}
        results.append(seal(dict(fictional=True,record_id=person["record_id"],person=person,historical=reference,
            synthetic_hazard=modern,equivalent_horizons=horizons,modern_combination_math_only=policies,
            comparison_warning="Synthetic hazard parameters deliberately differ; no historical values adopted as contemporary defaults")))
    reproduction=reproduce_published_examples(package);comparison=formulation_comparisons()
    for name,value in (("reference_fixtures",results),("published_reproduction",reproduction),("formulation_comparison",comparison)):
        (output/(name+".json")).write_text(canonical(value)+"\n",encoding="utf-8")
    export_results([r["historical"] for r in results],csv_path=output/"fictional_reference.csv",json_path=output/"fictional_reference.json")
    fields=("example_id","published_value","reproduced_value","absolute_difference","relative_difference","declared_tolerance","within_declared_tolerance","inputs_completely_specified","full_algorithm_reproduced","unresolved_implementation_detail")
    with (output/"published_reproduction.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows({k:r.get(k) for k in fields} for r in reproduction["rows"])
    return dict(fictional_people=len(results),published_rows=len(reproduction["rows"]),full_algorithm_matches=0,
        partial_comparisons=reproduction["partial_comparisons"],partial_within_printed_resolution=reproduction["partial_within_printed_resolution"],
        partial_outside_resolution=[r["example_id"] for r in reproduction["rows"] if r["within_declared_tolerance"] is False],
        no_clinical_validation=True,parameter_package_hash=package.record()["scientific_hash"])


def main():
    parser=argparse.ArgumentParser(description="Historical reproduction only, not clinical validation")
    parser.add_argument("--output",type=Path,default=ROOT/"abm/output/individual_risk_m2a/reproduction")
    args=parser.parse_args();print(json.dumps(run(args.output)))

if __name__=="__main__":main()
