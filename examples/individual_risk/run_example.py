"""Run invented arithmetic examples. Not for clinical use."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tb_individual_risk import IndividualTBRiskModel,export_results,write_records
from synthetic_parameters import synthetic_parameters,fictional_people


def main():
    parser=argparse.ArgumentParser(description="INVENTED SYNTHETIC demonstration, not clinical")
    parser.add_argument("--output",type=Path,default=ROOT/"abm/output/individual_risk_m1")
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    model=IndividualTBRiskModel(synthetic_parameters());people=fictional_people()
    write_records(people,csv_path=args.output/"fictional_people.csv",json_path=args.output/"fictional_people.json")
    comparison=model.compare_policies(people)
    for policy,result in comparison.items():
        export_results(result,csv_path=args.output/(policy+".csv"),json_path=args.output/(policy+".json"))
    print("INVENTED SYNTHETIC DATA AND PARAMETERS; NOT FOR CLINICAL USE")
    for row in comparison["highest_applicable_rr"]:
        print(row["record_id"],"rank=",row["rank"],"absolute risk=",row["central_absolute_risk"],"status=",row["status"])
    print("Full precision CSV/JSON written to",args.output)

if __name__=="__main__":main()
