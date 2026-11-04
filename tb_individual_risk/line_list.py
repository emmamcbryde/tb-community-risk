"""Local Python/DataFrame/CSV/JSON adapters. Invalid rows remain present."""
import csv
import io
import json
from pathlib import Path
from .contracts import ValidationError,Observation,PersonRecord
from .provenance import canonical


def records(source):
    if isinstance(source,dict):return [source]
    if hasattr(source,"to_dict") and hasattr(source,"columns"):
        return source.to_dict(orient="records")
    if isinstance(source,(list,tuple)):return list(source)
    raise ValidationError("line_list","Supply dictionaries/DataFrame, or read_csv/read_json for local files.")


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as stream:
        reader=csv.DictReader(stream)
        if not reader.fieldnames or len(reader.fieldnames)!=len(set(reader.fieldnames)):
            raise ValidationError("CSV","Missing or duplicate header.")
        result=[]
        for row in reader:
            if None in row:
                result.append({"record_id":row.get("record_id"),"_malformed_CSV_row":row});continue
            identifier_kind=row.pop("__identifier_kind",None)
            mask=row.pop("__present_fields",None)
            if mask is not None:
                try:
                    present=json.loads(mask)
                    if not isinstance(present,list) or not all(isinstance(k,str) for k in present) or len(present)!=len(set(present)) or not set(present)<=set(row) or any(v not in (None,"") for k,v in row.items() if k not in present):raise ValueError("Invalid source field mask")
                    row={k:v for k,v in row.items() if k in present}
                except (ValueError,TypeError):
                    result.append({"record_id":row.get("record_id"),"_malformed_CSV_presence_mask":mask});continue
            if identifier_kind in {"integer","observation","unknown"} and "record_id" in row:
                try:
                    if identifier_kind=="unknown" and row["record_id"] not in (None,"","null"):raise ValueError("Nonempty unknown ID")
                    row["record_id"]=json.loads(row["record_id"]) if identifier_kind!="unknown" else None
                    if identifier_kind=="observation" and not isinstance(row["record_id"],dict):raise ValueError("Invalid identifier observation")
                    if identifier_kind=="integer" and (not isinstance(row["record_id"],int) or isinstance(row["record_id"],bool)):raise ValueError("Invalid integer identifier")
                except (ValueError,TypeError):
                    result.append({"record_id":None,"_invalid_identifier_encoding":row});continue
            elif identifier_kind not in {None,"string","omitted"}:
                result.append({"record_id":None,"_invalid_identifier_kind":identifier_kind});continue
            parsed={}
            for name,value in row.items():
                if name in {"record_id","schema_id"}:parsed[name]=value;continue
                if value is None or value.strip()=="":parsed[name]=None;continue
                if value.strip().lower() in {"true","false"}:
                    parsed[name]=value.strip().lower()=="true";continue
                try:parsed[name]=json.loads(value)
                except (ValueError,TypeError):parsed[name]=value
            result.append(parsed)
        return result


def read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        value=json.load(stream,parse_constant=lambda token: (_ for _ in ()).throw(ValidationError("JSON","Non-finite JSON constant: "+token)))
    return records(value)


def write_records(rows,*,csv_path,json_path):
    rows=records(rows)
    # Object-valued cells are explicit JSON; identifiers remain literal strings.
    names=["record_id"]+sorted(set().union(*(set(r) for r in rows))-{"record_id"}) if rows else ["record_id"]
    with Path(csv_path).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=names+["__present_fields","__identifier_kind"]);writer.writeheader()
        for row in rows:
            value={k:(v if k in {"record_id","schema_id"} else canonical(v)) for k,v in row.items()}
            identifier=row.get("record_id")
            kind="omitted" if "record_id" not in row else ("unknown" if identifier is None else "observation" if isinstance(identifier,dict) else "integer" if isinstance(identifier,int) and not isinstance(identifier,bool) else "string")
            value["__identifier_kind"]=kind
            if kind in {"integer","observation"}:value["record_id"]=canonical(identifier)
            value["__present_fields"]=canonical(sorted(row))
            writer.writerow(value)
    Path(json_path).write_text(canonical(rows)+"\n",encoding="utf-8")


def export_results(results,*,csv_path,json_path):
    rows=records(results)
    fields=("record_id","rank","formulation","parameter_package","source_version","central_absolute_risk","lower_missing_scenario","upper_missing_scenario",
            "tbi_probability","conditional_progression_risk","effective_rr","expected_effective_rr","combination_policy",
            "horizon_years","status","ranking_uncertain","overlapping_missing_information_ranges","warning_status",
            "principal_risk_drivers","missing_critical_variables","input_completeness","warnings","scientific_hash","full_result_json")
    with Path(csv_path).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for row in rows:
            value={k:row.get(k) for k in fields}
            for k,v in list(value.items()):
                if isinstance(v,(dict,list,tuple)):value[k]=canonical(v)
            value["full_result_json"]=canonical(row);writer.writerow(value)
    Path(json_path).write_text(canonical(rows)+"\n",encoding="utf-8")
