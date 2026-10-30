"""Offline schema/hash checks for a completed empirical diagnostic export."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from jsonschema import Draft202012Validator
from engine.country_tb.data_common import REPO_ROOT,profile_hash,digest
from engine.country_tb.validation_data import require,load_validation_snapshot,load_mortality_validation_snapshot


def check(paths):
    s=json.loads((REPO_ROOT/'schemas/country_tb/v1/empirical_validation.schema.json').read_text(encoding='utf-8'))
    Draft202012Validator.check_schema(s);families={d['properties']['schema_id']['const']:name for name,d in s['$defs'].items()
        if 'properties' in d and 'schema_id' in d['properties']}
    counts={k:0 for k in families};hashes=set();occurrences=0;checksums={}
    def visit(v):
        nonlocal occurrences
        if isinstance(v,dict):
            if 'profile_hash' in v:
                require(v['profile_hash']==profile_hash(v),'Scientific hash mismatch in export.')
                hashes.add(v['profile_hash']);occurrences+=1
            family=v.get('schema_id')
            if family in families:
                Draft202012Validator({'$ref':'#/$defs/'+families[family],'$defs':s['$defs']}).validate(v)
                counts[family]+=1
            for a in v.values():visit(a)
        elif isinstance(v,list):
            for a in v:visit(a)
    for p in paths:
        p=Path(p);b=p.read_bytes();checksums[str(p.relative_to(REPO_ROOT))]=digest(b);visit(json.loads(b.decode('utf-8')))
    return {'passed':True,'family_occurrences':counts,'unique_scientific_hashes':len(hashes),
        'hash_occurrences':occurrences,'file_checksums':checksums,'network_access':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('files',nargs='+',type=Path);a=p.parse_args()
    print(json.dumps(check([f.resolve() for f in a.files]),sort_keys=True))


if __name__=='__main__':main()
