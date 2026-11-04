"""Bounded offline M2A tests, reproduction, compilation and protected-state audit."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
ALLOWED_M1_EDITS={"tb_individual_risk/__init__.py","tb_individual_risk/contracts.py","tb_individual_risk/model.py","tb_individual_risk/parameters.py","tb_individual_risk/provenance.py","tb_individual_risk/risk_factors.py","tb_individual_risk/ranking.py","tb_individual_risk/line_list.py","docs/individual_risk_model/README.md"}


def git(path,*args):return subprocess.check_output(["git","-C",str(path),*args])
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def protection_check():
    base=json.loads((ROOT/"docs/individual_risk_model/protected_start_m2a.json").read_text())
    frozen=json.loads((ROOT/"docs/individual_risk_model/protected_start.json").read_text())
    for name,h in frozen["tracked_files"].items():
        assert sha(ROOT/name)==h,"Frozen release file changed: "+name
    for name,h in base["tracked_files"].items():
        if name not in ALLOWED_M1_EDITS:assert sha(ROOT/name)==h,"Unexpected retained-file change: "+name
    refs=dict(line.split(" ",1) for line in git(ROOT,"for-each-ref","--format=%(refname) %(objectname)").decode().splitlines())
    allowed={"refs/heads/feature/individual-tb-risk-platform-v1","refs/remotes/origin/feature/individual-tb-risk-platform-v1"}
    assert {k:v for k,v in refs.items() if k not in allowed}=={k:v for k,v in base["refs"].items() if k not in allowed},"Protected refs changed"
    for name,expected in base["protected_worktrees"].items():
        p=Path(expected["path"])
        changed=set(git(p,"diff","HEAD","--name-only","-z").decode().split("\0"))|set(git(p,"ls-files","--others","--exclude-standard","-z").decode().split("\0"));changed.discard("")
        actual=dict(path=str(p),head=git(p,"rev-parse","HEAD").decode().strip(),branch=git(p,"branch","--show-current").decode().strip(),status=git(p,"status","--porcelain").decode(),index_diff=hashlib.sha256(git(p,"diff","--cached","--binary")).hexdigest(),worktree_diff=hashlib.sha256(git(p,"diff","--binary")).hexdigest(),changed_files={n:sha(p/n) for n in sorted(changed) if (p/n).is_file()})
        assert actual==expected,"Protected worktree changed: "+name
    return dict(passed=True,frozen_release_files=len(frozen["tracked_files"]),protected_worktrees=len(base["protected_worktrees"]),retained_files_except_authorized_edits=len(base["tracked_files"])-len(ALLOWED_M1_EDITS),frozen_peeled_tag=git(ROOT,"rev-parse","sa-health-apy-he-v1.0.0^{}").decode().strip())


def main():
    output=ROOT/"abm/output/individual_risk_m2a";output.mkdir(parents=True,exist_ok=True)
    files=set((ROOT/"tb_individual_risk").rglob("*.py"))|set((ROOT/"tb_individual_risk/reference_data").glob("*.json"))|set((ROOT/"examples/individual_risk").glob("*.py"))|set((ROOT/"tests").glob("test_individual_risk_m*.py"))|{Path(__file__),ROOT/"scripts/reproduce_menzies_2008.py"}
    before={str(p.relative_to(ROOT)).replace("\\","/"):sha(p) for p in sorted(files)}
    python_files=[p for p in files if p.suffix==".py"]
    for p in python_files:compile(p.read_bytes(),str(p),"exec")
    commands=[]
    for name,cmd in [("focused_tests",[sys.executable,"-m","pytest","tests/test_individual_risk_m1.py","tests/test_individual_risk_m2a_audit.py","tests/test_individual_risk_m2a_reference.py","-q","--durations=8"]),("reproduction",[sys.executable,"scripts/reproduce_menzies_2008.py"])]:
        start=time.perf_counter();r=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);elapsed=time.perf_counter()-start
        log=output/(name+".txt");log.write_bytes(r.stdout);text=r.stdout.decode("utf-8",errors="replace")
        commands.append(dict(command=subprocess.list2cmdline(cmd),exit_code=r.returncode,elapsed_seconds=elapsed,output=text,log_sha256=sha(log)))
        if r.returncode:print(text);raise SystemExit(r.returncode)
    assert before=={str(p.relative_to(ROOT)).replace("\\","/"):sha(p) for p in sorted(files)}
    from tb_individual_risk.packages import PINS,load_parameter_package
    from tb_individual_risk.provenance import scientific_hash
    from tb_individual_risk.reference_comparison import reproduce_published_examples
    package=load_parameter_package("menzies_2008_published_reference_v1")
    report=reproduce_published_examples(package)
    assert report==reproduce_published_examples(package)
    hashes=set();occurrences=0
    def verify(value):
        nonlocal occurrences
        if isinstance(value,dict):
            if "scientific_hash" in value:
                payload=dict(value);h=payload.pop("scientific_hash");assert scientific_hash(payload)==h,value.get("schema_id");hashes.add(h);occurrences+=1
            for v in value.values():verify(v)
        elif isinstance(value,list):
            for v in value:verify(v)
    for p in (output/"reproduction").glob("*.json"):verify(json.loads(p.read_text(encoding="utf-8")))
    protection=protection_check();subprocess.run(["git","diff","--check"],cwd=ROOT,check=True)
    receipt=dict(schema_id="individual_tb.m2a_validation_receipt.v1",recorded_at_utc=datetime.now(timezone.utc).isoformat(),starting_head="8b533f81a6a1435bed11f10db52c6d8d916224a1",commands=commands,compiled_python_files=len(python_files),source_checksums=before,source_files_unchanged=True,package_file_checksums=PINS,scientific_hash_occurrences=occurrences,unique_scientific_hashes=len(hashes),protection=protection,reproduction_summary={k:report[k] for k in ("full_algorithm_matches","partial_comparisons","partial_within_printed_resolution")},clinical_validation=False,ordinary_prediction_network=False)
    (output/"validation_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k not in {"source_checksums","commands"}}))
    for command in commands:print(command["command"],"elapsed_seconds=",command["elapsed_seconds"]);print(command["output"])

if __name__=="__main__":main()
