"""Bounded offline M1 validation; generated evidence stays in abm/output."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def git(path,*args):
    return subprocess.check_output(["git","-C",str(path),*args]).decode()


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def protection_check():
    base=json.loads((ROOT/"docs/individual_risk_model/protected_start.json").read_text(encoding="utf-8"))
    for name,value in base["tracked_files"].items():
        if sha(ROOT/name)!=value:raise AssertionError("Frozen existing file changed: "+name)
    refs=dict(x.split(" ",1) for x in git(ROOT,"for-each-ref","--format=%(refname) %(objectname)").splitlines())
    allowed={"refs/heads/feature/individual-tb-risk-platform-v1","refs/remotes/origin/feature/individual-tb-risk-platform-v1"}
    if {k:v for k,v in refs.items() if k not in allowed}!={k:v for k,v in base["refs"].items() if k not in allowed}:
        raise AssertionError("Protected refs changed")
    for name,expected in base["protected_worktrees"].items():
        p=Path(expected["path"])
        changed=set(git(p,"diff","HEAD","--name-only","-z").split("\0"))|set(git(p,"ls-files","--others","--exclude-standard","-z").split("\0"));changed.discard("")
        actual=dict(path=str(p),head=git(p,"rev-parse","HEAD").strip(),branch=git(p,"branch","--show-current").strip(),status=git(p,"status","--porcelain"),
            index_diff=hashlib.sha256(subprocess.check_output(["git","-C",str(p),"diff","--cached","--binary"])).hexdigest(),
            worktree_diff=hashlib.sha256(subprocess.check_output(["git","-C",str(p),"diff","--binary"])).hexdigest(),
            changed_files={n:sha(p/n) for n in sorted(changed) if (p/n).is_file()})
        if actual!=expected:raise AssertionError("Protected worktree changed: "+name)
    assert git(ROOT,"rev-parse","sa-health-apy-he-v1.0.0^{}").strip()==base["starting_head"]
    return dict(passed=True,existing_frozen_files=len(base["tracked_files"]),protected_worktrees=len(base["protected_worktrees"]),all_other_refs_unchanged=True)


def main():
    folder=ROOT/"abm/output/individual_risk_m1";folder.mkdir(parents=True,exist_ok=True)
    files=sorted(list((ROOT/"tb_individual_risk").glob("*.py"))+list((ROOT/"examples/individual_risk").glob("*.py"))+[ROOT/"tests/test_individual_risk_m1.py",Path(__file__)])
    before={str(p.relative_to(ROOT)).replace("\\","/"):sha(p) for p in files}
    for p in files:compile(p.read_bytes(),str(p),"exec")
    commands=[]
    for name,cmd in [("focused_tests",[sys.executable,"-m","pytest","tests/test_individual_risk_m1.py","-q","--durations=8"]),
                     ("synthetic_example",[sys.executable,"examples/individual_risk/run_example.py"])]:
        start=time.perf_counter();result=subprocess.run(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);elapsed=time.perf_counter()-start
        log=folder/(name+".txt");log.write_bytes(result.stdout)
        record=dict(command=subprocess.list2cmdline(cmd),exit_code=result.returncode,elapsed_seconds=elapsed,
                    log=str(log.relative_to(ROOT)).replace("\\","/"),log_sha256=sha(log),output=result.stdout.decode("utf-8",errors="replace"))
        commands.append(record)
        if result.returncode:print(record["output"]);raise SystemExit(result.returncode)
    assert before=={str(p.relative_to(ROOT)).replace("\\","/"):sha(p) for p in files}
    protection=protection_check()
    subprocess.run(["git","diff","--check"],cwd=ROOT,check=True)
    from tb_individual_risk import IndividualTBRiskModel,read_csv,read_json
    from tb_individual_risk.provenance import scientific_hash
    from examples.individual_risk.synthetic_parameters import synthetic_parameters,fictional_people
    model=IndividualTBRiskModel(synthetic_parameters());people=fictional_people()
    assert read_csv(folder/"fictional_people.csv")==read_json(folder/"fictional_people.json")==people
    ranked=json.loads((folder/"highest_applicable_rr.json").read_text(encoding="utf-8"))
    for row in ranked:
        content=dict(row);content.pop("scientific_hash")
        assert scientific_hash(content)==row["scientific_hash"]
    comparison={policy:[dict(record_id=r["record_id"],rank=r["rank"],absolute_risk=r["central_absolute_risk"],effective_rr=r.get("effective_rr")) for r in rows]
                for policy,rows in model.compare_policies(people).items()}
    receipt=dict(schema_id="individual_tb.m1_validation_receipt.v1",recorded_at_utc=datetime.now(timezone.utc).isoformat(),
        starting_commit="03cc16e4a52a55e10019dab16bd4f599571c4b90",commands=commands,compiled_files=len(files),source_checksums=before,
        source_files_unchanged_during_checks=True,protection=protection,csv_json_roundtrip=True,all_result_hashes_verified=True,
        example_people=len(people),example_resolved=sum(r["central_absolute_risk"] is not None for r in ranked),comparison=comparison,
        parameter_hash=model.parameters.record()["scientific_hash"],clinical_use_approved=False,
        limitations=["All example parameters invented; no clinical central parameter set","No external validation or TSTin3D reproduction claim","Population generator contract only","No new UI, interventions, economics or transmission"])
    (folder/"validation_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in receipt.items() if k not in {"source_checksums","comparison","commands"}}))
    for row in commands:print(row["command"],"elapsed_seconds=",row["elapsed_seconds"]);print(row["output"])

if __name__=="__main__":main()
