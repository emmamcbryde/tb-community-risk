"""Focused Milestone 2A data and forensic checks; no epidemiological execution."""
from copy import deepcopy
from decimal import Decimal
import ast
import json
from pathlib import Path
import socket
import subprocess

import pytest

from engine.country_tb.data_common import (
    DataValidationError, REPO_ROOT, canonical, csv_bytes, decimal, digest,
    profile_hash, publish_immutable,
)
from engine.country_tb.incidence import (
    WHO_2025, COLUMNS, _manifest_hash, load_who_snapshot, country_index, national_profile,
    scaled_profile, user_series_profile, compare_snapshots, RATE_UNIT,
)
from engine.country_tb.demography import (
    DEMOGRAPHIC_KEY, normalize_records, canonical_sexes, combined_rows,
    aggregate_display, demographic_profile, resolve_source,
)
from engine.country_tb.duplicate_audit import audit_legacy
from engine.country_tb.data_profile import load_candidates, baseline_data_profile, preview_files


@pytest.fixture(scope="session")
def who():
    return load_who_snapshot(WHO_2025)


@pytest.fixture(scope="session")
def demographic_sources():
    return load_candidates()


def record(sex="female", age=0, upper=None, population="10", **updates):
    return {"source_id": "test", "source_vintage": "v1", "geography_id": "who:VNM",
            "geography_level": "national", "year": 2024, "estimate_status": "estimate",
            "projection_scenario": "historical", "sex": sex, "age_lower": age,
            "age_upper_exclusive": upper, "unit": "persons", "population_variant": "resident",
            "population": population, **updates}


def candidate(category="un_wpp", country="VNM", records=None, identifier=None):
    rows = records if records is not None else [record("female", population="12"), record("male", population="8")]
    rows = [{**r, "geography_id": "who:" + country} for r in rows]
    return {"snapshot_id": identifier or category, "rows": rows,
            "input_source_hash": "0"*64, "transformed_file_hash": "1"*64,
            "manifest": {"source_id": category, "source_category": category, "source_vintage": "v1",
                "review_status": "reviewed", "manifest_hash": "2"*64, "schema_version": "country_tb.demography_snapshot.v1",
                "transformation_version": "country_tb_demography_import_1.0.0", "transformations": ["test input"],
                "source_url": "urn:test", "import_date": "2026-10-27", "terms": {"status": "test"},
                "citation": "Test only", "projection_boundary": None}}


def who_bundle(tmp_path, who, rows=None, mutation=None, missing_column=False):
    manifest = deepcopy(who["manifest"])
    source_rows = deepcopy(rows if rows is not None else [r for r in who["rows"] if r["iso3"] == "AUS"])
    columns = COLUMNS[:-1] if missing_column else COLUMNS
    payload = csv_bytes([{c:r[c] for c in columns} for r in source_rows], columns)
    identifier = "who-gtb2025-incidence-" + digest(payload)[:12]
    manifest.update(snapshotId=identifier, subsetIso3=None)
    manifest["dataFile"] = {"filename": identifier+".csv", "sha256": digest(payload), "rows": len(source_rows)}
    manifest["coverage"] = {"countriesAndAreas":len({r["iso3"] for r in source_rows}), "withEstimates":1,
                            "yearRange":[min(int(r["year"]) for r in source_rows),max(int(r["year"]) for r in source_rows)]}
    if mutation:
        mutation(manifest)
    manifest["analyticalHash"] = _manifest_hash(manifest)
    (tmp_path/(identifier+".csv")).write_bytes(payload)
    (tmp_path/(identifier+"_manifest.json")).write_bytes(canonical(manifest))
    return identifier


def test_original_1212_decomposition_and_preserved_source():
    source = (REPO_ROOT/"data/population_age_latest.csv").read_bytes()
    assert digest(source) == "e51101fa3f5054b9aeb11999d38f6b39a3b0ab4625daa00c5c7fd52e77e535d4"
    audit, summary, repaired = audit_legacy(source)
    assert (summary["source_rows"],summary["unique_legacy_keys"],summary["duplicate_groups"],
            summary["involved_rows"],summary["excess_rows"]) == (21917,20705,1111,2323,1212)
    assert len(summary["per_country_excess"]) == 11
    assert summary["per_country_excess"]["VCT"] == 202
    assert all(v == 101 for k,v in summary["per_country_excess"].items() if k != "VCT")
    assert summary["categories"][0]["excess_rows"] == 1212
    assert all(c["excess_rows"] == 0 for c in summary["categories"][1:])
    assert all(r["max_population_difference"] == "0" for r in audit)
    assert repaired == (REPO_ROOT/"data/demography/legacy_audit/legacy_unique_quarantine.csv").read_bytes()
    saved = json.loads((REPO_ROOT/"data/demography/legacy_audit/duplicate_summary.json").read_text())
    assert summary == saved


def test_legacy_resolution_does_not_depend_on_file_order():
    rows = [{"location":"Test","iso_code":"TST","year":"2021","age":str(a),"population":str(a+1)}
            for a in range(3) for _ in range(2)]
    columns=tuple(rows[0])
    a=audit_legacy(csv_bytes(rows,columns))[2]
    rows.reverse()
    assert a == audit_legacy(csv_bytes(rows,columns))[2]
    rows[0]["population"]="999"
    with pytest.raises(DataValidationError):
        audit_legacy(csv_bytes(rows,columns))


def test_complete_key_and_exact_duplicate_removal():
    assert DEMOGRAPHIC_KEY == ("source_id","source_vintage","geography_id","geography_level","year",
        "estimate_status","projection_scenario","sex","age_lower","age_upper_exclusive","unit","population_variant")
    rows,count=normalize_records([record(),record()])
    assert len(rows)==1 and count==1
    with pytest.raises(DataValidationError,match="Conflicting"):
        normalize_records([record(),record(population="11")])
    a=record(source_checksum="0"*64);b=record(source_checksum="1"*64)
    with pytest.raises(DataValidationError):
        normalize_records([a,b])


def test_sex_omission_explains_toy_incomplete_keys_without_losing_records():
    rows=[record("female",population="12"),record("male",population="8"),record("both",population="20")]
    assert len({(r["geography_id"],r["year"],r["age_lower"]) for r in rows})==1
    normalized,count=normalize_records(rows)
    assert len(normalized)==3 and count==0
    result=canonical_sexes(rows,absolute_tolerance="0",relative_tolerance="0")
    assert {r["sex"] for r in result["rows"]}=={"female","male"}
    assert sum(decimal(r["population"]) for r in result["rows"])==20
    assert len(result["raw_records"])==3
    assert combined_rows(result["rows"])[0]["population"]=="20"
    assert result["published_both_checks"][0]["difference"]=="0"
    with pytest.raises(DataValidationError):
        combined_rows(rows)


def test_both_sex_reconciliation_tolerance_and_conflict():
    rows=[record("female",population="12"),record("male",population="8"),record("both",population="19")]
    assert canonical_sexes(rows)["published_both_checks"][0]["difference"]=="1"
    rows[-1]["population"]="18"
    with pytest.raises(DataValidationError):
        canonical_sexes(rows)


def test_both_only_and_unknown_are_visible_and_never_split():
    result=canonical_sexes([record("both",population="20")])
    assert result["rows"][0]["sex"]=="both" and not result["sex_specific_available"]
    assert "sex_specific_modelling_unavailable" in result["warnings"]
    assert "sex_unknown" in canonical_sexes([record("unknown")])["warnings"]
    with pytest.raises(DataValidationError):
        canonical_sexes([record("unknown"),record("male")])


@pytest.mark.parametrize("field,value",[("source_vintage","v2"),("estimate_status","projection"),
                                      ("projection_scenario","high"),("year",2025)])
def test_other_key_dimensions_are_preserved(field,value):
    a=record();b=record(**{field:value})
    if field in {"estimate_status","projection_scenario"}:
        b.update(estimate_status="projection",projection_scenario="high")
    result,count=normalize_records([a,b])
    assert len(result)==2 and count==0


def test_overlapping_age_bands_fail_and_missing_age_warns():
    with pytest.raises(DataValidationError,match="Overlapping"):
        normalize_records([record(age=0,upper=5),record(age=4,upper=None)])
    result=canonical_sexes([record("both",age=1,upper=None)])
    assert "incomplete_age_coverage" in result["warnings"]
    assert result["missing_age_intervals"]==[{"sex":"both","age_lower":0,"age_upper_exclusive":1}]


def test_single_ages_and_open_tail_aggregate_exactly():
    rows=[record("female",age=a,upper=a+1,population="1") for a in range(100)]
    rows.append(record("female",age=100,upper=None,population="2"))
    result=aggregate_display(rows)
    assert [r["population"] for r in result]==["5","10","30","20","37"]
    assert sum(decimal(r["population"]) for r in result)==102
    assert rows[-1]["age_upper_exclusive"] is None
    with pytest.raises(DataValidationError):
        aggregate_display([record("both",age=0,upper=None)])


@pytest.mark.parametrize("bad",["-1","NaN","Infinity",True,None])
def test_invalid_population_fails(bad):
    with pytest.raises(DataValidationError):
        normalize_records([record(population=bad)])


def test_row_order_and_display_are_not_scientific_changes():
    rows=[record("female"),record("male")]
    a=demographic_profile("who:VNM",[candidate(records=rows)])
    b=demographic_profile("who:VNM",[candidate(records=rows[::-1])])
    assert a["profile_hash"]==b["profile_hash"]
    b["display"]={"name":"A changed translation"}
    assert profile_hash(b)==a["profile_hash"]
    for field,value in [("population","11"),("sex","unknown"),("age_lower",1),("year",2025),
                         ("source_vintage","v2"),("projection_scenario","high")]:
        changed=deepcopy(a);changed["rows"][0][field]=value
        assert profile_hash(changed)!=a["profile_hash"],field
    restored=json.loads(canonical(a))
    assert restored["provenance"]==a["provenance"] and profile_hash(restored)==a["profile_hash"]


def test_source_hierarchies_and_review_and_ties():
    aus=[candidate("un_wpp","AUS"),candidate("abs","AUS"),candidate("user","AUS")]
    assert resolve_source("who:AUS",aus)[0]["manifest"]["source_category"]=="user"
    assert resolve_source("who:AUS",aus[:-1])[0]["manifest"]["source_category"]=="abs"
    vnm=[candidate("fallback"),candidate("un_wpp"),candidate("user")]
    assert resolve_source("who:VNM",vnm[:-1])[0]["manifest"]["source_category"]=="un_wpp"
    assert resolve_source("who:VNM",vnm)[0]["manifest"]["source_category"]=="user"
    vnm[-1]["manifest"]["review_status"]="not_reviewed"
    assert resolve_source("who:VNM",vnm)[0]["manifest"]["source_category"]=="un_wpp"
    with pytest.raises(DataValidationError,match="Equally ranked"):
        resolve_source("who:VNM",[candidate(identifier="one"),candidate(identifier="two")])
    assert resolve_source("who:VNM",[candidate(identifier="one"),candidate(identifier="two")],
                          preferred_snapshot_id="two")[0]["snapshot_id"]=="two"


def test_baseline_year_exact_or_explicit_reviewed_proxy():
    source=candidate(records=[record("both",year=2023)])
    with pytest.raises(DataValidationError):
        demographic_profile("who:VNM",[source])
    rule={"method":"nearest","review_status":"reviewed","rationale":"Explicit test proxy","maximum_year_displacement":1}
    result=demographic_profile("who:VNM",[source],year_rule=rule)
    assert result["baseline_year"]==2024 and result["rows"][0]["year"]==2023
    assert result["source_resolution"]["baseline_resolution"]["year_displacement"]==-1
    assert "baseline_uses_explicit_nearest_year_proxy" in result["warnings"]
    with pytest.raises(DataValidationError):
        demographic_profile("who:VNM",[source],year_rule={**rule,"review_status":"not_reviewed"})


def test_medium_default_and_alternatives_preserved():
    rows=[record("both",year=2025,estimate_status="projection",projection_scenario=scenario,population=count)
          for scenario,count in (("medium","20"),("high","30"),("low","10"))]
    result=demographic_profile("who:VNM",[candidate(records=rows)],year=2025)
    assert result["combined_population"]=="20"
    assert result["inactive_alternatives"][0]["scenarios"]==["high","low","medium"]
    assert demographic_profile("who:VNM",[candidate(records=rows)],year=2025,scenario="high")["combined_population"]=="30"


def test_immutable_atomic_import_is_idempotent_and_never_appends(tmp_path):
    files={"population.csv":b"header\nvalue\n","manifest.json":b'{"schema":"test"}'}
    target=publish_immutable(tmp_path,"test-snapshot",files)
    before={p.name:p.read_bytes() for p in target.iterdir()}
    assert publish_immutable(tmp_path,"test-snapshot",files)==target
    assert {p.name:p.read_bytes() for p in target.iterdir()}==before
    with pytest.raises(DataValidationError):
        publish_immutable(tmp_path,"test-snapshot",{"population.csv":b"different"})
    assert before["population.csv"]==(target/"population.csv").read_bytes()


def test_pinned_who_and_country_selector(who):
    assert (len(who["rows"]),who["report_year"],who["final_observation_year"])==(5347,2025,2024)
    index=country_index(who)
    assert len(index)==217
    assert next(r for r in index if r["country_id"]=="PRK")["completeness"]=="unavailable"
    assert next(r for r in index if r["country_id"]=="ANT")["available_incidence_years"]==list(range(2000,2010))
    ids={r["country_id"] for r in index}
    assert {"ANT","CUW","SXM","GBR","BMU"}<=ids
    assert country_index({**who,"rows":list(reversed(who["rows"]))})==index


def test_who_missing_manifest_fixture_and_checksum_failure(tmp_path,who):
    with pytest.raises(DataValidationError):
        load_who_snapshot(WHO_2025,tmp_path)
    with pytest.raises(DataValidationError):
        load_who_snapshot("who-gtb2025-incidence-fixture-1709dd393fec",tmp_path)
    identifier=who_bundle(tmp_path,who)
    path=tmp_path/(identifier+".csv")
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(DataValidationError,match="checksum"):
        load_who_snapshot(identifier,tmp_path)


@pytest.mark.parametrize("mutation",[
    lambda m:m.update(contractVersion="unknown"),
    lambda m:m.update(importerVersion="unknown"),
    lambda m:m.update(snapshotKind="test_fixture"),
    lambda m:m.update(isCompleteDataset=False),
    lambda m:m["coverage"].update(yearRange=[2000,2025]),
])
def test_who_manifest_semantics(tmp_path,who,mutation):
    identifier=who_bundle(tmp_path,who,mutation=mutation)
    with pytest.raises(DataValidationError):
        load_who_snapshot(identifier,tmp_path)


def test_who_manifest_tampering_and_columns(tmp_path,who):
    identifier=who_bundle(tmp_path,who,missing_column=True)
    with pytest.raises(DataValidationError,match="columns"):
        load_who_snapshot(identifier,tmp_path)
    path=tmp_path/(identifier+"_manifest.json")
    m=json.loads(path.read_text());m["terms"]["status"]="tampered";path.write_bytes(canonical(m))
    with pytest.raises(DataValidationError,match="Manifest analytical"):
        load_who_snapshot(identifier,tmp_path)


@pytest.mark.parametrize("field,bad",[("incidence_per_100k","-1"),("incidence_per_100k","NaN"),
    ("incidence_per_100k_hi","Infinity"),("incidence_per_100k_lo","999999"),("population","0"),("iso3","bad")])
def test_who_invalid_values(tmp_path,who,field,bad):
    rows=deepcopy([r for r in who["rows"] if r["iso3"]=="AUS"])
    rows[-1][field]=bad
    identifier=who_bundle(tmp_path,who,rows)
    with pytest.raises(DataValidationError):
        load_who_snapshot(identifier,tmp_path)


def test_who_duplicate_country_year(tmp_path,who):
    rows=deepcopy([r for r in who["rows"] if r["iso3"]=="AUS"])
    rows.append(deepcopy(rows[-1]))
    identifier=who_bundle(tmp_path,who,rows)
    with pytest.raises(DataValidationError,match="Duplicate"):
        load_who_snapshot(identifier,tmp_path)


def test_no_network_during_load(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError("Network forbidden")
    monkeypatch.setattr(socket,"socket",forbidden)
    assert load_who_snapshot(WHO_2025)["final_observation_year"]==2024


def test_incidence_scaling_preserves_parent_and_hashes(who):
    parent=national_profile(who,"VNM");before=deepcopy(parent)
    result=scaled_profile(parent,"iso3166-2:VN-01",2,rationale="Test scale",provenance={"source":"test"},review_status="not_reviewed")
    assert result["mode"]=="subnational_scaled"
    assert decimal(result["rows"][-1]["central"])==2*decimal(parent["rows"][-1]["central"])
    assert parent==before and result["national_parent"]==parent
    assert result["profile_hash"]!=parent["profile_hash"]
    assert json.loads(canonical(result))["transformation"]==result["transformation"]
    with pytest.raises(DataValidationError):
        scaled_profile(parent,"iso3166-2:VN-01",{2024:2},rationale="Test",provenance={"source":"x"},review_status="reviewed")
    changed=deepcopy(parent);changed["rows"][-1]["upper"]="999"
    assert profile_hash(changed)!=parent["profile_hash"]


def test_user_series_preserves_notifications_numerators_and_gaps():
    rows=[{"year":year,"original_rate":"10","numerator":"2","denominator":"20000",
           "lower":"8","upper":"12","unit":RATE_UNIT,"quantity":"notifications"} for year in (2022,2024)]
    profile=user_series_profile("subnational_user_series","custom:test",rows,provenance={"source":"test"},review_status="reviewed")
    assert profile["rows"][0]["quantity"]=="notifications" and profile["missing_years"]==[2023]
    assert profile["rows"][0]["numerator"]=="2"
    conversion={"method":"ascertainment_fraction","fraction":".5","review_status":"reviewed","source":"test"}
    adjusted=user_series_profile("custom_population_series","custom:test",rows,provenance={"source":"test"},
                                 review_status="reviewed",conversion=conversion)
    assert adjusted["rows"][0]["central"]=="20" and adjusted["rows"][0]["original_rate"]=="10"
    assert adjusted["rows"][0]["quantity"]=="estimated_tb_incidence"
    with pytest.raises(DataValidationError):
        user_series_profile("custom_population_series","custom:test",rows+rows[:1],provenance={"source":"test"},review_status="reviewed")


def test_snapshot_comparison_is_pure_and_reports_bounds_and_country_changes(who):
    before=deepcopy(who)
    newer=deepcopy(who);newer["snapshot_id"]="test-new"
    newer["rows"]=[r for r in newer["rows"] if r["iso3"]!="ANT"]
    newer["rows"][0]["incidence_per_100k_hi"]="999"
    result=compare_snapshots(who,newer)
    assert result["countries_removed"]==["ANT"]
    assert result["historical_revisions"][0]["changes"]["incidence_per_100k_hi"]["new"]=="999"
    assert who==before


def test_real_baselines_sexes_sources_counts_and_warnings(who,demographic_sources):
    expected={"AUS":("abs_erp_jun2024","27204809"),"VNM":("un_wpp2024_medium","100987682"),
              "PRK":("un_wpp2024_medium","26498835")}
    for country,(source,total) in expected.items():
        profile=baseline_data_profile(country,who,demographic_sources);d=profile["demography"]
        assert d["provenance"]["source_id"]==source and d["combined_population"]==total
        assert {r["sex"] for r in d["rows"]}=={"female","male"}
        assert all(r["year"]==2024 for r in d["rows"])
        assert len(d["rows"])==202 and len(d["combined_rows"])==101 and len(d["display_bands"])==10
        assert sum(decimal(r["population"]) for r in d["display_bands"])==decimal(total)
        assert max(r["age_lower"] for r in d["rows"] if r["age_upper_exclusive"] is None)==100
        files=preview_files(profile)
        assert {"incidence.csv","population_by_sex_age.csv","combined_rows.csv","display_bands.csv","profile.json","README.md"}<=files.keys()
        if country!="AUS":
            assert all(r["estimate_status"]=="projection" for r in d["rows"])
            assert "baseline_2024_is_source_projection" in d["warnings"]
        if country=="PRK":
            assert "no_published_incidence_estimates" in profile["warnings"]
    assert baseline_data_profile("VNM",who,demographic_sources)["profile_hash"]==baseline_data_profile("VNM",who,list(reversed(demographic_sources)))["profile_hash"]


def test_actual_medium_2025_and_no_unavailable_extrapolation(demographic_sources):
    p=demographic_profile("who:VNM",demographic_sources,year=2025)
    assert all(r["projection_scenario"]=="medium" and r["year"]==2025 for r in p["rows"])
    with pytest.raises(DataValidationError):
        demographic_profile("who:VNM",demographic_sources,year=2031)


def test_no_deployed_imports_or_prohibited_execution_dependencies():
    for directory in ("app","ui","general_pages","general_ui","explicit_tbi_ui"):
        folder=REPO_ROOT/directory
        if folder.exists():
            for path in folder.rglob("*.py"):
                assert "engine.country_tb" not in path.read_text(encoding="utf-8-sig"),path
    banned=("streamlit","engine.apy","engine.static","engine.dynamic","engine.starsim","engine.profiles.engine_mapping")
    for path in (REPO_ROOT/"engine/country_tb").glob("*.py"):
        tree=ast.parse(path.read_text())
        for node in ast.walk(tree):
            imports=[a.name for a in node.names] if isinstance(node,ast.Import) else [node.module or ""] if isinstance(node,ast.ImportFrom) else []
            assert not any(name.startswith(banned) for name in imports),path


def test_protected_fingerprints_and_frozen_artifacts():
    receipt=REPO_ROOT/"abm/output/country_tb_m2a_audit/protected_start.json"
    if not receipt.is_file():
        pytest.skip("Session protection receipt not present; source/content isolation tests still run.")
    saved=json.loads(receipt.read_text())
    for name in ("gtbreport2025","tb-community-risk","tb-community-risk-no-transmission","tb-community-risk-tbi-calibration"):
        directory=REPO_ROOT.parent/name
        def git(*args):return subprocess.check_output(["git","-C",str(directory),*args],stderr=subprocess.DEVNULL)
        assert git("rev-parse","HEAD").decode().strip()==saved[name]["head"]
        assert digest(git("diff","HEAD","--binary"))==saved[name]["tracked_diff_sha256"]
        untracked={s:digest((directory/s).read_bytes()) for s in git("ls-files","--others","--exclude-standard").decode().splitlines() if (directory/s).is_file()}
        assert digest(json.dumps(untracked,sort_keys=True).encode())==saved[name]["untracked_hash"]
    changed=subprocess.check_output(["git","diff","--name-only","694851fda0394affef6731c44ff215302e83cfc7"]).decode().splitlines()
    assert not any(p.startswith(("abm/","validation/","reports/","app/","ui/","engine/apy/")) for p in changed)


def test_offline_future_who_publish_preserves_prior_and_is_idempotent(tmp_path,who):
    from engine.country_tb.incidence import freeze_who_bundle
    identifier=who_bundle(tmp_path,who)
    manifest=json.loads((tmp_path/(identifier+"_manifest.json")).read_text())
    payload=(tmp_path/(identifier+".csv")).read_bytes()
    registry=tmp_path/"registry"
    first=freeze_who_bundle(manifest,payload,registry)
    assert freeze_who_bundle(manifest,payload,registry)==first
    assert load_who_snapshot(identifier,registry)["source_data_hash"]==digest(payload)
    assert (tmp_path/(identifier+".csv")).read_bytes()==payload


def test_user_demography_import_idempotence_vintage_and_both_only(tmp_path):
    from engine.country_tb.demography_import import freeze_user_demography
    from engine.country_tb.demography import CSV_COLUMNS,load_demography_snapshot
    rows=[{k:v for k,v in record("both").items() if k in CSV_COLUMNS}]
    rows[0]["age_upper_exclusive"]=""
    source=tmp_path/"user.csv";source.write_bytes(csv_bytes(rows,CSV_COLUMNS))
    metadata={"source_id":"user_stats","source_vintage":"v1","source_url":"urn:test:authority",
              "terms":{"status":"user_permission"},"citation":"Authoritative supplied test counts",
              "review_status":"reviewed"}
    out=tmp_path/"snapshots"
    first=freeze_user_demography(source,metadata,out,import_date="2026-10-27")
    assert freeze_user_demography(source,metadata,out,import_date="2026-10-27")==first
    second=freeze_user_demography(source,{**metadata,"source_vintage":"v2"},out,import_date="2026-10-27")
    assert second!=first and first.is_dir()
    loaded=load_demography_snapshot(first.name,out)
    p=demographic_profile("who:VNM",[loaded])
    assert not p["sex_specific_available"] and p["combined_population"]=="10"


def test_demography_checksum_and_metadata_tampering(tmp_path):
    from engine.country_tb.demography import load_demography_snapshot
    import shutil
    source=REPO_ROOT/"data/demography/snapshots/abs-erp-jun2024-07ca11bcb092"
    target=tmp_path/source.name;shutil.copytree(source,target)
    data=target/"population.csv.gz";data.write_bytes(data.read_bytes()+b"tamper")
    with pytest.raises(DataValidationError,match="checksum"):
        load_demography_snapshot(source.name,tmp_path)
    data.write_bytes((source/"population.csv.gz").read_bytes())
    manifest=target/"manifest.json";m=json.loads(manifest.read_text());m["source_vintage"]="tampered"
    manifest.write_bytes(canonical(m))
    with pytest.raises(DataValidationError,match="manifest checksum"):
        load_demography_snapshot(source.name,tmp_path)


def test_demography_load_is_offline(monkeypatch):
    from engine.country_tb.demography import load_demography_snapshot
    def forbidden(*args,**kwargs):raise AssertionError("Network forbidden")
    monkeypatch.setattr(socket,"socket",forbidden)
    assert len(load_demography_snapshot("abs-erp-jun2024-07ca11bcb092")["rows"])==303


def test_executable_schema_and_restore_original_snapshot(who,demographic_sources):
    from jsonschema import Draft202012Validator
    from engine.country_tb.data_profile import restore_baseline_data
    schema=json.loads((REPO_ROOT/"schemas/country_tb/v1/data_layer.schema.json").read_text())
    Draft202012Validator.check_schema(schema)
    validator=Draft202012Validator(schema)
    p=baseline_data_profile("AUS",who,demographic_sources)
    for obj in (p,p["incidence"],p["demography"]):
        validator.validate(obj)
    restored=restore_baseline_data(canonical(p))
    assert restored==p
    p["demography"]["rows"][0]["population"]="999"
    with pytest.raises(DataValidationError,match="checksum"):
        restore_baseline_data(canonical(p))
