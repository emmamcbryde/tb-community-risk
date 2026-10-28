"""Offline developer preview: synthetic calibration and unreviewed WHO shapes.

Country alpha is an explicitly invented supplied-scale diagnostic. Synthetic
TBI anchors and assay parameters are software-validation data, not evidence.
Frozen country populations retain their source identity. No production imports.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.country_tb.calibration_common import EXCLUSIONS, artifact
from engine.country_tb.catalytic_calibration import calibration_policy, calibrate_scale
from engine.country_tb.catalytic_core import evaluate_age
from engine.country_tb.data_profile import load_candidates
from engine.country_tb.data_common import (
    DataValidationError, canonical, csv_bytes, digest, publish_immutable,
)
from engine.country_tb.demography import demographic_profile
from engine.country_tb.incidence import WHO_2025, load_who_snapshot, national_profile, user_series_profile
from engine.country_tb.pressure_shape import apply_backcast, shape_from_incidence
from engine.country_tb.tbi_baseline import baseline_tbi_profile, preview_files, uncertainty_scenario_set
from engine.country_tb.tbi_observations import assay_operator, observations_profile

REVIEW = {"status": "reviewed", "source_id": "synthetic:M2C",
          "rationale": "Invented software-validation fixture; not an empirical scientific review."}


def synthetic_example():
    geo = "custom:synthetic-m2c"
    with localcontext() as context:
        context.prec = 80
        incidence = [{"year": year, "original_rate": str(100 * (Decimal("-.015") * (year-2000)).exp()),
            "lower": str(50 * (Decimal("-.015") * (year-2000)).exp()),
            "upper": str(150 * (Decimal("-.015") * (year-2000)).exp()),
            "quantity": "estimated_tb_incidence", "unit": "cases_per_100000_person_year"} for year in range(2000, 2025)]
    source = user_series_profile("custom_population_series", geo, incidence,
                                provenance={"source_id": "synthetic:M2C"}, review_status="reviewed")
    shape = shape_from_incidence(source, review=REVIEW)
    history = apply_backcast(shape, start_time="1900")
    assay = assay_operator(assay_class="synthetic_identity", name="Invented perfect observation",
        positivity_threshold="synthetic_only", sensitivity=1, specificity=1,
        review=REVIEW, provenance={"source_id": "synthetic:M2C"})
    records = []
    for age, sex in ((5, "female"), (20, "male"), (65, "female")):
        prevalence = evaluate_age(history, alpha=".013", age=age, sex=sex)["any_tbi"]
        records.append({"observation_id": f"synthetic:{sex}:{age}", "age_lower": str(age),
            "age_upper_exclusive": str(age+1), "sex": sex, "geography_id": geo, "geography_level": "custom",
            "observation_policy": "diagnostic_synthetic_anchor", "survey_year": 2024, "survey_time": "2024", "prevalence": prevalence,
            "numerator": None, "denominator": None, "sampling_frame": "invented exact-age anchors",
            "survey_design": {"kind": "synthetic_only"}, "source": {"source_id": "synthetic:M2C"},
            "review": REVIEW, "assay_operator": assay, "uncertainty": {"kind": "exact"},
            "weight_source": {"source_id": "synthetic:M2C"},
            "strata": [{"age_lower": str(age), "age_upper_exclusive": None, "sex": sex,
                        "age_policy": "exact_age", "weight": "1"}]})
    observations = observations_profile(records, review=REVIEW)
    policy = calibration_policy(alpha_upper="1", observation_weights={r["observation_id"]: "1" for r in records},
        residual_tolerance=".000001", review=REVIEW, method="weighted_least_squares")
    calibration = calibrate_scale(history, policy, observations=observations)
    rows = [{"source_id": "synthetic", "source_vintage": "M2C", "geography_id": geo, "geography_level": "custom",
        "year": 2024, "estimate_status": "estimate", "projection_scenario": "historical", "sex": sex,
        "age_lower": age, "age_upper_exclusive": age+1 if age < 100 else None,
        "population": str(count if age < 100 else tail), "unit": "persons", "population_variant": "invented_residents"}
        for sex, count, tail in (("female", 1000, 10), ("male", 900, 8)) for age in range(101)]
    candidate = {"snapshot_id": "synthetic-demo-m2c", "rows": rows, "input_source_hash": "0"*64,
        "transformed_file_hash": "1"*64, "manifest": {"source_id": "synthetic", "source_category": "user",
        "source_vintage": "M2C", "review_status": "reviewed", "manifest_hash": "2"*64,
        "schema_version": "country_tb.demography_snapshot.v1", "transformation_version": "synthetic_fixture",
        "transformations": ["Invented separate sex-specific counts, ages 0-99 and an unestimated 100+ tail"],
        "source_url": "urn:synthetic:M2C", "import_date": "2026-10-28", "terms": {"status": "invented_software_fixture"},
        "citation": "Invented validation data; not a country population", "projection_boundary": None}}
    demo = demographic_profile(geo, [candidate])
    return shape, history, observations, policy, calibration, baseline_tbi_profile(demo, calibration)


def shape_files(shape):
    source = {r["year"]: r for r in shape["source_incidence"]["rows"]}
    diagnostics = {r["year"]: r for r in shape["diagnostics"]}
    rows = []
    for row in shape["rows"]:
        year = row["year"]
        rows.append({"year": year, "raw_incidence": source[year]["central"],
            "who_lower": source[year]["lower"], "who_upper": source[year]["upper"],
            "smoothed_incidence": row["fitted_incidence"], "pressure_shape_proxy": row["shape"],
            "log_residual": diagnostics[year]["log_residual"],
            "local_log_slope_per_year": diagnostics[year]["local_log_slope_per_year"],
            "covid_period": diagnostics[year]["covid_period"], "truncated_edge": diagnostics[year]["truncated_edge"],
            "normalization_year": shape["normalization"]["year"],
            "model_exclusions": ";".join(EXCLUSIONS), "interpretation": shape["interpretation"]})
    return {"shape.json": canonical(shape) + b"\n", "shape_diagnostics.csv": csv_bytes(rows, tuple(rows[0])),
        "README.md": ("# Pressure-shape diagnostic\n\n"
            "Relative infection-pressure shape proxy; no infection scale or country TBI estimate.\n\n"
            "Review status: " + shape["review"]["status"] + ". WHO bounds are preserved, not pressure intervals.\n\n"
            "Exclusions: " + "; ".join(EXCLUSIONS) + ".\n").encode()}


def country_diagnostics(source, demography, *, supplied_alpha='.005'):
    """Invented scale .005 is a software scenario, never a country estimate."""
    review = {'status': 'reviewed', 'source_id': 'diagnostic:M2C:user-authorized-scenario',
        'rationale': 'Computational diagnostic selection only; no empirical scientific approval of country pressure or scale.'}
    profiles = {}
    for bound in ('lower','central','upper'):
        q = shape_from_incidence(source, review=review, source_bound=bound)
        for backcast in ('constant_earliest','bounded_backward_trend'):
            parameters = None if backcast == 'constant_earliest' else {
                'slope_per_year': '-.01', 'minimum_shape': '0', 'maximum_shape': '100', 'review': review}
            h = apply_backcast(q, start_time='1900', method=backcast, parameters=parameters)
            policy = calibration_policy(alpha_upper='1', observation_weights={}, residual_tolerance='0', review=review)
            fit = calibrate_scale(h, policy, supplied_alpha=supplied_alpha,
                scale_provenance={'source_id': 'synthetic:diagnostic-supplied-scale',
                    'diagnostic_only': True, 'not_a_country_tbi_estimate': True})
            name = 'who_'+bound+'__'+backcast
            profiles[name] = baseline_tbi_profile(demography, fit, baseline_time='2024.5', scenario_id=name)
    return uncertainty_scenario_set(profiles, review=review)


def synthetic_assay_anchor_scenarios(profile, history):
    """Invented prevalence/uncertainty and assay parameters, with no pseudo-counts."""
    profiles = {}
    for name,se,sp,prevalence in [('central','.85','.96','.2'), ('assay_low','.75','.98','.2'),
            ('assay_high','.95','.92','.2'), ('anchor_lower','.85','.96','.15'), ('anchor_upper','.85','.96','.25')]:
        op = assay_operator(assay_class='IGRA', name='Invented IGRA diagnostic operator',
            positivity_threshold='synthetic_only', sensitivity=se, specificity=sp, review=REVIEW,
            provenance={'source_id':'synthetic:assay-parameters-not-evidence'})
        record = {'observation_id':'synthetic:overall', 'geography_id':history['geography_id'],
            'geography_level':'custom','observation_policy':'diagnostic_synthetic_anchor',
            'survey_year':2024,'survey_time':'2024','age_lower':'20','age_upper_exclusive':'21','sex':'both',
            'numerator':None,'denominator':None,'prevalence':prevalence,'review':REVIEW,
            'source':{'source_id':'synthetic:anchor-not-survey'},'sampling_frame':'invented exact age 20',
            'survey_design':{'kind':'synthetic_only'},'uncertainty':{'kind':'interval','lower':'.15','upper':'.25','method':'invented diagnostic range'},
            'assay_operator':op,'weight_source':{'source_id':'synthetic:equal-sex-weights'},
            'strata':[{'age_lower':'20','age_upper_exclusive':None,'age_policy':'exact_age','sex':s,'weight':'1'}
                for s in ('female','male')]}
        observations = observations_profile([record],review=REVIEW)
        policy = calibration_policy(alpha_upper='1',observation_weights={'synthetic:overall':'1'},
            residual_tolerance='.000001',review=REVIEW)
        fit = calibrate_scale(history,policy,observations=observations)
        profiles[name] = baseline_tbi_profile(profile['source_demography'],fit,scenario_id='synthetic_'+name)
    return uncertainty_scenario_set(profiles,review=REVIEW)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "abm/output/country_tb_m2c_preview")
    parser.add_argument("--shape-countries", nargs="*", default=["AUS", "VNM", "PRK"])
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if not output.is_relative_to((ROOT / "abm/output").resolve()):
        raise DataValidationError("Developer exports must stay under abm/output.")
    shape, history, observations, policy, calibration, profile = synthetic_example()
    files = preview_files(profile)
    files.update({"pressure_shape.json": canonical(shape)+b"\n", "pressure_history.json": canonical(history)+b"\n",
                  "observations.json": canonical(observations)+b"\n", "calibration_policy.json": canonical(policy)+b"\n",
                  "calibration.json": canonical(calibration)+b"\n"})
    receipt = artifact("country_tb.calibration_export_receipt.v1", profile_hash_reference=profile["profile_hash"],
                      files={name: digest(payload) for name, payload in sorted(files.items())})
    files["manifest.json"] = canonical(receipt)+b"\n"
    print(publish_immutable(output, "synthetic-calibration-" + profile["profile_hash"][:12], files))
    synthetic_scenarios = synthetic_assay_anchor_scenarios(profile, history)
    print(publish_immutable(output, 'synthetic-scenarios-'+synthetic_scenarios['profile_hash'][:12],
        {'scenarios.json':canonical(synthetic_scenarios)+b'\n'}))
    who = load_who_snapshot(WHO_2025) if args.shape_countries else None
    candidates = load_candidates() if args.shape_countries else None
    for country in args.shape_countries:
        source = national_profile(who, country)
        try:
            q = shape_from_incidence(source, review={"status": "not_reviewed", "source_id": WHO_2025,
                "rationale": "Developer preview only; no scientific shape approval assigned."})
            files = shape_files(q)
            identity = "who-" + country + "-shape-" + q["profile_hash"][:12]
        except DataValidationError as exc:
            unavailable = artifact("country_tb.pressure_shape_unavailable.v1", geography_id="who:"+country,
                source_incidence_hash=source["profile_hash"], reason=str(exc))
            files = {"unavailable.json": canonical(unavailable)+b"\n"}
            identity = "who-" + country + "-unavailable-" + unavailable["profile_hash"][:12]
        print(publish_immutable(output, identity, files))
        if country != 'PRK':
            scenarios = country_diagnostics(source, demographic_profile('who:'+country,candidates))
            central = next(r['profile'] for r in scenarios['rows'] if r['scenario_id']=='who_central__constant_earliest')
            files = preview_files(central)
            files['scenarios.json'] = canonical(scenarios)+b'\n'
            print(publish_immutable(output,'diagnostic-'+country+'-'+scenarios['profile_hash'][:12],files))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
