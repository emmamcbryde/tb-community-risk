"""Bounded offline AppTest evidence and3E performance receipt.

This records rendered interaction,not a browser screenshot or external validation.
Generated files are confined to abm/output. No model defaults are approved.
"""
from pathlib import Path
import argparse
from datetime import datetime,timezone
import json
import sys
from time import perf_counter

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'abm/output/country_tb_m3e_rendered');p.add_argument('--prospective',action='store_true');args=p.parse_args()
    output=args.output.resolve()
    if not output.is_relative_to((ROOT/'abm/output').resolve()):raise ValueError('Outputs must remain under abm/output.')
    start=perf_counter()
    import streamlit
    from streamlit.testing.v1 import AppTest
    from country_model_review_ui.state import PAGES,working_example,SOURCE_MODES,material_hash
    from country_model_review_ui.model import execute_user_run,select_population
    from country_model_review_ui.gates import prepare_diagnostic
    from engine.country_tb.data_common import canonical,digest,publish_immutable,DataValidationError
    app=AppTest.from_file(str(ROOT/'country_model_review_ui/app.py'),default_timeout=60).run()
    assert not app.exception and not app.error
    timings={'first_app_load_seconds':perf_counter()-start}
    country=app.selectbox(key='edit_country');vietnam=next(x for x in country.options if x.endswith('(VNM)'))
    start=perf_counter();country.set_value(vietnam).run();app.button(key='load_population').click().run();assert not app.exception and not app.error
    timings['first_Viet_Nam_frozen_demography_selection_seconds']=perf_counter()-start
    app.button(key='restore_working_example').click().run();assert not app.exception
    app.sidebar.radio[0].set_value(PAGES[4]).run();app.button(key='run_user').click().run()
    assert not app.exception and not app.error
    result=app.session_state['review_result'];timings.update(app.session_state['review_timings'])
    evidence=[]
    for page in PAGES:
        start=perf_counter();app.sidebar.radio[0].set_value(page).run();elapsed=perf_counter()-start
        assert not app.exception
        errors=[x.value for x in app.error]
        assert not errors or all('incidence inversion is disabled' in e for e in errors),errors
        evidence.append(dict(page=page,elapsed_seconds=elapsed,headers=[x.value for x in app.header],
            captions=[x.value for x in app.caption],warnings=[x.value for x in app.warning],expected_blocked_messages=errors,
            metrics=[dict(label=x.label,value=x.value) for x in app.metric],
            tables=[dict(columns=list(x.value.columns),rows=len(x.value)) for x in app.dataframe]))
        if page==PAGES[5]:timings['results_page_wall_seconds']=elapsed
    config=working_example();blocked_config=working_example();blocked_config['tbi']['mode']=SOURCE_MODES[1]
    context=select_population(blocked_config);start=perf_counter()
    try:prepare_diagnostic(blocked_config,context)
    except DataValidationError as exc:blocked=str(exc)
    else:raise AssertionError('Default incidence diagnostic unexpectedly enabled')
    timings['default_incidence_gate_refusal_seconds']=perf_counter()-start
    prospective=None
    if args.prospective:
        config['natural']['horizon']='One-year prospective (bounded)';start=perf_counter();prospective,t=execute_user_run(config)
        timings['optional_one_year_run_wall_seconds']=perf_counter()-start;timings['optional_one_year_forward_seconds']=t['deterministic_forward_seconds']
        assert prospective['projection']['prospective'] is not None,prospective['projection'].get('prospective_unavailable_reason')
    receipt=dict(schema_id='country_tb.review_ui_rendered_receipt.v1',streamlit_version=streamlit.__version__,
        runtime_basis='wall_clock; cold frozen demography vs warmed working example explicitly separated',
        production_available=False,externally_validated=False,browser_verified=False,
        verification_type='Streamlit_AppTest_not_live_browser',timings=timings,rendered_pages=evidence,
        default_incidence_gate=blocked,working_example='Synthetic adult cohort,not Australian estimates',
        example_configuration_hash=material_hash(result['configuration']),example_result_hash=result['profile_hash'],
        source_integrity='immutable source loaders; no network requests',
        files={str(p.relative_to(ROOT)).replace('\\','/'):digest(p.read_bytes()) for base in ['country_model_review_ui','tests/country_tb']
            for p in (ROOT/base).rglob('*.py') if base=='country_model_review_ui' or 'review_ui' in p.name})
    files={'rendered_receipt.json':canonical(receipt)+b'\n','working_example.json':canonical(result)+b'\n'}
    if prospective is not None:files['one_year_example.json']=canonical(prospective)+b'\n'
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    folder=publish_immutable(output,'review-ui-'+stamp,files)
    print(json.dumps(dict(directory=str(folder),timings=timings,rendered_pages=len(evidence),browser_verified=False,
        synthetic_instantaneous_rate=result['projection']['instantaneous']['full_population_rate_per100000PY'],
        one_year_status=None if prospective is None else prospective['projection']['status'])))


if __name__=='__main__':main()
