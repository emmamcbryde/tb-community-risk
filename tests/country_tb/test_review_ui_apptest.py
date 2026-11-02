"""Rendered Streamlit1.32+ interactions; separate from browser verification."""
from copy import deepcopy
from pathlib import Path
import ast
import json
import time
import pytest
from streamlit.testing.v1 import AppTest
from streamlit.source_util import get_pages

from country_model_review_ui.state import PAGES,SOURCE_MODES,working_example,default_config,TBI_COLUMNS,csv_text,material_hash
from country_model_review_ui.model import catalog
from engine.country_tb.data_common import REPO_ROOT
from engine.country_tb.risk_factor_catalog import catalogue_records

ENTRY=REPO_ROOT/'country_model_review_ui/app.py'


def app(config=None):
    a=AppTest.from_file(str(ENTRY),default_timeout=45)
    if config is not None:a.session_state['review_config']=deepcopy(config)
    a.run();assert not a.exception
    return a


def page(a,name):
    a.sidebar.radio[0].set_value(name).run();assert not a.exception;return a


def run(a):
    page(a,PAGES[4]);a.button(key='run_user').click().run();assert not a.exception and not a.error
    return a.session_state['review_result']


def all_text(elements):return '\n'.join(str(x.value) for x in elements)


def test_only_isolated_navigation_and_file_discovery():
    a=app();assert a.sidebar.radio[0].options==list(PAGES)
    discovered=get_pages(str(ENTRY));assert len(discovered)==1
    assert all(Path(v['script_path']).parent==ENTRY.parent for v in discovered.values())
    assert not (ENTRY.parent/'pages').exists()
    for p in ENTRY.parent.rglob('*.py'):
        tree=ast.parse(p.read_text(encoding='utf8'))
        for n in ast.walk(tree):
            if isinstance(n,ast.ImportFrom) and n.module:
                assert not n.module.startswith(('app.','apps.','engine.apy','explicit_tbi_ui','pages.','starsim','summer'))
    assert any(m.value=='27,204,809' for m in a.metric)
    assert any(m.value=='6,588' for m in a.metric)


def test_all_seven_pages_render_without_unexpected_errors():
    a=app(working_example())
    for name in PAGES:
        page(a,name)
        assert not a.error or all('incidence inversion is disabled' in x.value for x in a.error)


def test_restore_working_example_and_full_run():
    a=app();a.button(key='restore_working_example').click().run()
    assert a.session_state['review_config']['example'] and not a.exception
    r=run(a);assert float(r['baseline']['summary']['components']['any_tbi']['prevalence_full'])==pytest.approx(.2)
    page(a,PAGES[5]);assert not a.exception
    assert any(m.label=='Population retained' and m.value=='1,000' for m in a.metric)
    assert 'not externally validated' in all_text(a.warning)
    assert any(m.value.startswith('48.18') for m in a.metric)


def test_overall_user_entry_total_only_rendered():
    c=working_example();c['tbi'].update(derive=False,overall_percent='')
    a=app(c);page(a,PAGES[1]);a.text_input(key='edit_Overall prevalence (%)').input('12.3456').run()
    a.button(key='build_TBI').click().run();assert not a.exception and not a.error
    preview=a.session_state['review_TBI_preview']['baseline']
    assert all(r['any_tbi']=='0.123456' and r['recent'] is None for r in preview['rows'])
    run(a);page(a,PAGES[5]);assert 'recent/remote composition' in all_text(a.warning)


@pytest.mark.parametrize('bysex',[False,True])
def test_age_and_age_by_sex_csv_through_rendered_input(bysex):
    c=working_example();c['tbi'].update(derive=False,input_kind='Age-by-sex total profile' if bysex else 'Age-specific total profile')
    base=dict(subgroup_id='main',source='Invented render fixture',review_status='diagnostic',geography_id='custom:synthetic_review',baseline_year=2024,
        allocation_policy='piecewise_constant_within_band_v1',recent='',remote_only='',never_infected='',prior_remote_plus_recent='')
    rows=[dict(base,record_id=s,age_lower=0,age_upper_exclusive='',sex=s,any_tbi=p) for s,p in [('female','.1'),('male','.3')]] if bysex else [
        dict(base,record_id='young',age_lower=0,age_upper_exclusive=45,sex='both',any_tbi='.1'),
        dict(base,record_id='older',age_lower=45,age_upper_exclusive='',sex='both',any_tbi='.3')]
    a=app(c);page(a,PAGES[1]);a.text_area(key='edit_TBI profile CSV').input(csv_text(rows,TBI_COLUMNS)).run()
    a.button(key='build_TBI').click().run();assert not a.error and not a.exception
    r=run(a);assert float(r['baseline']['summary']['components']['any_tbi']['prevalence_full'])==pytest.approx(.2)


def test_direct_recent_remote_rendered():
    c=working_example();c['tbi'].update(derive=False,input_kind='Direct recent and remote')
    a=app(c);page(a,PAGES[1]);a.text_input(key='edit_Recent (%)').input('5').run()
    a.text_input(key='edit_Remote-only (%)').input('15').run();a.text_input(key='edit_Prior remote + recent subset (%) — optional').input('2').run()
    r=run(a);assert r['baseline']['summary']['components']['recent']['complete_count']=='50'
    assert r['baseline']['summary']['components']['prior_remote_plus_recent']['complete_count']=='20'


def test_synthetic_subgroup_example_rendered():
    a=app();a.button(key='restore_subgroup_example').click().run();r=run(a)
    assert len(r['baseline']['subgroups'])==3 and r['baseline']['summary']['population']=='1000'
    assert float(r['baseline']['summary']['components']['any_tbi']['prevalence_full'])==pytest.approx(.1)
    page(a,PAGES[5]);assert not a.exception
    assert any('invented' in str(df.value).lower() for df in a.dataframe)


def test_incomplete_subgroups_visible_in_results():
    c=working_example(mixture=True);c['population']['subgroups'][1]['source_mode']='Unresolved'
    a=app(c);r=run(a);page(a,PAGES[5])
    assert r['baseline']['summary']['unresolved_population']=='200'
    assert 'Incomplete mixture' in all_text(a.error)
    assert any(m.value=='Unavailable' for m in a.metric)


def test_incidence_fork_and_direct_navigation_blocked():
    a=app(working_example());page(a,PAGES[1]);a.radio(key='edit_TBI_mode').set_value(SOURCE_MODES[1]).run()
    assert 'Experimental diagnostic' in all_text(a.warning) and 'blocked' in all_text(a.error).lower()
    page(a,PAGES[4]);assert a.button(key='run_incidence').disabled
    assert 'review_result' not in a.session_state
    page(a,PAGES[3]);a.checkbox(key='edit_Explicitly enable scientific incidence-to-TBI diagnostic').check().run()
    page(a,PAGES[4]);assert a.button(key='run_incidence').disabled
    assert not a.exception


def test_child_no_extrapolation_and_rr74_labels_rendered():
    c=working_example();a=app(c);page(a,PAGES[3])
    a.selectbox(key='edit_Remote-infection tail').set_value('No extrapolation').run()
    assert 'never zero hazard' in all_text(a.caption)
    assert any('Pediatric' in str(x.value) for x in a.dataframe)
    page(a,PAGES[2]);assert '94.31%' in all_text(a.warning)
    assert 'not automatically applied below15' in all_text(a.info)


def test_requested_association_and_achieved_phi_displayed():
    c=working_example();definition=catalogue_records(catalog())['diabetes']['exposure_definition']
    c['risk'].update(no_factors=False,legacy_effects=True,records=[dict(record_id='diabetes',subgroup_id='all',factor_id='diabetes',age_lower=15,
        age_upper_exclusive=None,sex='both',state='value',prevalence='.2',allocation_policy='piecewise_constant_within_supplied_band_v1',
        supplied_exposure_definition=definition,source='Synthetic fixture',review_status='diagnostic',reason='Synthetic')])
    a=app(c);r=run(a);page(a,PAGES[5])
    assert any({'Requested_c','Achieved_phi','Feasible_lower','Feasible_upper'}<=set(d.value.columns) for d in a.dataframe)
    page(a,PAGES[2]);a.text_input(key='edit_Normalized association c (−1 to +1)').input('.5').run()
    assert any('Preview' in x.value or 'preview' in x.value for x in a.caption)
    assert any('Results stale' in x.value for x in a.warning)


def test_material_edits_stale_and_hide_old_result():
    a=app(working_example());r=run(a);old=r['profile_hash'];page(a,PAGES[1])
    a.text_input(key='edit_Overall prevalence (%)').input('30').run();assert 'Results stale' in all_text(a.warning)
    page(a,PAGES[5]);assert any('Results stale' in x.value for x in a.warning)
    assert not a.metric and a.session_state['review_result']['profile_hash']==old


def test_no_eager_exports_and_precision_download_preparation(monkeypatch):
    a=app(working_example());run(a)
    import country_model_review_ui.views.results as view
    original=view.prepare_exports
    monkeypatch.setattr(view,'prepare_exports',lambda *args:pytest.fail('Export eagerly generated'))
    page(a,PAGES[5]);assert 'review_exports' not in a.session_state
    monkeypatch.setattr(view,'prepare_exports',original);a.button(key='prepare_review_exports').click().run()
    assert not a.exception
    data=a.session_state['review_exports'];value=json.loads(data['json'])
    assert value['run']['configuration_hash']==material_hash(a.session_state['review_config']) and not value['production_available']
    assert 'exact' in data['json'].decode('utf8')


def test_evidence_comparisons_and_validation_only_labels():
    a=app(working_example());page(a,PAGES[6]);df=a.dataframe[0].value
    assert list(df['RMSE (percentage points)'])==[8.25,2.71,1.18,22.14]
    assert 'not external validation' in all_text(a.warning)
    run(a);page(a,PAGES[5]);assert 'Validation only,not calibration' in all_text(a.warning)
    assert any('Observed WHO' in str(c.proto) for c in a.get('arrow_vega_lite_chart'))


def test_duplicate_and_remove_subgroup_callbacks():
    a=app(working_example(mixture=True));a.button(key='duplicate_subgroup').click().run();assert not a.exception
    groups=a.session_state['review_config']['population']['subgroups'];assert len(groups)==4 and groups[-1]['population']==''
    duplicate=groups[-1]['subgroup_id'];a.selectbox(key='edit_selected_group').set_value(duplicate).run()
    a.button(key='remove_subgroup').click().run();assert not a.exception
    assert len(a.session_state['review_config']['population']['subgroups'])==3
