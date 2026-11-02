"""Launch: python -m streamlit run country_model_review_ui/app.py

This subdirectory has no pages/ folder. Navigation is an explicit isolated
allow-list; root legacy pages and deployed entrypoints are never imported.
"""
from pathlib import Path
from time import perf_counter
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import streamlit as st
from engine.country_tb.data_common import DataValidationError
from country_model_review_ui import APP_VERSION
from country_model_review_ui.state import PAGES,initialize,invalidate_exports,stale
from country_model_review_ui.views import population,baseline,risk,natural,run,results,evidence

VIEWS=dict(zip(PAGES,(population.render,baseline.render,risk.render,natural.render,run.render,results.render,evidence.render)))


def main():
    start=perf_counter();st.set_page_config(page_title='Country TB · Scientific review',page_icon='🔬',layout='wide')
    initialize(st.session_state)
    st.markdown('''<style>
    .block-container{padding-top:2rem;max-width:1400px} h1{letter-spacing:-.03em}
    [data-testid="stSidebar"]{background:#eff3f6} .stButton>button[kind="primary"]{background:#6247aa;border-color:#6247aa}
    .review-kicker{color:#177e89;font-size:.8rem;font-weight:700;letter-spacing:.12em;margin-bottom:.4rem}
    </style>''',unsafe_allow_html=True)
    st.markdown('<div class="review-kicker">COUNTRY / SUBNATIONAL TB · INTERNAL RESEARCH WORKSPACE</div>',unsafe_allow_html=True)
    st.title('Scientific review workspace')
    st.caption('Inspect data,assumptions and unresolved calculations. Not a production application or an externally validated model. For planning and sequencing,not denying care.')
    with st.sidebar:
        st.subheader('Review workflow');page=st.radio('Review pages',PAGES,key='review_page',label_visibility='collapsed')
        st.caption('Seven isolated review pages. No legacy workflows or economics pages.')
        st.markdown('---');st.caption('Observed · User supplied · Model derived · Diagnostic sensitivity · Unresolved · Calculation blocked')
        st.caption(APP_VERSION)
    status=st.empty();config=st.session_state['review_config']
    try:VIEWS[page](config)
    except (DataValidationError,ValueError,KeyError,TypeError) as exc:
        st.error('Calculation blocked — '+str(exc))
        st.caption('Inputs are preserved. Correct the source/schema/assumption issue and retry; no fallback or clipping was applied.')
    invalidate_exports(st.session_state)
    if stale(st.session_state):status.warning('Results stale — material inputs changed. Rerun before interpreting prior results.')
    elif st.session_state.get('review_result'):status.info('Current internal review run · scientific approval pending')
    else:status.caption('No accepted diagnostic run for these inputs.')
    if 'first_app_render_seconds' not in st.session_state['review_timings']:
        st.session_state['review_timings']['first_app_render_seconds']=perf_counter()-start


if __name__=='__main__':main()
