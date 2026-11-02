import pandas as pd
import streamlit as st
from engine.country_tb.review_parameters_3d import evidence_registry,progression_sets
from country_model_review_ui.forward import RR74_WARNING


def render(config):
    st.header('Evidence and unresolved decisions')
    st.warning('Milestone3D was internal scientific evidence synthesis and model checking,not external validation or scientific sign-off.')
    st.subheader('Held-out empirical positivity comparisons')
    st.dataframe(pd.DataFrame([
        {'Diagnostic dataset':'Viet Nam contacts2020','RMSE (percentage points)':8.25},
        {'Diagnostic dataset':'Viet Nam contacts2021','RMSE (percentage points)':2.71},
        {'Diagnostic dataset':'US-born NHANES','RMSE (percentage points)':1.18},
        {'Diagnostic dataset':'2020 contact fit transported to2021','RMSE (percentage points)':22.14}]),hide_index=True,use_container_width=True)
    st.write('The Viet Nam samples are selected household contacts in different campaign districts,not national probability samples or a longitudinal national cohort. Assay,exposure and selection differ. These comparisons do not validate national Viet Nam inference.')
    st.write('Australian premigration child applications are not a homogeneous Australian resident population. Mixed assays,pediatric uncertainty,selection and residence/exposure histories prevent a national latent-profile claim. US NHANES offers probability-sample age/nativity estimates,with specificity and open-age limitations.')
    st.subheader('Selectable progression assumptions')
    st.write('The adult four-point curve is sensitivity-only: focal modelled infection at age20,historical cohort inputs,unresolved net-hazard transport. The unverified2-year5% knot remains excluded. Pediatric curves remain external-required; adult hazards are not scaled into child curves.')
    st.write('No extrapolation means unresolved beyond the evidence window. Fixed and final-segment tails are explicit sensitivities,not approved central late-reactivation hazards. Reinfection resets the most-recent-event clock centrally; prior remote history is auditable without another person count.')
    rows=[]
    for s in progression_sets():
        rows.append({'Set':s['parameter_set_id'],'Cumulative risk knots':'; '.join(r['time_years']+'y: '+f"{float(r['cumulative_risk'])*100:.1f}%" for r in s['cumulative_risk_knots'] if r['time_years']!='0'),
            'Status':'Sensitivity only; no central approval','Population':s['applicable_population']})
    st.dataframe(pd.DataFrame(rows),hide_index=True,use_container_width=True)
    st.subheader('Risk evidence and definition compatibility');st.warning(RR74_WARNING)
    registry=evidence_registry()
    st.dataframe(pd.DataFrame([{'Factor':r['factor_id'],'Legacy definition':r['exact_legacy_definition'],'Review status':r['status'],
        'Effect measure':r['effect_measure'],'Evidence population / limitations':r['extraction_notes'],'Conditional progression multiplier approved':False,
        'Source':r['url']} for r in registry['factor_rows']]),hide_index=True,use_container_width=True)
    st.write('Disease associations can include acquisition,ascertainment and confounding; they are not automatically conditional-TBI hazard ratios. Fully calcified granulomas differ from fibrotic old-TB lesions. HIV treatment/CD4 and transplant organ/regimen matter. No universal cap is introduced.')
    st.subheader('Mortality,assay and observation decisions')
    st.write('ABS2022–2024 life-table values are a validation snapshot,not an adopted competing-mortality default. Country/sex/subgroup source licensing,all-cause versus TB-death scope and open ages need review. Australian mortality is not transferred to another country.')
    st.write('TST,QFT-GIT and QFT-Plus are distinct. Adult pooled test performance is not a latent-viability gold standard and does not automatically transport to children,BCG settings or immunosuppression. One observation cannot jointly identify assay performance and infection-pressure scale.')
    st.write('The operator predicts first episodes after baseline; WHO incidence includes new/recurrent episodes and may include exposure acquired elsewhere. Prior TB,imported disease,recurrence and historical depletion need compatible data. No fractions are invented or subtracted.')
    st.subheader('Scientific-review questions')
    st.markdown('- Which adult,child and late-tail estimands are transportable?\n- Which factor definitions and conditional effects are defensible?\n- Which mortality,assay and migration inputs match the population?\n- Which episode/ascertainment target is compatible?\n- Which unresolved-population and numerical thresholds should be approved?')
    st.caption('Static inversion is a diagnostic/initialization tool,not a substitute for dynamic transmission. No Summer/Starsim workstream is connected.')
    st.markdown('[3D reviewer document](https://github.com/emmamcbryde/tb-community-risk/blob/feature/country-subnational-tb-model-v1/docs/country_subnational_model/milestone3d_review.md) · [Evidence and sources](https://github.com/emmamcbryde/tb-community-risk/blob/feature/country-subnational-tb-model-v1/docs/country_subnational_model/milestone3d_evidence.md)')
