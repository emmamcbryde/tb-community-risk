import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
from scipy.optimize import minimize_scalar
import scipy
from engine.dynamic.exec_dynamic import run_dynamic_model
from engine.infection_backcast import (
    calc_ari_from_incidence,
    infection_prob_by_age_split,
)

# =====================================================
# Calibration-controlled constants (NOT user-facing)
# =====================================================
BASELINE_DIAG_MONTHS = 12.0  # mean time to treatment pre-intervention
INCIDENCE_FLOOR = 0.1  # per 100k, avoids zero pathologies
ARI_FLOOR = 1e-6  # avoids ARI numerical issues


# =====================================================
# Default fallback age distribution
# =====================================================
def default_age_distribution():
    ages = list(range(0, 101))
    prop = np.array([1 / 101] * 101)
    return pd.DataFrame({"AgeGroup": ages, "Proportion": prop})


# =====================================================
# Load country-specific OWID population structure
# =====================================================
def load_population_data(
    country_code="AUS", file_path="data/population_age_latest.csv"
):

    df = pd.read_csv(file_path)
    df_country = df[df["iso_code"].str.upper() == country_code.upper()]

    if df_country.empty:
        st.warning(f"No population structure for {country_code}. Using default global.")
        return default_age_distribution(), pd.DataFrame(
            {"age": range(0, 101), "population": [1] * 101}
        )

    total_pop_country = df_country["population"].sum()
    df_country = df_country.copy()
    df_country["Proportion"] = df_country["population"] / total_pop_country

    bin_edges = list(range(0, 105, 5))
    bin_labels = [
        f"{bin_edges[i]}–{bin_edges[i+1]-1}" for i in range(len(bin_edges) - 1)
    ]
    bin_labels.append("100+")

    df_country["AgeBin"] = pd.cut(
        df_country["age"], bins=bin_edges + [200], labels=bin_labels, right=False
    )

    df_age_groups = df_country.groupby("AgeBin", as_index=False)[
        ["population", "Proportion"]
    ].sum()
    return df_age_groups, df_country


# =====================================================
# Incidence history builder
# =====================================================
def build_incidence_history(
    hist_pattern, user_incidence, years_back_needed, uploaded_inc_df=None
):
    """
    Returns dict: year offset (0=now, -1 last year, ...) -> incidence per 100k
    Includes flooring + smoothing for stability.
    """
    years_back = list(range(0, years_back_needed + 1))

    if hist_pattern == "Constant":
        inc_hist = {-k: float(user_incidence) for k in years_back}

    elif hist_pattern == "Falling 3%/year":
        inc_hist = {-k: float(user_incidence) * (1.03**k) for k in years_back}

    elif hist_pattern == "Rising 3%/year":
        inc_hist = {-k: float(user_incidence) * (0.97**k) for k in years_back}

    elif hist_pattern == "Upload CSV (year, incidence)":
        if uploaded_inc_df is None:
            st.warning("Incidence CSV not uploaded yet – using constant incidence.")
            inc_hist = {-k: float(user_incidence) for k in years_back}
        else:
            df = uploaded_inc_df.copy().sort_values("year")
            if not {"year", "incidence"}.issubset(df.columns):
                st.warning(
                    "Incidence CSV missing required columns – using constant incidence."
                )
                inc_hist = {-k: float(user_incidence) for k in years_back}
            else:
                years = df["year"].values.astype(int)
                incs = df["incidence"].values.astype(float)

                # floor to avoid zeros
                incs = np.maximum(incs, INCIDENCE_FLOOR)

                year_min = int(years[0])
                year_max = int(years[-1])

                inc_min = float(np.min(incs))
                inc_max = float(np.max(incs))

                # estimate geometric trend ignoring tiny values
                ratios = []
                for i in range(1, len(incs)):
                    if incs[i - 1] > INCIDENCE_FLOOR and incs[i] > INCIDENCE_FLOOR:
                        ratios.append(incs[i] / incs[i - 1])
                trend = float(np.exp(np.mean(np.log(ratios)))) if ratios else 1.0

                inc_map = dict(zip(years, incs))

                # Treat "now" as the last year in the CSV
                ref_year = year_max

                inc_hist = {}
                for k in years_back:
                    target_year = ref_year - k

                    if year_min <= target_year <= year_max:
                        nearest = min(
                            inc_map.keys(), key=lambda y: abs(y - target_year)
                        )
                        inc_hist[-k] = float(inc_map[nearest])

                    elif target_year > year_max:
                        # forward extrapolation (if needed)
                        j = target_year - year_max
                        extrap = incs[-1] * (trend**j)
                        inc_hist[-k] = float(min(extrap, inc_max))

                    else:
                        # backward extrapolation
                        j = year_min - target_year
                        extrap = incs[0] * (trend ** (-j))
                        inc_hist[-k] = float(max(extrap, inc_min))

    else:
        inc_hist = {-k: float(user_incidence) for k in years_back}

    # smooth with 3-year moving average
    inc_series = pd.Series(inc_hist).sort_index()
    inc_series = inc_series.rolling(window=3, center=True, min_periods=1).mean()

    # final floor
    inc_hist = {
        int(k): float(max(v, INCIDENCE_FLOOR)) for k, v in inc_series.to_dict().items()
    }
    return inc_hist


def compute_ltbi_from_inc_hist(ages, inc_hist, shift_years=0):
    """
    Compute LTBI probabilities at a reference time shift_years in the past.
    shift_years=0 => now
    shift_years=20 => treat incidence at -20 as 'now' for LTBI calculation
    """
    max_age = int(max(ages))
    min_key = int(min(inc_hist.keys()))

    inc_ref = {}
    for k in range(0, max_age + 1):
        src_key = -(k + shift_years)
        inc_ref[-k] = float(
            inc_hist.get(src_key, inc_hist.get(min_key, INCIDENCE_FLOOR))
        )

    ari_hist = calc_ari_from_incidence(inc_ref)
    ari_hist = {t: max(float(a), ARI_FLOOR) for t, a in ari_hist.items()}

    ltbi_ever, ltbi_recent, ltbi_remote = infection_prob_by_age_split(ages, ari_hist)
    return ltbi_ever, ltbi_recent, ltbi_remote


def calibrate_beta_to_history(
    age_counts,
    ages,
    inc_hist,
    calib_years,
    risk_inputs,
    pre_det_months,
    delta_pre,
    beta_bounds=(0.0, 50.0),
):
    """
    Fit beta to match observed incidence history over the last calib_years (excluding year 0).
    Observed series used: offsets -calib_years .. -1
    """

    total_pop = float(sum(age_counts.values()))

    # observed incidence per 100k (oldest -> newest)
    obs = np.array([inc_hist[-k] for k in range(calib_years, 0, -1)], dtype=float)

    # LTBI state at start of calibration window
    ltbi_ever0, ltbi_recent0, _ = compute_ltbi_from_inc_hist(
        ages, inc_hist, shift_years=calib_years
    )

    # incidence at the start year (used to seed I0)
    inc0 = float(inc_hist.get(-calib_years, obs[0]))

    base_params = {
        "age_counts": age_counts,
        "ltbi_ever": ltbi_ever0,
        "ltbi_recent": ltbi_recent0,
        "initial_incidence_per_100k": inc0,
        "pre_det_months": pre_det_months,
        "delta_pre": delta_pre,
        "delta_post": delta_pre,
        "ltbi_coverage": 0.0,
        "rollout_years": 0,
        "treatment_method": "None",
        "testing_method": "None",
    }
    base_params.update(risk_inputs)

    def objective(beta):
        p = dict(base_params)
        p["beta"] = float(beta)

        sim = run_dynamic_model(p, years=calib_years, intervention=False)
        model_inc_per100k = sim["annual_incidence"] * 100000.0 / total_pop

        # annual_incidence length = calib_years; align to obs length = calib_years
        err = model_inc_per100k - obs
        rmse = float(np.sqrt(np.mean(err**2)))
        return rmse

    res = minimize_scalar(objective, bounds=beta_bounds, method="bounded")
    beta_hat = float(res.x)
    rmse = float(res.fun)

    # re-run to return fitted series for plotting
    p_final = dict(base_params)
    p_final["beta"] = beta_hat
    sim_fit = run_dynamic_model(p_final, years=calib_years, intervention=False)
    fit_inc_per100k = sim_fit["annual_incidence"] * 100000.0 / total_pop

    return beta_hat, rmse, obs, fit_inc_per100k


# =====================================================
# Main Dynamic Model UI
# =====================================================
def render_dynamic_ui():

    st.header("📈 Dynamic LTBI → TB Model (with calibration)")
    st.info(
        "This model estimates TB incidence based on LTBI, risk factors, and interventions.\n\n"
        "• Baseline = no new LTBI test-and-treat\n"
        "• Intervention = selected LTBI testing and treatment\n"
        "• Outputs are annualised and shown as rates (per 100,000) and counts\n"
        "• Transmission (β) is calibrated to historical incidence"
    )

    # --------------------------------------------------
    # Core inputs
    # --------------------------------------------------
    population = st.sidebar.number_input("Population size", min_value=50, value=10000)
    user_incidence = st.sidebar.number_input(
        "Baseline annual incidence (per 100k)", 0, 500, 30
    )
    time_horizon = st.sidebar.slider("Time horizon (years)", 1, 30, 20)

    calib_years = st.sidebar.slider(
        "Calibration window (years)",
        5,
        30,
        15,
        help="β is fitted so that modelled incidence matches the historical incidence pattern over this many past years.",
    )

    # --------------------------------------------------
    # Risk factors
    # --------------------------------------------------
    smoker_pct = st.sidebar.slider(
        "Smoker (%)",
        0,
        100,
        30,
        help="Proportion of the population with current tobacco smoking.",
    )
    alcohol_pct = st.sidebar.slider(
        "Excess alcohol use (%)",
        0,
        100,
        15,
        help="Proportion of the population with excess alcohol consumption.",
    )
    diabetes_pct = st.sidebar.slider(
        "Diabetes (%)",
        0,
        100,
        10,
        help="Proportion of the population with diagnosed diabetes mellitus.",
    )
    renal_pct = st.sidebar.slider(
        "Renal impairment (%)",
        0,
        100,
        5,
        help="Proportion of the population with moderate–severe chronic kidney disease.",
    )
    HIV_treated_pct = st.sidebar.slider(
        "HIV Treated with antiretrovirals (%)",
        0,
        100,
        3,
        help="Proportion of the population with HIV under treatment.",
    )
    HIV_untreated_pct = st.sidebar.slider(
        "HIV Untreated (%)",
        0,
        100,
        3,
        help="Proportion of the population with HIV not under treatment.",
    )

    risk_inputs = {
        "smoker_pct": smoker_pct,
        "alcohol_pct": alcohol_pct,
        "diabetes_pct": diabetes_pct,
        "renal_pct": renal_pct,
        "HIV_treated_pct": HIV_treated_pct,
        "HIV_untreated_pct": HIV_untreated_pct,
    }

    # --------------------------------------------------
    # LTBI Test & Treat (pulse)
    # --------------------------------------------------
    testing_method = st.sidebar.selectbox("Testing method", ["TST", "IGRA", "None"])
    treatment_method = st.sidebar.selectbox(
        "Treatment regimen", ["1HP", "3HP", "4R", "6H", "9H", "None"]
    )

    ltbi_coverage = (
        st.sidebar.slider(
            "LTBI Test & Treat total coverage (%)",
            0,
            100,
            50,
            help="Effective fraction of the population that will successfully complete LTBI treatment over the rollout period.",
        )
        / 100.0
    )

    rollout_years = st.sidebar.slider("Rollout duration (years)", 1, 10, 5)

    # --------------------------------------------------
    # Diagnosis improvement: percent reduction slider (intervention only)
    # --------------------------------------------------
    st.sidebar.subheader("Diagnosis improvement (intervention)")
    diag_reduction_pct = st.sidebar.slider(
        "Percent reduction in time before treatment (%)",
        0,
        100,
        50,
        help="Reduces the mean time to diagnosis/treatment relative to baseline.",
    )

    pre_det_months = BASELINE_DIAG_MONTHS
    post_det_months = max(
        BASELINE_DIAG_MONTHS * (1.0 - diag_reduction_pct / 100.0), 0.1
    )

    delta_pre = 12.0 / pre_det_months
    delta_post = 12.0 / post_det_months

    # --------------------------------------------------
    # Historical incidence pattern
    # --------------------------------------------------
    st.sidebar.subheader("Historical Incidence Pattern")
    st.sidebar.markdown(
        "**Historical TB incidence (for LTBI back-calculation & calibration)**  \n"
        "Historical incidence determines LTBI by age, and is also used to calibrate β."
    )

    hist_pattern = st.sidebar.selectbox(
        "Choose pattern:",
        [
            "Constant",
            "Falling 3%/year",
            "Rising 3%/year",
            "Upload CSV (year, incidence)",
        ],
    )

    uploaded_inc_df = None
    if hist_pattern == "Upload CSV (year, incidence)":
        inc_file = st.sidebar.file_uploader("Upload incidence CSV", type="csv")
        if inc_file:
            try:
                tmp = pd.read_csv(inc_file)
                if {"year", "incidence"}.issubset(tmp.columns):
                    uploaded_inc_df = tmp.sort_values("year")
                    st.sidebar.success("Incidence history loaded.")
                else:
                    st.sidebar.error("CSV must contain columns: year, incidence")
            except Exception as e:
                st.sidebar.error(f"Could not read file: {e}")

    # --------------------------------------------------
    # Age distribution
    # --------------------------------------------------
    st.sidebar.subheader("Age Distribution")
    age_method = st.sidebar.radio(
        "Choose method:",
        ["Country ISO code (recommended)", "Upload custom CSV", "Default global"],
    )

    if age_method == "Country ISO code (recommended)":
        country = st.sidebar.text_input("ISO3 code", "AUS")
        age_df_display, df_country = load_population_data(country)

    elif age_method == "Upload custom CSV":
        file = st.sidebar.file_uploader(
            "Upload CSV with AgeGroup,Proportion", type="csv"
        )
        if file:
            df = pd.read_csv(file)
            if {"AgeGroup", "Proportion"}.issubset(df.columns):
                age_df_display = df
                df_country = pd.DataFrame(
                    {
                        "age": df["AgeGroup"].astype(int),
                        "population": df["Proportion"] * population,
                    }
                )
            else:
                st.error("CSV must include AgeGroup and Proportion. Using default.")
                age_df_display = default_age_distribution()
                df_country = pd.DataFrame(
                    {"age": range(0, 101), "population": [population / 101] * 101}
                )
        else:
            age_df_display = default_age_distribution()
            df_country = pd.DataFrame(
                {"age": range(0, 101), "population": [population / 101] * 101}
            )

    else:
        age_df_display = default_age_distribution()
        df_country = pd.DataFrame(
            {"age": range(0, 101), "population": [population / 101] * 101}
        )

    st.subheader("📊 Age Distribution (5-year bins)")
    st.dataframe(age_df_display)

    total_pop_country = df_country["population"].sum()
    age_counts = {
        int(row["age"]): population * (row["population"] / total_pop_country)
        for _, row in df_country.iterrows()
    }
    ages = sorted(age_counts.keys())
    max_age = int(max(ages))

    # Build incidence history deep enough for calibration + LTBI
    years_back_needed = max_age + calib_years + 5
    inc_hist = build_incidence_history(
        hist_pattern, user_incidence, years_back_needed, uploaded_inc_df=uploaded_inc_df
    )

    # LTBI now (for forward simulation)
    ltbi_ever, ltbi_recent, ltbi_remote = compute_ltbi_from_inc_hist(
        ages, inc_hist, shift_years=0
    )

    # --------------------------------------------------
    # LTBI by age stacked chart (≤ 60 years)
    # --------------------------------------------------
    ltbi_age_df = pd.DataFrame(
        {
            "Age": ages,
            "LTBI_recent": 100 * pd.Series(ltbi_recent),
            "LTBI_remote": 100 * (pd.Series(ltbi_ever) - pd.Series(ltbi_recent)),
        }
    )
    ltbi_age_df = ltbi_age_df[ltbi_age_df["Age"] <= 60]
    ltbi_age_df = ltbi_age_df.melt(id_vars="Age", var_name="Type", value_name="Percent")

    st.subheader("📉 LTBI Prevalence by Age (stacked %, ages 0–60)")
    chart = (
        alt.Chart(ltbi_age_df)
        .mark_area()
        .encode(
            x="Age:Q",
            y="Percent:Q",
            color="Type:N",
            tooltip=["Age", "Type", alt.Tooltip("Percent:Q", format=".1f")],
        )
    )
    st.altair_chart(chart, width="stretch")
    st.caption(
        "This chart shows estimated LTBI prevalence by age (≤60 years).\n\n"
        "• Recent LTBI = infection acquired in the last 5 years\n"
        "• Remote LTBI = infection acquired more than 5 years ago\n"
        "• The stacked height equals total LTBI prevalence at each age\n\n"
        "The shape reflects the selected historical TB incidence trajectory."
    )

    # --------------------------------------------------
    # Run simulations
    # --------------------------------------------------
    if st.sidebar.button("Run Dynamic Simulation"):

        st.info("Calibrating β to historical incidence…")

        beta_hat, rmse, obs_inc, fit_inc = calibrate_beta_to_history(
            age_counts=age_counts,
            ages=ages,
            inc_hist=inc_hist,
            calib_years=calib_years,
            risk_inputs=risk_inputs,
            pre_det_months=pre_det_months,
            delta_pre=delta_pre,
            beta_bounds=(0.0, 50.0),
        )

        st.success(
            f"Calibration complete: β = {beta_hat:.2f} (RMSE={rmse:.2f} per 100k)"
        )

        # Plot fit "back in time"
        if (
            hist_pattern == "Upload CSV (year, incidence)"
            and uploaded_inc_df is not None
        ):
            ref_year = int(uploaded_inc_df["year"].max())
            years_axis = [ref_year - k for k in range(calib_years, 0, -1)]
            x_title = "Calendar year"
        else:
            years_axis = list(range(-calib_years, 0))  # -calib_years .. -1
            x_title = "Years before present (relative)"

        df_cal = pd.DataFrame(
            {
                "Year": years_axis,
                "Observed": obs_inc,
                "Model (calibrated)": fit_inc,
            }
        ).melt(id_vars="Year", var_name="Series", value_name="Incidence_per100k")

        st.subheader("🧪 Calibration fit to historical incidence")
        st.caption(
            "Observed incidence vs modelled incidence over the calibration window."
        )
        cal_chart = (
            alt.Chart(df_cal)
            .mark_line()
            .encode(
                x=alt.X("Year:Q", title=x_title),
                y=alt.Y("Incidence_per100k:Q", title="Incidence per 100,000 per year"),
                color="Series:N",
                tooltip=[
                    "Year",
                    "Series",
                    alt.Tooltip("Incidence_per100k:Q", format=".1f"),
                ],
            )
        )
        st.altair_chart(cal_chart, use_container_width=True)

        # -------------------------
        # Build baseline + intervention params at "now"
        # -------------------------
        params_base = {}
        params_int = {}

        for p in (params_base, params_int):
            p["beta"] = beta_hat
            p.update(risk_inputs)

            p["ltbi_ever"] = ltbi_ever
            p["ltbi_recent"] = ltbi_recent
            p["age_counts"] = age_counts

            p["pre_det_months"] = pre_det_months
            p["delta_pre"] = delta_pre
            p["delta_post"] = delta_post

            p["initial_incidence_per_100k"] = float(user_incidence)

        # baseline = no intervention
        params_base["treatment_method"] = "None"
        params_base["testing_method"] = "None"
        params_base["ltbi_coverage"] = 0.0
        params_base["rollout_years"] = 0

        # intervention
        params_int["treatment_method"] = treatment_method
        params_int["testing_method"] = testing_method
        params_int["ltbi_coverage"] = ltbi_coverage
        params_int["rollout_years"] = rollout_years

        try:
            baseline = run_dynamic_model(
                params_base, years=time_horizon, intervention=False
            )
            intervention = run_dynamic_model(
                params_int, years=time_horizon, intervention=True
            )

            total_pop = float(sum(age_counts.values()))
            years_out = baseline["annual_incidence_time"]  # 0..time_horizon-1

            base_cases = baseline["annual_incidence"]
            int_cases = intervention["annual_incidence"]

            df_out = pd.DataFrame(
                {
                    "Year": years_out,
                    "Baseline_inc_count": base_cases,
                    "Intervention_inc_count": int_cases,
                    "Baseline_inc_per100k": base_cases * 100000.0 / total_pop,
                    "Intervention_inc_per100k": int_cases * 100000.0 / total_pop,
                }
            )

            df_out["Cases_averted_count"] = (
                df_out["Baseline_inc_count"] - df_out["Intervention_inc_count"]
            )
            df_out["Cases_averted_per100k"] = (
                df_out["Baseline_inc_per100k"] - df_out["Intervention_inc_per100k"]
            )

            df_out = df_out.round(1)

            st.success("Simulation complete.")

            st.markdown(
                "**Annual TB incidence**  \n"
                "Shown as both a rate (per 100,000) and as the absolute number of cases "
                "in the selected population."
            )

            st.subheader("📈 Annual Incidence per 100,000")
            st.line_chart(
                df_out.set_index("Year")[
                    ["Baseline_inc_per100k", "Intervention_inc_per100k"]
                ]
            )

            st.subheader("📈 Annual Incidence (counts)")
            st.line_chart(
                df_out.set_index("Year")[
                    ["Baseline_inc_count", "Intervention_inc_count"]
                ]
            )

            st.subheader("🔍 Cases Averted")
            st.write(df_out[["Year", "Cases_averted_count", "Cases_averted_per100k"]])

        except Exception as e:
            st.error(f"Dynamic model failed: {e}")
