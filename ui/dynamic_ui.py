import streamlit as st
import pandas as pd
import numpy as np
import altair as alt

# SciPy is recommended for calibration. If unavailable, RW refinement falls back to constant beta series.
try:
    from scipy.optimize import minimize_scalar, minimize

    SCIPY_AVAILABLE = True
except Exception:
    minimize_scalar = None
    minimize = None
    SCIPY_AVAILABLE = False

from engine.dynamic.exec_dynamic import run_dynamic_model
from engine.infection_backcast import (
    calc_ari_from_incidence,
    infection_prob_by_age_split,
)

# =====================================================
# Hard-coded calibration + display configuration
# =====================================================
BASELINE_DIAG_MONTHS = 12.0
INCIDENCE_FLOOR = 0.1
ARI_FLOOR = 1e-6

CALIB_YEARS_FIT = 20  # fit window length
CALIB_YEARS_SHOW = 10  # show only last 10 years

# Random-walk beta calibration (hard-coded; no UI controls)
BETA_RW_PCT = 10
BETA_RW_WEIGHT = 0.005  # you found this is needed given objective scaling
BETA_BOUNDS = (0.01, 50.0)

# Wider ARI adjustment bounds (improves constant/falling fits)
ARI_ADJ_BOUNDS = (0.05, 5.0)
ARI_ADJ_GRID_POINTS = 21


# =====================================================
# Streamlit/Altair compatibility helper
# =====================================================
def show_altair(chart):
    try:
        st.altair_chart(chart, width="stretch")  # newer Streamlit
    except TypeError:
        st.altair_chart(chart, use_container_width=True)  # older Streamlit


# =====================================================
# Utility: stable hashes for session invalidation
# =====================================================
# =====================================================
# Default fallback age distribution
# =====================================================
def default_age_distribution():
    ages = list(range(0, 101))
    prop = np.array([1 / 101] * 101)
    return pd.DataFrame({"AgeGroup": ages, "Proportion": prop})


# =====================================================
# Load country-specific population structure (OWID CSV)
# =====================================================
def load_population_data(
    country_code="AUS", file_path="data/population_age_latest.csv"
):
    df = pd.read_csv(file_path)
    mask = df["iso_code"].str.upper() == country_code.upper()
    df_country = df.loc[mask].copy()

    if df_country.empty:
        st.warning(f"No population structure for {country_code}. Using default global.")
        return (
            default_age_distribution(),
            pd.DataFrame({"age": range(0, 101), "population": [1] * 101}),
        )

    total_pop_country = df_country["population"].sum()
    df_country["Proportion"] = df_country["population"] / total_pop_country

    # Build 5-year bins for display
    bin_edges = list(range(0, 105, 5))
    bin_labels = [
        f"{bin_edges[i]}–{bin_edges[i+1]-1}" for i in range(len(bin_edges) - 1)
    ]
    bin_labels.append("100+")

    df_country["AgeBin"] = pd.cut(
        df_country["age"], bins=bin_edges + [200], labels=bin_labels, right=False
    )

    # explicit observed= to silence pandas FutureWarning
    df_age_groups = df_country.groupby("AgeBin", as_index=False, observed=False)[
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
    Returns dict: offset -> incidence per 100k, where 0 = present, -1 = last year, ...
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

                # avoid zeros
                incs = np.maximum(incs, INCIDENCE_FLOOR)

                year_min = int(years[0])
                year_max = int(years[-1])

                inc_min = float(np.min(incs))
                inc_max = float(np.max(incs))

                # geometric mean trend
                ratios = []
                for i in range(1, len(incs)):
                    if incs[i - 1] > INCIDENCE_FLOOR and incs[i] > INCIDENCE_FLOOR:
                        ratios.append(incs[i] / incs[i - 1])
                trend = float(np.exp(np.mean(np.log(ratios)))) if ratios else 1.0

                inc_map = dict(zip(years, incs))

                # define "present" as last year in CSV
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
                        # forward extrapolation
                        j = target_year - year_max
                        extrap = incs[-1] * (trend**j)
                        inc_hist[-k] = float(min(extrap, inc_max))
                    else:
                        # backward extrapolation; never below min observed
                        j = year_min - target_year
                        extrap = incs[0] * (trend ** (-j))
                        inc_hist[-k] = float(max(extrap, inc_min))
    else:
        inc_hist = {-k: float(user_incidence) for k in years_back}

    # Smooth (3-year moving average)
    inc_series = pd.Series(inc_hist).sort_index()
    inc_series = inc_series.rolling(window=3, center=True, min_periods=1).mean()

    # Final floor
    inc_hist = {
        int(k): float(max(v, INCIDENCE_FLOOR)) for k, v in inc_series.to_dict().items()
    }
    return inc_hist


def compute_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
    """
    Compute LTBI probabilities for a reference time shift_years in the past.
    ari_adjustment scales incidence→ARI mapping, and therefore LTBI.
    """
    max_age = int(max(ages))
    min_key = int(min(inc_hist.keys()))

    inc_ref = {}
    for k in range(0, max_age + 1):
        src_key = -(k + shift_years)
        inc_ref[-k] = float(
            inc_hist.get(src_key, inc_hist.get(min_key, INCIDENCE_FLOOR))
        )

    # calc_ari_from_incidence supports adjustment= in your updated engine;
    # keep a fallback for older versions.
    try:
        ari_hist = calc_ari_from_incidence(inc_ref, adjustment=float(ari_adjustment))
    except TypeError:
        ari_hist = calc_ari_from_incidence(inc_ref)

    ari_hist = {t: max(float(a), ARI_FLOOR) for t, a in ari_hist.items()}
    ltbi_ever, ltbi_recent, ltbi_remote = infection_prob_by_age_split(ages, ari_hist)
    return ltbi_ever, ltbi_recent, ltbi_remote


def minimize_scalar_bounded_grid(fn, lo, hi, n=41, refine_steps=4):
    """SciPy-free bounded minimiser: coarse-to-fine grid search."""
    best_x = None
    best_f = float("inf")

    for _ in range(refine_steps):
        xs = np.linspace(lo, hi, n)
        fs = np.array([fn(x) for x in xs], dtype=float)

        idx = int(np.argmin(fs))
        best_x = float(xs[idx])
        best_f = float(fs[idx])

        if idx == 0:
            lo, hi = float(xs[0]), float(xs[1])
        elif idx == n - 1:
            lo, hi = float(xs[-2]), float(xs[-1])
        else:
            lo, hi = float(xs[idx - 1]), float(xs[idx + 1])

    return best_x, best_f


def calibrate_beta_and_ltbi_scale(
    age_counts,
    ages,
    inc_hist,
    calib_years,
    risk_inputs,
    pre_det_months,
    delta_pre,
    beta_bounds=(0.0, 50.0),
    adj_bounds=(0.2, 2.0),
    adj_grid_points=13,
):
    """
    Jointly calibrate:
      - beta (transmission)
      - ari_adjustment (scales incidence→ARI→LTBI mapping)

    Uses nested optimisation:
      - grid over ari_adjustment
      - within each, optimise beta
    """
    total_pop = float(sum(age_counts.values()))

    # Observed incidence per 100k over calibration window (oldest -> newest)
    obs = np.array([inc_hist[-k] for k in range(calib_years, 0, -1)], dtype=float)

    # Seed incidence at start of window (used only to initialise I0)
    inc0 = float(inc_hist.get(-calib_years, obs[0]))

    best = {"rmse": float("inf"), "beta": None, "adj": None, "fit": None}

    adj_values = np.linspace(adj_bounds[0], adj_bounds[1], adj_grid_points)

    for adj in adj_values:
        # LTBI at start of calibration window using this adjustment
        ltbi_ever0, ltbi_recent0, _ = compute_ltbi_from_inc_hist(
            ages, inc_hist, shift_years=calib_years, ari_adjustment=float(adj)
        )

        base_params = {
            "age_counts": age_counts,
            "ltbi_ever": ltbi_ever0,
            "ltbi_recent": ltbi_recent0,
            "initial_incidence_per_100k": inc0,
            "pre_det_months": float(pre_det_months),
            "delta_pre": float(delta_pre),
            "delta_post": float(delta_pre),
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
            model_inc_per100k = (
                np.array(sim["annual_incidence"], dtype=float) * 100000.0 / total_pop
            )
            err = model_inc_per100k - obs
            return float(np.sqrt(np.mean(err**2)))

        if SCIPY_AVAILABLE:
            res = minimize_scalar(objective, bounds=beta_bounds, method="bounded")
            beta_hat = float(res.x)
            rmse = float(res.fun)
        else:
            beta_hat, rmse = minimize_scalar_bounded_grid(
                objective, beta_bounds[0], beta_bounds[1]
            )

        if rmse < best["rmse"]:
            p_final = dict(base_params)
            p_final["beta"] = beta_hat
            sim_fit = run_dynamic_model(p_final, years=calib_years, intervention=False)
            fit_inc_per100k = (
                np.array(sim_fit["annual_incidence"], dtype=float)
                * 100000.0
                / total_pop
            )

            best.update(
                {
                    "rmse": rmse,
                    "beta": beta_hat,
                    "adj": float(adj),
                    "fit": fit_inc_per100k,
                }
            )

    return best["beta"], best["adj"], best["rmse"], obs, best["fit"]


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
        "• Transmission (β) and LTBI scale are calibrated to historical incidence"
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
        help="β and LTBI scale are fitted so that modelled incidence matches the historical incidence pattern "
        "over this many past years.",
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
    # Diagnosis improvement (intervention only)
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

    # Scale to user-selected population
    total_pop_country = float(df_country["population"].sum())
    age_counts = {
        int(row["age"]): float(population)
        * (float(row["population"]) / total_pop_country)
        for _, row in df_country.iterrows()
    }
    ages = sorted(age_counts.keys())
    max_age = int(max(ages))

    # Build incidence history deep enough for calibration + LTBI
    years_back_needed = max_age + calib_years + 5
    inc_hist = build_incidence_history(
        hist_pattern, user_incidence, years_back_needed, uploaded_inc_df=uploaded_inc_df
    )

    st.info(
        "LTBI-by-age will be displayed after calibration (click **Run Dynamic Simulation**)."
    )

    # --------------------------------------------------
    # Run simulations
    # --------------------------------------------------
    if st.sidebar.button("Run Dynamic Simulation"):
        st.info("Calibrating β and LTBI scale to historical incidence…")

        beta_hat, ari_adj_hat, rmse, obs_inc, fit_inc = calibrate_beta_and_ltbi_scale(
            age_counts=age_counts,
            ages=ages,
            inc_hist=inc_hist,
            calib_years=calib_years,
            risk_inputs=risk_inputs,
            pre_det_months=pre_det_months,
            delta_pre=delta_pre,
            beta_bounds=(0.0, 50.0),
            adj_bounds=(0.2, 2.0),
            adj_grid_points=13,
        )

        st.success(
            f"Calibration complete: β={beta_hat:.2f}, ARI adjustment={ari_adj_hat:.2f} "
            f"(RMSE={rmse:.2f} per 100k)"
        )

        # Recompute LTBI-at-now using calibrated adjustment
        ltbi_ever_cal, ltbi_recent_cal, _ = compute_ltbi_from_inc_hist(
            ages, inc_hist, shift_years=0, ari_adjustment=ari_adj_hat
        )

        # --- LTBI BY AGE STACKED CHART (≤ 60 YEARS), CALIBRATED ---
        st.subheader("📉 LTBI Prevalence by Age (after calibration, ages 0–60)")

        ltbi_age_df = pd.DataFrame(
            {
                "Age": ages,
                "LTBI_recent": 100 * pd.Series(ltbi_recent_cal),
                "LTBI_remote": 100
                * (pd.Series(ltbi_ever_cal) - pd.Series(ltbi_recent_cal)),
            }
        )
        ltbi_age_df = ltbi_age_df[ltbi_age_df["Age"] <= 60].copy()
        ltbi_age_df = ltbi_age_df.melt(
            id_vars="Age", var_name="Type", value_name="Percent"
        )

        ltbi_chart = (
            alt.Chart(ltbi_age_df)
            .mark_area()
            .encode(
                x="Age:Q",
                y="Percent:Q",
                color="Type:N",
                tooltip=["Age", "Type", alt.Tooltip("Percent:Q", format=".1f")],
            )
        )
        show_altair(ltbi_chart)

        st.caption(
            "Calibrated LTBI by age (≤60 years).\n\n"
            "• Recent LTBI = infection acquired in the last 5 years\n"
            "• Remote LTBI = infection acquired more than 5 years ago\n"
            "• The stacked height equals total LTBI prevalence at each age\n\n"
            "This uses the calibrated ARI adjustment."
        )

        # Calibration fit plot (relative years: -calib_years .. -1)
        years_axis = list(range(-calib_years, 0))
        df_cal = pd.DataFrame(
            {"Year": years_axis, "Observed": obs_inc, "Model fit": fit_inc}
        ).melt(id_vars="Year", var_name="Series", value_name="Incidence_per100k")

        st.subheader("🧪 Calibration fit to historical incidence")
        cal_chart = (
            alt.Chart(df_cal)
            .mark_line()
            .encode(
                x=alt.X("Year:Q", title="Years relative to present (past < 0)"),
                y=alt.Y("Incidence_per100k:Q", title="Incidence per 100,000 per year"),
                color="Series:N",
                tooltip=[
                    "Year",
                    "Series",
                    alt.Tooltip("Incidence_per100k:Q", format=".1f"),
                ],
            )
        )
        show_altair(cal_chart)

        # -------------------------
        # Build baseline + intervention params at "now"
        # -------------------------
        params_base = {}
        params_int = {}

        for p in (params_base, params_int):
            p["beta"] = beta_hat
            p.update(risk_inputs)

            p["ltbi_ever"] = ltbi_ever_cal
            p["ltbi_recent"] = ltbi_recent_cal
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

            base_cases = np.array(baseline["annual_incidence"], dtype=float)
            int_cases = np.array(intervention["annual_incidence"], dtype=float)

            years_out = baseline.get("annual_incidence_time")
            if years_out is None:
                years_out = np.arange(0, len(base_cases))

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

            # Combined “past fit + future projection” on one chart
            df_future = df_out[
                ["Year", "Baseline_inc_per100k", "Intervention_inc_per100k"]
            ].copy()
            df_future = df_future.melt(
                id_vars="Year", var_name="Series", value_name="Incidence_per100k"
            )
            df_future["Series"] = df_future["Series"].replace(
                {
                    "Baseline_inc_per100k": "Baseline (projected)",
                    "Intervention_inc_per100k": "Intervention (projected)",
                }
            )

            df_all = pd.concat([df_cal, df_future], ignore_index=True)

            rule = (
                alt.Chart(pd.DataFrame({"Year": [0.0]}))
                .mark_rule(strokeDash=[4, 4])
                .encode(x="Year:Q")
            )
            combo = (
                alt.Chart(df_all)
                .mark_line()
                .encode(
                    x=alt.X(
                        "Year:Q",
                        title="Years relative to present (past < 0, future ≥ 0)",
                    ),
                    y=alt.Y(
                        "Incidence_per100k:Q", title="Incidence per 100,000 per year"
                    ),
                    color="Series:N",
                    tooltip=[
                        "Year",
                        "Series",
                        alt.Tooltip("Incidence_per100k:Q", format=".1f"),
                    ],
                )
            )

            st.subheader("📈 Incidence: fitted history + projection")
            show_altair(combo + rule)

        except Exception as e:
            st.error(f"Dynamic model failed: {e}")
