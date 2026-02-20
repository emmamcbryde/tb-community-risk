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
CALIB_YEARS_SHOW = 10  # show only last N years on charts

# Random-walk beta calibration (hard-coded; no UI controls)
BETA_RW_PCT = 10
BETA_RW_WEIGHT = 0.005  # penalty weight (relative to data misfit)
BETA_BOUNDS = (0.01, 50.0)

# ARI adjustment bounds (scales incidence -> ARI -> LTBI)
ARI_ADJ_BOUNDS = (0.05, 5.0)
ARI_ADJ_GRID_POINTS = 21


# =====================================================
# Streamlit/Altair compatibility helper
# =====================================================
def show_altair(chart):
    """Compatible with both older and newer Streamlit versions."""
    try:
        st.altair_chart(chart, width="stretch")  # newer Streamlit
    except TypeError:
        st.altair_chart(chart, use_container_width=True)  # older Streamlit


def show_df(df):
    """Compatible dataframe width helper."""
    try:
        st.dataframe(df, width="stretch")
    except TypeError:
        st.dataframe(df, use_container_width=True)


# =====================================================
# Utility: stable hashes for session invalidation
# =====================================================
def hash_df(df, cols=None):
    if df is None:
        return None
    if cols is not None:
        df = df[cols].copy()
    return int(pd.util.hash_pandas_object(df, index=False).sum())


def clear_calibration():
    for k in [
        "cal_done",
        "cal_sig",
        "cal_beta_forward",
        "cal_beta_series",
        "cal_ari_adj",
        "cal_rmse_rw",
        "cal_obs_inc",
        "cal_fit_inc",
        "cal_ref_year",
        "cal_ltbi_ever_now",
        "cal_ltbi_recent_now",
        "cal_state_now",
    ]:
        st.session_state.pop(k, None)
    clear_simulation()


def clear_simulation():
    for k in ["sim_sig", "sim_df_future"]:
        st.session_state.pop(k, None)


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

    total_pop_country = float(df_country["population"].sum())
    df_country["Proportion"] = df_country["population"] / total_pop_country

    # 5-year bins for display
    bin_edges = list(range(0, 105, 5))
    bin_labels = [
        f"{bin_edges[i]}–{bin_edges[i+1]-1}" for i in range(len(bin_edges) - 1)
    ]
    bin_labels.append("100+")

    df_country["AgeBin"] = pd.cut(
        df_country["age"], bins=bin_edges + [200], labels=bin_labels, right=False
    )

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
                    "Incidence CSV must have columns: year, incidence. Using constant."
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

                # geometric trend
                ratios = []
                for i in range(1, len(incs)):
                    if incs[i - 1] > INCIDENCE_FLOOR and incs[i] > INCIDENCE_FLOOR:
                        ratios.append(incs[i] / incs[i - 1])
                trend = float(np.exp(np.mean(np.log(ratios)))) if ratios else 1.0

                inc_map = dict(zip(years, incs))
                ref_year = year_max  # "present" = last year in CSV

                inc_hist = {}
                for k in years_back:
                    target_year = ref_year - k

                    if year_min <= target_year <= year_max:
                        nearest = min(
                            inc_map.keys(), key=lambda y: abs(y - target_year)
                        )
                        inc_hist[-k] = float(inc_map[nearest])

                    elif target_year > year_max:
                        # forward extrapolation (missing recent years)
                        j = target_year - year_max
                        extrap = incs[-1] * (trend**j)
                        inc_hist[-k] = float(min(extrap, inc_max))

                    else:
                        # backward extrapolation; never below minimum observed
                        j = year_min - target_year
                        extrap = incs[0] * (trend ** (-j))
                        inc_hist[-k] = float(max(extrap, inc_min))
    else:
        inc_hist = {-k: float(user_incidence) for k in years_back}

    # smooth with 3-year moving average
    inc_series = pd.Series(inc_hist).sort_index()
    inc_series = inc_series.rolling(window=3, center=True, min_periods=1).mean()

    inc_hist = {
        int(k): float(max(v, INCIDENCE_FLOOR)) for k, v in inc_series.to_dict().items()
    }
    return inc_hist


def compute_ltbi_from_inc_hist(ages, inc_hist, shift_years=0, ari_adjustment=1.0):
    """
    Compute LTBI probabilities at a reference time shift_years in the past.
    ari_adjustment scales incidence→ARI mapping (Stýblo-style), and therefore LTBI.
    """
    max_age = int(max(ages))
    min_key = int(min(inc_hist.keys()))

    inc_ref = {}
    for k in range(0, max_age + 1):
        src_key = -(k + shift_years)
        inc_ref[-k] = float(
            inc_hist.get(src_key, inc_hist.get(min_key, INCIDENCE_FLOOR))
        )

    try:
        ari_hist = calc_ari_from_incidence(inc_ref, adjustment=float(ari_adjustment))
    except TypeError:
        ari_hist = calc_ari_from_incidence(inc_ref)

    ari_hist = {t: max(float(a), ARI_FLOOR) for t, a in ari_hist.items()}
    ltbi_ever, ltbi_recent, ltbi_remote = infection_prob_by_age_split(ages, ari_hist)
    return ltbi_ever, ltbi_recent, ltbi_remote


def calibrate_beta_and_ltbi_scale(
    age_counts,
    ages,
    inc_hist,
    calib_years,
    risk_inputs,
    pre_det_months,
    delta_pre,
    beta_bounds=BETA_BOUNDS,
    adj_bounds=ARI_ADJ_BOUNDS,
    adj_grid_points=ARI_ADJ_GRID_POINTS,
):
    """
    Coarse calibration:
      - grid search ARI adjustment (scales incidence -> ARI -> LTBI),
      - for each, optimise a constant beta to fit incidence over `calib_years`.
    Returns:
      beta_hat, ari_adj_hat, rmse_hat, obs_inc_per100k (oldest -> newest)
    """
    total_pop = float(sum(age_counts.values()))

    # Observed incidence per 100k over calibration window (oldest -> newest)
    obs = np.array([inc_hist[-k] for k in range(calib_years, 0, -1)], dtype=float)

    # Seed incidence at start of window
    inc0 = float(inc_hist.get(-calib_years, obs[0]))

    best = {"rmse": float("inf"), "beta": None, "adj": None}
    adj_values = np.linspace(adj_bounds[0], adj_bounds[1], adj_grid_points)
    beta_min, beta_max = float(beta_bounds[0]), float(beta_bounds[1])

    for adj in adj_values:
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

        def rmse_for_beta(beta):
            p = dict(base_params)
            p["beta"] = float(beta)
            sim = run_dynamic_model(p, years=int(calib_years), intervention=False)
            model_inc_per100k = (
                np.array(sim["annual_incidence"], dtype=float) * 100000.0 / total_pop
            )
            err = model_inc_per100k - obs
            return float(np.sqrt(np.mean(err**2)))

        if SCIPY_AVAILABLE and minimize_scalar is not None:
            res = minimize_scalar(
                rmse_for_beta, bounds=(beta_min, beta_max), method="bounded"
            )
            beta_hat = float(res.x)
            rmse = float(res.fun)
        else:
            grid = np.linspace(beta_min, beta_max, 31)
            vals = [rmse_for_beta(b) for b in grid]
            j = int(np.argmin(vals))
            beta_hat = float(grid[j])
            rmse = float(vals[j])

        if rmse < best["rmse"]:
            best.update({"rmse": rmse, "beta": beta_hat, "adj": float(adj)})

    return best["beta"], best["adj"], best["rmse"], obs


def refine_beta_random_walk(
    age_counts,
    ages,
    inc_hist,
    calib_years,
    risk_inputs,
    pre_det_months,
    delta_pre,
    ari_adjustment,
    beta_init,
    beta_bounds=BETA_BOUNDS,
):
    """
    Refine beta from a scalar into a smooth beta(t) over the calibration window.
    Returns:
      beta_series_hat (length calib_years)
    """
    beta_min, beta_max = float(beta_bounds[0]), float(beta_bounds[1])
    beta_init = float(np.clip(beta_init, beta_min, beta_max))

    # If SciPy isn't available, fall back to constant beta series
    if not (SCIPY_AVAILABLE and minimize is not None) or calib_years < 2:
        return np.full(calib_years, beta_init, dtype=float)

    total_pop = float(sum(age_counts.values()))
    obs = np.array([inc_hist[-k] for k in range(calib_years, 0, -1)], dtype=float)
    obs_scale = max(float(np.mean(obs)), 1.0)

    inc0 = float(inc_hist.get(-calib_years, obs[0]))
    sigma_rw = float(np.log(1.0 + BETA_RW_PCT / 100.0))

    ltbi_ever0, ltbi_recent0, _ = compute_ltbi_from_inc_hist(
        ages, inc_hist, shift_years=calib_years, ari_adjustment=float(ari_adjustment)
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

    x0 = np.full(calib_years, np.log(beta_init), dtype=float)
    bounds = [(np.log(beta_min), np.log(beta_max))] * calib_years

    def simulate_per100k_from_x(x):
        beta_series = np.exp(x)
        p = dict(base_params)
        p["beta_series"] = np.asarray(beta_series, dtype=float)
        p["beta"] = float(beta_series[-1])  # scalar fallback
        sim = run_dynamic_model(p, years=int(calib_years), intervention=False)
        pred = np.array(sim["annual_incidence"], dtype=float) * 100000.0 / total_pop
        return pred

    def objective(x):
        pred = simulate_per100k_from_x(x)
        data = float(np.mean(((pred - obs) / obs_scale) ** 2))
        dx = np.diff(x)
        smooth = float(np.mean((dx / sigma_rw) ** 2))
        return data + BETA_RW_WEIGHT * smooth

    res = minimize(
        objective, x0, method="L-BFGS-B", bounds=bounds, options={"maxiter": 120}
    )
    x_hat = np.asarray(res.x, dtype=float)
    beta_series_hat = np.exp(x_hat)

    return beta_series_hat


def _get_annual(sim):
    if "annual_incidence" in sim:
        y = np.array(sim["annual_incidence"], dtype=float)
    elif "incidence" in sim:
        y = np.array(sim["incidence"], dtype=float)
    else:
        raise KeyError(f"Simulation output missing incidence keys: {list(sim.keys())}")

    t = sim.get("annual_incidence_time", None)
    if t is None:
        t = np.arange(0, len(y), dtype=float)
    else:
        t = np.array(t, dtype=float)

    return t, y


# =====================================================
# Main UI
# =====================================================
def render_dynamic_ui():
    st.header("📈 Dynamic LTBI → TB Model")

    st.info(
        "Workflow:\n"
        "1) **Calibrate** (fits β(t) and ARI adjustment; shows backcast + LTBI by age today)\n"
        "2) Choose intervention settings\n"
        "3) **Simulate** (projects baseline vs intervention)\n"
    )

    # -------------------------
    # Inputs that affect calibration
    # -------------------------
    st.sidebar.subheader("Core inputs")
    population = st.sidebar.number_input("Population size", min_value=50, value=10000)
    user_incidence = st.sidebar.number_input(
        "Baseline annual incidence (per 100k)", 0, 500, 30
    )
    time_horizon = st.sidebar.slider("Projection horizon (years)", 1, 30, 20)

    st.sidebar.subheader("Risk factors")
    smoker_pct = st.sidebar.slider("Smoker (%)", 0, 100, 30)
    alcohol_pct = st.sidebar.slider("Excess alcohol use (%)", 0, 100, 15)
    diabetes_pct = st.sidebar.slider("Diabetes (%)", 0, 100, 10)
    renal_pct = st.sidebar.slider("Renal impairment (%)", 0, 100, 5)
    HIV_treated_pct = st.sidebar.slider("HIV treated (%)", 0, 100, 3)
    HIV_untreated_pct = st.sidebar.slider("HIV untreated (%)", 0, 100, 3)

    risk_inputs = {
        "smoker_pct": smoker_pct,
        "alcohol_pct": alcohol_pct,
        "diabetes_pct": diabetes_pct,
        "renal_pct": renal_pct,
        "HIV_treated_pct": HIV_treated_pct,
        "HIV_untreated_pct": HIV_untreated_pct,
    }

    st.sidebar.subheader("Historical incidence pattern")
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
    ref_year = None
    if hist_pattern == "Upload CSV (year, incidence)":
        inc_file = st.sidebar.file_uploader("Upload incidence CSV", type="csv")
        if inc_file:
            tmp = pd.read_csv(inc_file)
            if {"year", "incidence"}.issubset(tmp.columns):
                uploaded_inc_df = tmp.sort_values("year")
                ref_year = int(uploaded_inc_df["year"].max())
                st.sidebar.success("Incidence history loaded.")
            else:
                st.sidebar.error("CSV must contain columns: year, incidence")

    st.sidebar.subheader("Age distribution")
    age_method = st.sidebar.radio(
        "Choose method:",
        ["Country ISO code (recommended)", "Upload custom CSV", "Default global"],
    )

    age_df_display = None
    df_country = None
    country = None
    age_upload_hash = None

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
                age_upload_hash = hash_df(df[["AgeGroup", "Proportion"]])
                age_df_display = df
                df_country = pd.DataFrame(
                    {
                        "age": df["AgeGroup"].astype(int),
                        "population": df["Proportion"] * population,
                    }
                )
            else:
                st.sidebar.error(
                    "CSV must include AgeGroup and Proportion. Using default."
                )
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

    # Show age distribution
    st.subheader("📊 Age Distribution (5-year bins)")
    show_df(age_df_display)

    # Scale age distribution to chosen population
    total_pop_country = float(df_country["population"].sum())
    age_counts = {
        int(row["age"]): float(population)
        * (float(row["population"]) / total_pop_country)
        for _, row in df_country.iterrows()
    }
    ages = sorted(age_counts.keys())
    max_age = int(max(ages))
    total_pop = float(sum(age_counts.values()))

    # Build incidence history deep enough for calibration + LTBI backcast
    years_back_needed = max_age + CALIB_YEARS_FIT + 5
    inc_hist = build_incidence_history(
        hist_pattern,
        user_incidence,
        years_back_needed,
        uploaded_inc_df=uploaded_inc_df,
    )

    # -------------------------
    # Invalidate calibration if upstream inputs changed
    # -------------------------
    inc_hash = (
        hash_df(uploaded_inc_df, cols=["year", "incidence"])
        if uploaded_inc_df is not None
        else None
    )

    cal_sig = (
        int(population),
        float(user_incidence),
        hist_pattern,
        inc_hash,
        age_method,
        (country or ""),
        age_upload_hash,
        tuple(sorted((k, float(v)) for k, v in risk_inputs.items())),
    )

    if st.session_state.get("cal_sig", None) != cal_sig:
        clear_calibration()

    # -------------------------
    # Calibration button
    # -------------------------
    st.sidebar.markdown("---")
    calibrate_clicked = st.sidebar.button("1) Calibrate")

    if calibrate_clicked:
        st.info("Running calibration…")

        pre_det_months = BASELINE_DIAG_MONTHS
        delta_pre = 12.0 / pre_det_months

        # Fit beta (scalar) + ARI adjustment
        beta_hat, ari_adj_hat, _, obs_inc = calibrate_beta_and_ltbi_scale(
            age_counts=age_counts,
            ages=ages,
            inc_hist=inc_hist,
            calib_years=CALIB_YEARS_FIT,
            risk_inputs=risk_inputs,
            pre_det_months=pre_det_months,
            delta_pre=delta_pre,
        )

        # Refine beta into a random-walk series
        beta_series_hat = refine_beta_random_walk(
            age_counts=age_counts,
            ages=ages,
            inc_hist=inc_hist,
            calib_years=CALIB_YEARS_FIT,
            risk_inputs=risk_inputs,
            pre_det_months=pre_det_months,
            delta_pre=delta_pre,
            ari_adjustment=ari_adj_hat,
            beta_init=beta_hat,
        )
        beta_series_hat = np.asarray(beta_series_hat, dtype=float)
        beta_forward = float(beta_series_hat[-1])

        # Run ONE consistent backcast with the fitted beta_series to:
        #  - compute fit_inc (per 100k)
        #  - capture the end-of-calibration state to start projections (avoids a jump)
        inc0 = float(inc_hist.get(-CALIB_YEARS_FIT, obs_inc[0]))
        ltbi_ever0, ltbi_recent0, _ = compute_ltbi_from_inc_hist(
            ages,
            inc_hist,
            shift_years=CALIB_YEARS_FIT,
            ari_adjustment=float(ari_adj_hat),
        )

        p_cal = {
            "beta": beta_forward,  # scalar fallback
            "beta_series": beta_series_hat,  # used during calibration years
            "age_counts": age_counts,
            "ltbi_ever": ltbi_ever0,
            "ltbi_recent": ltbi_recent0,
            "initial_incidence_per_100k": float(inc0),
            "pre_det_months": float(pre_det_months),
            "delta_pre": float(delta_pre),
            "delta_post": float(delta_pre),
            "ltbi_coverage": 0.0,
            "rollout_years": 0,
            "treatment_method": "None",
            "testing_method": "None",
        }
        p_cal.update(risk_inputs)

        sim_cal = run_dynamic_model(
            p_cal, years=int(CALIB_YEARS_FIT), intervention=False
        )
        fit_inc = (
            np.array(sim_cal["annual_incidence"], dtype=float) * 100000.0 / total_pop
        )
        rmse_rw = float(np.sqrt(np.mean((fit_inc - np.asarray(obs_inc)) ** 2)))

        # End-of-calibration state (present): THIS is what we will start projections from.
        state_now = {
            "S": float(np.asarray(sim_cal["S"], dtype=float)[-1]),
            "L_fast": float(np.asarray(sim_cal["L_fast"], dtype=float)[-1]),
            "L_slow": float(np.asarray(sim_cal["L_slow"], dtype=float)[-1]),
            "I": float(np.asarray(sim_cal["I"], dtype=float)[-1]),
            "R": float(np.asarray(sim_cal["R"], dtype=float)[-1]),
        }

        # LTBI by age today (calibrated ARI adjustment) - for display only
        ltbi_ever_now, ltbi_recent_now, _ = compute_ltbi_from_inc_hist(
            ages,
            inc_hist,
            shift_years=0,
            ari_adjustment=float(ari_adj_hat),
        )

        st.session_state["cal_done"] = True
        st.session_state["cal_sig"] = cal_sig
        st.session_state["cal_beta_forward"] = beta_forward
        st.session_state["cal_beta_series"] = beta_series_hat
        st.session_state["cal_ari_adj"] = float(ari_adj_hat)
        st.session_state["cal_rmse_rw"] = rmse_rw
        st.session_state["cal_obs_inc"] = np.asarray(obs_inc, dtype=float)
        st.session_state["cal_fit_inc"] = np.asarray(fit_inc, dtype=float)
        st.session_state["cal_ref_year"] = ref_year
        st.session_state["cal_ltbi_ever_now"] = ltbi_ever_now
        st.session_state["cal_ltbi_recent_now"] = ltbi_recent_now
        st.session_state["cal_state_now"] = state_now

        clear_simulation()

    # -------------------------
    # Display calibration outputs (persist after reruns)
    # -------------------------
    if st.session_state.get("cal_done", False):
        beta_forward = float(st.session_state["cal_beta_forward"])
        beta_series_hat = np.asarray(st.session_state["cal_beta_series"], dtype=float)
        ari_adj_hat = float(st.session_state["cal_ari_adj"])
        rmse_rw = float(st.session_state["cal_rmse_rw"])
        obs_inc = np.asarray(st.session_state["cal_obs_inc"], dtype=float)
        fit_inc = np.asarray(st.session_state["cal_fit_inc"], dtype=float)
        ref_year_used = st.session_state.get("cal_ref_year", None)

        st.subheader("🧪 Calibration results")

        st.success(
            f"Calibrated over {CALIB_YEARS_FIT} years (showing last {CALIB_YEARS_SHOW}). "
            f"β(t) range {beta_series_hat.min():.2f}–{beta_series_hat.max():.2f}. "
            f"Using β={beta_forward:.2f} for projections. "
            f"ARI adjustment={ari_adj_hat:.2f}. RMSE={rmse_rw:.2f} per 100k."
        )

        show_years = min(CALIB_YEARS_SHOW, CALIB_YEARS_FIT)
        obs_show = obs_inc[-show_years:]
        fit_show = fit_inc[-show_years:]

        if ref_year_used is not None:
            years_axis = list(
                range(int(ref_year_used) - show_years, int(ref_year_used))
            )
            x_title = "Calendar year"
        else:
            years_axis = list(range(-show_years, 0))
            x_title = "Years before present (relative)"

        df_cal = pd.DataFrame(
            {"Year": years_axis, "Observed": obs_show, "Model fit": fit_show}
        ).melt(id_vars="Year", var_name="Series", value_name="Incidence_per100k")

        cal_chart = (
            alt.Chart(df_cal)
            .mark_line()
            .encode(
                x=alt.X("Year:Q", title=x_title),
                y=alt.Y("Incidence_per100k:Q", title="Incidence per 100,000 per year"),
                color="Series:N",
                tooltip=[
                    alt.Tooltip("Year:Q", format=".0f"),
                    "Series:N",
                    alt.Tooltip("Incidence_per100k:Q", format=".1f"),
                ],
            )
        )
        show_altair(cal_chart)

        # LTBI by age today (calibrated)
        st.subheader("📉 LTBI prevalence by age today (after calibration, ages 0–60)")

        ltbi_ever_now = st.session_state["cal_ltbi_ever_now"]
        ltbi_recent_now = st.session_state["cal_ltbi_recent_now"]

        ltbi_age_df = pd.DataFrame(
            {
                "Age": ages,
                "LTBI_recent": 100 * pd.Series(ltbi_recent_now),
                "LTBI_remote": 100
                * (pd.Series(ltbi_ever_now) - pd.Series(ltbi_recent_now)),
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
                x=alt.X("Age:Q", title="Age"),
                y=alt.Y("Percent:Q", title="Percent"),
                color=alt.Color("Type:N", title=None),
                tooltip=[
                    alt.Tooltip("Age:Q", format=".0f"),
                    "Type:N",
                    alt.Tooltip("Percent:Q", format=".1f"),
                ],
            )
        )
        show_altair(ltbi_chart)

    # -------------------------
    # Intervention controls + simulate button
    # -------------------------
    st.sidebar.markdown("---")
    st.sidebar.subheader("Intervention settings (for simulation)")

    testing_method = st.sidebar.selectbox("Testing method", ["TST", "IGRA", "None"])
    treatment_method = st.sidebar.selectbox(
        "Treatment regimen", ["1HP", "3HP", "4R", "6H", "9H", "None"]
    )
    ltbi_coverage = (
        st.sidebar.slider("LTBI Test & Treat total coverage (%)", 0, 100, 50) / 100.0
    )
    rollout_years = st.sidebar.slider("Rollout duration (years)", 1, 10, 5)

    st.sidebar.subheader("Diagnosis improvement (intervention)")
    diag_reduction_pct = st.sidebar.slider(
        "Percent reduction in time before treatment (%)", 0, 100, 50
    )

    simulate_clicked = st.sidebar.button("2) Simulate")

    # simulation invalidation (if intervention settings changed)
    sim_sig = (
        st.session_state.get("cal_sig", None),
        int(time_horizon),
        testing_method,
        treatment_method,
        float(ltbi_coverage),
        int(rollout_years),
        float(diag_reduction_pct),
    )
    if st.session_state.get("sim_sig", None) != sim_sig:
        clear_simulation()

    if simulate_clicked:
        if not st.session_state.get("cal_done", False):
            st.error("Please run **1) Calibrate** first.")
        else:
            st.info("Running baseline and intervention projections…")

            beta_forward = float(st.session_state["cal_beta_forward"])
            state_now = st.session_state["cal_state_now"]

            pre_det_months = BASELINE_DIAG_MONTHS
            post_det_months = max(
                BASELINE_DIAG_MONTHS * (1.0 - diag_reduction_pct / 100.0), 0.1
            )
            delta_pre = 12.0 / pre_det_months
            delta_post = 12.0 / post_det_months

            params_base = {}
            params_int = {}

            for p in (params_base, params_int):
                p["beta"] = float(beta_forward)
                p.update(risk_inputs)
                p["age_counts"] = age_counts

                # START FROM THE END-OF-CALIBRATION STATE (prevents a visible jump at t=0)
                p["initial_state"] = dict(state_now)

                # diagnosis params
                p["pre_det_months"] = float(pre_det_months)
                p["delta_pre"] = float(delta_pre)
                p["delta_post"] = float(delta_post)

            # baseline = no intervention
            params_base["treatment_method"] = "None"
            params_base["testing_method"] = "None"
            params_base["ltbi_coverage"] = 0.0
            params_base["rollout_years"] = 0

            # intervention
            params_int["treatment_method"] = treatment_method
            params_int["testing_method"] = testing_method
            params_int["ltbi_coverage"] = float(ltbi_coverage)
            params_int["rollout_years"] = int(rollout_years)

            sim_base = run_dynamic_model(
                params_base, years=int(time_horizon), intervention=False
            )
            sim_int = run_dynamic_model(
                params_int, years=int(time_horizon), intervention=True
            )

            t_out, base_cases = _get_annual(sim_base)
            _, int_cases = _get_annual(sim_int)

            df_future = pd.DataFrame(
                {
                    "Year": t_out.astype(float),
                    "Baseline_inc_per100k": base_cases * 100000.0 / total_pop,
                    "Intervention_inc_per100k": int_cases * 100000.0 / total_pop,
                    "Baseline_inc_count": base_cases,
                    "Intervention_inc_count": int_cases,
                }
            )
            df_future["Cases_averted_count"] = (
                df_future["Baseline_inc_count"] - df_future["Intervention_inc_count"]
            )
            df_future["Cases_averted_per100k"] = (
                df_future["Baseline_inc_per100k"]
                - df_future["Intervention_inc_per100k"]
            )
            df_future = df_future.round(1)

            st.session_state["sim_sig"] = sim_sig
            st.session_state["sim_df_future"] = df_future

    # -------------------------
    # Display simulation outputs (persist after reruns)
    # -------------------------
    if st.session_state.get("sim_df_future", None) is not None and st.session_state.get(
        "cal_done", False
    ):
        df_future = st.session_state["sim_df_future"].copy()

        # Build combined plot: last CALIB_YEARS_SHOW of past fit + future baseline/intervention
        obs_inc = np.asarray(st.session_state["cal_obs_inc"], dtype=float)
        fit_inc = np.asarray(st.session_state["cal_fit_inc"], dtype=float)

        show_years = min(CALIB_YEARS_SHOW, CALIB_YEARS_FIT)
        obs_show = obs_inc[-show_years:]
        fit_show = fit_inc[-show_years:]
        years_past = np.arange(-show_years, 0).astype(float)

        df_past = pd.DataFrame(
            {
                "Year": np.concatenate([years_past, years_past]),
                "Series": ["Observed (past)"] * show_years
                + ["Model fit (past)"] * show_years,
                "Incidence_per100k": np.concatenate([obs_show, fit_show]),
            }
        )

        df_future_long = df_future.melt(
            id_vars="Year",
            value_vars=["Baseline_inc_per100k", "Intervention_inc_per100k"],
            var_name="Series",
            value_name="Incidence_per100k",
        )
        df_future_long["Series"] = df_future_long["Series"].replace(
            {
                "Baseline_inc_per100k": "Baseline (projected)",
                "Intervention_inc_per100k": "Intervention (projected)",
            }
        )

        df_all = pd.concat([df_past, df_future_long], ignore_index=True)
        df_all["Incidence_count"] = (
            df_all["Incidence_per100k"] * total_pop / 100000.0
        ).round(1)

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
                    "Year:Q", title="Years relative to present (past < 0, future ≥ 0)"
                ),
                y=alt.Y("Incidence_per100k:Q", title="Incidence per 100,000 per year"),
                color=alt.Color("Series:N", title=None),
                tooltip=[
                    alt.Tooltip("Year:Q", format=".0f"),
                    "Series:N",
                    alt.Tooltip("Incidence_per100k:Q", format=".1f"),
                    alt.Tooltip(
                        "Incidence_count:Q", format=".1f", title="Incidence (count)"
                    ),
                ],
            )
        )

        st.subheader("📈 Incidence: backcast fit + projected baseline vs intervention")
        show_altair(combo + rule)

        # Optional continuity check (useful when stress testing)
        last_fit = float(fit_inc[-1])
        first_proj = float(df_future["Baseline_inc_per100k"].iloc[0])
        st.caption(
            f"Continuity check: last fitted year (−1) = {last_fit:.1f} per 100k; "
            f"first projected year (0→1) baseline = {first_proj:.1f} per 100k."
        )

        st.subheader("📋 Projected annual outcomes (future)")
        show_df(
            df_future[
                [
                    "Year",
                    "Baseline_inc_per100k",
                    "Intervention_inc_per100k",
                    "Baseline_inc_count",
                    "Intervention_inc_count",
                    "Cases_averted_per100k",
                    "Cases_averted_count",
                ]
            ]
        )
