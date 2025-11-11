import pandas as pd
from pathlib import Path


def load_all_parameters(file_path="data/parameters.xlsx"):
    """
    Load all sheets from the parameters Excel file into a dictionary of DataFrames.
    """
    file = Path(file_path)
    if not file.exists():
        raise FileNotFoundError(f"Parameter file not found: {file.resolve()}")

    xls = pd.ExcelFile(file)
    data = {}
    for sheet in xls.sheet_names:
        try:
            df = xls.parse(sheet)
            data[sheet.lower()] = df
        except Exception as e:
            print(f"⚠️ Skipping sheet '{sheet}': {e}")
    print(f"✅ Loaded sheets: {list(data.keys())}")
    return data


def extract_core_parameters(file_path="data/parameters.xlsx"):
    """
    Extracts key TB model parameters from a complex parameters.xlsx file.
    Returns a structured dictionary suitable for simulation input.
    """

    xls = pd.ExcelFile(file_path)
    params = {}

    # --- 1. Cascade of Care ---
    if "cascade_of_care" in xls.sheet_names:
        df = xls.parse("cascade_of_care")
        if "p" in df.columns:
            params["cascade"] = df.set_index("p").to_dict(orient="index")

    # --- 2. Mortality ---
    if "mortality" in xls.sheet_names:
        df = xls.parse("mortality")
        if {"Age", "Prob"}.issubset(df.columns):
            params["mortality"] = df.groupby("Age")["Prob"].mean().to_dict()

    # --- 3. TB mortality ---
    if "tb_mortality" in xls.sheet_names:
        df = xls.parse("tb_mortality")
        if {"age", "Prob"}.issubset(df.columns):
            params["tb_mortality"] = df.groupby("age")["Prob"].mean().to_dict()

    # --- 4. Reactivation rates (RRates or rradjrates) ---
    if "rrates" in [s.lower() for s in xls.sheet_names]:
        df = xls.parse("RRates")
        keep = [c for c in ["aaa", "ysa", "Rate"] if c in df.columns]
        df = df[keep].dropna()
        params["reactivation"] = df.groupby("aaa")["Rate"].mean().to_dict()

    elif "rradjrates" in [s.lower() for s in xls.sheet_names]:
        df = xls.parse("rradjrates")
        keep = [c for c in ["aaa", "rate"] if c in df.columns]
        df = df[keep].dropna()
        params["reactivation"] = df.groupby("aaa")["rate"].mean().to_dict()

    # --- 5. SAE and SAE mortality ---
    if "sae_rate" in xls.sheet_names:
        df = xls.parse("sae_rate")
        if {"Age", "Rate"}.issubset(df.columns):
            params["sae_rate"] = df.groupby("Age")["Rate"].mean().to_dict()

    if "sae_mortality" in xls.sheet_names:
        df = xls.parse("sae_mortality")
        if {"Age", "Rate"}.issubset(df.columns):
            params["sae_mortality"] = df.groupby("Age")["Rate"].mean().to_dict()

    # --- 6. Community risk multipliers ---
    if "community_risk" in xls.sheet_names:
        df = xls.parse("community_risk")
        if {"factor", "risk_multiplier"}.issubset(df.columns):
            params["community_risk"] = df.set_index("factor")["risk_multiplier"].to_dict()

    # --- 7. Comorbidities ---
    if "comorbidities" in xls.sheet_names:
        df = xls.parse("comorbidities")
        if {"Condition", "Relative_Risk"}.issubset(df.columns):
            params["comorbidities"] = df.set_index("Condition")["Relative_Risk"].to_dict()

    # --- 8. Utilities ---
    if "utilities" in xls.sheet_names:
        df = xls.parse("utilities")
        if {"health_state", "utility"}.issubset(df.columns):
            params["utilities"] = df.set_index("health_state")["utility"].to_dict()

    print("✅ Extracted key parameter groups:", list(params.keys()))
    return params


def summarize_core_parameters(params):
    """
    Summarize extracted key parameters in a human-readable dictionary.
    """

    summary = {
        "num_cascade_params": len(params.get("cascade", {})),
        "num_age_specific_mortality": len(params.get("mortality", {})),
        "num_tb_mortality": len(params.get("tb_mortality", {})),
        "num_reactivation_groups": len(params.get("reactivation", {})),
        "num_community_risk_factors": len(params.get("community_risk", {})),
        "num_comorbidities": len(params.get("comorbidities", {})),
    }

    print("📊 Parameter summary:")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    return summary
