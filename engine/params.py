import pandas as pd
from pathlib import Path


def load_all_parameters(directory="data"):
    """
    Load all CSV parameter files in the directory into a dictionary of DataFrames.
    Keys are derived from the CSV filenames.
    """
    params = {}
    for fname in Path(directory).glob("*.csv"):
        key = fname.stem.lower()
        try:
            params[key] = pd.read_csv(fname)
        except Exception as e:
            raise RuntimeError(f"Failed to read {fname}: {e}")
    return params



def extract_core_parameters(directory="data"):
    """
    Extract key TB model parameters from multiple CSV files.
    This version contains **NO EXCEL DEPENDENCIES** and matches your renamed files.
    """

    # Load raw CSV files into a dict of DataFrames
    csv_files = {f.stem.lower(): pd.read_csv(f) for f in Path(directory).glob("*.csv")}
    params = {}

    # --- 1. Cascade of Care ---
    # Using cascade.csv
    if "cascade" in csv_files:
        df = csv_files["cascade"]
        if "p" in df.columns:
            params["cascade"] = df.set_index("p").to_dict(orient="index")

    # --- 2. Mortality ---
    if "mortality" in csv_files:
        df = csv_files["mortality"]
        if {"Age", "Prob"}.issubset(df.columns):
            params["mortality"] = df.groupby("Age")["Prob"].mean().to_dict()

    # --- 3. TB Mortality ---
    if "tb_mortality" in csv_files:
        df = csv_files["tb_mortality"]
        if {"age", "Prob"}.issubset(df.columns):
            params["tb_mortality"] = df.groupby("age")["Prob"].mean().to_dict()

    # --- 4. Reactivation rates ---
    # Handles rrates.csv or rradjrates.csv
    if "rrates" in csv_files:
        df = csv_files["rrates"]
        if {"aaa", "Rate"}.issubset(df.columns):
            params["reactivation"] = df.groupby("aaa")["Rate"].mean().to_dict()

    elif "rradjrates" in csv_files:
        df = csv_files["rradjrates"]
        if {"aaa", "rate"}.issubset(df.columns):
            params["reactivation"] = df.groupby("aaa")["rate"].mean().to_dict()

    # --- 5. SAE rates ---
    if "sae_rate" in csv_files:
        df = csv_files["sae_rate"]
        if {"Age", "Rate"}.issubset(df.columns):
            params["sae_rate"] = df.groupby("Age")["Rate"].mean().to_dict()

    if "sae_mortality" in csv_files:
        df = csv_files["sae_mortality"]
        if {"Age", "Rate"}.issubset(df.columns):
            params["sae_mortality"] = df.groupby("Age")["Rate"].mean().to_dict()

    # --- 6. Community risk multipliers ---
    if "community_risk" in csv_files:
        df = csv_files["community_risk"]
        if {"factor", "risk_multiplier"}.issubset(df.columns):
            params["community_risk"] = df.set_index("factor")["risk_multiplier"].to_dict()

    # --- 7. Comorbidities ---
    if "comorbidities" in csv_files:
        df = csv_files["comorbidities"]
        if {"Condition", "Relative_Risk"}.issubset(df.columns):
            params["comorbidities"] = df.set_index("Condition")["Relative_Risk"].to_dict()

    # --- 8. Utilities ---
    if "utilities" in csv_files:
        df = csv_files["utilities"]
        if {"health_state", "utility"}.issubset(df.columns):
            params["utilities"] = df.set_index("health_state")["utility"].to_dict()

    print("✅ Extracted parameter groups:", list(params.keys()))
    return params



def summarize_core_parameters(params):
    """
    Summarize the extracted key parameters.
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
