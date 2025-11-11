import pandas as pd
from pathlib import Path

def load_parameters(file_path="data/parameters.xlsx"):
    """
    Load all parameter tables from the Excel workbook.
    Returns a dictionary of DataFrames, one per sheet.
    """

    file = Path(file_path)
    if not file.exists():
        raise FileNotFoundError(f"Parameter file not found: {file.resolve()}")

    xls = pd.ExcelFile(file)
    params = {}
    for sheet in xls.sheet_names:
        try:
            params[sheet.lower()] = xls.parse(sheet)
        except Exception as e:
            print(f"⚠️ Could not read sheet '{sheet}': {e}")
    return params


def get_default_values():
    """
    Provide a fallback minimal parameter set for testing.
    """
    return {
        "reactivation_rate": 0.001,
        "mortality_rate": 0.0007,
        "tb_mortality_rate": 0.05,
        "treatment_efficacy": 0.8,
        "testing_sensitivity": 0.85,
        "testing_specificity": 0.95,
    }
