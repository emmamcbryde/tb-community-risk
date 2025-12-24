import pandas as pd
from pathlib import Path

def load_static_parameters(directory="data/params_static"):
    params = {}
    for f in Path(directory).glob("*.csv"):
        params[f.stem] = pd.read_csv(f)
    return params
