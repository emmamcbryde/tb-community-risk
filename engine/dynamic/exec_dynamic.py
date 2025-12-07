from .dynamic_model import simulate_dynamic

def run_dynamic_model(params, years, intervention=True):
    """
    Thin wrapper around simulate_dynamic so the UI can pass
    intervention=True/False to distinguish baseline vs intervention runs.
    """
    return simulate_dynamic(params, years, intervention=intervention)

