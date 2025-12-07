from .dynamic_model import simulate_dynamic

def run_dynamic_model(params, years):
    results = simulate_dynamic(params, years)
    return results
