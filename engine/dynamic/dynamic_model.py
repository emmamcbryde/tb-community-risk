import numpy as np

def simulate_dynamic(params, years):
    # TODO: place dynamic ODE logic here
    t = np.arange(years)
    incidence = np.zeros(years)
    
    return {
        "t": t,
        "incidence": incidence,
    }
