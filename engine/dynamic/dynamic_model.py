import numpy as np
import pandas as pd

def simulate_dynamic(params, years, initial_conditions=None):
    """
    Simulate dynamic TB model using ODEs, with Fast and Slow Latent TB compartments.

    Parameters:
    - params: Dictionary containing model parameters (beta, sigma_fast, sigma_slow, gamma, etc.)
    - years: Number of years to run the simulation
    - initial_conditions: Initial values for S, L_fast, L_slow, I, and R

    Returns:
    - Dictionary with time series for incidence and compartment values
    """

    # Define model parameters from the CSV files
    beta = params['community_risk']['beta']  # Transmission rate
    sigma_fast = params['reactivation']['sigma_fast']  # Fast progression from latent to active
    sigma_slow = params['reactivation']['sigma_slow']  # Slow progression from latent to active
    gamma = params['interventions']['gamma']  # Recovery/treatment rate

    # Set initial conditions (default if not provided)
    if initial_conditions is None:
        initial_conditions = {
            'S': 0.99,  # 99% susceptible
            'L_fast': 0.01,  # 1% in fast latent TB
            'L_slow': 0.0,   # No initial slow latent TB
            'I': 0.0,   # No initial infected individuals
            'R': 0.0    # No initial recovered individuals
        }

    # Time steps
    t = np.arange(0, years, 1)  # simulation time in years

    # Initialize arrays to store results
    S = np.zeros(len(t))
    L_fast = np.zeros(len(t))
    L_slow = np.zeros(len(t))
    I = np.zeros(len(t))
    R = np.zeros(len(t))

    # Initial values
    S[0], L_fast[0], L_slow[0], I[0], R[0] = initial_conditions['S'], initial_conditions['L_fast'], initial_conditions['L_slow'], initial_conditions['I'], initial_conditions['R']

    # Simulate the system using Euler's method (simple numerical integration)
    for i in range(1, len(t)):
        # Calculate the force of infection
        lambda_t = beta * I[i-1] / (S[i-1] + L_fast[i-1] + L_slow[i-1] + I[i-1] + R[i-1])  # Force of infection

        # ODEs with fast and slow latent compartments
        dS = -lambda_t * S[i-1]
        dL_fast = lambda_t * S[i-1] - sigma_fast * L_fast[i-1] - (1/5) * L_fast[i-1]  # Fast latent progression
        dL_slow = (1/5) * L_fast[i-1] - sigma_slow * L_slow[i-1]  # Slow latent progression
        dI = sigma_fast * L_fast[i-1] + sigma_slow * L_slow[i-1] - gamma * I[i-1]
        dR = gamma * I[i-1]

        # Update compartments
        S[i] = S[i-1] + dS
        L_fast[i] = L_fast[i-1] + dL_fast
        L_slow[i] = L_slow[i-1] + dL_slow
        I[i] = I[i-1] + dI
        R[i] = R[i-1] + dR

    # Output results
    results = {
        "time": t,
        "S": S,
        "L_fast": L_fast,
        "L_slow": L_slow,
        "I": I,
        "R": R,
        "incidence": I  # Active cases
    }

    return results
