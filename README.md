TB COMMUNITY RISK MODEL
Static and Dynamic LTBI → TB Models
=====================================================

OVERVIEW
-----------------------------------------------------
This repository contains two related tuberculosis (TB) models
with a shared Streamlit interface:

1. Dynamic model
   - Continuous-time compartmental model
   - Sub-annual integration (dt = 0.1 years)
   - Annualised outputs
   - Includes transmission feedback

2. Static model
   - Multi-year algebraic approximation
   - No ODE solver
   - Mirrors dynamic model inputs and outputs
   - Useful for rapid scenario exploration

Both models support:
- Country-based or custom age distributions
- Historical TB incidence patterns
- LTBI back-calculation by age
- Risk factor prevalence inputs
- LTBI test-and-treat interventions (pulse style)
- Baseline vs intervention comparison
- Outputs as:
    • Incidence per 100,000
    • Incidence counts
    • Cases averted


RUNNING LOCALLY
-----------------------------------------------------

1) Create a virtual environment

   python -m venv .venv

   Windows:
   .venv\Scripts\activate

   macOS/Linux:
   source .venv/bin/activate

2) Install required packages

   pip install -r requirements.txt

3) Run the application

   streamlit run ui/app.py

The app will open in your browser at:
   http://localhost:8501

Python must be installed on the local machine to run locally.


RUNNING ON STREAMLIT CLOUD
-----------------------------------------------------
To deploy without requiring Python on users' machines:

1) Push the repository to GitHub
2) Create a new app on Streamlit Cloud
3) Set the entry point to:
      ui/app.py
4) Ensure:
      - requirements.txt is in the repo root
      - data files (e.g., data/population_age_latest.csv) are committed

Users accessing the hosted app only need a web browser.


HISTORICAL INCIDENCE INPUT
-----------------------------------------------------
The LTBI age distribution depends on historical TB incidence.

Available options:
- Constant incidence
- Falling 3% per year
- Rising 3% per year
- Upload CSV (year, incidence)

CSV format:
   year,incidence
   2010,40
   2011,38
   2012,36

If historical data do not cover all required years,
the model extrapolates forward/backward using geometric trends.

For numerical stability:
- Incidence is floored at a small positive value
- Extremely low or zero incidence values are adjusted
  to prevent pathological LTBI estimates


AGE DISTRIBUTION INPUT
-----------------------------------------------------
Options:
- Country ISO3 code (uses OWID population file)
- Upload custom CSV
- Default global distribution

Custom CSV format:
   AgeGroup,Proportion
   0,0.012
   1,0.012
   ...

Proportions must sum to 1.


RISK FACTORS
-----------------------------------------------------
Risk factors modify progression from LTBI → active TB.
They are combined multiplicatively.

Included risk factors:
- Smoking
- Excess alcohol use
- Diabetes
- Renal impairment
- HIV (treated)
- HIV (untreated)

Relative risk defaults:
- Smoking: 1.5
- Alcohol: 2.0
- Diabetes: 3.0
- Renal disease: 2.5
- HIV treated: 4.0
- HIV untreated: 10.0

Risk multiplier calculation:
For each factor:
   Effective RR = (1 - p) * 1 + p * RR

Total risk multiplier:
   Product of all effective RRs


LTBI TEST-AND-TREAT INTERVENTION
-----------------------------------------------------
Treatment removes individuals from L_fast and L_slow
and returns them to the susceptible compartment.

Coverage input represents:
   The fraction of the total population that will
   successfully complete LTBI treatment.

In the dynamic model:
   Coverage is distributed across rollout years.

In the static model:
   Coverage is applied per year in a pulse-style approximation.


OUTPUTS
-----------------------------------------------------
Both models produce:

- Annual incidence per 100,000
- Annual incidence (counts)
- Cases averted (baseline – intervention)
- LTBI prevalence by age (stacked recent vs remote)


NUMERICAL STABILITY
-----------------------------------------------------
To prevent instability in LTBI back-calculation:

- Incidence is floored at a small positive value
- ARI is prevented from reaching exactly zero
- Sub-annual time step (dt = 0.1 years) is used
  for smooth integration


STATIC VS DYNAMIC MODEL
-----------------------------------------------------
Dynamic model:
- Includes transmission feedback
- Continuous-time approximation
- More realistic epidemic behaviour

Static model:
- No transmission feedback
- Multi-year algebraic stepping
- Faster, approximate projection
- Useful for quick comparisons


TROUBLESHOOTING
-----------------------------------------------------
If the app does not open locally:

- Try a different browser (Edge vs Chrome)
- Confirm http://localhost:8501 is accessible
- Check antivirus or endpoint protection
- On managed work devices, local web apps may be restricted

If using Streamlit Cloud:
- Confirm your organisation allows access to Streamlit domains


CONTACT / SUPPORT
-----------------------------------------------------
If issues persist, provide:
- Operating system
- Browser and version
- Error message (if any)
- Whether issue occurs in multiple browsers

