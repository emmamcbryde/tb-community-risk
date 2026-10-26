# Explicit-TBI application

From the repository root, run:

```powershell
python -m streamlit run explicit_tbi_ui\app.py
```

The application provides **Set up**, **Run analysis**, **Results**, and
**Evidence and technical information**. The main Set up page contains the known
recent and remote-only TBI prevalence inputs, as percentages of total population.

Milestone 3G adds **Progression risk among people with TBI** on Set up. The
working default is the provisional **Risk-enriched population — 8×** scenario.
Standard uses 1×; Custom requires a positive finite multiplier. This multiplies
future progression hazard, not prevalence or final probability. The 8× value
is a working scenario assumption, not calibrated to 20 cases or independently
estimated. Existing risk-factor OR multiplication remains inactive.

Results include baseline-group contributions and an analytical 1× versus selected
risk comparison using otherwise identical assumptions. At the working defaults,
the 20-year comparator is about 2.22 cases at 1× and 16.38 at 8×. Neither new
infections during follow-up nor endogenous transmission are modelled. See the
current [milestone status](../docs/explicit_recent_remote_tbi_v2_status.md).

The former root-level `explicit_tbi_app.py` entry point has been removed. Its
location next to the legacy `pages/` directory caused automatic legacy multipage
discovery. This directory deliberately has no `pages/` subdirectory. Do not add
legacy page links here or run an entry point alongside the repository's legacy
pages. No repository-wide Streamlit configuration is changed.

Python 3.12.4 / Streamlit 1.32.0 was the observed local environment at the start
of Milestone 3F. The repository pins Streamlit 1.39.0; both versions are checked
for the isolated app. The frozen `streamlit_app.py` uses `st.navigation`, which
requires a supported newer Streamlit version, and retains its original navigation.

All session and widget keys use the `explicit_tbi_ui_v1:` namespace. Unmarked,
foreign, retired or incompatible state is rejected on startup, including stale
economics. The legacy economics page also stops before loading anything when
explicit-app keys are present. Explicit results retain an unavailable,
deliberately non-costable event-ledger marker. No economic outputs are produced.

Presentation helpers live in `app/explicit_tbi_formatting.py`. Percentages use
one decimal, with a meaningful second decimal permitted for baseline active TB.
Ordinary results use at most two decimal places and thousands separators. Whole
people are integers. Fractional years remain deliberately supported. Percentage
widgets commit back to model proportions only on a user edit. JSON and CSV
provenance downloads keep full precision. Hazards, residuals, hashes and timing
diagnostics belong under Evidence and technical information.

The unused future currency contract formats whole Australian dollars with
thousands separators and a consistent leading minus, for example `AUD 95,052`,
`-AUD 92,370`, and `AUD 0`. Half-dollar ties round away from zero. A future,
authorised arithmetic ICER would use the same whole-dollar rule per DALY. These
formatters do not calculate costs or create economic results.

Run rendered isolation and formatting checks:

```powershell
python -m pytest tests/test_explicit_tbi_ui_isolation.py tests/test_explicit_tbi_streamlit_app.py
python scripts/validate_explicit_tbi_ui.py
```

The live protocol audit starts its own localhost server, verifies the page
registry and rendered UI messages, measures cold/warm clicks and dispatch counts,
and rejects a direct Economics route. Generated evidence goes under
`abm/output/milestone3f/`. `--streamlit-target` optionally selects an isolated
Streamlit validation install; `--frozen` checks the original application's
rendered navigation; `--apptest` runs regression tests under that selected version.
