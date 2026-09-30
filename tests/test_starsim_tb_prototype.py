"""Tests for the Starsim TB natural-history feasibility prototype (Milestone 1).

Software-verification tests only: passing them shows the code does what the
specification says, not that the model or its demonstration values are
scientifically valid.

Monte Carlo comparisons use a prespecified tolerance of 4 binomial standard
errors (two-sided false-failure probability about 6e-5 per comparison), plus the
known discretisation error where the discrete scheme is not exact. Seeds are
fixed, so each test is also deterministic.
"""

from __future__ import annotations

import importlib.util
import itertools
import json
import math
import unittest

import numpy as np

from engine.starsim_tb import analytics as an
from engine.starsim_tb.parameters import (
    BIOLOGICAL_STATES,
    PrototypeParameters,
    demonstration_parameters,
    from_recent_mean_residence,
    validate_initial_state_proportions,
)

HAS_STARSIM = importlib.util.find_spec("starsim") is not None
Z = 4.0


def cohort(state: str, n: int = 20_000, years: float = 20.0, **changes) -> PrototypeParameters:
    """Closed cohort starting entirely in ``state`` with no transmission."""
    props = {k: 0.0 for k in BIOLOGICAL_STATES}
    props[state] = 1.0
    return demonstration_parameters(
        n_agents=n,
        initial_state_proportions=props,
        start_year=0.0,
        stop_year=years,
        network_mean_degree=0.0,
        transmission_hazard_per_contact_per_year=0.0,
        **changes,
    )


# Pure helpers (no Starsim required) ------------------------------------------------------

class HazardConversionTests(unittest.TestCase):
    def test_exact_conversion(self) -> None:
        self.assertAlmostEqual(an.hazard_to_probability(0.01, 1 / 12), 1 - math.exp(-0.01 / 12), places=15)
        self.assertAlmostEqual(an.hazard_to_probability(2.0, 1.0), 1 - math.exp(-2.0), places=15)
        self.assertEqual(an.hazard_to_probability(0.0, 1 / 12), 0.0)
        np.testing.assert_allclose(an.hazard_to_probability(np.array([0.1, 1.0]), 0.5), 1 - np.exp(-np.array([0.05, 0.5])))

    def test_round_trip(self) -> None:
        for h in (1e-6, 0.001, 0.19, 5.0):
            self.assertAlmostEqual(an.probability_to_hazard(an.hazard_to_probability(h, 1 / 12), 1 / 12), h, delta=1e-9 * max(h, 1))

    def test_linear_approximation_is_diagnostic_only_and_differs(self) -> None:
        exact = an.hazard_to_probability(2.0, 1.0)
        approx = an.linear_probability_approximation(2.0, 1.0)
        self.assertGreater(approx, 1.0)  # h dt can exceed one: not a probability
        self.assertLess(exact, 1.0)
        # For the demonstration monthly hazards the two agree to second order only.
        h, dt = 0.2, 1 / 12
        self.assertAlmostEqual(an.linear_probability_approximation(h, dt) - an.hazard_to_probability(h, dt), (h * dt) ** 2 / 2, delta=1e-5)

    def test_negative_inputs_rejected(self) -> None:
        with self.assertRaises(ValueError):
            an.hazard_to_probability(-0.1, 1.0)


class CompetingRiskHelperTests(unittest.TestCase):
    def test_step_probabilities_are_order_invariant(self) -> None:
        hazards = {"active": 0.01, "remote": 0.19, "other": 0.3}
        reference = an.competing_step_probabilities(hazards, 1 / 12)
        for perm in itertools.permutations(hazards):
            probs = an.competing_step_probabilities({k: hazards[k] for k in perm}, 1 / 12)
            for key in hazards:
                self.assertEqual(probs[key], reference[key])
        self.assertAlmostEqual(sum(reference.values()), 1 - math.exp(-0.5 / 12), places=15)

    def test_naive_sequential_bernoulli_is_order_dependent(self) -> None:
        # Diagnostic contrast: applying independent Bernoullis in sequence biases towards the first.
        pa, pb = an.hazard_to_probability(0.5, 1.0), an.hazard_to_probability(1.5, 1.0)
        a_first = (pa, (1 - pa) * pb)
        b_first = ((1 - pb) * pa, pb)
        self.assertNotAlmostEqual(a_first[0], b_first[0], places=3)

    @unittest.skipUnless(HAS_STARSIM, "starsim not installed")
    def test_assignment_identical_under_any_hazard_ordering(self) -> None:
        from engine.starsim_tb.prototype import assign_competing_outcomes

        rng = np.random.default_rng(1)
        u_event, u_which = rng.random(200_000), rng.random(200_000)
        hazards = {"active": 0.4, "remote": 3.0, "third": 1.0}
        reference = assign_competing_outcomes(u_event, u_which, hazards, 0.25)
        for perm in itertools.permutations(hazards):
            out = assign_competing_outcomes(u_event, u_which, {k: hazards[k] for k in perm}, 0.25)
            np.testing.assert_array_equal(out, reference)
        expected = an.competing_step_probabilities(hazards, 0.25)
        n = len(u_event)
        for key, p in expected.items():
            self.assertAlmostEqual(np.mean(reference == key), p, delta=an.binomial_tolerance(p, n, Z))
        self.assertAlmostEqual(np.mean(reference == ""), 1 - sum(expected.values()), delta=an.binomial_tolerance(0.5, n, Z))


class AnalyticFormulaTests(unittest.TestCase):
    h_f, gamma, h_r = 0.01, 0.19, 0.001

    def test_closed_form_matches_numerical_integration(self) -> None:
        for t in (0.5, 5.0, 20.0, 60.0):
            closed = an.recent_cumulative_active(self.h_f, self.gamma, self.h_r, t)
            numeric = an.recent_cumulative_active_numerical(self.h_f, self.gamma, self.h_r, t)
            self.assertAlmostEqual(closed, numeric, delta=1e-10)

    def test_reduces_to_recent_only_formula_when_remote_hazard_zero(self) -> None:
        t = 20.0
        self.assertAlmostEqual(
            an.recent_cumulative_active(self.h_f, self.gamma, 0.0, t),
            self.h_f / (self.h_f + self.gamma) * (1 - math.exp(-(self.h_f + self.gamma) * t)),
            places=14,
        )

    def test_probabilities_are_coherent(self) -> None:
        t = 20.0
        recent = an.recent_still_recent(self.h_f, self.gamma, t)
        ever_remote = an.recent_ever_remote(self.h_f, self.gamma, t)
        direct = an.recent_active_while_recent(self.h_f, self.gamma, t)
        self.assertAlmostEqual(recent + ever_remote + direct, 1.0, places=14)
        self.assertLess(an.recent_active_via_remote(self.h_f, self.gamma, self.h_r, t), ever_remote)

    def test_equal_hazard_limit(self) -> None:
        # H == h_r uses the limiting form; compare with a hair either side.
        mid = an.recent_active_via_remote(0.1, 0.2, 0.3, 10.0)
        near = an.recent_active_via_remote(0.1, 0.2, 0.3 + 1e-7, 10.0)
        self.assertAlmostEqual(mid, near, places=6)

    def test_monthly_discrete_chain_is_close_to_continuous_time(self) -> None:
        # Check 10 (deterministic part): the discrete scheme converges as dt shrinks.
        t = 20.0
        target = an.recent_cumulative_active(self.h_f, self.gamma, self.h_r, t)
        errors = {}
        for steps_per_year in (12, 52, 365):
            dist = an.discrete_chain_distribution([1, 0, 0], self.h_f, self.gamma, self.h_r, 1 / steps_per_year, 20 * steps_per_year)
            errors[steps_per_year] = abs(dist[2] - target)
            # Recent-only quantities are exact at grid points for any dt.
            self.assertAlmostEqual(dist[0], an.recent_still_recent(self.h_f, self.gamma, t), places=12)
        self.assertLess(errors[12], 1e-4)
        self.assertLess(errors[52], errors[12])
        self.assertLess(errors[365], errors[52])

    def test_expected_person_time(self) -> None:
        t = 20.0
        self.assertAlmostEqual(an.expected_person_years_recent(self.h_f, self.gamma, t), (1 - math.exp(-4.0)) / 0.2, places=12)
        total = (
            an.expected_person_years_recent(self.h_f, self.gamma, t)
            + an.expected_person_years_remote(self.h_f, self.gamma, self.h_r, t)
        )
        self.assertLess(total, t)


class ParameterContractTests(unittest.TestCase):
    def test_demonstration_values(self) -> None:
        p = demonstration_parameters()
        self.assertEqual(p.n_agents, 10_000)
        self.assertAlmostEqual(p.timestep_years, 1 / 12)
        self.assertEqual(p.recent_progression_hazard_per_year, 0.01)
        self.assertEqual(p.remote_progression_hazard_per_year, 0.001)
        self.assertAlmostEqual(p.recent_mean_residence_years, 5.0, places=12)
        self.assertEqual(p.n_timesteps, 240)
        self.assertAlmostEqual(from_recent_mean_residence(5.0).recent_to_remote_hazard_per_year, 0.19, places=12)

    def test_initial_proportions_validate_and_sum_to_one(self) -> None:
        # Check 15.
        p = demonstration_parameters()
        self.assertAlmostEqual(math.fsum(p.initial_state_proportions.values()), 1.0, places=12)
        validate_initial_state_proportions(p.initial_state_proportions)
        bad = [
            {"susceptible": 0.5, "recent_infection": 0.2, "remote_infection": 0.2, "active_pulmonary_tb": 0.2},
            {"susceptible": 1.0, "recent_infection": 0.0, "remote_infection": 0.0},
            {"susceptible": 1.1, "recent_infection": -0.1, "remote_infection": 0.0, "active_pulmonary_tb": 0.0},
            {"susceptible": 0.9, "recent_infection": 0.1, "remote_infection": 0.0, "active_pulmonary_tb": 0.0, "subclinical": 0.0},
        ]
        for props in bad:
            with self.subTest(props=props), self.assertRaises(ValueError):
                demonstration_parameters(initial_state_proportions=props)

    def test_other_invalid_parameters_rejected(self) -> None:
        for changes in (
            {"recent_progression_hazard_per_year": -0.01},
            {"remote_progression_hazard_per_year": math.nan},
            {"stop_year": 2000.0},
            {"timestep_years": 0.07},  # 20 years is not a whole number of steps
            {"n_agents": 0},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                demonstration_parameters(**changes)

    def test_canonical_serialisation_and_hash(self) -> None:
        a = demonstration_parameters()
        props = dict(reversed(list(a.initial_state_proportions.items())))
        b = demonstration_parameters(initial_state_proportions=props)
        self.assertEqual(a.canonical_json(), b.canonical_json())
        self.assertEqual(a.parameter_hash(), b.parameter_hash())
        self.assertNotEqual(a.parameter_hash(), a.with_updates(rand_seed=a.rand_seed + 1).parameter_hash())
        self.assertEqual(json.loads(a.canonical_json())["parameters"]["recent_progression_hazard_per_year"], 0.01)


# Simulation tests -------------------------------------------------------------------------

@unittest.skipUnless(HAS_STARSIM, "starsim not installed")
class SimulationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from engine.starsim_tb import prototype

        cls.proto = prototype
        cls.demo = prototype.run_prototype(demonstration_parameters())

    def run_pars(self, pars: PrototypeParameters, **kwargs):
        return self.proto.run_prototype(pars, **kwargs)

    # 1, 2 ---------------------------------------------------------------------------------
    def test_state_exclusivity_and_conservation_every_timestep(self) -> None:
        for run in (self.demo, self.run_pars(cohort("recent_infection", n=5_000))):
            res = run.integrity.results
            self.assertTrue(np.all(res.n_state_violations[:] == 0))
            self.assertTrue(np.all(res.population_balance[:] == 0))
            self.assertTrue(np.all(res.n_in_exactly_one_state[:] == run.pars.n_agents))
            d = run.disease
            stocks = sum(np.asarray(d.results[f"n_{k}"][:]) for k in BIOLOGICAL_STATES)
            self.assertTrue(np.all(stocks == run.pars.n_agents))
            membership = sum(getattr(d, k).raw.astype(int) for k in BIOLOGICAL_STATES)
            self.assertTrue(np.all(membership == 1))

    def test_stock_changes_reconcile_with_flows(self) -> None:
        r = self.demo.disease.results
        ds = np.diff(np.asarray(r.n_susceptible[:]))
        drecent = np.diff(np.asarray(r.n_recent_infection[:]))
        dremote = np.diff(np.asarray(r.n_remote_infection[:]))
        dactive = np.diff(np.asarray(r.n_active_pulmonary_tb[:]))
        new_inf, r2r = np.asarray(r.new_infections[1:]), np.asarray(r.new_recent_to_remote[1:])
        a_rec, a_rem = np.asarray(r.new_active_from_recent[1:]), np.asarray(r.new_active_from_remote[1:])
        np.testing.assert_array_equal(ds, -new_inf)
        np.testing.assert_array_equal(drecent, new_inf - r2r - a_rec)
        np.testing.assert_array_equal(dremote, r2r - a_rem)
        np.testing.assert_array_equal(dactive, a_rec + a_rem)
        self.assertEqual(r.new_infections[0], 0)
        self.assertEqual(r.new_active[0], 0)

    # 3 ------------------------------------------------------------------------------------
    def test_zero_transmission_gives_no_new_infections(self) -> None:
        run = self.run_pars(demonstration_parameters(transmission_hazard_per_contact_per_year=0.0))
        r = run.disease.results
        self.assertEqual(int(np.sum(r.new_infections[:])), 0)
        self.assertEqual(len(set(np.asarray(r.n_susceptible[:]).tolist())), 1)

    # 4 ------------------------------------------------------------------------------------
    def test_zero_progression_gives_no_new_active(self) -> None:
        run = self.run_pars(
            demonstration_parameters(recent_progression_hazard_per_year=0.0, remote_progression_hazard_per_year=0.0)
        )
        r = run.disease.results
        self.assertEqual(int(np.sum(r.new_active[:])), 0)
        self.assertEqual(len(set(np.asarray(r.n_active_pulmonary_tb[:]).tolist())), 1)
        self.assertGreater(int(np.sum(r.new_infections[:])), 0)  # the initial active cases still transmit

    # 5 ------------------------------------------------------------------------------------
    def test_zero_ageing_hazard_keeps_recent_recent_unless_progressing(self) -> None:
        run = self.run_pars(cohort("recent_infection", n=10_000, recent_to_remote_hazard_per_year=0.0))
        r = run.disease.results
        self.assertTrue(np.all(np.asarray(r.n_remote_infection[:]) == 0))
        self.assertEqual(int(r.cum_recent_to_remote[-1]), 0)
        np.testing.assert_array_equal(
            np.asarray(r.n_recent_infection[:]) + np.asarray(r.n_active_pulmonary_tb[:]), run.pars.n_agents
        )
        p = an.remote_cumulative_active(0.01, 20.0)  # recent progression is the only exit
        self.assertAlmostEqual(r.n_active_pulmonary_tb[-1] / run.pars.n_agents, p, delta=an.binomial_tolerance(p, run.pars.n_agents, Z))

    # 6 ------------------------------------------------------------------------------------
    def test_very_high_ageing_hazard_moves_survivors_to_remote_quickly(self) -> None:
        n = 20_000
        run = self.run_pars(cohort("recent_infection", n=n, years=2.0, recent_to_remote_hazard_per_year=1_000.0))
        r = run.disease.results
        self.assertEqual(int(r.n_recent_infection[1]), 0)  # exp(-1000/12) is effectively zero
        p_direct = 0.01 / 1000.01
        self.assertLessEqual(int(r.cum_active_from_recent[-1]), n * p_direct + Z * math.sqrt(n * p_direct) + 1)
        self.assertGreater(r.n_remote_infection[1], 0.999 * n)

    # 7 ------------------------------------------------------------------------------------
    def test_remote_progression_lower_than_recent(self) -> None:
        pars = demonstration_parameters()
        self.assertLess(pars.remote_progression_hazard_per_year, pars.recent_progression_hazard_per_year)
        n = 40_000
        props = {"susceptible": 0.0, "recent_infection": 0.5, "remote_infection": 0.5, "active_pulmonary_tb": 0.0}
        run = self.run_pars(cohort("recent_infection", n=n).with_updates(initial_state_proportions=props))
        s = self.proto.summary(run)
        rate_recent = s["rate_per_100k_person_years"]["active_from_recent"]
        rate_remote = s["rate_per_100k_person_years"]["active_from_remote"]
        self.assertGreater(rate_recent, 3 * rate_remote)
        # Empirical incidence per person-year matches each input hazard (Poisson, 4 SE).
        for events, py, hazard in (
            (s["cumflow_total"]["active_from_recent"], s["persontime_years_total"]["recent_infection"], 0.01),
            (s["cumflow_total"]["active_from_remote"], s["persontime_years_total"]["remote_infection"], 0.001),
        ):
            expected = hazard * py
            self.assertLess(abs(events - expected), Z * math.sqrt(expected) + 0.01 * expected)

    # 8 ------------------------------------------------------------------------------------
    def test_remote_cumulative_progression_matches_analytic(self) -> None:
        for h_r in (0.001, 0.05):  # demonstration value, then a larger one for statistical power
            n = 20_000
            run = self.run_pars(cohort("remote_infection", n=n, remote_progression_hazard_per_year=h_r))
            r = run.disease.results
            times = np.arange(len(r.n_active_pulmonary_tb)) * run.pars.timestep_years
            for ti in (12, 60, 240):
                p = an.remote_cumulative_active(h_r, times[ti])
                with self.subTest(h_r=h_r, years=times[ti]):
                    self.assertAlmostEqual(r.cum_active_from_remote[ti] / n, p, delta=an.binomial_tolerance(p, n, Z))

    def test_remote_progression_unbiased_across_independent_seeds(self) -> None:
        # Runs sharing a seed share Starsim's common random numbers, so single-seed
        # comparisons in different cohorts are correlated; pool independent seeds here.
        n, years, h_r, seeds = 20_000, 10.0, 0.05, range(1, 9)
        finals = []
        for seed in seeds:
            run = self.run_pars(cohort("remote_infection", n=n, years=years, remote_progression_hazard_per_year=h_r, rand_seed=seed))
            finals.append(run.disease.results.cum_active_from_remote[-1] / n)
        p = an.remote_cumulative_active(h_r, years)
        self.assertAlmostEqual(float(np.mean(finals)), p, delta=an.binomial_tolerance(p, n * len(finals), Z))
        # Replicate spread is consistent with binomial sampling (not degenerate, not inflated).
        sd = float(np.std(finals, ddof=1))
        self.assertGreater(sd, 0.3 * math.sqrt(p * (1 - p) / n))
        self.assertLess(sd, 2.0 * math.sqrt(p * (1 - p) / n))

    # 9 ------------------------------------------------------------------------------------
    def test_recent_competing_outcomes_match_analytic(self) -> None:
        h_f, gamma = 0.01, 0.19
        n, seeds = 20_000, (11, 12, 13)
        pooled = {"active": 0, "recent": 0, "remote": 0}
        for seed in seeds:
            run = self.run_pars(cohort("recent_infection", n=n, remote_progression_hazard_per_year=0.0, rand_seed=seed))
            r = run.disease.results
            pooled["active"] += int(r.n_active_pulmonary_tb[-1])
            pooled["recent"] += int(r.n_recent_infection[-1])
            pooled["remote"] += int(r.n_remote_infection[-1])
        total = n * len(seeds)
        t = 20.0
        expected = {
            "active": an.recent_active_while_recent(h_f, gamma, t),
            "recent": an.recent_still_recent(h_f, gamma, t),
            "remote": an.recent_ever_remote(h_f, gamma, t),
        }
        for key, p in expected.items():
            with self.subTest(outcome=key):
                self.assertAlmostEqual(pooled[key] / total, p, delta=an.binomial_tolerance(p, total, Z))

    def test_recent_with_later_remote_progression_matches_analytic(self) -> None:
        n = 40_000
        run = self.run_pars(cohort("recent_infection", n=n))
        r = run.disease.results
        p = an.recent_cumulative_active(0.01, 0.19, 0.001, 20.0)
        discretisation = abs(an.discrete_chain_distribution([1, 0, 0], 0.01, 0.19, 0.001, 1 / 12, 240)[2] - p)
        self.assertAlmostEqual(r.n_active_pulmonary_tb[-1] / n, p, delta=an.binomial_tolerance(p, n, Z) + discretisation)
        # Person-time in recent infection against its expectation (trapezoid output).
        expected_py = n * an.expected_person_years_recent(0.01, 0.19, 20.0)
        self.assertLess(abs(np.sum(r.person_years_recent[:]) - expected_py) / expected_py, 0.02)

    # 10 -----------------------------------------------------------------------------------
    def test_monthly_and_finer_timesteps_consistent(self) -> None:
        n, years = 20_000, 10.0
        finals = {}
        for dt in (1 / 12, 1 / 48):
            run = self.run_pars(cohort("recent_infection", n=n, years=years, timestep_years=dt))
            r = run.disease.results
            finals[dt] = {k: int(r[f"n_{k}"][-1]) / n for k in ("recent_infection", "remote_infection", "active_pulmonary_tb")}
        for key in finals[1 / 12]:
            a, b = finals[1 / 12][key], finals[1 / 48][key]
            pooled = (a + b) / 2
            with self.subTest(state=key):
                self.assertLess(abs(a - b), Z * math.sqrt(2 * pooled * (1 - pooled) / n) + 1e-4)

    # 11, 12 -------------------------------------------------------------------------------
    def test_fixed_seed_reproduces_event_history_and_outputs(self) -> None:
        again = self.run_pars(demonstration_parameters())
        self.assertEqual(self.proto.results_table(self.demo), self.proto.results_table(again))
        for arr in ("ti_infected", "ti_remote", "ti_active", "active_origin"):
            np.testing.assert_array_equal(getattr(self.demo.disease, arr).raw, getattr(again.disease, arr).raw)
        self.assertEqual(self.demo.metadata["parameter_hash_sha256"], again.metadata["parameter_hash_sha256"])

    def test_different_seed_changes_stochastic_outcomes(self) -> None:
        other = self.run_pars(demonstration_parameters(rand_seed=demonstration_parameters().rand_seed + 1))
        self.assertNotEqual(self.proto.results_table(self.demo), self.proto.results_table(other))
        self.assertFalse(np.array_equal(self.demo.disease.ti_active.raw, other.disease.ti_active.raw, equal_nan=True))

    # 13, 14 -------------------------------------------------------------------------------
    def test_only_small_minority_of_recent_infections_progress(self) -> None:
        n = 20_000
        run = self.run_pars(cohort("recent_infection", n=n, years=40.0))
        r = run.disease.results
        from_recent = int(r.cum_active_from_recent[-1]) / n
        lifetime_direct = 0.01 / 0.2  # h_f / (h_f + gamma) = 5 %
        self.assertAlmostEqual(from_recent, an.recent_active_while_recent(0.01, 0.19, 40.0), delta=an.binomial_tolerance(lifetime_direct, n, Z))
        self.assertLess(from_recent, 0.10)

    def test_latent_infection_is_not_a_conveyor_to_active_tb(self) -> None:
        n = 20_000
        run = self.run_pars(cohort("recent_infection", n=n))
        r = run.disease.results
        self.assertLess(r.n_active_pulmonary_tb[-1], 0.10 * n)
        self.assertGreater(r.n_remote_infection[-1], 0.80 * n)
        self.assertGreater(r.cum_recent_to_remote[-1], 10 * r.cum_active_from_recent[-1])

    # Order invariance within the simulation ----------------------------------------------
    def test_simulation_is_invariant_to_hazard_listing_order(self) -> None:
        # The module passes its recent exits to assign_competing_outcomes; listing them
        # in the opposite order must give an identical event history.
        proto = self.proto
        original = proto.assign_competing_outcomes

        def reversed_order(u_event, u_which, hazards, dt):
            return original(u_event, u_which, dict(reversed(list(hazards.items()))), dt)

        pars = cohort("recent_infection", n=5_000, years=10.0, recent_to_remote_hazard_per_year=0.5, recent_progression_hazard_per_year=0.3)
        base = self.run_pars(pars)
        proto.assign_competing_outcomes = reversed_order
        try:
            flipped = self.run_pars(pars)
        finally:
            proto.assign_competing_outcomes = original
        self.assertEqual(proto.results_table(base), proto.results_table(flipped))
        np.testing.assert_array_equal(base.disease.ti_active.raw, flipped.disease.ti_active.raw)

    # 15 (simulation side) -----------------------------------------------------------------
    def test_initial_state_counts_match_proportions(self) -> None:
        n = 100_000
        pars = demonstration_parameters(n_agents=n, stop_year=2001.0)
        run = self.run_pars(pars)
        for key, p in pars.initial_state_proportions.items():
            with self.subTest(state=key):
                self.assertAlmostEqual(run.disease.initial_counts[key] / n, p, delta=an.binomial_tolerance(p, n, Z))
                self.assertEqual(run.disease.initial_counts[key], int(run.disease.results[f"n_{key}"][0]))

    # Transmission -------------------------------------------------------------------------
    def test_only_active_tb_transmits(self) -> None:
        # Many recent and remote infections, no active TB and no progression: nobody is infectious.
        props = {"susceptible": 0.5, "recent_infection": 0.25, "remote_infection": 0.25, "active_pulmonary_tb": 0.0}
        run = self.run_pars(
            demonstration_parameters(
                initial_state_proportions=props,
                transmission_hazard_per_contact_per_year=50.0,
                recent_progression_hazard_per_year=0.0,
                remote_progression_hazard_per_year=0.0,
            )
        )
        self.assertEqual(int(np.sum(run.disease.results.new_infections[:])), 0)
        # rel_trans alone would allow transmission; the infectious alias is what restricts it.
        self.assertTrue(np.all(run.disease.infectious.raw == run.disease.active_pulmonary_tb.raw))

    def test_transmission_enters_recent_and_reconciles(self) -> None:
        import starsim as ss

        pars = demonstration_parameters(
            transmission_hazard_per_contact_per_year=2.0,
            recent_progression_hazard_per_year=0.0,
            remote_progression_hazard_per_year=0.0,
            recent_to_remote_hazard_per_year=0.0,
            stop_year=2005.0,
        )
        run = self.run_pars(pars, analyzers=[ss.infection_log()])
        d, r = run.disease, run.disease.results
        new = int(np.sum(r.new_infections[:]))
        self.assertGreater(new, 100)
        self.assertEqual(new, int(r.cum_infections[-1]))
        self.assertEqual(new, int(r.n_susceptible[0] - r.n_susceptible[-1]))
        # Everyone infected during the run is in recent infection (no exits enabled).
        infected_in_run = (d.ti_infected.raw > 0) & np.isfinite(d.ti_infected.raw)
        self.assertEqual(int(np.count_nonzero(infected_in_run)), new)
        self.assertTrue(np.all(d.recent_infection.raw[infected_in_run]))
        self.assertEqual(int(r.n_recent_infection[-1] - r.n_recent_infection[0]), new)

        log = run.sim.analyzers.infection_log.logs[self.proto.MODULE_NAME].to_df()
        self.assertEqual(len(log), new)  # no seeds are logged after initialisation
        self.assertEqual(log["target"].nunique(), new)  # nobody infected twice
        sources = log["source"].to_numpy(dtype=np.int64)
        targets = log["target"].to_numpy(dtype=np.int64)
        self.assertTrue(np.all(d.active_pulmonary_tb.raw[sources]))
        self.assertTrue(np.all(d.ti_active.raw[sources] <= d.ti_infected.raw[targets]))
        # Initially infected agents keep their original infection time.
        initially_infected = d.ti_infected.raw == 0
        self.assertEqual(int(np.count_nonzero(initially_infected)), d.initial_counts["recent_infection"] + d.initial_counts["remote_infection"])

    def test_metadata_records_reproducibility_fields(self) -> None:
        meta = self.demo.metadata
        for key in ("python_version", "starsim_version", "platform", "git_commit", "parameter_hash_sha256",
                    "rand_seed", "start_year", "stop_year", "timestep_years", "network", "runtime_seconds",
                    "python_build_platform", "host_machine"):
            self.assertIn(key, meta)
        self.assertEqual(meta["starsim_version"], "3.6.1")
        self.assertEqual(meta["parameter_hash_sha256"], self.demo.pars.parameter_hash())
        json.dumps(meta, allow_nan=False)
        table = self.proto.results_table(self.demo)
        self.assertTrue(all(k.split("_")[0] in {"time", "stock", "flow", "cumflow", "stockprop", "persontime"} for k in table))


if __name__ == "__main__":
    unittest.main()
