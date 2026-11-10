#!/usr/bin/env python3
"""
Tests of the clipped concentrations of the rate laws (ADR 0004): every
rate law reads max(c, 0), so a species driven below 0 by a solver
overshoot is neither consumed further nor produced backwards, and the
ODE right hand side stays exactly stoichiometry @ rates (mass balance).

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np
import pytest
from scipy.integrate import solve_ivp

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics.beam import ConstantBeam
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.ode_functions_solver import (
    DEPLETION_LIMIT,
    derive,
    jacobian,
    reaction_rates,
)

if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.utils.sim_types import SimMatrices

# ------------------------------- Constants --------------------------------- #
# Depleted concentration [µmol/L], far below the depletion threshold
DEPLETED = 1e-3 * DEPLETION_LIMIT * 1e6

STATES = ["positive", "negative", "depleted", "mixed"]


def toy_env(config: dict[str, tp.Any]) -> SimEnv:
    return SimEnv.from_dict(config, log=False)  # type: ignore[arg-type]


@pytest.fixture
def dimerization_env() -> SimEnv:
    """2 A -> B (with A produced by irradiation)."""
    return toy_env(
        {
            "parameters": {"pH": 7.0},
            "concentrations": {"initial": {"A": 1e-6, "B": 0.0}},
            "reactions": {
                "k_reaction": [{"reaction": "2 A -> B", "k_value": 1e9}],
                "radiolytic": {"A": 1.0},
            },
        }
    )


@pytest.fixture
def balance_env() -> SimEnv:
    """A + B -> C - D (negative product D), C -> B + A, 2 A -> E and an
    enzymatic S -> P. Conserved: A + C + 2 E, B + C, S + P."""
    return toy_env(
        {
            "parameters": {"pH": 7.0},
            "concentrations": {
                "initial": {
                    "A": 1e-3,
                    "B": 2e-3,
                    "C": 1e-4,
                    "D": 1e-6,
                    "E": 0.0,
                    "S": 1e-3,
                    "P": 0.0,
                    "enz": 1e-6,
                }
            },
            "reactions": {
                "k_reaction": [
                    {"reaction": "A + B -> C - D", "k_value": 1e9},
                    {"reaction": "C -> B + A", "k_value": 3e2},
                    {"reaction": "2 A -> E", "k_value": 2e8},
                ],
                "enzymatic": [
                    {"reaction": "S -- enz >> P", "k_value": 1e3, "k_micha": 1e-3},
                ],
                "radiolytic": {"A": 1.0, "S": 0.5},
            },
        }
    )


# --------------------------------------------------------------------------- #
#                                  Overshoot                                  #
# --------------------------------------------------------------------------- #
def test_self_reacting_species_below_zero_is_not_consumed(
    dimerization_env: SimEnv,
) -> None:
    env = dimerization_env
    mat = env.build_all_matrices()
    a, b = env.species["A"].index, env.species["B"].index

    y = np.zeros(mat.n_species - mat.state.n_cst)
    y[a] = -1e-3  # µmol/L
    dydt = derive(0.0, y, None, mat)
    assert dydt[a] == 0.0  # was -2 k A^2 < 0: further consumption
    assert dydt[b] == 0.0

    # at A > 0 it is consumed as usual
    y[a] = 1e-3
    assert derive(0.0, y, None, mat)[a] < 0.0


def test_run_started_slightly_negative_stays_bounded(dimerization_env: SimEnv) -> None:
    env = dimerization_env
    mat = env.build_all_matrices()
    a, b = env.species["A"].index, env.species["B"].index
    beam = ConstantBeam(dose_rate=1e3, max_dose=1.0)

    y0 = np.zeros(mat.n_species - mat.state.n_cst)
    y0[a] = -1e-3  # µmol/L

    sol = solve_ivp(
        lambda t, y: derive(t, y, beam, mat),
        (0.0, 1e-3),
        y0,
        method="LSODA",
        jac=lambda t, y: jacobian(t, y, beam, mat),
        rtol=1e-8,
        atol=1e-12,
    )
    assert sol.success
    assert np.isfinite(sol.y).all()
    # A never goes below its starting overshoot, B is only ever produced
    assert sol.y[a].min() >= y0[a] - 1e-9
    assert sol.y[b].min() >= -1e-9
    # irradiation brings A back above 0
    assert sol.y[a, -1] > 0.0


def test_reverse_reaction_not_run_backwards(balance_env: SimEnv) -> None:
    """C -> B + A with C < 0 must not produce A and B (nor consume C)."""
    env = balance_env
    mat = env.build_all_matrices()
    y = np.zeros(mat.n_species - mat.state.n_cst)
    y[env.species["C"].index] = -1.0
    y[env.species["D"].index] = 1.0
    np.testing.assert_array_equal(derive(0.0, y, None, mat), 0.0)


# --------------------------------------------------------------------------- #
#                                Mass balance                                 #
# --------------------------------------------------------------------------- #
def _states(env: SimEnv, mat: SimMatrices) -> dict[str, NDArray[np.float64]]:
    """Solver states [µmol/L]: positive, with negative species, with
    near-depleted species, and both."""
    n_dyn = mat.n_species - mat.state.n_cst
    idx = {sp.label: sp.index for sp in env.species.iter_dyn_species()}
    positive = np.random.default_rng(7).lognormal(0.0, 1.0, n_dyn)

    negative = positive.copy()
    for label in ("A", "S"):
        negative[idx[label]] *= -1.0

    depleted = positive.copy()
    for label in ("B", "D", "P", "enz"):
        depleted[idx[label]] = DEPLETED

    mixed = negative.copy()
    mixed[idx["D"]] = DEPLETED
    mixed[idx["enz"]] = -1e-3
    return {
        "positive": positive,
        "negative": negative,
        "depleted": depleted,
        "mixed": mixed,
    }


@pytest.mark.parametrize("state", STATES)
def test_derive_is_stoichiometry_times_full_rate_vector(
    balance_env: SimEnv, state: str
) -> None:
    env = balance_env
    mat = env.build_all_matrices()
    beam = ConstantBeam(dose_rate=1e3, max_dose=10.0)
    y = _states(env, mat)[state]

    rates = reaction_rates(0.0, mat.state.to_full_molar(y), mat, beam)
    # same column order as ReactionRates.stacked(): k, radiolytic, enzymatic
    stoichiometry = np.hstack(
        [
            mat.k_reactions.stoichiometry,
            mat.radiolytic.stoichiometry,
            mat.enzymatic.stoichiometry,
        ]
    )
    expected = (stoichiometry @ rates.stacked())[: y.size] * 1e6
    np.testing.assert_allclose(
        derive(0.0, y, beam, mat),
        expected,
        rtol=1e-12,
        atol=1e-12 * np.abs(expected).max(),
    )


@pytest.mark.parametrize("state", STATES)
def test_conserved_combinations_have_zero_derivative(
    balance_env: SimEnv, state: str
) -> None:
    env = balance_env
    mat = env.build_all_matrices()
    y = _states(env, mat)[state]
    idx = {sp.label: sp.index for sp in env.species.iter_dyn_species()}
    # no beam: the radiolytic production does not conserve anything
    dydt = derive(0.0, y, None, mat)
    assert np.abs(dydt).max() > 0.0  # the state does evolve (not vacuous)

    for weights in (
        {"A": 1.0, "C": 1.0, "E": 2.0},
        {"B": 1.0, "C": 1.0},
        {"S": 1.0, "P": 1.0},
    ):
        terms = [w * dydt[idx[label]] for label, w in weights.items()]
        assert abs(sum(terms)) <= 1e-12 * sum(abs(t) for t in terms), weights
