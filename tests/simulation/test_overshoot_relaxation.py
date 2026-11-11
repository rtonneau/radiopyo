#!/usr/bin/env python3
"""
Tests of the overshoot relaxation of the species below 0 (ADR 0007).

A dynamic species j with c_j < 0 gets, on its own row, the linear
continuation from c = 0+: dc_j/dt += J_jj(0+) * c_j, with J_jj(0+) the
diagonal Jacobian entry computed with c_j read at 0+. Every other row still
reads c_j as 0 (no reaction runs backwards), and a state without species
below 0 gives exactly the same derivative as before. The Jacobian of the
term is checked against finite differences in test_jacobian.py.

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
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.ode_functions_solver import (
    DEPLETION_LIMIT,
    derive,
    jacobian,
    overshoot_slopes,
    reaction_rates,
    species_rates,
)

if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.utils.sim_types import SimMatrices


def toy_env(config: dict[str, tp.Any]) -> SimEnv:
    return SimEnv.from_dict(config, log=False)  # type: ignore[arg-type]


@pytest.fixture
def decay_env() -> SimEnv:
    """A -> B, fast first-order decay (k = 1e9 1/s); B is not consumed."""
    return toy_env(
        {
            "parameters": {"pH": 7.0},
            "concentrations": {"initial": {"A": 1e-6, "B": 0.0}},
            "reactions": {
                "k_reaction": [{"reaction": "A -> B", "k_value": 1e9}],
                "radiolytic": {"A": 1.0},
            },
        }
    )


@pytest.fixture
def mixed_env() -> SimEnv:
    """Several kinds of consumption: 2nd order A + B -> C - D (negative
    product D), 2 A -> E, C -> B + A, and an enzymatic S -> P."""
    return toy_env(
        {
            "parameters": {"pH": 7.0},
            "concentrations": {
                "initial": {"A": 1e-3, "B": 2e-3, "C": 1e-4, "D": 1e-6, "S": 1e-3, "enz": 1e-6}
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


def _state(mat: SimMatrices, seed: int = 3) -> NDArray[np.float64]:
    n_dyn = mat.n_species - mat.state.n_cst
    return np.random.default_rng(seed).lognormal(0.0, 1.0, n_dyn)


def _neg(env: SimEnv, mat: SimMatrices, label: str) -> NDArray[np.float64]:
    y = np.zeros(mat.n_species - mat.state.n_cst)
    y[env.species[label].index] = -1.0
    return y


# --------------------------------------------------------------------------- #
#                          Single stiff first-order decay                     #
# --------------------------------------------------------------------------- #
def test_stiff_decay_below_zero_is_restored(decay_env: SimEnv) -> None:
    env = decay_env
    mat = env.build_all_matrices()
    a, b = env.species["A"].index, env.species["B"].index
    k = 1e9

    y = np.zeros(mat.n_species - mat.state.n_cst)
    y[a] = -1e-3  # µmol/L
    dydt = derive(0.0, y, None, mat)
    assert dydt[a] == pytest.approx(-k * y[a], rel=1e-12)  # > 0: back to 0
    assert dydt[a] > 0.0
    assert dydt[b] == 0.0  # B is not produced from a species below 0

    # Integrated, the overshoot decays to 0 like exp(-k t)
    sol = solve_ivp(
        lambda t, yy: derive(t, yy, None, mat),
        (0.0, 1e-7),
        y,
        method="LSODA",
        jac=lambda t, yy: jacobian(t, yy, None, mat),
        rtol=1e-8,
        atol=1e-14,
    )
    assert sol.success
    assert abs(sol.y[a, -1]) < 1e-6 * abs(y[a])  # exp(-100) of the overshoot
    assert sol.y[a].max() <= 1e-9  # never crosses above 0
    assert sol.y[b].max() <= 1e-9  # nothing produced


def test_slope_value_and_sign(decay_env: SimEnv) -> None:
    env = decay_env
    mat = env.build_all_matrices()
    a = env.species["A"].index
    y_full = mat.state.to_full_molar(_neg(env, mat, "A"))
    slopes = overshoot_slopes(y_full, np.array([a]), mat)
    np.testing.assert_allclose(slopes, [-1e9], rtol=1e-12)


# --------------------------------------------------------------------------- #
#                       Species with no consumption stays put                 #
# --------------------------------------------------------------------------- #
def test_species_without_consumption_stays_put(decay_env: SimEnv) -> None:
    env = decay_env
    mat = env.build_all_matrices()
    b = env.species["B"].index

    y = np.zeros(mat.n_species - mat.state.n_cst)
    y[b] = -1.0
    np.testing.assert_array_equal(derive(0.0, y, None, mat), 0.0)

    y_full = mat.state.to_full_molar(y)
    np.testing.assert_array_equal(overshoot_slopes(y_full, np.array([b]), mat), 0.0)
    # same with a species produced by the others: B follows S @ rates only
    y[env.species["A"].index] = 1.0
    produced = derive(0.0, y, None, mat)
    assert produced[b] == pytest.approx(1e9 * 1.0 * 1.0, rel=1e-12)  # k A, µmol/L/s


def test_2a_repeated_reactant_has_no_slope() -> None:
    """2 A -> E: d(rate)/dA = 2 k A = 0 at 0+, so A < 0 stays put."""
    env = toy_env(
        {
            "parameters": {"pH": 7.0},
            "concentrations": {"initial": {"A": 1e-6, "E": 0.0}},
            "reactions": {"k_reaction": [{"reaction": "2 A -> E", "k_value": 1e9}]},
        }
    )
    mat = env.build_all_matrices()
    y = np.zeros(mat.n_species - mat.state.n_cst)
    y[env.species["A"].index] = -1.0
    np.testing.assert_array_equal(derive(0.0, y, None, mat), 0.0)


# --------------------------------------------------------------------------- #
#                      Other rows and non-negative states unchanged           #
# --------------------------------------------------------------------------- #
def test_state_without_negative_species_is_s_times_rates(mixed_env: SimEnv) -> None:
    mat = mixed_env.build_all_matrices()
    y = _state(mat)
    y_full = mat.state.to_full_molar(y)
    expected = species_rates(reaction_rates(0.0, y_full, mat, None), mat)
    np.testing.assert_array_equal(
        derive(0.0, y, None, mat), mat.state.to_dynamic_micromolar(expected)
    )


@pytest.mark.parametrize("labels", [("A",), ("A", "S"), ("D", "enz"), ("B", "C", "D")])
def test_other_rows_are_unchanged(mixed_env: SimEnv, labels: tuple[str, ...]) -> None:
    """Only the rows of the negative species differ from the derivative of
    the clipped state: mass balance of the other species is unaffected."""
    env = mixed_env
    mat = env.build_all_matrices()
    y = _state(mat)
    negative = [env.species[label].index for label in labels]
    y[negative] *= -1.0

    dydt = derive(0.0, y, None, mat)
    clipped = derive(0.0, np.maximum(y, 0.0), None, mat)
    others = np.setdiff1d(np.arange(y.size), negative)
    np.testing.assert_array_equal(dydt[others], clipped[others])

    # the negative rows: clipped derivative + J_jj(0+) c_j
    y_full = mat.state.to_full_molar(y)
    slopes = overshoot_slopes(y_full, np.array(negative), mat)
    np.testing.assert_allclose(
        dydt[negative], clipped[negative] + slopes * y[negative], rtol=1e-12, atol=0.0
    )


def test_negative_product_species_relaxes_against_the_switch(mixed_env: SimEnv) -> None:
    """D < 0 is a negative product of A + B -> C - D: its row is
    -k A B s(D) (0 at D <= 0) and the slope is -k A B / eps."""
    env = mixed_env
    mat = env.build_all_matrices()
    a, b, d = (env.species[n].index for n in ("A", "B", "D"))
    y = _state(mat)
    y[d] = -1e-3
    y_full = mat.state.to_full_molar(y)
    slope = overshoot_slopes(y_full, np.array([d]), mat)[0]
    assert slope == pytest.approx(-1e9 * y_full[a] * y_full[b] / DEPLETION_LIMIT, rel=1e-12)
    assert derive(0.0, y, None, mat)[d] == pytest.approx(slope * y_full[d] * 1e6, rel=1e-12)
    assert derive(0.0, y, None, mat)[d] > 0.0  # restoring


def test_enzymatic_substrate_slope(mixed_env: SimEnv) -> None:
    env = mixed_env
    mat = env.build_all_matrices()
    s, enz = env.species["S"].index, env.species["enz"].index
    y = _state(mat)
    y[s] = -1e-3
    y_full = mat.state.to_full_molar(y)
    k, km = mat.enzymatic.k_values[0], mat.enzymatic.k_michaelis[0]
    slope = overshoot_slopes(y_full, np.array([s]), mat)[0]
    assert slope == pytest.approx(-k * y_full[enz] / km, rel=1e-12)
