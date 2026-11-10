#!/usr/bin/env python3
"""
Tests of the depletion switch of enzymatic reactions: an enzymatic
reaction with a negative product is slowed down by c / (c + DEPLETION_LIMIT)
as that product runs out (same rule as the k-reactions, ADR 0005).

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.ode_functions_solver import (
    DEPLETION_LIMIT,
    depletion_switch,
    reaction_rates,
)
from radiopyo.simulation.unit_cell import UnitCell

K_VALUE = 1e3
K_MICHA = 1e-3


@pytest.fixture
def toy_env() -> SimEnv:
    """S --enz>> P - D: the enzymatic reaction consumes D."""
    config = {
        "parameters": {"pH": 7.0},
        "concentrations": {
            "initial": {"S": 1e-3, "enz": 1e-6, "D": 1e-6},
        },
        "reactions": {
            "enzymatic": [
                {"reaction": "S -- enz >> P - D", "k_value": K_VALUE, "k_micha": K_MICHA},
            ],
            "radiolytic": {"S": 1.0},
        },
    }
    return SimEnv.from_dict(config, log=False)  # type: ignore[arg-type]


def test_negative_products_array(toy_env: SimEnv) -> None:
    mat = toy_env.build_all_matrices()
    ez = mat.enzymatic
    d = toy_env.species["D"].index
    assert ez.negative_products.shape == (mat.n_species, 1)
    assert ez.negative_products.dtype == np.bool_
    assert not ez.negative_products.flags.writeable
    assert ez.negative_products.sum() == 1
    assert ez.negative_products[d, 0]
    assert ez.stoichiometry[d, 0] == -1  # D is consumed
    assert ez.stoichiometry[toy_env.species["P"].index, 0] == 1


@pytest.mark.parametrize(
    ("d", "switch"),
    [
        (0.0, 0.0),
        (-1e-9, 0.0),
        (DEPLETION_LIMIT, 0.5),
        (2 * DEPLETION_LIMIT, 2 / 3),
        (1e-10, 1e-10 / (1e-10 + DEPLETION_LIMIT)),
        (1e-6, 1e-6 / (1e-6 + DEPLETION_LIMIT)),
    ],
)
def test_enzymatic_depletion_switch(toy_env: SimEnv, d: float, switch: float) -> None:
    mat = toy_env.build_all_matrices()
    y_full = np.full(mat.n_species, 1e-4)
    y_full[toy_env.species["D"].index] = d
    s, e = 1e-4, 1e-4

    assert depletion_switch(y_full, mat.enzymatic)[0] == pytest.approx(switch, rel=1e-14)
    rate = reaction_rates(0.0, y_full, mat, None).enzymatic[0]
    expected = switch * K_VALUE * e * s / (K_MICHA + s)
    assert rate == pytest.approx(expected, rel=1e-14)


def test_enzymatic_rate_halved_at_the_limit_unchanged_at_1e6(toy_env: SimEnv) -> None:
    mat = toy_env.build_all_matrices()
    d = toy_env.species["D"].index
    y_full = np.full(mat.n_species, 1e-4)
    y_full[d] = 1.0  # fully supplied: the plain Michaelis-Menten rate
    full = reaction_rates(0.0, y_full, mat, None).enzymatic[0]
    y_full[d] = DEPLETION_LIMIT
    assert reaction_rates(0.0, y_full, mat, None).enzymatic[0] == pytest.approx(
        full / 2, rel=1e-6
    )
    y_full[d] = 1e-6
    assert reaction_rates(0.0, y_full, mat, None).enzymatic[0] == pytest.approx(
        full, rel=1e-5
    )


def test_enzymatic_depletion_switch_vectorised(toy_env: SimEnv) -> None:
    mat = toy_env.build_all_matrices()
    states = np.full((mat.n_species, 3), 1e-4)
    states[toy_env.species["D"].index] = [0.0, DEPLETION_LIMIT, 1e-6]
    rates = reaction_rates(np.zeros(3), states, mat, None).enzymatic
    assert rates[0, 0] == 0.0
    assert rates[0, 1] == pytest.approx(0.5 * rates[0, 2], rel=1e-5)
    assert rates[0, 2] > 0.0
    for col in range(3):
        assert rates[0, col] == reaction_rates(
            0.0, states[:, col], mat, None
        ).enzymatic[0]


def test_bundled_catalase_has_no_negative_product() -> None:
    cell = UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    ez = cell.mat.enzymatic
    assert ez.n_reactions > 0
    assert not ez.negative_products.any()
    assert {f.name for f in dataclasses.fields(ez)} >= {"negative_products"}
