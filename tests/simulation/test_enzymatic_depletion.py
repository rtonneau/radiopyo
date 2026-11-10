#!/usr/bin/env python3
"""
Tests of the negative-product block of enzymatic reactions: an enzymatic
reaction with a negative product does not occur while that product is
depleted (same rule as the k-reactions).

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
    negative_product_block,
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
    ("d", "blocked"),
    [(0.0, True), (DEPLETION_LIMIT, True), (2 * DEPLETION_LIMIT, False), (1e-6, False)],
)
def test_enzymatic_block_on_off(toy_env: SimEnv, d: float, blocked: bool) -> None:
    mat = toy_env.build_all_matrices()
    y_full = np.full(mat.n_species, 1e-4)
    y_full[toy_env.species["D"].index] = d
    s, e = 1e-4, 1e-4

    assert bool(negative_product_block(y_full, mat.enzymatic)[0]) == blocked
    rate = reaction_rates(0.0, y_full, mat, None).enzymatic[0]
    expected = 0.0 if blocked else K_VALUE * e * s / (K_MICHA + s)
    assert rate == pytest.approx(expected, rel=1e-14)


def test_enzymatic_block_vectorised(toy_env: SimEnv) -> None:
    mat = toy_env.build_all_matrices()
    states = np.full((mat.n_species, 2), 1e-4)
    states[toy_env.species["D"].index] = [0.0, 1e-6]
    rates = reaction_rates(np.zeros(2), states, mat, None).enzymatic
    assert rates[0, 0] == 0.0
    assert rates[0, 1] > 0.0


def test_bundled_catalase_has_no_negative_product() -> None:
    cell = UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    ez = cell.mat.enzymatic
    assert ez.n_reactions > 0
    assert not ez.negative_products.any()
    assert {f.name for f in dataclasses.fields(ez)} >= {"negative_products"}
