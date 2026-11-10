#!/usr/bin/env python3
"""
Tests of the analytic Jacobian of the ODE right hand side: compared with
the central finite-difference Jacobian of ``derive`` on the bundled
configs (beam on/off, enzymes, fixed species) and on toy
environments (repeated reactant, active negative-product blocks).

The Jacobian is piecewise (ADR 0003): the finite differences are taken
with a relative step far smaller than the distance of the states to the
depletion threshold, so both sides agree on every switch.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.physics.beam import ConstantBeam
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.ode_functions_solver import (
    DEPLETION_LIMIT,
    derive,
    jacobian,
    negative_product_block,
)
from radiopyo.simulation.unit_cell import UnitCell

if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.physics.beam import ParticleBeam
    from radiopyo.utils.sim_types import SimMatrices

# ------------------------------- Fixtures ---------------------------------- #
CONFIGS = {
    "labarbe": rp.LABARBE_CONFIG_FILE,
    "water": rp.PURE_WATER_CONFIG_FILE,
}

# The beam (constant, 1e7 Gy/s up to 10 Gy) is ON during [0, 1e-6] s.
T_BEAM_ON = 5e-7
T_BEAM_OFF = 1e-3

# Depletion threshold in the solver unit [µmol/L]
LIMIT_UMOL = DEPLETION_LIMIT * 1e6
# A concentration far below the threshold [µmol/L]
DEPLETED = 1e-3 * LIMIT_UMOL


@pytest.fixture(scope="module", params=list(CONFIGS))
def cell(request: pytest.FixtureRequest) -> UnitCell:
    out = UnitCell.from_toml(CONFIGS[request.param])
    if len(out.beam_list) == 0:
        out.set_beam(ConstantBeam(dose_rate=1e7, max_dose=10.0), LET=0.201)
    return out


def toy_env(config: dict[str, tp.Any]) -> SimEnv:
    return SimEnv.from_dict(config, log=False)  # type: ignore[arg-type]


@pytest.fixture
def mass_action_env() -> SimEnv:
    """Repeated reactant (2 A), a fixed species (F) and a negative product."""
    return toy_env(
        {
            "parameters": {"pH": 7.0},
            "concentrations": {
                "fixed": {"F": 1e-3},
                "initial": {"A": 1e-3, "B": 2e-3, "D": 1e-6},
            },
            "reactions": {
                "k_reaction": [
                    {"reaction": "A + B -> C - D", "k_value": 1e9},
                    {"reaction": "2 A + F -> E", "k_value": 2e8},
                    {"reaction": "C -> B + A", "k_value": 3e2},
                ],
                "radiolytic": {"A": 1.0},
            },
        }
    )


@pytest.fixture
def enzymatic_env() -> SimEnv:
    """S --enz>> P - D: the enzymatic reaction consumes D."""
    return toy_env(
        {
            "parameters": {"pH": 7.0},
            "concentrations": {"initial": {"S": 1e-3, "enz": 1e-6, "D": 1e-6}},
            "reactions": {
                "enzymatic": [
                    {"reaction": "S -- enz >> P - D", "k_value": 1e3, "k_micha": 1e-3},
                ],
                "radiolytic": {"S": 1.0},
            },
        }
    )


# --------------------------------- Helpers --------------------------------- #
def random_state(mat: SimMatrices, seed: int = 0) -> NDArray[np.float64]:
    """Random positive solver state [µmol/L], far from the depletion
    threshold (a few decades around 1 µmol/L)."""
    rng = np.random.default_rng(seed)
    n_dyn = mat.n_species - mat.state.n_cst
    return rng.lognormal(mean=0.0, sigma=1.5, size=n_dyn)


class FDJacobian(tp.NamedTuple):
    """Finite-difference Jacobian and the size of its rounding noise."""

    jac: NDArray[np.float64]
    noise: NDArray[np.float64]


def fd_jacobian(
    y: NDArray[np.float64],
    beam: tp.Optional[ParticleBeam],
    mat: SimMatrices,
    t: float,
    rel_step: float = 1e-3,
) -> FDJacobian:
    """Central finite-difference Jacobian of derive (relative step), with
    the rounding noise eps * |f| / h of every entry."""
    jac = np.zeros((y.size, y.size))
    noise = np.zeros((y.size, y.size))
    for j in range(y.size):
        h = rel_step * abs(y[j])
        up, down = y.copy(), y.copy()
        up[j] += h
        down[j] -= h
        f_up = derive(t, up, beam, mat)
        f_down = derive(t, down, beam, mat)
        jac[:, j] = (f_up - f_down) / (2 * h)
        noise[:, j] = 10 * np.finfo(float).eps * (np.abs(f_up) + np.abs(f_down)) / (2 * h)
    return FDJacobian(jac, noise)


def assert_jacobian_matches(
    analytic: NDArray[np.float64],
    numeric: FDJacobian,
    rtol: float = 1e-6,
) -> None:
    """Entry-wise match: the tolerance is rtol times the largest entry of
    the column, plus the rounding noise of the finite difference."""
    jac, noise = numeric
    assert analytic.shape == jac.shape
    tolerance = rtol * np.abs(jac).max(axis=0, keepdims=True) + noise
    excess = np.abs(analytic - jac) - tolerance
    assert (excess <= 0).all(), (
        f"{(excess > 0).sum()} entries differ, worst at "
        f"{np.unravel_index(excess.argmax(), excess.shape)}"
    )


# --------------------------------------------------------------------------- #
#                      Bundled configs vs finite differences                  #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("t", [T_BEAM_ON, T_BEAM_OFF], ids=["beam_on", "beam_off"])
def test_matches_finite_differences(cell: UnitCell, t: float) -> None:
    y = random_state(cell.mat)
    analytic = jacobian(t, y, cell.beam, cell.mat)
    numeric = fd_jacobian(y, cell.beam, cell.mat, t)
    assert_jacobian_matches(analytic, numeric)


def test_shape_and_pure(cell: UnitCell) -> None:
    mat = cell.mat
    y = random_state(mat)
    y_before = y.copy()
    n_dyn = mat.n_species - mat.state.n_cst
    jac = jacobian(T_BEAM_ON, y, cell.beam, mat)
    assert jac.shape == (n_dyn, n_dyn)
    assert jac.dtype == np.float64
    np.testing.assert_array_equal(y, y_before)
    np.testing.assert_array_equal(jac, jacobian(T_BEAM_ON, y, cell.beam, mat))


def test_independent_of_beam(cell: UnitCell) -> None:
    y = random_state(cell.mat)
    on = jacobian(T_BEAM_ON, y, cell.beam, cell.mat)
    off = jacobian(T_BEAM_OFF, y, cell.beam, cell.mat)
    no_beam = jacobian(T_BEAM_ON, y, None, cell.mat)
    np.testing.assert_array_equal(on, off)
    np.testing.assert_array_equal(on, no_beam)


def test_labarbe_has_enzymes_and_fixed_species() -> None:
    """The finite-difference tests above exercise both on this config."""
    mat = UnitCell.from_toml(rp.LABARBE_CONFIG_FILE).mat
    assert mat.enzymatic.n_reactions > 0
    assert mat.state.n_cst > 0


def test_depleted_species_still_consumed(cell: UnitCell) -> None:
    """Every solved species far below the threshold: the ones still
    consumed (by the constant species) keep a non-zero row."""
    mat = cell.mat
    y = DEPLETED * random_state(mat)
    dydt = derive(T_BEAM_OFF, y, cell.beam, mat)
    consumed = np.flatnonzero(dydt < 0)
    assert consumed.size >= 2

    analytic = jacobian(T_BEAM_OFF, y, cell.beam, mat)
    numeric = fd_jacobian(y, cell.beam, mat, T_BEAM_OFF)
    assert_jacobian_matches(analytic, numeric)
    assert analytic[consumed].any()


# --------------------------------------------------------------------------- #
#                              Toy environments                               #
# --------------------------------------------------------------------------- #
def test_repeated_reactant_and_fixed_species(mass_action_env: SimEnv) -> None:
    env = mass_action_env
    mat = env.build_all_matrices()
    y = random_state(mat)
    analytic = jacobian(0.0, y, None, mat)
    assert_jacobian_matches(analytic, fd_jacobian(y, None, mat, 0.0))

    # 2 A + F -> E, rate = k * A^2 * F, so d rate/d A = 2 k A F; E is
    # produced by this reaction only.
    kmat = mat.k_reactions
    a, e, f = (env.species[n].index for n in ("A", "E", "F"))
    r = int(np.flatnonzero(kmat.stoichiometry[e] > 0)[0])
    assert list(kmat.reactant_index[r]).count(a) == 2
    y_full = mat.state.to_full_molar(y)
    expected = 2 * kmat.k_values[r] * y_full[a] * y_full[f]
    assert analytic[e, a] == pytest.approx(expected, rel=1e-12)
    # the fixed species is not solved: no column for it
    n_dyn = mat.n_species - mat.state.n_cst
    assert mat.state.n_cst > 0
    assert analytic.shape == (n_dyn, n_dyn)


def test_negative_product_block_active_k_reaction(mass_action_env: SimEnv) -> None:
    env = mass_action_env
    mat = env.build_all_matrices()
    a, b, c, d = (env.species[n].index for n in ("A", "B", "C", "D"))
    y = random_state(mat)

    free = jacobian(0.0, y, None, mat)
    assert free[c, a] != 0.0
    assert free[c, b] != 0.0  # A + B -> C - D is active

    y[d] = DEPLETED  # D depleted: A + B -> C - D is blocked
    assert negative_product_block(mat.state.to_full_molar(y), mat.k_reactions).any()
    blocked = jacobian(0.0, y, None, mat)
    assert_jacobian_matches(blocked, fd_jacobian(y, None, mat, 0.0))
    # C is produced by A + B only (blocked): only its own decay remains
    assert blocked[c, a] == 0.0
    assert blocked[c, b] == 0.0
    assert blocked[c, c] != 0.0


def test_negative_product_block_active_enzymatic(enzymatic_env: SimEnv) -> None:
    env = enzymatic_env
    mat = env.build_all_matrices()
    d, p, s, enz = (env.species[n].index for n in ("D", "P", "S", "enz"))

    y = random_state(mat)
    free = jacobian(0.0, y, None, mat)
    assert_jacobian_matches(free, fd_jacobian(y, None, mat, 0.0))
    assert free[p, s] != 0.0
    assert free[p, enz] != 0.0

    y[d] = DEPLETED
    assert negative_product_block(mat.state.to_full_molar(y), mat.enzymatic).all()
    blocked = jacobian(0.0, y, None, mat)
    assert_jacobian_matches(blocked, fd_jacobian(y, None, mat, 0.0))
    assert not blocked.any()  # the only reaction is blocked


def test_michaelis_menten_derivatives(enzymatic_env: SimEnv) -> None:
    env = enzymatic_env
    mat = env.build_all_matrices()
    y = random_state(mat)
    y_full = mat.state.to_full_molar(y)
    s, enz, p = (env.species[n].index for n in ("S", "enz", "P"))
    k, km = mat.enzymatic.k_values[0], mat.enzymatic.k_michaelis[0]
    jac = jacobian(0.0, y, None, mat)
    assert jac[p, enz] == pytest.approx(k * y_full[s] / (km + y_full[s]), rel=1e-12)
    assert jac[p, s] == pytest.approx(
        k * y_full[enz] * km / (km + y_full[s]) ** 2, rel=1e-12
    )


def test_depleted_consumed_species_keeps_row(mass_action_env: SimEnv) -> None:
    env = mass_action_env
    mat = env.build_all_matrices()
    b, c = env.species["B"].index, env.species["C"].index
    y = random_state(mat)
    # B is depleted and consumed by A + B -> C - D (C -> B + A, its only
    # source, is off since C is depleted too)
    y[[b, c]] = DEPLETED

    jac = jacobian(0.0, y, None, mat)
    assert jac[b].any()
    assert_jacobian_matches(jac, fd_jacobian(y, None, mat, 0.0))
