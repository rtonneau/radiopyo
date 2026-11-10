#!/usr/bin/env python3
"""
Tests of the analytic Jacobian of the ODE right hand side: compared with
the central finite-difference Jacobian of ``derive`` on the bundled
configs (beam on/off, enzymes, fixed species) and on toy
environments (repeated reactant, negative products).

The rates are smooth (depletion switch c / (c + eps), ADR 0005) except at
the clipping kink c = 0 (ADR 0004): the finite differences are taken on
states away from 0, with a relative step small enough for the curvature of
the switch near eps.

The bundled configs have no negative product: ``negative_cell`` is the
Labarbe config plus a k-reaction with two negative products and an
enzymatic reaction with one.

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
    depletion_switch,
    derive,
    jacobian,
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


NEGATIVE_CONFIG = """\
[includes]
toml = ["radiopyo/config_Labarbe.toml"]

[[reactions.k_reaction]]
reaction = "e_aq + O2 -> O2_r_minus - H2O2 - OH_minus"
k_value = 1e9

[[reactions.enzymatic]]
reaction = "H2O2 -- SOD >> O2 - O2_r_minus"
k_value = 1e5
k_micha = 1e-3

[concentrations.fixed]
SOD = 1e-7
"""


@pytest.fixture(scope="module")
def negative_cell(tmp_path_factory: pytest.TempPathFactory) -> UnitCell:
    """Labarbe config with negative products (k: H2O2 and OH_minus,
    enzymatic: O2_r_minus)."""
    path = tmp_path_factory.mktemp("negative") / "negative.toml"
    path.write_text(NEGATIVE_CONFIG, encoding="utf-8")
    out = UnitCell.from_toml(path)
    assert out.mat.k_reactions.negative_products.sum() == 2
    assert out.mat.enzymatic.negative_products.sum() == 1
    return out


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


def test_negative_product_switch_k_reaction(mass_action_env: SimEnv) -> None:
    env = mass_action_env
    mat = env.build_all_matrices()
    a, b, c, d = (env.species[n].index for n in ("A", "B", "C", "D"))
    y = random_state(mat)

    free = jacobian(0.0, y, None, mat)
    assert free[c, a] != 0.0
    assert free[c, b] != 0.0  # A + B -> C - D is active
    # the switch of D (~1 umol/L >> eps) is 1 up to 1e-6
    free_switch = depletion_switch(mat.state.to_full_molar(y), mat.k_reactions).min()
    assert free_switch == pytest.approx(1.0, abs=1e-5)

    # D far below eps: A + B -> C - D is almost off but the Jacobian keeps the
    # (tiny) smooth dependence, no longer a hard zero
    y[d] = DEPLETED
    y_full = mat.state.to_full_molar(y)
    switch = depletion_switch(y_full, mat.k_reactions).min()
    assert 0.0 < switch < 1e-2
    slowed = jacobian(0.0, y, None, mat)
    assert_jacobian_matches(slowed, fd_jacobian(y, None, mat, 0.0, rel_step=1e-4))
    assert slowed[c, a] == pytest.approx(switch * free[c, a] / free_switch, rel=1e-12)
    assert slowed[c, b] == pytest.approx(switch * free[c, b] / free_switch, rel=1e-12)
    assert slowed[c, c] != 0.0
    # product rule: the switch itself depends on D (dC/dt = ... + k A B s(D))
    k_ab = 1e9
    expected = (
        k_ab * y_full[a] * y_full[b] * DEPLETION_LIMIT / (y_full[d] + DEPLETION_LIMIT) ** 2
    )
    assert slowed[c, d] == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("factor", [0.1, 1.0, 3.0, 100.0])
def test_switch_reactant_and_repeated_negative_product(factor: float) -> None:
    """A species that is both a reactant and a negative product of the
    same reaction, and a negative product of stoichiometry 2 next to a
    second one: one switch factor per negative species, rate and Jacobian
    consistent."""
    env = toy_env(
        {
            "parameters": {"pH": 7.0},
            "concentrations": {"initial": {"A": 1e-3, "B": 1e-6, "D": 1e-6}},
            "reactions": {
                "k_reaction": [
                    {"reaction": "A + D -> C - D", "k_value": 1e9},
                    {"reaction": "A -> C - 2 D - B", "k_value": 1e3},
                ],
                "radiolytic": {"A": 1.0},
            },
        }
    )
    mat = env.build_all_matrices()
    d = env.species["D"].index
    y = random_state(mat)
    y[d] = factor * LIMIT_UMOL
    y_full = mat.state.to_full_molar(y)
    switch = depletion_switch(y_full, mat.k_reactions)
    s_d = y_full[d] / (y_full[d] + DEPLETION_LIMIT)
    b = env.species["B"].index
    s_b = y_full[b] / (y_full[b] + DEPLETION_LIMIT)
    r1, r2 = (r.get_index() for r in env.reactions.iter_all_kreactions())
    assert switch[r1] == pytest.approx(s_d, rel=1e-12)
    assert switch[r2] == pytest.approx(s_d * s_b, rel=1e-12)  # 2 D: one factor
    analytic = jacobian(0.0, y, None, mat)
    assert_jacobian_matches(analytic, fd_jacobian(y, None, mat, 0.0, rel_step=1e-4))


@pytest.mark.parametrize("factor", [0.1, 1.0, 3.0, 100.0])
def test_switch_derivative_near_the_limit_k_reaction(
    mass_action_env: SimEnv, factor: float
) -> None:
    env = mass_action_env
    mat = env.build_all_matrices()
    d = env.species["D"].index
    y = random_state(mat)
    y[d] = factor * LIMIT_UMOL
    analytic = jacobian(0.0, y, None, mat)
    assert_jacobian_matches(analytic, fd_jacobian(y, None, mat, 0.0, rel_step=1e-4))


@pytest.mark.parametrize("factor", [0.1, 1.0, 3.0, 100.0])
def test_switch_derivative_near_the_limit_enzymatic(
    enzymatic_env: SimEnv, factor: float
) -> None:
    env = enzymatic_env
    mat = env.build_all_matrices()
    d, p, s, enz = (env.species[n].index for n in ("D", "P", "S", "enz"))
    y = random_state(mat)
    y[d] = factor * LIMIT_UMOL
    y_full = mat.state.to_full_molar(y)
    jac = jacobian(0.0, y, None, mat)
    assert_jacobian_matches(jac, fd_jacobian(y, None, mat, 0.0, rel_step=1e-4))

    k, km = mat.enzymatic.k_values[0], mat.enzymatic.k_michaelis[0]
    law = k * y_full[enz] * y_full[s] / (km + y_full[s])
    switch = y_full[d] / (y_full[d] + DEPLETION_LIMIT)
    assert switch == pytest.approx(factor / (1 + factor), rel=1e-12)
    # rate = law * switch: dP/dD = law * eps / (D + eps)^2 (product rule)
    expected = law * DEPLETION_LIMIT / (y_full[d] + DEPLETION_LIMIT) ** 2
    assert jac[p, d] == pytest.approx(expected, rel=1e-12)
    assert jac[p, enz] == pytest.approx(switch * law / y_full[enz], rel=1e-12)


def test_michaelis_menten_derivatives(enzymatic_env: SimEnv) -> None:
    env = enzymatic_env
    mat = env.build_all_matrices()
    y = random_state(mat)
    y_full = mat.state.to_full_molar(y)
    s, enz, p = (env.species[n].index for n in ("S", "enz", "P"))
    k, km = mat.enzymatic.k_values[0], mat.enzymatic.k_michaelis[0]
    jac = jacobian(0.0, y, None, mat)
    # D (~1 umol/L) is far above eps: the depletion switch is 1 up to 1e-6
    switch = depletion_switch(y_full, mat.enzymatic)[0]
    assert switch == pytest.approx(1.0, abs=1e-5)
    assert jac[p, enz] == pytest.approx(
        switch * k * y_full[s] / (km + y_full[s]), rel=1e-12
    )
    assert jac[p, s] == pytest.approx(
        switch * k * y_full[enz] * km / (km + y_full[s]) ** 2, rel=1e-12
    )


def test_depleted_consumed_species_keeps_row(mass_action_env: SimEnv) -> None:
    env = mass_action_env
    mat = env.build_all_matrices()
    b, c = env.species["B"].index, env.species["C"].index
    y = random_state(mat)
    # B is depleted and consumed by A + B -> C - D
    y[[b, c]] = DEPLETED

    jac = jacobian(0.0, y, None, mat)
    assert jac[b].any()
    assert_jacobian_matches(jac, fd_jacobian(y, None, mat, 0.0))


# --------------------------------------------------------------------------- #
#                     Clipped concentrations (ADR 0004)                       #
# --------------------------------------------------------------------------- #
def test_negative_species_contributes_zero(mass_action_env: SimEnv) -> None:
    env = mass_action_env
    mat = env.build_all_matrices()
    a, c = env.species["A"].index, env.species["C"].index

    y = random_state(mat)
    y[a] = -y[a]  # a solver overshoot below 0
    jac = jacobian(0.0, y, None, mat)
    assert not jac[:, a].any()  # no derivative depends on A (also 2 A + F -> E)
    assert jac[:, c].any()  # C -> B + A is unaffected
    assert_jacobian_matches(jac, fd_jacobian(y, None, mat, 0.0))


def test_zero_species_uses_positive_side(mass_action_env: SimEnv) -> None:
    env = mass_action_env
    mat = env.build_all_matrices()
    a, b, c, e = (env.species[n].index for n in ("A", "B", "C", "E"))

    y = random_state(mat)
    y[a] = 0.0
    y_full = mat.state.to_full_molar(y)
    jac = jacobian(0.0, y, None, mat)
    k_ab = 1e9  # A + B -> C - D: rate = k A B s(D), d rate / dA = k B s(D) at A = 0
    switch = depletion_switch(y_full, mat.k_reactions).min()
    assert jac[c, a] == pytest.approx(k_ab * y_full[b] * switch, rel=1e-12)
    assert jac[c, b] == 0.0  # d rate / dB = k A = 0
    assert jac[e, a] == 0.0  # 2 A + F -> E: d rate / dA = 2 k A F = 0


def test_negative_enzymatic_species_contribute_zero(enzymatic_env: SimEnv) -> None:
    env = enzymatic_env
    mat = env.build_all_matrices()
    p, s, enz = (env.species[n].index for n in ("P", "S", "enz"))
    k, km = mat.enzymatic.k_values[0], mat.enzymatic.k_michaelis[0]

    y = random_state(mat)
    y[s] = -y[s]
    jac = jacobian(0.0, y, None, mat)
    assert not jac[:, s].any()
    assert jac[p, enz] == 0.0  # clipped substrate = 0: no rate
    assert_jacobian_matches(jac, fd_jacobian(y, None, mat, 0.0))

    y = random_state(mat)
    y[enz] = -y[enz]
    jac = jacobian(0.0, y, None, mat)
    assert not jac[:, enz].any()
    assert jac[p, s] == 0.0  # clipped enzyme = 0: no rate
    assert_jacobian_matches(jac, fd_jacobian(y, None, mat, 0.0))

    # exactly 0: the derivative of the positive side
    y = random_state(mat)
    y[s] = 0.0
    jac = jacobian(0.0, y, None, mat)
    y_full = mat.state.to_full_molar(y)
    switch = depletion_switch(y_full, mat.enzymatic)[0]
    assert jac[p, s] == pytest.approx(switch * k * y_full[enz] / km, rel=1e-12)
    assert jac[p, enz] == 0.0


# --------------------------------------------------------------------------- #
#               Negative products on a bundled config (ADR 0005)              #
# --------------------------------------------------------------------------- #
def _negative_species(cell: UnitCell) -> list[int]:
    """Dynamic species that are negative products (k or enzymatic)."""
    mat = cell.mat
    n_dyn = mat.n_species - mat.state.n_cst
    flags = mat.k_reactions.negative_products.any(axis=1)
    flags |= mat.enzymatic.negative_products.any(axis=1)
    return [int(i) for i in np.flatnonzero(flags) if i < n_dyn]


@pytest.mark.parametrize("seed", range(3))
def test_negative_cell_random_positive_states(
    negative_cell: UnitCell, seed: int
) -> None:
    mat = negative_cell.mat
    y = random_state(mat, seed=seed)
    analytic = jacobian(T_BEAM_ON, y, negative_cell.beam, mat)
    numeric = fd_jacobian(y, negative_cell.beam, mat, T_BEAM_ON)
    assert_jacobian_matches(analytic, numeric)


@pytest.mark.parametrize("factor", [0.1, 1.0, 3.0, 100.0])
def test_negative_cell_negative_products_near_the_limit(
    negative_cell: UnitCell, factor: float
) -> None:
    mat = negative_cell.mat
    species = _negative_species(negative_cell)
    assert len(species) == 3  # H2O2, OH_minus, O2_r_minus
    y = random_state(mat, seed=7)
    y[species] = factor * LIMIT_UMOL * np.array([1.0, 0.5, 2.0])
    switches = depletion_switch(mat.state.to_full_molar(y), mat.k_reactions)
    assert 0.0 < switches.min() < 1.0
    analytic = jacobian(T_BEAM_OFF, y, negative_cell.beam, mat)
    numeric = fd_jacobian(y, negative_cell.beam, mat, T_BEAM_OFF, rel_step=1e-4)
    assert_jacobian_matches(analytic, numeric)


def test_negative_cell_with_clipped_species(negative_cell: UnitCell) -> None:
    mat = negative_cell.mat
    rng = np.random.default_rng(11)
    for _ in range(3):
        y = random_state(mat, seed=int(rng.integers(1000)))
        flip = rng.random(y.size) < 0.3
        flip[_negative_species(negative_cell)[0]] = True  # a negative product
        y[flip] *= -1.0
        analytic = jacobian(T_BEAM_ON, y, negative_cell.beam, mat)
        assert not analytic[:, flip].any()  # clipped: no dependence below 0
        numeric = fd_jacobian(y, negative_cell.beam, mat, T_BEAM_ON)
        assert_jacobian_matches(analytic, numeric)
