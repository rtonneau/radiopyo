#!/usr/bin/env python3
"""
Tests of the solver core: the frozen matrices (SimMatrices and its
parts), the single reaction rate function shared by derive() and
ODEResult.reconstruct_reactions(), and the µmol/L <-> mol/L state helper.

Expected values are rebuilt from the reaction objects of the SimEnv
(reactants, products, k values, ...), not from the matrices under test.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.chemistry.basic_reaction import Operation
from radiopyo.physics.beam import ConstantBeam
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.ode_functions_solver import (
    DEPLETION_LIMIT,
    depletion_switch,
    derive,
    dose_rates,
    mass_action_rates,
    michaelis_menten_rates,
    reaction_rates,
    species_rates,
)
from radiopyo.simulation.ode_result import DummyScipyOdeResult, ODEResult
from radiopyo.simulation.unit_cell import UnitCell
from radiopyo.utils.sim_types import (
    EzMatrices,
    GVector,
    KMatrices,
    SimMatrices,
    StateLayout,
)

# ------------------------------- Fixtures ---------------------------------- #
CONFIGS = {
    "labarbe": rp.LABARBE_CONFIG_FILE,
    "water": rp.PURE_WATER_CONFIG_FILE,
}

# The Labarbe config defines a constant beam of 1e7 Gy/s up to 10 Gy, i.e.
# ON during [0, 1e-6] s.
T_BEAM_ON = 5e-7
T_BEAM_OFF = 1e-3


@pytest.fixture(scope="module", params=list(CONFIGS))
def cell(request: pytest.FixtureRequest) -> UnitCell:
    """UnitCell of a bundled config (not to be modified by the tests)."""
    out = UnitCell.from_toml(CONFIGS[request.param])
    if len(out.beam_list) == 0:
        out.set_beam(ConstantBeam(dose_rate=1e7, max_dose=10.0), LET=0.201)
    return out


@pytest.fixture(scope="module")
def labarbe_cell() -> UnitCell:
    return UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)


def random_full_state(
    mat: SimMatrices,
    rng: np.random.Generator,
    n_t: int | None = None,
) -> np.ndarray:
    """Random positive full molar state [mol/L] (constant species at the
    end, as built by StateLayout)."""
    n_dyn = mat.n_species - mat.state.n_cst
    shape = (n_dyn,) if n_t is None else (n_dyn, n_t)
    y = rng.lognormal(mean=-12, sigma=3, size=shape) * 1e6  # µmol/L
    return mat.state.to_full_molar(y)


@pytest.fixture
def negative_product_env() -> SimEnv:
    """Small env with a negative product (D) and a fixed reactant (F)."""
    config = {
        "parameters": {"pH": 7.0},
        "concentrations": {
            "fixed": {"F": 1e-3},
            "initial": {"A": 1e-3, "B": 2e-3, "D": 1e-6},
        },
        "reactions": {
            "k_reaction": [
                {"reaction": "A + B -> C - D", "k_value": 1e9},
                {"reaction": "2 A + F -> E", "k_value": 2e8},
            ],
            "radiolytic": {"A": 1.0},
        },
    }
    return SimEnv.from_dict(config, log=False)  # type: ignore[arg-type]


# --------------------------------------------------------------------------- #
#                               Frozen matrices                               #
# --------------------------------------------------------------------------- #
def test_matrices_are_frozen_and_read_only(labarbe_cell: UnitCell) -> None:
    mat = labarbe_cell.mat
    with pytest.raises(dataclasses.FrozenInstanceError):
        mat.k_reactions = mat.k_reactions  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        mat.k_reactions.k_values = np.zeros(1)  # type: ignore[misc]
    for part in (mat.k_reactions, mat.radiolytic, mat.enzymatic, mat.state):
        for field in dataclasses.fields(part):
            array = getattr(part, field.name)
            assert not array.flags.writeable, f"{type(part).__name__}.{field.name}"
            with pytest.raises(ValueError):
                array[(0,) * array.ndim] = 1


def test_matrices_shapes(cell: UnitCell) -> None:
    env, mat = cell.env, cell.mat
    n_sp = len(env.species)
    n_kr = env.reactions.total_kreactions()
    n_g = len(env.reactions.radiolytic)
    n_ez = len(env.reactions.enzymatic)

    assert mat.n_species == n_sp
    assert mat.state.n_cst == env.species.len_cst_species()
    assert mat.k_reactions.n_reactions == n_kr
    assert mat.k_reactions.stoichiometry.shape == (n_sp, n_kr)
    assert mat.k_reactions.negative_products.shape == (n_sp, n_kr)
    assert mat.k_reactions.reactant_index.shape[0] == n_kr
    assert mat.radiolytic.stoichiometry.shape == (n_sp, n_g)
    assert mat.enzymatic.n_reactions == n_ez
    assert mat.enzymatic.stoichiometry.shape == (n_sp, n_ez)


def test_k_stoichiometry_matches_reactions(cell: UnitCell) -> None:
    env, kmat = cell.env, cell.mat.k_reactions
    for _, reaction, k_value in env.reactions.iter_all_kreactions_and_kvalues(cell._T):
        r = reaction.get_index()
        expected = np.zeros(len(env.species))
        for r_sp in reaction.reactants:
            if not env.species.is_cst_species(r_sp.species):
                expected[r_sp.species.index] -= r_sp.stoi
        for r_sp in reaction.products:
            if env.species.is_cst_species(r_sp.species):
                continue
            sign = 1 if r_sp.ope is Operation.ADD else -1
            expected[r_sp.species.index] += sign * r_sp.stoi
        np.testing.assert_array_equal(kmat.stoichiometry[:, r], expected, str(reaction))
        assert kmat.k_values[r] == k_value


def test_radiolytic_vector_matches_reactions(cell: UnitCell) -> None:
    env, gvec = cell.env, cell.mat.radiolytic
    for col, reaction in enumerate(env.reactions.radiolytic):
        assert gvec.yields[col] == reaction.kr()
        expected = np.zeros(len(env.species))
        if not env.species.is_cst_species(reaction.species()):
            expected[reaction.species().index] = 1.0
        np.testing.assert_array_equal(gvec.stoichiometry[:, col], expected)


# --------------------------------------------------------------------------- #
#                                 StateLayout                                 #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("cst_sp", [[], [1e-3, 55.5]], ids=["no-fixed", "2-fixed"])
@pytest.mark.parametrize("n_t", [None, 4], ids=["1-D", "2-D"])
def test_state_layout_round_trip(cst_sp: list, n_t: int | None) -> None:
    layout = StateLayout(cst_sp=np.array(cst_sp, dtype=np.float64))
    assert layout.n_cst == len(cst_sp)

    shape = (3,) if n_t is None else (3, n_t)
    y = np.arange(1.0, 1.0 + math.prod(shape)).reshape(shape)  # µmol/L
    full = layout.to_full_molar(y)

    assert full.shape == (3 + len(cst_sp), *shape[1:])
    np.testing.assert_array_equal(full[:3], y * 1.0e-6)
    expected_cst = np.array(cst_sp).reshape(-1, *([1] * (len(shape) - 1)))
    np.testing.assert_array_equal(
        full[3:], np.broadcast_to(expected_cst, full[3:].shape)
    )

    back = layout.to_dynamic_micromolar(full)
    assert back.shape == y.shape
    np.testing.assert_allclose(back, y, rtol=1e-15)


def test_derive_with_zero_fixed_species() -> None:
    """2 A -> B with no constant species at all."""
    k = 3.0e9
    kmat = KMatrices(
        k_values=np.array([k]),
        reactant_index=np.array([[0, 0]], dtype=np.intp),
        stoichiometry=np.array([[-2.0], [1.0]]),
        negative_products=np.zeros((2, 1), dtype=np.bool_),
    )
    mat = SimMatrices(
        k_reactions=kmat,
        radiolytic=GVector(yields=np.zeros(0), stoichiometry=np.zeros((2, 0))),
        enzymatic=EzMatrices(
            k_values=np.zeros(0),
            k_michaelis=np.zeros(0),
            substrate_index=np.zeros(0, dtype=np.intp),
            enzyme_index=np.zeros(0, dtype=np.intp),
            stoichiometry=np.zeros((2, 0)),
            negative_products=np.zeros((2, 0), dtype=np.bool_),
        ),
        state=StateLayout(cst_sp=np.zeros(0)),
    )
    y = np.array([2.0, 5.0])  # µmol/L
    rate = k * (2.0e-6) ** 2  # mol/L/s
    np.testing.assert_allclose(
        derive(0.0, y, None, mat),
        np.array([-2 * rate, rate]) * 1e6,
        rtol=1e-14,
    )


# --------------------------------------------------------------------------- #
#                               Reaction rates                                #
# --------------------------------------------------------------------------- #
def test_k_rates_follow_mass_action_law(cell: UnitCell) -> None:
    env, mat = cell.env, cell.mat
    y_full = random_full_state(mat, np.random.default_rng(0))
    rates = mass_action_rates(y_full, mat.k_reactions)

    for _, reaction, k_value in env.reactions.iter_all_kreactions_and_kvalues(cell._T):
        factors = [
            y_full[r_sp.species.index]
            for r_sp in reaction.reactants
            for _ in range(r_sp.stoi)
        ]
        expected = math.prod(factors) * k_value
        assert rates[reaction.get_index()] == pytest.approx(expected, rel=1e-14), str(
            reaction
        )


def test_michaelis_menten_rates(labarbe_cell: UnitCell) -> None:
    env, mat = labarbe_cell.env, labarbe_cell.mat
    assert len(env.reactions.enzymatic) > 0
    y_full = random_full_state(mat, np.random.default_rng(1))
    rates = michaelis_menten_rates(y_full, mat.enzymatic)

    for col, reaction in enumerate(env.reactions.enzymatic):
        enzyme = y_full[reaction.enzyme.index]
        substrate = y_full[reaction.substrate.index]
        expected = (
            reaction.k_value * enzyme * substrate / (reaction.k_micha + substrate)
        )
        assert rates[col] == pytest.approx(expected, rel=1e-14)


def test_dose_rates(labarbe_cell: UnitCell) -> None:
    beam = labarbe_cell.beam
    times = np.array([0.0, T_BEAM_ON, T_BEAM_OFF])
    expected = np.array([beam.at(t).dose_rate() for t in times])
    assert expected[1] > 0.0
    assert expected[2] == 0.0

    np.testing.assert_array_equal(dose_rates(times, beam), expected)
    assert dose_rates(T_BEAM_ON, beam) == expected[1]
    assert dose_rates(T_BEAM_ON, None) == 0.0
    np.testing.assert_array_equal(dose_rates(times, None), np.zeros(3))


@pytest.mark.parametrize("t", [T_BEAM_ON, T_BEAM_OFF], ids=["beam-on", "beam-off"])
def test_derive_is_stoichiometry_times_rates(cell: UnitCell, t: float) -> None:
    mat, beam = cell.mat, cell.beam
    y = (
        cell.env.initial_values()
        * np.random.default_rng(2).lognormal(0, 1, mat.n_species - mat.state.n_cst)
        + 1e-3
    )
    y_full = mat.state.to_full_molar(y)
    rates = reaction_rates(t, y_full, mat, beam)

    np.testing.assert_array_equal(
        rates.radiolytic, mat.radiolytic.yields * beam.at(t).dose_rate()
    )
    expected = (
        mat.k_reactions.stoichiometry @ rates.k_reactions
        + mat.radiolytic.stoichiometry @ rates.radiolytic
        + mat.enzymatic.stoichiometry @ rates.enzymatic
    )
    n_dyn = mat.n_species - mat.state.n_cst
    assert np.all(expected[n_dyn:] == 0.0)  # constant species do not move
    np.testing.assert_allclose(
        species_rates(rates, mat), expected, rtol=1e-12, atol=1e-300
    )

    dydt = derive(t, y, beam, mat)
    assert dydt.shape == y.shape
    scale = np.max(np.abs(expected))
    np.testing.assert_allclose(
        dydt, expected[:n_dyn] * 1e6, rtol=1e-12, atol=1e-12 * scale * 1e6
    )


def test_reconstruct_reactions_matches_rates_per_time_point(
    labarbe_cell: UnitCell,
) -> None:
    env, mat, beam = labarbe_cell.env, labarbe_cell.mat, labarbe_cell.beam
    n_t = 6
    rng = np.random.default_rng(3)
    times = np.array([1e-7, 3e-7, T_BEAM_ON, 9e-7, 2e-6, T_BEAM_OFF])
    full = random_full_state(mat, rng, n_t=n_t)
    n_dyn = mat.n_species - mat.state.n_cst
    y = full[:n_dyn] * 1e6  # µmol/L
    y[0, 1] = 0.0  # a depleted species

    res = ODEResult(
        env,
        beam,
        mat,
        labarbe_cell.beam_list.current_LET,
        labarbe_cell.T_ZERO,
        DummyScipyOdeResult(t=times, y=y),
    )
    df = res.reconstruct_reactions(recompute=True)

    # Column order: k-reactions, radiolytic, then enzymatic reactions
    expected_columns = (
        [r.as_label() for r in env.reactions.iter_all_kreactions()]
        + [r.as_label() for r in env.reactions.radiolytic]
        + [r.as_label() for r in env.reactions.enzymatic]
    )
    assert list(df.columns) == expected_columns
    np.testing.assert_array_equal(df.index.values, times)

    for i, t in enumerate(times):
        y_full = mat.state.to_full_molar(y[:, i])
        expected = reaction_rates(t, y_full, mat, beam).stacked()
        np.testing.assert_allclose(df.iloc[i].to_numpy(), expected, rtol=1e-14, atol=0)

    # Radiolytic columns are 0 once the beam is OFF.
    g_cols = [r.as_label() for r in env.reactions.radiolytic]
    assert np.all(df.loc[T_BEAM_OFF, g_cols] == 0.0)
    assert np.all(df.loc[T_BEAM_ON, g_cols] > 0.0)

    # Cached result is returned when not recomputing
    np.testing.assert_array_equal(res.reconstruct_reactions().to_numpy(), df.to_numpy())


def test_reconstruct_reactions_without_time_points(labarbe_cell: UnitCell) -> None:
    res = ODEResult(
        labarbe_cell.env, labarbe_cell.beam, labarbe_cell.mat, None, labarbe_cell.T_ZERO
    )
    df = res.reconstruct_reactions()
    assert df.shape[0] == 0
    assert df.shape[1] == labarbe_cell.env.reactions.total_kreactions() + len(
        labarbe_cell.env.reactions.radiolytic
    ) + len(labarbe_cell.env.reactions.enzymatic)


# --------------------------------------------------------------------------- #
#                    Negative products / depletion switch                     #
# --------------------------------------------------------------------------- #
def test_negative_product_depletion_switch(negative_product_env: SimEnv) -> None:
    env = negative_product_env
    mat = env.build_all_matrices()
    kmat = mat.k_reactions
    idx = {label: env.species[label].index for label in "ABDF"}
    r_neg = next(r for r in env.reactions.iter_all_kreactions() if r.has_reactant("B"))
    r_cst = next(r for r in env.reactions.iter_all_kreactions() if r.has_reactant("F"))

    assert kmat.negative_products[idx["D"], r_neg.get_index()]
    assert kmat.negative_products.sum() == 1
    assert kmat.stoichiometry[idx["D"], r_neg.get_index()] == -1  # D is consumed
    assert kmat.stoichiometry[idx["F"], r_cst.get_index()] == 0  # F is fixed

    def state(d: float) -> np.ndarray:
        y_full = np.full(mat.n_species, 1e-4)
        y_full[-mat.state.n_cst :] = mat.cst_sp
        y_full[idx["D"]] = d
        return y_full

    # Switch c / (c + eps): 0 at 0 (and below), 1/2 at eps, ~1 well above
    for d, switch in [
        (0.0, 0.0),
        (-1e-9, 0.0),
        (DEPLETION_LIMIT, 0.5),
        (2 * DEPLETION_LIMIT, 2 / 3),
        (1e-10, 1e-10 / (1e-10 + DEPLETION_LIMIT)),
        (1e-6, 1e-6 / (1e-6 + DEPLETION_LIMIT)),
    ]:
        y_full = state(d)
        factors = depletion_switch(y_full, kmat)
        assert factors[r_neg.get_index()] == pytest.approx(switch, rel=1e-14)
        assert factors[r_cst.get_index()] == 1.0  # no negative product
        rates = reaction_rates(0.0, y_full, mat, None).k_reactions
        expected_neg = switch * 1e9 * 1e-4 * 1e-4
        assert rates[r_neg.get_index()] == pytest.approx(expected_neg, rel=1e-14)
        # The other reaction (with a fixed reactant F) is never slowed down
        expected_cst = (
            2e8 * 1e-4 * 1e-4 * mat.cst_sp[idx["F"] - (mat.n_species - mat.state.n_cst)]
        )
        assert rates[r_cst.get_index()] == pytest.approx(expected_cst, rel=1e-14)

    # Almost no slowdown at 1e-6 mol/L (< 1e-5 relative), under 1 % at 1e-10
    assert depletion_switch(state(1e-6), kmat)[r_neg.get_index()] > 1 - 1e-5
    assert depletion_switch(state(1e-10), kmat)[r_neg.get_index()] > 0.99

    # Vectorised over time points: same factors column by column
    ds = (0.0, DEPLETION_LIMIT, 1e-6)
    states = np.stack([state(d) for d in ds], axis=1)
    np.testing.assert_allclose(
        depletion_switch(states, kmat)[r_neg.get_index()],
        [0.0, 0.5, 1e-6 / (1e-6 + DEPLETION_LIMIT)],
        rtol=1e-14,
    )
    np.testing.assert_array_equal(depletion_switch(states, kmat)[r_cst.get_index()], 1.0)
    vec = reaction_rates(np.zeros(3), states, mat, None).k_reactions
    for col, d in enumerate(ds):
        np.testing.assert_array_equal(
            vec[:, col], reaction_rates(0.0, state(d), mat, None).k_reactions
        )


# --------------------------------------------------------------------------- #
#                     Clipped concentrations (ADR 0004)                       #
# --------------------------------------------------------------------------- #
def _overshoot(
    y_full: np.ndarray, mat: SimMatrices, rng: np.random.Generator
) -> np.ndarray:
    """Copy of a full state with a third of the dynamic species pushed below
    0 (as a solver overshoot would do); constant species are untouched."""
    out = y_full.copy()
    n_dyn = mat.n_species - mat.state.n_cst
    flip = rng.random(n_dyn) < 0.33
    out[:n_dyn][flip] *= -1.0
    return out


def test_rates_read_clipped_concentrations(cell: UnitCell) -> None:
    mat, beam = cell.mat, cell.beam
    rng = np.random.default_rng(3)
    y_neg = _overshoot(random_full_state(mat, rng), mat, rng)
    assert (y_neg < 0).any()
    y_clip = np.maximum(y_neg, 0.0)

    mass_action = mass_action_rates(y_neg, mat.k_reactions)
    np.testing.assert_array_equal(
        mass_action, mass_action_rates(y_clip, mat.k_reactions)
    )
    assert (mass_action >= 0.0).all()
    michaelis = michaelis_menten_rates(y_neg, mat.enzymatic)
    np.testing.assert_array_equal(
        michaelis, michaelis_menten_rates(y_clip, mat.enzymatic)
    )
    assert (michaelis >= 0.0).all()

    rates_neg = reaction_rates(T_BEAM_ON, y_neg, mat, beam)
    rates_clip = reaction_rates(T_BEAM_ON, y_clip, mat, beam)
    np.testing.assert_array_equal(rates_neg.stacked(), rates_clip.stacked())

    # Same for several states at once (n_species, n_t)
    states = np.stack([y_neg, y_clip], axis=1)
    both = reaction_rates(np.array([T_BEAM_ON, T_BEAM_ON]), states, mat, beam)
    np.testing.assert_array_equal(both.stacked()[:, 0], both.stacked()[:, 1])
    np.testing.assert_array_equal(both.stacked()[:, 0], rates_neg.stacked())


def test_michaelis_menten_negative_enzyme_or_substrate(labarbe_cell: UnitCell) -> None:
    mat = labarbe_cell.mat
    ezmat = mat.enzymatic
    y_full = random_full_state(mat, np.random.default_rng(4))
    for index in (ezmat.enzyme_index, ezmat.substrate_index):
        y_neg = y_full.copy()
        y_neg[index] *= -1.0
        np.testing.assert_array_equal(michaelis_menten_rates(y_neg, ezmat), 0.0)
