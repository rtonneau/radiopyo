#!/usr/bin/env python3
"""
Tests of the gross production / consumption integrals of ODEResult
(integrate_production_of, integrate_consumption_of), equilibrium reactions
included: an equilibrium counts as production when its net integral
produces the species and as consumption otherwise, never both.
Consumption integrals are negative.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.physics.beam import SinglePulseBeam
from radiopyo.simulation.ode_result import ODEResult
from radiopyo.simulation.unit_cell import UnitCell


# ------------------------------- Fixtures ---------------------------------- #
@pytest.fixture(scope="module")
def result() -> ODEResult:
    cell = UnitCell.from_toml(rp.PURE_WATER_CONFIG_FILE)
    cell.env.species["O2"].set_initial_cc(65e-6)
    cell.set_beam(SinglePulseBeam(dose_rate=1e3, max_dose=1.0), LET=1.17)
    return cell.prepare_single_pulse(t_span=(0.0, 1.0)).run()


# --------------------------------- Tests ----------------------------------- #
@pytest.mark.parametrize("species", ["H2O2", "O2", "OH_minus"])
def test_consumption_with_equilibrium_does_not_raise(result, species):
    """Several equilibria involve OH_minus: dropping the producing ones
    while iterating over them raised under pandas 3."""
    cons = result.integrate_consumption_of(species)
    assert np.all(np.isfinite(cons.to_numpy()))
    # consumption integrals are negative by convention
    assert np.all(cons.to_numpy() <= 0.0)


@pytest.mark.parametrize("species", ["H2O2", "OH_minus"])
def test_equilibrium_is_either_production_or_consumption(result, species):
    prod = result.integrate_production_of(species)
    cons = result.integrate_consumption_of(species)
    eq_labels = set(result.compute_yields_of_eq_reactions_by_species(species).columns)
    assert eq_labels, "the test needs at least one equilibrium"
    assert not (eq_labels & set(prod.index) & set(cons.index))


@pytest.mark.parametrize("species", ["H2O2", "O2"])
def test_production_minus_consumption_is_the_net_change(result, species):
    df = result.to_pandas()[species]
    net = df.iloc[-1] - df.iloc[0]
    # The initial burst is added to y0 before the first stored point, so
    # it counts in the production but not in the stored difference.
    burst = sum(
        r.compute_initial_cc_before(result.t_zero_pulse, dr) * 1e6
        for r in result.env.radiolytic_reaction_involving(species)
        for dr in result.burst_dose_rates or []
    )
    # consumption integrals are negative: their sum adds to the production
    gross = (
        result.integrate_production_of(species).sum()
        + result.integrate_consumption_of(species).sum()
    )
    assert gross == pytest.approx(net + burst, rel=2e-2, abs=1e-3)
