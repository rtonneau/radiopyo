#!/usr/bin/env python3
"""
Tests for radiopyo.chemistry.species module

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 21/08/2024
"""
import pytest
import numpy as np

import radiopyo.chemistry.species as sp
from radiopyo.chemistry.utils.units import Temperature


# ----------------------------- BasicSpecies tests ---------------------------- #
@pytest.fixture
def raw_H2O2():
    """
    GEnerate a unique BasicSpecies instance.
    """
    return sp.BasicSpecies("H2O2")


def test_raw_species_label(raw_H2O2):
    assert raw_H2O2.label == "H2O2"


def test_raw_species_eq(raw_H2O2):
    other = sp.BasicSpecies("H2O2")
    assert raw_H2O2 == other
    assert raw_H2O2 == "H2O2"
    other = sp.BasicSpecies("H2O3")
    assert raw_H2O2 != other
    assert raw_H2O2 != "H2O3"


def test_raw_species_eq_unsupported_type(raw_H2O2):
    # BasicSpecies no longer supports `in` (former "contains" test): the
    # equality with str / BasicSpecies is the remaining comparison API.
    # Foreign types: __eq__ returns NotImplemented -> Python falls back
    # to identity, i.e. not equal (no exception).
    assert (raw_H2O2 == 3.2) is False
    assert raw_H2O2 != 3.2
    assert raw_H2O2.__eq__(3.2) is NotImplemented


def test_raw_species_hash(raw_H2O2):
    assert hash(raw_H2O2) == hash(sp.BasicSpecies("H2O2"))
    assert len({raw_H2O2, sp.BasicSpecies("H2O2"), sp.BasicSpecies("H2")}) == 2


def test_raw_species_string(raw_H2O2):
    assert str(raw_H2O2) == "H2O2"


# ---------------------------- SimSpecies tests ---------------------------- #
@pytest.fixture
def H2():
    """
    Generate a unique SimSpecies instance.
    """
    return sp.SimSpecies(sp.BasicSpecies("H2"))


@pytest.fixture
def H2O2():
    """
    Generate a unique SimSpecies instance.
    """
    return sp.SimSpecies(sp.BasicSpecies("H2O2"), index=4, cc_init=25.4)


def test_simspecies_init(H2O2, H2O, H2):
    assert H2.label == "H2"
    assert H2.index == -1
    assert H2.cc_init is None  # empty initial concentration
    assert H2O2.label == "H2O2"
    assert H2O2.index == 4
    assert H2O2.cc_init == np.polynomial.Polynomial(25.4)
    assert H2O.label == "H2O"
    assert H2O.index == -1
    assert H2O.cc_init == np.polynomial.Polynomial(
        [55.50, 6.075e-3, 4.110e-4, 1.496e-6, 2.619e-9]
    )


def test_initial_cc(H2O2, H2O, H2):
    assert H2O2.initial_cc() == 25.4
    assert H2O2.initial_cc(Temperature(298.15)) == 25.4
    assert H2.initial_cc() == 0


def test_set_cc_init(H2):
    H2.set_initial_cc(42)
    assert H2.cc_init == np.polynomial.Polynomial(42)
    assert H2.initial_cc() == 42
    H2.set_initial_cc([1, 2, 1, 1])
    assert H2.cc_init == np.polynomial.Polynomial([1, 2, 1, 1])
    assert H2.initial_cc(0.0) == 1.0
    assert H2.initial_cc() == 16301.0  # Then default T = 25 °C


def test_default_water(H2O):
    assert sp.H2O.initial_cc(0) == 55.5
    assert sp.H2O.initial_cc() == H2O.initial_cc()
    assert sp.get_H2O_cc() == H2O.initial_cc()
