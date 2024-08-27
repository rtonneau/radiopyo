#!/usr/bin/env python3
"""
Provides common fixtures for chemistry tests.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 27/08/2024
"""
import pytest
import numpy as np
from radiopyo.chemistry.species_collection import ListSpecies
from radiopyo.chemistry.species import RawSpecies, ChemSpecies


@pytest.fixture
def raw_species():
    """
    Generate a list of RawSpecies instances.
    """
    species_list = ["OH_minus", "H3O_plus", "H_plus", "OH", "ROO", "dyn"]
    return [RawSpecies(label) for label in species_list]


@pytest.fixture
def H2O():
    """
    Generate a unique ChemSpecies instance.
    """
    return ChemSpecies(
        label="H2O",
        cc_init=[55.50, 6.075e-3, 4.110e-4, 1.496e-6, 2.619e-9],
    )


@pytest.fixture
def species(raw_species):
    """
    Generate a list of ChemSpecies instances.
    """
    return [ChemSpecies(label=elt.label) for elt in raw_species]


@pytest.fixture
def empty_list_species():
    """Provides an empty ListSpecies"""
    return ListSpecies()


@pytest.fixture
def non_existing_species():
    """Provides an empty ListSpecies"""
    return ChemSpecies("non_existing")


@pytest.fixture
def dyn_species(species):
    """Provides the list of dynamical species used as test"""
    return species


@pytest.fixture
def cst_species(H2O):
    """Provides the list of constant species used as test"""
    return [H2O, ChemSpecies("cst_species")]


@pytest.fixture
def list_species(dyn_species, cst_species):
    """
    Provides a ListSpecies populated with with 3 dyn species and 1 cst
    species:
        dyn species:
            OH_minus, H3O_plus, H_plus

        cst species:
            H2O, cst_species
    """
    list_species = ListSpecies()
    # First add cst species so that it is unordered by default
    for sp in cst_species:
        list_species.add_cst_species(sp.label, 0)
    for sp in dyn_species:
        list_species.add_dyn_species(sp.label)
    return list_species
