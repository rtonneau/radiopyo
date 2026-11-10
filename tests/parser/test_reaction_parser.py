#!/usr/bin/env python3
"""
Tests of the reaction string parser: enzyme names of enzymatic reactions.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.species import BasicSpecies
from radiopyo.parser import RPARSER


@pytest.mark.parametrize(
    "reaction, enzyme, substrate",
    [
        ("2H2O2 -- catalase >> O2 + 2 H2O", "catalase", "H2O2"),
        ("O2_r_minus -- SOD >> H2O2", "SOD", "O2_r_minus"),
        ("H2O2 -- GPx >> H2O", "GPx", "H2O2"),
        ("H2O2 -- Cat2 >> H2O", "Cat2", "H2O2"),
        ("A -- E >> B", "E", "A"),
    ],
)
def test_enzyme_name_is_a_plain_species_named_after_the_text(
    reaction: str, enzyme: str, substrate: str
) -> None:
    gen = RPARSER.parse_reaction(reaction)
    assert type(gen["enzyme"]) is BasicSpecies
    assert gen["enzyme"].label == enzyme
    assert list(gen["reactants"]) == [substrate]


def test_enzyme_with_capitals_builds_a_michaelis_menten() -> None:
    from radiopyo.chemistry.michaelis import MichaelisMenten

    gen = RPARSER.parse_reaction("O2_r_minus -- SOD >> H2O2")
    gen["constants"] = {"k_value": 1e9, "k_micha": 1.0}
    reaction = MichaelisMenten.from_generic(gen)
    assert reaction.enzyme.label == "SOD"
    assert str(reaction) == "O2_r_minus -- SOD >> H2O2"
