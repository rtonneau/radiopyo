import pytest
import numpy as np
from numpy.polynomial import Polynomial

from radiopyo.chemistry import k_reactions as kr
from radiopyo.chemistry.species import BasicSpecies, SimSpecies
from radiopyo.chemistry.utils.units import Temperature
from tests.chemistry.test_basic_chemistry import TestBasicReaction

param_kvalues = pytest.mark.parametrize(
    "k_value, temp, result",
    [
        [[13.339, -2.220e3, 7.333e5, -1.065e8], None, 13303373557.26],
        [[13.339, -2.220e3, 7.333e5, -1.065e8], 25, 13303373557.26],
        [
            [13.339, -2.220e3, 7.333e5, -1.065e8],
            Temperature(298.15, "kelvin"),
            13303373557.26,
        ],
        [
            Polynomial([13.339, -2.220e3, 7.333e5, -1.065e8]),
            Temperature(25),
            13303373557.26,
        ],
        [
            13303373557.26,
            Temperature(298.15, "kelvin"),
            13303373557.26,
        ],
    ],
)


class TestKReaction(TestBasicReaction):

    # Update methods from TestBasicReaction to match the KReaction
    @pytest.fixture
    def empty_reaction(self):
        return kr.KReaction(1e8)

    @pytest.fixture
    def reaction(self, empty_reaction, reactants, products):
        for reactant in reactants:
            empty_reaction.add_as_reactant(reactant)
        for product in products:
            empty_reaction.add_as_product(product)
        empty_reaction.index = 0
        return empty_reaction

    # End overrides

    @pytest.mark.parametrize(
        "k_val, result",
        [
            (None, None),
            (4, Polynomial(np.log10(4))),
        ],
    )
    def test_init(self, k_val, result):
        r = kr.KReaction(k_val)
        assert r._k_value == result

    def test_k_value_error(self):
        r = kr.KReaction()
        with pytest.raises(ValueError):
            r.k_value()

    @param_kvalues
    def test_k_value(self, k_value, temp, result):
        r = kr.KReaction(k_value)
        assert pytest.approx(r.k_value(temp), abs=0.1) == result

    @param_kvalues
    def test_set_k_value(self, empty_reaction, k_value, temp, result):
        empty_reaction.set_k_value(k_value)
        assert pytest.approx(empty_reaction.k_value(temp), abs=0.1) == result

    def test_as_label(self, reaction):
        assert reaction.as_label() == "KR_0"

    def test_as_label_error(self, empty_reaction):
        with pytest.raises(KeyError):
            empty_reaction.as_label()

    def test_raw_k_value(self):
        assert kr.KReaction().raw_k_value() is None
        assert kr.KReaction(4).raw_k_value() == pytest.approx(4)
        coefs = [13.339, -2.220e3, 7.333e5, -1.065e8]
        assert list(kr.KReaction(coefs).raw_k_value()) == coefs
