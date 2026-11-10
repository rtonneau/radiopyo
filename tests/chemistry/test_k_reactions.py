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


# ------------------------------- Arrhenius --------------------------------- #
# k = A exp(-Ea / (R T)), Ea in kJ/mol, T in K, R = 8.314462618 J/mol/K.
# Hand computation: A = 1e10, Ea = 20 kJ/mol, T = 298.15 K
#   Ea / (R T) = 20000 / 2478.957 = 8.06791  ->  k = 1e10 * 3.1344e-4 = 3.1344e6
@pytest.mark.parametrize(
    "temp, expected",
    [
        (None, 3.134379609536337e6),  # default temperature: 25 degC
        (25, 3.134379609536337e6),
        (Temperature(298.15, "kelvin"), 3.134379609536337e6),
        (37, 4.282711680029348e6),  # 310.15 K
    ],
)
def test_arrhenius_k_value_kj_per_mol(temp, expected):
    arr = kr.Arrhenius(A=1e10, Ea=20.0)
    assert arr.k_value(temp) == pytest.approx(expected, rel=1e-9)


def test_arrhenius_zero_activation_energy_is_a_constant():
    assert kr.Arrhenius(A=2e9, Ea=0.0).k_value(60) == pytest.approx(2e9)


def test_arrhenius_from_generic_and_to_dict_roundtrip():
    from radiopyo.parser import RPARSER

    gen = RPARSER.parse_reaction("H_r + OH_r -> H2O")
    gen["constants"] = {"A": 1e10, "Ea": 20.0}
    arr = kr.Arrhenius.from_generic(gen)
    assert arr.k_value(25) == pytest.approx(3.134379609536337e6, rel=1e-9)
    assert arr.to_dict()["A"] == 1e10
    assert arr.to_dict()["Ea"] == 20.0


def test_arrhenius_in_a_sim_env_follows_the_simulation_temperature():
    from radiopyo.simulation.environment_build import arrhenius_from_dict

    (arr,), _ = arrhenius_from_dict(
        [{"reaction": "H_r + OH_r -> H2O", "A": 1e10, "Ea": 20.0}]
    )
    assert arr.k_value(37) > arr.k_value(25)
    assert arr.k_value(37) == pytest.approx(4.282711680029348e6, rel=1e-9)
