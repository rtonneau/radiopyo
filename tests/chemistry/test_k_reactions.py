import pytest
import numpy as np
from numpy.polynomial import Polynomial

from radiopyo.chemistry import k_reactions as kr
from radiopyo.chemistry.species import ChemSpecies
from tests.chemistry.test_basic_chemistry import TestBasicReaction


class TestKReaction(TestBasicReaction):

    # Update methods from TestBasicReaction to match the KReaction
    @pytest.fixture
    def empty_reaction(self):
        return kr.KReaction(1e8)

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
