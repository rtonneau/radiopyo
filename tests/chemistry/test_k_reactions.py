import pytest
import numpy as np


from radiopyo.chemistry import k_reactions as kr
from radiopyo.chemistry.species import ChemSpecies
from tests.chemistry.test_basic_chemistry import TestBasicReaction


class TestKReaction(TestBasicReaction):

    # Update methods from TestBasicReaction to match the KReaction
    @pytest.fixture
    def empty_reaction(self):
        return kr.KReaction(1e8)

    # End overrides
