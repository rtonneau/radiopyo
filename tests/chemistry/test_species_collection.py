import pytest

from radiopyo.chemistry.species_collection import ListSpecies


class TestListSpecies:
    @pytest.fixture
    def list_species(self):
        return ListSpecies()
