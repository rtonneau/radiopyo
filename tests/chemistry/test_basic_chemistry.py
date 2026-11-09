from typing import List
import pytest

from radiopyo.chemistry import basic_reaction
from radiopyo.chemistry.species import BasicSpecies, SimSpecies


def make_sp(label: str) -> SimSpecies:
    return SimSpecies(BasicSpecies(label))


class DummyReaction(basic_reaction.BasicReaction):
    """
    Dummy class to test the BasicReaction class which is abstract.
    """

    def make_hash(self, method="sha256"):
        return ""

    def as_label(self):
        return ""

    def to_dict(self):
        return {}

    @classmethod
    def from_generic(cls, generic, list_sp=None):
        raise NotImplementedError


class TestBasicReaction:
    @pytest.fixture
    def fake_species(self):
        return make_sp("FakeSpecies")

    @pytest.fixture
    def products(self):
        return [make_sp("HO2"), make_sp("H2O")]

    @pytest.fixture
    def reactants(self):
        return [make_sp("OH"), make_sp("H2O2")]

    @pytest.fixture
    def empty_reaction(self):
        return DummyReaction()

    @pytest.fixture
    def reaction(self, empty_reaction, reactants, products):
        for reactant in reactants:
            empty_reaction.add_as_reactant(reactant)
        for product in products:
            empty_reaction.add_as_product(product)
        return empty_reaction

    def test_init(self, empty_reaction):
        assert isinstance(empty_reaction.reactants, basic_reaction.ReactionSpeciesCollection)
        assert len(empty_reaction.reactants) == 0
        assert isinstance(empty_reaction.products, basic_reaction.ReactionSpeciesCollection)
        assert len(empty_reaction.products) == 0
        assert empty_reaction.index is None

    def test_add_as_product(self, empty_reaction):
        products = [make_sp("HO2"), make_sp("H2O")]
        # test single entry
        empty_reaction.add_as_product(products[0])
        assert list(empty_reaction.products.species_col) == [products[0].label]
        assert empty_reaction.products[products[0]].species is products[0]
        assert empty_reaction.products.stoichiometry_of(products[0]) == 1
        # Adding the same species twice used to increment the stoichiometry
        # but now raises: stoichiometry is given explicitly via `stoi`.
        with pytest.raises(KeyError):
            empty_reaction.add_as_product(products[0])
        assert empty_reaction.products.stoichiometry_of(products[0]) == 1

        # test multiple entries
        empty_reaction.add_as_product(products[1], stoi=2)
        assert list(empty_reaction.products.species_col) == [
            products[0].label,
            products[1].label,
        ]
        assert empty_reaction.products.stoichiometry_of(products[0]) == 1
        assert empty_reaction.products.stoichiometry_of(products[1]) == 2

    def test_add_as_reactant(self, empty_reaction):
        reactants = [make_sp("OH"), make_sp("H2O2")]
        # test single entry
        empty_reaction.add_as_reactant(reactants[0])
        assert list(empty_reaction.reactants.species_col) == [reactants[0].label]
        assert empty_reaction.reactants[reactants[0]].species is reactants[0]
        assert empty_reaction.reactants.stoichiometry_of(reactants[0]) == 1
        with pytest.raises(KeyError):
            empty_reaction.add_as_reactant(reactants[0])
        assert empty_reaction.reactants.stoichiometry_of(reactants[0]) == 1

        # test multiple entries
        empty_reaction.add_as_reactant(reactants[1], stoi=2)
        assert list(empty_reaction.reactants.species_col) == [
            reactants[0].label,
            reactants[1].label,
        ]
        assert empty_reaction.reactants.stoichiometry_of(reactants[0]) == 1
        assert empty_reaction.reactants.stoichiometry_of(reactants[1]) == 2

    def test_reactants_as_string(self, reaction, reactants):
        assert reaction.reactants.as_string() == " + ".join(
            [sp.label for sp in reactants]
        )

    def test_products_as_string(self, reaction, products):
        assert reaction.products.as_string() == " + ".join(
            [sp.label for sp in products]
        )

    def test_has_reactant(self, reaction, reactants, fake_species):
        for sp in reactants:
            assert reaction.has_reactant(sp) is True
            assert reaction.has_reactant(sp) == reaction.has_reactant(sp.label)
        assert reaction.has_reactant(fake_species.label) is False

    def test_has_product(self, reaction, products, fake_species):
        for sp in products:
            assert reaction.has_product(sp) is True
            assert reaction.has_product(sp) == reaction.has_product(sp.label)
        assert reaction.has_product(fake_species.label) is False

    def test_number_of_products(self, reaction, products):
        assert reaction.number_of_products() == len(products)

    def test_number_of_reactants(self, reaction, reactants):
        assert reaction.number_of_reactants() == len(reactants)

    def test_get_index(self, empty_reaction):
        assert empty_reaction.get_index() == -1

    def test_set_index(self, reaction):
        reaction.set_index(4)
        assert reaction.get_index() == 4

    def test_tot_number_reactants(self, reaction, reactants):
        assert reaction.tot_number_reactants() == len(reactants)
        reaction.add_as_reactant(make_sp("extra"), stoi=3)
        assert reaction.tot_number_reactants() == len(reactants) + 3

    def test_tot_number_products(self, reaction, products):
        assert reaction.tot_number_products() == len(products)
        reaction.add_as_product(make_sp("extra"), stoi=3)
        assert reaction.tot_number_products() == len(products) + 3

    def test_get_stoi_of_reactant(self, reaction, reactants, fake_species):
        # test for the reactants
        assert reaction.get_stoi_of(reactants[0].label) == 1
        reaction.add_as_reactant(make_sp("extra"), stoi=2)
        assert reaction.get_stoi_of("extra") == 2
        assert reaction.get_stoi_of(reactants[0]) == 1

        # test for a fake entry
        with pytest.raises(KeyError):
            reaction.get_stoi_of(fake_species.label)

    def test_get_stoi_of_product(self, reaction, products, fake_species):
        # test for the products
        assert reaction.get_stoi_of(products[0].label) == 1
        reaction.add_as_product(make_sp("extra"), stoi=2)
        assert reaction.get_stoi_of("extra") == 2
        assert reaction.get_stoi_of(products[0]) == 1

        # test for a fake entry
        with pytest.raises(KeyError):
            reaction.get_stoi_of(fake_species.label)

    def test_contains(self, reaction, reactants, products, fake_species):
        for sp in reactants:
            assert sp.label in reaction
            assert sp.raw_species in reaction
        for sp in products:
            assert sp.label in reaction
            assert sp.raw_species in reaction
        assert fake_species.label not in reaction
        with pytest.raises(NotImplementedError):
            assert 3.2 in reaction
        # SimSpecies are not accepted by BasicReaction.__contains__ (only
        # str / BasicSpecies): pin the current behaviour.
        with pytest.raises(NotImplementedError):
            assert (reactants + products)[0] in reaction

    def test_iter_products(self, reaction):
        for product, stoi in reaction.products.iter_stoi():
            assert product.label in reaction.products
            assert reaction.products.stoichiometry_of(product.label) == stoi

    def test_iter_reactants(self, reaction):
        for reactant, stoi in reaction.reactants.iter_stoi():
            assert reactant.label in reaction.reactants
            assert reaction.reactants.stoichiometry_of(reactant.label) == stoi

    def test_iter_species(self, reaction):
        for sp in reaction.iter_species():
            assert isinstance(sp, basic_reaction.ReactionSpecies)
            assert sp.label in reaction
            assert sp.species.raw_species in reaction
