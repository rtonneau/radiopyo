from typing import List
import pytest

from radiopyo.chemistry import basic_chemistry
from radiopyo.chemistry.species import ChemSpecies

class DummyReaction(basic_chemistry.BasicReaction):
    def make_hash(self, method="sha256"):
        return ""
    
    def as_label(self):
        return ""
    
    @classmethod
    def from_kwargs(cls, list_sp=None, **kwargs):
        return cls()
    

class TestBasicReaction():
    @pytest.fixture
    def fake_species(self):
        return ChemSpecies("FakeSpecies")
    
    @pytest.fixture
    def products(self):
        return [ChemSpecies("HO2"), ChemSpecies("H2O")]
    
    @pytest.fixture
    def reactants(self):
        return [ChemSpecies("OH"), ChemSpecies("H2O2")]
    
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
        assert isinstance(empty_reaction.reactants, dict)
        assert len(empty_reaction.reactants) == 0
        assert isinstance(empty_reaction.products, dict)
        assert len(empty_reaction.products) == 0
        assert isinstance(empty_reaction.stoi_products, dict)
        assert len(empty_reaction.stoi_products) == 0
        assert isinstance(empty_reaction.stoi_reactants, dict)
        assert len(empty_reaction.stoi_reactants) == 0
        assert empty_reaction.index is None

    def test_add_as_product(self, empty_reaction):
        products = [ChemSpecies("HO2"), ChemSpecies("H2O")]
        # test single entry
        empty_reaction.add_as_product(products[0])
        assert empty_reaction.products == {products[0].label: products[0]}
        assert empty_reaction.stoi_products == {products[0].label: 1}
        empty_reaction.add_as_product(products[0])
        assert empty_reaction.products == {products[0].label: products[0]}
        assert empty_reaction.stoi_products == {products[0].label: 2}

        #test multiple entries
        empty_reaction.add_as_product(products[1], multiple=2)
        assert empty_reaction.products == {products[0].label: products[0], products[1].label: products[1]}
        assert empty_reaction.stoi_products == {products[0].label: 2, products[1].label: 2}

    def test_add_as_reactant(self, empty_reaction):
        reactants = [ChemSpecies("OH"), ChemSpecies("H2O2")]
        # test single entry
        empty_reaction.add_as_reactant(reactants[0])
        assert empty_reaction.reactants == {reactants[0].label: reactants[0]}
        assert empty_reaction.stoi_reactants == {reactants[0].label: 1}
        empty_reaction.add_as_reactant(reactants[0])
        assert empty_reaction.reactants == {reactants[0].label: reactants[0]}
        assert empty_reaction.stoi_reactants == {reactants[0].label: 2}

        #test multiple entries
        empty_reaction.add_as_reactant(reactants[1], multiple=2)
        assert empty_reaction.reactants == {reactants[0].label: reactants[0], reactants[1].label: reactants[1]}
        assert empty_reaction.stoi_reactants == {reactants[0].label: 2, reactants[1].label: 2}

    def test_reactants_as_string(self, reaction, reactants):
        assert reaction.reactants_as_string() == " + ".join([sp.label for sp in reactants])

    def test_products_as_string(self, reaction, products):
        assert reaction.products_as_string() == " + ".join([sp.label for sp in products])

    def test_has_reactant(self, reaction, reactants, fake_species):
        for sp in reactants:
            assert reaction.has_reactant(sp) == 1
            assert reaction.has_reactant(sp) == reaction.has_reactant(sp.label)
        assert reaction.has_reactant(fake_species.label) == 0

    def test_has_product(self, reaction, products, fake_species):
        for sp in products:
            assert reaction.has_product(sp) == 1
            assert reaction.has_product(sp) == reaction.has_product(sp.label)
        assert reaction.has_product(fake_species.label) == 0

    def test_number_of_products(self, reaction, products):
        assert reaction.number_of_products() == len(products)

    def test_number_of_reactants(self, reaction, reactants):
        assert reaction.number_of_reactants() == len(reactants)

    def test_get_index(self, reaction):
        assert reaction.get_index() == -1

    def test_set_index(self, reaction):
        reaction.set_index(4)
        assert reaction.get_index() == 4

    def test_tot_number_reactants(self, reaction, reactants):
        assert reaction.tot_number_reactants() == len(reactants)
        reaction.add_as_reactant(reactants[0], multiple=3)
        assert reaction.tot_number_reactants() == len(reactants) + 3

    def test_tot_number_products(self, reaction, products):
        assert reaction.tot_number_products() == len(products)
        reaction.add_as_product(products[0], multiple=3)
        assert reaction.tot_number_products() == len(products) + 3

    def test_get_stoi_of_reactant(self, reaction, reactants, fake_species):
        # test for the reactants
        assert reaction.get_stoi_of(reactants[0].label) == 1
        reaction.add_as_reactant(reactants[0])
        assert reaction.get_stoi_of(reactants[0]) == 2
        
        #test for a fake entry
        with pytest.raises(KeyError):
            reaction.get_stoi_of(fake_species.label)

    def test_get_stoi_of_product(self, reaction, products, fake_species):
        #test for the products
        assert reaction.get_stoi_of(products[0].label) == 1
        reaction.add_as_product(products[0])
        assert reaction.get_stoi_of(products[0]) == 2

        #test for a fake entry
        with pytest.raises(KeyError):
            reaction.get_stoi_of(fake_species.label)

    def test_contains(self, reaction, reactants, products, fake_species):
        for sp in reactants:
            assert sp.label in reaction
            assert sp in reaction
        for sp in products:
            assert sp.label in reaction
            assert sp in reaction
        assert fake_species.label not in reaction
        with pytest.raises(NotImplementedError):
            assert 3.2 in reaction

    def test_iter_products(self, reaction):
        for product, stoi in reaction.iter_products():
            assert product.label in reaction.products
            assert reaction.stoi_products[product.label] == stoi

    def test_iter_reactants(self, reaction):
        for reactant, stoi in reaction.iter_reactants():
            assert reactant.label in reaction.reactants
            assert reaction.stoi_reactants[reactant.label] == stoi

    def test_iter_species(self, reaction):
        for sp in reaction.iter_species():
            assert isinstance(sp, ChemSpecies)
            assert sp in reaction
