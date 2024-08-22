import pytest

from radiopyo.chemistry import radiolytic
from radiopyo.chemistry.species import ChemSpecies
from tests.chemistry.test_basic_chemistry import TestBasicReaction

@pytest.mark.parametrize("ge, kr", [(2.8, 2.9020e-7), (0.62, 6.4258e-8), (0.47, 4.8712e-8), (0.73,7.5659e-8)])
class TestConversion:
    def test_ge_to_kr(self, ge, kr):
        assert radiolytic.ge_to_kr(ge) == pytest.approx(kr, rel=1e-3)

    def test_kr_to_ge(self, ge, kr):
        assert radiolytic.kr_to_ge(kr) == pytest.approx(ge, rel=1e-3)


class TestRadiolyticReaction(TestBasicReaction):
    @pytest.fixture
    def products(self):
        return [ChemSpecies("e_aq")]

    @pytest.fixture
    def reactants(self):
        return []
    

    @pytest.fixture
    def empty_reaction(self):
        return radiolytic.RadiolyticReaction(radiolytic.ge_to_kr(1.5))

    @pytest.fixture
    def reaction(self, empty_reaction, products):
        for product in products:
            empty_reaction.add_as_product(product)
        return empty_reaction
    
    # Update methods from TestBasicReaction to match the RadiolyticReaction
    def test_add_as_reactant(self, empty_reaction, fake_species):
        print("coucou", fake_species, type(empty_reaction))
        with pytest.raises(NotImplementedError):
            empty_reaction.add_as_reactant(fake_species)

    def test_reactants_as_string(self, reaction):
        assert reaction.reactants_as_string() == ""

    def test_tot_number_reactants(self, reaction):
        assert reaction.tot_number_reactants() == 0

    def test_get_stoi_of_reactant(self, reaction, reactants, fake_species):
        with pytest.raises(KeyError): 
            reaction.get_stoi_of(fake_species)

    def test_from_fict(self, products):
        reaction = radiolytic.RadiolyticReaction.from_dict(products[0].label, 1.5)
        assert reaction._reaction_cst == radiolytic.ge_to_kr(1.5)
        assert len(reaction.products) == 1
        assert products[0].label in reaction.products
        assert reaction.stoi_products[products[0].label] == 1
        assert len(reaction.reactants) == 0

    def test_kr(self, reaction): 
        assert reaction.kr() == pytest.approx(radiolytic.ge_to_kr(1.5), rel=1e-4)
    
    def test_ge(self, reaction):
        assert reaction.ge() == pytest.approx(1.5, rel=1e-4)

    def test_label(self, reaction, products):
        assert reaction.label() == products[0].label
    
    def test_species(self, reaction, products):
        assert reaction.species() == products[0]

    def test_as_label(self, reaction, products):
        assert reaction.as_label() == f"G_{products[0].label}"

    @pytest.mark.parametrize("dr, t0, res", [(0.1, 1e-7, 1.55e-15), (1, 5e-6, 7.77e-13), (100, 1e-6, 1.55e-11)])
    def test_compute_initial_cc_before(self, reaction, dr, t0, res):
        assert reaction.compute_initial_cc_before(dr, t0) == pytest.approx(res, rel=1e-4)

    @pytest.mark.parametrize("ge, kr", [(2.8, 2.9020e-7), (0.62, 6.4258e-8), (0.47, 4.8712e-8), (0.73,7.5659e-8)])
    def test_update_from_ge(self, reaction, ge, kr):
        reaction.update_from_ge(ge)
        assert reaction.kr() == pytest.approx(kr, rel=1e-4)



    
