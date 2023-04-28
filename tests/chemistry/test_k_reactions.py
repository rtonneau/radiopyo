import pytest

from radiopyo.chemistry.k_reactions import KReaction


@pytest.fixture
def k_reaction_dict():
    return {'reactants': ['e_aq', 'e_aq'],
            'products': ["H2", "OH_minus", "OH_minus"],
            'k_value': 1.1e10
            }


def test_creation_from_dict(k_reaction_dict):
    reaction = KReaction.from_dict(**k_reaction_dict)

    # Test reactants
    assert len(reaction.reactants) == 1
    assert reaction.reactants[0] == 'e_aq'
    assert reaction.stoi_reactants[0] == 2

    # Test products
    assert len(reaction.products) == 2
    assert reaction.products[0] == 'H2'
    assert reaction.stoi_products[0] == 1
    assert reaction.stoi_products[1] == 2

    # Test k value
    assert 1.1e10 == pytest.approx(1.1e10)
