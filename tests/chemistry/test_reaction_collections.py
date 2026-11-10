"""
Tests for the reaction collections (KReactionCollection, BiReaction),
the shared reaction helpers and the MichaelisMenten / RadiolyticReaction
label and hash behaviour.
"""
import hashlib

import numpy as np
import pytest

from radiopyo.chemistry import k_reactions as kr
from radiopyo.chemistry.acid_base import AcidBase
from radiopyo.chemistry.e_reactions import EReaction
from radiopyo.chemistry.michaelis import MichaelisMenten
from radiopyo.chemistry.radiolytic import RadiolyticReaction
from radiopyo.chemistry.species import BasicSpecies
from radiopyo.chemistry.species_collection import ListSpecies
from radiopyo.parser import RPARSER


def _kr(reaction: str, k: float, index: int) -> kr.KReaction:
    r = kr.KReaction.from_generic(
        {**RPARSER.parse_reaction(reaction), "constants": {"k_value": k}}
    )
    r.index = index
    return r


# ------------------------- KReactionCollection ----------------------------- #
@pytest.fixture
def collection():
    coll = kr.KReactionCollection()
    coll.add_reaction(_kr("A -> B", 1.0, 2), label="c")
    coll.add_reaction(_kr("B -> C", 2.0, 0), label="a")
    coll.add_reaction(_kr("C -> D", 3.0, 1), label="b")
    return coll


def _assert_lookups(coll: kr.KReactionCollection) -> None:
    for label, reaction in coll.iter_labels_and_reactions():
        assert coll.get_reaction(label) is reaction


def test_collection_lookup_after_sort_by_label(collection):
    collection.sort_by_label()
    assert [label for label, _ in collection.iter_labels_and_reactions()] == [
        "a",
        "b",
        "c",
    ]
    _assert_lookups(collection)
    assert str(collection.get_reaction("a")) == "B -> C"


def test_collection_lookup_after_sort_by_reaction_index(collection):
    collection.sort_by_reaction_index()
    assert [r.index for r in collection] == [0, 1, 2]
    _assert_lookups(collection)
    assert str(collection.get_reaction("c")) == "A -> B"


def test_collection_lookup_after_update_label(collection):
    reaction = collection.get_reaction("a")
    collection.update_label("a", "z")
    # Reaction keeps its position
    assert [label for label, _ in collection.iter_labels_and_reactions()] == [
        "c",
        "z",
        "b",
    ]
    assert collection.get_reaction("z") is reaction
    with pytest.raises(KeyError):
        collection.get_reaction("a")
    _assert_lookups(collection)
    # Adding after a rename does not clash with the renamed position
    collection.add_reaction(_kr("D -> E", 4.0, 3), label="a")
    _assert_lookups(collection)
    assert len(collection) == 4


def test_collection_update_label_errors(collection):
    with pytest.raises(KeyError):
        collection.update_label("unknown", "x")
    with pytest.raises(KeyError):
        collection.update_label("a", "b")


# ------------------------------ BiReactions -------------------------------- #
@pytest.fixture
def acid_base():
    gen = RPARSER.parse_reaction("H2O2 -> H3O_plus + HO2_minus")
    gen["constants"]["k_forward"] = 0.0355
    gen["constants"]["k_backward"] = [13.339, -2.220e3, 7.333e5, -1.065e8]
    return AcidBase(gen, pKa=11.65, list_sp=ListSpecies())


def test_acid_base_kreactions_own_k(acid_base):
    for T in (None, 10.0, 40.0):
        assert acid_base.forward.k_value(T) == acid_base.k_forward(T)
        assert acid_base.backward.k_value(T) == acid_base.k_backward(T)
    assert acid_base.forward.k_value() == pytest.approx(0.0355)
    assert acid_base.backward.k_value() == pytest.approx(13303373557.26, abs=0.1)
    assert not hasattr(acid_base, "_k_forward")
    assert not hasattr(acid_base, "_k_backward")


def test_acid_base_to_dict(acid_base):
    out = acid_base.to_dict()
    assert out["k_forward"] == pytest.approx(0.0355)
    assert out["k_backward"] == [13.339, -2.220e3, 7.333e5, -1.065e8]


def test_acid_base_get_couple(acid_base):
    couple = acid_base.get_couple()
    assert couple.species.label == "H2O2"
    assert couple.conjugate.label == "HO2_minus"


def test_base_get_couple():
    gen = RPARSER.parse_reaction("NH3 + H2O -> NH4_plus + OH_minus")
    gen["constants"]["k_forward"] = 1.0
    gen["constants"]["k_backward"] = 2.0
    couple = AcidBase(gen, pKa=9.25, list_sp=ListSpecies()).get_couple()
    assert {couple.species.label, couple.conjugate.label} == {"NH3", "NH4_plus"}


def test_ereaction_kreactions_own_k():
    gen = RPARSER.parse_reaction("A + B -> C")
    gen["constants"]["k_forward"] = 1e8
    eq = EReaction(gen, constant=1e4, list_sp=ListSpecies())
    assert eq.forward.k_value() == pytest.approx(1e8)
    assert eq.backward.k_value() == pytest.approx(1e4)
    assert eq.k_forward() == eq.forward.k_value()
    assert eq.k_backward() == eq.backward.k_value()
    assert eq.k_equilibrium() == pytest.approx(1e4)


# --------------------------- Shared helpers -------------------------------- #
def test_from_generic_species_order():
    r = _kr("2 A + B -> C + 3 D", 4.0, 0)
    assert [sp.label for sp in r.reactants] == ["A", "B"]
    assert [sp.label for sp in r.products] == ["C", "D"]
    # n·k convention (ADR 0001) only for a single reactant
    single = _kr("2 A -> B", 4.0, 1)
    assert single.k_value() == pytest.approx(2.0)
    assert single.raw_k_value() == pytest.approx(4.0)


# --------------------------- MichaelisMenten ------------------------------- #
@pytest.fixture
def michaelis():
    gen = RPARSER.parse_reaction("2H2O2 -- catalase >> O2 + 2H2O")
    gen["constants"]["k_value"] = 132400000.0
    gen["constants"]["k_micha"] = 1.1
    return MichaelisMenten.from_generic(gen)


def test_michaelis_repr_closed(michaelis):
    assert repr(michaelis) == (
        "MichaelisMenten('2H2O2 -- catalase >> O2 + 2H2O', "
        "k_value=132400000.0, k_micha=1.1)"
    )


def test_michaelis_hash_unchanged(michaelis):
    legacy = (
        "MichaelisMenten('2H2O2 -- catalase >> O2 + 2H2O', "
        "k_value=132400000.0, k_micha=1.1"
    )
    assert michaelis.make_hash() == hashlib.sha256(legacy.encode()).hexdigest()


def test_michaelis_single_substrate():
    # The enzymatic grammar refuses 2 substrates: build the generic by hand
    gen = RPARSER.parse_reaction("H2O2 + O2 -> H2O")
    gen["enzyme"] = BasicSpecies("catalase")
    gen["constants"]["k_value"] = 1.0
    gen["constants"]["k_micha"] = 1.0
    with pytest.raises(ValueError):
        MichaelisMenten.from_generic(gen)


def test_michaelis_as_label(michaelis):
    with pytest.raises(KeyError):
        michaelis.as_label()
    michaelis.index = 3
    assert michaelis.as_label() == "EZ_3"


# -------------------------- RadiolyticReaction ----------------------------- #
def test_radiolytic_label_attribute_and_call():
    r = RadiolyticReaction.from_product(BasicSpecies("e_aq"), 2.6)
    assert r.label == "e_aq"
    assert r.label() == "e_aq"
    assert type(r.label()) is str
    assert str(r) == repr(r) == "RadiolyticReaction('e_aq')"
    assert r.as_label() == "G_e_aq"
    assert np.isclose(r.ge(), 2.6)


def test_michaelis_to_dict_roundtrip_keeps_k_value(michaelis):
    """2H2O2 -- catalase: to_dict used to write n*k (the k doubled when reloaded)."""
    d = michaelis.to_dict()
    assert d["k_value"] == pytest.approx(132400000.0)
    gen = RPARSER.parse_reaction(d["reaction"])
    gen["constants"] = {"k_value": d["k_value"], "k_micha": d["k_micha"]}
    again = MichaelisMenten.from_generic(gen)
    assert again.k_value == pytest.approx(michaelis.k_value)
    assert again.to_dict() == d
