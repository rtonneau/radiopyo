import pytest

from random import shuffle

from radiopyo.chemistry.species_collection import ListSpecies
from radiopyo.chemistry.species import BasicSpecies, SimSpecies


@pytest.fixture
def ordered_list_species(dyn_species, cst_species):
    """
    Provides an ordered version of the ListSpecies.
    """
    list_species = ListSpecies()
    for sp in dyn_species:
        list_species.add_dyn_species(sp.raw_species)
    for sp in cst_species:
        list_species.add_cst_species(sp.raw_species, 0)
    list_species.sorted = True
    return list_species


class TestListSpecies_empty:
    def test_init(self, empty_list_species):
        assert len(empty_list_species.dyn_sp) == 0
        assert len(empty_list_species.cst_sp) == 0
        assert empty_list_species.next_index == 0
        assert empty_list_species.sorted is False

    def test_sort_empty(self, empty_list_species):
        out = empty_list_species.sort(inplace=False)
        assert len(out) == 0
        assert out.next_index == 0
        assert out.sorted is True

    def test_set_order_empty(self, empty_list_species):
        empty_list_species.set_order([], [])
        assert empty_list_species.next_index == 0
        assert empty_list_species.sorted is True

    def test_add_dyn_species(self, empty_list_species):
        # From empty list
        empty_list_species.add_dyn_species(BasicSpecies("H2O2"))
        assert len(empty_list_species.dyn_sp) == 1
        assert empty_list_species.next_index == 1
        assert len(empty_list_species.cst_sp) == 0
        assert empty_list_species.sorted is False

    def test_add_dyn_species_twice(self, empty_list_species):
        # From empty list
        empty_list_species.add_dyn_species(BasicSpecies("H2O2"))
        # From non empty list
        with pytest.raises(KeyError):
            empty_list_species.add_dyn_species(BasicSpecies("H2O2"))

    def test_add_cst_species(self, empty_list_species):
        empty_list_species.add_cst_species(BasicSpecies("H2O"), 0)
        assert len(empty_list_species.cst_sp) == 1
        assert empty_list_species.next_index == 1
        assert len(empty_list_species.dyn_sp) == 0
        assert empty_list_species.sorted is False

    def test_add_dyn_species_existing_cst_species(self, empty_list_species):
        empty_list_species.add_cst_species(BasicSpecies("H2O"), 0)
        with pytest.raises(KeyError):
            empty_list_species.add_dyn_species(BasicSpecies("H2O"))

    def test_add_cst_species_existing_dyn_species(self, empty_list_species):
        empty_list_species.add_dyn_species(BasicSpecies("H2O"))
        with pytest.raises(KeyError):
            empty_list_species.add_cst_species(BasicSpecies("H2O"), 0)


class TestListSpecies:

    def test_contains(
        self, list_species, dyn_species, cst_species, non_existing_species
    ):
        # test for dynamical species
        for sp in dyn_species:
            assert sp in list_species
            assert sp.label in list_species
        # test for constant species
        for sp in cst_species:
            assert sp in list_species
            assert sp.label in list_species
        # test for non existing species
        assert non_existing_species not in list_species

    def test_iter(self, list_species):
        for elt in list_species:
            assert isinstance(elt, SimSpecies)
            assert elt in list_species

    def test_len(self, list_species, dyn_species, cst_species):
        assert len(list_species) == len(dyn_species) + len(cst_species)

    def test_getitem_with_existing_dyn(self, list_species, dyn_species):
        # test for dynamical species
        assert isinstance(list_species[dyn_species[0].label], SimSpecies)

    def test_getitem_with_existing_cst(self, list_species, cst_species):
        # test for constant species
        assert isinstance(list_species[cst_species[0].label], SimSpecies)

    def test_getitem_with_non_existing_species(
        self, list_species, non_existing_species
    ):
        # test for non existing species
        with pytest.raises(KeyError):
            list_species[non_existing_species.label]

    def test_len_cst_species(self, list_species, cst_species):
        assert list_species.len_cst_species() == len(cst_species)

    def test_len_dyn_species(self, list_species, dyn_species):
        assert list_species.len_dyn_species() == len(dyn_species)

    def test_get_with_existing_dyn(self, list_species, dyn_species):
        assert list_species.get(dyn_species[0].label).label == dyn_species[0].label

    def test_get_with_existing_cst(self, list_species, cst_species):
        assert list_species.get(cst_species[0].label).label == cst_species[0].label

    def test_get_non_existing_with_no_default(self, list_species, non_existing_species):
        assert list_species.get(non_existing_species.label) is None

    def test_get_non_existing_with_default_str(
        self, list_species, non_existing_species
    ):
        # The default is returned untouched (here, the str itself)
        assert (
            list_species.get(non_existing_species.label, non_existing_species.label)
            == non_existing_species.label
        )

    def test_get_non_existing_with_default_SimSpecies(
        self, list_species, non_existing_species
    ):
        assert (
            list_species.get(non_existing_species.label, non_existing_species)
            is non_existing_species
        )

    def test_remove_species_with_existing_dyn(self, list_species, dyn_species):
        assert dyn_species[0].label in list_species, "Species not in initial list"
        list_species.remove_species(dyn_species[0].label)
        assert dyn_species[0].label not in list_species

    def test_remove_species_with_existing_cst(self, list_species, cst_species):
        assert cst_species[0].label in list_species, "Species not in initial list"
        list_species.remove_species(cst_species[0].label)
        assert cst_species[0].label not in list_species

    def test_remove_species_with_non_existing_species(
        self, list_species, non_existing_species
    ):
        assert non_existing_species.label not in list_species, "Species in initial list"
        list_species.remove_species(non_existing_species.label)
        assert non_existing_species.label not in list_species

    def test_get_or_create_dyn_species_with_existing_dyn(
        self, list_species, dyn_species
    ):
        assert (
            list_species.get_or_create_dyn_species(dyn_species[0].raw_species).label
            == dyn_species[0].label
        )

    def test_get_or_create_dyn_species_with_existing_cst(
        self, list_species, cst_species
    ):
        assert (
            list_species.get_or_create_dyn_species(cst_species[0].raw_species).label
            == cst_species[0].label
        )

    def test_get_or_create_dyn_species_with_non_existing_species(
        self, list_species, non_existing_species
    ):
        assert (
            list_species.get_or_create_dyn_species(non_existing_species.raw_species).label
            == non_existing_species.label
        )

    def test_get_or_create_cst_species_with_existing_dyn(
        self, list_species, dyn_species
    ):
        assert (
            list_species.get_or_create_cst_species(dyn_species[0].raw_species).label
            == dyn_species[0].label
        )

    def test_get_or_create_cst_species_with_existing_cst(
        self, list_species, cst_species
    ):
        assert (
            list_species.get_or_create_cst_species(cst_species[0].raw_species).label
            == cst_species[0].label
        )

    def test_get_or_create_cst_species_with_non_existing_species(
        self, list_species, non_existing_species
    ):
        assert (
            list_species.get_or_create_cst_species(non_existing_species.raw_species).label
            == non_existing_species.label
        )

    def test_sorted_iter_cst_species(self, ordered_list_species):
        _prev = len(ordered_list_species.dyn_sp) - 1
        for elt in ordered_list_species.sorted_iter_cst_species():
            assert elt.index > _prev
            _prev = elt.index

    def test_sorted_iter_dyn_species(self, ordered_list_species):
        _prev = -1
        for elt in ordered_list_species.sorted_iter_dyn_species():
            assert elt.index > _prev
            _prev = elt.index

    def test_iter_cst_species(self, list_species):
        for elt in list_species.iter_cst_species():
            assert isinstance(elt, SimSpecies)
            assert elt.label in list_species.cst_sp

    def test_iter_dyn_species(self, list_species):
        for elt in list_species.iter_dyn_species():
            assert isinstance(elt, SimSpecies)
            assert elt.label in list_species.dyn_sp

    def test_is_dyn_species_with_existing_dyn(self, list_species, dyn_species):
        assert list_species.is_dyn_species(dyn_species[0].label)
        assert list_species.is_dyn_species(dyn_species[0])

    def test_is_dyn_species_with_existing_cst(self, list_species, cst_species):
        assert not list_species.is_dyn_species(cst_species[0].label)

    def test_is_dyn_species_with_non_existing_species(
        self, list_species, non_existing_species
    ):
        assert not list_species.is_dyn_species(non_existing_species.label)

    def test_is_cst_species_with_existing_cst(self, list_species, cst_species):
        assert list_species.is_cst_species(cst_species[0].label)
        assert list_species.is_cst_species(cst_species[0])

    def test_is_cst_species_with_existing_dyn(self, list_species, dyn_species):
        assert not list_species.is_cst_species(dyn_species[0].label)

    def test_is_cst_species_with_non_existing_species(
        self, list_species, non_existing_species
    ):
        assert not list_species.is_cst_species(non_existing_species.label)

    def test_set_order(self, list_species, dyn_species, cst_species):
        shuffle(dyn_species)
        shuffle(cst_species)
        list_species.set_order(
            [elt.label for elt in dyn_species],
            [elt.label for elt in cst_species],
        )
        for from_ref, elt in zip(dyn_species, list_species.sorted_iter_dyn_species()):
            assert elt.label == from_ref.label
        for from_ref, elt in zip(cst_species, list_species.sorted_iter_cst_species()):
            assert elt.label == from_ref.label
        assert list_species.sorted is True

    def test_set_order_wrong_arguments(self, list_species, dyn_species, cst_species):
        with pytest.raises(ValueError):
            list_species.set_order([elt.label for elt in dyn_species], [])
        with pytest.raises(ValueError):
            list_species.set_order([], [elt.label for elt in cst_species])

    def test_copy_with_no_change(self, list_species):
        copy = list_species.copy()
        assert copy is not list_species
        # SimSpecies are re-created by copy(): same labels, distinct objects
        assert list(copy.dyn_sp) == list(list_species.dyn_sp)
        assert list(copy.cst_sp) == list(list_species.cst_sp)
        for label, sim_sp in copy.dyn_sp.items():
            assert sim_sp is not list_species.dyn_sp[label]
            assert sim_sp.raw_species is list_species.dyn_sp[label].raw_species
        assert copy.sorted == list_species.sorted
        assert copy.next_index == list_species.next_index

    def test_copy_keeps_indices(self, list_species):
        list_species.sort(inplace=True)
        copy = list_species.copy()
        for sim_sp in list_species:
            assert copy[sim_sp.label].index == sim_sp.index
            assert copy[sim_sp.label].initial_cc() == sim_sp.initial_cc()

    def test_copy_with_changes(self, list_species, non_existing_species):
        copy = list_species.copy()
        assert copy is not list_species
        copy.add_dyn_species(non_existing_species.raw_species)
        assert list(copy.dyn_sp) != list(list_species.dyn_sp)
        assert list(copy.cst_sp) == list(list_species.cst_sp)
        assert copy.sorted == list_species.sorted
        assert copy.next_index != list_species.next_index

    def test_sort_with_not_inplace_no_starter(self, list_species):
        ls = list_species.sort(inplace=False, starter=None)
        assert id(ls) != id(list_species)
        for sp in ls.iter_dyn_species():
            assert sp.index < ls.len_dyn_species()
            assert ls.dyn_sp[sp.label] is not list_species.dyn_sp[sp.label]
        for sp in ls.iter_cst_species():
            assert sp.index >= ls.len_dyn_species()
            assert ls.cst_sp[sp.label] is not list_species.cst_sp[sp.label]

    def test_sort_with_inplace_no_starter(self, list_species):
        list_species.sort(inplace=True, starter=None)
        for sp in list_species.iter_dyn_species():
            assert sp.index < list_species.len_dyn_species()
        for sp in list_species.iter_cst_species():
            assert sp.index >= list_species.len_dyn_species()

    def test_sort_with_not_inplace_with_starter(self, list_species):
        dyn_sp = [sp.label for sp in list_species.sorted_iter_dyn_species()]
        ref = dyn_sp.copy()
        while dyn_sp == ref:
            shuffle(dyn_sp)

        ls = list_species.sort(inplace=False, starter=dyn_sp[:2])
        assert ls.dyn_sp[dyn_sp[0]].index == 0
        assert ls.dyn_sp[dyn_sp[1]].index == 1

    def test_sort_with_not_inplace_and_starter_error(self, list_species):
        dyn_sp = [sp.label for sp in list_species.sorted_iter_dyn_species()]
        # Try to set a ConstantSpecies as starter
        dyn_sp[0] = next(iter(list_species.cst_sp))
        with pytest.raises(TypeError):
            list_species.sort(inplace=True, starter=dyn_sp[:2])
