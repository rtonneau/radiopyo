#!/usr/bin/env python3
"""
Module containing the ListSpecies class definition.

The ListSpecies class is a collection (i.e. container) of species.
Its main goal is to centralize species creation so that it is possible
to easily attribute unique id and loops over species.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 24/08/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import operator as op
import typing as tp
from contextlib import suppress
from itertools import chain

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry import exceptions as radiopyo_exceptions
from radiopyo.chemistry.species import ChemSpecies

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    pass

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List = []

# ---------------------------------- Logging -------------------------------- #

# --------------------------------------------------------------------------- #
#                            DECORATOR DEFINITIONS                            #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #


class ListSpecies:
    """
    Collection (i.e. container) of species. Main goal is to centralized
    species creation so that it is possible to easily attribute unique
    id and loop over species.

    Attributes:
        dyn_sp: Dict[str, ChemSpecies]
            Dictionary storing dyn species
        cst_sp: Dict[str, ChemSpecies]
            Dictionary storing cst species
        next_index: int
            Index for the next Species (Cst or Dyn) creation
        sorted: bool
            Whether or not indexes of dyn species and cst species are
            contiguously separated. For simulations
            (n DynSpecies & m cstSpecies), indexes should belike:
                DynSpecies: "O2"->0; "e_aq"->1; ... ;"H2O2"->n-1
                CstSpecies: "H2O"->n; ... ; "catalase"->n+m-1

    Notes:
        For now it feels smarter to keep 'next_index' as an auto
        incremented variable instead of just computing the global length
        to assign new index to species. That way, even upon deletion,
        no index will be attributed twice!
    """

    dyn_sp: tp.Dict[str, ChemSpecies]  # Dict like store for Chem species
    cst_sp: tp.Dict[str, ChemSpecies]
    next_index: int
    sorted: bool

    def __init__(self) -> None:
        self.dyn_sp, self.cst_sp = {}, {}
        self.next_index = 0
        self.sorted = False

    def __contains__(self, key: str | ChemSpecies) -> bool:
        """
        Checks if a given key is in the ListSpecies instance.

        Args:
            key (str | ChemSpecies): The key to be searched. It can be
                either a string or a ChemSpecies instance.

        Returns:
            bool: True if the key is found in the ListSpecies instance,
            False otherwise.
        """
        label = key.label if isinstance(key, (ChemSpecies,)) else key

        return label in self.dyn_sp or label in self.cst_sp

    def __iter__(self) -> tp.Iterator[ChemSpecies]:
        """
        Returns an iterator over all the species (both Cst and Dyn) in
        the ListSpecies instance.

        Returns:
            Iterator(ChemSpecies)
        """
        return chain(self.dyn_sp.values(), self.cst_sp.values())

    def __len__(self) -> int:
        """
        Returns the total number of species in the ListSpecies instance.

        Returns:
            int: The total number of species.
        """
        return self.len_cst_species() + self.len_dyn_species()

    def __getitem__(self, item: str) -> ChemSpecies:
        """
        Returns the ChemSpecies instance associated with the given item.

        Args:
            item (str): The label of the species to be retrieved.

        Returns:
            ChemSpecies: The ChemSpecies instance associated with the
                given item.

        Raises:
            KeyError: If the given item is not found in the ListSpecies
                instance.
        """
        with suppress(KeyError):
            return self.dyn_sp[item]
        with suppress(KeyError):
            return self.cst_sp[item]
        raise KeyError(item)

    def len_cst_species(self) -> int:
        """
        Returns the number of constant species in the ListSpecies
        instance.

        Returns:
            int: The number of constant species.
        """
        return len(self.cst_sp)

    def len_dyn_species(self) -> int:
        """
        Returns the number of dynamic species in the ListSpecies
        instance.

        Returns:
            int: The number of dynamic species.
        """
        return len(self.dyn_sp)

    def get(
        self,
        label: str,
        default: tp.Optional[str | ChemSpecies] = None,
    ) -> tp.Optional[ChemSpecies]:
        """
        Getter method emulating the dict one.

        Args:
            label (str): The label of the species to be retrieved.
            default (ChemSpecies, optional): The default value to return
                if the species does not exist. Defaults to None.

        Returns:
            ChemSpecies: The ChemSpecies instance associated with the
                given label, or the default value if not found.
        """
        with suppress(KeyError):
            return self.__getitem__(label)
        if isinstance(default, str):
            return ChemSpecies(default)
        return default

    def remove_species(self, label: str) -> None:
        """
        Removes a species from the ListSpecies instance by its label.

        Args:
            label (str): The label of the species to be removed.

        Returns:
            None
        """
        if label in self.dyn_sp:
            self.dyn_sp.pop(label)
            return

        if label in self.cst_sp:
            self.cst_sp.pop(label)
            return

    def add_dyn_species(
        self, label: str, cc_init: tp.Optional[tp.List[float] | float] = None
    ) -> None:
        """
        Adds a dynamic species to the ListSpecies instance.

        Args:
            label (str): The label of the dynamic species.
            cc_init (tp.Optional[tp.List[float] | float], optional):
                The initial concentration coefficient of the dynamic
                species. Defaults to None.

        Returns:
            None
        """
        if label in self:
            raise KeyError(f"{label} is already recorded")

        self.dyn_sp[label] = ChemSpecies(label, self.next_index, cc_init)
        self.next_index += 1

    def add_cst_species(self, label: str, cc: float) -> None:
        """
        Adds a constant species to the ListSpecies instance.

        Args:
            label (str): The label of the constant species.
            cc (float): The concentration coefficient of the
                constant species.

        Returns:
            None
        """
        if label in self:
            raise KeyError(f"{label} is already recorded")
        self.cst_sp[label] = ChemSpecies(label, self.next_index, cc)
        self.next_index += 1

    def get_or_create_dyn_species(
        self,
        label: str,
    ) -> ChemSpecies:
        """
        Retrieves or creates a dynamic species from the ListSpecies
        instance.

        Args:
            label (str): The label of the dynamic species.

        Returns:
            ChemSpecies: The retrieved or created dynamic species.
        """
        if (out := self.get(label)) is not None:
            return out

        with suppress(radiopyo_exceptions.IsConstantSpeciesError):
            self.add_dyn_species(label)
        return self.__getitem__(label)

    def get_or_create_cst_species(
        self,
        label: str,
    ) -> ChemSpecies:
        """
        Retrieves or creates a constant species from the ListSpecies
        instance.

        Args:
            label (str): The label of the constant species.

        Returns:
            ChemSpecies: The retrieved or created constant species.
        """
        if (out := self.get(label)) is not None:
            return out
        with suppress(radiopyo_exceptions.IsSimSpeciesError):
            self.add_cst_species(label, 0.0)
        return self.__getitem__(label)

    def sorted_iter_cst_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Get an iterator over all the CstSpecies (ordered).

        Returns:
            Iterator(ChemSpecies)
        """
        return map(
            op.itemgetter(1), sorted(self.cst_sp.items(), key=lambda x: x[1].index)
        )

    def sorted_iter_dyn_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Get an iterator over all the dynamic species (ordered).

        Returns:
            Iterator(ChemSpecies)
        """
        return map(
            op.itemgetter(1), sorted(self.dyn_sp.items(), key=lambda x: x[1].index)
        )

    def iter_cst_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Get an iterator over all the constant species (not ordered).

        Returns:
            Iterator(ChemSpecies)
        """
        return iter(self.cst_sp.values())

    def iter_dyn_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Get an iterator over all the DynSpecies (not ordered).

        Returns:
            Iterator(ChemSpecies)
        """
        return iter(self.dyn_sp.values())

    def set_order(self, dyn_sp: tp.List, cst_sp: tp.List) -> None:
        """
        Set the order of the species in the ListSpecies.

        Args:
            dyn_sp (List): A list of dynamic species labels.
            cst_sp (List): A list of constant species labels.

        Returns:
            None
        """
        if len(dyn_sp) != len(self.dyn_sp) or len(cst_sp) != len(self.cst_sp):

            raise ValueError(
                "Unable to sort species, different numbers\n"
                f"Number of dyn sp read: {len(self.dyn_sp)}\n"
                f"Number of cst sp read: {len(self.cst_sp)}\n"
                f"Length of order list: {len(dyn_sp)} - {len(cst_sp)}"
            )
        for idx, label in enumerate(dyn_sp):
            self[label].index = idx
        for idx, label in enumerate(cst_sp, start=len(dyn_sp)):
            self[label].index = idx
        self.next_index = idx + 1
        self.sorted = True

    def sort(
        self, inplace: bool = True, starter: tp.Optional[tp.List[str]] = None
    ) -> tp.Optional[ListSpecies]:
        """
        Sort the species in the ListSpecies.
        For n DynSpecies & m cstSpecies, indexes should be like:
            - DynSpecies: "O2"->0; "e_aq"->1; ... ;"H2O2"->n-1
            - CstSpecies: "H2O"->n; ... ; "catalase"->n+m-1
        Args:
            inplace (bool, optional): Whether to sort the species in
                place. Defaults to True.
            starter (List[str], optional): A list of dynamic species
                labels to start the sorting with. Defaults to None.

        Returns:
            ListSpecies: The sorted ListSpecies if inplace is False,
                otherwise None.
        """

        out = copy.deepcopy(self) if not inplace else self
        if starter is None:
            starter = []
        all_dyn = {elt.label for elt in out.iter_dyn_species()}

        for idx, label in enumerate(starter):
            if self.is_cst_species(label):
                raise TypeError(
                    f"{label} is defined as constant. Not able to use "
                    "it in an A/B equilibrium"
                )
            out[label].index = idx

        for idx, label in enumerate(all_dyn - set(starter), start=len(starter)):
            out[label].index = idx

        for idx, species in enumerate(
            out.iter_cst_species(), start=out.len_dyn_species()
        ):
            out[species.label].index = idx

        out.next_index = idx + 1
        out.sorted = True
        return out if not inplace else None

    def is_cst_species(self, key: str | ChemSpecies) -> bool:
        """
        Checks if a given key is a constant species.

        Args:
            key (str | ChemSpecies): The key to be searched. It can be
                either a string or a ChemSpecies instance.

        Returns:
            bool: True if the key is found in the constant species,
                False otherwise.
        """
        if isinstance(key, ChemSpecies):
            key = key.label
        return key in self.cst_sp

    def is_dyn_species(self, key: str | ChemSpecies) -> bool:
        """
        Checks if a given key is a dynamic species.

        Args:
            key (str | ChemSpecies): The key to be searched. It can be
                either a string or a ChemSpecies instance.

        Returns:
            bool: True if the key is found in the dynamic species,
                False otherwise.
        """
        # If CstSpecies from another sim, should still be tested for this one.
        if isinstance(key, ChemSpecies):
            key = key.label
        return key in self.dyn_sp
