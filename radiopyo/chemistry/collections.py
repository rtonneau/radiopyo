#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 24/08/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import copy
import operator as op
import typing as tp
from contextlib import suppress
from functools import wraps
from itertools import chain

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry import exceptions as radiopyo_exceptions
from radiopyo.chemistry.acid_base import ABCouple
from radiopyo.chemistry.species import ChemSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = []

# -------------------------------------- Logging ------------------------------------- #

# ------------------------------------------------------------------------------------ #
#                                DECORATOR DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class ListSpecies():
    """
    Collection (i.e. container) of species. Main goal is to centralized species
    creation so that it is possible to easily attribute unique id
    and loops over species.

    Attributes:
        dyn_sp: Dict[str, ChemSpecies]
            Dictionary storing dyn species
        cst_sp: Dict[str, ChemSpecies]
            Dictionary storing cst species
        next_index: int
            Index for the next Species (Cst or Dyn) creation
        sorted: bool
            Whether or not indexes of dyn species and cst species are contiguously
            separated. For simulations (n DynSpecies & m cstSpecies), indexes should be 
            like:
                DynSpecies: "O2"->0; "e_aq"->1; ... ;"H2O2"->n-1
                CstSpecies: "H2O"->n; ... ; "catalase"->n+m-1

    N.B.: For now it feels smarter to keep 'next_index' as an auto incremented variable
    instead of just computing the global length to assign new index to species. That
    way, even upon deletion, no index will be attributed twice!
    """
    dyn_sp: tp.Dict[str, ChemSpecies]  # Dict like store for Chem species
    cst_sp: tp.Dict[str, ChemSpecies]
    next_index: int
    sorted: bool

    def __init__(self) -> None:
        self.dyn_sp, self.cst_sp = {}, {}
        self.next_index = 0
        self.sorted = False

    def __contains__(self, key: str) -> bool:
        if key in self.dyn_sp or key in self.cst_sp:
            return True
        return any(isinstance(elt, ABCouple) and key in elt for elt in self)

    def __iter__(self) -> tp.Iterator[ChemSpecies]:
        return chain(self.dyn_sp.values(), self.cst_sp.values())

    def __len__(self) -> int:
        """Total number of species => Cst + Dyn"""
        return self.len_cst_species() + self.len_dyn_species()

    def __getitem__(self, item: str) -> ChemSpecies:
        """
        First look in dyn, then cst.

        Raises:
            KeyError

        Returns:
            None
        """
        with suppress(KeyError):
            return self.dyn_sp[item]
        return self.cst_sp[item]

    def len_cst_species(self) -> int:
        """
        Compute the number of CstSpecies defined.

        Returns:
            int
        """
        return len(self.cst_sp)

    def len_dyn_species(self) -> int:
        """
        Compute the number of DynSpecies defined

        Returns:
            int
        """
        return len(self.dyn_sp)

    def _convert_to_CstSpecies(self, item: str | ChemSpecies) -> None:
        """Convert a ChemSpecies into a CstSpecies.

        Args:
            item (str | ChemSpecies): species to convert

        Raises:
            KeyError: is item does not exists in CstSpecies dictionary

        Returns:
            None
        """
        if isinstance(item, str):
            item = self.__getitem__(item)
        if isinstance(item, ChemSpecies):
            return
        self.dyn_sp.pop(item.label)
        self._add_species(ChemSpecies(
            item.label, item.index, item.initial_cc()))
        self.sorted = False

    def get(self,
            label: str,
            default: tp.Optional[ChemSpecies] = None,
            ) -> tp.Optional[ChemSpecies]:
        """
        Getter method emulating the dict one.

        Args:
            label: str
                species name
            default: SimSpecies, optional:
                default value to returns if the species does not exists
                (default is None)

        """
        with suppress(KeyError):
            return self.__getitem__(label)
        return default

    def add_acid_base(self,
                      acid: str,
                      base: str,
                      pKa: float,
                      dynamical: bool = True,
                      ) -> None:
        """ 
        """
        if acid in self:
            raise KeyError(f"Acid species: {acid} already defined outside of AcidBase")
        if base in self:
            raise KeyError(f"Base species: {base} already defined outside of AcidBase")

        ab = ABCouple(acid, base, pKa, self.next_index)
        self.next_index += 1
        if dynamical:
            self.dyn_sp[ab.label] = ab
        else:
            self.cst_sp[ab.label] = ab

    def add_dyn_species(self, label: str) -> None:
        """
        Method use to add a new DynSpecies. The nex_index attributes is assigned to the
        species and increased.

        TODO:
            - Convert from CstSpecies to DynSpecies? Not sure it is a good idea though

        Args:
            label: str
                name of the new species

        Raises:
            IsConstantSpeciesError: if species exists and is defined as constant.

        Returns:
            None
        """
        if label in self.cst_sp:
            raise radiopyo_exceptions.IsConstantSpeciesError(
                f"{label} is already recorded as CstSpecies")
        if label not in self.dyn_sp:
            self.dyn_sp[label] = ChemSpecies(label, self.next_index)
            self.next_index += 1

    def add_cst_species(self, label: str, cc: float) -> None:
        """
        Method use to add a new CstSpecies. The nex_index attributes is assigned to the
        species and increased.

        Args:
            label: str
                name of the new species

        Raises:
            IsSimSpeciesError: if species exists and is defined as constant.

        Returns:
            None
        """
        if label in self.cst_sp:
            raise radiopyo_exceptions.IsSimSpeciesError(
                f"{label} is already recorded as SimSpecies")
        if label not in self.dyn_sp:
            self.cst_sp[label] = ChemSpecies(label, self.next_index, cc)
            self.next_index += 1

    def get_or_create_dyn_species(self,
                                  label: str,
                                  ) -> ChemSpecies:
        """
        Method use to retrieve a species. If it does not exist, the subsequent
        DynSpecies is created via a call to ListSpecies.add_dyn_species. If the species 
        does exist, it is return no matter what kind it is (Cst or Dyn).

        Args:
            label: str
                name of the new species

        Returns:
            DynSpecies | CstSpecies
        """
        with suppress(KeyError):
            return self.dyn_sp[label]
        with suppress(radiopyo_exceptions.IsConstantSpeciesError):
            self.add_dyn_species(label)
        return self.__getitem__(label)

    def get_or_create_cst_species(self,
                                  label: str,
                                  ) -> ChemSpecies:
        """
        Method use to retrieve a species. If it does not exist, the subsequent
        CstSpecies is created via a call to ListSpecies.add_cst_species. If the species 
        does exist, it is return no matter what kind it is (Cst or Dyn).

        Args:
            label: str
                name of the new species

        Returns:
            DynSpecies | CstSpecies
        """
        with suppress(KeyError):
            return self.cst_sp[label]
        with suppress(radiopyo_exceptions.IsSimSpeciesError):
            self.add_cst_species(label, 0.0)
        return self.__getitem__(label)

    def sorted_iter_cst_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Get an iterator over all the CstSpecies (ordered).

        Returns:
            Iterator(CstSpecies) 
        """
        return map(op.itemgetter(1),
                   sorted(self.cst_sp.items(), key=lambda x: x[1].index))

    def sorted_iter_dyn_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Get an iterator over all the DynSpecies (ordered).

        Returns:
            Iterator(DynSpecies) 
        """
        return map(op.itemgetter(1),
                   sorted(self.dyn_sp.items(), key=lambda x: x[1].index))

    def iter_cst_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Get an iterator over all the CstSpecies (not ordered).

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

    def sort(self, inplace: bool = True) -> tp.Optional[ListSpecies]:
        """
        Sort all SimSpecies according to the following: 
        For n DynSpecies & m cstSpecies, indexes should be like:
            - DynSpecies: "O2"->0; "e_aq"->1; ... ;"H2O2"->n-1
            - CstSpecies: "H2O"->n; ... ; "catalase"->n+m-1

        Args:
            inplace: bool, optional
                Replace the existing index or return a deepcopy of the ListSpecies 
                instance (default is True)

        Returns:
            None if inplace = True
            ListSpecies if inplace = False
        """
        out = copy.deepcopy(self) if not inplace else self
        for idx, key in enumerate(chain(out.dyn_sp.keys(), out.cst_sp.keys())):
            out[key].index = idx
        out.next_index = idx + 1
        out.sorted = True
        return out if not inplace else None

    def is_dyn_species(self, key: str | ChemSpecies) -> bool:
        """ """
        # If CstSpecies from another sim, should still be tested for this one.
        if isinstance(key, ChemSpecies):
            key = key.label
        return key in self.dyn_sp
