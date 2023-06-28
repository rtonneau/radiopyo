#!/usr/bin/env python3
"""
The species submodule contains all class definitions related to species definitions. It
also provides a collection to ease their handling.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import copy
import operator as op
from contextlib import suppress
from itertools import chain
from typing import (
    Dict,
    Iterable,
    Iterator,
    Optional,
    Protocol,
)

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry import exceptions as radiopyo_exceptions

# ----------------------------- Type checking ONLY imports---------------------------- #
# ------------------------------------------------------------------------------------ #
__all__ = [
    "DynSpecies",
    "RawSpecies",
    "CstSpecies",
]

# ------------------------------------------------------------------------------------ #
#                                 PROTOCOL DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


class SimSpecies(Protocol):
    """Protocol => used for type hints"""
    label: str
    index: int

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def initial_cc(self) -> float: ...


# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class RawSpecies():
    """ Basic class for species. """
    label: str

    def __init__(self, label: str):
        self.label = label

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.label == other
        if isinstance(other, RawSpecies):
            return self.label == other.label
        raise NotImplementedError(
            f"Equality between RawSpecies and {other}")

    def __str__(self) -> str:
        return self.label


class DynSpecies(RawSpecies):
    """Class used for species whose cc can vary during the simulation."""
    cc_init: float  # Store initial value
    index: int

    def __init__(self, label: str, index: int, cc_init: float = 0) -> None:
        super().__init__(label)
        self.index = index
        self.cc_init = cc_init

    def __str__(self) -> str:
        return f"DynSpecies('{self.label}')"

    def __repr__(self) -> str:
        return self.__str__()

    def __copy__(self) -> DynSpecies:
        return DynSpecies(self.label, self.index, self.cc_init)

    def initial_cc(self) -> float:
        return self.cc_init

    def set_initial_cc(self, value: float) -> None:
        self.cc_init = float(value)


class CstSpecies(RawSpecies):
    """Class used for species whose cc cannot vary during the simulation."""
    cc_value: float
    index: int

    def __init__(self, label: str, index: int, cc: float) -> None:
        super().__init__(label)
        self.index = index
        self.cc_value = cc

    def __str__(self) -> str:
        return (f"CstSpecies('{self.label}', cc='{self.cc_value}')")

    def __repr__(self) -> str:
        return self.__str__()

    def __copy__(self) -> CstSpecies:
        return CstSpecies(self.label, self.index, self.cc_value)

    def initial_cc(self) -> float:
        return self.cc_value


class ListSpecies():
    """
    Collection (i.e. container) of species. Main goal is to centralized species
    creation so that it is possible to easily attribute unique id
    and loops over species.

    Attributes:
        dyn_sp: Dict[str, DynSpecies]
            Dictionary storing dyn species
        cst_sp: Dict[str, CstSpecies]
            Dictionary storing cst species
        next_index: int
            Index for the next Species (Cst or Dyn) creation
        sorted: bool
            Whether or not indexes of dyn species and cst species are contiguously
            separated. For simulations (n DynSpecies & m cstSpecies), indexes should be 
            like:
                DynSpecies: "O2"->0; "e_aq"->1; ... ;"H2O2"->n-1
                CstSpecies: "H2O"->n; ... ; "catalase"->n+m-1

    """
    dyn_sp: Dict[str, DynSpecies]  # Dict like store for dyn species
    cst_sp: Dict[str, CstSpecies]
    next_index: int
    sorted: bool

    def __init__(self) -> None:
        self.dyn_sp, self.cst_sp = {}, {}
        self.next_index = 0
        self.sorted = False

    def __contains__(self, key: str) -> bool:
        return key in self.dyn_sp or key in self.cst_sp

    def __iter__(self) -> Iterator[SimSpecies]:
        return chain(self.dyn_sp.values(), self.cst_sp.values())

    def __len__(self) -> int:
        """Total number of species => Cst + Dyn"""
        return self.len_cst_species() + self.len_dyn_species()

    def __getitem__(self, item: str) -> SimSpecies:
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

    def _convert_to_CstSpecies(self, item: str | SimSpecies) -> None:
        """Convert a CstSpecies into a DynSpecies.

        Args:
            item (str | SimSpecies): species to convert

        Raises:
            KeyError: is item does not exists in CstSpecies dictionary

        Returns:
            None
        """
        if isinstance(item, str):
            item = self.__getitem__(item)
        if isinstance(item, CstSpecies):
            return
        self.dyn_sp.pop(item.label)
        self._add_species(CstSpecies(
            item.label, item.index, item.initial_cc()))
        self.sorted = False

    def _add_species(self, species: DynSpecies | CstSpecies) -> None:
        """
        Generic method to add a new DynSpecies or CstSpecies. 

        Args:
            species: DynSpecies | CstSpecies

        Raises:
            TypeError
        """
        if isinstance(species, DynSpecies) and species.label not in self.dyn_sp:
            self.dyn_sp[species.label] = species
        elif isinstance(species, CstSpecies) and species.label not in self.cst_sp:
            self.cst_sp[species.label] = species
        else:
            raise TypeError(f"Unknown species type: {type(species)}")

    def get(self,
            label: str,
            default: Optional[SimSpecies] = None,
            ) -> Optional[SimSpecies]:
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
            self.dyn_sp[label] = DynSpecies(label, self.next_index)
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
            self.cst_sp[label] = CstSpecies(label, self.next_index, cc)
            self.next_index += 1

    def get_or_create_dyn_species(self,
                                  label: str,
                                  ) -> SimSpecies:
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
                                  ) -> SimSpecies:
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

    def sorted_iter_cst_species(self) -> Iterable[CstSpecies]:
        """
        Get an iterator over all the CstSpecies (ordered).

        Returns:
            Iterator(CstSpecies) 
        """
        return map(op.itemgetter(1),
                   sorted(self.cst_sp.items(), key=lambda x: x[1].index))

    def sorted_iter_dyn_species(self) -> Iterable[DynSpecies]:
        """
        Get an iterator over all the DynSpecies (ordered).

        Returns:
            Iterator(DynSpecies) 
        """
        return map(op.itemgetter(1),
                   sorted(self.dyn_sp.items(), key=lambda x: x[1].index))

    def iter_cst_species(self) -> Iterable[CstSpecies]:
        """
        Get an iterator over all the CstSpecies (not ordered).

        Returns:
            Iterator(CstSpecies) 
        """
        return iter(self.cst_sp.values())

    def iter_dyn_species(self) -> Iterable[DynSpecies]:
        """
        Get an iterator over all the DynSpecies (not ordered).

        Returns:
            Iterator(DynSpecies) 
        """
        return iter(self.dyn_sp.values())

    def sort(self, inplace: bool = True) -> Optional[ListSpecies]:
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

    def is_dyn_species(self, key: str | SimSpecies) -> bool:
        """ """
        # If CstSpecies from another sim, should still be tested for this one.
        if isinstance(key, (CstSpecies, DynSpecies)):
            key = key.label
        return key in self.dyn_sp
