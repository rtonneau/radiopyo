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
    dyn_AB: tp.Dict[str, ABCouple]
    next_index: int
    sorted: bool

    def __init__(self) -> None:
        self.dyn_sp, self.cst_sp = {}, {}
        self.dyn_AB = {}
        self.next_index = 0
        self.sorted = False

    def __contains__(self, key: str | ChemSpecies | ABCouple) -> bool:
        label = key.label if isinstance(key, (ChemSpecies, ABCouple)) else key

        if "/" in label:  # it is an A/B
            return label in self.dyn_AB

        # It is not an A/B.
        # 1) Check if def as dyn_sp or cst_sp
        if label in self.dyn_sp or label in self.cst_sp:
            return True
        # 2) if not, then check if member of an A/B couple
        return any(label in ab_couple for ab_couple in self.iter_AB())

    def __iter__(self) -> tp.Iterator[ChemSpecies]:
        return chain(self.dyn_AB.values(),
                     self.dyn_sp.values(),
                     self.cst_sp.values())

    def __len__(self) -> int:
        """Total number of species => Cst + Dyn"""
        return self.len_cst_species() +\
            self.len_dyn_species() +\
            self.len_dyn_AB()

    def __getitem__(self, item: str) -> ChemSpecies | ABCouple:
        """
        First look in dyn, then cst.

        Raises:
            KeyError

        Returns:
            ChemSpecies | ABCouple
        """
        with suppress(KeyError):
            return self.dyn_sp[item]
        with suppress(KeyError):
            return self.cst_sp[item]
        with suppress(KeyError):
            return self.dyn_AB[item]
        for ab_reaction in self.iter_AB():
            with suppress(KeyError):
                return ab_reaction[item]
        raise KeyError(item)

    def developed_len(self) -> int:
        """
        """
        return self.len_species() * ChemSpecies.LENGTH_INDEX + \
            self.len_AB() * ABCouple.LENGTH_INDEX

    def len_dyn_AB(self) -> int:
        """ 
        """
        return len(self.dyn_AB)

    def len_AB(self) -> int:
        """
        """
        return self.len_dyn_AB()

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

    def len_species(self) -> int:
        """ 
        """
        return self.len_cst_species() + self.len_dyn_species()

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
                      pKa: float
                      ) -> None:
        """
        Rem: There is no reason to have a constant_AB variable.  
        """
        match (self.get(acid, default=None), self.get(base, default=None)):
            # Neither of acid nor base are pre-existing, easy.
            case (None, None):
                ab = ABCouple(acid,
                              base,
                              pKa,
                              self.next_index,
                              )
            # Acid is existing but base is not, then check if cst. If not, delete it.
            case (_acid, None) if isinstance(_acid, ChemSpecies):
                if not self.is_dyn_species(_acid) or self.is_AB_species(_acid):
                    # If cst, then error instead
                    raise radiopyo_exceptions.ChemistryError(
                        "Unable to define an AB couple with a cst species or a species "
                        "already involved in another A/B in it"
                    )
                self.remove_species(acid)
                ab = ABCouple(acid,
                              base,
                              pKa,
                              self.next_index,
                              )
            # Base is existing but acid is not, then check if cst. If not, delete it.
            case (None, _base) if isinstance(_base, ChemSpecies):
                if not self.is_dyn_species(_base) or self.is_AB_species(_base):
                    # If cst, then error instead
                    raise radiopyo_exceptions.ChemistryError(
                        "Unable to define an AB couple with a cst species or a species "
                        "already involved in another A/B in it"
                    )
                self.remove_species(base)
                ab = ABCouple(acid,
                              base,
                              pKa,
                              self.next_index,
                              )
            # Both base and acid are already existing, then check if cst.
            # If not, delete them both.
            case (_acid, _base) if (isinstance(_acid, ChemSpecies) and
                                    isinstance(_base, ChemSpecies)):
                if not self.is_dyn_species(_acid) or self.is_AB_species(_acid):
                    # If cst, then error instead
                    raise radiopyo_exceptions.ChemistryError(
                        "Unable to define an AB couple with a cst species or a species "
                        "already involved in another A/B in it"
                    )
                if not self.is_dyn_species(_base) or self.is_AB_species(_base):
                    # If cst, then error instead
                    raise radiopyo_exceptions.ChemistryError(
                        "Unable to define an AB couple with a cst species or a species "
                        "already involved in another A/B in it"
                    )
                self.remove_species(base)
                self.remove_species(acid)
                ab = ABCouple(acid,
                              base,
                              pKa,
                              self.next_index,
                              )
            # Default case, just in case...
            case _:
                raise radiopyo_exceptions.ChemistryError("Unable to add the A/B couple:"
                                                         f" {acid}/{base} (pKa={pKa})")
        if ab.label not in self.dyn_AB:
            self.dyn_AB[ab.label] = ab
            self.next_index += ABCouple.LENGTH_INDEX  # an A/B couple needs 2 locations

    def remove_species(self, label: str) -> None:
        """ 
        """
        if "/" in label:
            raise NotImplementedError("Unable to remove and Acid/Base link")

        if label in self.dyn_sp:
            self.dyn_sp.pop(label)
            return

        if label in self.cst_sp:
            self.cst_sp.pop(label)
            return

    def _convert_to_CstSpecies(self, item: str | ChemSpecies) -> None:
        """Convert a CstSpecies into a DynSpecies.

        Args:
            item (str | SimSpecies): species to convert

        Raises:
            KeyError: is item does not exists in CstSpecies dictionary
            ChemistryError: item is an AB couple.
        Returns:
            None
        """
        if isinstance(item, str):
            item = self.__getitem__(item)
            if isinstance(item, ABCouple):
                raise radiopyo_exceptions.ChemistryError(
                    f"Conversion of A/B: {item} to a cst Species not allowed."
                )
        if self.is_dyn_species(item):
            return
        self.dyn_sp.pop(item.label)
        self.add_cst_species(item.label, item.initial_cc())
        self.sorted = False

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
        if label not in self:
            self.dyn_sp[label] = ChemSpecies(label, self.next_index)
            self.next_index += ChemSpecies.LENGTH_INDEX

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
                f"{label} is already recorded as CstSpecies")
        if label not in self:
            self.cst_sp[label] = ChemSpecies(label, self.next_index, cc)
            self.next_index += ChemSpecies.LENGTH_INDEX

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
            ChemSpecies
        """
        if (out := self.get(label)) is not None:
            return out

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
            ChemSpecies
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

    def iter_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Get an iterator over all the Species
        (not ordered but all dyn first then cst)

        Returns:
            Iterator(ChemSpecies) 
        """
        return chain(self.iter_dyn_species(), self.iter_cst_species())

    def iter_AB(self) -> tp.Iterator[ABCouple]:
        """
        Get an iterator over all the AB couples
        (not ordered but all dyn first then cst).

        Returns:
            Iterator(ABCouple) 
        """
        return iter(self.dyn_AB.values())

    def iter_sim_species(self) -> tp.Iterator[ABCouple | ChemSpecies]:
        """ 
        """
        return chain(self.iter_AB(), self.iter_dyn_species())

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
        for idx_ab, ab_couple in enumerate(out.iter_AB()):
            out[ab_couple.label].index = idx_ab*ABCouple.LENGTH_INDEX
        for idx, species in enumerate(out.iter_species(),
                                      start=self.len_AB()*ABCouple.LENGTH_INDEX):
            out[species.label].index = idx*ChemSpecies.LENGTH_INDEX

        out.next_index = idx + 1
        out.sorted = True
        return out if not inplace else None

    def is_dyn_species(self, key: str | ChemSpecies) -> bool:
        """ """
        # If CstSpecies from another sim, should still be tested for this one.
        if isinstance(key, ChemSpecies):
            key = key.label
        return key in self.dyn_sp

    def is_AB_species(self, key: str | ABCouple | ChemSpecies) -> bool:
        """ """
        if isinstance(key, (ChemSpecies, ABCouple)):
            key = key.label
        if key in self.dyn_AB:
            return True
        return any(key in reaction for reaction in self.iter_AB())
