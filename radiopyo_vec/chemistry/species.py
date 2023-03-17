from __future__ import annotations

from itertools import chain
from typing import (
    Dict,
    Iterable,
    Iterator,
    Protocol,
)

from radiopyo_vec.chemistry import exceptions as radiopyo_exceptions

__all__ = [
    "DynSpecies",
    "RawSpecies",
    "CstSpecies",
]

# ------------------------------------------------------------------------------------ #
#                                 PROTOCOL DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


class SimSpecies(Protocol):
    label: str
    index: int

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...


# ------------------------------------------------------------------------------------ #
#                                  SPECIES DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


class RawSpecies():
    """"""
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
    """"""
    index: int

    def __init__(self, label: str, index: int) -> None:
        super().__init__(label)
        self.index = index

    def __str__(self) -> str:
        return f"SimSpecies<{self.label}, index={self.index}>"

    def __repr__(self) -> str:
        return self.__str__()


class CstSpecies(RawSpecies):
    """"""
    cc_value: float
    index: int

    def __init__(self, label: str, index: int, cc: float) -> None:
        super().__init__(label)
        self.index = index
        self.cc_value = cc

    def __str__(self) -> str:
        return (f"CstSpecies<{self.label}, "
                f"cc={self.cc_value:.1g}, "
                f"index= {self.index}"
                f">")

    def __repr__(self) -> str:
        return self.__str__()


class ListSpecies():
    dyn_sp: Dict[str, DynSpecies]
    cst_sp: Dict[str, CstSpecies]
    index: int

    def __init__(self) -> None:
        self.dyn_sp, self.cst_sp = {}, {}
        self.index = 0

    def __contains__(self, key: str) -> bool:
        return key in self.dyn_sp or key in self.cst_sp

    def __iter__(self) -> Iterator[SimSpecies]:
        return chain(self.cst_sp.values(), self.dyn_sp.values())

    def __len__(self) -> int:
        return len(self.dyn_sp) + len(self.cst_sp)

    def get(self, label: str) -> SimSpecies:
        try:
            return self.dyn_sp[label]
        except KeyError:
            return self.cst_sp[label]

    def _add_species(self, species: DynSpecies | CstSpecies) -> None:
        if isinstance(species, DynSpecies) and species not in self.dyn_sp:
            self.dyn_sp[species.label] = species
        elif isinstance(species, CstSpecies) and species not in self.cst_sp:
            self.cst_sp[species.label] = species

    def add_dyn_species(self, label: str) -> None:
        if label in self.cst_sp:
            raise radiopyo_exceptions.IsConstantSpeciesError(
                f"{label} is already recorded as CstSpecies")
        if label not in self.dyn_sp:
            self.dyn_sp[label] = DynSpecies(label, self.index)
            self.index += 1

    def add_cst_species(self, label: str, cc: float) -> None:
        if label in self.cst_sp:
            raise radiopyo_exceptions.IsSimSpeciesError(
                f"{label} is already recorded as SimSpecies")
        if label not in self.dyn_sp:
            self.cst_sp[label] = CstSpecies(label, self.index, cc)
            self.index += 1

    def get_or_create_dyn_species(self,
                                  label: str,
                                  ) -> SimSpecies:
        if label not in self.dyn_sp:
            try:
                self.add_dyn_species(label)
            except radiopyo_exceptions.IsConstantSpeciesError:
                return self.cst_sp[label]
        return self.dyn_sp[label]

    def get_or_create_cst_species(self,
                                  label: str,
                                  ) -> SimSpecies:
        if label not in self.cst_sp:
            try:
                self.add_cst_species(label, 0.0)
            except radiopyo_exceptions.IsSimSpeciesError:
                return self.dyn_sp[label]
        return self.cst_sp[label]

    def iter_cst_species(self) -> Iterable[CstSpecies]:
        return iter(self.cst_sp.values())

    def iter_dyn_species(self) -> Iterable[DynSpecies]:
        return iter(self.dyn_sp.values())
