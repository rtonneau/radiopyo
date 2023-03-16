from __future__ import annotations

from itertools import chain
from typing import (
    Dict,
    Iterable,
)

__all__ = [
    "SimSpecies",
    "RawSpecies",
    "CstSpecies",
]

# ------------------------------------------------------------------------------------ #
#                                 PROTOCOL DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


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

    def as_str(self) -> str:
        """"""
        return self.label


class SimSpecies(RawSpecies):
    """"""
    index: int

    def __init__(self, label: str, index: int) -> None:
        super().__init__(label)
        self.index = index

    def __str__(self) -> str:
        return f"SimSpecies<{self.label}>"

    def __repr__(self) -> str:
        return self.__str__()


class CstSpecies(SimSpecies):
    """"""
    cc_value: float

    def __init__(self, label: str, index: int, cc: float) -> None:
        super().__init__(label, index)
        self.cc_value = cc

    def __str__(self) -> str:
        return f"CstSpecies<{self.label}, {self.cc_value:.1g}>"

    def __repr__(self) -> str:
        return self.__str__()


class ListSpecies():
    sim_sp: Dict[str, SimSpecies]
    cst_sp: Dict[str, CstSpecies]
    index: int

    def __init__(self) -> None:
        self.sim_sp, self.cst_sp = {}, {}
        self.index = 0

    def __contains__(self, key: str) -> bool:
        return key in self.sim_sp or key in self.cst_sp

    def __iter__(self) -> Iterable[SimSpecies|CstSpecies]:
        return iter(chain(self.sim_sp.values(), self.cst_sp.values()))

    def get(self, label: str) -> SimSpecies | CstSpecies:
        try:
            return self.sim_sp[label]
        except KeyError:
            return self.cst_sp[label]

    def total_species(self) -> int:
        return len(self.sim_sp) + len(self.cst_sp)

    def _add_species(self, species: SimSpecies | CstSpecies) -> None:
        if isinstance(species, SimSpecies) and species not in self.sim_sp:
            self.sim_sp[species.label] = species
        elif isinstance(species, CstSpecies) and species not in self.cst_sp:
            self.cst_sp[species.label] = species

    def add_sim_species(self, label: str) -> None:
        if label in self.cst_sp:
            raise KeyError(
                f"Species {label} is already recorded as CstSpecies")
        if label not in self.sim_sp:
            self.sim_sp[label] = SimSpecies(label, self.index)
            self.index += 1

    def add_cst_species(self, label: str, cc: float) -> None:
        if label in self.cst_sp:
            raise KeyError(
                f"Species {label} is already recorded as SimSpecies")
        if label not in self.sim_sp:
            self.cst_sp[label] = CstSpecies(label, self.index, cc)
            self.index += 1

    def get_or_create_sim_species(self, label: str) -> SimSpecies:
        if label not in self.sim_sp:
            self.add_sim_species(label)
        return self.sim_sp[label]

    def get_or_create_cst_species(self, label: str) -> CstSpecies:
        if label not in self.cst_sp:
            self.add_cst_species(label, 0.0)
        return self.cst_sp[label]

    def iter_cst_species(self) -> Iterable[CstSpecies]:
        return iter(self.cst_sp.values())

    def iter_sim_species(self) -> Iterable[SimSpecies]:
        return iter(self.sim_sp.values())
