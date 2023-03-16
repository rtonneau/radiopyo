# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterator,
    Optional,
)

if TYPE_CHECKING:
    from radiopyo.chemistry.acid_base import ABPartition

__all__ = [
    "MappedSpecies",
    "SnapSpecies",
]


@dataclass(order=False)
class SnapSpecies:
    label: str
    cc: float = 0
    index: Optional[int] = None
    ab_derive: float = 1

    def update_cc(self, value: float) -> None:
        self.cc = value

    def update_acid_from_partition(self, partition: ABPartition) -> None:
        self.cc = partition.acid
        self.ab_derive = partition.derive_acid

    def update_base_from_partition(self, partition: ABPartition) -> None:
        self.cc = partition.base
        self.ab_derive = partition.derive_base


class MappedSpecies():
    species: Dict[str, SnapSpecies]
    number_tracked_species: int

    def __init__(self, number_tracked_species: int = 0) -> None:
        self.species = {}
        self.number_tracked_species = number_tracked_species

    def __getitem__(self, key: str) -> SnapSpecies:
        return self.species[key]

    def __iter__(self) -> Iterator[SnapSpecies]:
        return iter(self.species.values())

    def __contains__(self, key: str) -> bool:
        return key in self.species

    def __str__(self) -> str:
        out = ["MappedSpecies< "]
        out += [f"{str(elt)}, " for elt in self.species.values()]
        out.append(" >")
        return "".join(out)

    def __repr__(self) -> str:
        return self.__str__()

    def add_untracked(self, sp: SnapSpecies) -> None:
        self.species[sp.label] = sp

    def add_tracked(self, sp: SnapSpecies) -> None:
        self.species[sp.label] = sp
        self.number_tracked_species += 1

    def clone(self) -> MappedSpecies:
        out = MappedSpecies(self.number_tracked_species)
        #out.species = deepcopy(self.species)
        out.species = self.species.copy()
        return out
