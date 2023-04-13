from __future__ import annotations

from contextlib import suppress
from itertools import chain
from typing import (
    Dict,
    Iterable,
    Iterator,
    Optional,
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
    def initial_cc(self) -> float: ...


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
    cc_init: float  # Store initial value
    index: int

    def __init__(self, label: str, index: int, cc_init: float = 0) -> None:
        super().__init__(label)
        self.index = index
        self.cc_init = cc_init

    def __str__(self) -> str:
        return f"SimSpecies<{self.label}, index={self.index}>"

    def __repr__(self) -> str:
        return self.__str__()

    def __copy__(self) -> DynSpecies:
        return DynSpecies(self.label, self.index, self.cc_init)

    def initial_cc(self) -> float:
        return self.cc_init

    def set_initial_cc(self, value: float) -> None:
        self.cc_init = float(value)


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

    def __copy__(self) -> CstSpecies:
        return CstSpecies(self.label, self.index, self.cc_value)

    def initial_cc(self) -> float:
        return self.cc_value


class ListSpecies():
    dyn_sp: Dict[str, DynSpecies]
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
        return self.len_cst_species() + self.len_dyn_species()

    def __getitem__(self, item: str) -> SimSpecies:
        with suppress(KeyError):
            return self.dyn_sp[item]
        return self.cst_sp[item]

    def __deepcopy__(self) -> ListSpecies:
        out = ListSpecies()
        out.next_index = self.next_index
        for label, dyn_sp in self.dyn_sp.items():
            out.dyn_sp[label] = dyn_sp.__copy__()
        for label, cst_sp in self.cst_sp.items():
            out.cst_sp[label] = cst_sp.__copy__()
        return out

    def len_cst_species(self) -> int:
        return len(self.cst_sp)

    def len_dyn_species(self) -> int:
        return len(self.dyn_sp)

    def _convert_to_CstSpecies(self, item: str | SimSpecies) -> None:
        if isinstance(item, str):
            item = self.__getitem__(item)
        if isinstance(item, CstSpecies):
            return
        self.dyn_sp.pop(item.label)
        self._add_species(CstSpecies(
            item.label, item.index, item.initial_cc()))
        self.sorted = False

    def _add_species(self, species: DynSpecies | CstSpecies) -> None:
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
        with suppress(KeyError):
            return self.__getitem__(label)
        return default

    def add_dyn_species(self, label: str) -> None:
        if label in self.cst_sp:
            raise radiopyo_exceptions.IsConstantSpeciesError(
                f"{label} is already recorded as CstSpecies")
        if label not in self.dyn_sp:
            self.dyn_sp[label] = DynSpecies(label, self.next_index)
            self.next_index += 1

    def add_cst_species(self, label: str, cc: float) -> None:
        if label in self.cst_sp:
            raise radiopyo_exceptions.IsSimSpeciesError(
                f"{label} is already recorded as SimSpecies")
        if label not in self.dyn_sp:
            self.cst_sp[label] = CstSpecies(label, self.next_index, cc)
            self.next_index += 1

    def get_or_create_dyn_species(self,
                                  label: str,
                                  ) -> SimSpecies:
        with suppress(KeyError):
            return self.dyn_sp[label]
        with suppress(radiopyo_exceptions.IsConstantSpeciesError):
            self.add_dyn_species(label)
        return self.__getitem__(label)

    def get_or_create_cst_species(self,
                                  label: str,
                                  ) -> SimSpecies:
        with suppress(KeyError):
            return self.cst_sp[label]
        with suppress(radiopyo_exceptions.IsSimSpeciesError):
            self.add_cst_species(label, 0.0)
        return self.__getitem__(label)

    def iter_cst_species(self) -> Iterable[CstSpecies]:
        return iter(self.cst_sp.values())

    def iter_dyn_species(self) -> Iterable[DynSpecies]:
        return iter(self.dyn_sp.values())

    def sort(self, inplace: bool = True) -> Optional[ListSpecies]:
        out = self.__deepcopy__() if not inplace else self
        for idx, key in enumerate(chain(out.dyn_sp.keys(), out.cst_sp.keys())):
            out[key].index = idx
        out.next_index = idx + 1
        out.sorted = True
        return out if not inplace else None
