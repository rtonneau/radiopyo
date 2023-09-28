#!/usr/bin/env python3
"""
The species submodule contains all class definitions related to species definitions. It
also provides a collection to ease their handling.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ------------------------------------------------------------------------------------ #
__all__ = [
    "ChemSpecies",
    "RawSpecies",
]


# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class RawSpecies():
    """ Basic class for species. """
    label: str

    def __init__(self, label: str):
        self.label = label

    def __contains__(self, other: object) -> bool:
        return self.__eq__(other)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.label == other
        if isinstance(other, RawSpecies):
            return self.label == other.label
        raise NotImplementedError(
            f"Equality between RawSpecies and {other}")

    def __str__(self) -> str:
        return self.label


class ChemSpecies(RawSpecies):
    """Class used for species with cc (cst or dyn)."""
    LENGTH_INDEX: int = 1

    cc_init: float  # Store initial value as [mol/l]
    index: int

    def __init__(self,
                 label: str,
                 index: int = -1,
                 cc_init: float = -1) -> None:
        super().__init__(label)
        self.index = index
        self.cc_init = cc_init

    def __repr__(self) -> str:
        return f"ChemSpecies('{self.label}', {self.index}, {self.cc_init})"

    def __copy__(self) -> ChemSpecies:
        return ChemSpecies(self.label, self.index, self.cc_init)

    @property
    def cc(self) -> float:
        return self.cc_init

    def initial_cc(self) -> float:
        return self.cc_init if self.cc_init > 0 else 0

    def set_initial_cc(self, value: float) -> None:
        self.cc_init = float(value)
