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

import numpy as np

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ------------------------------------------------------------------------------------ #
__all__ = [
    "ChemSpecies",
    "RawSpecies",
    "H2O",
]


# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class RawSpecies:
    """Basic class for species."""

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
        raise NotImplementedError(f"Equality between RawSpecies and {other}")

    def __str__(self) -> str:
        return self.label


class ChemSpecies(RawSpecies):
    """Class used for species with cc (cst or dyn)."""

    _DEFAULT_T: float = 25

    cc_init: np.polynomial.Polynomial  # Store initial value as [mol/l]
    index: int

    def __init__(
        self,
        label: str,
        index: int = -1,
        cc_init: tp.Optional[float | tp.List[float]] = None,
    ) -> None:
        super().__init__(label)
        self.index = index
        cc_init = cc_init if cc_init is not None else -1
        self.cc_init = np.polynomial.Polynomial(cc_init)

    def __repr__(self) -> str:
        return f"ChemSpecies('{self.label}', {self.index}, {self.cc_init})"

    def __copy__(self) -> ChemSpecies:
        return ChemSpecies(self.label, self.index, self.cc_init.coef)

    def initial_cc(self, T: tp.Optional[float] = None) -> float:
        T = T if T is not None else self._DEFAULT_T  # Default value
        return self.cc_init(T) if self.cc_init(T) > 0 else 0

    def set_initial_cc(self, value: float | tp.List[float]) -> None:
        self.cc_init = np.polynomial.Polynomial(value)


# Define H2O here
H2O = ChemSpecies("H2O", cc_init=[55.50, 6.075e-3, 4.110e-4, 1.496e-6, 2.619e-9])
