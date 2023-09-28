#!/usr/bin/env python3
"""
The acid_base submodule contains all class definitions related to acid-base equilibrium.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.species import ChemSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

__all__ = [
    "ABPartition",
    "ABCouple",
    "ABDict",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class ABPartition(tp.NamedTuple):
    """
    Acid-Base equilibrium results are stored in a NamedTuple. It is computed
    from the total concentration of the Acid-Base couple. The class provides several
    attributes for a more detailed usage. 

    Attributes:
        A: float
            Base concentration
        HA: float
            Acid concentration
        dA: float
            Derivative of the base concentration with respect to the total concentration
            of the acid/base couple
        dHA: float
            Derivative of the acid concentration with respect to the total concentration
            of the acid/base couple        
    """
    A: float
    HA: float
    dA: float
    dHA: float

    @property
    def acid(self) -> float: return self.HA

    @property
    def base(self) -> float: return self.A

    @property
    def derive_acid(self) -> float: return self.dHA

    @property
    def derive_base(self) -> float: return self.dA


# ---------------------------------------------------------------------------- #
#                         ACID BASE REACTION DEFINITION                        #
# ---------------------------------------------------------------------------- #

class ABCouple(ChemSpecies):
    """ Class """
    LENGTH_INDEX: int = 2

    acid: ChemSpecies
    base: ChemSpecies
    pKa: float

    def __init__(self,
                 acid: ChemSpecies | str,
                 base: ChemSpecies | str,
                 pKa: float,
                 index: int = -1):
        self.pKa = pKa
        self.acid = ChemSpecies(acid, index=index) if isinstance(acid, str) else acid
        self.base = ChemSpecies(base, index=index+1) if isinstance(base, str) else base
        super().__init__(self.__str__(),
                         index=index,
                         cc_init=0)

    def __contains__(self, key: str) -> bool:
        return self.acid == key or self.base == key

    def __getitem__(self, item: str) -> ChemSpecies:
        match item:
            case self.acid.label:
                return self.acid
            case self.base.label:
                return self.base
            case _:
                raise KeyError(f"Unknown {item} for A/B couple: {self}")

    @property
    def ka(self) -> float:
        return 10**(-self.pKa)

    def compute_acid_partition(self, cc_tot: NDArray, cc_H_plus: NDArray) -> NDArray:
        """
        """
        return cc_tot / (1.0 + self.ka / cc_H_plus)

    def compute_base_partition(self, cc_tot: NDArray, cc_H_plus: NDArray) -> NDArray:
        """
        """
        return cc_tot / (1.0 + cc_H_plus / self.ka)

    def compute_partition(self, cc_tot: float, cc_H_plus: float) -> ABPartition:
        """
        Method to compute the Acid/base partition.

        Args:
            cc_tot: float
                Total concentration of species i.e. [acid] + [base]
            cc_H_plus: float
                Concentration of species H_plus

        Returns:
            ABPartition
        """
        return ABPartition(
            cc_tot / (1.0 + cc_H_plus / self.ka),  # Base
            cc_tot / (1.0 + self.ka / cc_H_plus),  # Acid
            1.0 / (1.0 + cc_H_plus / self.ka),     # dBase / dCt
            1.0 / (1.0 + self.ka / cc_H_plus),     # dAcid / dCt
        )

    def __str__(self) -> str:
        """
        rem: pKa should not be included in order to easily find similar AB couples.
        """
        return f"{self.acid}/{self.base}"

    def __repr__(self) -> str:
        return ("ABCouple("
                f"acid='{self.acid}', "
                f"base='{self.base}', "
                f"pKa={round(self.pKa, 5)}, "
                f"index={self.index}"
                ")")


class ABDict(tp.TypedDict):
    """
    Typed-Dict used to parse config files like a breeze.
    """
    acid: str
    base: str
    pKa: float
