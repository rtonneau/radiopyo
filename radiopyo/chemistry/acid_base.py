#!/usr/bin/env python3
"""
The acid_base submodule contains all class definitions related to acid-base equilibrium.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import hashlib
import typing as tp

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.species import ListSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.species import SimSpecies

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


class ABCouple():
    """
    Acid/Base couple definition. An A/B couple is composed of 2 chemical species, one
    base and one acid. The relative concentrations of both species is defined by the pKa
    constant of the A/B reaction. The relative concentrations can be computed based on
    the total concentration of both species and the concentration of H_plus ion.

    Attributes:
        acid: CstSpecies | DynSpecies
        base: CstSpecies | DynSpecies
        pKa: float
        ka: float (computed from pKa)

    Usage:
        >>> from radiopyo.species import ListSpecies
        >>> list_of_species = ListSpecies()
        >>> acid = list_of_species.get_or_create_dyn_species("OH_r")
        >>> base = list_of_species.get_or_create_dyn_species("O_r_minus")
        >>> ab = ABCouple(acid, base, pKa=11.9)
        >>> partition = ab.compute_partition(cc_tot, cc_H_plus)

    """
    acid: SimSpecies
    base: SimSpecies
    pKa: float

    def __init__(self, acid: SimSpecies, base: SimSpecies, pKa: float):
        self.pKa = pKa
        self.acid = acid
        self.base = base

    @property
    def ka(self) -> float:
        return 10**(-self.pKa)

    def make_hash(self, method: str = "sha256") -> str:
        """
        """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(self).encode())
        return hasher.hexdigest()

    @classmethod
    def from_dict(cls,
                  acid: str,
                  base: str,
                  pKa: float,
                  list_sp: tp.Optional[ListSpecies] = None,
                  ) -> ABCouple:
        """Class method to conveniently create a new Acid/Base couple from dictionary
        destructuring (from input file parsing). The method also require an argument 
        "list_sp" which keeps track of the currently existing species so that new
        species can be created if necessary (assuming dynamic species).

        Args:
            acid: str
                acid species name
            base: str
                base species name
            pKa: float
                pKa of the reaction
            list_sp: ListSpecies
                Collection of species already involved in the simulation. This
                collection will grow if acid and/or base do not exist yet.

        Returns:
            ABCouple
        """
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()
        return cls(acid=list_sp.get_or_create_dyn_species(acid.strip()),
                   base=list_sp.get_or_create_dyn_species(base.strip()),
                   pKa=pKa,
                   )

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
        return f"{self.acid}/{self.base})"

    def __repr__(self) -> str:
        return f"ABCouple('{self.__str__()}', pKa={round(self.pKa, 5)})"


class ABDict(tp.TypedDict):
    """
    Typed-Dict used to parse config files like a breeze.
    """
    acid: str
    base: str
    pKa: float
