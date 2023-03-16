
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    List,
)

import numpy as np
from scipy.constants import Avogadro, elementary_charge  # type: ignore[import]

from radiopyo.chemistry.base_chemistry import BaseChemicalReaction, ReactionSpecies

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.utils.collections import MappedSpecies

__all__ = ["ge_to_kr",
           "RadiolyticReaction",
           "umol_per_joule_to_ge"
           ]

# Ge units => radical / 100eV / incident particle
# Kr units => mol/l/Gy


def umol_per_joule_to_ge(value: float) -> float:
    return value * 1e-6 * Avogadro / elementary_charge * 100


def ge_to_kr(ge: float) -> float:
    d: float = 1.0  # Solvent density [kg/l]
    return ge * d / elementary_charge / 100.0 / Avogadro


def kr_to_ge(kr: float) -> float:
    d: float = 1.0  # Solvent density [kg/l]
    return kr / d * elementary_charge * 100.0 * Avogadro


class RadiolyticReaction(BaseChemicalReaction):
    """
    """
    _reaction_cst: float

    def __init__(self, species: List[ReactionSpecies], reaction_cst: float):
        if len(species) > 1:
            raise TypeError(
                f"A radiolytic reaction can contain only one product! ({species})")
        super().__init__(species)
        self._reaction_cst = reaction_cst

    def label(self) -> str:
        return self.products[0].label

    def kr(self) -> float:
        return self._reaction_cst

    def ge(self) -> float:
        return kr_to_ge(self._reaction_cst)

    @classmethod
    def new_from_ge(cls, species: str, ge: float) -> RadiolyticReaction:
        """new_from_ge _summary_

        Args:
            species (str): _description_
            ge (float): _description_

        Returns:
            RadiolyticReaction: _description_
        """
        return cls([ReactionSpecies.new_product(species),], ge_to_kr(ge))

    # ----- Implement Protocol IsChemicalReaction from base_chemistry.py ----- #
    def compute_reaction(self,
                         dose_rate: float,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        return self._rate_vec * self.kr() * dose_rate

    def compute_jacobian(self,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        return np.zeros(species.number_tracked_species)
    # ------------------------------------------------------------------------ #

    def add_as_reactant(self, sp: str) -> None:
        NotImplementedError("Not implemented for Radiolytic Reactions")

    def __str__(self) -> str:
        return f"Radiolytic<{self.products[0].label}, Ge:{self.ge():.2f}>"

    def __repr__(self) -> str:
        return self.__str__()
