
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
)

from scipy.constants import Avogadro, elementary_charge  # type: ignore[import]

from .base_chemistry import BaseChemicalReaction

if TYPE_CHECKING:
    from radiopyo_vec.chemistry.species import ListSpecies, SimSpecies

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

    def __init__(self, reaction_cst: float):
        super().__init__()
        self._reaction_cst = reaction_cst

    def label(self) -> str:
        return list(self.products.keys())[0]

    def kr(self) -> float:
        return self._reaction_cst

    def ge(self) -> float:
        return kr_to_ge(self._reaction_cst)

    @classmethod
    def new_from_dict(cls,
                      product: str,
                      ge: float,
                      list_sp: ListSpecies,
                      ) -> RadiolyticReaction:
        """new_from_ge _summary_

        Args:
            species (str): _description_
            ge (float): _description_

        Returns:
            RadiolyticReaction: _description_
        """
        new = cls(ge_to_kr(ge))
        new.add_as_product(list_sp.get_or_create_dyn_species(product))
        return new

    # ------------------------------------------------------------------------ #

    def add_as_reactant(self, sp: SimSpecies) -> None:
        NotImplementedError("Not implemented for Radiolytic Reactions")

    def __str__(self) -> str:
        return f"Radiolytic<{self.label()}, Ge:{self.ge():.2f}>"

    def __repr__(self) -> str:
        return self.__str__()
