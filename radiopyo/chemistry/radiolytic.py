
#!/usr/bin/env python3
"""
The radiolytic submodule contains all class and function definitions to include
reactions involving High Energy Particles with the background species => Production of
radicals and charged species. These reactions are defined using so called "G-values" ->
radiolytic yields i.e.  [#radical / 100eV / incident particle].

/!\ Here, the solvent density is assumed to be the same as Water ==> 1 kg/l

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import hashlib
import typing as tp

from scipy.constants import Avogadro, elementary_charge  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.species import ListSpecies

from .base_chemistry import BaseChemicalReaction

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.species import SimSpecies

# ------------------------------------------------------------------------------------ #

__all__ = ["ge_to_kr",
           "RadiolyticReaction",
           "umol_per_joule_to_ge",
           "RadiolyticReactionDict",
           ]

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #
# Ge units => radical / 100eV / incident particle
# Kr units => mol/l/Gy

# Define some functions for unit conversion:


def umol_per_joule_to_ge(value: float) -> float:
    return value * 1e-6 * Avogadro / elementary_charge * 100


def ge_to_kr(ge: float) -> float:
    d: float = 1.0  # Solvent density [kg/l]
    return ge * d / elementary_charge / 100.0 / Avogadro


def kr_to_ge(kr: float) -> float:
    d: float = 1.0  # Solvent density [kg/l]
    return kr / d * elementary_charge * 100.0 * Avogadro

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class RadiolyticReaction(BaseChemicalReaction):
    """
    A radiolytic reaction is constructed from the BaseChemicalReaction. Here no reactant
    should be involved. Moreover, the production rate only depends on the nature and
    energy of incident radiation ==> LET and whether or not beam is ON.
    ==> Only the list of products and the reaction cst (G value) are then relevant.

    Attributes:
        _reaction_cst: float
            The G value of the reaction in [mol/l/Gy]
    """
    _reaction_cst: float

    def __init__(self, reaction_cst: float):
        super().__init__()
        self._reaction_cst = reaction_cst  # [mol/l/Gy]

    def make_hash(self, method: str = "sha256") -> str:
        """
        The G value is not included in the hash. The aim is to have a single hash for
        two radiolytic reactions of the same species even for different Gvalue. The LET
        or actual Gvalue should be given sideway.
        (Therefore it limits the number of entries in a database)
        """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(self).encode())
        return hasher.hexdigest()

    def species(self) -> SimSpecies:
        return self.products[self.label()]

    def label(self) -> str:
        return list(self.products.keys())[0]

    def kr(self) -> float:
        return self._reaction_cst

    def ge(self) -> float:
        return kr_to_ge(self._reaction_cst)

    def update_from_ge(self, value: float) -> None:
        self._reaction_cst = ge_to_kr(value)

    @classmethod
    def from_dict(cls,
                  product: str,
                  ge: float,
                  list_sp: tp.Optional[ListSpecies] = None,
                  ) -> RadiolyticReaction:
        """Create a RadiolyticReaction from a dictionary (via destructuring)

        Args:
            product: str
                species to be produced
            ge: float
                G Value in [radical / 100eV / incident particle]
            list_sp: ListSpecies
                Collection of species already involved in the simulation. This
                collection will grow if some products do not exist yet.
        Returns:
            RadiolyticReaction
        """
        new = cls(ge_to_kr(ge))
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()
        new.add_as_product(list_sp.get_or_create_dyn_species(product.strip()))
        return new

    # ------------------------------------------------------------------------ #

    def add_as_reactant(self, sp: SimSpecies, multiple: int = 1) -> None:
        """No reactant should be involved"""
        NotImplementedError("Not implemented for Radiolytic Reactions")

    def as_label(self) -> str:
        """
        This method is used for the data mining. It should be implemented by child
        classes in order to easily define short labels for each reactions based on the
        reaction type and the index of reaction. For instance, it can be used as column 
        name in a pandas DataFrame.

        Returns:
            str

        Raises:
            KeyError
                If index is None.
        """
        if self.index is not None:
            return f"G_{self.species().label}"
        raise KeyError("Unknown reaction index. Cannot make a label.")

    def __str__(self) -> str:
        return (f"RadiolyticReaction('{self.label()}')")

    def __repr__(self) -> str:
        return self.__str__()


RadiolyticReactionDict = tp.Dict[str, float]
