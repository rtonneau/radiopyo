#!/usr/bin/env python3
"""
The michaelis submodule contains all class definitions related to enzymatic reactions
treatment via the Michaelis-Menten kinetic.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    List,
    TypedDict,
)

# --------------------------------- Internal imports --------------------------------- #
from .base_chemistry import BaseChemicalReaction

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from .species import ListSpecies, SimSpecies

# ------------------------------------------------------------------------------------ #

__all__ = [
    "MichaelisMenten",
]

# ------------------------------------------------------------------------------------ #
#                                    CLASS DEFINITION                                  #
# ------------------------------------------------------------------------------------ #


class MichaelisMenten(BaseChemicalReaction):
    """ 
    Class implementing the Michaelis-Menten kinetic. Child class of BaseChemicalReaction
    with few additions. In the Michaelis-Menten kinetic the reaction rate is defined by:
        k_value * [enzyme] * [substrate] / (k_micha + [substrate])

    This kinetic is really convenient as it ensures a positive denominator during
    calculation (k_micha > 0).

    Attributes:
        k_value: float
            reaction constant
        k_micha: float
            Michaelis-Menten constant
        substrate: SimSpecies
            Species used as substrate
        enzyme: SimSpecies
            Species used as enzyme
    """
    k_value: float
    k_micha: float
    substrate: SimSpecies
    enzyme: SimSpecies

    def __init__(self,
                 k_value: float,
                 k_micha: float,
                 ):
        super().__init__()
        self.k_value = k_value
        self.k_micha = k_micha

    @classmethod
    def from_dict(cls,
                  enzyme: str,
                  substrate: List[str] | str,
                  products: List[str],
                  k_value: float,
                  k_micha: float,
                  list_sp: ListSpecies,
                  ) -> MichaelisMenten:
        """Create a MichaelisMenten from a dictionary (via destructuring)

        Args:
            enzyme: str
                Name of the species used as enzyme
            substrate: List[str]
                Name of the species used as substrate. For some reaction, more than one
                molecule of substrate are required, therefore a list is expected. For
                instance if 2 molecule of H2 are needed: ["H2", "H2"]
            products: List[str]
                List of products
            k_value: float
                k value of the reaction
            k_micha: float
                Michaelis-Menten constant of the reaction
            list_sp: ListSpecies
                Collection of species already involved in the simulation. This
                collection will grow if some reactants and/or products do not exist yet.

        Returns:
            MichaelisMenten

        Raises:
            TypeError:
                if more than one substrate "kind" if defined.
        """
        new = cls(k_value, k_micha)

        new.enzyme = list_sp.get_or_create_cst_species(enzyme.strip())
        if isinstance(substrate, str):
            substrate = [substrate.strip()]
        elif len(substrate) > 1 and len({elt.strip() for elt in substrate}) > 1:
            raise TypeError(
                f"Substrate in Michaelis Menten reaction"
                f" must all be the same if provide as list. Found: {substrate}"
            )

        new.substrate = list_sp.get_or_create_dyn_species(substrate[0].strip())

        for sp in substrate:
            new.add_as_reactant(list_sp.get_or_create_dyn_species(sp.strip()))

        for sp in products:
            new.add_as_product(list_sp.get_or_create_dyn_species(sp.strip()))
        return new

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
            return f"EZ_{self.index}"
        raise KeyError("Unknown reaction index. Cannot make a label.")

    def __str__(self) -> str:
        out = ""
        for idx, (elt, stoi) in enumerate(self.iter_reactants()):
            if idx > 0:
                out += " + "
            if stoi > 1:
                out += f"{stoi} {elt.label}"
            else:
                out += f"{elt.label}"
        out += " -> "
        for idx, (elt, stoi) in enumerate(self.iter_products()):
            if idx > 0:
                out += " + "
            if stoi > 1:
                out += f"{stoi} {elt.label}"
            else:
                out += f"{elt.label}"
        return out

    def __repr__(self) -> str:
        out = ["KReaction<Reactants:["]
        for (elt, stoi) in self.iter_reactants():
            out.append(f"{elt.label},"*stoi)
        out[-1] = out[-1].strip(",")
        out.append("], Products:[")
        for (elt, stoi) in self.iter_products():
            out.append(f"{elt.label},"*stoi)
        out[-1] = out[-1].strip(",")
        out.append(f"], kValue: {self.k_value:.2e}")
        out.append(f", index: {self.index}")
        out.append(">")
        return "".join(out)

class MichaelisMentenDict(TypedDict):
    enzyme: str
    substrate: List[str] | str
    products: List[str]
    k_value: float
    k_micha: float