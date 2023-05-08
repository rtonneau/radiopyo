#!/usr/bin/env python3
"""
The k_reactions submodule contains all class definitions related to basic chemical
reactions defined via a reaction constant.

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
    from .species import ListSpecies

# ------------------------------------------------------------------------------------ #

__all__ = [
    "KReaction",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #
class KReaction(BaseChemicalReaction):
    """
    Almost all the logic is already included in the BaseChemicalReaction class.
    Add the reaction rate as attribute.
    For reactions where only one reactant is involved e.g.
        2 H_r -> H2
    in that case, 2k is usually given/measured instead of k (reaction constant).
    User should simply pass the 2k value at construction and then it will be accounted
    for during calculation of reaction rate.

    Attributes:
        k_value: float
            Reaction rate
    """
    k_value: float

    def __init__(self,
                 k_value: float,
                 ):
        super().__init__()
        self.k_value = k_value

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
            return f"KR_{self.index}"
        raise KeyError(
            "Unknown reaction index. Cannot make a label without index.")

    @classmethod
    def from_dict(cls,
                  reactants: List[str],
                  products: List[str],
                  k_value: float,
                  list_sp: ListSpecies,
                  ) -> KReaction:
        """Create a KReaction from a dictionary (via destructuring)

        Args:
            reactants: List[str]
                List of reactants
            products: List[str]
                List of products
            k_value: float
                k value of the reactions
            list_sp: ListSpecies
                Collection of species already involved in the simulation. This
                collection will grow if some reactants and/or products do not exist yet.

        Returns:
            KReaction
        """
        new = cls(k_value)
        for sp in reactants:
            sim_sp = list_sp.get_or_create_dyn_species(sp.strip())
            new.add_as_reactant(sim_sp)

        # Case of n*A -> B, then n*k is measured/given
        if len(new.reactants) == 1:
            # Compute k: divide by reactant's stoichio
            new.k_value /= new.stoi_reactants[sp]

        for sp in products:
            sim_sp = list_sp.get_or_create_dyn_species(sp.strip())
            new.add_as_product(sim_sp)
        return new

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

class KReactionDict(TypedDict):
    reactants: List[str]
    products: List[str]
    k_value: float