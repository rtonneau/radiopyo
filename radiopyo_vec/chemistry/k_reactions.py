from __future__ import annotations

from functools import reduce
from operator import mul
from typing import (
    TYPE_CHECKING,
    List,
)

import numpy as np

from .base_chemistry import BaseChemicalReaction

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo_vec.utils.collections import MappedSpecies

    from .species import ListSpecies

__all__ = [
    "KReaction",
]

# ---------------------------------------------------------------------------- #
#                             K REACTION DEFINITION                            #
# ---------------------------------------------------------------------------- #


class KReaction(BaseChemicalReaction):
    k_value: float

    def __init__(self,
                 k_value: float,
                 ):
        super().__init__()
        self.k_value = k_value

    @classmethod
    def from_dict(cls,
                  reactants: List[str],
                  products: List[str],
                  k_value: float,
                  list_sp: ListSpecies,
                  ) -> KReaction:
        """Create a KReaction from a dictionary

        Args:
            reactants (List[str]): List of reactants.
            products (List[str]): List of products.
            k_value (float): k value of the reactions.

        Returns:
            KReaction
        """
        new = cls(k_value)
        for sp in reactants:
            sim_sp = list_sp.get_or_create_dyn_species(sp)
            new.add_as_reactant(sim_sp)

        # Case of n*A -> B, then n*k is measured/given
        if len(new.reactants) == 1:
            # Compute k: divide by reactant's stoichio
            new.k_value /= new.stoi_reactants[sp]

        for sp in products:
            sim_sp = list_sp.get_or_create_dyn_species(sp)
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
        out.append(f"], kValue: {self.k_value:.2e}>")
        return "".join(out)
