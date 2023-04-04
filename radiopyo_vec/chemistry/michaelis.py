from __future__ import annotations

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

    from .species import ListSpecies, SimSpecies

__all__ = [
    "MichaelisMenten",
]

# ---------------------------------------------------------------------------- #
#                             K REACTION DEFINITION                            #
# ---------------------------------------------------------------------------- #


class MichaelisMenten(BaseChemicalReaction):
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
        """Create a KReaction from a dictionary

        Args:
            reactants (List[str]): List of reactants.
            products (List[str]): List of products.
            k_value (float): k value of the reactions.

        Returns:
            KReaction
        """
        new = cls(k_value, k_micha)

        new.enzyme = list_sp.get_or_create_cst_species(enzyme)
        if isinstance(substrate, str):
            substrate = [substrate]
        elif len(substrate) > 1 and len(set(substrate)) > 1:
            raise TypeError(
                f"Substrate in Michaelis Menten reaction"
                f" must all be the same if provide as list. Found: {substrate}"
            )

        new.substrate = list_sp.get_or_create_dyn_species(substrate[0])

        for sp in substrate:
            new.add_as_reactant(list_sp.get_or_create_dyn_species(sp))

        for sp in products:
            new.add_as_product(list_sp.get_or_create_dyn_species(sp))
        return new

    def as_label(self) -> str:
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
