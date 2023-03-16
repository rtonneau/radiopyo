from __future__ import annotations

from functools import reduce
from itertools import repeat
from operator import mul
from typing import (
    TYPE_CHECKING,
    List,
)

import numpy as np

from radiopyo_vec.chemistry.base_chemistry import (
    BaseChemicalReaction,
)

from radiopyo_vec.chemistry.species import ListSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo_vec.utils.collections import MappedSpecies

__all__ = [
    "ReactionSpeciesEnum",
    "ReactionSpecies",
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
            sim_sp = list_sp.get_or_create_sim_species(sp)
            new.add_as_reactant(sim_sp)

        # Case of n*A -> B, then n*k is measured/given
        if len(new.reactants) == 1:
            # Compute k: divide by reactant's stoichio
            new.k_value /= new.stoi_reactants[0]

        for sp in products:
            sim_sp = list_sp.get_or_create_sim_species(sp)
            new.add_as_product(sim_sp)
        return new

    def compute_reaction(self,
                         dose_rate: float,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        """Compute the reaction rate of a reaction"""

        res = reduce(mul, [species[sp.label].cc**stoi
                           for sp, stoi in self.iter_reactants()])

        return self._rate_vec * res * self.k_value

    def compute_jacobian(self,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        """_summary_

        Args:
            species   (Dict[str, float]): _description_
            indexes   (Dict[str, int]): _description_
            ab_derive (Dict[str, float]): _description_

        Returns:
            NDArray[np.float64]: _description_
        """
        out = np.zeros(species.number_tracked_species)
        # First loop is the derivation species
        for sp, stoi in self.iter_reactants():
            idx = species[sp.label].index
            if idx is None:  # idx is None (not a tracked species then)
                continue
            # Actual derive
            out[idx] = self.k_value * stoi * \
                np.power(species[sp.label].cc, stoi-1)
            # d_R/dCt => (d_R/dsp) * (dsp/dCt)
            out[idx] *= species[sp.label].ab_derive

            # Loop on the other species
            for sp2, stoi2 in filter(lambda x: x[0] != sp, self.iter_reactants()):
                out[idx] *= (species[sp2.label].cc ** stoi2)
        return out

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
