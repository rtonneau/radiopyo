from __future__ import annotations

from functools import reduce
from operator import mul
from typing import (
    TYPE_CHECKING,
    List,
)

import numpy as np

from radiopyo.chemistry.base_chemistry import (
    BaseChemicalReaction,
    ReactionSpecies,
    ReactionSpeciesEnum,
)

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.utils.collections import MappedSpecies

__all__ = [
    "MichaelisMenten",
]

# ---------------------------------------------------------------------------- #
#                             K REACTION DEFINITION                            #
# ---------------------------------------------------------------------------- #


class MichaelisMenten(BaseChemicalReaction):
    k_value: float
    k_micha: float
    _substrate: ReactionSpecies
    _enzyme: ReactionSpecies

    def __init__(self,
                 species: List[ReactionSpecies],
                 stoichio: List[int],
                 k_value: float,
                 k_micha: float,
                 ):
        super().__init__(species)
        self.k_value = k_value
        self.k_micha = k_micha
        self.stoichio = stoichio

    @classmethod
    def new_empty(cls) -> MichaelisMenten:
        """new_empty _summary_

        Returns:
            MichaelisMenten: _description_
        """
        return cls([], [], 0.0, 0.0,)

    @classmethod
    def from_dict(cls,
                  enzyme: str,
                  substrate: List[str] | str,
                  products: List[str],
                  k_value: float,
                  k_micha: float,
                  ) -> MichaelisMenten:
        """Create a KReaction from a dictionary

        Args:
            reactants (List[str]): List of reactants.
            products (List[str]): List of products.
            k_value (float): k value of the reactions.

        Returns:
            KReaction
        """
        new = cls.new_empty()
        new.set_k_value(k_value)
        new.set_k_micha(k_micha)
        new._enzyme = ReactionSpecies.new_enzyme(enzyme)
        if isinstance(substrate, str):
            substrate = [substrate]
        elif len(substrate) > 1 and len(set(substrate)) > 1:
            raise TypeError(
                f"Substrate in Michaelis Menten reaction"
                f" must all be the same if provide as list. Found: {substrate}"
            )

        new._substrate = ReactionSpecies.new_reactant(substrate[0])

        for sp in substrate:
            new.add_as_reactant(sp)

        for sp in products:
            new.add_as_product(sp)
        return new

    def set_k_value(self, k_value: float) -> None:
        """set_k_value _summary_

        Args:
            k_value (float): _description_
        """
        self.k_value = float(k_value)

    def set_k_micha(self, k_micha: float) -> None:
        """set_k_value _summary_

        Args:
            k_value (float): _description_
        """
        self.k_micha = float(k_micha)

    def compute_reaction(self,
                         dose_rate: float,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        """Compute the reaction rate of an enzymatic reaction.
        Reaction rate: R = k_value * [enzyme] * [substrate] / (k_m + [substrate])
        Enzyme => UntrackedSpecies.
        """
        res = species[self._enzyme.label].cc * \
            species[self._substrate.label].cc
        res /= (self.k_micha + species[self._substrate.label].cc)
        return self._rate_vec * self.k_value * res

    def compute_jacobian(self,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        """
        d_R/d_subs = k_val * [enzyme] * ([subs]+k_m - [subs]) / ([subs]+k_m)**2
        d_R/d_subs = k_val * [enzyme] * k_m / ([subs]+k_m)**2

        Args:
            species   (MappedSpecies): _description_

        Returns:
            NDArray[np.float64]: _description_
        """
        res = self.k_value * self.k_micha * species[self._enzyme.label].cc
        res /= (species[self._substrate.label].cc+self.k_micha)
        res /= (species[self._substrate.label].cc+self.k_micha)
        # d_R/dCt => (d_R/dsp) * (dsp/dCt)
        return self._rate_vec * species[self._substrate.label].ab_derive * res

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
