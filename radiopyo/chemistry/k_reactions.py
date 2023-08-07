#!/usr/bin/env python3
"""
The k_reactions submodule contains all class definitions related to basic chemical
reactions defined via a reaction constant.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import hashlib
import typing as tp

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.utils import RPARSER

from .base_chemistry import BaseChemicalReaction
from .species import ListSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.parser.reaction_parser import GenericReaction

# ------------------------------------------------------------------------------------ #

__all__ = [
    "KReaction",
    "KReactionDict",
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

    def raw_k_value(self) -> float:
        """
        Produce the k_value as provided by config files.
        """
        # Case of n*A -> B, then n*k is measured/given
        if len(self.reactants) == 1:
            sp = list(self.reactants.keys())[0]
            return self.k_value * self.stoi_reactants[sp]
        return self.k_value

    def make_hash(self, method: str = "sha256") -> str:
        """
        """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(self).encode())
        return hasher.hexdigest()

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
    def from_generic(cls,
                     generic: GenericReaction,
                     k_value: float,
                     list_sp: tp.Optional[ListSpecies] = None,) -> KReaction:
        """
        Create a KReaction from a 'GenericReaction', Typed-Dict structure produced by
        radiopyo module parser.
        TODO: k_value should probably be contained in the GenericReaction instance.
        """

        reactants: tp.List[str] = []
        products: tp.List[str] = []

        for sp, (_, stoi) in generic["reactants"].items():
            reactants += [sp]*stoi
        for sp, (_, stoi) in generic["products"].items():
            products += [sp]*stoi
        return KReaction.from_dict(reactants, products, k_value, list_sp)

    @classmethod
    def from_kwargs(cls,
                    list_sp: tp.Optional[ListSpecies] = None,
                    **kwargs: str | tp.List[str] | float, ) -> KReaction:
        """ 
        This method combines the two possible definitions of a KReaction in a config
        file.
        TODO: remove the dict approach to keep only the single string reaction.
        """
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        match kwargs:
            case {"reactants": _, "products": _, "k_value": _}\
                    if len(kwargs) == 3:
                return KReaction.from_dict(
                    list_sp=list_sp,
                    **kwargs,  # type:ignore[arg-type]
                )
            case {"reaction": r, "k_value": k} if len(kwargs) == 2:
                gen = RPARSER.parse_reaction(r)  # type:ignore[arg-type]
                return KReaction.from_generic(
                    gen,
                    k,  # type:ignore[arg-type]
                    list_sp=list_sp,
                )

        raise ValueError(f"Wrong KReaction syntax {kwargs}")

    @classmethod
    def from_dict(cls,
                  reactants: tp.List[str],
                  products: tp.List[str],
                  k_value: float,
                  list_sp: tp.Optional[ListSpecies] = None,
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
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

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
        out = []
        for sp, stoi in self.iter_reactants():
            _ = "" if stoi <= 1 else f"{stoi}"
            out.append(_+f"{sp.label}")
            out.append("+")
        out[-1] = "->"
        for sp, stoi in self.iter_products():
            _ = "" if stoi <= 1 else f"{stoi}"
            out.append(_+f"{sp.label}")
            out.append("+")
        return " ".join(out[:-1])

    def __repr__(self) -> str:
        return (f"KReaction('{self.__str__()}', "
                f"k_value={round(self.k_value,5)})")


class KReactionDict(tp.TypedDict, total=False):
    """ """
    reactants: tp.Optional[tp.List[str]]
    products: tp.Optional[tp.List[str]]
    reaction: tp.Optional[str]
    constants: tp.Optional[tp.Dict[str, float]]
    k_value: float
