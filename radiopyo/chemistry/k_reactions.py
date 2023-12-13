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

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.utils import RPARSER

from .basic_chemistry import BasicReaction
from .collections import ListSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.parser.reaction_parser import GenericReaction

# ------------------------------------------------------------------------------------ #

__all__ = [
    "KReaction",
    "KReactionDict",
    "KReactionCollection",
    "KValType"
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #

KValType: tp.TypeAlias = tp.Union[float,
                                  tp.List[float],
                                  np.polynomial.Polynomial,
                                  tp.List[np.polynomial.Polynomial],
                                  ]


class KReaction(BasicReaction):
    """
    Almost all the logic is already included in the BaseChemicalReaction class.
    Add the reaction rate as attribute.
    For reactions where only one reactant is involved e.g.
        2 H_r -> H2
    in that case, 2k is usually given/measured instead of k (reaction constant).
    User should simply pass the 2k value at construction and then it will be accounted
    for during calculation of reaction rate.

    Attributes:
        k_value: List[np.polynomial.Polynomial, ]
            Reaction rate
    """
    k_value: tp.List[np.polynomial.Polynomial]

    def __init__(self,
                 k_value: KValType,
                 ):
        super().__init__()
        self.k_value = []
        if isinstance(k_value, (tuple, list)):
            if len(k_value) == 0:
                raise ValueError("Len 0 for k_value passed as a collection")
            for elt in k_value:
                if not isinstance(elt, np.polynomial.Polynomial):
                    elt = np.polynomial.Polynomial(elt)
                self.k_value.append(elt)
        elif isinstance(k_value, np.polynomial.Polynomial):
            self.k_value.append(k_value)
        else:
            self.k_value = [np.polynomial.Polynomial(k_value), ]

    def make_hash(self, method: str = "sha256") -> str:
        """
        """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(self).encode())
        return hasher.hexdigest()

    def as_label(self) -> tp.List[str]:
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
            return [f"KR_{self.index}",]
        raise KeyError(
            "Unknown reaction index. Cannot make a label without index.")

    def raw_k_value(self) -> float | tp.List[float]:
        """
        Produce the k_value as provided by config files.
        """
        if len(self.k_value) > 1:
            return self.k_value.coef

        # Case of n*A -> B, then n*k is measured/given
        if len(self.reactants) == 1:
            sp = next(iter(self.reactants.keys()))
            return self.k_value(25) * self.stoi_reactants[sp]
        return self.k_value(25)

    @classmethod
    def from_generic(cls,
                     generic: GenericReaction,
                     list_sp: tp.Optional[ListSpecies] = None,) -> tp.Self:
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
        return cls.from_dict(reactants,
                             products,
                             generic["constants"]["k_value"],
                             list_sp)

    @classmethod
    def from_kwargs(cls,
                    list_sp: tp.Optional[ListSpecies] = None,
                    **kwargs: str | tp.List[str] | float | tp.List[float],
                    ) -> tp.Self:
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
                return cls.from_dict(
                    list_sp=list_sp,
                    **kwargs,  # type:ignore[arg-type]
                )
            case {"reaction": r, "k_value": k} if len(kwargs) == 2:
                gen = RPARSER.parse_reaction(r)  # type:ignore[arg-type]
                gen["constants"]["k_value"] = k  # type:ignore[assignment]
                return cls.from_generic(
                    gen,
                    list_sp=list_sp,
                )

        raise ValueError(f"Wrong KReaction syntax {kwargs}")

    @classmethod
    def from_dict(cls,
                  reactants: tp.List[str],
                  products: tp.List[str],
                  k_value: float | tp.List[float],
                  list_sp: tp.Optional[ListSpecies] = None,
                  ) -> tp.Self:
        """Create a KReaction from a dictionary (via destructuring)

        Args:
            reactants: List[str]
                List of reactants
            products: List[str]
                List of products
            k_value: float | tp.List[float]
                k value of the reactions or the coefficients of the poly version
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
        if len(new.reactants) == 1 and len(new.k_value) > 1:
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
                f"k_value={self.raw_k_value()})")


class KReactionDict(tp.TypedDict, total=False):
    """ """
    reactants: tp.Optional[tp.List[str]]
    products: tp.Optional[tp.List[str]]
    reaction: tp.Optional[str]
    constants: tp.Optional[tp.Dict[str, float]]
    k_value: float | tp.List[float]


class KReactionCollection():
    """
    Attributes:
        reactions: tp.Dict[str, KReaction]
            All KReactions belonging to the collection.
    """
    reactions: tp.Dict[str, KReaction]

    def __init__(self) -> None:
        self.reactions = {}

    def __len__(self) -> int:
        return len(self.reactions)

    def __iter__(self) -> tp.Iterator[KReaction]:
        return iter(self.reactions.values())

    def items(self) -> tp.Iterator[tp.Tuple[str, KReaction]]:
        return iter(self.reactions.items())

    def add_reaction(self,
                     label: str,
                     generic: GenericReaction,
                     list_sp: tp.Optional[ListSpecies] = None,
                     ) -> None:
        """
        Push a new Reaction in the Collection.
        """
        self.reactions[label] = KReaction.from_generic(generic,
                                                       list_sp=list_sp, )
