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
from radiopyo.parser import RPARSER

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
    "KValType",
    "polynomial_factory",
]


# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #
def polynomial_factory(
    value: float | tp.List[float] | np.polynomial.Polynomial,
) -> np.polynomial.Polynomial:
    """ """
    if isinstance(value, np.polynomial.Polynomial):
        return value
    if isinstance(value, (list, tuple, np.ndarray)):
        if len(value) == 1:
            _coef = [np.log10(value[0])]
        elif len(value) > 1:
            _coef = value
        else:
            raise ValueError("Cannot make a polynomial with no coeff!")
    else:
        _coef = [np.log10(value)]

    return np.polynomial.Polynomial(_coef)


# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #

KValType: tp.TypeAlias = tp.Union[float, tp.List[float], np.polynomial.Polynomial]


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

    _DEFAULT_T: float = 25
    T_inv: float
    _k_value: tp.Optional[np.polynomial.Polynomial]

    def __init__(
        self,
        k_value: tp.Optional[KValType] = None,
        T: tp.Optional[float] = None,
    ):
        super().__init__()
        self.T_inv = 1.0 / T if T is not None else 1.0 / self._DEFAULT_T
        self._k_value = polynomial_factory(k_value) if k_value is not None else None

    def k_value(self, T: tp.Optional[float] = None) -> float:
        if self._k_value is None:
            raise ValueError(
                f"Unable to compute reaction constant. No data provided for reaction: {self}"
            )
        T_inv = 1.0 / T if T is not None else self.T_inv
        return np.power(10, self._k_value(T_inv))

    def set_k_value(self, value: KValType) -> None:
        """ """
        self._k_value = polynomial_factory(value)

    def set_T(self, value: float) -> None:
        self.T_inv = 1.0 / value

    def get_T(self) -> float:
        return 1.0 / self.T_inv

    def make_hash(self, method: str = "sha256") -> str:
        """ """
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
            return [
                f"KR_{self.index}",
            ]
        raise KeyError("Unknown reaction index. Cannot make a label without index.")

    def raw_k_value(self) -> float | tp.List[float]:
        """
        Produce the k_value as provided by config files.
        """
        if self._k_value is None:
            raise ValueError(
                f"Unable to compute reaction constant. No data provided for reaction: {self}"
            )
        if len(self._k_value) > 1:
            return self._k_value.coef

        # Case of n*A -> B, then n*k is measured/given
        if len(self.reactants) == 1:
            sp = next(iter(self.reactants.keys()))
            return self.k_value() * self.stoi_reactants[sp]
        return self.k_value()

    @classmethod
    def from_generic(
        cls,
        generic: GenericReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Create a KReaction from a 'GenericReaction', Typed-Dict structure produced by
        radiopyo module parser.
        TODO: k_value should probably be contained in the GenericReaction instance.
        """

        reactants: tp.List[str] = []
        products: tp.List[str] = []

        for sp, (_, stoi) in generic["reactants"].items():
            reactants += [sp] * stoi
        for sp, (_, stoi) in generic["products"].items():
            products += [sp] * stoi
        return cls.from_dict(
            reactants, products, generic["constants"]["k_value"], list_sp
        )

    @classmethod
    def from_kwargs(
        cls,
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
            case {"reactants": _, "products": _, "k_value": _} if len(kwargs) == 3:
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
    def from_dict(
        cls,
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
        if isinstance(k_value, (float, int)):
            k_value = [k_value]

        new = cls(k_value)
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        for sp in reactants:
            sim_sp = list_sp.get_or_create_dyn_species(sp.strip())
            new.add_as_reactant(sim_sp)

        # Case of n*A -> B, then n*k is measured/given
        if len(new.reactants) == 1 and len(k_value) == 1:
            # Compute k: divide by reactant's stoichio
            new.set_k_value(k_value[0] / new.stoi_reactants[sp])

        for sp in products:
            sim_sp = list_sp.get_or_create_dyn_species(sp.strip())
            new.add_as_product(sim_sp)
        return new

    def __str__(self) -> str:
        out = []
        for sp, stoi in self.iter_reactants():
            _ = "" if stoi <= 1 else f"{stoi}"
            out.append(_ + f"{sp.label}")
            out.append("+")
        out[-1] = "->"
        for sp, stoi in self.iter_products():
            _ = "" if stoi <= 1 else f"{stoi}"
            out.append(_ + f"{sp.label}")
            out.append("+")
        return " ".join(out[:-1])

    def __repr__(self) -> str:
        return f"KReaction('{self.__str__()}', k_value={self.raw_k_value()})"


class KReactionDict(tp.TypedDict, total=False):
    """ """

    reactants: tp.Optional[tp.List[str]]
    products: tp.Optional[tp.List[str]]
    reaction: tp.Optional[str]
    constants: tp.Optional[tp.Dict[str, float]]
    k_value: float | tp.List[float]


class KReactionCollection:
    """
    Attributes:
        reactions: tp.Dict[str, KReaction]
            All KReactions belonging to the collection.
    """

    _reactions: tp.List[tp.Tuple[str, KReaction]]
    _indexer: tp.Dict[str, int]

    def __init__(self) -> None:
        self._reactions = []
        self._indexer = {}

    def __len__(self) -> int:
        return len(self._reactions)

    def __iter__(self) -> tp.Iterator[tp.Tuple[str, KReaction]]:
        return iter(self._reactions)

    @property
    def last(self) -> KReaction:
        return self._reactions[-1][1]

    def iter_reactions_and_kvalues(
        self,
        T: tp.Optional[float] = None,
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """ """
        for label, reaction in self:
            yield label, reaction, reaction.k_value(T)

    def _get_reaction_index(self, label: str) -> int:
        return self._indexer[label]

    def get_reaction(self, label: str) -> KReaction:
        if label not in self._indexer:
            raise KeyError(f"No reaction labeled {label} in the reaction collection")
        return self._reactions[self._get_reaction_index(label)][1]

    def add_reaction(
        self,
        label: str,
        generic: GenericReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> None:
        """
        Push a new Reaction in the Collection.
        """
        if label in self._indexer:
            raise KeyError(
                f"Already existing reaction with label {label} in the collection"
            )
        self._indexer[label] = len(self)
        self._reactions.append(
            (
                label,
                KReaction.from_generic(
                    generic,
                    list_sp=list_sp,
                ),
            )
        )
