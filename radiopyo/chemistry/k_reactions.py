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
import operator
import typing as tp
import uuid
from abc import ABC, abstractmethod

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.parser import RPARSER

from .basic_chemistry import BasicReaction
from .collections import ListSpecies
from .utils.units import Temperature

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
        _coef = [-np.inf] if value <= 0 else [np.log10(value)]

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
        T: Temperature in °C (optional, default:25).
    """

    T: Temperature
    _k_value: tp.Optional[np.polynomial.Polynomial]

    def __init__(
        self,
        k_value: tp.Optional[KValType] = None,
        T: tp.Optional[float] = None,
    ):
        super().__init__()
        self.T = Temperature(T, unit="celsius")
        self._k_value = polynomial_factory(k_value) if k_value is not None else None

    def k_value(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Args:
            T: Temperature in °C
        """
        if self._k_value is None:
            raise ValueError(
                f"Unable to compute reaction constant. No data provided for reaction: {self}"
            )
        _T = Temperature(T) if T is not None else Temperature()
        return np.power(10, self._k_value(1.0 / _T.kelvin))

    def set_k_value(self, value: KValType) -> None:
        """ """
        self._k_value = polynomial_factory(value)

    def set_T(self, value: float, unit: str = "celsius") -> None:
        self.T = Temperature(value, unit)

    def get_T(self, unit: str = "celsius") -> float:
        return self.T.get(unit)

    def make_hash(self, method: str = "sha256") -> str:
        """ """
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
        if self.index is not None and self.index >= 0:
            return f"KR_{self.index}"

        raise KeyError("Unknown reaction index. Cannot make a label without index.")

    def raw_k_value(self) -> float | tp.List[float] | None:
        """
        Produce the k_value as provided by config files.
        """
        if self._k_value is None:
            return None
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
            reactants,
            products,
            generic["constants"].get("k_value", None),
            list_sp,
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
        k_value: None | float | tp.List[float],
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
        if k_value is not None and len(new.reactants) == 1 and len(k_value) == 1:
            # Compute k: divide by reactant's stoichio
            new.set_k_value(k_value[0] / new.stoi_reactants[sp])

        for sp in products:
            sim_sp = list_sp.get_or_create_dyn_species(sp.strip())
            new.add_as_product(sim_sp)
        return new

    def __str__(self) -> str:
        out = []
        out.append(self.reactants_as_string())
        out.append("->")
        out.append(self.products_as_string())
        return " ".join(out)

    def __repr__(self) -> str:
        return f"KReaction('{self.__str__()}', k_value={self.raw_k_value()})"


class Arrhenius(KReaction):
    """
    Args:
        A: float
            pre-exponential factor [L/mol/s]
        Ea: float
            Activation Energy [kJ/mol]
    """

    _pre_exp: float
    _activation_energy: float

    def __init__(
        self,
        A: float,
        Ea: float,
        T: tp.Optional[float] = None,
    ):
        """ """
        super().__init__(k_value=None, T=T)
        self._pre_exp = A
        self._activation_energy = Ea

    def k_value(self, T: Temperature | float | None = None) -> float:
        """
        Args:
            T: float | Temperature, optional (default=25°C)
                Temperature value in °C
        """
        T = Temperature(T)
        return self._pre_exp * np.exp(-self._activation_energy / T.kelvin)

    def raw_k_value(self) -> tp.List[float]:
        """
        Produce the k_value as provided by config files.
        """
        return [self._pre_exp, self._activation_energy]

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
            reactants,
            products,
            [
                generic["constants"]["A"],  # type: ignore [list-item]
                generic["constants"]["Ea"],  # type: ignore [list-item]
            ],
            list_sp,
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
            case {"reactants": _, "products": _, "A": _, "Ea": _} if len(kwargs) == 4:
                return cls.from_dict(
                    list_sp=list_sp,
                    **kwargs,  # type:ignore[arg-type]
                )
            case {"reaction": r, "A": _, "Ea": _} if len(kwargs) == 3:
                gen = RPARSER.parse_reaction(r)  # type: ignore[arg-type]
                gen["constants"]["A"] = kwargs["A"]  # type: ignore[assignment]
                gen["constants"]["Ea"] = kwargs["Ea"]  # type: ignore[assignment]

                return cls.from_generic(
                    gen,
                    list_sp=list_sp,
                )

        raise ValueError(f"Wrong Arrhenius Reaction syntax {kwargs}")

    @classmethod
    def from_dict(
        cls,
        reactants: tp.List[str],
        products: tp.List[str],
        k_value: None | float | tp.List[float],
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """Create a KReaction from a dictionary (via destructuring)

        Args:
            reactants: List[str]
                List of reactants
            products: List[str]
                List of products
            k_value: tp.List[float]
                k value of the reactions or the coefficients of the poly version
            list_sp: ListSpecies
                Collection of species already involved in the simulation. This
                collection will grow if some reactants and/or products do not exist yet.

        Returns:
            Arrhenius
        """
        if k_value is None or isinstance(k_value, float):
            raise TypeError
        new = cls(k_value[0], k_value[1])

        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        for sp in reactants:
            sim_sp = list_sp.get_or_create_dyn_species(sp.strip())
            new.add_as_reactant(sim_sp)

        for sp in products:
            sim_sp = list_sp.get_or_create_dyn_species(sp.strip())
            new.add_as_product(sim_sp)
        return new


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

    def __iter__(self) -> tp.Iterator[KReaction]:
        return iter(map(operator.itemgetter(1), self._reactions))

    def __getitem__(self, index: int) -> KReaction:
        return self._reactions[index][1]

    @property
    def last(self) -> KReaction:
        return self._reactions[-1][1]

    def merge_with(
        self,
        other: tp.Self,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> None:
        """ """
        for label, reaction in other.iter_labels_and_reactions():
            self.add_reaction(reaction, label=label, list_sp=list_sp)

    def update_label(self, current: str, new: str) -> tp.Self:
        """ """
        try:
            idx = self._indexer[current]
        except KeyError as err:
            raise KeyError(
                f"Reaction label: {current} not in KReactionCollection"
            ) from err
        if new in self._indexer:
            raise KeyError(f"Reaction label: {new} already in KReactionCollection")

        del self._indexer[current]
        self._indexer[new] = idx
        reaction = self._reactions[idx][1]
        self._reactions.append((new, reaction))
        del self._reactions[idx]
        return self

    def sort_by_label(self) -> tp.Self:
        """
        May be useful to sort the list of reactions so that iteration is not random.
        """
        self._reactions = sorted(self._reactions, key=operator.itemgetter(0))
        return self

    def sort_by_reaction_index(self) -> tp.Self:
        """
        May be useful to sort the list of reactions so that iteration is not random.
        """
        self._reactions = sorted(self._reactions, key=lambda x: x[1].get_index())
        return self

    def iter_labels_and_reactions(self) -> tp.Iterator[tp.Tuple[str, KReaction]]:
        """ """
        return iter(self._reactions)

    def iter_reactions_and_kvalues(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """
        Args:
            T: Temperature object (optional, default:25).
        """
        for label, reaction in self.iter_labels_and_reactions():
            yield (
                label,
                reaction,
                reaction.k_value(T),
            )

    def _get_reaction_index(self, label: str) -> int:
        return self._indexer[label]

    def get_reaction(self, label: str) -> KReaction:
        if label not in self._indexer:
            raise KeyError(f"No reaction labeled {label} in the reaction collection")
        return self._reactions[self._get_reaction_index(label)][1]

    def add_reaction(
        self,
        reaction: GenericReaction | KReaction,
        label: tp.Optional[str] = None,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Push a new Reaction in the Collection.
        """
        if label is None:
            label = str(uuid.uuid4())  # Make a unique ID
        if label in self._indexer:
            raise KeyError(
                f"Already existing reaction with label {label} in the collection"
            )
        self._indexer[label] = len(self)
        if isinstance(reaction, dict):
            reaction = KReaction.from_generic(
                reaction,
                list_sp=list_sp,
            )
        self._reactions.append((label, reaction))
        return self


class BiReaction(KReactionCollection, ABC):
    """ """

    _k_backward: tp.Optional[np.polynomial.Polynomial]
    _k_forward: tp.Optional[np.polynomial.Polynomial]

    @property
    def forward(self) -> KReaction:
        """ """
        return self._reactions[0][1]

    @property
    def backward(self) -> KReaction:
        """ """
        return self._reactions[1][1]

    def add_reaction(
        self,
        reaction: GenericReaction | KReaction,
        label: tp.Optional[str] = None,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        raise NotImplementedError(
            "Not possible to add more reaction to an equilibrium collection"
        )

    def iter_reactions_and_kvalues(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        yield "forward", self.get_reaction("forward"), self.k_forward(T)
        yield "backward", self.get_reaction("backward"), self.k_backward(T)

    def __str__(self) -> str:
        """ """
        out = []
        out.append(self.forward.reactants_as_string())
        out.append("<->")
        out.append(self.forward.products_as_string())
        return " ".join(out)

    def add_forward_reaction(
        self,
        reaction: GenericReaction | KReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        super().add_reaction(reaction, label="forward", list_sp=list_sp)
        return self

    def add_backward_reaction(
        self,
        reaction: GenericReaction | KReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        super().add_reaction(reaction, label="backward", list_sp=list_sp)
        return self

    @abstractmethod
    def k_equilibrium(self, T: tp.Optional[Temperature | float] = None) -> float:
        """ """

    @abstractmethod
    def k_forward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """ """

    @abstractmethod
    def k_backward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """ """
