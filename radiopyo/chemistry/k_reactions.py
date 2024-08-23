#!/usr/bin/env python3
"""
The k_reactions submodule contains all class definitions related to basic chemical
reactions defined via a reaction constant.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import hashlib
import operator
import typing as tp
import uuid
from abc import ABC, abstractmethod

import numpy as np

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser import RPARSER

from .basic_reaction import BasicReaction
from .species_collection import ListSpecies
from .utils.units import Temperature

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.parser.reaction_parser import GenericReaction

# --------------------------------------------------------------------------- #

__all__ = [
    "KReaction",
    "KReactionDict",
    "KReactionCollection",
    "KValType",
    "polynomial_factory",
]


# --------------------------------------------------------------------------- #
#                              FUNCTION DEFINITIONS                           #
# --------------------------------------------------------------------------- #
def polynomial_factory(
    value: float | tp.List[float] | np.polynomial.Polynomial,
) -> np.polynomial.Polynomial:
    """
    A factory function to create a numpy polynomial from a given value.

    Args:
        value (float | tp.List[float] | np.polynomial.Polynomial):
            The value to be converted into a numpy polynomial. It can be
            a float, a list of floats or a numpy polynomial.

    Returns:
        np.polynomial.Polynomial:
            A numpy polynomial created from the given value.

    Raises:
        ValueError:
            If the given value is not a valid input for creating a
            numpy polynomial.
    """
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


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #

KValType: tp.TypeAlias = tp.Union[float, tp.List[float], np.polynomial.Polynomial]


class KReaction(BasicReaction):
    """
    Almost all the logic is already included in the BaseChemicalReaction
    class. Add the reaction rate as attribute.For reactions where only
    one reactant is involved e.g.
        2 H_r -> H2
    in that case, 2k is usually given/measured instead of
    k (reaction constant). User should simply pass the 2k value at
    construction and then it will be accounted for during calculation of
    reaction rate.

    Attributes:
        k_value: List[np.polynomial.Polynomial, ]
            Reaction rate
        T: Temperature in °C (optional, default:25).
    """

    T: Temperature
    _k_value: tp.Optional[np.polynomial.Polynomial]
    _OPE_STRING: str = "->"

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
        Returns the reaction constant at a given temperature.

        Args:
            T (Temperature | float, optional): Temperature in °C.
                Defaults to None.

        Returns:
            float: The reaction constant at the given temperature.

        Raises:
            ValueError: If no constant reaction value is provided for
                the reaction.
        """
        if self._k_value is None:
            raise ValueError(
                "Unable to compute reaction constant. No data provided"
                f" for reaction: {self}"
            )
        _T = Temperature(T) if T is not None else Temperature()
        return np.power(10, self._k_value(1.0 / _T.kelvin))

    def set_k_value(self, value: KValType) -> None:
        """
        Set the reaction constant value.

        Args:
            value (KValType): The new reaction constant value.

        Returns:
            None
        """
        self._k_value = polynomial_factory(value)

    def set_T(self, value: float, unit: str = "celsius") -> None:
        """
        Sets the temperature value of the reaction.

        Args:
            value (float): The new temperature value.
            unit (str, optional): The unit of the temperature value.
                Defaults to "celsius".

        Returns:
            None
        """
        self.T = Temperature(value, unit)

    def get_T(self, unit: str = "celsius") -> float:
        """
        Returns the temperature value of the reaction in the specified
        unit.

        Args:
            unit (str, optional): The unit of the temperature value.
                Defaults to "celsius".

        Returns:
            float: The temperature value in the specified unit.
        """
        return self.T.get(unit)

    def make_hash(self, method: str = "sha256") -> str:
        """
        Produce a unique identifier (content based) for easy comparison.

        Args:
            method (str): Hash method to use (default: "sha256").

        Returns:
            str: A hash value representing the object.
        """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(self).encode())
        return hasher.hexdigest()

    def as_label(self) -> str:
        """
        This method generates a short label for the reaction based on
        its index.

        Returns:
            str: A string label in the format "KR_{index}" if index is
                valid.

        Raises:
            KeyError: If the index is None or invalid.
        """
        if self.index is not None and self.index >= 0:
            return f"KR_{self.index}"

        raise KeyError("Unknown reaction index. Cannot make a label without index.")

    def raw_k_value(self) -> float | tp.List[float] | None:
        """
        Produce the k_value as provided by config files.

        Returns:
            float | tp.List[float] | None: k_value as provided by
                config files.
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
        Create a KReaction from a 'GenericReaction', Typed-Dict
        structure produced by radiopyo module parser.

        Args:
            generic (GenericReaction): Typed-Dict structure holding the
                reaction metadata.
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.

        Returns:
            KReaction: An instance of KReaction.
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
        Create a KReaction from keyword arguments.

        Args:
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.
            **kwargs: Keyword arguments to create the instance.

        Returns:
            KReaction: An instance of KReaction.

        Raises:
            ValueError: If the keyword arguments do not match the
                expected syntax.
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
        """
        Create a KReaction from a dictionary (via destructuring).

        Args:
            reactants: List[str]
                List of reactants
            products: List[str]
                List of products
            k_value: float | tp.List[float]
                k value of the reactions or the coefficients of the
                    poly version
            list_sp: ListSpecies
                Collection of species already involved in the
                simulation. This collection will grow if some reactants
                and/or products do not exist yet.

        Returns:
            KReaction
                An instance of KReaction.
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
        """
        Returns a string representation of the reaction.

        This method combines the reactants, products and the operation
        string to create a string that represents the reaction.

        Returns:
            str: A string representation of the reaction.
        """
        out = []
        out.append(self.reactants_as_string())
        out.append(self._OPE_STRING)
        out.append(self.products_as_string())
        return " ".join(out)

    def __repr__(self) -> str:
        return f"KReaction('{self.__str__()}', k_value={self.raw_k_value()})"


class Arrhenius(KReaction):
    """
    Class representing an Arrhenius reaction.

    Args:
        A: float
            The pre-exponential factor of the reaction in [L/mol/s].
        Ea: float
            The activation energy of the reaction in [kJ/mol].

    Attributes:
        _pre_exp: float
            The pre-exponential factor of the reaction.
        _activation_energy: float
            The activation energy of the reaction.
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
        Produce the k_value of the Arrhenius reaction at a given
        temperature.

        Args:
            T: Temperature | float | None
                The temperature at which the k_value is calculated. If
                None, the default temperature is used.

        Returns:
            float
                The k_value of the Arrhenius reaction at the given
                temperature.
        """
        T = Temperature(T)
        return self._pre_exp * np.exp(-self._activation_energy / T.kelvin)

    def raw_k_value(self) -> tp.List[float]:
        """
        Produce the k_value as provided by config files.

        Returns:
            List[float]: A list containing the pre-exponential factor
            and the activation energy.
        """
        return [self._pre_exp, self._activation_energy]

    @classmethod
    def from_generic(
        cls,
        generic: GenericReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Create an Arrhenius from a 'GenericReaction', Typed-Dict
        structure produced by radiopyo module parser.

        Args:
            generic (GenericReaction): Typed-Dict structure holding the
                reaction metadata.
            list_sp (ListSpecies): Collection of species already
                involved in the simulation. Defaults to None.

        Returns:
            Arrhenius: An instance of Arrhenius.
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
        This method combines the two possible definitions of an
        Arrhenius in a config file.
        Args:
            list_sp (ListSpecies): Collection of species already
                involved in the simulation. Defaults to None.
            **kwargs: Keyword arguments to create the instance.
                - reactants (ListSpecies): Collection of reactants in
                    the reaction.
                - products (ListSpecies): Collection of products in
                    the reaction.
                - A (float): Pre-exponential factor.
                - Ea (float): Activation energy.
                - reaction (str): Reaction string.

        Returns:
            Arrhenius: An instance of Arrhenius.
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
        """
        Create an Arrhenius from a dictionary (via destructuring).

        Args:
            reactants (List[str]): List of reactants
            products (List[str]): List of products
            k_value (float | List[float]): k value of the reactions or
                the coefficients of the poly version
            list_sp (ListSpecies): Collection of species already
                involved in the simulation. This collection will grow if
                some reactants and/or products do not exist yet.

        Returns:
            Arrhenius: An instance of Arrhenius.
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


class CstReaction(KReaction):
    """ """

    def as_label(self) -> str:
        """
        This method generates a short label for the reaction based on
        its index.

        Returns:
            str: A string label in the format "CR_{index}" if
                index is valid.

        Raises:
            KeyError: If the index is None or invalid.
        """
        if self.index is not None and self.index >= 0:
            return f"CR_{self.index}"

        raise KeyError("Unknown reaction index. Cannot make a label without index.")

    @classmethod
    def from_generic(
        cls,
        generic: GenericReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Create a CstReaction from a 'GenericReaction', Typed-Dict
        structure produced by radiopyo module parser.

        Args:
            generic (GenericReaction): Typed-Dict structure holding the
                reaction metadata.
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.

        Returns:
            CstReaction: An instance of CstReaction.
        """
        products: tp.List[str] = []

        for sp, (_, stoi) in generic["products"].items():
            products += [sp] * stoi
        return cls.from_dict(
            [],
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
        Create a CstReaction from keyword arguments.

        Args:
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.
            **kwargs: Keyword arguments to create the instance.

        Returns:
            CstReaction: An instance of CstReaction.

        Raises:
            ValueError: If the keyword arguments do not match the
                expected syntax.
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

    def __str__(self) -> str:
        """
        Returns a string representation of the reaction.

        This method combines the reactants, products and the operation
        string to create a string that represents the reaction.

        Returns:
            str: A string representation of the reaction.
        """
        out = []
        out.append(self._OPE_STRING)
        out.append(self.products_as_string())
        return " ".join(out)

    def __repr__(self) -> str:
        """
        Returns a string representation of the CstReaction instance.

        The representation includes the reaction string and the k_value.

        Returns:
            str: A string representation of the CstReaction instance.
        """
        return f"CstReaction('{self.__str__()}', k_value={self.raw_k_value()})"

    def tot_number_reactants(self) -> int:
        """
        Count the total number of reactants involved in the reaction
        => Stoichiometry dependent! In other words, it counts the number
        of molecules/entities defined as reactants.

        Returns:
            int : The total number of reactants i.e. 0 as there is
                no reactant

        Notes:
            For a CstReaction, there is no reactant needed!
        """
        return 0


class KReactionDict(tp.TypedDict, total=False):
    """
    TypedDict structure holding the metadata of a KReaction like object.

    Attributes:
        reactants (tp.Optional[tp.List[str]]): List of reactants in
            the reaction.
        products (tp.Optional[tp.List[str]]): List of products in
            the reaction.
        reaction (tp.Optional[str]): String representation of
            the reaction.
        constants (tp.Optional[tp.Dict[str, float]]): Dictionary of
            constants in the reaction.
        k_value (float | tp.List[float]): Value of the reaction
            rate constant.
    """

    reactants: tp.Optional[tp.List[str]]
    products: tp.Optional[tp.List[str]]
    reaction: tp.Optional[str]
    constants: tp.Optional[tp.Dict[str, float]]
    k_value: float | tp.List[float]


class KReactionCollection:
    """
    A collection of KReaction objects.

    Attributes:
        _reactions (tp.List[tp.Tuple[str, KReaction]]):
            A list of tuples containing the unique label and the
            KReaction instance. The label is used to identify and
            manipulate the reaction instance.
        _indexer (tp.Dict[str, int]):
            A dictionary mapping each label to its index in the
            _reactions list.
    """

    _reactions: tp.List[tp.Tuple[str, KReaction]]
    _indexer: tp.Dict[str, int]

    def __init__(self) -> None:
        self._reactions = []
        self._indexer = {}

    def __len__(self) -> int:
        """
        Returns the total number of reactions in the collection.

        Returns:
            int: The number of reactions in the collection.
        """
        return len(self._reactions)

    def __iter__(self) -> tp.Iterator[KReaction]:
        """
        Returns an iterator over the KReaction instances in the
        collection.

        Yields:
            KReaction: The KReaction instance.
        """
        return iter(map(operator.itemgetter(1), self._reactions))

    def __getitem__(self, index: int) -> KReaction:
        """
        Returns the KReaction object at the specified index in the
        collection.

        Args:
            index (int): The index of the KReaction object to
                be returned.

        Returns:
            KReaction: The KReaction object at the specified index.
        """
        return self._reactions[index][1]

    @property
    def last(self) -> KReaction:
        """
        Returns the last KReaction instance in the collection.

        Returns:
            KReaction: The last KReaction instance.
        """
        return self._reactions[-1][1]

    def merge_with(
        self,
        other: tp.Self,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> None:
        """
        Merges the current KReactionCollection with another collection.

        Args:
            other (KReactionCollection): The other KReactionCollection
                to merge with.
            list_sp (Optional[ListSpecies]): The list of species to use
                for the merge. Defaults to None.

        Returns:
            None
        """
        for label, reaction in other.iter_labels_and_reactions():
            self.add_reaction(reaction, label=label, list_sp=list_sp)

    def update_label(self, current: str, new: str) -> tp.Self:
        """
        Updates the label of a reaction in the collection.

        Args:
            current (str): The current label of the reaction.
            new (str): The new label for the reaction.

        Returns:
            KReactionCollection: The updated KReactionCollection
                instance.

        Raises:
            KeyError: If the current label is not in the collection
                or the new label already exists.
        """
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
        Sorts the list of reactions by their label.

        Returns:
            KReactionCollection: The sorted KReactionCollection
                instance.
        """
        self._reactions = sorted(self._reactions, key=operator.itemgetter(0))
        return self

    def sort_by_reaction_index(self) -> tp.Self:
        """
        Sorts the list of reactions by their index.

        Returns:
            KReactionCollection: The sorted KReactionCollection
                instance.
        """
        self._reactions = sorted(self._reactions, key=lambda x: x[1].get_index())
        return self

    def iter_labels_and_reactions(self) -> tp.Iterator[tp.Tuple[str, KReaction]]:
        """
        Returns an iterator over the labels and reactions in the
        collection.

        Yields:
            Tuple[str, KReaction]: A tuple containing the label and
                the KReaction instance.
        """
        return iter(self._reactions)

    def iter_reactions_and_kvalues(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """
        Returns an iterator over the labels, reactions and their
        corresponding k values in the collection.

        Args:
            T: Temperature in °C (optional, default:25)

        Yields:
            Tuple[str, KReaction, float]: A tuple containing the label,
                the KReaction instance and its k value.
        """
        for label, reaction in self.iter_labels_and_reactions():
            yield (
                label,
                reaction,
                reaction.k_value(T),
            )

    def _get_reaction_index(self, label: str) -> int:
        """
        Returns the index of a reaction given its label.

        Args:
            label (str): The label of the reaction.

        Returns:
            int: The index of the reaction.
        """
        return self._indexer[label]

    def get_reaction(self, label: str) -> KReaction:
        """
        Retrieves a reaction from the collection by its label.

        Args:
            label (str): The label of the reaction to retrieve.

        Returns:
            KReaction: The reaction instance associated with the given
                label.

        Raises:
            KeyError: If the label is not found in the reaction
                collection.
        """
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
        Adds a new reaction to the collection.

        Args:
            reaction (GenericReaction | KReaction): The reaction to add.
            label (Optional[str]): The label for the reaction. If not
                provided, a unique ID will be generated.
            list_sp (Optional[ListSpecies]): The list of species to
                use for the reaction.

        Returns:
            Self: The updated collection instance.

        Raises:
            KeyError: If a reaction with the same label already
                exists in the collection.
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
    """
    Abstract base class representing a collection of 2 reactions.

    This class is the base class for equilibrium reactions. It provides
    abstract methods and properties that are common to the basic
    structure and methods for handling a couple of reactions.
    Subclasses should implement the specific behavior of the reaction.

    A BiReaction always defines a couple of forward and backward
    reactions.

    Attributes:
        _OPE_STRING (str): The string representing the operation of
            the reaction.
        _k_backward (Optional[np.polynomial.Polynomial]): The backward
            rate constant of the reaction.
        _k_forward (Optional[np.polynomial.Polynomial]): The forward
            rate constant of the reaction.
    """

    _OPE_STRING: str = "<->"
    _k_backward: tp.Optional[np.polynomial.Polynomial]
    _k_forward: tp.Optional[np.polynomial.Polynomial]

    def __contains__(self, species: object) -> bool:
        """
        Checks if a given species is involved in the reaction.

        Args:
            species (object): The species to be checked.

        Returns:
            bool: True if the species is in the reaction,
                False otherwise.
        """
        return species in self.forward or species in self.backward

    @property
    def forward(self) -> KReaction:
        """
        Property method to get the forward reaction of the BiReaction.

        Returns:
            KReaction: The forward reaction of the BiReaction.
        """
        return self.get_reaction("forward")

    @property
    def backward(self) -> KReaction:
        """
        Property method to get the backward reaction of the BiReaction.

        Returns:
            KReaction: The backward reaction of the BiReaction.
        """
        return self.get_reaction("backward")

    def add_reaction(
        self,
        reaction: GenericReaction | KReaction,
        label: tp.Optional[str] = None,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Add a reaction to the BiReaction.

        This method is abstract and should not be called. It always
        raises a NotImplementedError, because it is not possible to
        enlarge a BiReaction. A BiReaction always defines only a couple
        of forward and backward reactions.

        Args:
            reaction (GenericReaction | KReaction): The reaction to
                be added.
            label (str, optional): Label of the reaction.
                Defaults to None.
            list_sp (ListSpecies, optional): List of species.
                Defaults to None.

        Raises:
            NotImplementedError: Always raised.

        """
        raise NotImplementedError(
            "Not possible to add more reaction to an equilibrium collection"
        )

    def iter_reactions_and_kvalues(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """
        Returns an iterator over the labels, reactions and their
        corresponding k values in the collection.

        Args:
            T: Temperature in °C (optional, default:25)

        Yields:
            Tuple[str, KReaction, float]: A tuple containing the label,
                the KReaction instance and its k value.
        """
        yield "forward", self.get_reaction("forward"), self.k_forward(T)
        yield "backward", self.get_reaction("backward"), self.k_backward(T)

    def __str__(self) -> str:
        """
        Returns a string representation of the reaction.

        This method combines the reactants, products and the operation
        string to create a string that represents the reaction.

        Returns:
            str: A string representation of the reaction.
        """
        out = []
        out.append(self.forward.reactants_as_string())
        out.append(self._OPE_STRING)
        out.append(self.forward.products_as_string())
        return " ".join(out)

    def add_forward_reaction(
        self,
        reaction: GenericReaction | KReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Adds a reaction to the BiReaction as the forward reaction.

        Args:
            reaction (GenericReaction | KReaction): The reaction to be
                added as the forward reaction.
            list_sp (Optional[ListSpecies]): The list of species to use
                for the reaction. Defaults to None.

        Returns:
            Self: The updated BiReaction instance.
        """
        super().add_reaction(reaction, label="forward", list_sp=list_sp)
        self.forward._OPE_STRING = "<->"
        return self

    def add_backward_reaction(
        self,
        reaction: GenericReaction | KReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Adds a reaction to the BiReaction as the backward reaction.

        Args:
            reaction (GenericReaction | KReaction): The reaction to be
                added as the backward reaction.
            list_sp (Optional[ListSpecies]): The list of species to use
                for the reaction. Defaults to None.

        Returns:
            Self: The updated BiReaction instance.
        """
        super().add_reaction(reaction, label="backward", list_sp=list_sp)
        self.backward._OPE_STRING = "<->"
        return self

    @abstractmethod
    def k_equilibrium(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Abstract method to compute the equilibrium constant of the
        reaction.

        Args:
            T (Optional[Temperature | float]): The temperature in
                °C (optional, default:None)

        Returns:
            float: The equilibrium constant at the given temperature.
        """

    @abstractmethod
    def k_forward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Abstract method to compute the forward reaction constant of
        the reaction.

        Args:
            T (Optional[Temperature | float]): The temperature in
                °C (optional, default:None)

        Returns:
            float: The forward reaction constant at the given
                temperature.
        """

    @abstractmethod
    def k_backward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Abstract method to compute the backward reaction constant of
        the reaction.

        Args:
            T (Optional[Temperature | float]): The temperature in
                °C (optional, default:None)

        Returns:
            float: The backward reaction constant at the given
                temperature.
        """
