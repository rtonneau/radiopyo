#!/usr/bin/env python3
"""
The k_reactions submodule contains all class definitions related to basic chemical
reactions defined via a reaction constant.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import operator
import typing as tp
import uuid
from abc import ABC, abstractmethod

import numpy as np
from scipy.constants import R as GAS_CONSTANT  # type: ignore[import-untyped]

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.basic_reaction import BasicReaction
from radiopyo.chemistry.generics import revert_reaction
from radiopyo.chemistry.species_collection import ListSpecies
from radiopyo.chemistry.utils.factories import polynomial_factory
from radiopyo.chemistry.utils.units import Temperature

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.generics import GenericReaction

# --------------------------------------------------------------------------- #

__all__ = [
    "KReaction",
    "KReactionCollection",
    "KValType",
]


# --------------------------------------------------------------------------- #
#                              FUNCTION DEFINITIONS                           #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #

KValType: tp.TypeAlias = tp.Union[float, tp.List[float], np.polynomial.Polynomial]


class KReaction(BasicReaction):
    """
    Almost all the logic is already included in the BasicReaction
    class. Add the reaction rate constant as attribute. It is stored as
    a polynomial of log10(k) in 1/T [1/K]:
        log10(k) = a0 + a1 (1/T) + a2 (1/T)**2 + ...
    so a single value k is stored as [log10(k)].

    n·k convention (ADR 0001):
        For reactions with a single reactant species, e.g.
            2 H_r -> H2
        2k is usually given/measured instead of k (reaction constant).
        When such a reaction is built from a config ('from_generic'),
        the configured value is divided by the stoichiometry n of the
        reactant ONLY IF k is a single number (scalar k). Polynomial
        k(T), Arrhenius and equilibrium (acid/base, EReaction) k are used
        as given. 'raw_k_value()' gives back the configured (n·k) value.
        The constructor and 'set_k_value' never divide.

    Attributes:
        _k_value: np.polynomial.Polynomial | None
            log10 of the reaction rate constant, polynomial of 1/T[K].
    """

    _k_value: tp.Optional[np.polynomial.Polynomial]
    _OPE_STRING: str = "->"

    def __init__(
        self,
        k_value: tp.Optional[KValType] = None,
    ):
        super().__init__()
        self._k_value = (
            polynomial_factory(k_value, log10=True) if k_value is not None else None
        )

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
        _T = Temperature(T)
        return np.power(10, self._k_value(1.0 / _T.kelvin))

    @property
    def k_polynomial(self) -> tp.Optional[np.polynomial.Polynomial]:
        """
        The stored log10(k) polynomial (of 1/T [1/K]), None if no rate
        constant is defined.
        """
        return self._k_value

    def set_k_value(self, value: KValType) -> None:
        """
        Set the reaction constant value. A float is a k value, a list
        of floats or a polynomial are log10(k) coefficients (see class
        docstring). No n·k division is applied here.

        Args:
            value (KValType): The new reaction constant value.

        Returns:
            None
        """
        self._k_value = polynomial_factory(value, log10=True)

    def _apply_nk_convention(self) -> None:
        """
        n·k convention (ADR 0001): for n A -> ..., the configured scalar
        k is n·k, divide it by n. Polynomial k are kept as given.
        """
        if (
            self._k_value is not None
            and len(self.reactants) == 1
            and len(self._k_value) == 1
        ):
            reactant = next(iter(self.reactants))
            self.set_k_value(self.k_value() / reactant.stoi)

    def to_dict(self) -> tp.Dict[str, tp.Any]:
        """
        Convert the KReaction to a dict structure
        {
            "reaction": str,
            "k_value": float | List[float]
        }
        """
        return {
            "reaction": self.__str__(),
            "k_value": self.raw_k_value(),
        }

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
        return self._index_label("KR")

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
            return tp.cast("list[float]", self._k_value.coef)  # ndarray at runtime

        # n·k convention (ADR 0001): give back the configured n·k value
        if len(self.reactants) == 1:
            sp = next(iter(self.reactants))
            return self.k_value() * sp.stoi
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
        reaction = cls(generic["constants"].get("k_value", None))
        reaction._add_species_from_generic(generic, list_sp)
        reaction._apply_nk_convention()
        return reaction

    def __str__(self) -> str:
        """
        Returns a string representation of the reaction.

        This method combines the reactants, products and the operation
        string to create a string that represents the reaction.

        Returns:
            str: A string representation of the reaction.
        """
        out = []
        out.append(self.reactants.as_string())
        out.append(self._OPE_STRING)
        out.append(self.products.as_string())
        return " ".join(out)

    def __repr__(self) -> str:
        return f"KReaction('{self.__str__()}', k_value={self.raw_k_value()})"


class Arrhenius(KReaction):
    """
    Class representing an Arrhenius reaction:
        k(T) = A * exp(-Ea / (R * T))
    with T the simulation temperature in K and R the gas constant
    (8.314 J/mol/K). Ea is given in kJ/mol (converted to J/mol here).

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
    ):
        """ """
        super().__init__(k_value=None)
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
        return float(
            self._pre_exp
            * np.exp(-self._activation_energy * 1e3 / (GAS_CONSTANT * T.kelvin))
        )

    def raw_k_value(self) -> tp.List[float]:
        """
        Produce the k_value as provided by config files.

        Returns:
            List[float]: A list containing the pre-exponential factor
            and the activation energy.
        """
        return [self._pre_exp, self._activation_energy]

    def to_dict(self) -> tp.Dict[str, tp.Any]:
        """
        Convert the KReaction to a dict structure
        {
            "reaction": str,
            "k_value": float | List[float]
        }
        """
        return {
            "reaction": self.__str__(),
            "A": self._pre_exp,
            "Ea": self._activation_energy,
        }

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
        reaction = cls(
            generic["constants"]["A"],  # type: ignore[arg-type]
            generic["constants"]["Ea"],  # type: ignore[arg-type]
        )
        reaction._add_species_from_generic(generic, list_sp)
        return reaction


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
        return self._index_label("CR")

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
        reaction = cls(generic["constants"].get("k_value", None))
        reaction._add_species_from_generic(generic, list_sp, reactants=False)
        return reaction

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
        out.append(self.products.as_string())
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
            _reactions list. It is rebuilt whenever _reactions is
            reordered (see '_reindex').
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

    def _reindex(self) -> None:
        """Rebuild the label -> position mapping after a reordering."""
        self._indexer = {label: idx for idx, (label, _) in enumerate(self._reactions)}

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
        Updates the label of a reaction in the collection. The reaction
        keeps its position in the collection.

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

        self._reactions[idx] = (new, self._reactions[idx][1])
        del self._indexer[current]
        self._indexer[new] = idx
        return self

    def sort_by_label(self) -> tp.Self:
        """
        Sorts the list of reactions by their label.

        Returns:
            KReactionCollection: The sorted KReactionCollection
                instance.
        """
        self._reactions = sorted(self._reactions, key=operator.itemgetter(0))
        self._reindex()
        return self

    def sort_by_reaction_index(self) -> tp.Self:
        """
        Sorts the list of reactions by their index.

        Returns:
            KReactionCollection: The sorted KReactionCollection
                instance.
        """
        self._reactions = sorted(self._reactions, key=lambda x: x[1].get_index())
        self._reindex()
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
    reactions. Both are KReactions owning their rate constant: the
    forward/backward k are always read from them (no copy is kept by
    the BiReaction).

    Attributes:
        _OPE_STRING (str): The string representing the operation of
            the reaction.
    """

    _OPE_STRING: str = "<->"

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
        yield "forward", self.forward, self.k_forward(T)
        yield "backward", self.backward, self.k_backward(T)

    def __str__(self) -> str:
        """
        Returns a string representation of the reaction.

        This method combines the reactants, products and the operation
        string to create a string that represents the reaction.

        Returns:
            str: A string representation of the reaction.
        """
        out = []
        out.append(self.forward.reactants.as_string())
        out.append(self._OPE_STRING)
        out.append(self.forward.products.as_string())
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

    def _add_reaction_pair(
        self,
        reaction: GenericReaction,
        k_forward: KValType,
        k_backward: KValType,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> None:
        """
        Build the forward reaction from 'reaction' and the backward one
        from its reverse, and set their rate constants (log10
        polynomials, see KReaction). No n·k division is applied.

        Args:
            reaction (GenericReaction): The forward reaction.
            k_forward (KValType): The forward rate constant.
            k_backward (KValType): The backward rate constant.
            list_sp (Optional[ListSpecies]): The list of species to use
                for the reactions. Defaults to None.
        """
        self.add_forward_reaction(reaction, list_sp=list_sp)
        self.forward.set_k_value(k_forward)
        self.add_backward_reaction(revert_reaction(reaction), list_sp=list_sp)
        self.backward.set_k_value(k_backward)

    @abstractmethod
    def to_dict(self) -> tp.Dict[str, tp.Any]:
        """
        Convert the reaction to a dictionary.

        Returns:
            dict    A dictionary representation of the reaction.
        """

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

    def k_forward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Compute the forward reaction constant of the reaction at a given
        temperature (read from the forward KReaction).

        Args:
            T (Optional[Temperature | float]): The temperature in
                °C (optional, default:None)

        Returns:
            float: The forward reaction constant at the given
                temperature.
        """
        return self.forward.k_value(T)

    def k_backward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Compute the backward reaction constant of the reaction at a
        given temperature (read from the backward KReaction).

        Args:
            T (Optional[Temperature | float]): The temperature in
                °C (optional, default:None)

        Returns:
            float: The backward reaction constant at the given
                temperature.
        """
        return self.backward.k_value(T)
