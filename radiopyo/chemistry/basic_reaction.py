#!/usr/bin/env python3
"""
The base_chemistry submodule contains all class definitions and
protocols to build some more advanced chemistry ==> building blocks of
basic chemistry.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp
from abc import ABC, abstractmethod
from contextlib import suppress
from dataclasses import dataclass
from functools import reduce
from itertools import chain
from operator import add

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser import RPARSER

from .generics import GenericReaction, Operation, dict_to_generic_reaction
from .species import BasicSpecies, ChemicalSpecies, SimSpecies

if tp.TYPE_CHECKING:
    from .species_collection import ListSpecies

# ------------------------- Type checking ONLY imports----------------------- #
# --------------------------------------------------------------------------- #
__all__: tp.List[str] = [
    "BasicReaction",
    "ReactionSpecies",
    "ReactionSpeciesCollection",
]


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ReactionSpecies:
    """ """

    species: SimSpecies
    stoi: int
    ope: Operation

    def __init__(
        self,
        species: SimSpecies,
        stoi: tp.Optional[int] = None,
        ope: tp.Optional[Operation] = None,
    ):
        object.__setattr__(self, "species", species)
        object.__setattr__(self, "stoi", stoi if stoi is not None else 1)
        object.__setattr__(self, "ope", ope if ope is not None else Operation.ADD)

    @property
    def label(self) -> str:
        return self.species.label


@dataclass(init=False)
class ReactionSpeciesCollection:
    """"""

    species_col: tp.Dict[str, ReactionSpecies]
    _order: tp.List[str]

    def __init__(self) -> None:
        self.species_col = {}
        self._order = []

    def __len__(self) -> int:
        return len(self.species_col)

    def __iter__(self) -> tp.Iterator[ReactionSpecies]:
        """ """
        return iter((self.species_col[sp] for sp in self._order))

    def __contains__(
        self, key: str | BasicSpecies | SimSpecies | ReactionSpecies
    ) -> bool:
        if not isinstance(key, str):
            key = key.label
        return key in self.species_col

    def __getitem__(self, key: str | BasicSpecies | SimSpecies) -> ReactionSpecies:
        if not isinstance(key, str):
            key = key.label
        return self.species_col[key]

    def items(self) -> tp.Generator[tp.Tuple[str, ReactionSpecies], None, None]:
        yield from ((label, self.species_col[label]) for label in self._order)

    def iter_stoi(self) -> tp.Iterable[tp.Tuple[ReactionSpecies, int]]:
        for sp in self.__iter__():
            yield sp, sp.stoi

    def _add_species(
        self,
        species: SimSpecies,
        stoi: tp.Optional[int] = None,
        ope: tp.Optional[Operation] = None,
    ) -> None:
        if species.label in self._order:
            raise KeyError(f"{species.label} is already recorded in {self}")
        self.species_col[species.label] = ReactionSpecies(species, stoi, ope)
        self._order.append(species.label)

    def add_species(
        self,
        species: SimSpecies,
        stoi: tp.Optional[int] = None,
        ope: tp.Optional[Operation] = None,
    ) -> None:
        """ """
        self._add_species(species, stoi, ope)

    def stoichiometry_of(self, key: str | BasicSpecies | SimSpecies) -> int:
        if not isinstance(key, str):
            key = key.label
        return self.species_col[key].stoi

    def tot_number(self) -> int:
        """
        Count the total number of elements involved in the reaction
        => Stoichiometry dependent! In other words, it counts the number
        of molecules/entities.

        Returns:
            int
        """
        return 0 if len(self) == 0 else reduce(add, [sp.stoi for sp in self.__iter__()])

    def as_string(self) -> str:
        """
        Method to generate a string representation of the collection in
        the reaction.

        Returns:
            str
                A string representation of the species collection,
                with stoichiometry.
        """
        out: tp.List[str] = []
        for sp in self.__iter__():
            _ = "" if sp.stoi <= 1 else f"{sp.stoi}"
            if len(out) > 0 or sp.ope is Operation.SUBTRACT:
                out.append(str(sp.ope))
            out.append(_ + f"{sp.label}")
        return " ".join(out)


class BasicReaction(ABC):
    """
    Abstract base class for a basic chemical reaction.

    This class provides a foundation for representing a chemical
    reaction, including the reactants, products, stoichiometric
    coefficients, and reaction index. Represents very basic chemical
    reaction. The most simple reaction is 2A + B -> C + D where
        A and B are defined as reactants
        C and D are defined as products
    It is important to also keep track of the stoichiometry of each
    species. Chemical Reactions should have a unique index to easily
    refer to it.

    Attributes:
    reactants : tp.Dict[str, ChemSpecies]
        A dictionary of reactants in the reaction, where each key is a
        reactant name and each value is a ChemSpecies object
        representing the reactant.
    products : tp.Dict[str, ChemSpecies]
        A dictionary of products in the reaction, where each key is a
        product name and each value is a ChemSpecies object representing
        the product.
    stoi_reactants : tp.Dict[str, int]
        A dictionary of stoichiometric coefficients for reactants, where
        each key is a reactant name and each value is the corresponding
        stoichiometric coefficient.
    stoi_products : tp.Dict[str, int]
        A dictionary of stoichiometric coefficients for products, where
        each key is a product name and each value is the corresponding
        stoichiometric coefficient.
    index : tp.Optional[int]
        An optional index for the reaction, which can be used to
        identify the reaction in a collection of reactions.

    Notes:
        This class is intended to be subclassed by more specific
        reaction classes, such as elementary reactions or complex
        reactions.
    """

    reactants: ReactionSpeciesCollection
    products: ReactionSpeciesCollection
    index: tp.Optional[int]

    def __init__(self) -> None:
        self.reactants = ReactionSpeciesCollection()
        self.products = ReactionSpeciesCollection()
        self.index = None

    def __contains__(self, other: object) -> bool:
        """
        Check where the reaction involves a specific species
        Args:
            other: str | BasicSpecies

        Returns:
            bool

        Raises:
            NotImplementedError:
                if 'other' have any other type
        """
        if isinstance(other, (str, BasicSpecies)):
            return (other in self.reactants) or (other in self.products)
        raise NotImplementedError

    def __eq__(self, other: tp.Any) -> bool:
        """
        Check if two reactions are equal.

        Args:
            other: The other reaction to compare.

        Returns:
            bool: True if the reactions are equal, False otherwise.
        """
        if isinstance(other, str):
            r = RPARSER.parse_reaction(other)
        elif isinstance(other, dict):
            r = dict_to_generic_reaction(other)
        elif isinstance(other, BasicReaction):
            r = other.to_generic()
        else:
            raise TypeError(f"Cannot compare {self} with {type(other)}")

        for label, sp in r["reactants"].items():
            if (
                label not in self.reactants
                or sp[1] != self.reactants[label].stoi
                or sp[2] != self.reactants[label].ope
            ):
                return False
        for label, sp in r["products"].items():
            if (
                label not in self.products
                or sp[1] != self.products[label].stoi
                or sp[2] != self.products[label].ope
            ):
                return False
        return len(r["reactants"]) + len(r["products"]) == len(self.reactants) + len(
            self.products
        )

    @abstractmethod
    def to_dict(self) -> tp.Dict[str, tp.Any]:
        """
        Convert the reaction to a dictionary.

        Returns:
            dict    A dictionary representation of the reaction.
        """

    @abstractmethod
    def make_hash(self, method: str = "sha256") -> str:
        """
        Produce a unique identifier (content based) for easy comparison.

        Args:
            method (str): Hash method to use (default: "sha256").

        Returns:
            str: A hash value representing the object.
        """

    @abstractmethod
    def as_label(self) -> str:
        """
        This method is used for the data mining. It should be
        implemented by child classes in order to easily define short
        labels for each reactions based on the reaction type and the
        index of reaction. For instance, it can be used as column name
        in a pandas DataFrame.

        Returns:
            str
        """

    @classmethod
    @abstractmethod
    def from_generic(
        cls,
        generic: GenericReaction,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Create a Reaction from a 'GenericReaction', Typed-Dict
        structure produced by radiopyo module parser.

        Args:
            generic (GenericReaction): Typed-Dict structure holding the
                reaction metadata.
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.

        Returns:
            Self: An instance of Reaction.
        """

    def to_generic(self) -> GenericReaction:
        """
        Convert the reaction to a GenericReaction.

        Returns:
            GenericReaction: A GenericReaction representation of the
                reaction.
        """
        return GenericReaction(
            products={
                label: (sp.species.raw_species, sp.stoi, sp.ope)
                for label, sp in self.products.items()
            },
            reactants={
                label: (sp.species.raw_species, sp.stoi, sp.ope)
                for label, sp in self.reactants.items()
            },
            constants={},
            enzyme=None,
        )

    def get_stoi_of(self, species: SimSpecies | ChemicalSpecies | str) -> int:
        """
        Method to extract the stoichiometric coefficient of a species
        regardless if the species is a reactant or a product.

        Arguments:
            species: ChemSpecies | str

        Returns:
            float the stoichiometric coefficient of species

        Raises:
            KeyError:
                If species not part of the reaction.
        """
        with suppress(KeyError):
            return self.reactants.stoichiometry_of(species)
        with suppress(KeyError):
            return self.products.stoichiometry_of(species)

        raise KeyError(f"Species {species} not in reaction {self}")

    def tot_number_reactants(self) -> int:
        """
        Count the total number of reactants involved in the reaction
        => Stoichiometry dependent! In other words, it counts the number
        of molecules/entities defined as reactants.

        Returns:
            int
        """
        return self.reactants.tot_number()

    def tot_number_products(self) -> int:
        """
        Count the total number of reactants involved in the reaction
        => Stoichiometry dependent! In other words, it counts the number
        of molecules/entities defined as reactants.

        Returns:
            int
        """
        return self.products.tot_number()

    def iter_species(self) -> tp.Iterable[ReactionSpecies]:
        """
        Iterator over all the species involved in the reaction,
        including both reactants and products.

        Returns:
            Iterable[ChemSpecies]: An iterator over all the species
                involved in the reaction.
        """
        return chain(iter(self.reactants), iter(self.products))

    def get_index(self, default: int = -1) -> int:
        """
        Returns the index of the reaction if it exists, otherwise
        returns the default value.

        Args:
            default (int): The default value to return if the index does
            not exist. Defaults to -1.

        Returns:
            int: The index of the reaction or the default value.
        """
        return self.index if self.index is not None else default

    def set_index(self, index: int) -> None:
        """
        Sets the index of the reaction.

        Args:
            index (int): The index to be set.

        Returns:
            None
        """
        self.index = index

    def number_of_reactants(self) -> int:
        """
        Returns the total number of reactants involved in the reaction.

        Returns:
            int
        """
        return len(self.reactants)

    def number_of_products(self) -> int:
        """
        Returns the number of products involved in the reaction.

        Returns:
            int
                The number of products.
        """
        return len(self.products)

    def has_product(self, product: SimSpecies | ChemicalSpecies | str) -> bool:
        """
        Checks if the reaction involves a specific product.

        Args:
            product (ChemSpecies | str): The product to be checked. Can
            be either a ChemSpecies object or a string.

        Returns:
            bool: True if the reaction involves the product,
            False otherwise.
        """
        return product in self.products

    def has_reactant(self, reactant: SimSpecies | ChemicalSpecies | str) -> int:
        """
        Checks if a given reactant is involved in the reaction.

        Args:
            reactant (ChemSpecies | str): The reactant to be checked.
            Can be either a ChemSpecies object or a string.

        Returns:
            int: 1 if the reactant is part of the reaction, 0 otherwise.
        """
        return reactant in self.reactants

    def add_as_product(
        self,
        species: SimSpecies,
        stoi: tp.Optional[int] = None,
        ope: tp.Optional[Operation] = None,
    ) -> None:
        """
        Method to add a new product to the reaction by keeping the
        stoichiometry relevant.
            -> If the species already exists as product,
            increment its stoichiometry.
            -> If the species does not already exists as product,
            add it to the dict of products and set its stoichiometry
            to 1.

        Args:
            species: SimSpecies
            stoi (int, optional): The stoichiometric coefficient of
                the reactant.
            ope (Operation, optional): The operation to apply.

        Returns:
            None
        """
        self.products.add_species(species, stoi, ope)

    def add_as_reactant(
        self,
        species: SimSpecies,
        stoi: tp.Optional[int] = None,
        ope: tp.Optional[Operation] = None,
    ) -> None:
        """
        Method to conveniently add a new reactant to the reaction by
        keeping the stoichiometry relevant:
            -> If the species already exists as reactant,
            increment its stoichiometry
            -> If the species does not already exists as reactant,
            add it to the dict of reactants and set its stoichiometry
            to 1

        Args:
            species: SimSpecies
            stoi (int, optional): The stoichiometric coefficient of
                the reactant.
            ope (Operation, optional): The operation to apply.

        Returns:
            None
        """
        self.reactants.add_species(species, stoi, ope)
