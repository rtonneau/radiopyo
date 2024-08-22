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
from functools import reduce
from itertools import chain
from operator import add

# ----------------------------- Internal imports ---------------------------- #
from .species import ChemSpecies, RawSpecies

if tp.TYPE_CHECKING:
    from .collections import ListSpecies

# ------------------------- Type checking ONLY imports----------------------- #
# --------------------------------------------------------------------------- #
__all__ = [
    "BasicReaction",
    "IsChemicalReaction",
]

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class IsChemicalReaction(tp.Protocol):
    """
    Protocol for a chemical reaction.

    A class implementing this protocol should provide the following
    attributes:

    Attributes:
    reactants: tp.Dict[str, ChemSpecies]
        A dictionary of reactants in the reaction.
    products: tp.Dict[str, ChemSpecies]
        A dictionary of products in the reaction.
    stoi_reactants: tp.Dict[str, int]
        A dictionary of stoichiometric coefficients for reactants.
    stoi_products: tp.Dict[str, int]
        A dictionary of stoichiometric coefficients for products.
    index: tp.Optional[int]
        An optional index for the reaction.

    Methods:
    __contains__(other: object) -> bool
        Checks if a given object is in the reaction.

    Notes:
        This protocol is used to define the interface for a chemical
        reaction. Any class implementing this protocol should provide
        the required attributes and methods.
    """

    reactants: tp.Dict[str, ChemSpecies]
    products: tp.Dict[str, ChemSpecies]
    stoi_reactants: tp.Dict[str, int]
    stoi_products: tp.Dict[str, int]
    index: tp.Optional[int]

    def __contains__(self, other: object) -> bool:
        """"""

    def iter_species(self) -> tp.Iterable[ChemSpecies]:
        """Iter over species involved in the reaction"""

    def iter_reactants(self) -> tp.Iterable[tp.Tuple[ChemSpecies, int]]:
        """Iter over reactants only"""

    def iter_products(self) -> tp.Iterable[tp.Tuple[ChemSpecies, int]]:
        """Iter over products only"""

    def as_label(self) -> str:
        """ """

    def has_product(self, product: ChemSpecies | str) -> int:
        """"""

    def has_reactant(self, product: ChemSpecies | str) -> int:
        """"""

    def make_hash(self, method: str) -> str:
        """"""


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

    reactants: tp.Dict[str, ChemSpecies]
    products: tp.Dict[str, ChemSpecies]
    stoi_reactants: tp.Dict[str, int]
    stoi_products: tp.Dict[str, int]
    index: tp.Optional[int]

    def __init__(self) -> None:
        self.reactants = {}
        self.products = {}
        self.stoi_products = {}
        self.stoi_reactants = {}
        self.index = None

    def __contains__(self, other: object) -> bool:
        """
        Check where the reaction involves a specific species
        Args:
            other: str | RawSpecies

        Returns:
            bool

        Raises:
            NotImplementedError:
                if 'other' have any other type
        """
        if isinstance(other, (str, RawSpecies)):
            return len([elt for elt in self.iter_species() if elt == other]) > 0
        raise NotImplementedError

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
    def from_kwargs(
        cls,
        list_sp: tp.Optional[ListSpecies] = None,
        **kwargs: str | tp.List[str] | float | tp.List[float],
    ) -> tp.Self:
        """
        Abstract class method to create a new instance from keyword
        arguments.

        Args:
            list_sp (ListSpecies): Optional list of species.
                Defaults to None.
            **kwargs: Keyword arguments to create the instance.

        Returns:
            tp.Self: A new instance of the class.
        """

    def get_stoi_of(self, species: ChemSpecies | str) -> float:
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
        if isinstance(species, ChemSpecies):
            species = species.label
        if species in self.stoi_products:
            return self.stoi_products[species]
        if species in self.stoi_reactants:
            return self.stoi_reactants[species]
        raise KeyError(f"Species {species} not in reaction {self}")

    def tot_number_reactants(self) -> int:
        """
        Count the total number of reactants involved in the reaction
        => Stoichiometry dependent! In other words, it counts the number
        of molecules/entities defined as reactants.

        Returns:
            int
        """
        return (
            0
            if len(self.reactants) == 0
            else reduce(add, list(self.stoi_reactants.values()))
        )

    def tot_number_products(self) -> int:
        """
        Count the total number of reactants involved in the reaction
        => Stoichiometry dependent! In other words, it counts the number
        of molecules/entities defined as reactants.

        Returns:
            int
        """
        return (
            0
            if len(self.products) == 0
            else reduce(add, list(self.stoi_products.values()))
        )

    def iter_species(self) -> tp.Iterable[ChemSpecies]:
        """
        Iterator over all the species involved in the reaction,
        including both reactants and products.

        Returns:
            Iterable[ChemSpecies]: An iterator over all the species
                involved in the reaction.
        """
        return chain(self.reactants.values(), self.products.values())

    def iter_reactants(self) -> tp.Iterable[tp.Tuple[ChemSpecies, int]]:
        """
        Iter over reactants only, including the stoichiometry.

        Returns:
            Iterator: Tuple[SimSpecies, int]
                The second element is the stoichiometry of the reactant.
        """
        return [
            (sp, self.stoi_reactants[label]) for label, sp in self.reactants.items()
        ]

    def iter_products(self) -> tp.Iterable[tp.Tuple[ChemSpecies, int]]:
        """
        Iterator over products only, including the stoichiometry.

        Returns:
            Iterator: Tuple[ChemSpecies, int]
                The second element is the stoichiometry of the product.
        """
        return [(sp, self.stoi_products[label]) for label, sp in self.products.items()]

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
        return len(list(self.iter_reactants()))

    def number_of_products(self) -> int:
        """
        Returns the number of products involved in the reaction.

        Returns:
            int
                The number of products.
        """
        return len(list(self.iter_products()))

    def has_product(self, product: ChemSpecies | str) -> bool:
        """
        Checks if the reaction involves a specific product.

        Args:
            product (ChemSpecies | str): The product to be checked. Can
            be either a ChemSpecies object or a string.

        Returns:
            bool: True if the reaction involves the product,
            False otherwise.
        """
        if isinstance(product, str):
            return self.stoi_products.get(product, 0) > 0
        return self.stoi_products.get(product.label, 0) > 0

    def has_reactant(self, reactant: ChemSpecies | str) -> int:
        """
        Checks if a given reactant is involved in the reaction.

        Args:
            reactant (ChemSpecies | str): The reactant to be checked.
            Can be either a ChemSpecies object or a string.

        Returns:
            int: 1 if the reactant is part of the reaction, 0 otherwise.
        """
        if isinstance(reactant, str):
            return self.stoi_reactants.get(reactant, 0) > 0
        return self.stoi_reactants.get(reactant.label, 0) > 0

    def add_as_product(self, sp: ChemSpecies, multiple: int = 1) -> None:
        """
        Method to add a new product to the reaction by keeping the
        stoichiometry relevant.
            -> If the species already exists as product,
            increment its stoichiometry.
            -> If the species does not already exists as product,
            add it to the dict of products and set its stoichiometry
            to 1.

        Args:
            sp (ChemSpecies): The product to be added.
            multiple (int, optional): The stoichiometric coefficient of
                the product. Defaults to 1.

        Returns:
            None
        """
        if sp.label in self.products:
            self.stoi_products[sp.label] += multiple
        else:
            self.products[sp.label] = sp
            self.stoi_products[sp.label] = multiple

    def add_as_reactant(self, sp: ChemSpecies, multiple: int = 1) -> None:
        """
        Method to conveniently add a new reactant to the reaction by
        keeping the stoichiometry relevant:
            -> If the species already exists as reactant,
            increment its stoichiometry
            -> If the species does not already exists as reactant,
            add it to the dict of reactants and set its stoichiometry
            to 1

        Args:
            sp: SimSpecies
            multiple (int, optional): The stoichiometric coefficient of
                the reactant. Defaults to 1.

        Returns:
            None
        """
        if sp.label in self.reactants:
            self.stoi_reactants[sp.label] += multiple
        else:
            self.reactants[sp.label] = sp
            self.stoi_reactants[sp.label] = multiple

    def reactants_as_string(self) -> str:
        """
        Method to generate a string representation of the reactants in
        the reaction.
        Returns:

            str
                A string representation of the reactants,
                with stoichiometry.
        """
        out = []
        for sp, stoi in self.iter_reactants():
            _ = "" if stoi <= 1 else f"{stoi}"
            out.append(_ + f"{sp.label}")
            out.append("+")
        return " ".join(out[:-1])

    def products_as_string(self) -> str:
        """
        Method to generate a string representation of the products
        in the reaction.
        Returns:
            str
                A string representation of the products,
                with stoichiometry.
        """
        out = []
        for sp, stoi in self.iter_products():
            _ = "" if stoi <= 1 else f"{stoi}"
            out.append(_ + f"{sp.label}")
            out.append("+")
        return " ".join(out[:-1])
