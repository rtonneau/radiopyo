#!/usr/bin/env python3
"""
The base_chemistry submodule contains all class definitions and protocols to build some
more advanced chemistry ==> building blocks of basic chemistry.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from functools import reduce
from itertools import chain
from operator import add
from typing import (
    Dict,
    Iterable,
    Optional,
    Protocol,
    Tuple,
)

# --------------------------------- Internal imports --------------------------------- #
from .species import RawSpecies, SimSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
# ------------------------------------------------------------------------------------ #
__all__ = [
    "BaseChemicalReaction",
    "IsChemicalReaction",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class IsChemicalReaction(Protocol):
    """Protocol => used for type hints"""
    reactants: Dict[str, SimSpecies]
    products: Dict[str, SimSpecies]
    stoi_reactants: Dict[str, int]
    stoi_products: Dict[str, int]
    index: Optional[int]

    def __contains__(self, other: object) -> bool:
        """"""

    def iter_species(self) -> Iterable[SimSpecies]:
        """Iter over species involved in the reaction"""

    def iter_reactants(self) -> Iterable[Tuple[SimSpecies, int]]:
        """ Iter over reactants only"""

    def iter_products(self) -> Iterable[Tuple[SimSpecies, int]]:
        """Iter over products only"""

    def as_label(self) -> str:
        """ """

    def has_product(self, product: SimSpecies | str) -> int:
        """"""

    def has_reactant(self, product: SimSpecies | str) -> int:
        """"""

    def make_hash(self, method: str) -> str:
        """"""


class BaseChemicalReaction():
    """
    Class for very basic chemical reaction. The most simple reaction is
    2A + B -> C + D
    where
        A and B are defined as reactants
        C and D are defined as products
    It is important to also keep track of the stoichiometry of each species.
    Chemical Reactions should have a unique index to easily refer to it.
    Many useful methods can be defined based on this simple implementation.

    Attributes:
        reactant: Dict[str, SimSpecies]
        products: Dict[str, SimSpecies]
        stoi_reactants: Dict[str, int]
        stoi_products: Dict[str, int]
        index: int
    """
    reactants: Dict[str, SimSpecies]
    products: Dict[str, SimSpecies]
    stoi_reactants: Dict[str, int]
    stoi_products: Dict[str, int]
    index: Optional[int]

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

    def make_hash(self, method: str = "sha256") -> str:
        """
        """
        raise NotImplementedError

    def tot_number_reactants(self) -> int:
        """
        Count the total number of reactants involved in the reaction
        => Stoichiometry dependent! In other words, it counts the number of
        molecules/entities defined as reactants.

        Returns:
            int
        """
        return reduce(add, list(self.stoi_reactants.values()))

    def iter_species(self) -> Iterable[SimSpecies]:
        """Iter over ALL the species involved in the reaction

        Returns:
            Iterator: SimSpecies
        """
        return chain(self.reactants.values(), self.products.values())

    def iter_reactants(self) -> Iterable[Tuple[SimSpecies, int]]:
        """ Iter over reactants only, including the stoichiometry.

        Returns:
            Iterator: Tuple[SimSpecies, int]
                The second element is the stoichiometry of the reactant.
        """
        return [(sp, self.stoi_reactants[label])
                for label, sp in self.reactants.items()]

    def iter_products(self) -> Iterable[Tuple[SimSpecies, int]]:
        """Iter over products only, including the stoichiometry.
        Returns:
            Iterator: Tuple[SimSpecies, int]
                The second element is the stoichiometry of the product.
        """
        return [(sp, self.stoi_products[label])
                for label, sp in self.products.items()]

    def set_index(self, index: int) -> None:
        """ Set the index of the reaction.

        Args:
            index: int

        Returns:
            None

        """
        self.index = index

    def number_of_reactants(self) -> int:
        """Count the number of reactants involved, independently of the stoichiometry.
        Returns:
            int
        """
        return len(list(self.iter_reactants()))

    def number_of_products(self) -> int:
        """Count the number of reactants involved, independently of the stoichiometry.
        Returns:
            int
        """
        return len(list(self.iter_products()))

    def has_product(self, product: SimSpecies | str) -> bool:
        """Check if the reaction involves a specific product.

        Args:
            product: str | SimSpecies

        Returns:
            bool

        """
        if isinstance(product, str):
            return self.stoi_products.get(product, 0) > 0
        return self.stoi_products.get(product.label, 0) > 0

    def has_reactant(self, reactant: SimSpecies | str) -> int:
        """Check if the reaction involves a specific reactant.

        Args:
            reactant: str | SimSpecies

        Returns:
            bool

        """
        if isinstance(reactant, str):
            return self.stoi_reactants.get(reactant, 0) > 0
        return self.stoi_reactants.get(reactant.label, 0) > 0

    def add_as_product(self, sp: SimSpecies, multiple: int = 1) -> None:
        """
        Method to conveniently add a new product to the reaction by keeping the
        stoichiometry relevant:
            -> If the species already exists as product, increment its stoichiometry
            -> If the species does not already exists as product, add it to the dict of
            products and set its stoichiometry to 1

        Args:
            sp: SimSpecies

        Returns:
            None
        """
        if sp.label in self.products:
            self.stoi_products[sp.label] += multiple
        else:
            self.products[sp.label] = sp
            self.stoi_products[sp.label] = multiple

    def add_as_reactant(self, sp: SimSpecies, multiple: int = 1) -> None:
        """
        Method to conveniently add a new reactant to the reaction by keeping the
        stoichiometry relevant:
            -> If the species already exists as reactant, increment its stoichiometry
            -> If the species does not already exists as reactant, add it to the dict of
            reactants and set its stoichiometry to 1

        Args:
            sp: SimSpecies

        Returns:
            None
        """
        if sp.label in self.reactants:
            self.stoi_reactants[sp.label] += multiple
        else:
            self.reactants[sp.label] = sp
            self.stoi_reactants[sp.label] = multiple

    def as_label(self) -> str:
        """ 
        This method is used for the data mining. It should be implemented by child
        classes in order to easily define short labels for each reactions based on the
        reaction type and the index of reaction. For instance, it can be used as column 
        name in a pandas DataFrame.

        Returns:
            str
        """
        raise NotImplementedError
