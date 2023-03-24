from __future__ import annotations

from functools import reduce
from itertools import chain
from operator import add
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterable,
    List,
    Protocol,
    Tuple,
)

import numpy as np

import radiopyo_vec.chemistry.exceptions as radiopyo_exception

from .species import RawSpecies, SimSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray


__all__ = [
    "BaseChemicalReaction",
    "IsChemicalReaction",
]

# ---------------------------------------------------------------------------- #
#                          BASIC CHEMISTRY DEFINITIONS                         #
# ---------------------------------------------------------------------------- #


class IsChemicalReaction(Protocol):
    """ """
    reactants: List[SimSpecies]
    products: List[SimSpecies]
    stoi_reactants: List[int]
    stoi_products: List[int]
    reactants_idx: List[int | None]
    products_idx: List[int | None]

    def iter_species(self) -> Iterable[SimSpecies]:
        """Iter over species involved in the reaction"""

    def iter_reactants(self) -> Iterable[Tuple[RawSpecies, int]]:
        """ Iter over reactants only"""

    def iter_products(self) -> Iterable[Tuple[RawSpecies, int]]:
        """Iter over products only"""


class BaseChemicalReaction():
    """"""
    reactants: Dict[str, SimSpecies]
    products: Dict[str, SimSpecies]
    stoi_reactants: Dict[str, int]
    stoi_products: Dict[str, int]

    def __init__(self) -> None:
        self.reactants = {}
        self.products = {}
        self.stoi_products = {}
        self.stoi_reactants = {}

    def __contains__(self, other: object) -> bool:
        """"""
        if isinstance(other, (str, RawSpecies)):
            return len([elt for elt in self.iter_species() if elt == other]) > 0
        raise NotImplementedError

    def tot_number_reactants(self) -> int:
        return reduce(add, list(self.stoi_reactants.values()))

    def iter_species(self) -> Iterable[SimSpecies]:
        """Iter over species involved in the reaction"""
        return chain(self.reactants.values(), self.products.values())

    def iter_reactants(self) -> Iterable[Tuple[SimSpecies, int]]:
        """ Iter over reactants only"""
        return [(sp, self.stoi_reactants[label])
                for label, sp in self.reactants.items()]

    def iter_products(self) -> Iterable[Tuple[SimSpecies, int]]:
        """Iter over products only"""
        return [(sp, self.stoi_products[label])
                for label, sp in self.products.items()]

    def set_index(self, index: int) -> None:
        """"""
        self.index = index

    def number_of_reactants(self) -> float:
        """Count the number of reactants involved"""
        return len(list(self.iter_reactants()))

    def number_of_products(self) -> float:
        """Count the number of reactants involved"""
        return len(list(self.iter_products()))

    def is_product(self, product: SimSpecies | str) -> int:
        """"""
        if isinstance(product, str):
            return self.stoi_products.get(product, 0)
        return self.stoi_products.get(product.label, 0)

    def is_reactant(self, reactant: SimSpecies | str) -> int:
        """"""
        if isinstance(reactant, str):
            return self.stoi_reactants.get(reactant, 0)
        return self.stoi_reactants.get(reactant.label, 0)

    def _get_product_idx(self, product: SimSpecies | str) -> int:
        """"""
        key = product
        if isinstance(product, RawSpecies):
            key = product.label
        for idx, (sp, _) in enumerate(self.iter_products()):
            if sp.label == key:
                return idx
        raise radiopyo_exception.NotReactionProduct

    def _get_reactant_idx(self, reactant: SimSpecies | str) -> int:
        """"""
        key = reactant
        if isinstance(reactant, RawSpecies):
            key = reactant.label
        for idx, (sp, _) in enumerate(self.iter_reactants()):
            if sp.label == key:
                return idx
        raise radiopyo_exception.NotReactionReactant

    def add_as_product(self, sp: SimSpecies) -> None:
        """"""
        if sp.label in self.products:
            self.stoi_products[sp.label] += 1
        else:
            self.products[sp.label] = sp
            self.stoi_products[sp.label] = 1

    def add_as_reactant(self, sp: SimSpecies) -> None:
        """"""
        if sp.label in self.reactants:
            self.stoi_reactants[sp.label] += 1
        else:
            self.reactants[sp.label] = sp
            self.stoi_reactants[sp.label] = 1
