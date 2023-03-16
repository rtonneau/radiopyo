from __future__ import annotations

from enum import Enum
from itertools import chain
from typing import (
    TYPE_CHECKING,
    Iterable,
    List,
    NamedTuple,
    Protocol,
    Tuple,
    Type,
)

import numpy as np

from radiopyo_vec.chemistry.species import RawSpecies, SimSpecies

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo_vec.utils.collections import MappedSpecies

__all__ = [
    "ReactionRate",
    "ReactionRateEnum",
    "IsChemicalReaction",
    "ReactionSpecies",
    "BaseChemicalReaction",
    "ReactionSpeciesEnum",
    "ReactionType"
]

# ---------------------------------------------------------------------------- #
#                          SPECIFIC SPECIES DEFINITION                         #
# ---------------------------------------------------------------------------- #


class ReactionRateEnum(Enum):
    """"""
    Production = 1
    Consumption = 2


class ReactionType(Enum):
    """"""
    KReaction = 1
    Radiolytic = 2
    Michaelis = 3


class ReactionRate(NamedTuple):
    role: ReactionRateEnum
    reaction_index: int

    @classmethod
    def new_production(cls,
                       reaction_index: int
                       ) -> ReactionRate:
        return cls(ReactionRateEnum.Production, reaction_index)

    @classmethod
    def new_consumption(cls,
                        reaction_index: int,
                        ) -> ReactionRate:
        return cls(ReactionRateEnum.Consumption, reaction_index)


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
    _rate_vec: NDArray[np.float64]

    def init_rate_vec(self, species: MappedSpecies) -> None:
        """"""

    def compute_reaction(self,
                         dose_rate: float,
                         species: MappedSpecies
                         ) -> NDArray[np.float64]:
        """Compute the reaction rate of a reaction"""

    def compute_jacobian(self,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        """Compute the jacobian of the reaction rate"""

    def iter_species(self) -> Iterable[SimSpecies]:
        """Iter over species involved in the reaction"""

    def iter_reactants(self) -> Iterable[Tuple[RawSpecies, int]]:
        """ Iter over reactants only"""

    def iter_products(self) -> Iterable[Tuple[RawSpecies, int]]:
        """Iter over products only"""


class BaseChemicalReaction():
    """"""
    reactants: List[SimSpecies]
    products: List[SimSpecies]
    stoi_reactants: List[int]
    stoi_products: List[int]
    reactants_idx: List[int | None]
    products_idx: List[int | None]
    _rate_vec: NDArray[np.float64]

    def __init__(self) -> None:
        self.reactants = []
        self.products = []
        self.stoi_products = []
        self.stoi_reactants = []
        self.reactants_idx = []
        self.products_idx = []

    def __contains__(self, other: object) -> bool:
        """"""
        if isinstance(other, (str, RawSpecies)):
            return len([elt for elt in self.iter_species() if elt == other]) > 0
        raise NotImplementedError

    def init_rate_vec(self, species: MappedSpecies) -> None:
        self._rate_vec = np.zeros(species.number_tracked_species)
        for idx, stoi in zip(self.reactants_idx, self.stoi_reactants):
            if idx is None:  # No need to add it, related to an untracked species
                continue
            self._rate_vec[idx] -= stoi  # Reactants are consumed!
        for idx, stoi in zip(self.products_idx, self.stoi_products):
            if idx is None:  # No need to add it, related to an untracked species
                continue
            self._rate_vec[idx] += stoi  # Products are ... produced!

    def compute_reaction(self,
                         dose_rate: float,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        """Compute the reaction rate of a reaction"""
        raise NotImplementedError

    def compute_jacobian(self,
                         species: MappedSpecies,
                         ) -> NDArray[np.float64]:
        """Compute the jacobian of the reaction rate"""
        raise NotImplementedError

    def iter_species(self) -> Iterable[SimSpecies]:
        """Iter over species involved in the reaction"""
        return chain(self.reactants, self.products)

    def iter_reactants(self) -> Iterable[Tuple[SimSpecies, int]]:
        """ Iter over reactants only"""
        return zip(self.reactants, self.stoi_reactants)

    def iter_products(self) -> Iterable[Tuple[SimSpecies, int]]:
        """Iter over products only"""
        return zip(self.products, self.stoi_products)

    def set_index(self, index: int) -> None:
        """"""
        self.index = index

    def number_of_reactants(self) -> float:
        """Count the number of reactants involved"""
        return len(list(self.iter_reactants()))

    def number_of_products(self) -> float:
        """Count the number of reactants involved"""
        return len(list(self.iter_products()))

    def get_product_stoichio(self, product: str) -> int:
        """"""
        return self.stoi_products[self._get_product_idx(product)]

    def get_reactant_stoichio(self, reactant: str) -> int:
        """"""
        return self.stoi_reactants[self._get_reactant_idx(reactant)]

    def is_reactant(self, sp: str | Type[RawSpecies]) -> bool:
        """"""
        return len([elt for elt in self.iter_reactants() if elt == sp]) > 0

    def is_product(self, sp: str | Type[RawSpecies]) -> bool:
        """"""
        return len([elt for elt in self.iter_products() if elt == sp]) > 0

    def _get_product_idx(self, product: SimSpecies | str) -> int:
        """"""
        key = product
        if isinstance(product, SimSpecies):
            key = product.label
        for idx, (sp, _) in enumerate(self.iter_products()):
            if sp.label == key:
                return idx
        raise KeyError

    def _get_reactant_idx(self, reactant: SimSpecies | str) -> int:
        """"""
        key = reactant
        if isinstance(reactant, SimSpecies):
            key = reactant.label
        for idx, (sp, _) in enumerate(self.iter_reactants()):
            if sp.label == key:
                return idx
        raise KeyError

    def add_as_product(self, sp: SimSpecies) -> None:
        """"""
        try:
            index = self._get_product_idx(sp)
            self.stoi_products[index] += 1
        except KeyError:
            self.products.append(sp)
            self.stoi_products.append(1)

    def add_as_reactant(self, sp: SimSpecies) -> None:
        """"""
        try:
            index = self._get_reactant_idx(sp)
            self.stoi_reactants[index] += 1
        except KeyError:
            self.reactants.append(sp)
            self.stoi_reactants.append(1)
