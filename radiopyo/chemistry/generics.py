#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 04/04/2025
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import typing as tp
from enum import Enum

# ----------------------------- Internal imports ---------------------------- #
from .species import BasicSpecies

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    pass

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List[str] = [
    "GenericReaction",
    "Operation",
    "dict_to_generic_reaction",
    "revert_reaction",
]

# ---------------------------------- Logging -------------------------------- #

# --------------------------------------------------------------------------- #
#                            DECORATOR DEFINITIONS                            #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
def revert_reaction(reaction: GenericReaction) -> GenericReaction:
    """ """
    # Revert the reaction only if all products are related to addition
    # operation.
    for _, (_, _, ope) in reaction["products"].items():
        if ope != Operation.ADD:
            raise TypeError(
                f"Cannot revert reaction with subtraction operation: {reaction}"
            )
    return GenericReaction(
        products=reaction["reactants"].copy(),
        reactants=reaction["products"].copy(),
        constants=reaction["constants"].copy(),
        enzyme=copy.copy(reaction["enzyme"]),
    )


def dict_to_generic_reaction(other: dict) -> GenericReaction:
    return GenericReaction(
        products={k: tuple(v) for k, v in other.get("products", {}).items()},
        reactants={k: tuple(v) for k, v in other.get("reactants", {}).items()},
        constants=dict(other.get("constants", {}).items()),
        enzyme=other.get("enzyme"),
    )


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #
class Operation(Enum):
    ADD = 1
    SUBTRACT = 2

    def __str__(self) -> str:
        match self:
            case Operation.ADD:
                return "+"
            case Operation.SUBTRACT:
                return "-"


class GenericReaction(tp.TypedDict):
    products: tp.Dict[str, tp.Tuple[BasicSpecies, int, Operation]]
    reactants: tp.Dict[str, tp.Tuple[BasicSpecies, int, Operation]]
    constants: tp.Dict[str, float | tp.List[float]]
    enzyme: BasicSpecies | None
