#!/usr/bin/env python3
"""
The e_reactions submodule contains all class definitions related to chemical
reactions with equilibrium.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 13/12/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.parser.reaction_parser import revert_reaction

from .k_reactions import (
    BiReaction,
    KReaction,
    KReactionCollection,
    polynomial_factory,
)
from .utils.units import Temperature

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    import numpy as np

    from radiopyo.parser.reaction_parser import GenericReaction

    from .collections import ListSpecies

# ------------------------------------------------------------------------------------ #

__all__ = [
    "EReaction",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class EReaction(BiReaction):
    """

    The math for equilibrium constant are the following:
    K_equilibrium = K_forward / K_backward
    <=> log(K_e) = log(K_f) - log(K_b)

    For polynomial form, log(K) = a0 + a1 (1/T) + a2 (1/T)**2 + ...
    """

    _k_equilibrium: np.polynomial.Polynomial

    def __init__(
        self,
        reaction: GenericReaction,
        constant: tp.Optional[float | tp.List[float]] = None,
        list_sp: tp.Optional[ListSpecies] = None,
    ):
        """ """
        super().__init__()
        # Check whether enough information is provided to compute all K's.
        if constant is None:
            try:
                self._k_forward = polynomial_factory(reaction["constants"]["k_forward"])
                self._k_backward = polynomial_factory(
                    reaction["constants"]["k_backward"]
                )
                self._k_equilibrium = self._k_forward - self._k_backward
            except KeyError as err:
                raise ValueError(
                    "EReaction needs 2 out of the 3 (K_equilibrium, K_forward, K_backward)"
                ) from err

        elif (
            "k_backward" in reaction["constants"]
            and reaction["constants"]["k_backward"] is not None
        ):
            self._k_backward = polynomial_factory(reaction["constants"]["k_backward"])
            self._k_equilibrium = polynomial_factory(constant)
            self._k_forward = self._k_equilibrium + self._k_backward
        elif (
            "k_forward" in reaction["constants"]
            and reaction["constants"]["k_forward"] is not None
        ):
            self._k_forward = polynomial_factory(reaction["constants"]["k_forward"])
            self._k_equilibrium = polynomial_factory(constant)
            self._k_backward = self._k_forward - self._k_equilibrium
        else:
            raise ValueError(
                "EReaction needs 2 out of the 3 (K_equilibrium, K_forward, K_backward)"
            )

        # Need to add forward and backward reactions via the superclass
        # KReactionCollection as this method is not allowed for an EReaction.
        KReactionCollection.add_reaction(self, "forward", reaction, list_sp)
        self.forward.set_k_value(self._k_forward)  # type: ignore [arg-type]
        KReactionCollection.add_reaction(
            self, "backward", revert_reaction(reaction), list_sp
        )
        self.backward.set_k_value(self._k_backward)  # type: ignore [arg-type]

    def k_equilibrium(self, T: tp.Optional[float] = None) -> float:
        """
        Args:
            T: temperature in °C.
        """
        _T = Temperature(T, unit="celsius")
        return 10 ** self._k_equilibrium(1 / _T.kelvin)

    def k_forward(self, T: tp.Optional[float] = None) -> float:
        """
        Args:
            T: temperature in °C.
        """
        return self.forward.k_value(T)

    def k_backward(self, T: tp.Optional[float] = None) -> float:
        """
        Args:
            T: temperature in °C.
        """
        return self.backward.k_value(T)
