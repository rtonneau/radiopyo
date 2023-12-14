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

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from .k_reactions import KReaction, KReactionCollection

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.parser.reaction_parser import GenericReaction

    from .collections import ListSpecies

# ------------------------------------------------------------------------------------ #

__all__ = [
    "EReaction",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class EReaction(KReactionCollection):
    """
    The math for equilibrium constant are the following:
    K_equilibrium = K_forward / K_backward
    """
    _k_equilibrium: tp.Optional[np.polynomial.Polynomial]
    _k_backward: tp.Optional[np.polynomial.Polynomial]
    _k_forward: np.polynomial.Polynomial
    T_inv: float

    def __init__(self,
                 forward: GenericReaction,
                 backward: GenericReaction,
                 constant: tp.Optional[float | tp.List[float]] = None,
                 temperature: tp.Optional[float] = None,
                 list_sp: tp.Optional[ListSpecies] = None,
                 ):
        """
        """
        super().__init__()
        # Check whether enough information is provided to compute all K's.
        if constant is None:
            self._k_forward = np.polynomial.Polynomial(
                forward["constants"]["k_value"])
            self._k_backward = np.polynomial.Polynomial(
                backward["constants"]["k_value"])
            self._k_equilibrium = None
        elif ("k_value" not in forward["constants"] or
              forward["constants"]["k_value"] is None):
            self._k_backward = np.polynomial.Polynomial(
                backward["constants"]["k_value"])
            self._k_equilibrium = np.polynomial.Polynomial(constant)
            self._k_forward = self._k_equilibrium * self._k_backward
        elif ("k_value" not in backward["constants"] or
              backward["constants"]["k_value"] is None):
            self._k_forward = np.polynomial.Polynomial(
                forward["constants"]["k_value"])
            self._k_equilibrium = np.polynomial.Polynomial(constant)
            self._k_backward = None
        else:
            raise ValueError(
                "EReaction needs 2 out of the 3 (K_equilibrium, K_forward, K_backward)")

        # Need to add forward and backward reactions via the superclass
        # KReactionCollection as this method is not allowed for an EReaction.
        super().add_reaction("forward", forward, list_sp)
        super().add_reaction("backward", backward, list_sp)

    def add_reaction(self,
                     label: str,
                     generic: GenericReaction,
                     list_sp: ListSpecies | None = None,
                     ) -> None:
        raise NotImplementedError(
            "Not possible to add more reaction to an equilibrium collection")

    def k_equilibrium(self, T: tp.Optional[float] = None) -> float:
        """
        """
        T = T or self.T
        if self._k_equilibrium is not None:
            return self._k_equilibrium(T)
        if self._k_backward is not None:
            return self._k_forward(T) / self._k_backward(T)
        raise ValueError(
            "Not enough information about reaction constant to compute a value")

    def k_forward(self, T: tp.Optional[float] = None) -> float:
        """
        """
        T = T or self.T
        return self.k_forward(T)

    def k_backward(self, T: tp.Optional[float] = None) -> float:
        """ 
        """
        T = T or self.T
        if self._k_backward is not None:
            return self._k_backward(T)
        if self._k_equilibrium is not None:
            return self._k_forward(T) / self._k_equilibrium(T)
        raise ValueError(
            "Not enough information about reaction constant to compute a value")

    def iter_reactions_and_kvalues(self,
                                   T: tp.Optional[float] = None,
                                   ) -> tp.Iterator[str, KReaction, float]:
        """
        """
        T = T or self.T
        yield "forward", self.get_reaction("forward"), self.k_forward(T)
        yield "backward", self.get_reaction("backward"), self.k_backward(T)
