#!/usr/bin/env python3
"""
The acid_base submodule contains all class definitions related to acid-base equilibrium.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 14/12/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np

from radiopyo.parser.reaction_parser import GenericReaction, revert_reaction

# --------------------------------- Internal imports --------------------------------- #
from .k_reactions import (
    KReaction,
    KReactionCollection,
    polynomial_factory,
)
from .utils.units import Temperature

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from .collections import ListSpecies

__all__ = [
    "AcidBase",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class AcidBase(KReactionCollection):
    """
    This one is more tricky about temperature management as various polynomial forms
    are involved:
        pKa = a0 + a1 T + a2 T**2 + ...
        log(K) = a0 + a1 (1/T) + a2 (1/T)**2 + ...
    => pKa uses T while log(K) uses (1/T)
    Therefore, the calculation of reaction constants cannot be left to the respective
    KReactions as these are not aware of it.

    Rem:
    K_equilibrium = K_forward / K_backward
    <=> log(K_e) = log(K_f) - log(K_b)
    K_equilibrium = 10 ** (-pKa)
    """

    _pKa: np.polynomial.Polynomial
    _k_backward: tp.Optional[np.polynomial.Polynomial]
    _k_forward: tp.Optional[np.polynomial.Polynomial]

    def __init__(
        self,
        reaction: GenericReaction,
        pKa: float | tp.List[float],
        list_sp: tp.Optional[ListSpecies] = None,
    ):
        """ """
        super().__init__()
        self._pKa = np.polynomial.Polynomial(pKa)

        if (
            "k_backward" in reaction["constants"]
            and reaction["constants"]["k_backward"] is not None
        ):
            self._k_backward = polynomial_factory(reaction["constants"]["k_backward"])
            self._k_forward = None
        elif (
            "k_forward" in reaction["constants"]
            and reaction["constants"]["k_forward"] is not None
        ):
            self._k_forward = polynomial_factory(reaction["constants"]["k_forward"])
            self._k_backward = None
        else:
            raise ValueError("AcidBase needs at least one of K_forward or K_backward")

        # Need to add forward and backward reactions via the superclass
        # KReactionCollection as this method is not allowed for an EReaction.
        super().add_reaction("forward", reaction, list_sp)
        super().add_reaction("backward", revert_reaction(reaction), list_sp)

    @property
    def forward(self) -> KReaction:
        """ """
        return self._reactions[0][1]

    @property
    def backward(self) -> KReaction:
        """ """
        return self._reactions[1][1]

    def add_reaction(
        self,
        label: str,
        generic: GenericReaction,
        list_sp: ListSpecies | None = None,
    ) -> None:
        raise NotImplementedError(
            "Not possible to add more reaction to an equilibrium collection"
        )

    def k_equilibrium(self, T: tp.Optional[float] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)

        For pKa, polynomial form depends on T in °C.
        """
        _T = Temperature(T, unit="celsius")
        return 10 ** (-self._pKa(_T.celsius))

    def k_forward(self, T: tp.Optional[float] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        _T = Temperature(T, unit="celsius")
        if self._k_backward is not None:
            _val = self._k_backward(1.0 / _T.kelvin) - self._pKa(_T.celsius)
            return 10**_val
        if self._k_forward is None:
            raise ValueError("AcidBase: Both K_forward and K_backward are None")
        return 10 ** self._k_forward(1.0 / _T.kelvin)

    def k_backward(self, T: tp.Optional[float] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        _T = Temperature(T, unit="celsius")
        if self._k_forward is not None:
            _val = self._k_forward(1.0 / _T.kelvin) + self._pKa(_T.celsius)
            return 10**_val
        if self._k_backward is None:
            raise ValueError("AcidBase: Both K_forward and K_backward are None")
        return 10 ** self._k_backward(1.0 / _T.kelvin)

    def iter_reactions_and_kvalues(
        self,
        T: tp.Optional[float] = None,
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        yield "forward", self.get_reaction("forward"), self.k_forward(T)
        yield "backward", self.get_reaction("backward"), self.k_backward(T)
