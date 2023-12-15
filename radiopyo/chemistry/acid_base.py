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

# --------------------------------- Internal imports --------------------------------- #
from .k_reactions import KReaction, KReactionCollection

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.parser.reaction_parser import GenericReaction

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
        forward: GenericReaction,
        backward: GenericReaction,
        pKa: float | tp.List[float],
        list_sp: tp.Optional[ListSpecies] = None,
    ):
        """ """
        super().__init__()
        self._pKa = np.polynomial.Polynomial(pKa)

        if (
            "k_value" not in forward["constants"]
            or forward["constants"]["k_value"] is None
        ):
            self._k_backward = np.polynomial.Polynomial(
                backward["constants"]["k_value"]
            )
            self._k_forward = None
        elif (
            "k_value" not in backward["constants"]
            or backward["constants"]["k_value"] is None
        ):
            self._k_forward = np.polynomial.Polynomial(forward["constants"]["k_value"])
            self._k_backward = None
        else:
            raise ValueError("AcidBase needs at least one of K_forward or K_backward")

        # Need to add forward and backward reactions via the superclass
        # KReactionCollection as this method is not allowed for an EReaction.
        super().add_reaction("forward", forward, list_sp)
        super().add_reaction("backward", backward, list_sp)

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
        """ """
        T = T if T is not None else KReaction._DEFAULT_T
        return 10 ** (-self._pKa(T))

    def k_forward(self, T: tp.Optional[float] = None) -> float:
        """ """
        T = 1.0 / T if T is not None else 1.0 / KReaction._DEFAULT_T
        if self._k_backward is not None:
            _val = self._k_backward(T) - self._pKa(1.0 / T)
            return 10**_val
        if self._k_forward is None:
            raise ValueError("AcidBase: Both K_forward and K_backward are None")
        return 10 ** self._k_forward(T)

    def k_backward(self, T: tp.Optional[float] = None) -> float:
        """ """
        T = 1.0 / T if T is not None else 1.0 / KReaction._DEFAULT_T
        if self._k_forward is not None:
            _val = self._k_forward(T) + self._pKa(1.0 / T)
            return 10**_val
        if self._k_backward is None:
            raise ValueError("AcidBase: Both K_forward and K_backward are None")
        return 10 ** self._k_backward(T)

    def iter_reactions_and_kvalues(
        self,
        T: tp.Optional[float] = None,
    ) -> tp.Iterator[tp.Tuple[str, KReaction, float]]:
        """ """
        T = T if T is not None else KReaction._DEFAULT_T
        yield "forward", self.get_reaction("forward"), self.k_forward(T)
        yield "backward", self.get_reaction("backward"), self.k_backward(T)
