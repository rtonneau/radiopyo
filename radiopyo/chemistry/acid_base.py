#!/usr/bin/env python3
"""
The acid_base submodule contains all class definitions related to acid-base equilibrium.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 14/12/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp
from contextlib import suppress

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.parser import RPARSER
from radiopyo.parser.reaction_parser import GenericReaction, revert_reaction

from .collections import ListSpecies
from .k_reactions import (
    BiReaction,
    KReactionCollection,
    polynomial_factory,
)
from .utils.units import Temperature

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

__all__ = [
    "AcidBase",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class AcidBase(BiReaction):
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
        KReactionCollection.add_reaction(
            self, reaction, label="forward", list_sp=list_sp
        )
        KReactionCollection.add_reaction(
            self, revert_reaction(reaction), label="backward", list_sp=list_sp
        )

    @classmethod
    def from_kwargs(
        cls,
        list_sp: tp.Optional[ListSpecies] = None,
        **kwargs: str | tp.List[str] | float | tp.List[float],
    ) -> tp.Self:
        """
        This method combines the two possible definitions of a KReaction in a config
        file.
        TODO: remove the dict approach to keep only the single string reaction.
        """
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        try:
            gen = RPARSER.parse_reaction(kwargs["reaction"])  # type:ignore[arg-type]
            pKa = kwargs["pKa"]
        except KeyError as err:
            raise AttributeError(
                "Need a reaction string AND a pKa to create an Acid Base rel."
            ) from err

        if "k_value" in kwargs:
            raise AttributeError(
                "'k_value' is ambiguous for an AcidBase reaction. "
                "Use k_forward or k_backward instead."
            )

        with suppress(KeyError):
            gen["constants"]["k_forward"] = kwargs["k_forward"]  # type:ignore[assignment]
        with suppress(KeyError):
            gen["constants"]["k_backward"] = kwargs["k_backward"]  # type:ignore[assignment]

        return cls(gen, pKa=pKa, list_sp=list_sp)  # type:ignore[arg-type]

    def k_equilibrium(self, T: tp.Optional[Temperature] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)

        For pKa, polynomial form depends on T in °C.
        """
        _T = T if T is not None else Temperature()
        return 10 ** (-self._pKa(_T.celsius))

    def k_forward(self, T: tp.Optional[Temperature] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        _T = T if T is not None else Temperature()
        if self._k_backward is not None:
            _val = self._k_backward(1.0 / _T.kelvin) - self._pKa(_T.celsius)
            return 10**_val
        if self._k_forward is None:
            raise ValueError("AcidBase: Both K_forward and K_backward are None")
        return 10 ** self._k_forward(1.0 / _T.kelvin)

    def k_backward(self, T: tp.Optional[Temperature] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        _T = T if T is not None else Temperature()
        if self._k_forward is not None:
            _val = self._k_forward(1.0 / _T.kelvin) + self._pKa(_T.celsius)
            return 10**_val
        if self._k_backward is None:
            raise ValueError("AcidBase: Both K_forward and K_backward are None")
        return 10 ** self._k_backward(1.0 / _T.kelvin)
