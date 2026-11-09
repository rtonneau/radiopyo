#!/usr/bin/env python3
"""
The `e_reactions` submodule contains classes that define chemical
reactions with equilibrium. These reactions are defined using the
equilibrium approach, also known as the "chemical equilibrium" or
"thermodynamic equilibrium".

This submodule provides the necessary classes to work with chemical
reactions in which the concentrations of the reactants and products are
in equilibrium.

The classes defined here are:

- `EReaction`: This class represents a chemical reaction with
 equilibrium.It inherits from the `BasicReaction` class and provides
 additional methods to work with these types of reactions.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 13/12/2023
"""

# ---------------------------- External imports ----------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.generics import revert_reaction
from radiopyo.chemistry.k_reactions import (
    BiReaction,
)
from radiopyo.chemistry.utils.factories import polynomial_factory
from radiopyo.chemistry.utils.units import Temperature

# ------------------------ Type checking ONLY imports------------------------ #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.generics import GenericReaction
    from radiopyo.chemistry.species_collection import ListSpecies


# --------------------------------------------------------------------------- #

__all__ = [
    "EReaction",
]

# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #


class EReaction(BiReaction):
    """
    This class represents equilibrium reactions. It computes the
    equilibrium constant by the following formula:
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
        """
        Initializes an EReaction instance.

        This method is responsible for computing the equilibrium
        constant (K_equilibrium) from the given reaction and constant.

        Args:
            reaction (GenericReaction): The reaction to be used for
                computing the equilibrium constant.
            constant (Optional[float | List[float]], optional): The
                constant value for the reaction. Defaults to None.
            list_sp (Optional[ListSpecies], optional): The list of
                species involved in the reaction. Defaults to None.

        Returns:
            None
        """
        super().__init__()
        # Check whether enough information is provided to compute all K's.
        if constant is None:
            try:
                self._k_forward = polynomial_factory(
                    reaction["constants"]["k_forward"],
                    log10=True,
                )
                self._k_backward = polynomial_factory(
                    reaction["constants"]["k_backward"],
                    log10=True,
                )
                if not np.isinf(self._k_forward.coef[0]) and not np.isinf(
                    self._k_backward.coef[0]
                ):
                    self._k_equilibrium = self._k_forward - self._k_backward
            except KeyError as err:
                raise ValueError(
                    "EReaction needs 2 out of the 3 (K_equilibrium, K_forward, K_backward)"
                    f"{reaction} & constant: {constant}"
                ) from err

        elif (
            "k_backward" in reaction["constants"]
            and reaction["constants"]["k_backward"] is not None
        ):
            self._k_backward = polynomial_factory(
                reaction["constants"]["k_backward"],
                log10=True,
            )
            self._k_equilibrium = polynomial_factory(
                constant,
                log10=True,
            )
            self._k_forward = self._k_equilibrium + self._k_backward
        elif (
            "k_forward" in reaction["constants"]
            and reaction["constants"]["k_forward"] is not None
        ):
            self._k_forward = polynomial_factory(
                reaction["constants"]["k_forward"],
                log10=True,
            )
            self._k_equilibrium = polynomial_factory(
                constant,
                log10=True,
            )
            self._k_backward = self._k_forward - self._k_equilibrium
        else:
            raise ValueError(
                "EReaction needs 2 out of the 3 (K_equilibrium, K_forward, K_backward)"
                f"{reaction} & constant: {constant}"
            )

        self.add_forward_reaction(reaction, list_sp)
        self.forward.set_k_value(self._k_forward)

        self.add_backward_reaction(revert_reaction(reaction), list_sp)
        self.backward.set_k_value(self._k_backward)

    def k_equilibrium(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Compute the equilibrium constant of the reaction at a given
        temperature.

        Args:
            T (Optional[Temperature | float]): The temperature in °C
                (optional, default:None)

        Returns:
            float: The equilibrium constant at the given temperature.
        """
        _T = Temperature(T)
        return tp.cast(float, 10 ** self._k_equilibrium(1 / _T.kelvin))

    def k_forward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Compute the forward reaction constant of the reaction at a given
        temperature.

        Args:
            T (Optional[Temperature | float]): The temperature in
                °C (optional, default:None)

        Returns:
            float: The forward reaction constant at the given
                temperature.
        """
        return self.forward.k_value(T)

    def k_backward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Compute the backward reaction constant of the reaction at a
            given temperature.

        Args:
            T (Optional[Temperature | float]): The temperature in
                °C (optional, default:None)

        Returns:
            float: The backward reaction constant at the
                given temperature.
        """
        return self.backward.k_value(T)

    def to_dict(self) -> tp.Dict[str, tp.Any]:
        """
        Convert the reaction to a dictionary.

        Returns:
            dict    A dictionary representation of the reaction.
        """
        return {
            "reaction": self.__str__(),
            "k_forward": self.forward.raw_k_value(),
            "k_backward": self.backward.raw_k_value(),
        }
