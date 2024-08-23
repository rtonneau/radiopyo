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
from contextlib import suppress

import numpy as np

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser import RPARSER
from radiopyo.parser.reaction_parser import revert_reaction

from .species_collection import ListSpecies
from .k_reactions import (
    BiReaction,
    polynomial_factory,
)
from .utils.units import Temperature

# ------------------------ Type checking ONLY imports------------------------ #
if tp.TYPE_CHECKING:
    from radiopyo.parser.reaction_parser import GenericReaction


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
                self._k_forward = polynomial_factory(reaction["constants"]["k_forward"])
                self._k_backward = polynomial_factory(
                    reaction["constants"]["k_backward"]
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
                f"{reaction} & constant: {constant}"
            )

        self.add_forward_reaction(reaction, list_sp)
        self.forward.set_k_value(self._k_forward)  # type: ignore [arg-type]

        self.add_backward_reaction(revert_reaction(reaction), list_sp)
        self.backward.set_k_value(self._k_backward)  # type: ignore [arg-type]

    @classmethod
    def from_kwargs(
        cls,
        list_sp: tp.Optional[ListSpecies] = None,
        **kwargs: str | tp.List[str] | float | tp.List[float],
    ) -> tp.Self:
        """
        Create an EquilibriumReaction from keyword arguments.

        Args:
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.
            **kwargs: Keyword arguments to create the instance.
                - reaction (str): The reaction string.
                - k_forward (float): The forward reaction rate constant.
                - k_backward (float): The backward reaction rate
                    constant.
                - k_equilibrium (float): The equilibrium constant.

        Returns:
            EquilibriumReaction: An instance of EquilibriumReaction.

        Raises:
            AttributeError: If 'k_value' is provided, as it is
                ambiguous for an EquilibriumReaction.
            AttributeError: If no reaction string is provided.
        """
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        try:
            gen = RPARSER.parse_reaction(kwargs["reaction"])  # type:ignore[arg-type]
        except KeyError as err:
            raise AttributeError(
                "Need a reaction string to create an Equilibrium Reaction"
            ) from err

        if "k_value" in kwargs:
            raise AttributeError(
                "'k_value' is ambiguous for an Equilibrium Reaction. "
                "Use k_forward or k_backward instead."
            )

        with suppress(KeyError):
            gen["constants"]["k_forward"] = kwargs[
                "k_forward"
            ]  # type:ignore[assignment]
        with suppress(KeyError):
            gen["constants"]["k_backward"] = kwargs[
                "k_backward"
            ]  # type:ignore[assignment]

        constant = kwargs.get("k_equilibrium", None)

        # type:ignore[arg-type]
        return cls(gen, constant=constant, list_sp=list_sp)  # type:ignore[arg-type]

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
        _T = Temperature(T) if T is not None else Temperature(unit="celsius")
        return 10 ** self._k_equilibrium(1 / _T.kelvin)

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
