#!/usr/bin/env python3
"""
The species submodule contains all class definitions related to species definitions.
It also provides a collection to ease their handling.

Classes:
    ChemSpecies: Class used for species with cc (cst or dyn).
    RawSpecies: Basic class for species.
    H2O: A predefined species for water.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from .utils.units import Temperature

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ------------------------------------------------------------------------------------ #
__all__ = [
    "ChemSpecies",
    "RawSpecies",
    "H2O",
]


# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class RawSpecies:
    """
    Basic class for species.

    This class provides a basic structure for species, including a label and methods for equality and string representation.

    Attributes:
        label (str): The label of the species.

    Methods:
        __contains__(other): Checks if the given object is equal to the current species.
        __eq__(other): Checks if the given object is equal to the current species.
        __str__(): Returns a string representation of the object.
    """

    label: str

    def __init__(self, label: str):
        """
        Initializes a new instance of the RawSpecies class.

        Args:
            label (str): The label of the species.
        """
        self.label = label

    def __contains__(self, other: object) -> bool:
        """
        Checks if the given object is equal to the current species.

        Args:
            other (object): The object to be compared.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        return self.__eq__(other)

    def __eq__(self, other: object) -> bool:
        """
        Checks if the given object is equal to the current species.

        Args:
            other (object): The object to be compared.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        if isinstance(other, str):
            return self.label == other
        if isinstance(other, RawSpecies):
            return self.label == other.label
        raise NotImplementedError(f"Equality between RawSpecies and {other}")

    def __str__(self) -> str:
        """
        Returns a string representation of the object.

        Returns:
            str: The label of the object.
        """
        return self.label


class ChemSpecies(RawSpecies):
    """
    Class used to define species from concentration coefficient (cc).

    This class extends the RawSpecies class and adds attributes and methods related to the concentration coefficient.

    Attributes:
        label (str): The label of the species.
        cc_init (Polynomial): The initial concentration coefficient of the species.

    Methods:
        set_initial_cc(value): Sets the initial concentration coefficient of the species.
        ... (other methods specific to ChemSpecies)
    """

    cc_init: np.polynomial.Polynomial  # Store initial value as [mol/l]
    index: int

    def __init__(
        self,
        label: str,
        index: int = -1,
        cc_init: tp.Optional[float | tp.List[float]] = None,
    ) -> None:
        """
        Initializes a new instance of the ChemSpecies class.

        Args:
            label (str): The label of the species.
            index (int): The index of the species. Defaults to -1.
            cc_init (float | list[float], optional): The initial concentration of the species. Defaults to None.

        Returns:
            None
        """
        super().__init__(label)
        self.index = index
        cc_init = cc_init if cc_init is not None else -1
        self.cc_init = np.polynomial.Polynomial(cc_init)

    def __repr__(self) -> str:
        """
        Returns a string representation of the object.

        Returns:
            str: The representation of the object.
        """
        return f"ChemSpecies('{self.label}', {self.index}, {self.cc_init})"

    def __copy__(self) -> ChemSpecies:
        """
        Creates a deep copy of the current ChemSpecies instance.

        Returns:
            ChemSpecies: A new ChemSpecies instance with the same label, index, and initial concentration coefficients.
        """
        return ChemSpecies(self.label, self.index, self.cc_init.coef)

    def initial_cc(self, T: tp.Optional[float | Temperature] = None) -> float:
        """
        Returns the initial concentration of the species at a given temperature.

        Parameters:
            T (float | Temperature): Temperature in Celsius. If None, default temperature is used.

        Returns:
            float: The initial concentration of the species at the given temperature.
        """
        T = Temperature() if T is None else Temperature(T)
        return self.cc_init(T.celsius) if self.cc_init(T.celsius) > 0 else 0

    def set_initial_cc(self, value: float | tp.List[float]) -> None:
        """
        Sets the initial concentration coefficient of the species.

        Parameters:
            value (float | list[float]): The initial concentration coefficient(s) of the species.

        Returns:
            None
        """
        self.cc_init = np.polynomial.Polynomial(value)


H2O = ChemSpecies("H2O", cc_init=[55.50, 6.075e-3, 4.110e-4, 1.496e-6, 2.619e-9])


# Define H2O here
def get_H2O_cc(T: tp.Optional[float | Temperature] = None) -> float:
    """
    Returns the initial concentration of H2O at a given temperature.

    Parameters:
        T (float | Temperature): Temperature in Celsius. If None, default temperature is used.

    Returns:
        float: The initial concentration of H2O at the given temperature.
    """
    return H2O.initial_cc(T)
