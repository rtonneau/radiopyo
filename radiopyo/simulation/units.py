#!/usr/bin/env python3
"""
Unit converters between molar concentrations and gas pressures/mass units.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 31/07/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

# ----------------------------- Internal imports ---------------------------- #

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List = [
    "CO2Convert",
    "ConcentrationConvert",
    "N2Convert",
    "O2Convert",
    "SpeciesConvert",
]

# ---------------------------------- Logging -------------------------------- #

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #

class ConcentrationConvert(object):
    """ """

    HENRY: float  # = 1.3e-3 # units are M/atm, for O2
    mass_per_mole: float  # e.g. O2: 32 g / mole

    def __init__(self, henry_cst: float, mass: tp.Optional[float] = None) -> None:
        self.HENRY = henry_cst
        self.mass_per_mole = mass if mass is not None else -1.0

    def molar_to_atm(self, value: float | NDArray) -> float | NDArray:
        """ """
        return value / self.HENRY

    def atm_to_molar(self, value: float | NDArray) -> float | NDArray:
        """ """
        return value * self.HENRY

    def molar_to_partial(self, value: float | NDArray) -> float | NDArray:
        """ """
        return self.molar_to_atm(value) * 100

    def partial_to_molar(self, value: float | NDArray) -> float | NDArray:
        """ """
        return self.atm_to_molar(value / 100)

    def molar_to_mmHg(self, value: float | NDArray) -> float | NDArray:
        """ """
        return self.molar_to_atm(value) * 760

    def mmHg_to_molar(self, value: float | NDArray) -> float | NDArray:
        """ """
        return self.atm_to_molar(value / 760)

    def mg_per_l_to_molar(self, value: float | NDArray) -> float | NDArray:
        """ """
        if self.mass_per_mole <= 0.0:
            raise ValueError(
                f"Molar mass of species not defined (value={self.mass_per_mole})"
            )
        return value * 1e-3 / self.mass_per_mole


class SpeciesConvert(object):
    """ """

    convertor: ConcentrationConvert
    value: float | NDArray  # Concentration in Molar unit

    def __init__(self, value: float | NDArray) -> None:
        """value is a concentration in Molar unit"""
        self.value = value

    @classmethod
    def from_atm(cls: tp.Type[tp.Self], value: float | NDArray) -> tp.Self:
        try:
            return cls(cls.convertor.atm_to_molar(value))
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    @classmethod
    def from_molar(cls: tp.Type[tp.Self], value: float | NDArray) -> tp.Self:
        try:
            return cls(value)
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    @classmethod
    def from_umolar(cls: tp.Type[tp.Self], value: float | NDArray) -> tp.Self:
        try:
            return cls(value * 1e-6)
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    @classmethod
    def from_partial(cls: tp.Type[tp.Self], value: float | NDArray) -> tp.Self:
        try:
            return cls(cls.convertor.partial_to_molar(value))
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    @classmethod
    def from_mmHg(cls: tp.Type[tp.Self], value: float) -> tp.Self:
        try:
            return cls(cls.convertor.mmHg_to_molar(value))
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    @classmethod
    def from_mg_per_l(cls: tp.Type[tp.Self], value: float) -> tp.Self:
        try:
            return cls(cls.convertor.mg_per_l_to_molar(value))
        except AttributeError as err:
            raise TypeError(f"No convertor (with Henry cst) defined for {cls}") from err

    def to_atm(self) -> float | NDArray:
        """ """
        return self.convertor.molar_to_atm(self.value)

    def to_mmHg(self) -> float | NDArray:
        """ """
        return self.convertor.molar_to_mmHg(self.value)

    def to_partial(self) -> float | NDArray:
        """ """
        return self.convertor.molar_to_partial(self.value)

    def to_molar(self) -> float | NDArray:
        """ """
        return self.value

    def to_umolar(self) -> float | NDArray:
        """ """
        return self.value * 1e6


class O2Convert(SpeciesConvert):
    convertor = ConcentrationConvert(henry_cst=1.3e-3, mass=32.0)


class CO2Convert(SpeciesConvert):
    convertor = ConcentrationConvert(henry_cst=3.4e-2, mass=44.0)


class N2Convert(SpeciesConvert):
    convertor = ConcentrationConvert(henry_cst=6.1e-4, mass=28.0)
