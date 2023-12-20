#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 18/12/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = [
    "Temperature",
]

# -------------------------------------- Logging ------------------------------------- #

# ------------------------------------------------------------------------------------ #
#                                DECORATOR DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #
class Temperature(object):
    """ """

    _DEFAULT: float = 25  # Default temperature in °C
    _value: float  # Temperature in Kelvin

    def __init__(
        self, value: tp.Optional[float | Temperature] = None, unit: str = "celsius"
    ) -> None:
        """ """
        if isinstance(value, Temperature):
            value = value._value
            unit = "kelvin"
        elif value is None:
            value = self._DEFAULT
        unit = unit.lower()
        match unit:
            case "celsius":
                self._value = self.celsius_to_kelvin(value)
            case "kelvin":
                if value < 0.0:
                    raise ValueError(
                        f"Non-sense negative temperature (in kelvin): {value}"
                    )
                self._value = value
            case _:
                raise ValueError(f"Unknown temperature unit: {unit}")

    @property
    def kelvin(self) -> float:
        return self._value

    @property
    def celsius(self) -> float:
        return self.kelvin_to_celsius(self._value)

    @staticmethod
    def celsius_to_kelvin(value: float) -> float:
        if value < -273.15:
            raise ValueError(f"Non-sense negative temperature: {value}")
        return value + 273.15

    @staticmethod
    def kelvin_to_celsius(value: float) -> float:
        if value < 0.0:
            raise ValueError(f"Non-sense negative temperature: {value}")
        return value - 273.15

    def get(self, unit: str = "celsius") -> float:
        """ """
        match unit:
            case "celsius":
                return self.celsius
            case "kelvin":
                return self.kelvin
            case _:
                raise ValueError(f"Unknown temperature unit: {unit}")

    def copy(self) -> tp.Self:
        """ """
        return self.__class__(value=self._value, unit="kelvin")
