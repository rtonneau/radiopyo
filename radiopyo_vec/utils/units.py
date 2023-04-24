#!/usr/bin/env python3
"""
The unit submodule contains class definitions to handle scientific unit conversion.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from enum import Enum

__all__ = [
    "Units",
    "SimUnit",
]


class Units(float, Enum):
    NANO = 1e-9
    MICRO = 1e-6
    MILLI = 1e-3
    UNITY = 1
    KILO = 1e3
    MEGA = 1e6
    GIGA = 1e9


class SimUnit:
    unit: Units

    def __init__(self, unit: str | Units):
        self.set_base_unit(unit)

    @staticmethod
    def str_to_unit(unit: str) -> Units:
        for elt in Units:
            if unit.lower() == elt.name.lower():
                return elt
        raise KeyError(f"{unit} is not a valid Unit "
                       f"({[elt.name for elt in Units]}).")

    def convert_to(self, value: float, unit: Units) -> float:
        return value * self.unit.value / unit.value

    def convert_from(self, value: float, unit: Units) -> float:
        return value * unit.value / self.unit.value

    def set_base_unit(self, unit: str | Units) -> None:
        if isinstance(unit, str):
            unit = SimUnit.str_to_unit(unit)
        self.unit = unit

    def __str__(self) -> str:
        return f"SimUnit({self.unit})"

    def __repr__(self) -> str:
        return self.__str__()
