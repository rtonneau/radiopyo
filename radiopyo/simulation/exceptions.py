#!/usr/bin/env python3
"""
Definition of some Module level exceptions.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 17/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
# --------------------------------- Internal imports --------------------------------- #
# ----------------------------- Type checking ONLY imports---------------------------- #
# ------------------------------------------------------------------------------------ #

__all__ = [
    "RadioPyoError",
    "SimulationError",
    "NoBeamDefineError",
]

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class RadioPyoError(Exception):
    """Generic Module Related Exception"""


class SimulationError(RadioPyoError):
    """"""


class NoBeamDefinedError(SimulationError):
    """"""
