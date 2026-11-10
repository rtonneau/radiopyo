#!/usr/bin/env python3
"""
Definition of some Module level exceptions.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 17/03/2023
"""
# ----------------------------- External imports ---------------------------- #
# ----------------------------- Internal imports ---------------------------- #
from radiopyo.exceptions import RadioPyoError

# ------------------------ Type checking ONLY imports ----------------------- #
# --------------------------------------------------------------------------- #

__all__ = [
    "NoBeamDefinedError",
    "RadioPyoError",
    "SimulationError",
]

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class SimulationError(RadioPyoError):
    """"""


class NoBeamDefinedError(SimulationError):
    """"""
