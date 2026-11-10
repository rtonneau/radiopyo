#!/usr/bin/env python3
"""
Definition of some Physics (beam) level exceptions.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 08/05/2023
"""
# ----------------------------- External imports ---------------------------- #
# ----------------------------- Internal imports ---------------------------- #
from radiopyo.exceptions import RadioPyoError

# ------------------------- Type checking ONLY imports----------------------- #
# --------------------------------------------------------------------------- #

__all__ = [
    "BeamError",
    "ParameterError",
]


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #
class BeamError(RadioPyoError):
    """"""


class ParameterError(BeamError):
    """"""
