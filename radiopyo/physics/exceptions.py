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
    "ShortBeamOnError",
]


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #
class BeamError(RadioPyoError):
    """"""


class ParameterError(BeamError):
    """"""


class ShortBeamOnError(BeamError, ValueError):
    """
    A beam whose configured ON time (pulse on_time, micro or train
    on_time, or the single ON span of a constant / single pulse beam
    with a max dose) is shorter than MIN_BEAM_ON_TIME: it cannot be
    integrated correctly (the initial burst covers most of it).
    """
