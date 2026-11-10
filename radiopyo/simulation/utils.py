#!/usr/bin/env python3
"""
Backward-compatible facade: re-exports the content of ``units``,
``result_collection`` and ``persistence``.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 31/07/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.simulation.persistence import (
    AlreadyExistingGroup,
    GroupNotFound,
    H5Manager,
    HDFManagerException,
    SimCollection,
    SimIdentifier,
)
from radiopyo.simulation.result_collection import ResCollection
from radiopyo.simulation.units import (
    CO2Convert,
    ConcentrationConvert,
    N2Convert,
    O2Convert,
    SpeciesConvert,
)

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List = [
    "AlreadyExistingGroup",
    "CO2Convert",
    "ConcentrationConvert",
    "GroupNotFound",
    "H5Manager",
    "HDFManagerException",
    "N2Convert",
    "O2Convert",
    "ResCollection",
    "SimCollection",
    "SimIdentifier",
    "SpeciesConvert",
]
