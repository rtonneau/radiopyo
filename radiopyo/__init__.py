#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 01/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp
from pathlib import Path

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.physics.beam import ConstantBeam, PulsedBeam, SinglePulseBeam
from radiopyo.simulation.unit_cell import UnitCell
from radiopyo.utils import LABARBE_CONFIG_FILE, PACKAGE_DATA_PATH
from radiopyo.utils.logging import (
    disable_logging,
    enable_logging,
    logConfig_from_pyproject,
    set_logging_level,
)

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = [
    "ConstantBeam",
    "PulsedBeam",
    "SinglePulseBeam",
    "UnitCell",
    "disable_logging",
    "enable_logging",
    "set_logging_level",
    "PACKAGE_DATA_PATH",
    "LABARBE_CONFIG_FILE",
    "enable_logging",
    "set_logging_level",
    "enable_logging",
    "set_logging_level",
]


# -------------------------------------- Logging ------------------------------------- #
logConfig_from_pyproject(Path(__file__).parent.parent / "pyproject.toml")

# ------------------------------------------------------------------------------------ #
#                                DECORATOR DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #
