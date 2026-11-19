#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 01/03/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp
from pathlib import Path

# ----------------------------- Internal imports ---------------------------- #
from radiopyo import rust
from radiopyo.parser.exceptions import (
    ConfigError,
    ConfigSyntaxError,
    IncludeNotFoundError,
)
from radiopyo.parser.utils import bundled_configs
from radiopyo.physics.beam import (
    ConstantBeam,
    MacroPulsedBeam,
    PulsedBeam,
    SinglePulseBeam,
    min_beam_on_time,
    set_min_beam_on_time,
)
from radiopyo.simulation.unit_cell import UnitCell
from radiopyo.utils import (
    LABARBE_CONFIG_FILE,
    PACKAGE_DATA_PATH,
    PURE_WATER_CONFIG_FILE,
)
from radiopyo.utils.logging import (
    disable_logging,
    enable_logging,
    logConfig_from_pyproject,
    set_logging_level,
)
from radiopyo.utils.progress import progress_bar_enabled, set_progress_bar
from radiopyo.utils.version import VERSION

__version__ = VERSION

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    pass

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List = [
    "LABARBE_CONFIG_FILE",
    "PACKAGE_DATA_PATH",
    "PURE_WATER_CONFIG_FILE",
    "ConfigError",
    "ConfigSyntaxError",
    "ConstantBeam",
    "IncludeNotFoundError",
    "MacroPulsedBeam",
    "PulsedBeam",
    "SinglePulseBeam",
    "UnitCell",
    "__version__",
    "bundled_configs",
    "disable_logging",
    "enable_logging",
    "min_beam_on_time",
    "progress_bar_enabled",
    "rust",
    "set_logging_level",
    "set_min_beam_on_time",
    "set_progress_bar",
]


# ---------------------------------- Logging -------------------------------- #
logConfig_from_pyproject(Path(__file__).parent / "utils/default_log_config.toml")

# --------------------------------------------------------------------------- #
#                            DECORATOR DEFINITIONS                            #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#                              FUNCTION DEFINITIONS                           #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #
