#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 26/07/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp
from pathlib import Path

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = [
    "PACKAGE_DATA_PATH",
    "LABARBE_CONFIG_FILE",
]

# -------------------------------------- Logging ------------------------------------- #

# ------------------------------------------------------------------------------------ #
PACKAGE_DATA_PATH = Path(f"{Path(__file__).parent.parent}//data").absolute()
LABARBE_CONFIG_FILE = PACKAGE_DATA_PATH / "config_Labarbe.toml"
