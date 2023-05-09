#!/usr/bin/env python3
"""
The utils.logging submodule aims to load logging config from pyproject.toml file.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 08/05/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import logging
import logging.config
from typing import TYPE_CHECKING

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from pathlib import Path
# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "logConfig_from_pyproject",
    "InfoOnly",
    "InfoHandler",
]
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #

class InfoHandler(logging.StreamHandler):
    def __init__(self, *args, **kwargs) -> None: # type: ignore[no-untyped-def]
        logging.StreamHandler.__init__(self, *args, **kwargs)
        self.addFilter(InfoWarning())

class InfoWarning(logging.Filter):
    def filter(self, record:logging.LogRecord) -> bool:
        return record.levelno in [logging.INFO, logging.WARNING]
    
def logConfig_from_pyproject(file: str| Path) -> None:
    """Reads (same way as logging.config.ConfigFile) logging config from pyproject.toml
    file.

    Parameters
    ----------
    file : str| Path
        Path to pyproject.toml

    Raises
    ------
    KeyError
        If tool table not found in TOML file
    KeyError
        If logging not found in tool table
    """
    import toml  # type: ignore[import]

    toml_dict = toml.load(file)

    tool_table = toml_dict.get("tool", {})

    if not tool_table:
        raise KeyError(
            "Tool table not found in TOML file. "
            "See https://peps.python.org/pep-0518/#tool-table"
        )

    tool_table_logging = tool_table.get("logging", {})

    if not tool_table_logging:
        raise KeyError("Logging section not found in tool table. See documentation")

    logging.config.dictConfig(tool_table_logging)