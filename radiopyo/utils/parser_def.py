#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 23/06/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import logging
import os
import typing as tp

from .sim_types import default_ConfigDict

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    import pathlib

    from .sim_types import ConfigDict

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "FileParser",
    "list_available_parsers",
    "resolve_includes",
]

# -------------------------------------- Logging ------------------------------------- #
logger = logging.getLogger("radiopyo")

# ------------------------------------------------------------------------------------ #
#                                DECORATOR DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #


def list_available_parsers() -> tp.Dict[str, tp.Type[FileParser]]:
    return {cls.ext: cls for cls in FileParser.__subclasses__()}


def resolve_includes(includes: tp.Dict[str, tp.List[str]]) -> ConfigDict:
    """ """
    out = default_ConfigDict()
    if len(includes) == 0:
        return out
    available_parsers = list_available_parsers()
    for ext, f_list in includes.items():
        if ("."+ext) not in available_parsers:
            logger.warning(f"Unable to include .{ext} file. "
                           "No parser available")
            continue
        for f in f_list:
            logger.info(f"Resolving include: {f}")
            _ = available_parsers["."+ext](f).parse()
            out = _.update(out)
    return out
# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class FileParser(object):
    """
    Generic Class parsing files.
    """
    path: str | pathlib.Path
    ext: tp.ClassVar[str]

    def __init__(self, path: str | pathlib.Path):
        self.file_type()
        self.path = path
        if not os.path.exists(self.path):
            raise FileNotFoundError(self.path)

    def parse(self) -> ConfigDict:
        """ 
        Method parsing the file (self.path) and returning a ConfigDict.

        Returns
            ConfigDict: TypedDictionary use to build the simulation Environment.
        """
        raise NotImplementedError

    def file_type(self) -> str:
        """
        raises
            AttributeError: if ext class attribute is not defined by the user.

        returns
            str: file extension
        """
        try:
            return self.ext if self.ext.startswith(".") else "."+self.ext
        except AttributeError as e:
            raise AttributeError("You need to provide explicit file "
                                 "extension for user defined FileParser via its "
                                 "'ext' class attribute") from e
