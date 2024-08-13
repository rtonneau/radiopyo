#!/usr/bin/env python3
"""
This module contains the implementation of the abstract base class for file parsers
used in radiopyo.

The module provides the abstract base class :class:`FileParser` and the function
:func:`list_available_parsers` to list the available parsers.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 23/06/2023
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import logging
import os
import pathlib
import typing as tp
from abc import ABC, abstractmethod

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.parser.utils import combine_config

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.utils.sim_types import ConfigDict

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
    """
    Returns a dictionary containing all available parsers and their corresponding extensions.

    Args:
        None

    Returns:
        Dict[str, Type[FileParser]]: A dictionary with parser extensions as keys and parser classes as values.
    """
    return {cls.ext: cls for cls in FileParser.__subclasses__()}  # type: ignore [type-abstract]


def resolve_includes(config: ConfigDict) -> ConfigDict:
    """
    Resolves the includes in the given configuration dictionary.

    Args:
        config (ConfigDict): The configuration dictionary.

    Returns:
        ConfigDict: The updated configuration dictionary.

    Raises:
        FileNotFoundError: If the included file does not exist.

    """
    if len(config["includes"]) == 0:
        return config  # Do nothing
    available_parsers = list_available_parsers()
    includes = config["includes"].copy()
    config["includes"] = {}

    # Loop ovr all files in 'includes' section (sorted by extension)
    for ext, f_list in includes.items():
        # Check that it is a known extension
        if ("." + ext) not in available_parsers:
            logger.warning(f"Unable to include .{ext} file. No parser available")
            continue

        if not isinstance(f_list, (list, tuple)):
            f_list = [
                f_list,
            ]

        # Loop over all file with the same extension ('ext')
        for f in f_list:
            # If ask to include "built-in" configuration
            if f.startswith("radiopyo"):
                path = pathlib.Path(__file__).parent.parent
                path /= rf"data/{f.split('/')[-1].strip()}"
            else:
                path = pathlib.Path(f)

            path = path.absolute()
            if not path.exists():
                raise FileNotFoundError(path)
            logger.info(f"Resolving include: {path}")
            _ = resolve_includes(available_parsers["." + ext](path).parse())
            config = combine_config(_, config)

    return config


# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class FileParser(ABC):
    """
    FileParser is an abstract class designed to parse files into a dictionary.

    This class is meant to be subclassed by specific parsers for different file types.
    Each subclass should define `ext`, a class variable that specifies the file type
    the parser can handle, and `parse`, a method that parses the file and returns
    a dictionary.

    Parameters
    ----------
    path : str or pathlib.Path
        The path to the file to be parsed.

    Attributes
    ----------
    path : str or pathlib.Path
        The path to the file to be parsed.
    ext : classvar str
        The file extension that the parser can handle.

    Methods
    -------
    parse() -> dict
        Parse the file and return a dictionary.
    file_type() -> None
        Check that the file has the correct extension.

    """

    path: str | pathlib.Path
    ext: tp.ClassVar[str]

    def __init__(self, path: str | pathlib.Path):
        self.file_type()
        self.path = path
        if not os.path.exists(self.path):
            raise FileNotFoundError(self.path)

    @abstractmethod
    def parse(self) -> ConfigDict:
        """
        Method parsing the file (self.path) and returning a ConfigDict.

        Returns:
            ConfigDict: TypedDictionary use to build the simulation Environment.
        """

    def file_type(self) -> str:
        """
        raises:
            AttributeError: if ext class attribute is not defined by the user.

        returns:
            str: file extension
        """
        try:
            return self.ext if self.ext.startswith(".") else "." + self.ext
        except AttributeError as e:
            raise AttributeError(
                "You need to provide explicit file "
                "extension for user defined FileParser via its "
                "'ext' class attribute"
            ) from e
