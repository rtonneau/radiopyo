#!/usr/bin/env python3
"""
This module contains the implementation of the abstract base class for
file parsers used in radiopyo.

The module provides the abstract base class :class:`FileParser` and
the function
:func:`list_available_parsers` to list the available parsers.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 23/06/2023
"""
# ----------------------------- External imports----------------------------- #
from __future__ import annotations

import logging
import os
import pathlib
import typing as tp
from abc import ABC, abstractmethod

# ------------------------------ Internal imports --------------------------- #
from radiopyo.parser.utils import (
    PACKAGE_PREFIX,
    blank_config,
    combine_config,
    include_path,
)

# -------------------------- Type checking ONLY imports---------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.utils.sim_types import ConfigDict

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "FileParser",
    "list_available_parsers",
    "resolve_includes",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                           DECORATOR DEFINITIONS                             #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


def list_available_parsers() -> tp.Dict[str, tp.Type[FileParser]]:
    """
    Returns a dictionary containing all available parsers and their
    corresponding extensions.

    Args:
        None

    Returns:
        Dict[str, Type[FileParser]]: A dictionary with parser extensions
            as keys and parser classes as values.
    """
    return {cls.ext: cls for cls in FileParser.__subclasses__()}  # type: ignore [type-abstract]


def resolve_includes(
    config: ConfigDict, _chain: tp.Tuple[pathlib.Path, ...] = ()
) -> ConfigDict:
    """
    Resolves the includes in the given configuration dictionary.

    The included files are loaded (recursively) and merged in the order
    they are listed: when two of them define the same thing, the last one
    wins. The config itself has the highest precedence, it overrides
    everything coming from its includes (see combine_config for the
    details of the merge).

    Args:
        config (ConfigDict): The configuration dictionary.
        _chain (Tuple[pathlib.Path, ...]): internal, the files being resolved
            (used to detect circular includes).

    Returns:
        ConfigDict: The configuration with all the includes merged in,
            and none left in its 'includes' section. 'config' itself is
            not modified. Without any include, 'config' is returned.

    Raises:
        FileNotFoundError: If the included file does not exist.
        ValueError: If a file includes itself (directly or not).

    """
    includes = config["includes"]
    if len(includes) == 0:
        return config  # Do nothing
    available_parsers = list_available_parsers()

    # Lowest precedence first: each included file overrides the previous ones
    merged = blank_config()

    # Loop over all files in 'includes' section (sorted by extension)
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
            path = include_path(f)
            if not path.exists():
                raise FileNotFoundError(path)
            resolved = path.resolve()
            if resolved in _chain:
                cycle = " -> ".join(str(p) for p in (*_chain, resolved))
                raise ValueError(f"Circular include: {cycle}")
            logger.info(f"Resolving include: {path}")
            included = resolve_includes(
                available_parsers["." + ext](path).parse(), (*_chain, resolved)
            )
            merged = combine_config(merged, included)

    return combine_config(merged, config)


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #


class FileParser(ABC):
    """
    FileParser is an abstract class designed to parse files into a
    dictionary.

    This class is meant to be subclassed by specific parsers for
    different file types. Each subclass should define `ext`, a class
    variable that specifies the file type the parser can handle, and
    `parse`, a method that parses the file and returns a dictionary.

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
            ConfigDict: TypedDictionary use to build the simulation
                Environment.
        """

    def anchor_includes(
        self,
        includes: tp.Mapping[str, str | tp.Sequence[str]],
    ) -> tp.Dict[str, tp.List[str]]:
        """
        The 'includes' section of the file being parsed, with its relative
        paths anchored to the folder of the file: a relative path is
        looked for next to the file first, so that a config does not
        depend on the working directory. When there is nothing at that
        place the path is kept as written, and it will be looked for in
        the working directory (see resolve_includes).

        Entries referring to the package data ('radiopyo/...') and
        absolute paths are kept as written.

        Parsers should use it to fill the 'includes' section of the
        ConfigDict they produce.

        Args:
            includes: the section as read from the file, extensions as
                keys, one path or a list of paths as values.

        Returns:
            Dict[str, List[str]]
        """
        base = pathlib.Path(self.path).absolute().parent
        out: tp.Dict[str, tp.List[str]] = {}
        for ext, entries in includes.items():
            if isinstance(entries, str):
                entries = [entries]
            out[ext] = [self._anchored(entry, base) for entry in entries]
        return out

    @staticmethod
    def _anchored(entry: str, base: pathlib.Path) -> str:
        """'entry' as the path next to 'base' if it exists, 'entry' else."""
        entry = entry.strip()
        if entry.startswith(PACKAGE_PREFIX) or pathlib.Path(entry).is_absolute():
            return entry
        candidate = base / entry
        return str(candidate) if candidate.exists() else entry

    def file_type(self) -> str:
        """
        raises:
            AttributeError: if ext class attribute is not defined by
                the user.

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
