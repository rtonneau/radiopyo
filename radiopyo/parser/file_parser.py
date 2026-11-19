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
from radiopyo.parser.exceptions import (
    ConfigError,
    IncludeNotFoundError,
    describe_chain,
)
from radiopyo.parser.utils import (
    PACKAGE_PREFIX,
    blank_config,
    combine_config,
    include_suggestions,
    locate_include,
)

# -------------------------- Type checking ONLY imports---------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.utils.sim_types import ConfigDict

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "FileParser",
    "list_available_parsers",
    "load_config",
    "normalize_includes",
    "resolve_includes",
]

# Accepted forms of the [includes] section, quoted in the errors.
INCLUDES_FORMS = (
    "Write either 'includes = [\"a.toml\", ...]' or a table keyed by "
    "extension: '[includes]' then 'toml = [\"a.toml\", ...]'."
)

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
    config: ConfigDict,
    _chain: tp.Tuple[pathlib.Path, ...] = (),
    *,
    source: str | pathlib.Path | None = None,
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
            (used to detect circular includes and named in the errors).
        source: the file 'config' was read from, if known: errors name
            it, and relative includes are also looked for next to it.

    Returns:
        ConfigDict: The configuration with all the includes merged in,
            and none left in its 'includes' section. 'config' itself is
            not modified. Without any include, 'config' is returned.

    Raises:
        IncludeNotFoundError: If an included file does not exist (it is a
            FileNotFoundError); the message lists the places tried and
            close matches.
        ConfigError: If a file includes itself (directly or not), for an
            include of an unknown file type or a malformed section.

    """
    if source is not None and not _chain:
        _chain = (pathlib.Path(source).resolve(),)
    includes = config["includes"]
    if len(includes) == 0:
        return config  # Do nothing
    available_parsers = list_available_parsers()
    base = _chain[-1].parent if _chain else None

    # Lowest precedence first: each included file overrides the previous ones
    merged = blank_config()

    for ext, entry in normalize_includes(includes, _chain):
        if ("." + ext) not in available_parsers:
            raise ConfigError(
                f"cannot include '{entry}': no parser for '.{ext}' files "
                f"(known: {', '.join(sorted(available_parsers))}). "
                + INCLUDES_FORMS,
                chain=_chain,
            )
        path = _found_include(entry, base, _chain)
        resolved = path.resolve()
        if resolved in _chain:
            cycle = " -> ".join(str(p) for p in (*_chain, resolved))
            raise ConfigError(f"Circular include: {cycle}", chain=_chain)
        logger.info(f"Resolving include: {path}")
        included = resolve_includes(
            available_parsers["." + ext](path).parse(), (*_chain, resolved)
        )
        merged = combine_config(merged, included)

    return combine_config(merged, config)


def load_config(path: str | pathlib.Path) -> ConfigDict:
    """
    The config of a file, includes resolved, read by the parser of its
    extension (see list_available_parsers).

    Raises:
        ConfigError: for an unknown extension, or any error of the file
            or of its includes (see resolve_includes).
    """
    parsers = list_available_parsers()
    suffix = pathlib.Path(path).suffix
    if suffix not in parsers:
        raise ConfigError(
            f"no parser for '{suffix}' files (known: {', '.join(sorted(parsers))})",
            chain=(pathlib.Path(path),),
        )
    return resolve_includes(parsers[suffix](path).parse(), source=path)


def normalize_includes(
    includes: tp.Any, chain: tp.Sequence[pathlib.Path] = ()
) -> tp.List[tp.Tuple[str, str]]:
    """
    The [includes] section as a list of (extension without dot, entry),
    in the order written. Accepted forms: a table {extension: entry or
    list of entries}, a list of entries or a single entry (the extension
    is then the suffix of each entry).

    Args:
        includes: the section, as read.
        chain: the files being read (named in the errors).

    Raises:
        ConfigError: for anything else (an entry that is not a string, an
            entry without suffix in the list form).
    """
    out: tp.List[tp.Tuple[str, str]] = []
    if isinstance(includes, (str, list, tuple)):
        entries = [includes] if isinstance(includes, str) else list(includes)
        for entry in _check_entries(entries, "includes", chain):
            ext = pathlib.PurePath(entry.strip()).suffix.lstrip(".")
            if not ext:
                raise ConfigError(
                    f"cannot include '{entry}': no file extension. " + INCLUDES_FORMS,
                    chain=chain,
                )
            out.append((ext, entry))
        return out
    if not isinstance(includes, tp.Mapping):
        raise ConfigError(
            f"invalid [includes] section ({type(includes).__name__}). "
            + INCLUDES_FORMS,
            chain=chain,
        )
    for ext, entries in includes.items():
        if isinstance(entries, str):
            entries = [entries]
        out += [(ext, e) for e in _check_entries(entries, f"includes.{ext}", chain)]
    return out


def _check_entries(
    entries: tp.Any, where: str, chain: tp.Sequence[pathlib.Path]
) -> tp.List[str]:
    """'entries' if it is a list of strings, ConfigError else."""
    if not isinstance(entries, (list, tuple)) or not all(
        isinstance(e, str) for e in entries
    ):
        raise ConfigError(
            f"'{where}' must hold file paths (strings), got {entries!r}. "
            + INCLUDES_FORMS,
            chain=chain,
        )
    return list(entries)


def _found_include(
    entry: str, base: pathlib.Path | None, chain: tp.Sequence[pathlib.Path]
) -> pathlib.Path:
    """The file designated by 'entry' (see locate_include), or raise."""
    entry = entry.strip()
    path, tried = locate_include(entry, base)
    if path is None:
        lines = [f"included file '{entry}' not found. Looked for:"]
        lines += [f"  - {p}" for p in tried]
        suggestions = include_suggestions(entry, base)
        if suggestions:
            quoted = ", ".join(f"'{s}'" for s in suggestions)
            lines.append(f"Did you mean: {quoted}?")
        lines.append(
            f"Bundled configs are written '{PACKAGE_PREFIX}<file>' (see "
            "radiopyo.bundled_configs()); a relative path is relative to the "
            "including file."
        )
        raise IncludeNotFoundError("\n".join(lines), chain=chain)
    relative = not (
        entry.startswith(PACKAGE_PREFIX) or pathlib.Path(entry).is_absolute()
    )
    if relative and (base is None or path != tried[0]):
        where = "the folder of the including file" if base is None else str(base)
        logger.warning(
            f"CONFIG::include '{entry}' (in {describe_chain(chain)}) is not in "
            f"{where}: using {path}, found in the working directory. Relative "
            "includes are looked for next to the including file first: move "
            "the file there or write its absolute path."
        )
    return path


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

    def anchor_includes(self, includes: tp.Any) -> tp.Dict[str, tp.List[str]]:
        """
        The 'includes' section of the file being parsed, with its relative
        paths anchored to the folder of the file: a relative path is
        looked for next to the file first, so that a config does not
        depend on the working directory. When there is nothing at that
        place the path is kept as written, and it will be looked for in
        the working directory (see resolve_includes).

        Entries referring to the package data ('radiopyo/...') and
        absolute paths are kept as written. The section is a table keyed
        by extension or a list of paths (see normalize_includes).

        Parsers should use it to fill the 'includes' section of the
        ConfigDict they produce.

        Args:
            includes: the section as read from the file.

        Returns:
            Dict[str, List[str]]: extensions (without dot) as keys.

        Raises:
            ConfigError: for a malformed section or an extension without
                parser, naming the file.
        """
        base = pathlib.Path(self.path).absolute().parent
        chain = (pathlib.Path(self.path).absolute(),)
        parsers = list_available_parsers()
        out: tp.Dict[str, tp.List[str]] = {}
        for ext, entry in normalize_includes(includes, chain):
            if "." + ext not in parsers:
                raise ConfigError(
                    f"cannot include '{entry}': no parser for '.{ext}' files "
                    f"(known: {', '.join(sorted(parsers))}). " + INCLUDES_FORMS,
                    chain=chain,
                )
            out.setdefault(ext, []).append(self._anchored(entry, base))
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
