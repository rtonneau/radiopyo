#!/usr/bin/env python3
"""
Exceptions raised while reading a configuration (config files, their
[includes] and the entries they define). Every one of them names where
the faulty entry comes from: a config file, the chain of files including
it ("main.toml -> a.toml"), or "<dict>" for a config given as a dict
(ADR 0012).

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 19/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import pathlib
import typing as tp

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.exceptions import RadioPyoError

# ------------------------- Type checking ONLY imports----------------------- #
# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "ConfigError",
    "ConfigSyntaxError",
    "IncludeNotFoundError",
    "describe_chain",
]

# Origin of an entry that was not read from a file.
DICT_ORIGIN = "<dict>"


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
def describe_chain(chain: tp.Sequence[str | pathlib.Path]) -> str:
    """'main.toml -> a.toml' for the chain of files (outermost first)."""
    return " -> ".join(str(p) for p in chain) if chain else DICT_ORIGIN


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #
class ConfigError(RadioPyoError, ValueError):
    """
    An error in a configuration. The message starts with where it comes
    from.

    Attributes:
        message: what is wrong, without the origin.
        origin: where the faulty entry comes from: the file (or the
            include chain, outermost first, joined by ' -> '), a
            description of the file of each field, or '<dict>'.
    """

    def __init__(
        self,
        message: str,
        *,
        origin: str | None = None,
        chain: tp.Sequence[str | pathlib.Path] = (),
    ):
        self.message = message
        self.origin = origin if origin is not None else describe_chain(chain)
        super().__init__(f"In config {self.origin}: {message}")


class IncludeNotFoundError(ConfigError, FileNotFoundError):
    """An [includes] entry designating no existing file."""


class ConfigSyntaxError(ConfigError):
    """A config file that cannot be read (e.g. invalid TOML)."""
