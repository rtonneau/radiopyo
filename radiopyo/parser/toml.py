#!/usr/bin/env python3
"""
The utils.toml_parser submodule aims to load simulation config from
TOML files.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 09/05/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging
import re
import tomllib
import typing as tp
from contextlib import suppress

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser.exceptions import ConfigSyntaxError
from radiopyo.parser.file_parser import FileParser
from radiopyo.parser.utils import BACKSLASH_HINT, blank_config

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.utils.sim_types import ConfigDict

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = ["TOMLFileParser"]

# --------------------------------------------------------------------------- #
logger = logging.getLogger("radiopyo")

# Position given at the end of the tomllib messages
_POSITION = re.compile(r"\(at line (\d+), column (\d+)\)")


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
def syntax_error_message(err: Exception, raw: bytes) -> str:
    """
    The message of a TOML syntax error: the tomllib message, the faulty
    line, and a hint for a Windows path written with backslashes in a
    double-quoted string.
    """
    if isinstance(err, UnicodeDecodeError):
        return f"invalid TOML: the file is not UTF-8 text ({err})"
    message = f"invalid TOML: {err}"
    match = _POSITION.search(str(err))
    if match:
        lines = raw.decode("utf-8", errors="replace").splitlines()
        lineno, colno = int(match.group(1)), int(match.group(2))
        if 0 < lineno <= len(lines):
            line = lines[lineno - 1]
            message += f"\n  {lineno} | {line}\n  {' ' * len(str(lineno))} | "
            message += " " * max(colno - 1, 0) + "^"
            text = str(err)
            if "\\" in line and ("'\\'" in text or "hex value" in text or "escape" in text):
                message += "\n" + BACKSLASH_HINT
    return message


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #
class TOMLFileParser(FileParser):
    """
    Class parsing TOML files.
    """

    ext = ".toml"

    def parse(self) -> ConfigDict:
        """
        No Need for Lark grammar, "tomllib" works just fine =)

        The result only holds what is written in the file. In particular
        'parameters' has no default value: they are applied when the
        simulation environment is built, so that a file that does not
        mention a parameter does not erase the one of its includes.
        Relative paths of the [includes] section are anchored to the
        folder of the file (see FileParser.anchor_includes). Every
        reaction string is checked (see FileParser.check_reactions).

        Raises:
            ConfigSyntaxError: if the file is not valid TOML.
            ConfigError: for a malformed [includes] section or an
                invalid reaction.
        """
        with open(self.path, mode="rb") as fb:
            raw = fb.read()
        try:
            config = tomllib.loads(raw.decode("utf-8"))
        except (tomllib.TOMLDecodeError, UnicodeDecodeError) as err:
            raise ConfigSyntaxError(
                syntax_error_message(err, raw), chain=(self.absolute_path,)
            ) from err

        # Produces a bare dict with all fields initiated
        out = blank_config()
        out["includes"] = self.anchor_includes(config.get("includes", {}))
        out["parameters"] = config.get("parameters", {})
        out["beam"] = config.get("beam", {})
        with suppress(KeyError):  # in case "concentrations" does not exist
            out["concentrations"]["fixed"] = config["concentrations"].get("fixed", {})
            out["concentrations"]["initial"] = config["concentrations"].get(
                "initial", {}
            )

        with suppress(KeyError):  # in case "reactions" does not exist
            out["reactions"]["acid_base"] = config["reactions"].get("acid_base", [])
            out["reactions"]["radiolytic"] = config["reactions"].get("radiolytic", {})
            out["reactions"]["k_reaction"] = config["reactions"].get("k_reaction", [])
            out["reactions"]["enzymatic"] = config["reactions"].get("enzymatic", [])
            out["reactions"]["equilibrium"] = config["reactions"].get("equilibrium", [])
            out["reactions"]["arrhenius"] = config["reactions"].get("arrhenius", [])
            out["reactions"]["constant"] = config["reactions"].get("constant", [])
        self.check_reactions(out)
        return out
