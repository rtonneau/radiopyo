#!/usr/bin/env python3
"""
The utils.toml_parser submodule aims to load simulation config from TOML files.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 09/05/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import logging
import typing as tp
from contextlib import suppress

import tomllib

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.parser.file_parser import FileParser
from radiopyo.utils.sim_types import ConfigDict, default_ConfigDict

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = ["TOMLFileParser"]

# ------------------------------------------------------------------------------------ #
logger = logging.getLogger("radiopyo")

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class TOMLFileParser(FileParser):
    """
    Class parsing TOML files.
    """
    ext = ".toml"

    def parse(self) -> ConfigDict:
        """
        No Need for Lark grammar, "tomllib" works just fine =) 
        """
        with open(self.path, mode="rb") as fb:
            config = tomllib.load(fb)

        # Produces a bare dict with all fields initiated
        out = default_ConfigDict()
        out["includes"] = config.get("includes", {})
        out["bio_param"] = config.get("bio_param", {})
        out["beam"] = config.get("beam", {})
        with suppress(KeyError):  # in case "concentrations" does not exist
            out["concentrations"]["fixed"] = config["concentrations"].get(
                "fixed", {})
            out["concentrations"]["initial"] = config["concentrations"].get(
                "initial", {})

        with suppress(KeyError):  # in case "reactions" does not exist
            out["reactions"]["acid_base"] = config["reactions"].get(
                "acid_base", [])
            out["reactions"]["radiolytic"] = config["reactions"].get(
                "radiolytic", {})
            out["reactions"]["k_reaction"] = config["reactions"].get(
                "k_reaction", [])
            out["reactions"]["enzymatic"] = config["reactions"].get(
                "enzymatic", [])
        return out
