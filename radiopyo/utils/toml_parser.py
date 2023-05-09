#!/usr/bin/env python3
"""
The utils.toml_parser submodule aims to load simulation config from TOML files.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 09/05/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import logging

import tomllib

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.sim_types import ConfigDict, FileParser

# ----------------------------- Type checking ONLY imports---------------------------- #

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

    def parse(self) -> ConfigDict:
        """ 
        """
        with open(self.path, mode="rb") as fb:
            config = tomllib.load(fb)
        
        # Sanity checks
        if "concentrations" not in config:
            logger.warn("CONFIG::No concentration section found in config. It means, no"
                        " species with constant concentration and all initial"
                        " set to 0.")
            config["concentrations"] = {"fixed":{}, "initial":{}}

        if "reactions" not in config:
            logger.warning("CONFIG::No 'reactions' section in the provided "
                           "configuration")
            config["reactions"] = {"acid_base": [], "k_reaction": [],
                                   "enzymatic": [], "radiolytic": {}}
            
        return {"bio_param": config.get("bio_param", {}),
                
                "concentrations": {
                    "fixed": config["concentrations"].get("fixed", {}),
                    "initial": config["concentrations"].get("initial", {}),
                    },

                "reactions": {
                    "acid_base": config["reactions"].get("acid_base", []),
                    "k_reaction": config["reactions"].get("k_reaction", []),
                    "radiolytic": config["reactions"].get("radiolytic", {}),
                    "enzymatic": config["reactions"].get("enzymatic", []),
                    },

                "beam": config.get("beam", {}),
                }