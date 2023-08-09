#!/usr/bin/env python3
"""
The beam ron_parser contains all the logic for parse .ron input files.
TODO:
    - This module could use some heavy re-lifting
    - Want to keep RON file format?

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp
from functools import reduce
from pathlib import Path

import lark
from more_itertools import grouper

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.parser.file_parser import FileParser

# ------------------------------------------------------------------------------------ #
#                               TYPE CHECKING ONLY IMPORT                              #
# ------------------------------------------------------------------------------------ #
if tp.TYPE_CHECKING:
    from radiopyo.utils.sim_types import ConfigDict

# ------------------------------------------------------------------------------------ #
#                                      DECORATORS                                      #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class TreeToRon(lark.Transformer):
    """ """
    list = list
    tuple = tuple

    def null(self, items: tp.List) -> None:
        return None

    def true(self, items: tp.List) -> bool:
        return True

    def false(self, items: tp.List) -> bool:
        return False

    def start(self, items: tp.List) -> tp.Dict:
        return {k: v for (k, v) in grouper(items, n=2)}

    def CNAME(self, items: lark.Token) -> str:
        return str(items)

    def ESCAPED_STRING(self, items: lark.Token) -> str:
        return str(items)

    def key(self, items: tp.List) -> str:
        return str(items[0])

    def list_dict(self, items: tp.List) -> tp.Dict:
        return reduce(lambda a, b: a | b, items, {})

    def dict(self, items: tp.List) -> tp.Dict:
        return reduce(lambda a, b: a | b, items, {})

    def pair(self, items: tp.List) -> tp.Dict:
        return {items[0]: items[1]}

    def string(self, items: tp.List) -> str:
        return items[0].strip('"')

    def number(self, items: tp.List) -> float:
        return float(items[0])

    def SIGNED_NUMBER(self, items: lark.Token) -> float:
        return float(items)


class RonFileParser(FileParser):
    """
    Class parsing RON files.
    """
    ext = ".ron"
    items: dict

    def parse(self) -> ConfigDict:
        """ 
        """
        # Read the grammar
        with open(Path(__file__).parent/r"grammar//ron_grammar.lark", "r") as file:
            parser = lark.Lark(file.read())

        # Parses the config file
        with open(self.path, "r") as file:
            items = TreeToRon().transform(parser.parse(file.read()))

        # Output a ConfigDict
        return {"includes": items.get("includes", {}),
                "bio_param": items.get("bio_param", {}),

                "concentrations": {
                    "fixed": items.get("fixed_concentrations", {}),
                    "initial": items.get("initial_concentrations", {})
        },

            "reactions": {
                    "acid_base": items.get("acid_base", []),
                    "k_reaction": items.get("k_reactions", []),
                    "radiolytic": items.get("radiolytic", {}),
                    "enzymatic": items.get("enzymatic_reactions", []),
        },

            "beam": items.get("beam", {}),

        }
