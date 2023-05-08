#!/usr/bin/env python3
"""
The sim_types submodule contains class types hint (TypedDict) in order to help users to
interface to the module.
These classes, help to make file parser to read input.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 08/05/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import os
from typing import (
    TYPE_CHECKING,
    Dict,
    List,
    TypedDict,
)

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    import pathlib

    from .chemistry.acid_base import ABDict
    from .chemistry.k_reactions import KReactionDict
    from .chemistry.michaelis import MichaelisMentenDict
    from .chemistry.radiolytic import RadiolyticReactionDict

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "ConfigDict",
    "ReactionsDict",
    "ConcentrationDict",
    "FileParser",
]
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #

class ConcentrationDict(TypedDict):
    fixed: Dict[str, float]
    initial: Dict[str, float]

class ReactionsDict(TypedDict):
    radiolytic: RadiolyticReactionDict
    acid_base: List[ABDict]
    k_reaction: List[KReactionDict]
    enzymatic: List[MichaelisMentenDict]

class ConfigDict(TypedDict):
    bio_param: Dict[str, float]
    concentrations: ConcentrationDict
    reactions: ReactionsDict
    beam: Dict[str, Dict[str, float]]

class FileParser(object):
    """
    Generic Class parsing files.
    """
    path: str | pathlib.Path

    def __init__(self, path: str | pathlib.Path):
        self.path = path
        if not os.path.exists(self.path):
            raise FileNotFoundError(self.path)
    
    def parse(self) -> ConfigDict:
        """ 
        Method parsing the file (self.path) and returning a ConfigDict.

        Return:
            ConfigDict: TypedDictionary use to build the simulation Environment.
        """
        raise NotImplementedError