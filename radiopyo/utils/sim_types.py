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

import typing as tp

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.acid_base import ABDict
    from radiopyo.chemistry.k_reactions import KReactionDict
    from radiopyo.chemistry.michaelis import MichaelisMentenDict
    from radiopyo.chemistry.radiolytic import RadiolyticReactionDict

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "ConfigDict",
    "ReactionsDict",
    "ConcentrationDict",
    "default_ConfigDict",
]
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class ConcentrationDict(tp.TypedDict):
    fixed: tp.Dict[str, float]
    initial: tp.Dict[str, float]


class ReactionsDict(tp.TypedDict):
    radiolytic: RadiolyticReactionDict
    acid_base: tp.List[ABDict]
    k_reaction: tp.List[KReactionDict]
    enzymatic: tp.List[MichaelisMentenDict]


class ConfigDict(tp.TypedDict):
    includes: tp.Dict[str, tp.List[str]]
    bio_param: tp.Dict[str, float]
    concentrations: ConcentrationDict
    reactions: ReactionsDict
    beam: tp.Dict[str, tp.Dict[str, float]]


def default_ConfigDict() -> ConfigDict:
    """ """
    return {"includes": {},
            "bio_param": {},
            "concentrations": {
                "fixed": {},
                "initial": {},
    },
        "beam": {},
        "reactions": {
            "acid_base": [],
            "k_reaction": [],
            "enzymatic": [],
            "radiolytic": {},
    },
    }


def combine_config(root: ConfigDict, other: ConfigDict) -> ConfigDict:
    """" Keep root, erase entries from root when other has it and add entries from other
    not included in root"""
    # Easy part
    root["bio_param"] |= other["bio_param"]
    root["beam"] |= other["beam"]
    root["concentrations"]["fixed"] |= other["concentrations"]["fixed"]
    root["concentrations"]["initial"] |= other["concentrations"]["initial"]

    # Update reactions
    root["reactions"]["radiolytic"] |= other["reactions"]["radiolytic"]
    for _ in other["reactions"]["k_reaction"]:
        pass

    return root
