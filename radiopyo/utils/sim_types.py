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

from radiopyo.chemistry.acid_base import ABCouple
from radiopyo.chemistry.k_reactions import KReaction
from radiopyo.chemistry.michaelis import MichaelisMenten

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
    #  => Radiolytic
    root["reactions"]["radiolytic"] |= other["reactions"]["radiolytic"]

    # => Acid Base
    r_root = {}
    for idx, _ in enumerate(root["reactions"]["acid_base"]):
        _ = ABCouple.from_dict(**_)
        r_root[str(_)] = (_, idx)
    to_add = []
    for _ in other["reactions"]["acid_base"]:
        ar = ABCouple.from_dict(**_)
        if str(ar) in r_root:
            if r_root[str(ar)][0].pKa != ar.pKa:
                idx = r_root[str(ar)][1]
                root["reactions"]["acid_base"][idx]["pKa"] = ar.pKa
            continue
        to_add.append(_.copy())
    root["reactions"]["acid_base"] += to_add

    # => k_reactions
    r_root = {}
    for idx, _ in enumerate(root["reactions"]["k_reaction"]):
        _ = KReaction.from_dict(**_)
        r_root[str(_)] = (_, idx)
    to_add = []
    for _ in other["reactions"]["k_reaction"]:
        kr = KReaction.from_dict(**_)
        if str(kr) in r_root:
            if r_root[str(kr)][0].k_value != kr.k_value:
                idx = r_root[str(kr)][1]
                root["reactions"]["k_reaction"][idx]["k_value"] = kr.k_value
            continue
        to_add.append(_.copy())
    root["reactions"]["k_reaction"] += to_add

    # => Enzymatic
    r_root = {}
    for idx, _ in enumerate(root["reactions"]["enzymatic"]):
        _ = MichaelisMenten.from_dict(**_)
        r_root[str(_)] = (_, idx)
    to_add = []
    for _ in other["reactions"]["enzymatic"]:
        er = MichaelisMenten.from_dict(**_)
        if str(er) in r_root:
            if r_root[str(er)][0].k_value != er.k_value and \
               r_root[str(er)][0].k_micha != er.k_micha:
                idx = r_root[str(er)][1]
                root["reactions"]["enzymatic"][idx]["k_value"] = er.k_value
                root["reactions"]["enzymatic"][idx]["k_micha"] = er.k_micha
            continue
        to_add.append(_.copy())
    root["reactions"]["enzymatic"] += to_add

    return root
