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
from radiopyo.sim_env import kreaction_from_dict, michaelis_from_dict

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "KwargsRun",
    "ConfigDict",
    "ReactionsDict",
    "ConcentrationDict",
    "default_ConfigDict",
]
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class KwargsRun(tp.TypedDict):
    t_span: tp.Tuple[float, float]
    t_eval: tp.Optional[NDArray[np.float64]]
    use_jac: bool
    method: str
    atol: tp.Optional[float]
    rtol: tp.Optional[float]
    max_step: tp.Optional[float]
    y0: tp.Optional[NDArray[np.float64]]


class ConcentrationDict(tp.TypedDict):
    fixed: tp.Dict[str, float]
    initial: tp.Dict[str, float]


class ReactionsDict(tp.TypedDict):
    radiolytic: tp.Dict
    acid_base: tp.List[tp.Dict]
    k_reaction: tp.List[tp.Dict]
    enzymatic: tp.List[tp.Dict]


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
            idx = r_root[str(ar)][1]
            root["reactions"]["acid_base"][idx]["pKa"] = ar.pKa
            continue
        to_add.append(_.copy())
    root["reactions"]["acid_base"] += to_add

    # => k_reactions
    r_root = {}
    for idx, _ in enumerate(root["reactions"]["k_reaction"]):
        _ = kreaction_from_dict(_, sim_sp=None)[0][0]
        r_root[str(_)] = (_, idx)
    to_add = []
    for _ in other["reactions"]["k_reaction"]:
        kr = kreaction_from_dict(_, sim_sp=None)[0][0]
        if str(kr) in r_root:
            idx = r_root[str(kr)][1]
            root["reactions"]["k_reaction"][idx]["k_value"] = kr.raw_k_value()
            continue
    to_add.append(_.copy())
    root["reactions"]["k_reaction"] += to_add

    # => Enzymatic
    r_root = {}
    for idx, _ in enumerate(root["reactions"]["enzymatic"]):
        _ = michaelis_from_dict(_, sim_sp=None)[0][0]
        r_root[str(_)] = (_, idx)
    to_add = []
    for _ in other["reactions"]["enzymatic"]:
        er = michaelis_from_dict(_, sim_sp=None)[0][0]
        if str(er) in r_root:
            idx = r_root[str(er)][1]
            root["reactions"]["enzymatic"][idx]["k_value"] = er.raw_k_value()
            root["reactions"]["enzymatic"][idx]["k_micha"] = er.k_micha
            continue
        to_add.append(_.copy())
    root["reactions"]["enzymatic"] += to_add

    return root


def dict_compare(d1: tp.Dict, d2: tp.Dict) -> tp.Tuple[tp.Set, tp.Set, tp.Dict, tp.Set]:
    """ 
    From:
    https://stackoverflow.com/questions/4527942
    """
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    removed = d1_keys - d2_keys
    added = d2_keys - d1_keys
    modified = {o: (d1[o], d2[o]) for o in shared_keys if d1[o] != d2[o]}
    same = {o for o in shared_keys if d1[o] == d2[o]}
    return added, removed, modified, same
