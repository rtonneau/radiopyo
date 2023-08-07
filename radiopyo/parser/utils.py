#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 26/07/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.acid_base import ABCouple
from radiopyo.physics.beam import BeamCollection
from radiopyo.simulation.environment_build import (
    kreaction_from_dict,
    michaelis_from_dict,
)

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.utils.sim_types import ConfigDict

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = []

# -------------------------------------- Logging ------------------------------------- #

# ------------------------------------------------------------------------------------ #
#                                DECORATOR DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


def combine_config(root: ConfigDict, other: ConfigDict) -> ConfigDict:
    """" 
    Function for config files assembly. It allows to define the environment in several
    independent config file. Convenient for reuse. 
    General behavior:
        Keep 'root', erase entries from 'root' when 'other' has it and add entries from
        'other' not included in 'root'
    """
    # Easy part
    root["bio_param"] |= other["bio_param"]
    root["concentrations"]["fixed"] |= other["concentrations"]["fixed"]
    root["concentrations"]["initial"] |= other["concentrations"]["initial"]

    # Update beam definition
    # Make sure to keep only the latest beam definition if multiple beam types
    root_beams = BeamCollection.from_config_dict(root["beam"], verbose=False)
    other_beams = BeamCollection.from_config_dict(other["beam"], verbose=False)
    if "default" in root_beams and len(other_beams) > 0:
        root_beams.remove_beam("default")

    root["beam"] = root_beams.merge(other_beams, verbose=False).as_dict()

    # Update reactions
    #  => Radiolytic
    root["reactions"]["radiolytic"] |= other["reactions"]["radiolytic"]

    # For the following, Reactions are built from Dict and then their unique identifier
    # (string version) are used to compare them. First, all reaction from 'root' are
    # read and stored in a Dict whose keys are the strings. Then each entries of 'other'
    # are parsed to extract their string. if it exists in the dict, their constants are
    # compared and updated accordingly.

    ## => Acid Base ##
    r_root = {}
    for idx, _ in enumerate(root["reactions"]["acid_base"]):
        _ = ABCouple.from_dict(**_)
        r_root[str(_)] = (_, idx)
    to_add = []  # List of all reaction to add to root afterward
    for _ in other["reactions"]["acid_base"]:
        ar = ABCouple.from_dict(**_)
        if str(ar) in r_root:
            idx = r_root[str(ar)][1]
            root["reactions"]["acid_base"][idx]["pKa"] = ar.pKa
            continue
        to_add.append(_.copy())
    root["reactions"]["acid_base"] += to_add

    ## => k_reactions ##
    # First, let's make a set of the kreaction presents in root. Each one should have a
    # unique ID ==> Make kReaction and then use the str version of it!
    r_root = {}
    for idx, _ in enumerate(root["reactions"]["k_reaction"]):
        _ = kreaction_from_dict(_, sim_sp=None)[0][0]
        r_root[str(_)] = (_, idx)

    # Loop over all KReaction in other and check whether it is also in root or not.
    to_add = []  # List of all reaction to add to root afterward
    for _ in other["reactions"]["k_reaction"]:
        kr = kreaction_from_dict(_, sim_sp=None)[0][0]
        # Same reaction in root => just make sur the kvalue comes from 'other'
        if str(kr) in r_root:
            idx = r_root[str(kr)][1]
            # raw_k_value is necessary to avoid multiple divide for reactions like
            # xA -> B + C (in this case, xk should be provided)
            root["reactions"]["k_reaction"][idx]["k_value"] = kr.raw_k_value()
            continue
        # KReaction from 'other' not in 'root'. Make a copy of the dict version.
        to_add.append(_.copy())
    # All shared reactions are updated at this point. Now just add all new reactions
    # from 'other' --> stored in 'to_add'
    root["reactions"]["k_reaction"] += to_add

    ## => Enzymatic ##
    r_root = {}
    for idx, _ in enumerate(root["reactions"]["enzymatic"]):
        _ = michaelis_from_dict(_, sim_sp=None)[0][0]
        r_root[str(_)] = (_, idx)
    to_add = []  # List of all reaction to add to root afterward
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


def dict_compare(d1: tp.Dict, d2: tp.Dict) -> tp.Dict:
    """ 
    Dictionary comparison, insightful.
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
    return {"added": added,
            "removed": removed,
            "modified": modified,
            "same": same}
