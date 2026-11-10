#!/usr/bin/env python3
"""
Tools used to assemble a simulation config out of several config files
(the ``[includes]`` mechanism): resolution of the include paths and merge
of two ConfigDict.

The merge only works on the dict structures (reactions are identified by
a key computed from their string, no reaction object is built), so that
this module does not depend on the simulation layer.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 26/07/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import logging
import pathlib
import typing as tp

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser import RPARSER
from radiopyo.physics.beam import BeamCollection
from radiopyo.utils import PACKAGE_DATA_PATH
from radiopyo.utils.sim_types import deep_update, default_ConfigDict

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.utils.sim_types import ConfigDict

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List[str] = [
    "PACKAGE_PREFIX",
    "ReactionKey",
    "blank_config",
    "combine_config",
    "include_path",
    "reaction_key",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# An include entry starting with this exact prefix refers to the data
# shipped with the package (PACKAGE_DATA_PATH).
PACKAGE_PREFIX = "radiopyo/"

# Identity of a reaction, whatever the order of its species: reactants
# {(label, stoichiometry)}, products {(label, stoichiometry, operation)},
# enzyme label (None if not enzymatic).
ReactionKey: tp.TypeAlias = tp.Tuple[
    tp.FrozenSet[tp.Tuple[str, int]],
    tp.FrozenSet[tp.Tuple[str, int, str]],
    tp.Optional[str],
]

# Sections of the config holding a list of reactions (one dict each)
_REACTION_LISTS: tp.Tuple[str, ...] = (
    "acid_base",
    "equilibrium",
    "k_reaction",
    "constant",
    "arrhenius",
    "enzymatic",
)

# An equilibrium is defined by 2 out of these 3 constants.
_EQUILIBRIUM_CONSTANTS: tp.Tuple[str, ...] = (
    "k_forward",
    "k_backward",
    "k_equilibrium",
)

# --------------------------------------------------------------------------- #
#                            DECORATOR DEFINITIONS                            #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


def blank_config() -> ConfigDict:
    """
    A ConfigDict with all its sections but nothing defined in it. In
    particular 'parameters' is empty (no default value): a config file
    only holds what the user wrote, defaults are applied when the
    simulation environment is built. Otherwise the defaults of an
    including file would erase the values of the files it includes.
    """
    config = default_ConfigDict()
    config["parameters"] = {}
    return config


def include_path(entry: str) -> pathlib.Path:
    """
    Absolute path designated by an entry of the [includes] section.

    - An entry starting with exactly 'radiopyo/' refers to the data
      shipped with the package: 'radiopyo/<file>' is
      PACKAGE_DATA_PATH / '<file>'.
    - Any other entry is a path, absolute or relative to the current
      working directory. (Relative paths written in a config file are
      anchored to that file when it is parsed, see
      FileParser.anchor_includes, so that the working directory is only
      a fallback.)
    """
    entry = entry.strip()
    if entry.startswith(PACKAGE_PREFIX):
        return PACKAGE_DATA_PATH / entry[len(PACKAGE_PREFIX) :]
    return pathlib.Path(entry).absolute()


def reaction_key(reaction: str) -> ReactionKey:
    """
    Identity of a reaction written as a string: the same for 'A + B -> C',
    'B + A -> C' and 'A+B->C', but not for 'C -> A + B'. Only the reaction
    string is parsed, no reaction object is built.

    Args:
        reaction: reaction string, as written in a config file.

    Raises:
        lark.LarkError: if the string is not a valid reaction.
    """
    generic = RPARSER.parse_reaction(reaction)
    enzyme = generic["enzyme"]
    return (
        frozenset(
            (label, stoi) for label, (_, stoi, _) in generic["reactants"].items()
        ),
        frozenset(
            (label, stoi, ope.name)
            for label, (_, stoi, ope) in generic["products"].items()
        ),
        None if enzyme is None else enzyme.label,
    )


def combine_config(root: ConfigDict, other: ConfigDict) -> ConfigDict:
    """
    Merge two config dicts; 'other' has the precedence over 'root'.
    Neither of them is modified, a new config is returned.

    - parameters, fixed/initial concentrations and radiolytic yields:
      the entries of 'other' replace those of 'root' with the same name,
      the others are added.
    - reactions (all the sections holding a list of reactions): two
      reactions are the same if they have the same species on each side
      whatever their order (see reaction_key), the arrow direction
      matters. A reaction of 'other' also in 'root' updates it in place
      field by field: the fields written in 'other' replace those of
      'root' (the reaction string too), the ones not written are kept.
      So a rate constant can be tweaked alone (e.g. only 'k_forward' of
      an acid/base keeps its 'k_backward' and 'pKa'). The exception is
      an equilibrium giving 2 or more of k_forward, k_backward and
      k_equilibrium: it is then fully defined by 'other', the constants
      of 'root' are dropped. New reactions are appended after those of
      'root'.
    - beam: see BeamCollection.merge. A 'default' beam of 'root' is
      dropped if 'other' defines any beam.
    - includes: not resolved here, the result has none.

    Args:
        root: lowest precedence config
        other: highest precedence config

    Returns:
        ConfigDict: the merged config.
    """
    out = _complete(root)
    new = _complete(other)

    out["includes"] = {}
    out["parameters"] |= new["parameters"]
    out["concentrations"]["fixed"] |= new["concentrations"]["fixed"]
    out["concentrations"]["initial"] |= new["concentrations"]["initial"]

    # Beams: make sure to keep only the latest beam definition if
    # multiple beam types
    root_beams = BeamCollection.from_config_dict(out["beam"], verbose=False)
    other_beams = BeamCollection.from_config_dict(new["beam"], verbose=False)
    if "default" in root_beams and len(other_beams) > 0:
        root_beams.remove_beam("default")
    out["beam"] = root_beams.merge(other_beams, verbose=False).as_dict()

    # Reactions
    out["reactions"]["radiolytic"] |= new["reactions"]["radiolytic"]
    reactions = tp.cast("tp.Dict[str, tp.Any]", out["reactions"])
    new_reactions = tp.cast("tp.Dict[str, tp.Any]", new["reactions"])
    for section in _REACTION_LISTS:
        reactions[section] = _merge_reaction_list(
            section,
            _as_list(reactions[section]),
            _as_list(new_reactions[section]),
        )
    return out


def _complete(config: ConfigDict) -> ConfigDict:
    """A deep copy of 'config' with all the sections of a ConfigDict."""
    return deep_update(blank_config(), copy.deepcopy(config))  # type: ignore[type-var]


def _as_list(reactions: tp.List[tp.Dict] | tp.Dict) -> tp.List[tp.Dict]:
    """A single reaction given as a table is a list of one reaction."""
    if isinstance(reactions, dict):
        return [reactions] if reactions else []
    return reactions


def _merge_reaction_list(
    section: str,
    base: tp.List[tp.Dict],
    override: tp.List[tp.Dict],
) -> tp.List[tp.Dict]:
    """
    Merge the reactions of 'override' into those of 'base' (see
    combine_config). The dicts of 'base' are updated in place.
    """
    position = {reaction_key(r["reaction"]): idx for idx, r in enumerate(base)}
    merged = list(base)
    for entry in override:
        idx = position.get(reaction_key(entry["reaction"]))
        if idx is None:  # a new reaction
            merged.append(entry)
            continue
        logger.debug(
            f"CONFIG::{section}: '{merged[idx]['reaction']}' redefined by "
            f"'{entry['reaction']}'"
        )
        if section == "equilibrium":
            given = [c for c in _EQUILIBRIUM_CONSTANTS if c in entry]
            if len(given) >= 2:  # fully defined by the new entry
                for constant in _EQUILIBRIUM_CONSTANTS:
                    merged[idx].pop(constant, None)
        merged[idx] |= entry
    return merged


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #
