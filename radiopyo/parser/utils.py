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
import difflib
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
    "bundled_configs",
    "combine_config",
    "disabled_constant",
    "drop_disabled_reactions",
    "include_path",
    "include_suggestions",
    "locate_include",
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

# Scalar rate constants of each reaction section: one of them equal to 0
# disables the reaction (the whole pair for acid_base and equilibrium),
# see drop_disabled_reactions and ADR 0011.
_RATE_CONSTANTS: tp.Dict[str, tp.Tuple[str, ...]] = {
    "acid_base": ("k_forward", "k_backward"),
    "equilibrium": _EQUILIBRIUM_CONSTANTS,
    "k_reaction": ("k_value",),
    "constant": ("k_value",),
    "arrhenius": ("A",),
    "enzymatic": ("k_value",),
}

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


def bundled_configs() -> tp.Dict[str, pathlib.Path]:
    """
    The TOML files shipped with the package (in PACKAGE_DATA_PATH, sub
    folders included), as {include entry: absolute path}, sorted by
    entry. The include entry is what to write in an [includes] section,
    e.g. 'radiopyo/advanced_water_radiolysis.toml'.
    """
    found = {
        PACKAGE_PREFIX + f.relative_to(PACKAGE_DATA_PATH).as_posix(): f.absolute()
        for f in PACKAGE_DATA_PATH.rglob("*.toml")
    }
    return dict(sorted(found.items()))


def locate_include(
    entry: str, base: pathlib.Path | None = None
) -> tp.Tuple[pathlib.Path | None, tp.List[pathlib.Path]]:
    """
    Where an [includes] entry is looked for, in order, and the first of
    those places holding a file.

    - 'radiopyo/<file>': PACKAGE_DATA_PATH / '<file>' only.
    - an absolute path: that path only.
    - a relative path: next to the including file ('base', its folder,
      when known), then in the current working directory.

    Args:
        entry: the entry, as written in the config.
        base: folder of the including file, None if unknown.

    Returns:
        (the path found or None, the list of the places tried)
    """
    entry = entry.strip()
    if entry.startswith(PACKAGE_PREFIX) or pathlib.Path(entry).is_absolute():
        tried = [include_path(entry)]
    else:
        tried = [] if base is None else [(base / entry).absolute()]
        cwd = include_path(entry)
        if cwd not in tried:
            tried.append(cwd)
    found = next((p for p in tried if p.is_file()), None)
    return found, tried


def include_suggestions(entry: str, base: pathlib.Path | None = None) -> tp.List[str]:
    """
    Entries close to 'entry' that would designate an existing file: the
    bundled configs (whatever the folder written before their name, e.g.
    'radiopyo/data/x.toml' or a bare 'x.toml' gives 'radiopyo/x.toml')
    and the config files next to the including file (folder 'base').
    """
    entry = entry.strip().replace("\\", "/")
    name = pathlib.PurePosixPath(entry).name
    bundled = bundled_configs()
    out = [key for key, path in bundled.items() if path.name == name]
    if out:  # the very file name of a bundled config: nothing else to say
        return out

    candidates = list(bundled)
    if base is not None and base.is_dir():
        candidates += sorted(
            f.relative_to(base).as_posix()
            for f in base.rglob("*.toml")
            if len(f.relative_to(base).parts) <= 3
        )
    for match in difflib.get_close_matches(entry, candidates, n=3, cutoff=0.6):
        if match not in out:
            out.append(match)
    # A close name in another folder (e.g. a typo in a bundled file name)
    names = {pathlib.PurePosixPath(c).name: c for c in reversed(candidates)}
    for match in difflib.get_close_matches(name, list(names), n=3, cutoff=0.75):
        if names[match] not in out:
            out.append(names[match])
    return out[:5]


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


def disabled_constant(section: str, reaction: tp.Mapping[str, tp.Any]) -> str | None:
    """
    Name of the first scalar rate constant of `reaction` (a reaction dict
    of the config section `section`) equal to 0, None if there is none.
    A polynomial constant (a list) is never zero: its coefficients are
    those of log10(k).
    """
    for name in _RATE_CONSTANTS.get(section, ()):
        value = reaction.get(name)
        is_scalar = isinstance(value, (int, float)) and not isinstance(value, bool)
        if is_scalar and value == 0:
            return name
    return None


def drop_disabled_reactions(config: ConfigDict) -> ConfigDict:
    """
    The config without its disabled reactions: a reaction whose scalar
    rate constant is 0 (k_value, Arrhenius A, or k_forward / k_backward /
    k_equilibrium of an acid/base or equilibrium pair, which drops the
    whole pair) is removed, with one warning each (ADR 0011). Applied on
    the merged config, so a reaction of an included file can be disabled
    by overriding its rate constant with 0. 'config' is not modified;
    calling it again on the result removes nothing.
    """
    reactions = tp.cast("tp.Dict[str, tp.Any]", config.get("reactions", {}))
    kept: tp.Dict[str, tp.Any] = {}
    dropped = False
    for section, entries in reactions.items():
        if section not in _RATE_CONSTANTS:
            kept[section] = entries
            continue
        out = []
        for entry in _as_list(entries):
            name = disabled_constant(section, entry)
            if name is None:
                out.append(entry)
                continue
            dropped = True
            what = "pair " if section in ("acid_base", "equilibrium") else ""
            logger.warning(
                f"CONFIG::[[reactions.{section}]] '{entry.get('reaction')}' has "
                f"{name} = 0: this reaction {what}is disabled (removed from the "
                "simulation)."
            )
        kept[section] = out
    if not dropped:
        return config
    new = tp.cast("tp.Dict[str, tp.Any]", dict(config))
    new["reactions"] = kept
    return tp.cast("ConfigDict", new)


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #
