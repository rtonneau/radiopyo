#!/usr/bin/env python3
"""
Where each entry of a merged config comes from (ADR 0012).

A ConfigOrigins table sits next to a ConfigDict (never inside it: the
ConfigDict is also sent as is to the Rust backend). It records, for each
entry (a parameter, a fixed or initial concentration, a radiolytic yield,
a reaction, a beam) the file that defined it and, field by field, the
file that wrote each of its fields. It is merged with the same
precedence as the configs themselves (see combine_config), so that an
error raised while building the simulation from the merged config can
name the file(s) the faulty entry comes from (see config_entry).

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 19/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import contextlib
import logging
import typing as tp
from dataclasses import dataclass, field

from lark.exceptions import LarkError

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser.exceptions import DICT_ORIGIN, ConfigError
from radiopyo.parser.utils import (
    EQUILIBRIUM_CONSTANTS,
    REACTION_LISTS,
    combine_config,
    reaction_key,
)

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.utils.sim_types import ConfigDict

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "ConfigOrigins",
    "EntryOrigin",
    "combine_with_origins",
    "config_entry",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# Identity of an entry: (section, id). The section is e.g. "parameters",
# "concentrations.fixed", "reactions.radiolytic", "reactions.k_reaction"
# or "beam"; the id is the name, the reaction_key or the beam label.
EntryKey: tp.TypeAlias = tp.Tuple[str, tp.Hashable]

# Sections whose entries are plain "name = value" pairs
_NAMED_SECTIONS: tp.Tuple[tp.Tuple[str, ...], ...] = (
    ("parameters",),
    ("concentrations", "fixed"),
    ("concentrations", "initial"),
    ("reactions", "radiolytic"),
)


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
def _error_text(err: BaseException) -> str:
    """The message of 'err', with its type when the message is terse."""
    text = str(err).strip()
    if isinstance(err, KeyError):
        return f"missing key {text}"
    if not text:
        return type(err).__name__
    return text if isinstance(err, (ValueError, TypeError)) else (
        f"{type(err).__name__}: {text}"
    )


@contextlib.contextmanager
def config_entry(
    what: str, origin: str | tp.Callable[[], str]
) -> tp.Iterator[None]:
    """
    Turn an error raised in the block into a ConfigError naming the entry
    being built ('what', e.g. "[[reactions.k_reaction]] 'A + B -> C'") and
    where it comes from ('origin', see ConfigOrigins.describe; a callable
    is only called on error). The original error is its __cause__. A
    ConfigError is never wrapped twice.
    """
    try:
        yield
    except ConfigError:
        raise
    except Exception as err:
        where = origin() if callable(origin) else origin
        raise ConfigError(f"{what}: {_error_text(err)}", origin=where) from err


def combine_with_origins(
    root: ConfigDict,
    root_origins: ConfigOrigins,
    other: ConfigDict,
    other_origins: ConfigOrigins,
) -> tp.Tuple[ConfigDict, ConfigOrigins]:
    """
    combine_config(root, other) and the origins of the result. A beam
    that cannot be built while merging names its origins.
    """
    origins = root_origins.combined(other_origins)
    return combine_config(root, other, beam_origin=origins.beam), origins


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #
@dataclass
class EntryOrigin:
    """
    Origin of one config entry.

    Attributes:
        defined: the file (or '<dict>') where the entry first appears.
        fields: for each field of the entry, the file that wrote it last.
        kind: the beam type, for a beam (another type replaces it).
    """

    defined: str
    fields: tp.Dict[str, str] = field(default_factory=dict)
    kind: str = ""

    def describe(self) -> str:
        """'a.toml', or 'a.toml (k_value from main.toml)' when some fields
        were written by another file."""
        others: tp.Dict[str, tp.List[str]] = {}
        for name, origin in self.fields.items():
            if origin != self.defined:
                others.setdefault(origin, []).append(name)
        if not others:
            return self.defined
        written = "; ".join(
            f"{', '.join(names)} from {origin}" for origin, names in others.items()
        )
        return f"{self.defined} ({written})"


class ConfigOrigins:
    """
    Origins of the entries of a config (see the module docstring). An
    entry missing from the table comes from '<dict>'.
    """

    def __init__(self) -> None:
        self._entries: tp.Dict[EntryKey, EntryOrigin] = {}

    def __len__(self) -> int:
        return len(self._entries)

    def __contains__(self, key: EntryKey) -> bool:
        return key in self._entries

    def __getitem__(self, key: EntryKey) -> EntryOrigin:
        return self._entries[key]

    @classmethod
    def of(cls, config: ConfigDict, origin: str) -> ConfigOrigins:
        """
        Every entry of 'config' (a single file, includes not resolved)
        and each of its fields come from 'origin'.
        """
        out = cls()
        data = tp.cast("tp.Mapping[str, tp.Any]", config)
        for path in _NAMED_SECTIONS:
            section: tp.Any = data
            for part in path:
                section = section.get(part, {}) if isinstance(section, dict) else {}
            if not isinstance(section, dict):
                continue
            for name in section:
                key = (".".join(path), str(name).strip())
                out._entries[key] = EntryOrigin(origin)

        reactions = data.get("reactions", {})
        for name in REACTION_LISTS:
            entries = reactions.get(name, []) if isinstance(reactions, dict) else []
            for entry in [entries] if isinstance(entries, dict) else entries:
                rkey = _reaction_id(entry)
                if rkey is None:
                    continue
                out._entries[(f"reactions.{name}", rkey)] = EntryOrigin(
                    origin, dict.fromkeys(entry, origin)
                )

        beams = data.get("beam", {})
        for beam_type, entries in beams.items() if isinstance(beams, dict) else ():
            for entry in [entries] if isinstance(entries, dict) else entries:
                if not isinstance(entry, dict):
                    continue
                label = str(entry.get("label", "default"))
                out._entries[("beam", label)] = EntryOrigin(
                    origin, dict.fromkeys(entry, origin), kind=str(beam_type).lower()
                )
        return out

    def combined(self, other: ConfigOrigins) -> ConfigOrigins:
        """
        The origins of combine_config(root, other) when 'self' holds those
        of root: same precedence rules, field by field. Neither table is
        modified.
        """
        out = ConfigOrigins()
        out._entries = {
            key: EntryOrigin(e.defined, dict(e.fields), e.kind)
            for key, e in self._entries.items()
        }
        new_beams = any(key[0] == "beam" for key in other._entries)
        if new_beams:  # a 'default' beam of root is dropped (combine_config)
            out._entries.pop(("beam", "default"), None)

        for key, entry in other._entries.items():
            current = out._entries.get(key)
            if current is None or (key[0] == "beam" and current.kind != entry.kind):
                out._entries[key] = EntryOrigin(
                    entry.defined, dict(entry.fields), entry.kind
                )
                continue
            if not entry.fields:  # a name = value entry: replaced
                out._entries[key] = EntryOrigin(entry.defined)
                continue
            if key[0] == "reactions.equilibrium":
                given = [c for c in EQUILIBRIUM_CONSTANTS if c in entry.fields]
                if len(given) >= 2:  # fully redefined: old constants dropped
                    for constant in EQUILIBRIUM_CONSTANTS:
                        current.fields.pop(constant, None)
            current.fields |= entry.fields
        return out

    # ------------------------------ Queries -------------------------------- #
    def describe(self, section: str, name: tp.Hashable) -> str:
        """Origin of the entry 'name' of 'section' (see EntryOrigin)."""
        if isinstance(name, str):
            name = name.strip()
        entry = self._entries.get((section, name))
        return DICT_ORIGIN if entry is None else entry.describe()

    def reaction(self, section: str, reaction: tp.Any) -> str:
        """Origin of a reaction (its dict or its string) of a section of
        [reactions] ("k_reaction", "acid_base", ...)."""
        rkey = _reaction_id(reaction if isinstance(reaction, dict) else {
            "reaction": reaction
        })
        if rkey is None:
            return DICT_ORIGIN
        return self.describe(f"reactions.{section}", rkey)

    def beam(self, label: str) -> str:
        """Origin of the beam named 'label' ('default' when unnamed)."""
        return self.describe("beam", label)


def _reaction_id(entry: tp.Any) -> tp.Hashable | None:
    """The reaction_key of a reaction dict, None if it has none."""
    if not isinstance(entry, dict) or not isinstance(entry.get("reaction"), str):
        return None
    try:
        return reaction_key(entry["reaction"])
    except (LarkError, TypeError):  # reported by the builder
        return None
