#!/usr/bin/env python3
"""
The beam_collection submodule contains the container used to store
several beam definitions (only one being used at a time).

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging
import typing as tp
from contextlib import suppress

# ----------------------------- Internal imports ---------------------------- #
# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from .beam import ParticleBeam

# -------------------------- Expose importable stuffs ----------------------- #
__all__ = [
    "BeamCollection",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class BeamCollection(object):
    """
    Container class for efficient beam storage so that multiple beams
    can be define (even though only one at a time can be use).
    """

    # Beam definitions are stored in this dict
    _beam: tp.Dict[str, ParticleBeam]
    # Keyword to access the currently used beam from _beam
    _current: str
    # LET to use with the corresponding beams
    _let: tp.Dict[str, float | None]

    def __init__(self) -> None:
        self._beam = {}
        self._let = {}
        self._current = ""

    def __contains__(self, key: str) -> bool:
        return key in self._beam

    def __iter__(self) -> tp.Iterator[tp.Tuple[str, ParticleBeam, float | None]]:
        for name, beam in self._beam.items():
            yield name, beam, self._let[name]

    def __len__(self) -> int:
        return len(self._beam)

    def __getitem__(self, key: str) -> tp.Tuple[ParticleBeam, float | None]:
        return self._beam[key], self._let[key]

    def __str__(self) -> str:
        out = ["BeamCollection:\n"]
        for name, beam, let in self:
            out.append(f"Beam {name} [{let}keV/µm] -> {beam}\n")
        return "".join(out)

    @property
    def current(self) -> ParticleBeam:
        return self._beam[self._current]

    @property
    def current_name(self) -> str:
        return self._current

    @property
    def current_LET(self) -> float:
        if self._let[self._current] is not None:
            return self._let[self._current]  # type: ignore [return-value]
        raise ValueError("LET not defined for current beam")

    def remove_beam(self, key: str) -> None:
        """Remove a beam entry.
        Args:
            key: str
                name of the beam to delete
        Raises:
            KeyError
                No beam with name: 'key'
        """
        if key not in self:
            raise KeyError(f"No beam definition named: {key}")
        del self._let[key]
        del self._beam[key]

    def update_let(self, let: float, key: tp.Optional[str] = None) -> None:
        """Modifies the LET of a beam
        Args:
            let: float
                LET value
            key: str, optional (default=None)
                name of the beam to update. Default is the current one.
        Raises:
            KeyError
                Beam with name 'key' not found
        """
        if key is None:
            key = self._current
        elif key not in self:
            raise KeyError(f"No beam definition named: {key}")
        self._let[key] = let

    def add_beam(
        self,
        beam: ParticleBeam,
        key: str = "default",
        let: tp.Optional[float] = None,
        use: bool = False,
    ) -> None:
        """
        Pushes a new beam into the container.
        Args:
            beam: ParticleBeam
                beam to push in
            key: str, optional (default="default")
                Name of the beam
            let: float, optional (default=None)
                LET of the beam
            use: bool, optional (default=False)
                Defines the newly pushed beam as the current one
        """
        if key is None or len(key) == 0:
            key = "default"
        # First remove any previous definition with the same 'key'
        if key in self:
            self.remove_beam(key)
        self._beam[key] = beam
        self._let[key] = let
        if use or len(self._current) == 0:
            self._current = key

    def use_beam(self, key: str) -> tp.Tuple[ParticleBeam, float | None]:
        """
        Defines which beam to use
        Args:
            key: str
                Beam name
        Raises:
            KeyError
                Beam with name 'key' not found
        Returns:
            ParticleBeam
            float
                The beam LET
        """
        if key not in self._beam:
            raise KeyError(f"No beam definition named: {key}")
        self._current = key
        return self.current, self.current_LET

    @staticmethod
    def build_beam_from_dict(
        beam_type: str,
        config: tp.Dict[str, float],
        LET: tp.Optional[float] = None,
        verbose: bool = True,
    ) -> tp.Tuple[str, ParticleBeam, float | None]:
        """
        Build a beam out of a dictionary. The caller's dictionary is
        left untouched.

        Raises:
            NotImplementedError
                beam_type not recognized
        """
        # Imported here: beam.py re-exports BeamCollection (circular import)
        from .beam import BEAM_TYPES

        config = dict(config)
        # if LET not provided, check if present in config.
        if LET is None:
            LET = config.get("LET")
        elif "LET" not in config:
            config["LET"] = LET  # So that LET is displayed in the bottom log

        label: str
        if "label" not in config:
            msg = (
                "BEAM:: Assigning 'default' name to beam definition "
                f"{beam_type} -> {config}"
            )
            if verbose:
                logger.info(msg)
            label = "default"
        else:
            label = str(config["label"])
            msg = f"BEAM:: Found beam '{label}' as {beam_type} -> {config}"
            if verbose:
                logger.info(msg)

        try:
            beam_class = BEAM_TYPES[beam_type.lower()]
        except KeyError:
            raise NotImplementedError(f"Unknown beam type: {beam_type}") from None
        beam = beam_class.from_dict(config)
        return label, beam, LET

    def add_beam_from_dict(
        self,
        beam_type: str,
        config: tp.Dict[str, float],
        LET: tp.Optional[float] = None,
        verbose: bool = True,
    ) -> None:
        """
        Add new beam defined in a dictionary.
        Args:
            beam_type: str
                type of the beam to add
            config: dict
                Dictionary containing values defining the beam
            LET: float, optional
                LET of the beam
            verbose: bool, optional (default=True)
                Enable or disable logging

        Raises:
            NotImplementedError
                beam_type not recognized
        """
        label, beam, LET = BeamCollection.build_beam_from_dict(
            beam_type, config, LET, verbose
        )

        if verbose and label in self:
            msg = (
                f"BEAM:: Beam definition: {label} already defined. "
                "Erasing previous definition"
            )
            logger.info(msg)

        self.add_beam(beam=beam, key=label, let=LET)

    @classmethod
    def from_config_dict(
        cls,
        param: tp.Dict[str, tp.List[tp.Dict]],
        verbose: bool = True,
    ) -> BeamCollection:
        """
        Create a BeamCollection out of a Dict containing multiple beam
        definitions.Param dictionary keywords should be valid beam type
        names. Values can be either be a single Dict-like beam or a
        list of it.

        Args:
            param: Dict
                Dict structure whose keywords are valid beam type names.
            verbose: bool, optional (default=True)
        """
        beams = cls()

        for beam_type, beam_config in param.items():
            # If only one beam, not in a list, put it in a list!
            if isinstance(beam_config, dict):
                beam_config = [
                    beam_config,
                ]
            for config in beam_config:
                config = dict(config)  # never modify the caller's dict
                LET: tp.Optional[float] = None
                with suppress(KeyError):
                    LET = config.pop("LET")

                beams.add_beam_from_dict(beam_type, config, LET, verbose)
        return beams

    def as_dict(self) -> tp.Dict[str, tp.List[tp.Dict]]:
        """
        Extract all the beam in a dictionary whose structure comply
        with the method 'from_config_dict'
        """
        out: tp.Dict[str, tp.List[tp.Dict]] = {}
        for name, beam, LET in self:
            x = beam.as_dict()
            if name != "default":
                x["label"] = name
            if LET is not None:
                x["LET"] = self._let[name]
            beam_type = beam.beam_type().lower()
            if beam_type not in out:
                out[beam_type] = [
                    x,
                ]
            else:
                out[beam_type].append(x)
        return out

    def merge(
        self,
        other: BeamCollection,
        verbose: bool = True,
    ) -> BeamCollection:
        """
        Merge two BeamCollection. Do not modify the current definition
        but return a new one.

        Args:
            other: BeamCollection
                collection to merge with
            verbose: bool, optional (default=True)
                Flag for logging

        Returns:
            BeamCollection
        """
        s: tp.Dict[str, tp.Tuple[tp.Dict, str]] = {}
        for name, _beam, let in self:
            d = _beam.as_dict()
            if let is not None:
                d["LET"] = let
            d["label"] = name
            s[name] = (d, _beam.beam_type().lower())

        for beam_type, beams in other.as_dict().items():
            for beam in beams:
                name = beam.get("label", "default")
                if name not in s:
                    s[name] = (beam, beam_type)
                    continue
                if s[name][1] == beam_type:
                    s[name][0].update(beam)
                    continue
                s[name] = (beam, beam_type)

        out = BeamCollection()
        for _, (beam, beam_type) in s.items():
            out.add_beam_from_dict(beam_type, beam, verbose=verbose)

        return out
