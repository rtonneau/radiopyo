#!/usr/bin/env python3
"""
Small parser to be able to parse reaction def from strings like:
xA + yB -> iC + kD
xA + yB --enzyme>> iC + kD
Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 29/06/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging
import typing as tp
import warnings
from pathlib import Path

import lark
import numpy as np

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.basic_reaction import GenericReaction, Operation
from radiopyo.chemistry.species import (
    Atom,
    BasicSpecies,
    HydratedElectron,
    Molecule,
)
from radiopyo.chemistry.utils.atomic_weights import Atomix

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    pass

# ------------------------- Expose importable stuffs ------------------------ #
__all__: tp.List = []

# ----------------------------------- Logging ------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                              DECORATOR DEFINITIONS                          #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #


class TreeToReaction(lark.Transformer):
    """Lark Class handling the conversion.
    See: https://lark-parser.readthedocs.io/en/latest/visitors.html
    For some doc and examples.
    """

    def start(
        self,
        items: tp.List,
    ) -> GenericReaction | Atom | Molecule:
        return items[0]

    def constant_reaction(self, items: tp.List) -> GenericReaction:
        """ """
        return GenericReaction(
            reactants={},
            products=items[0],
            constants={},
            enzyme=None,
        )

    def equilibrium_reaction(self, items: tp.List) -> GenericReaction:
        return GenericReaction(
            reactants=items[0],
            products=items[1],
            constants={},
            enzyme=None,
        )

    def enzymatic_reaction(self, items: tp.List) -> GenericReaction:
        return GenericReaction(
            reactants={items[0][0].label: items[0]},
            products=items[2],
            enzyme=items[1],
            constants={},
        )

    def k_reaction(self, items: tp.List) -> GenericReaction:
        return GenericReaction(
            reactants=items[0],
            products=items[1],
            constants={},
            enzyme=None,
        )

    def reactants(self, items: tp.List) -> tp.Dict:
        out: tp.Dict[str, tp.Tuple[BasicSpecies, int, Operation]] = {}
        for species, stoi, ope in items:
            label = species.label
            out[label] = (species, out.get(label, (0, 0))[1] + stoi, ope)
        return out

    def products(self, items: tp.List) -> tp.Dict:
        out: tp.Dict[str, tp.Tuple[BasicSpecies, int, Operation]] = {}
        for species, stoi, ope in items:
            label = species.label
            sign = 1 if ope is Operation.ADD else -1
            out[label] = (
                species,
                out.get(label, (0, 0))[1] + stoi * sign,  # Summing stoi
                Operation.ADD,  # by default, change it later.
            )
        for key, value in out.items():
            # if the sum of all stoichimetric coeff is < 0, then it is a subtraction.
            if value[1] < 0:
                out[key] = (value[0], abs(value[1]), Operation.SUBTRACT)
        # Remove all species with stoichiometry = 0
        return {key: value for key, value in out.items() if value[1] > 0}

    def signed_species(self, items: tp.List) -> tp.Tuple[BasicSpecies, int, Operation]:
        return (items[1][0], items[1][1], items[0])

    def species(self, items: tp.List) -> tp.Tuple[BasicSpecies, int, Operation]:
        """
        Parse a species from the given items.

        Args:
            items (tp.List): List of items to parse.

        Returns:
            tp.Tuple[BasicSpecies, int]: Parsed species and its stoichiometry.
        """
        if len(items) > 1:
            stoi = int(items[0])
            kwargs = items[1]
        else:
            stoi = 1
            kwargs = items[0]

        if isinstance(kwargs, str):
            return (BasicSpecies(kwargs), stoi, Operation.ADD)
        if "name" in kwargs:
            if kwargs["name"] == "e_aq":
                return (HydratedElectron(), stoi, Operation.ADD)
            return (Atom(**kwargs), stoi, Operation.ADD)
        if "atoms" in kwargs:  # It is a Molecule
            return (Molecule.from_atom_dict(**kwargs), stoi, Operation.ADD)
        raise ValueError(f"Unknown species: {kwargs}")

    def enzyme(self, items: tp.List) -> BasicSpecies:
        return BasicSpecies(items[0])

    def ion(self, items: tp.List) -> tp.Dict:
        return items[0]

    def anion(self, items: tp.List) -> tp.Dict:
        if isinstance(items[0], lark.Token) and items[0].type == "ELECTRON":
            return {"name": "e_aq", "charge": -1, "radical": False, "stoi": 1}
        items[0]["charge"] -= len(items) - 1
        return items[0]

    def cation(self, items: tp.List) -> tp.Dict:
        items[0]["charge"] += len(items) - 1
        return items[0]

    def radical(self, items: tp.List) -> tp.Dict:
        items[0].update({"radical": True})
        return items[0]

    def atom(self, items: tp.List) -> tp.Dict:
        return {"name": "".join(items), "radical": False, "charge": 0}

    def sub_mol(self, items: tp.List) -> tp.Tuple[tp.Dict, int]:
        (stoi,) = items[1:2] or [1]
        return (items[0], stoi)

    def molecule(self, items: tp.List) -> tp.Dict:
        return {"atoms": items, "radical": False, "charge": 0}

    def multiple(self, items: tp.List) -> int:
        return int(items[0])

    def stoi(self, items: tp.List) -> int:
        return int(items[0])

    def SIGN(self, sign: lark.lexer.Token) -> Operation:
        if sign == "+":
            return Operation.ADD
        return Operation.SUBTRACT

    def CAPITAL_LETTER(self, items: tp.List) -> str:
        return str(items[0])

    def LOWER_LETTER(self, items: tp.List) -> str:
        return str(items[0])

    def element(self, items: tp.List) -> str:
        return str(items[0])


class ReactionParser(object):
    """
    Actual Parser using the Lark grammar/Transformer.
    """

    GRAMMAR: Path = Path(__file__).parent / r"grammar//reaction_grammar.lark"

    parser: lark.Lark
    atomix: Atomix

    def __init__(self) -> None:
        with open(self.GRAMMAR, "r") as f:
            self.parser = lark.Lark(f.read())
        self.atomix = Atomix()

    def _parse(self, value: str) -> GenericReaction | Atom | Molecule:
        """ """
        # Parse and then transform the output.
        tree = self.parser.parse(value.strip())
        return TreeToReaction().transform(tree)

    def parse_reaction(self, reaction: str) -> GenericReaction:
        """
        Parse Reaction String.
        - K Reactions syntax: aA + bB -> cC + dD
        - Enzymatic reactions syntax: aA --ENZYME>> cC + dD

        Rem: Capital letters stand for Species, lower letters for
        stoichiometry.

        Raises:
            TypeError
                If cannot parse to a dict.

        Returns:
            GenericReaction
                TypedDict structure holding the reaction metadata.
        """
        r = self._parse(reaction)
        if not isinstance(r, dict):
            raise TypeError(f"Unable to parse: '{reaction}' as a chemical reaction")
        for _, (species, _, _) in r["reactants"].items():
            if isinstance(species, (Atom, Molecule)):
                try:
                    self._set_mass(species)
                except KeyError:
                    logger.warning(
                        f"Unable to compute mass of {species} "
                        "Setting it to default mass, i.e. 0"
                    )
                    species.mass = 0
        for _, (species, _, _) in r["products"].items():
            if isinstance(species, (Atom, Molecule)):
                try:
                    self._set_mass(species)
                except KeyError:
                    logger.warning(
                        f"Unable to compute mass of {species} "
                        "Setting it to default mass, i.e. 0"
                    )
                    species.mass = 0
        if isinstance(r["enzyme"], (Atom, Molecule)):
            try:
                self._set_mass(r["enzyme"])
            except KeyError:
                logger.warning(
                    f"Unable to compute mass of {species} "
                    "Setting it to default mass, i.e. 0"
                )
        return r

    def parse_species(self, species: str) -> BasicSpecies | Atom | Molecule:
        """
        Method to Parse a single species Name.

        Returns:
            BasicSpecies | Atom | Molecule
        """
        _ = self._parse(species)
        if not (isinstance(_, tuple) and len(_) == 3):
            raise TypeError(f"Unable to parse: '{species}' as a chemical species")
        sp = _[0]
        if isinstance(sp, (Atom, Molecule)):
            self._set_mass(sp)
            return sp
        if isinstance(sp, BasicSpecies):
            return sp
        raise TypeError(f"Unable to parse: '{species}' as a chemical species")

    def _set_mass(self, species: Atom | Molecule) -> None:
        """
        Try to get the mass of all Atoms.

        Raises:
            KeyError
                if unknown species
                i.e. not in 'data/_data_atomic_weights.txt'
        """
        # Not usable yet. Issues with some molecule's name and the
        # "parse_species" method used in several places.
        return
        if isinstance(species, Atom):
            species.mass = self.atomix.get(species.raw_label())
            return

        # species is a molecule
        for atom, _ in species:
            atom.mass = self.atomix.get(atom.raw_label())
        return
