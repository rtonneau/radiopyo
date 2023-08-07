#!/usr/bin/env python3
"""
Small parser to be able to parse reaction def from strings like:
xA + yB -> iC + kD
xA + yB --enzyme>> iC + kD
Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 29/06/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import logging
import typing as tp
from functools import reduce
from pathlib import Path

import lark
import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.utils.atomic_weights import Atomix

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from collections.abc import Iterator

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = []

# -------------------------------------- Logging ------------------------------------- #
logger = logging.getLogger("radiopyo")

# ------------------------------------------------------------------------------------ #
#                                DECORATOR DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class GenericReaction(tp.TypedDict):
    products: tp.Dict[str, tp.Tuple[Atom | Molecule, int]]
    reactants: tp.Dict[str, tp.Tuple[Atom | Molecule, int]]
    constants: tp.Dict[str, float]
    enzyme: BasicSpecies | Atom | Molecule | None


class KwargsAtom(tp.TypedDict):
    name: str
    charge: int
    radical: bool
    mass: float


class BasicSpecies(object):
    """ The most simple building block. Basically just a named thing."""
    name: str

    def __init__(self,
                 name: str,
                 ) -> None:
        self.name = name

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.__str__()

    def label(self) -> str:
        return self.__str__()


class ChemicalSpecies(BasicSpecies):
    """
    More complex entity. Contains charge, radical and mass.
    """
    charge: int
    radical: bool
    _mass: tp.Optional[float]

    def __init__(self,
                 name: str,
                 charge: int = 0,
                 radical: bool = False,
                 mass: tp.Optional[float] = None
                 ) -> None:
        super().__init__(name)
        self.charge = charge
        self.radical = radical
        self._mass = mass

    def __str__(self) -> str:
        out = [super().__str__(),]
        if self.radical:
            out.append("_r")
        if self.charge > 0:
            out += ["_plus"]*self.charge
        elif self.charge < 0:
            out += ["_minus"]*np.abs(self.charge)
        return "".join(out)

    @property
    def mass(self) -> float:
        if self._mass is None:
            raise ValueError(f"No mass defined for {self.__str__()}")
        return self._mass

    @mass.setter
    def mass(self, value: float) -> None:
        self._mass = value

    def _latex_radical_charge(self) -> str:
        """ For LateX rendering of radical and charge only.
        -> Shared by Atoms and Molecules.
        """
        out = ["^{"]
        if self.radical:
            out.append(r"\bullet")
        charge = f"{self.charge}" if abs(self.charge) > 1 else ""
        if self.charge > 0:
            out += [f"{charge}+"]
        elif self.charge < 0:
            out += [f"{charge}-"]
        out.append("}")
        return "".join(out) if len(out) > 2 else ""

    def to_latex(self) -> str:
        """Produces string formatted for LateX rendering (e.g. by MatplotLib)"""
        out = ["$", self.name, self._latex_radical_charge(), "$"]
        return "".join(out)


class HydratedElectron(ChemicalSpecies):
    """
    Hydrated electrons have their own class definition. No mass for them.
    """

    def __init__(self) -> None:
        super().__init__(name="e_aq",
                         charge=-1,
                         radical=False,
                         mass=0)

    def __str__(self) -> str:
        return self.name

    def to_latex(self) -> str:
        return "$e_{aq}$"


class Atom (ChemicalSpecies):
    """ Atom definition """

    def raw_label(self) -> str:
        """Just a getter for the name (does not include charge or radical)"""
        return self.name


class Molecule(ChemicalSpecies):
    """ Molecule definition -> Collection of Atoms.
     To create a Molecule out of atoms, they cannot be ions nor radical. It is the
     molecule which can be charged or radical"""
    atoms: tp.List[tp.Tuple[Atom, int]]

    def __init__(self,
                 atoms: tp.List[tp.Tuple[Atom, int]],
                 charge: int = 0,
                 radical: bool = False
                 ):
        super().__init__("", charge, radical)
        if len(atoms) == 2 and isinstance(atoms[0], Atom):
            self.atoms = [atoms,]
        else:
            self.atoms = atoms

        # Let's do some checks
        for atom, _ in self.atoms:
            if atom.charge != 0:
                raise ValueError(f"{atom} must be neutral (and not {atom.charge}) "
                                 "to compose a molecule")
            if atom.radical is not False:
                raise ValueError(f"{atom} must not be a radical "
                                 "to compose a molecule")

    def __iter__(self) -> Iterator[tp.Tuple[Atom, int]]:
        """Let's easily iterate of (atoms, stoichiometry)"""
        return iter(self.atoms)

    def __str__(self) -> str:
        out = []
        for atom, stoi in self.atoms:
            out.append(f"{atom.label()}{stoi if stoi>1 else ''}")
        if self.radical:
            out.append("_r")
        if self.charge > 0:
            out += ["_plus"]*self.charge
        elif self.charge < 0:
            out += ["_minus"]*np.abs(self.charge)
        return "".join(out)

    def to_latex(self) -> str:
        out = ["$"]
        for (atom, stoi) in self.atoms:
            out.append(f"{atom.to_latex().strip('$')}{f'_{stoi}' if stoi>1 else ''}")
        out.append(self._latex_radical_charge())
        out.append("$")
        return "".join(out)

    @classmethod
    def from_atom_dict(cls,
                       atoms: tp.List[tp.Tuple[KwargsAtom, int]],
                       charge: int = 0,
                       radical: bool = False,
                       ) -> Molecule:
        return cls([(Atom(**elt), stoi) for elt, stoi in atoms],
                   charge,
                   radical)

    @property
    def mass(self) -> float:
        if self._mass is not None:
            return self._mass
        return reduce(lambda tot, atom: tot+atom[1]*atom[0].mass, self.atoms, 0.0)

    @mass.setter
    def mass(self, value: float) -> None:
        self._mass = value


class TreeToReaction(lark.Transformer):
    """Lark Class handling the conversion.
    See: https://lark-parser.readthedocs.io/en/latest/visitors.html
    For some doc and examples.
    """

    def start(self, items: tp.List) -> GenericReaction | Atom | Molecule:
        return items[0]

    def ereaction(self, items: tp.List) -> GenericReaction:
        return GenericReaction(reactants={items[0][0].label(): items[0]},
                               products=items[2],
                               enzyme=items[1],
                               constants={},
                               )

    def kreaction(self, items: tp.List) -> GenericReaction:
        return GenericReaction(reactants=items[0],
                               products=items[1],
                               constants={},
                               enzyme=None)

    def reactants(self, items: tp.List) -> tp.List:
        return items[0]

    def products(self, items: tp.List) -> tp.List:
        return items[0]

    def species_list(self, items: tp.List[tp.Tuple[ChemicalSpecies, int]]) -> tp.Dict:
        """"""
        out: tp.Dict[str, tp.Tuple[ChemicalSpecies, int]] = {}
        for species, stoi in items:
            label = species.label()
            out[label] = (species, out.get(label, (0, 0))[1]+stoi)
        return out

    def species(self, items: tp.List) -> tp.Tuple[BasicSpecies, int]:

        if len(items) > 1:
            multiple = items[0]
            kwargs = items[1]
        else:
            multiple = 1
            kwargs = items[0]

        if isinstance(kwargs, str):
            return (BasicSpecies(kwargs), multiple)
        if "name" in kwargs:
            if kwargs["name"] == "e_aq":
                return (HydratedElectron(), multiple)
            return (Atom(**kwargs), multiple)
        if "atoms" in kwargs:  # It is a Molecule
            return (Molecule.from_atom_dict(**kwargs), multiple)
        raise ValueError(f"Unknown species: {kwargs}")

    def enzyme(self, items: tp.List) -> str:
        return items[0]

    def ion(self, items: tp.List) -> tp.Dict:
        return items[0]

    def anion(self, items: tp.List) -> tp.Dict:
        if isinstance(items[0], lark.Token) and items[0].type == "ELECTRON":
            return {"name": "e_aq", "charge": -1, "radical": False, "stoi": 1}
        items[0]["charge"] -= (len(items)-1)
        return items[0]

    def cation(self, items: tp.List) -> tp.Dict:
        items[0]["charge"] += (len(items)-1)
        return items[0]

    def radical(self, items: tp.List) -> tp.Dict:
        items[0].update({"radical": True})
        return items[0]

    def atom(self, items: tp.List) -> tp.Dict:
        return {"name": "".join(items), "radical": False, "charge": 0}

    def sub_mol(self, items: tp.List) -> tp.Tuple[tp.Dict, int]:
        stoi, = items[1:2] or [1]
        return (items[0], stoi)

    def molecule(self, items: tp.List) -> tp.Dict:
        return {"atoms": items, "radical": False, "charge": 0}

    def multiple(self, items: tp.List) -> int:
        return int(items[0])

    def stoi(self, items: tp.List) -> int:
        return int(items[0])

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
    GRAMMAR: Path = Path(__file__).parent/r"grammar//reaction_grammar.lark"

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

        Rem: Capital letters stand for Species, lower letters for stoichiometry.

        Raises:
            TypeError
                If cannot parse to a dict.

        Returns:
            GenericReaction
                TypedDict structure holding the reaction metadata.
        """
        r = self._parse(reaction)
        if not isinstance(r, dict):
            raise TypeError(
                f"Unable to parse: '{reaction}' as a chemical reaction")
        for _, (species, _) in r["reactants"].items():
            if isinstance(species, (Atom, Molecule)):
                try:
                    self._set_mass(species)
                except KeyError:
                    logger.warning(f"Unable to compute mass of {species}. "
                                   "Setting it to default mass, i.e. 0")
                    species.mass = 0
        for _, (species, _) in r["products"].items():
            if isinstance(species, (Atom, Molecule)):
                try:
                    self._set_mass(species)
                except KeyError:
                    logger.warning(f"Unable to compute mass of {species}. "
                                   "Setting it to default mass, i.e. 0")
                    species.mass = 0
        if isinstance(r["enzyme"], (Atom, Molecule)):
            try:
                self._set_mass(r["enzyme"])
            except KeyError:
                logger.warning(f"Unable to compute mass of {species}. "
                               "Setting it to default mass, i.e. 0")
                species.mass = 0
        return r

    def parse_species(self, species: str) -> BasicSpecies | Atom | Molecule:
        """
        Method to Parse a single species Name.

        Returns:
            BasicSpecies | Atom | Molecule        
        """
        _ = self._parse(species)
        if isinstance(_, (Atom, Molecule)):
            self._set_mass(_)
            return _
        if isinstance(_, BasicSpecies):
            return _
        raise TypeError(f"Unable to parse: '{species}' as a chemical species")

    def _set_mass(self, species: Atom | Molecule) -> None:
        """ 
        Try to get the mass of all Atoms.

        Raises:
            KeyError
                if unknown species (i.e. not in 'data/_data_atomic_weights.txt')
        """
        if isinstance(species, Atom):
            try:
                species.mass = self.atomix.get(species.raw_label())
                return
            except KeyError as e:
                # Let's add an exception for Carbon Centered molecules named R, ROOH,...
                if "R" not in species.raw_label():
                    msg = f"Unable to compute mass of: {species}"
                    raise KeyError(msg) from e
                species.mass = 0
                return

        for atom, _ in species:
            try:
                atom.mass = self.atomix.get(atom.raw_label())
            except KeyError as e:
                # Let's add an exception for Carbon Centered molecules named R, ROOH,...
                if "R" not in atom.raw_label():
                    msg = f"Unable to compute mass of: {atom} for molecule: {species}"
                    raise KeyError(msg) from e
                atom.mass = np.inf
                return
        return
