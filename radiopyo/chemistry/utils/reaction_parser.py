#!/usr/bin/env python3
"""
Small parser to be able to parse reaction like:
xA + yB -> iC + kD
xA + yB --enzyme>> iC + kD
Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 29/06/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp

import lark
import numpy as np

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    pass

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


class KwargsAtom(tp.TypedDict):
    name: str
    multiple: int
    charge: int
    radical: bool
    mass: float


class Atom(object):
    """
    """
    name: str
    multiple: int
    charge: int
    radical: bool
    mass: float

    def __init__(self,
                 name: str,
                 multiple: int = 1,
                 charge: int = 0,
                 radical: bool = False,
                 mass: float = 0
                 ) -> None:
        self.name = name
        self.multiple = multiple
        self.charge = charge
        self.radical = radical
        self.mass = mass

    def __str__(self) -> str:
        out = [f"{self.multiple}",
               self.name] if self.multiple > 1 else [self.name]
        if self.radical:
            out.append("_r")
        if self.charge > 0:
            out += ["_plus"]*self.charge
        elif self.charge < 0:
            out += ["_minus"]*self.charge
        return "".join(out)

    def __repr__(self) -> str:
        return self.__str__()


class ChemicalElement(Atom):
    def __init__(self, name: str, multiple: int = 1) -> None:
        super().__init__(name, multiple=multiple, charge=0, mass=0)

    def __str__(self) -> str:
        if self.multiple == 1:
            return self.name
        return f"{self.multiple}{self.name}"


class HydratedElectron(ChemicalElement):

    def __init__(self, multiple: int = 1) -> None:
        Atom.__init__(self, name="e_aq", multiple=multiple, charge=-1, mass=0)


class Molecule(object):
    """
    """
    atoms: tp.List[tp.Tuple[Atom, int]]
    multiple: int

    def __init__(self,
                 atoms: tp.List[tp.Tuple[Atom, int]],
                 multiple: int = 1,
                 charge: int = 0,
                 radical: bool = False
                 ):
        if len(atoms) == 2 and isinstance(atoms[0], Atom):
            self.atoms = [atoms,]
        else:
            self.atoms = atoms

        self.multiple = multiple
        self.radical = radical
        self.charge = charge

    @classmethod
    def from_atom_dict(cls,
                       atoms: tp.List[tp.Tuple[KwargsAtom, int]],
                       multiple: int = 1,
                       charge: int = 0,
                       radical: bool = False,
                       ) -> Molecule:
        return cls([(Atom(**elt), stoi) for elt, stoi in atoms],
                   multiple,
                   charge,
                   radical)

    def __str__(self) -> str:
        out = []
        if self.multiple > 1:
            out.append(str(self.multiple))
        for elt, stoi in self.atoms:
            out.append(f"{elt}")
            if stoi > 1:
                out[-1] += f"{stoi}"
        if self.radical:
            out.append("_r")
        if self.charge > 0:
            out += ["_plus"]*self.charge
        elif self.charge < 0:
            out += ["_minus"]*np.abs(self.charge)
        return "".join(out)

    def __repr__(self) -> str:
        return self.__str__()


class TreeToReaction(lark.Transformer):
    """ """
    
    species_list = list

    def start(self, items: tp.List) -> tp.Dict:
        return items[0]

    def ereaction(self, items: tp.List) -> tp.Dict:
        return {"substrate": items[0]} | items[1] | items[2]

    def kreaction(self, items: tp.List) -> tp.Dict:
        return items[0] | items[1]

    def reactants(self, items: tp.List) -> tp.Dict:
        return {"reactants": items[0]}

    def products(self, items: tp.List) -> tp.Dict:
        return {"products": items[0]}

    def species(self, items: tp.List) -> Molecule | Atom:
        if len(items) > 1:
            items[1]["multiple"] = items[0]
            kwargs = items[1]
        else:
            kwargs = items[0]

        if isinstance(kwargs, str):
            return ChemicalElement(kwargs)
        if "name" in kwargs:
            if kwargs["name"] == "e_aq":
                return HydratedElectron(kwargs["multiple"])
            return Atom(**kwargs)
        if "atoms" in kwargs:  # It is a Molecule
            return Molecule.from_atom_dict(**kwargs)
        raise ValueError(f"Unknown species: {kwargs}")

    def enzyme(self, items: tp.List) -> tp.Dict:
        return {"enzyme": items[0]}

    def ion(self, items: tp.List) -> tp.Dict:
        return items[0]

    def anion(self, items: tp.List) -> tp.Dict:
        if isinstance(items[0], lark.Token) and items[0].type == "ELECTRON":
            return {"name": "e_aq", "charge": -1, "radical": False, "stoi": 1}
        items[0]["charge"] -= 1
        return items[0]

    def cation(self, items: tp.List) -> tp.Dict:
        items[0]["charge"] += 1
        return items[0]

    def radical(self, items: tp.List) -> tp.Dict:
        items[0].update({"radical": True})
        return items[0]

    def atom(self, items: tp.List) -> tp.Dict:
        return {"name": items[0], "radical": False, "charge": 0}

    def sub_mol(self, items: tp.List) -> tp.Tuple[tp.Dict, int]:
        stoi, = items[1:2] or [1]
        return (items[0], stoi)

    def molecule(self, items: tp.List) -> tp.Dict:
        return {"atoms": items, "radical": False, "charge": 0, "multiple": 1}

    def multiple(self, items: tp.List) -> int:
        return int(items[0])

    def stoi(self, items: tp.List) -> int:
        return int(items[0])

    def CAPITAL_LETTER(self, items: tp.List) -> str:
        return str(items[0])

    def LOWER_LETTER(self, items: tp.List) -> str:
        return str(items[0])

    def element(self, items: tp.List) -> str:
        return items[0]

    def WORD(self, items: tp.List) -> str:
        return str(items)
