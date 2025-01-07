#!/usr/bin/env python3
"""
The species submodule contains all class definitions related to species
definitions.
It also provides a collection to ease their handling.

Classes:
    ChemSpecies: Class used for species with cc (cst or dyn).
    RawSpecies: Basic class for species.
    H2O: A predefined species for water.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""
# ----------------------------- External imports -----------------------------#
from __future__ import annotations

import typing as tp
from functools import reduce

import numpy as np

# ---------------------------- Internal imports ----------------------------- #
from radiopyo.chemistry.utils.factories import polynomial_factory
from radiopyo.chemistry.utils.units import Temperature

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from collections.abc import Iterator

# --------------------------------------------------------------------------- #
__all__ = [
    "H2O",
    "Atom",
    "BasicSpecies",
    "ChemicalSpecies",
    "HydratedElectron",
    "Molecule",
    "SimSpecies",
]


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #


class GenericReaction(tp.TypedDict):
    products: tp.Dict[str, tp.Tuple[ChemicalSpecies, int]]
    reactants: tp.Dict[str, tp.Tuple[ChemicalSpecies, int]]
    constants: tp.Dict[str, float | tp.List[float]]
    enzyme: ChemicalSpecies | None


class KwargsAtom(tp.TypedDict):
    name: str
    charge: int
    radical: bool
    mass: float


class BasicSpecies(object):
    """
    The most simple building block. Basically just a named 'thing'.
    """

    name: str

    def __init__(
        self,
        name: str,
    ) -> None:
        self.name = name

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.__str__()

    @property
    def label(self) -> str:
        return self.__str__()

    def __eq__(self, other: object) -> bool:
        """
        Checks if the given object is equal to the current species.

        Args:
            other (object): The object to be compared.

        Returns:
            bool: True if the objects are equal, False otherwise.
        """
        if isinstance(other, str):
            return self.label == other
        if isinstance(other, BasicSpecies):
            return self.label == other.label
        raise NotImplementedError(f"Equality between RawSpecies and {other}")


class ChemicalSpecies(BasicSpecies):
    """
    More complex entity. Contains charge, radical and mass.
    """

    charge: int
    radical: bool
    _mass: tp.Optional[float]

    def __init__(
        self,
        name: str,
        charge: int = 0,
        radical: bool = False,
        mass: tp.Optional[float] = None,
    ) -> None:
        super().__init__(name)
        self.charge = charge
        self.radical = radical
        self._mass = mass

    def __str__(self) -> str:
        out = [
            super().__str__(),
        ]
        if self.radical:
            out.append("_r")
        if self.charge > 0:
            out += ["_plus"] * self.charge
        elif self.charge < 0:
            out += ["_minus"] * np.abs(self.charge)
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
        """For LateX rendering of radical and charge only.
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
        """
        Produces string formatted for LateX rendering
        (e.g. by MatplotLib)
        """
        out = ["$", self.name, self._latex_radical_charge(), "$"]
        return "".join(out)

    def initial_cc(self, T: tp.Optional[float | Temperature] = None) -> float:
        """
        Returns the initial concentration of the species at a given
            temperature.

        Args:
            T (float | Temperature): Temperature in Celsius. If None,
                default temperature is used.

        Returns:
            float: The initial concentration of the species at the given
                temperature.
        """
        T = Temperature(T)
        return self.cc_init(T.celsius) if self.cc_init(T.celsius) > 0 else 0

    def set_initial_cc(
        self, value: float | tp.List[float] | np.polynomial.Polynomial
    ) -> None:
        """
        Sets the initial concentration coefficient of the species.

        Args:
            value (float | list[float]): The initial concentration
                coefficient(s) of the species.

        Returns:
            None
        """
        self.cc_init = polynomial_factory(value, log10=False)


class HydratedElectron(ChemicalSpecies):
    """
    Hydrated electrons have their own class definition.
    No mass for them.
    """

    def __init__(
        self,
    ) -> None:
        super().__init__(
            name="e_aq",
            charge=-1,
            radical=False,
            mass=0,
        )

    def __str__(self) -> str:
        return self.name

    def to_latex(self) -> str:
        return "$e_{aq}$"


class Atom(ChemicalSpecies):
    """Atom definition"""

    def raw_label(self) -> str:
        """
        Just a getter for the name (does not include charge or radical)
        """
        return self.name


class Molecule(ChemicalSpecies):
    """
    Molecule definition -> Collection of Atoms.
    To create a Molecule out of atoms, they cannot be ions nor radical.
    It is the molecule which can be charged or radical.
    """

    atoms: tp.List[tp.Tuple[Atom, int]]

    def __init__(
        self,
        atoms: tp.List[tp.Tuple[Atom, int]],
        charge: int = 0,
        radical: bool = False,
    ):
        super().__init__("", charge, radical)
        if len(atoms) == 2 and isinstance(atoms[0], Atom):
            self.atoms = [
                atoms,
            ]
        else:
            self.atoms = atoms

        # Let's do some checks
        for atom, _ in self.atoms:
            if atom.charge != 0:
                raise ValueError(
                    f"{atom} must be neutral (and not {atom.charge}) "
                    "to compose a molecule"
                )
            if atom.radical is not False:
                raise ValueError(f"{atom} must not be a radical to compose a molecule")

    def __iter__(self) -> Iterator[tp.Tuple[Atom, int]]:
        """Let's easily iterate of (atoms, stoichiometry)"""
        return iter(self.atoms)

    def __str__(self) -> str:
        out = []
        for atom, stoi in self.atoms:
            out.append(f"{atom.label}{stoi if stoi>1 else ''}")
        if self.radical:
            out.append("_r")
        if self.charge > 0:
            out += ["_plus"] * self.charge
        elif self.charge < 0:
            out += ["_minus"] * np.abs(self.charge)
        return "".join(out)

    def to_latex(self) -> str:
        out = ["$"]
        for atom, stoi in self.atoms:
            out.append(f"{atom.to_latex().strip('$')}{f'_{stoi}' if stoi>1 else ''}")
        out.append(self._latex_radical_charge())
        out.append("$")
        return "".join(out)

    @classmethod
    def from_atom_dict(
        cls,
        atoms: tp.List[tp.Tuple[KwargsAtom, int]],
        charge: int = 0,
        radical: bool = False,
    ) -> Molecule:
        return cls([(Atom(**elt), stoi) for elt, stoi in atoms], charge, radical)

    @property
    def mass(self) -> float:
        if self._mass is not None:
            return self._mass
        return reduce(lambda tot, atom: tot + atom[1] * atom[0].mass, self.atoms, 0.0)

    @mass.setter
    def mass(self, value: float) -> None:
        self._mass = value


class SimSpecies:
    """
    Class used to define species from concentration coefficient (cc).

    This class extends the RawSpecies class and adds attributes and
    methods related to the concentration coefficient.

    Attributes:
        label (str): The label of the species.
            cc_init (Polynomial): The initial concentration coefficient
            of the species.

    Methods:
        set_initial_cc(value): Sets the initial concentration
            coefficient of the species.
        ... (other methods specific to ChemSpecies)
    """

    species: ChemicalSpecies
    cc_init: np.polynomial.Polynomial  # Store initial value as [mol/l]
    index: int

    def __init__(
        self,
        raw_species: ChemicalSpecies,
        index: int = -1,
        cc_init: tp.Optional[tp.List[float] | float | np.polynomial.Polynomial] = None,
    ) -> None:
        """ """
        self.raw_species = raw_species
        self.index = index
        cc_init = cc_init if cc_init is not None else -1
        self.cc_init = polynomial_factory(cc_init, log10=False)

    def __repr__(self) -> str:
        """
        Returns a string representation of the object.

        Returns:
            str: The representation of the object.
        """
        return f"SimSpecies('{self.species}', {self.index}, {self.cc_init})"

    def __copy__(self) -> SimSpecies:
        """
        Creates a deep copy of the current SimSpecies instance.

        Returns:
            SimSpecies: A new SimSpecies instance with the same label,
                index, and initial concentration coefficients.
        """
        return SimSpecies(self.species, self.index, self.cc_init.coef)

    @property
    def label(self) -> str:
        return self.species.label

    def initial_cc(self, T: tp.Optional[float | Temperature] = None) -> float:
        """
        Returns the initial concentration of the species at a given
            temperature.

        Args:
            T (float | Temperature): Temperature in Celsius. If None,
                default temperature is used.

        Returns:
            float: The initial concentration of the species at the given
                temperature.
        """
        T = Temperature(T)
        return self.cc_init(T.celsius) if self.cc_init(T.celsius) > 0 else 0

    def set_initial_cc(
        self, value: float | tp.List[float] | np.polynomial.Polynomial
    ) -> None:
        """
        Sets the initial concentration coefficient of the species.

        Args:
            value (float | list[float]): The initial concentration
                coefficient(s) of the species.

        Returns:
            None
        """
        self.cc_init = polynomial_factory(value, log10=False)


H2O = SimSpecies(
    Molecule([(Atom("H", 0, False, 1.00784), 2), (Atom("O", 0, False, 15.999), 1)]),
    cc_init=[55.50, 6.075e-3, 4.110e-4, 1.496e-6, 2.619e-9],
)


# Define H2O here
def get_H2O_cc(T: tp.Optional[float | Temperature] = None) -> float:
    """
    Returns the initial concentration of H2O at a given temperature.

    Args:
        T (float | Temperature): Temperature in Celsius. If None,
            default temperature is used.

    Returns:
        float: The initial concentration of H2O at the given temperature.
    """
    return H2O.initial_cc(T)
