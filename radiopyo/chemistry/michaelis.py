#!/usr/bin/env python3
"""
The michaelis submodule:
The michaelis submodule contains all class definitions related to
enzymatic reactions treatment via the Michaelis-Menten kinetic.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ------------------------------ External imports --------------------------- #
from __future__ import annotations

import hashlib
import typing as tp

from radiopyo.chemistry.basic_chemistry import BasicReaction
from radiopyo.chemistry.collections import ListSpecies

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser import RPARSER

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.species import ChemSpecies
    from radiopyo.parser.reaction_parser import GenericReaction
# --------------------------------------------------------------------------- #

__all__ = [
    "MichaelisMenten",
]

# --------------------------------------------------------------------------- #
#                                CLASS DEFINITION                             #
# --------------------------------------------------------------------------- #


class MichaelisMenten(BasicReaction):
    """
    Class implementing the Michaelis-Menten kinetic. Child class of
    BaseChemicalReaction with few additions. In the Michaelis-Menten
    kinetic the reaction rate is defined by:
        k_value * [enzyme] * [substrate] / (k_micha + [substrate])

    This kinetic is really convenient as it ensures a positive
    denominator during calculation (k_micha > 0).

    Attributes:
        k_value: float
            reaction constant
        k_micha: float
            Michaelis-Menten constant
        substrate: SimSpecies
            Species used as substrate
        enzyme: SimSpecies
            Species used as enzyme
    """

    k_value: float
    k_micha: float
    substrate: ChemSpecies
    enzyme: ChemSpecies

    def __init__(
        self,
        k_value: float,
        k_micha: float,
    ):
        super().__init__()
        self.k_value = k_value
        self.k_micha = k_micha

    def raw_k_value(self) -> float:
        """
        Produce the k_value as provided by config files.

        Returns:
            float: k_value as provided by config files.
        """
        # Case of n*A -> B, then n*k is measured/given
        if len(self.reactants) == 1:
            sp = next(iter(self.reactants.keys()))
            return self.k_value * self.stoi_reactants[sp]
        return self.k_value

    def make_hash(self, method: str = "sha256") -> str:
        """
        Produce a unique identifier (content based) for easy comparison.

        Args:
            method (str): Hash method to use (default: "sha256").

        Returns:
            str: A hash value representing the object.
        """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(self).encode())
        return hasher.hexdigest()

    @classmethod
    def from_generic(
        cls,
        generic: GenericReaction,
        k_value: float,
        k_micha: float,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Create a MichaelisMenten from a 'GenericReaction', Typed-Dict
        structure produced by radiopyo module parser.

        Args:
            generic (GenericReaction): Typed-Dict structure holding the
                reaction metadata.
            k_value (float): k value of the reaction.
            k_micha (float): Michaelis-Menten constant of the reaction.
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.

        Returns:
            MichaelisMenten: An instance of MichaelisMenten.
        """
        if generic["enzyme"] is None:
            raise ValueError(f"No enzyme found for enzymatic reaction: {generic}")

        if isinstance(generic["enzyme"], str):
            enzyme: str = generic["enzyme"]
        else:
            enzyme = generic["enzyme"].label()

        substrate: tp.List[str] = []
        products: tp.List[str] = []

        for sp, (_, stoi) in generic["reactants"].items():
            substrate += [sp] * stoi
        for sp, (_, stoi) in generic["products"].items():
            products += [sp] * stoi
        return cls.from_dict(
            enzyme=enzyme,
            substrate=substrate,
            products=products,
            k_value=k_value,
            k_micha=k_micha,
            list_sp=list_sp,
        )

    @classmethod
    def from_kwargs(
        cls,
        list_sp: tp.Optional[ListSpecies] = None,
        **kwargs: str | tp.List[str] | float | tp.List[float],
    ) -> tp.Self:
        """
        Creates a MichaelisMenten reaction from keyword arguments.

        Args:
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.
            **kwargs: Keyword arguments to create the instance.

        Returns:
            MichaelisMenten: An instance of MichaelisMenten.

        Raises:
            ValueError: If the keyword arguments do not match the
                expected syntax.
        """
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        match kwargs:
            case {
                "enzyme": _,
                "substrate": _,
                "products": _,
                "k_value": _,
                "k_micha": _,
            } if len(kwargs) == 5:
                return cls.from_dict(
                    list_sp=list_sp,
                    **kwargs,  # type:ignore[arg-type]
                )
            case {"reaction": r, "k_value": k_value, "k_micha": k_micha} if len(
                kwargs
            ) == 3:
                gen = RPARSER.parse_reaction(r)  # type:ignore[arg-type]
                return cls.from_generic(
                    generic=gen,
                    k_value=k_value,  # type:ignore[arg-type]
                    k_micha=k_micha,  # type:ignore[arg-type]
                    list_sp=list_sp,
                )
        raise ValueError(f"Wrong KReaction syntax {kwargs}")

    @classmethod
    def from_dict(
        cls,
        enzyme: str,
        substrate: tp.List[str] | str,
        products: tp.List[str],
        k_value: float,
        k_micha: float,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Creates a MichaelisMenten instance from a dictionary.

        Args:
            enzyme (str): The enzyme name.
            substrate (List[str] | str): The substrate name(s).
            products (List[str]): The product name(s).
            k_value (float): The k value of the reaction.
            k_micha (float): The k_micha value of the reaction.
            list_sp (ListSpecies, optional): The list of species.
                Defaults to None.

        Returns:
            MichaelisMenten: A MichaelisMenten instance.
        """
        new = cls(k_value, k_micha)
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        new.enzyme = list_sp.get_or_create_cst_species(enzyme.strip())
        if isinstance(substrate, str):
            substrate = [substrate.strip()]
        elif len(substrate) > 1 and len({elt.strip() for elt in substrate}) > 1:
            raise TypeError(
                f"Substrate in Michaelis Menten reaction"
                f" must all be the same if provide as list. Found: {substrate}"
            )

        new.substrate = list_sp.get_or_create_dyn_species(substrate[0].strip())

        for sp in substrate:
            new.add_as_reactant(list_sp.get_or_create_dyn_species(sp.strip()))

        for sp in products:
            new.add_as_product(list_sp.get_or_create_dyn_species(sp.strip()))
        return new

    def as_label(self) -> str:
        """
        Returns a string label for the reaction based on its index.

        The label is in the format "EZ_{index}" if the index is valid.

        Raises:
            KeyError
                If the index is None or invalid.

        Returns:
            str
                A string label for the reaction.
        """
        if self.index is not None and self.index >= 0:
            return f"EZ_{self.index}"

        raise KeyError("Unknown reaction index. Cannot make a label.")

    def __str__(self) -> str:
        out = []
        out.append(self.reactants_as_string())
        out.append("--")
        out.append(f"{self.enzyme.label}")
        out.append(">>")
        out.append(self.products_as_string())
        return " ".join(out)

    def __repr__(self) -> str:
        """
        Returns a string representation of the MichaelisMenten instance.
        """
        return (
            f"MichaelisMenten('{self.__str__()}', "
            f"k_value={round(self.k_value, 5)}, "
            f"k_micha={round(self.k_micha, 5)}"
        )


class MichaelisMentenDict(tp.TypedDict):
    enzyme: str
    substrate: tp.List[str] | str
    products: tp.List[str]
    k_value: float
    k_micha: float
