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

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.basic_reaction import BasicReaction
from radiopyo.chemistry.species_collection import ListSpecies

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.basic_reaction import GenericReaction
    from radiopyo.chemistry.species import SimSpecies
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
    enzyme: SimSpecies

    def __init__(
        self,
        k_value: float,
        k_micha: float,
    ):
        super().__init__()
        self.k_value = k_value
        self.k_micha = k_micha

    @property
    def substrate(self) -> SimSpecies:
        return next(iter(self.reactants)).species

    def raw_k_value(self) -> float:
        """
        Produce the k_value as provided by config files.

        Returns:
            float: k_value as provided by config files.
        """
        # Case of n*A -> B, then n*k is measured/given
        if len(self.reactants) == 1:
            sp = next(iter(self.reactants))
            return self.k_value * sp.stoi
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
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> tp.Self:
        """
        Create a MichaelisMenten from a 'GenericReaction', Typed-Dict
        structure produced by radiopyo module parser.

        Args:
            generic (GenericReaction): Typed-Dict structure holding the
                reaction metadata.
            list_sp (ListSpecies, optional): Collection of species
                already involved in the simulation. Defaults to None.

        Returns:
            MichaelisMenten: An instance of MichaelisMenten.
        """
        if generic["enzyme"] is None:
            raise ValueError(f"No enzyme found for enzymatic reaction: {generic}")
        reaction = cls(
            generic["constants"]["k_value"],  # type:ignore[arg-type]
            generic["constants"]["k_micha"],  # type:ignore[arg-type]
        )
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

        reaction.enzyme = list_sp.get_or_create_dyn_species(generic["enzyme"])

        for _, (sp, stoi, ope) in generic["reactants"].items():
            if reaction.tot_number_reactants() < 1:
                substrate = list_sp.get_or_create_dyn_species(sp)
                reaction.add_as_reactant(substrate, stoi, ope)
                continue
            raise ValueError(
                f"MichaelisMenten reaction must have only one substrate. Found: {generic}"
            )

        for _, (sp, stoi, ope) in generic["products"].items():
            sim_sp = list_sp.get_or_create_dyn_species(sp)
            reaction.add_as_product(sim_sp, stoi, ope)
        return reaction

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
