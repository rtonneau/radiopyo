#!/usr/bin/env python3
"""
The radiolytic submodule contains all class and function definitions to
include reactions involving High Energy Particles with the background
species => Production of radicals and charged species. These reactions
are defined using so called "G-values" -> radiolytic yields i.e.
[#radical / 100eV / incident particle].

/!\ Here, the solvent density is assumed to be the same as Water,
thus a constant concentration of 1 kg/l.

Classes:
    - RadiolyticReaction: defines a simple reaction involving a High
        Energy Particle

Functions:
    - ge_to_kr: converts a G-value to K-rate constant
    - kr_to_ge: converts a K-rate constant to a G-value

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import hashlib
import typing as tp
from typing import List

from scipy.constants import Avogadro, elementary_charge  # type: ignore[import]

from radiopyo.chemistry.species import SimSpecies

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.species_collection import ListSpecies

from .basic_reaction import BasicReaction

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    pass
# --------------------------------------------------------------------------- #

__all__ = [
    "ge_to_kr",
    "RadiolyticReaction",
    "umol_per_joule_to_ge",
    "RadiolyticReactionDict",
]

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
# Ge units => radical / 100eV / incident particle
# Kr units => mol/l/Gy

# Define some functions for unit conversion:


def umol_per_joule_to_ge(value: float) -> float:
    """
    Converts a given value from umol/J to Ge units.

    Args:
        value (float): The value to be converted.

    Returns:
        float: The converted value in Ge units.

    Note:
        Unused and probably wrong function.
    """
    return value * 1e-6 * Avogadro / elementary_charge * 100


def ge_to_kr(ge: float) -> float:
    """
    Convert a G-value to a K-rate constant.

    Args:
        ge (float): The G-value to convert.

    Returns:
        float: The K-rate constant.
    """
    d: float = 1.0  # Solvent density [kg/l]
    return ge * d / elementary_charge / 100.0 / Avogadro


def kr_to_ge(kr: float) -> float:
    """
    Converts a given K-rate constant to a G-value.

    Args:
        kr (float): The K-rate constant to be converted.

    Returns:
        float: The converted G-value.
    """
    d: float = 1.0  # Solvent density [kg/l]
    return kr / d * elementary_charge * 100.0 * Avogadro


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #


class RadiolyticReaction(BasicReaction):
    """
    This class represents a radiolytic reaction. It is a specific type
    of chemical reaction where the production rate only depends on the
    nature and energy of the incident radiation. The reaction constant
    (G value) uniquely defines the rate of the reaction.

    The reaction consists of a list of products and a reaction constant
    (G value). A radiolytic reaction is constructed from the
    BaseChemicalReaction. Here no reactant should be involved.
    Moreover, the production rate only depends on the nature and energy
    of incident radiation ==> LET and whether or not beam is ON.
    ==> Only the list of products and the reaction cst (G value) are
    then relevant.

    Attributes:
        _reaction_cst (float): Reaction constant (G value)[mol/l/Gy].
    """

    _reaction_cst: float
    _OPE_STRING: str = "->"

    def __init__(self, reaction_cst: float):
        super().__init__()
        self._reaction_cst = reaction_cst  # [mol/l/Gy]

    def make_hash(self, method: str = "sha256") -> str:
        """
        Generate a hash for the reaction. The G value is not included
        in the hash because two radiolytic reactions of the same species
        can have different G values. The hash is intended to be used to
        identify the reaction uniquely, without the need to store the
        G value. When the reaction is to be used, the actual G value
        should be provided. The hash is computed from the reaction's
        index and the list of products.

        Args:
            method (str): The hashing algorithm to use. Defaults to
                "sha256".

        Returns:
            str: The hash of the instance as a hexadecimal string.

        Raises:
            ValueError: If the hashing algorithm is not supported.
        """
        hasher = hashlib.new(method, usedforsecurity=False)
        hasher.update(repr(self).encode())
        return hasher.hexdigest()

    def as_label(self) -> str:
        """
        Returns a string label for the reaction based on its index.

        Returns:
            str
                A string label in the format "G_{index}" if index is
                valid.

        Raises:
            KeyError
                If the index is None or invalid.
        """
        return f"G_{self.species().label}"

    @classmethod
    def from_kwargs(
        cls,
        list_sp: ListSpecies | None = None,
        **kwargs: str | List[str] | float | List[float],
    ) -> tp.Self:
        raise NotImplementedError

    def species(self) -> SimSpecies:
        """
        Returns the ChemSpecies instance associated with this reaction.

        Returns:
            ChemSpecies
                The ChemSpecies instance.
        """
        return self.products[self.label()].species

    def label(self) -> str:
        """
        Method to generate a string representation of the label of the
        reaction.

        Returns:
            str
                A string label of the reaction.
        """
        sp = next(iter(self.products))
        return sp.label

    def kr(self) -> float:
        """
        Returns the reaction constant.

        Returns:
            float
                The reaction constant.
        """
        return self._reaction_cst

    def ge(self) -> float:
        """
        Returns the reaction constant in Ge units.

        Returns:
            float
                The reaction constant in Ge units.
        """
        return kr_to_ge(self._reaction_cst)

    def update_from_ge(self, value: float) -> None:
        """
        Updates the reaction constant from a given G-value.

        Args:
            value (float): The G-value to convert.

        Returns:
            None
        """
        self._reaction_cst = ge_to_kr(value)

    def compute_initial_cc_before(self, t0: float, dose_rate: float) -> float:
        """
        Computes the production of the species by the beam for the first
        instants of beam ON i.e. when time <= t0.

        Arguments:
            t0: float
                Initial time of the Beam ON phase
            dose_rate: float
                Instantaneous dose rate of the related beam.

        Returns:
            float: initial concentration in mol/l
        """
        return self.kr() * dose_rate * t0

    @classmethod
    def from_dict(
        cls,
        product: str,
        ge: float,
        list_sp: tp.Optional[ListSpecies] = None,
    ) -> RadiolyticReaction:
        """
        Create a RadiolyticReaction from a dictionary
        (via destructuring)

        Args:
            product: str
                species to be produced
            ge: float
                G Value in [radical / 100eV / incident particle]
            list_sp: ListSpecies
                Collection of species already involved in the
                simulation. This collection will grow if some products
                do not exist yet.

        Returns:
            RadiolyticReaction
        """
        new = cls(ge_to_kr(ge))
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()
        new.add_as_product(list_sp.get_or_create_dyn_species(product.strip()))
        return new

    # ----------------------------------------------------------------------- #

    def add_as_reactant(self, sp: SimSpecies, multiple: int = 1) -> None:
        """No reactant should be involved"""
        raise NotImplementedError("Not implemented for Radiolytic Reactions")

    def __str__(self) -> str:
        """
        Returns a string representation of the RadiolyticReaction
        instance.
        """
        return f"RadiolyticReaction('{self.label()}')"

    def __repr__(self) -> str:
        """
        Returns a string representation of the object.

        Returns:
            str: A string representation of the object.
        """
        return self.__str__()


RadiolyticReactionDict = tp.Dict[str, float]
