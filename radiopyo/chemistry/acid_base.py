#!/usr/bin/env python3
"""
The acid_base submodule contains all class definitions related to
acid-base equilibrium.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 14/12/2023
"""
# ---------------------------- External imports ----------------------------- #
from __future__ import annotations

import typing as tp
from contextlib import suppress

import numpy as np

# ---------------------------- Internal imports ----------------------------- #
from radiopyo.chemistry.k_reactions import (
    BiReaction,
    KReaction,
)
from radiopyo.chemistry.species_collection import ListSpecies
from radiopyo.chemistry.utils.factories import polynomial_factory
from radiopyo.chemistry.utils.units import Temperature
from radiopyo.parser import RPARSER
from radiopyo.parser.reaction_parser import GenericReaction, revert_reaction

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.species import ChemSpecies


__all__ = [
    "AcidBase",
    "ABCouple",
]

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class ABCouple(tp.NamedTuple):
    """
    Represents an acid-base couple.

    Attributes:
        species (ChemSpecies): The acid or base species.
        conjugate (ChemSpecies): The conjugate base of the acid or
            conjugate acid of the base.
        pKa (float): The acid dissociation constant (pKa) of the couple.
        k_forward (float): The forward reaction rate constant
            of the couple.
        k_backward (float): The backward reaction rate constant
            of the couple.
        temperature (float | None): The temperature at which the
            couple's properties are defined (in Kelvin). If None, the
            temperature is not specified.

    Note:
        This class is a NamedTuple, which means it is immutable and its
        attributes can be accessed by name or by index.
    """

    species: ChemSpecies
    conjugate: ChemSpecies
    pKa: float
    k_forward: float
    k_backward: float
    temperature: None | float


class AcidBase(BiReaction):
    """
    Represents an acid-base reaction.

    This class is a subclass of BiReaction and manages the
    temperature-dependent reaction constants for acid-base reactions.

    Attributes:
        _pKa (np.polynomial.Polynomial): The polynomial form of the pKa.
        _k_backward (np.polynomial.Polynomial): The polynomial form of
            the backward reaction rate constant.
        _k_forward (np.polynomial.Polynomial): The polynomial form of
            the forward reaction rate constant.
        _HA (ChemSpecies): The acid species.
        _A (ChemSpecies): The conjugate base species.

    Notes:
        The reaction constants are calculated using the polynomial forms
        of the pKa and reaction rate constants:
            pKa = a0 + a1 T + a2 T**2 + ...
            log(K) = a0 + a1 (1/T) + a2 (1/T)**2 + ...
        Therefore, the calculation of reaction constants cannot be left
        to the respective KReactions as these are not aware of it.
        The temperature is assumed to be in Kelvin.
    """

    _k_backward: np.polynomial.Polynomial
    _k_forward: np.polynomial.Polynomial
    _pKa: np.polynomial.Polynomial
    _A: ChemSpecies
    _HA: ChemSpecies

    def __init__(
        self,
        reaction: GenericReaction,
        pKa: float | tp.List[float],
        list_sp: tp.Optional[ListSpecies] = None,
    ):
        """ """
        super().__init__()
        self._pKa = np.polynomial.Polynomial(pKa)

        if reaction["constants"].get("k_backward", None) is None:
            raise ValueError("For A/B reaction, k_backward must be defined")
        if reaction["constants"].get("k_forward", None) is None:
            raise ValueError("For A/B reaction, k_forward must be defined")

        self._k_backward = polynomial_factory(
            reaction["constants"]["k_backward"],
            log10=True,
        )
        self._k_forward = polynomial_factory(
            reaction["constants"]["k_forward"],
            log10=True,
        )

        # Need to add forward and backward reactions via the superclass
        self.add_forward_reaction(reaction, list_sp=list_sp)
        self.add_backward_reaction(revert_reaction(reaction), list_sp=list_sp)

        _ = self.get_couple()
        self._HA = _.species
        self._A = _.conjugate

    @property
    def acid(self) -> ChemSpecies:
        """
        A property to return the acid species.

        Returns:
            ChemSpecies: The acid species.
        """
        return self._HA

    @property
    def base(self) -> ChemSpecies:
        """
        A property to return the base species.

        Returns:
            ChemSpecies: The base species.
        """
        return self._A

    @classmethod
    def from_kwargs(
        cls,
        list_sp: tp.Optional[ListSpecies] = None,
        **kwargs: str | tp.List[str] | float | tp.List[float],
    ) -> tp.Self:
        """
        A class method to create an AcidBase reaction from keyword
        arguments. This method combines the two possible definitions of
        a KReaction in a config file.

        Args:
            list_sp (ListSpecies): A collection of species already
                involved in the simulation.
            **kwargs: Keyword arguments containing the reaction string
                and pKa value.

        Returns:
            AcidBase: An instance of the AcidBase class.

        Raises:
            AttributeError: If the reaction string or pKa value is
                missing.
            AttributeError: If 'k_value' is ambiguous for an AcidBase
                reaction.
        """
        # Allow shallow creation (outside of a simulation env)
        list_sp = list_sp if list_sp is not None else ListSpecies()

        try:
            gen = RPARSER.parse_reaction(kwargs["reaction"])  # type:ignore[arg-type]
            pKa = kwargs["pKa"]
        except KeyError as err:
            raise AttributeError(
                "Need a reaction string AND a pKa to create an Acid Base rel."
            ) from err

        if "k_value" in kwargs:
            raise AttributeError(
                "'k_value' is ambiguous for an AcidBase reaction. "
                "Use k_forward or k_backward instead."
            )

        with suppress(KeyError):
            gen["constants"]["k_forward"] = kwargs[
                "k_forward"
            ]  # type:ignore[assignment]
        with suppress(KeyError):
            gen["constants"]["k_backward"] = kwargs[
                "k_backward"
            ]  # type:ignore[assignment]

        return cls(gen, pKa=pKa, list_sp=list_sp)  # type:ignore[arg-type]

    def pKa(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Returns the acid dissociation constant (pKa) of the acid-base
        reaction at a given temperature.

        Args:
            T (Temperature | float, optional): Temperature in °C.
                Defaults to None.

        Returns:
            float: The pKa value at the given temperature.
        """
        _T = Temperature(T) if T is not None else Temperature()
        return self._pKa(_T.celsius)

    def k_equilibrium(
        self,
        T: tp.Optional[Temperature | float] = None,
    ) -> float:
        """
        Returns the acid-base equilibrium constant at a given
        temperature.

        Args:
            T (Temperature | float, optional): Temperature in °C.
                Defaults to None.

        Returns:
            float: The equilibrium constant at the given temperature.
        """
        _T = Temperature(T) if T is not None else Temperature()
        return 10 ** (-self._pKa(_T.celsius))

    def k_forward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Returns the forward reaction constant at a given temperature.

        Args:
            T (Temperature | float, optional): Temperature in °C.
                Defaults to None.

        Returns:
            float: The forward reaction constant at the
                given temperature.
        """
        _T = Temperature(T) if T is not None else Temperature()
        if self.forward.has_reactant("H2O"):
            return 10 ** self._k_forward(1.0 / _T.kelvin)  # / 55.5
        return 10 ** self._k_forward(1.0 / _T.kelvin)

    def k_backward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Returns the backward reaction constant at a given temperature.

        Args:
            T (Temperature | float, optional): Temperature in °C.
                Defaults to None.

        Returns:
            float: The backward reaction constant at the
                given temperature.
        """

        _T = Temperature(T) if T is not None else Temperature()
        if self.backward.has_reactant("H2O"):
            return 10 ** self._k_backward(1.0 / _T.kelvin)  # / 55.5
        return 10 ** self._k_backward(1.0 / _T.kelvin)

    def is_acid(self) -> bool:
        """
        Checks whether the reaction involved H3O+ or H+.

        Returns:
            bool: True if the reaction is an acid, False otherwise.
        """
        return bool("H3O_plus" in self.forward or "H_plus" in self.forward)

    def is_base(self) -> bool:
        """
        Checks whether the reaction involved OH-.

        Returns:
            bool: True if the reaction is a base, False otherwise.
        """
        return "OH_minus" in self.forward

    def get_deprotonation_reaction(self) -> KReaction:
        """
        Returns the deprotonation reaction of the current reaction.

        The deprotonation reaction is determined based on the nature of
        the reaction. If the reaction involves H3O+ or H+, it returns
        the reaction that produces H3O+ or H+. If the reaction involves
        OH-, it returns the reaction that has OH- as a reactant.

        Returns:
            KReaction: The deprotonation reaction.

        Raises:
            ValueError: If the reaction is neither an acid nor a base.
        """
        if self.is_acid():
            if self.forward.has_product("H3O_plus") or self.forward.has_product(
                "H_plus"
            ):
                return self.forward
            return self.backward

        if self.is_base():
            if self.forward.has_reactant("OH_minus"):
                return self.forward
            return self.backward

        raise ValueError(f"Unable to identify Acid or Base nature of {self}")

    def get_couple(self, T: tp.Optional[Temperature | float] = None) -> ABCouple:
        """
        Returns the acid/base couple of the current reaction.

        The acid/base couple is determined based on the nature of the
        reaction. If the reaction involves H3O or H+, it returns the
        reaction that produces H3O+ or H+. If the reaction involves OH-,
        it returns the reaction that has OH- as a reactant.

        Args:
            T: The temperature in Celsius. If None, the temperature is
                not considered.

        Returns:
            ABCouple: The acid/base couple.

        Raises:
            ValueError: If the reaction is neither an acid nor a base.

        Notes:
            To refactor
        """
        r = self.get_deprotonation_reaction()

        if self.is_acid():
            _l = len(r.reactants)
            if _l > 2:
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. "
                    "The reaction has more "
                    f"than 2 different reactants ({_l} reactants)"
                )
            if _l == 2 and not r.has_reactant("H2O"):
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. "
                    "The reaction has more "
                    "than 1 reactant but none of them is H2O."
                )

            species = None
            for sp, _ in r.iter_reactants():
                if sp == "H2O":
                    continue
                species = sp
            else:
                species = sp if species is None else species

            _l = len(r.products)
            if _l > 2:
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. "
                    "The reaction has more "
                    f"than 2 different products ({_l} products)"
                )
            if (_l == 2) and (
                not (r.has_product("H3O_plus") or r.has_product("H_plus"))
            ):
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. "
                    "The reaction has more than 1 product but none "
                    "of them is H3O_plus or H_plus."
                )

            conjugate = None
            for sp, _ in r.iter_products():
                if sp == "H3O_plus" or sp == "H_plus":
                    continue
                conjugate = sp
            else:
                conjugate = sp if conjugate is None else conjugate
        else:
            _l = len(r.products)
            if _l > 2:
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. "
                    "The reaction has more "
                    f"than 2 different products ({_l} products)"
                )
            if _l == 2 and not r.has_product("H2O"):
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. "
                    "The reaction has more "
                    "than 1 product but none of them is H2O."
                )

            conjugate = None
            for sp, _ in r.iter_products():
                if sp == "H2O":
                    continue
                conjugate = sp
            else:
                conjugate = sp if conjugate is None else conjugate

            _l = len(r.reactants)
            if _l > 2:
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. "
                    "The reaction has more "
                    f"than 2 different reactants ({_l} reactants)"
                )
            if _l == 2 and not r.has_reactant("OH_minus"):
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. "
                    "The reaction has more "
                    "than 1 reactant but none of them is OH_minus."
                )

            species = None
            for sp, _ in r.iter_reactants():
                if sp == "OH_minus":
                    continue
                species = sp
            else:
                species = sp if species is None else species

        return ABCouple(
            species=species,
            conjugate=conjugate,
            pKa=self.pKa(T),
            k_forward=self.k_forward(T),
            k_backward=self.k_backward(T),
            temperature=(Temperature(T).celsius if T is not None else None),
        )
