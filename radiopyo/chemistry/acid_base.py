#!/usr/bin/env python3
"""
The acid_base submodule contains all class definitions related to acid-base equilibrium.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 14/12/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import typing as tp
from contextlib import suppress

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.parser import RPARSER
from radiopyo.parser.reaction_parser import GenericReaction, revert_reaction

from .collections import ListSpecies
from .k_reactions import (
    BiReaction,
    KReaction,
    polynomial_factory,
)
from .utils.units import Temperature

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from .species import ChemSpecies


__all__ = [
    "AcidBase",
    "ABCouple",
]

# ------------------------------------------------------------------------------------ #
#                                   CLASS DEFINITIONS                                  #
# ------------------------------------------------------------------------------------ #


class ABCouple(tp.NamedTuple):
    """
    Attributes:
        species: ChemSpecies
            Acid or Base
        conjugate: ChemSpecies
            Conjugate base from Acid or Conjugate acid from Base.
        pKa: float
        k_forward: float
        k_backward: float
        temperature: float | None
    """

    species: ChemSpecies
    conjugate: ChemSpecies
    pKa: float
    k_forward: float
    k_backward: float
    temperature: None | float


class AcidBase(BiReaction):
    """
    This one is more tricky about temperature management as various polynomial forms
    are involved:
        pKa = a0 + a1 T + a2 T**2 + ...
        log(K) = a0 + a1 (1/T) + a2 (1/T)**2 + ...
    => pKa uses T while log(K) uses (1/T)
    Therefore, the calculation of reaction constants cannot be left to the respective
    KReactions as these are not aware of it.

    Rem:
    K_equilibrium = K_forward / K_backward
    <=> log(K_e) = log(K_f) - log(K_b)
    K_equilibrium = 10 ** (-pKa)
    """

    _k_backward: np.polynomial.Polynomial
    _k_forward: np.polynomial.Polynomial
    _pKa: np.polynomial.Polynomial

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

        self._k_backward = polynomial_factory(reaction["constants"]["k_backward"])
        self._k_forward = polynomial_factory(reaction["constants"]["k_forward"])

        # Need to add forward and backward reactions via the superclass
        self.add_forward_reaction(reaction, list_sp=list_sp)
        self.add_backward_reaction(revert_reaction(reaction), list_sp=list_sp)

    @classmethod
    def from_kwargs(
        cls,
        list_sp: tp.Optional[ListSpecies] = None,
        **kwargs: str | tp.List[str] | float | tp.List[float],
    ) -> tp.Self:
        """
        This method combines the two possible definitions of a KReaction in a config
        file.
        TODO: remove the dict approach to keep only the single string reaction.
        """
        # Allow shallow creation (outside of a simulation env)
        if list_sp is None:
            list_sp = ListSpecies()

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
        """ """
        _T = Temperature(T) if T is not None else Temperature()
        return self._pKa(_T.celsius)

    def k_equilibrium(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)

        For pKa, polynomial form depends on T in °C.
        """
        _T = Temperature(T) if T is not None else Temperature()
        return 10 ** (-self._pKa(_T.celsius))

    def k_forward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        _T = Temperature(T) if T is not None else Temperature()
        return 10 ** self._k_forward(1.0 / _T.kelvin)

    def k_backward(self, T: tp.Optional[Temperature | float] = None) -> float:
        """
        Args:
            T: Temperature in °C (optional, default:25)
        """
        _T = Temperature(T) if T is not None else Temperature()
        return 10 ** self._k_backward(1.0 / _T.kelvin)

    def is_acid(self) -> bool:
        if "H3O_plus" in self.forward or "H_plus" in self.forward:
            return True
        return False

    def is_base(self) -> bool:
        if "OH_minus" in self.forward:
            return True
        return False

    def get_deprotonation_reaction(self) -> KReaction:
        """ """
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
        Ugly function but working function
        ==> TODO: Refactor
        """
        r = self.get_deprotonation_reaction()

        if self.is_acid():
            _l = len(r.reactants)
            if _l > 2:
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. The reaction has more "
                    f"than 2 different reactants ({_l} reactants)"
                )
            if _l == 2 and not r.has_reactant("H2O"):
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. The reaction has more "
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
                    f"Unable to identify acid/base couple of {r}. The reaction has more "
                    f"than 2 different products ({_l} products)"
                )
            if _l == 2 and not (r.has_product("H3O_plus") or r.has_product("H_plus")):
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. The reaction has more "
                    "than 1 product but none of them is H3O_plus or H_plus."
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
                    f"Unable to identify acid/base couple of {r}. The reaction has more "
                    f"than 2 different products ({_l} products)"
                )
            if _l == 2 and not r.has_product("H2O"):
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. The reaction has more "
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
                    f"Unable to identify acid/base couple of {r}. The reaction has more "
                    f"than 2 different reactants ({_l} reactants)"
                )
            if _l == 2 and not r.has_reactant("OH_minus"):
                raise ValueError(
                    f"Unable to identify acid/base couple of {r}. The reaction has more "
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
