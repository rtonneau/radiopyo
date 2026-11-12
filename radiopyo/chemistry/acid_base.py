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

import numpy as np

# ---------------------------- Internal imports ----------------------------- #
from radiopyo.chemistry.k_reactions import (
    BiReaction,
    KReaction,
)
from radiopyo.chemistry.utils.factories import polynomial_factory
from radiopyo.chemistry.utils.units import Temperature

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.chemistry.basic_reaction import (
        ReactionSpecies,
        ReactionSpeciesCollection,
    )
    from radiopyo.chemistry.generics import GenericReaction
    from radiopyo.chemistry.species import SimSpecies
    from radiopyo.chemistry.species_collection import ListSpecies


__all__ = [
    "ABCouple",
    "AcidBase",
]

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class ABCouple(tp.NamedTuple):
    """
    Represents an acid-base couple.

    Attributes:
        species (SimSpecies): The acid or base species.
        conjugate (SimSpecies): The conjugate base of the acid or
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

    species: SimSpecies
    conjugate: SimSpecies
    pKa: float
    k_forward: float
    k_backward: float
    temperature: float | None


class AcidBase(BiReaction):
    """
    Represents an acid-base reaction.

    This class is a subclass of BiReaction and manages the
    temperature-dependent reaction constants for acid-base reactions.
    The forward and backward rate constants are owned by the forward
    and backward KReactions.

    Attributes:
        _pKa (np.polynomial.Polynomial): The polynomial form of the pKa.
        _HA (SimSpecies): The acid species.
        _A (SimSpecies): The conjugate base species.

    Notes:
        The constants are calculated using the polynomial forms of the
        pKa and of the reaction rate constants:
            pKa = a0 + a1 T + a2 T**2 + ...      (T in °C)
            log(K) = a0 + a1 (1/T) + a2 (1/T)**2 + ...   (T in Kelvin)
        The configured k_forward/k_backward are used as given (no n·k
        division, see ADR 0001).
    """

    _pKa: np.polynomial.Polynomial
    _A: SimSpecies
    _HA: SimSpecies

    def __init__(
        self,
        reaction: GenericReaction,
        pKa: float | tp.List[float],
        list_sp: tp.Optional[ListSpecies] = None,
    ):
        """ """
        super().__init__()
        self._pKa = np.polynomial.Polynomial(tp.cast("list[float]", pKa))

        if reaction["constants"].get("k_backward", None) is None:
            raise ValueError("For A/B reaction, k_backward must be defined")
        if reaction["constants"].get("k_forward", None) is None:
            raise ValueError("For A/B reaction, k_forward must be defined")

        self._add_reaction_pair(
            reaction,
            k_forward=polynomial_factory(
                reaction["constants"]["k_forward"], log10=True
            ),
            k_backward=polynomial_factory(
                reaction["constants"]["k_backward"], log10=True
            ),
            list_sp=list_sp,
        )

        _ = self.get_couple()
        self._HA = _.species
        self._A = _.conjugate

    @property
    def acid(self) -> SimSpecies:
        """
        A property to return the acid species.

        Returns:
            SimSpecies: The acid species.
        """
        return self._HA

    @property
    def base(self) -> SimSpecies:
        """
        A property to return the base species.

        Returns:
            SimSpecies: The base species.
        """
        return self._A

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
        _T = Temperature(T)
        return tp.cast(float, self._pKa(_T.celsius))

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
        _T = Temperature(T)
        return tp.cast(float, 10 ** (-self._pKa(_T.celsius)))

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
            ValueError: If the reaction is neither an acid nor a base, or
                if the couple cannot be identified.
        """
        r = self.get_deprotonation_reaction()

        if self.is_acid():
            # HA (+ H2O) -> A (+ H3O_plus | H_plus)
            species = self._couple_member(r, r.reactants, "reactant", ("H2O",))
            conjugate = self._couple_member(
                r, r.products, "product", ("H3O_plus", "H_plus")
            )
        else:
            # HA (+ OH_minus) -> A (+ H2O)
            conjugate = self._couple_member(r, r.products, "product", ("H2O",))
            species = self._couple_member(r, r.reactants, "reactant", ("OH_minus",))

        return ABCouple(
            species=species.species,
            conjugate=conjugate.species,
            pKa=self.pKa(T),
            k_forward=self.k_forward(T),
            k_backward=self.k_backward(T),
            temperature=(Temperature(T).celsius if T is not None else None),
        )

    @staticmethod
    def _couple_member(
        r: KReaction,
        collection: ReactionSpeciesCollection,
        kind: str,
        partners: tp.Tuple[str, ...],
    ) -> ReactionSpecies:
        """
        Find the member of the acid/base couple on one side of the
        deprotonation reaction 'r': the side holds at most 2 species,
        and if 2, one of them is a partner (H2O, H3O_plus, ...). The
        member is the (last) non-partner species.

        Args:
            r (KReaction): The deprotonation reaction (error messages).
            collection (ReactionSpeciesCollection): r.reactants or
                r.products.
            kind (str): "reactant" or "product" (error messages).
            partners (Tuple[str, ...]): Labels of the partner species.

        Returns:
            ReactionSpecies: The couple member.

        Raises:
            ValueError: If the couple member cannot be identified.
        """
        _l = len(collection)
        if _l > 2:
            raise ValueError(
                f"Unable to identify acid/base couple of {r}. "
                "The reaction has more "
                f"than 2 different {kind}s ({_l} {kind}s)"
            )
        if _l == 2 and not any(label in collection for label in partners):
            raise ValueError(
                f"Unable to identify acid/base couple of {r}. "
                f"The reaction has more than 1 {kind} but none "
                f"of them is {' or '.join(partners)}."
            )
        members = [sp for sp in collection if sp.label not in partners]
        return members[-1] if members else list(collection)[-1]

    def to_dict(self) -> tp.Dict[str, tp.Any]:
        """
        Convert the reaction to a dictionary.

        Returns:
            dict    A dictionary representation of the reaction.
        """
        return {
            "reaction": self.__str__(),
            "pKa": self.pKa(),
            "k_forward": self._raw_k(self.forward, self.k_forward),
            "k_backward": self._raw_k(self.backward, self.k_backward),
        }

    @staticmethod
    def _raw_k(
        reaction: KReaction,
        k_at: tp.Callable[[], float],
    ) -> float | tp.List[float]:
        """
        The rate constant as configured: the log10 polynomial
        coefficients if temperature dependent, else the k value.
        """
        poly = reaction.k_polynomial
        if poly is not None and len(poly) > 1:
            return tp.cast("list[float]", poly.coef.tolist())
        return k_at()
