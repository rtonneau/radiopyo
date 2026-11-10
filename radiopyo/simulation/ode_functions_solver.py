#!/usr/bin/env python3
"""
The ode_solver submodule contains all functions used to resolve ODE
system.

Every reaction rate is computed by the single pure function
``reaction_rates`` (used by both the ODE right hand side ``derive`` and
``ODEResult.reconstruct_reactions``). It is split in small steps that can
be reused when differentiating the system:
    1. mass action law of the k-reactions        -> mass_action_rates
    2. negative products blocking k-reactions     -> negative_product_block
    3. Michaelis-Menten law of enzymatic reactions -> michaelis_menten_rates
    4. radiolytic production (G values)           -> radiolytic_rates
The time derivative of the species is then the sum, over the reaction
families, of stoichiometry @ rates (see ``species_rates``).

Units: the full molar state and all rates are in [mol/L] and [mol/L/s].
Only ``derive`` deals with the [µmol/L] solver state, through
``SimMatrices.state`` (StateLayout).

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np

# ----------------------------- Internal imports ---------------------------- #
# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.physics.beam import ParticleBeam
    from radiopyo.utils.sim_types import (
        EzMatrices,
        GVector,
        KMatrices,
        SimMatrices,
    )

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "DEPLETION_LIMIT",
    "ReactionRates",
    "ScipyODEResult",
    "derive",
    "dose_rates",
    "mass_action_rates",
    "michaelis_menten_rates",
    "negative_product_block",
    "radiolytic_rates",
    "reaction_rates",
    "species_rates",
    "zero_depleted_consumption",
]

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #

# A species whose concentration is at or below this value [mol/L] is
# considered depleted.
DEPLETION_LIMIT: float = 1e-12


def dose_rates(
    t: float | NDArray[np.float64],
    beam: tp.Optional[ParticleBeam],
) -> NDArray[np.float64]:
    """Dose rate [Gy/s] of the beam at time(s) t [s].

    Args:
        t: one time (scalar) or an array of times.
        beam: ParticleBeam, or None for no irradiation (dose rate 0).

    Returns:
        NDArray[np.float64]: same shape as t (0-d array for a scalar t).
    """
    if not isinstance(t, np.ndarray):
        if beam is None:
            return np.asarray(0.0, dtype=np.float64)
        return np.asarray(beam.at(t).dose_rate(), dtype=np.float64)

    if beam is None:
        return np.zeros(t.shape, dtype=np.float64)
    out = [beam.at(t_val).dose_rate() for t_val in t.ravel()]
    return np.array(out, dtype=np.float64).reshape(t.shape)


def mass_action_rates(
    y_full: NDArray[np.float64],
    kmat: KMatrices,
) -> NDArray[np.float64]:
    """Mass action law of every k-reaction, ignoring negative products.

    Args:
        y_full: full molar state [mol/L], (n_species,) or (n_species, n_t).
        kmat: KMatrices

    Returns:
        NDArray[np.float64]: rates [mol/L/s], (n_kr,) or (n_kr, n_t).
    """
    # Unused reactant slots point to this extra row of 1.0
    ones = np.ones((1, *y_full.shape[1:]), dtype=np.float64)
    y_ext = np.concatenate([y_full, ones], axis=0)

    # (n_kr, max_reactants[, n_t]): concentration of each reactant molecule
    factors = y_ext[kmat.reactant_index]
    k_values = kmat.k_values if y_full.ndim == 1 else kmat.k_values[:, np.newaxis]
    return np.prod(factors, axis=1) * k_values


def negative_product_block(
    y_full: NDArray[np.float64],
    kmat: KMatrices,
) -> NDArray[np.bool_]:
    """Which k-reactions must not occur: the ones having at least one
    negative product that is depleted (concentration <= DEPLETION_LIMIT).

    Args:
        y_full: full molar state [mol/L], (n_species,) or (n_species, n_t).
        kmat: KMatrices

    Returns:
        NDArray[np.bool_]: True for blocked reactions, (n_kr,) or
        (n_kr, n_t).
    """
    depleted = y_full <= DEPLETION_LIMIT
    return kmat.negative_products.T @ depleted


def michaelis_menten_rates(
    y_full: NDArray[np.float64],
    ezmat: EzMatrices,
) -> NDArray[np.float64]:
    """Michaelis-Menten law of every enzymatic reaction.

    Args:
        y_full: full molar state [mol/L], (n_species,) or (n_species, n_t).
        ezmat: EzMatrices

    Returns:
        NDArray[np.float64]: rates [mol/L/s], (n_ez,) or (n_ez, n_t).
    """
    enzyme = y_full[ezmat.enzyme_index]
    substrate = y_full[ezmat.substrate_index]
    k_values, k_michaelis = ezmat.k_values, ezmat.k_michaelis
    if y_full.ndim > 1:
        k_values = k_values[:, np.newaxis]
        k_michaelis = k_michaelis[:, np.newaxis]
    return k_values * enzyme * substrate / (k_michaelis + substrate)


def radiolytic_rates(
    dose_rate: NDArray[np.float64],
    gvec: GVector,
) -> NDArray[np.float64]:
    """Radiolytic production rate of every radiolytic reaction.

    Args:
        dose_rate: dose rate(s) [Gy/s], 0-d or (n_t,) (see dose_rates).
        gvec: GVector

    Returns:
        NDArray[np.float64]: rates [mol/L/s], (n_g,) or (n_g, n_t).
    """
    return np.multiply.outer(gvec.yields, dose_rate)


def reaction_rates(
    t: float | NDArray[np.float64],
    y_full: NDArray[np.float64],
    mat: SimMatrices,
    beam: tp.Optional[ParticleBeam],
) -> ReactionRates:
    """Rate of every reaction, for one or several states.

    Pure function: same inputs, same output; nothing is modified.

    Args:
        t: time [s], a scalar for one state, or (n_t,) for n_t states.
        y_full: full molar state [mol/L], (n_species,) or (n_species, n_t)
            (see StateLayout.to_full_molar).
        mat: SimMatrices
        beam: ParticleBeam, or None for no irradiation.

    Returns:
        ReactionRates [mol/L/s]
    """
    k_rates = mass_action_rates(y_full, mat.k_reactions)
    k_rates = np.where(negative_product_block(y_full, mat.k_reactions), 0.0, k_rates)
    return ReactionRates(
        k_reactions=k_rates,
        radiolytic=radiolytic_rates(dose_rates(t, beam), mat.radiolytic),
        enzymatic=michaelis_menten_rates(y_full, mat.enzymatic),
    )


def species_rates(
    rates: ReactionRates,
    mat: SimMatrices,
) -> NDArray[np.float64]:
    """Time derivative of the full molar state [mol/L/s]: sum over the
    reaction families of stoichiometry @ rates. Rows of the constant
    species are 0.
    """
    dydt = mat.k_reactions.stoichiometry @ rates.k_reactions
    dydt += mat.enzymatic.stoichiometry @ rates.enzymatic
    dydt += mat.radiolytic.stoichiometry @ rates.radiolytic
    return dydt


def zero_depleted_consumption(
    y_full: NDArray[np.float64],
    dydt: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Set to 0 the (negative) derivative of every species already below
    DEPLETION_LIMIT, so that it cannot go further below 0."""
    return np.where((y_full < DEPLETION_LIMIT) & (dydt < 0), 0.0, dydt)


def derive(
    t: float,
    y: NDArray[np.float64],
    beam: tp.Optional[ParticleBeam],
    mat: SimMatrices,
    force_positive: bool = False,
) -> NDArray[np.float64]:
    """Compute the time derivative of species concentration (right hand
    side of the ODE system given to scipy.integrate.solve_ivp).
    See the RPM_extended_doc.pdf file for a more detailed explanation.

    Args:
        t: time [s]
        y: (n_dyn,) concentrations of the dynamical species [µmol/L]
        beam: ParticleBeam, or None for no irradiation.
        mat: SimMatrices
        force_positive: if True, depleted species are not consumed
            further (see zero_depleted_consumption).

    Returns:
        NDArray[np.float64]: (n_dyn,) dy/dt [µmol/L/s]
    """
    y_full = mat.state.to_full_molar(y)
    dydt = species_rates(reaction_rates(t, y_full, mat, beam), mat)
    if force_positive:
        dydt = zero_depleted_consumption(y_full, dydt)
    return mat.state.to_dynamic_micromolar(dydt)


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class ReactionRates(tp.NamedTuple):
    """Rates [mol/L/s] of every reaction, by reaction family. Each array
    is (n_reactions,) for one state, or (n_reactions, n_t) for n_t states.
    """

    k_reactions: NDArray[np.float64]
    radiolytic: NDArray[np.float64]
    enzymatic: NDArray[np.float64]

    def stacked(self) -> NDArray[np.float64]:
        """All rates in one array, ordered k-reactions, radiolytic then
        enzymatic (the column order of ODEResult.reconstruct_reactions;
        NOT the global reaction index order k -> enzymatic -> radiolytic).
        """
        return np.concatenate([self.k_reactions, self.radiolytic, self.enzymatic])


class ScipyODEResult(tp.Protocol):
    """
    Protocol mapping the output of scipy.integrate.solve_ivp
    """

    t: NDArray[np.float64]
    y: NDArray[np.float64]
    sol: tp.Optional[tp.Any]
    t_events: tp.Optional[tp.List[NDArray[np.float64]]]
    y_events: tp.Optional[tp.List[NDArray[np.float64]]]
    nfev: int
    njev: int
    nlu: int
    status: int
    message: str
    success: bool
