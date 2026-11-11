#!/usr/bin/env python3
"""
The ode_solver submodule contains all functions used to resolve ODE
system.

Every reaction rate is computed by the single pure function
``reaction_rates`` (used by both the ODE right hand side ``derive`` and
``ODEResult.reconstruct_reactions``). It is split in small steps that can
be reused when differentiating the system:
    1. mass action law of the k-reactions        -> mass_action_rates
    2. Michaelis-Menten law of enzymatic reactions -> michaelis_menten_rates
    3. radiolytic production (G values)           -> radiolytic_rates
    4. depletion switch of the negative products   -> depletion_switch
       (k-reactions and enzymatic reactions)
The time derivative of the species is then the sum, over the reaction
families, of stoichiometry @ rates (see ``species_rates``).

``jacobian`` is the analytic Jacobian of ``derive``: the same sum with the
rates replaced by their derivatives with respect to the concentrations
(mass_action_derivatives, michaelis_menten_derivatives,
depletion_switch_derivatives).

A k or enzymatic reaction with negative products (written ``- species`` on
the product side) is slowed down while one of them runs out: its rate is
multiplied, for each negative product, by the smooth depletion switch
c / (c + eps), with c the clipped concentration and eps = DEPLETION_LIMIT
(ADR 0005). The factor is 1/2 at c = eps and above 0.99 for c > 100 eps. It
replaces the former hard block at DEPLETION_LIMIT: ``derive`` is
continuous, and ``jacobian`` is no longer piecewise (it includes the
derivative of the switch, eps / (c + eps)**2, by the product rule).

Every rate law reads the clipped concentration max(c, 0) (ADR 0004): a
species driven below 0 by a solver overshoot counts as absent, so it is
not consumed further and cannot run a reaction backwards. The Jacobian is
only piecewise at the clipping kink c = 0 (ADR 0004).

Overshoot relaxation (ADR 0007): clipping alone leaves a stiff species
sitting near 0 with slope 0 below 0 and -k' above, a kink on which LSODA
stalls. A dynamic species j with c_j < 0 therefore gets, on its own row
only, the linear continuation from c = 0+:

    dc_j/dt = (stoichiometry @ rates)_j + J_jj(0+) * c_j

with J_jj(0+) the diagonal Jacobian entry computed with c_j read at 0+
(all the other species keep their value, a species below 0 reading as 0).
Every other row still reads c_j as 0, so no reaction runs backwards and
``derive`` is stoichiometry @ rates for a state without species below 0.
Only the columns of the negative species are computed (``overshoot_slopes``),
never a full Jacobian. ``jacobian`` adds the exact derivative of the term:
J_jj(0+) on the diagonal and c_j * d J_jj(0+) / d c_k on the row. The
reaction rates themselves (``reaction_rates``, hence
``ODEResult.reconstruct_reactions``) are unchanged.

Units: the full molar state and all rates are in [mol/L] and [mol/L/s].
Only ``derive`` deals with the [µmol/L] solver state, through
``SimMatrices.state`` (StateLayout).

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import dataclasses
import typing as tp

import numpy as np

# ----------------------------- Internal imports ---------------------------- #
# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.physics.beam import DoseRateSource
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
    "depletion_switch",
    "depletion_switch_derivatives",
    "derive",
    "dose_rates",
    "jacobian",
    "mass_action_derivatives",
    "mass_action_rates",
    "michaelis_menten_derivatives",
    "michaelis_menten_rates",
    "overshoot_slopes",
    "radiolytic_rates",
    "reaction_rates",
    "species_rates",
]

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #

# Scale eps [mol/L] of the depletion switch c / (c + eps): a negative product
# at this concentration halves the rate of its reaction.
DEPLETION_LIMIT: float = 1e-12


def dose_rates(
    t: float | NDArray[np.float64],
    beam: tp.Optional[DoseRateSource],
) -> NDArray[np.float64]:
    """Dose rate [Gy/s] of the beam at time(s) t [s].

    Args:
        t: one time (scalar) or an array of times.
        beam: ParticleBeam (or any DoseRateSource), or None for no
            irradiation (dose rate 0).

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
    """Mass action law of every k-reaction, without the depletion switch of
    the negative products (see ``depletion_switch``). Reads the clipped concentrations max(c, 0) (ADR 0004).

    Args:
        y_full: full molar state [mol/L], (n_species,) or (n_species, n_t).
        kmat: KMatrices

    Returns:
        NDArray[np.float64]: rates [mol/L/s], (n_kr,) or (n_kr, n_t).
    """
    # Unused reactant slots point to this extra row of 1.0
    ones = np.ones((1, *y_full.shape[1:]), dtype=np.float64)
    y_ext = np.concatenate([np.maximum(y_full, 0.0), ones], axis=0)

    # (n_kr, max_reactants[, n_t]): concentration of each reactant molecule
    factors = y_ext[kmat.reactant_index]
    k_values = kmat.k_values if y_full.ndim == 1 else kmat.k_values[:, np.newaxis]
    return np.prod(factors, axis=1) * k_values


def depletion_switch(
    y_full: NDArray[np.float64],
    mat: KMatrices | EzMatrices,
) -> NDArray[np.float64]:
    """Depletion switch of every reaction: the product, over its negative
    products, of c / (c + ε) with c = max(c, 0) the clipped concentration
    and ε = DEPLETION_LIMIT (ADR 0005). 1 for a reaction without negative
    product, tends to 0 as a negative product runs out.

    Args:
        y_full: full molar state [mol/L], (n_species,) or (n_species, n_t).
        mat: KMatrices or EzMatrices (both hold ``negative_products``).

    Returns:
        NDArray[np.float64]: factor in [0, 1] multiplying the rates,
        (n_r,) or (n_r, n_t).
    """
    negative = mat.negative_products
    out = np.ones((negative.shape[1], *y_full.shape[1:]), dtype=np.float64)
    columns = np.flatnonzero(negative.any(axis=0))
    if columns.size == 0:
        return out
    clipped = np.maximum(y_full, 0.0)
    per_species = clipped / (clipped + DEPLETION_LIMIT)
    for col in columns:
        out[col] = np.prod(per_species[negative[:, col]], axis=0)
    return out


def depletion_switch_derivatives(
    y_full: NDArray[np.float64],
    mat: KMatrices | EzMatrices,
) -> NDArray[np.float64]:
    """Derivative of the depletion switch of every reaction with respect to
    every species (product rule over the negative products).

    d/dc of c / (c + ε) is ε / (c + ε)² on the clipped concentration: 0 for
    c < 0 (ADR 0004); at exactly c = 0 the derivative of the c > 0 side
    (1 / ε) is used.

    Args:
        y_full: full molar state [mol/L], (n_species,).
        mat: KMatrices or EzMatrices (both hold ``negative_products``).

    Returns:
        NDArray[np.float64]: (n_r, n_species) d switch / d concentration
        [L/mol].
    """
    negative = mat.negative_products
    out = np.zeros((negative.shape[1], y_full.shape[0]), dtype=np.float64)
    columns = np.flatnonzero(negative.any(axis=0))
    if columns.size == 0:
        return out
    clipped = np.maximum(y_full, 0.0)
    per_species = clipped / (clipped + DEPLETION_LIMIT)
    slope = np.where(
        y_full < 0.0, 0.0, DEPLETION_LIMIT / (clipped + DEPLETION_LIMIT) ** 2
    )
    for col in columns:
        species = np.flatnonzero(negative[:, col])
        for pos, j in enumerate(species):
            others = np.prod(np.delete(per_species[species], pos))
            out[col, j] = others * slope[j]
    return out


def michaelis_menten_rates(
    y_full: NDArray[np.float64],
    ezmat: EzMatrices,
) -> NDArray[np.float64]:
    """Michaelis-Menten law of every enzymatic reaction, without the
    depletion switch of the negative products (see ``depletion_switch``).
    Reads the clipped concentrations max(c, 0) of the enzyme and the
    substrate (ADR 0004).

    Args:
        y_full: full molar state [mol/L], (n_species,) or (n_species, n_t).
        ezmat: EzMatrices

    Returns:
        NDArray[np.float64]: rates [mol/L/s], (n_ez,) or (n_ez, n_t).
    """
    enzyme = np.maximum(y_full[ezmat.enzyme_index], 0.0)
    substrate = np.maximum(y_full[ezmat.substrate_index], 0.0)
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
    beam: tp.Optional[DoseRateSource],
) -> ReactionRates:
    """Rate of every reaction, for one or several states.

    Pure function: same inputs, same output; nothing is modified.

    Args:
        t: time [s], a scalar for one state, or (n_t,) for n_t states.
        y_full: full molar state [mol/L], (n_species,) or (n_species, n_t)
            (see StateLayout.to_full_molar).
        mat: SimMatrices
        beam: ParticleBeam (or any DoseRateSource), or None for no irradiation.

    Returns:
        ReactionRates [mol/L/s]
    """
    return ReactionRates(
        k_reactions=mass_action_rates(y_full, mat.k_reactions)
        * depletion_switch(y_full, mat.k_reactions),
        radiolytic=radiolytic_rates(dose_rates(t, beam), mat.radiolytic),
        enzymatic=michaelis_menten_rates(y_full, mat.enzymatic)
        * depletion_switch(y_full, mat.enzymatic),
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


def mass_action_derivatives(
    y_full: NDArray[np.float64],
    kmat: KMatrices,
) -> NDArray[np.float64]:
    """Derivative of the mass action rate of every k-reaction with respect
    to every species, without the depletion switch of the negative products
    (see ``depletion_switch_derivatives``).

    Product rule on the reactant slots: d/dy_s of k * prod(slots) is k times
    the product of the other slots, summed over the slots holding s (a
    repeated reactant, 2 A, fills two slots).

    The rate reads the clipped concentrations max(c, 0) (ADR 0004): a
    species with c < 0 contributes 0 to every derivative it enters; at
    exactly c = 0 the derivative of the c > 0 side is used.

    Args:
        y_full: full molar state [mol/L], (n_species,).
        kmat: KMatrices

    Returns:
        NDArray[np.float64]: (n_kr, n_species) d rate / d concentration
        [1/s].
    """
    n_species = y_full.shape[0]
    n_kr, n_slots = kmat.reactant_index.shape
    # unused slots point to the extra 1.0
    y_ext = np.append(np.maximum(y_full, 0.0), 1.0)
    factors = y_ext[kmat.reactant_index]  # (n_kr, n_slots)

    rows = np.arange(n_kr)
    out = np.zeros((n_kr, n_species + 1), dtype=np.float64)
    for slot in range(n_slots):
        others = np.prod(np.delete(factors, slot, axis=1), axis=1) * kmat.k_values
        np.add.at(out, (rows, kmat.reactant_index[:, slot]), others)
    out[:, :n_species][:, y_full < 0.0] = 0.0  # clipped: no dependence below 0
    return out[:, :n_species]  # drop the extra column of the unused slots


def michaelis_menten_derivatives(
    y_full: NDArray[np.float64],
    ezmat: EzMatrices,
) -> NDArray[np.float64]:
    """Derivative of the Michaelis-Menten rate of every enzymatic reaction
    with respect to every species, without the depletion switch of the
    negative products (see ``depletion_switch_derivatives``).

    The rate reads the clipped concentrations max(c, 0) (ADR 0004): a
    species with c < 0 contributes 0; at exactly c = 0 the derivative of
    the c > 0 side is used.

    Args:
        y_full: full molar state [mol/L], (n_species,).
        ezmat: EzMatrices

    Returns:
        NDArray[np.float64]: (n_ez, n_species) d rate / d concentration
        [1/s].
    """
    raw_enzyme = y_full[ezmat.enzyme_index]
    raw_substrate = y_full[ezmat.substrate_index]
    enzyme = np.maximum(raw_enzyme, 0.0)
    substrate = np.maximum(raw_substrate, 0.0)
    denominator = ezmat.k_michaelis + substrate
    d_enzyme = np.where(raw_enzyme < 0.0, 0.0, ezmat.k_values * substrate / denominator)
    d_substrate = np.where(
        raw_substrate < 0.0,
        0.0,
        ezmat.k_values * enzyme * ezmat.k_michaelis / denominator**2,
    )

    rows = np.arange(ezmat.n_reactions)
    out = np.zeros((ezmat.n_reactions, y_full.shape[0]), dtype=np.float64)
    np.add.at(out, (rows, ezmat.enzyme_index), d_enzyme)
    np.add.at(out, (rows, ezmat.substrate_index), d_substrate)
    return out


def overshoot_slopes(
    y_full: NDArray[np.float64],
    negative: NDArray[np.intp],
    mat: SimMatrices,
) -> NDArray[np.float64]:
    """Overshoot relaxation slope J_jj(0+) [1/s] of the species below 0
    (ADR 0007): the diagonal Jacobian entry of species j with c_j read at
    0+ (every other species keeps its value, a species below 0 reading as
    0). It is 0 for a species no reaction consumes. Only the reactions
    that change species j are evaluated, one species at a time.

    Args:
        y_full: full molar state [mol/L], (n_species,).
        negative: indices of the dynamic species below 0.
        mat: SimMatrices

    Returns:
        NDArray[np.float64]: (len(negative),) J_jj(0+) [1/s].
    """
    return _overshoot(y_full, negative, mat, derivatives=False)[0]


def _overshoot(
    y_full: NDArray[np.float64],
    negative: NDArray[np.intp],
    mat: SimMatrices,
    derivatives: bool,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Slopes J_jj(0+) of the species ``negative`` and, if ``derivatives``,
    their gradients (len(negative), n_species): d J_jj(0+) / d c_k for k != j
    (0 for k = j: c_j is read at 0+; 0 for a k below 0, ADR 0004).
    """
    slopes = np.zeros(negative.size, dtype=np.float64)
    gradients = np.zeros((negative.size, y_full.shape[0]), dtype=np.float64)
    for pos, j in enumerate(negative):
        y0 = y_full.copy()
        y0[j] = 0.0
        for terms in (
            _k_relaxation_terms(y0, int(j), mat.k_reactions, derivatives),
            _ez_relaxation_terms(y0, int(j), mat.enzymatic, derivatives),
        ):
            if terms is None:
                continue
            slopes[pos] += terms[0]
            if terms[1] is not None:
                gradients[pos] += terms[1]
        gradients[pos, j] = 0.0
    return slopes, gradients


def _k_relaxation_terms(
    y0: NDArray[np.float64],
    j: int,
    kmat: KMatrices,
    derivatives: bool,
) -> tuple[float, tp.Optional[NDArray[np.float64]]] | None:
    """k-reactions part of ``_overshoot`` for species j: None when no
    k-reaction changes it. y0 is the state with c_j = 0."""
    cols = np.flatnonzero(kmat.stoichiometry[j])
    if cols.size == 0:
        return None
    sub = dataclasses.replace(
        kmat,
        k_values=kmat.k_values[cols],
        reactant_index=kmat.reactant_index[cols],
        stoichiometry=kmat.stoichiometry[:, cols],
        negative_products=kmat.negative_products[:, cols],
    )
    # d law / d c_j at c_j = 0+: only a reactant held once (k * product of the
    # other slots); held twice (2 A) it is 0, not held it is 0.
    slots = sub.reactant_index
    once = np.count_nonzero(slots == j, axis=1) == 1
    d_law = dataclasses.replace(
        sub,
        k_values=np.where(once, sub.k_values, 0.0),
        reactant_index=np.where(slots == j, y0.shape[0], slots),
    )
    return _relaxation_sum(
        y0,
        j,
        sub,
        law=mass_action_rates(y0, sub),
        law_j=mass_action_rates(y0, d_law),
        law_d=mass_action_derivatives(y0, sub) if derivatives else None,
        law_j_d=mass_action_derivatives(y0, d_law) if derivatives else None,
    )


def _ez_relaxation_terms(
    y0: NDArray[np.float64],
    j: int,
    ezmat: EzMatrices,
    derivatives: bool,
) -> tuple[float, tp.Optional[NDArray[np.float64]]] | None:
    """Enzymatic part of ``_overshoot`` for species j: None when no
    enzymatic reaction changes it. y0 is the state with c_j = 0."""
    cols = np.flatnonzero(ezmat.stoichiometry[j])
    if cols.size == 0:
        return None
    sub = dataclasses.replace(
        ezmat,
        k_values=ezmat.k_values[cols],
        k_michaelis=ezmat.k_michaelis[cols],
        substrate_index=ezmat.substrate_index[cols],
        enzyme_index=ezmat.enzyme_index[cols],
        stoichiometry=ezmat.stoichiometry[:, cols],
        negative_products=ezmat.negative_products[:, cols],
    )
    # law = k E S / (K + S) on the clipped E and S
    k, km = sub.k_values, sub.k_michaelis
    enzyme_slot, substrate_slot = sub.enzyme_index, sub.substrate_index
    enzyme = np.maximum(y0[enzyme_slot], 0.0)
    substrate = np.maximum(y0[substrate_slot], 0.0)
    denominator = km + substrate
    d_enzyme = k * substrate / denominator  # d law / d E
    d_substrate = k * enzyme * km / denominator**2  # d law / d S
    is_enzyme = enzyme_slot == j
    is_substrate = substrate_slot == j
    law_j = np.where(is_enzyme, d_enzyme, 0.0) + np.where(is_substrate, d_substrate, 0.0)

    law_d = law_j_d = None
    if derivatives:
        rows = np.arange(cols.size)
        d_es = k * km / denominator**2  # d2 law / dE dS
        d_ss = -2.0 * k * enzyme * km / denominator**3  # d2 law / dS2
        law_d = np.zeros((cols.size, y0.shape[0]), dtype=np.float64)
        np.add.at(law_d, (rows, enzyme_slot), d_enzyme)
        np.add.at(law_d, (rows, substrate_slot), d_substrate)
        law_j_d = np.zeros_like(law_d)
        np.add.at(law_j_d, (rows, enzyme_slot), np.where(is_substrate, d_es, 0.0))
        np.add.at(
            law_j_d,
            (rows, substrate_slot),
            np.where(is_enzyme, d_es, 0.0) + np.where(is_substrate, d_ss, 0.0),
        )
        law_d[:, y0 < 0.0] = 0.0  # clipped species: no dependence below 0
        law_j_d[:, y0 < 0.0] = 0.0
    return _relaxation_sum(
        y0,
        j,
        sub,
        law=k * enzyme * substrate / denominator,
        law_j=law_j,
        law_d=law_d,
        law_j_d=law_j_d,
    )


def _relaxation_sum(
    y0: NDArray[np.float64],
    j: int,
    sub: KMatrices | EzMatrices,
    law: NDArray[np.float64],
    law_j: NDArray[np.float64],
    law_d: tp.Optional[NDArray[np.float64]],
    law_j_d: tp.Optional[NDArray[np.float64]],
) -> tuple[float, tp.Optional[NDArray[np.float64]]]:
    """J_jj(0+) and its gradient from the reactions of ``sub`` (those that
    change species j), at the state y0 with c_j = 0.

    Rate of reaction i is law_i * switch_i. Its derivative with respect to
    c_j at 0+ is (product rule), with U_i the switch without the factor of
    species j (U_i is the whole switch if j is not a negative product):
      * j not a negative product of i:  d law_i / d c_j * U_i
      * j a negative product of i:      law_i * U_i / DEPLETION_LIMIT
        (the factor c_j / (c_j + eps) is 0 at 0+, its slope is 1/eps).
    So J_jj(0+) = sum_i S_ji * T_i * U_i with T_i the first or second
    expression without U_i; the gradient follows the product rule on
    T_i * U_i (``law_d`` and ``law_j_d`` are the derivatives of law and of
    d law / d c_j with respect to every species).

    Args:
        y0: (n_species,) full molar state with c_j = 0.
        j: species index.
        sub: the reactions of the family with S_ji != 0.
        law: (n_sub,) rate laws without switch.
        law_j: (n_sub,) d law / d c_j at 0+.
        law_d: (n_sub, n_species) d law / d c_k, or None (no gradient).
        law_j_d: (n_sub, n_species) d (d law / d c_j) / d c_k, or None.

    Returns:
        (J_jj(0+), its (n_species,) gradient or None).
    """
    flag = sub.negative_products[j]
    without_j = sub.negative_products.copy()
    without_j[j] = False
    reduced = dataclasses.replace(sub, negative_products=without_j)
    switch = depletion_switch(y0, reduced)
    weight = sub.stoichiometry[j]
    effective = np.where(flag, law / DEPLETION_LIMIT, law_j)
    slope = float(np.sum(weight * effective * switch))
    if law_d is None or law_j_d is None:
        return slope, None
    effective_d = np.where(flag[:, np.newaxis], law_d / DEPLETION_LIMIT, law_j_d)
    switch_d = depletion_switch_derivatives(y0, reduced)
    gradient = (weight * switch) @ effective_d + (weight * effective) @ switch_d
    return slope, gradient


def jacobian(
    t: float,
    y: NDArray[np.float64],
    beam: tp.Optional[DoseRateSource],
    mat: SimMatrices,
) -> NDArray[np.float64]:
    """Analytic Jacobian d(dy/dt)/dy of ``derive`` (the ``jac`` argument of
    scipy.integrate.solve_ivp).

    Pure function. The rate of a k or enzymatic reaction is the product of
    its rate law and of its depletion switch (ADR 0005), so every entry
    follows the product rule: d(law)/dc * switch + law * d(switch)/dc.
    It follows the clipped concentrations of ``derive`` (ADR 0004): a
    species below 0 contributes 0 to every derivative it enters, plus the
    derivative of its own overshoot relaxation term (ADR 0007, see
    ``overshoot_slopes``). Radiolytic
    reactions (they do not depend on the concentrations) and constant
    species (they are not solved) contribute nothing.

    The units cancel: the concentrations and their derivative are both
    scaled by 1e6, so the Jacobian is the same in [µmol/L] and [mol/L].

    Args:
        t: time [s]
        y: (n_dyn,) concentrations of the dynamical species [µmol/L]
        beam: ParticleBeam (or any DoseRateSource), or None for no irradiation.
        mat: SimMatrices

    Returns:
        NDArray[np.float64]: (n_dyn, n_dyn) [1/s], entry (i, j) is
        d (dy_i/dt) / d y_j.
    """
    y_full = mat.state.to_full_molar(y)
    n_dyn = y_full.shape[0] - mat.state.n_cst
    kmat, ezmat = mat.k_reactions, mat.enzymatic

    d_k = _product_rule(
        mass_action_rates(y_full, kmat),
        mass_action_derivatives(y_full, kmat),
        depletion_switch(y_full, kmat),
        depletion_switch_derivatives(y_full, kmat),
    )
    d_ez = _product_rule(
        michaelis_menten_rates(y_full, ezmat),
        michaelis_menten_derivatives(y_full, ezmat),
        depletion_switch(y_full, ezmat),
        depletion_switch_derivatives(y_full, ezmat),
    )

    jac = kmat.stoichiometry[:n_dyn] @ d_k[:, :n_dyn]
    jac += ezmat.stoichiometry[:n_dyn] @ d_ez[:, :n_dyn]

    # Overshoot relaxation (ADR 0007): the row of a species j below 0 gets
    # J_jj(0+) * c_j. d/dc_j is J_jj(0+) (the clipped column j is 0 above),
    # d/dc_k is c_j * d J_jj(0+) / d c_k.
    negative = np.flatnonzero(y < 0.0)
    if negative.size:
        slopes, gradients = _overshoot(y_full, negative, mat, derivatives=True)
        jac[negative] += y_full[negative, np.newaxis] * gradients[:, :n_dyn]
        jac[negative, negative] += slopes
    return jac


def _product_rule(
    law: NDArray[np.float64],
    law_derivatives: NDArray[np.float64],
    switch: NDArray[np.float64],
    switch_derivatives: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Derivatives (n_r, n_species) of the rates law * switch."""
    return (
        law_derivatives * switch[:, np.newaxis]
        + law[:, np.newaxis] * switch_derivatives
    )


def derive(
    t: float,
    y: NDArray[np.float64],
    beam: tp.Optional[DoseRateSource],
    mat: SimMatrices,
) -> NDArray[np.float64]:
    """Compute the time derivative of species concentration (right hand
    side of the ODE system given to scipy.integrate.solve_ivp).
    See the RPM_extended_doc.pdf file for a more detailed explanation.

    A dynamic species below 0 gets, on its own row, the overshoot
    relaxation term J_jj(0+) * c_j (ADR 0007, see the module docstring).

    Args:
        t: time [s]
        y: (n_dyn,) concentrations of the dynamical species [µmol/L]
        beam: ParticleBeam (or any DoseRateSource), or None for no irradiation.
        mat: SimMatrices

    Returns:
        NDArray[np.float64]: (n_dyn,) dy/dt [µmol/L/s]
    """
    y_full = mat.state.to_full_molar(y)
    dydt = species_rates(reaction_rates(t, y_full, mat, beam), mat)
    negative = np.flatnonzero(y < 0.0)
    if negative.size:
        dydt[negative] += overshoot_slopes(y_full, negative, mat) * y_full[negative]
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
