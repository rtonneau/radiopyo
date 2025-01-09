#!/usr/bin/env python3
"""
The ode_solver submodule contains all functions used to resolve ODE
system.

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

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "jacobian",
    "derive",
    "ScipyODEResult",
]

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #

global EPSILON
EPSILON = 1e-25  # µmol/L


def jacobian(
    t: float,
    y: NDArray[np.float64],
    beam: tp.Optional[ParticleBeam],
    H2O_index: int,
    H3O_plus_index: int,
    OH_minus_index: int,
    k_water: float,
    cst_sp: NDArray[np.float64],
    Ki_mat: tp.Tuple[
        NDArray[np.float64],
        NDArray[np.int8],
        NDArray[np.float32],
        NDArray[np.float32],
    ],
    Gi_mat: NDArray[np.float64],
    AB_mat: tp.Tuple[NDArray[np.int8], NDArray[np.float64]],
    Ez_mat: tp.Optional[
        tp.Tuple[
            NDArray[np.float32],
            NDArray[np.float32],
            NDArray[np.float64],
            NDArray[np.int8],
        ]
    ] = None,
) -> NDArray[np.float64]:
    """Computes the jacobian dy/ds_i where s_i are species[i]
    Returns:
        NDArray[np.float64]: jacobian (square) matrix
                 |   ds_0        ds_1    ....      ds_j
        ---------+---------------------------------------
            ds_0 | ds_0/ds_0   ds_0/ds_1 ....    ds_0/ds_j
            ds_1 | ds_1/ds_0   ds_1/ds_1 ....    ds_1/ds_j
            .... |    ....        ....   ....      ....
            ds_i | ds_i/ds_0   ds_i/ds_1 ....    ds_i/ds_j
    """

    y = np.append(y, cst_sp)

    # K Reactions component
    res = np.nan_to_num(np.einsum("j, jikl -> ikl ", y, Ki_mat[3]), nan=1)
    res = np.prod(res, axis=1).T * Ki_mat[0]
    res = np.matmul(Ki_mat[1], res.T)

    # Enzymatic Component
    if Ez_mat is not None:
        # En_subs, En_enzy, En_kmki => Enzi[0], Enzi[1], Enzi[2]
        # En_kmki[0,:] => Km & kmki[1,:] => k_val
        E = np.matmul(y, Ez_mat[1])
        S = np.matmul(y, Ez_mat[0]) + Ez_mat[2][0, :]
        JEnzi = np.prod(Ez_mat[2], axis=0) * E / S / S  # Vector...
        JEnzi = np.multiply(Ez_mat[0], JEnzi[None, :])
        JEnzi = np.matmul(Ez_mat[3], JEnzi.T)
        res += JEnzi

    # Take into account Acid/Base?

    ##
    return res[: -len(cst_sp), : -len(cst_sp)]


def derive(
    t: float,
    y: NDArray[np.float64],
    beam: tp.Optional[ParticleBeam],
    cst_sp: NDArray[np.float64],
    Ki_mat: tp.Tuple[
        NDArray[np.float64],
        NDArray[np.int8],
        NDArray[np.float32],
        NDArray[np.float32],
    ],
    Gi_mat: NDArray[np.float64],
    Ez_mat: tp.Optional[
        tp.Tuple[
            NDArray[np.float32],
            NDArray[np.float32],
            NDArray[np.float64],
            NDArray[np.int8],
        ]
    ] = None,
) -> NDArray[np.float64]:
    """Compute the time derivative of species concentration
    See the RPM_extended_doc.pdf file for a more detailed explanation.
    """

    # Force negative cc values to 0 + convert to [mol/L]
    y = y * 1.0e-6  # µmol/L => mol/l

    # Construct the cc vector with all relevant species
    # ==> Dynamical species (already in y)
    # ==> Add constant species
    y = np.append(y, cst_sp)

    ###########################################################################
    # Compute evolution rate of the species (contained in vector y)
    ###########################################################################

    # Compute the reaction rate of all k_reactions
    Rki = np.nan_to_num(np.einsum("j, jik -> ik ", y, Ki_mat[2]), nan=1)
    Rki = np.prod(Rki, axis=1) * Ki_mat[0]  # All reaction rates are here

    # Compute dydt, first only with chemistry:
    dydt = np.matmul(Ki_mat[1], Rki.T)

    # Add Enzymatic term
    if Ez_mat is not None:
        # En_subs, En_enzy, En_kmki => Enzi[0], Enzi[1], Enzi[2]
        # Extract Enzyme cc of each enzymatic reactions.
        E = np.matmul(y, Ez_mat[1])
        # Extract substrate cc of each enzymatic reactions.
        S = np.matmul(y, Ez_mat[0])
        # Compute the reaction rate of all enzymatic reactions.
        MM = Ez_mat[2][1, :] * E * S / (Ez_mat[2][0, :] + S)
        if np.isnan(MM).any():
            print("NaN found. E={E}; S={S}; MM={MM}")

        # Add this term to the total
        dydt += np.matmul(Ez_mat[3], MM.T)

    ###########################################################################
    if beam is not None:
        # Add the Time dependent Radiolytic term
        dydt += Gi_mat * beam.at(t).dose_rate()

    dydt[np.where((y < 0) & (dydt < 0))] = 0

    # ==> dynamical species
    # ==> convert back to [µmol/l]
    dydt = dydt[: -len(cst_sp)]

    # print(np.where((y < 0) & (dydt < 0)))

    return dydt * 1e6  # mol/l => µmol/l


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


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
