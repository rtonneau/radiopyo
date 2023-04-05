# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Optional,
    Tuple,
)

import numpy as np
import pandas as pd  # type: ignore[import]

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo_vec.physics.beam import ParticleBeam


# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #
def jacobian(t: float,
             y: NDArray[np.float64],
             beam: ParticleBeam,
             cst_sp: NDArray[np.float64],
             AB_mat: Tuple[NDArray[np.int8],
                           NDArray[np.int8],
                           NDArray[np.float64]],
             Ki_mat: Tuple[NDArray[np.float64],
                           NDArray[np.int8],
                           NDArray[np.float32],
                           NDArray[np.float32]],
             Gi_mat: NDArray[np.float64],
             Ez_mat: Optional[Tuple[NDArray[np.float32],
                                    NDArray[np.float32],
                                    NDArray[np.float64],
                                    NDArray[np.int8],
                                    ]] = None,
             ) -> NDArray[np.float64]:
    y = y.clip(min=0.0,) * 1e-6
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

    # Take into account Acid/Base
    Jab = AB_mat[0] + AB_mat[2].sum(axis=1)
    res = np.multiply(res, Jab[None, :])

    return res[:-len(cst_sp), :-len(cst_sp)]


def derive(t: float,
           y: NDArray[np.float64],
           beam: ParticleBeam,
           cst_sp: NDArray[np.float64],
           AB_mat: Tuple[NDArray[np.int8],
                         NDArray[np.int8],
                         NDArray[np.float64]],
           Ki_mat: Tuple[NDArray[np.float64],
                         NDArray[np.int8],
                         NDArray[np.float32],
                         NDArray[np.float32]],
           Gi_mat: NDArray[np.float64],
           Ez_mat: Optional[Tuple[NDArray[np.float32],
                                  NDArray[np.float32],
                                  NDArray[np.float64],
                                  NDArray[np.int8],
                                  ]] = None,
           ) -> NDArray[np.float64]:

    # Force negative cc values to 0
    y = y.clip(min=0.0,) * 1e-6
    y = np.append(y, cst_sp)
    # Build u vector => [1, y(t), y(t),...] n times y(t) with n = max
    # stoichiometry of reactants specie s (max_stoi).
    # max_stoi, Ri, Riy => ki_mat
    Rki = np.nan_to_num(np.einsum("j, jik -> ik ", y, Ki_mat[2]), nan=1)

    # Compute dydt, first only with k reactions term:
    dydt = np.matmul(Ki_mat[1], (np.prod(Rki, axis=1) * Ki_mat[0]).T)

    # Add Radiolytic term ==> Need deposition rate!
    # => Only term function of time... probably need refactorization => TODO
    dydt += Gi_mat * beam.at(t).dose_rate()

    # Add Enzymatic term
    if Ez_mat is not None:
        # En_subs, En_enzy, En_kmki => Enzi[0], Enzi[1], Enzi[2]
        # Extract substrate cc of each enzymatic reaction.
        E = np.matmul(y, Ez_mat[1])
        S = np.matmul(y, Ez_mat[0])
        MM = Ez_mat[2][1, :] * E * S / (Ez_mat[2][0, :] + S)
        dydt += np.matmul(Ez_mat[3], MM.T)

    # Last step ==> Handling acid/base equilibrium
    # => Force relation ship between total dydt for A/B Partners
    # ABy, ABi, ABo => ab_mat
    dydt = AB_mat[0]*dydt + np.matmul(AB_mat[2], np.matmul(dydt, AB_mat[1]).T)
    return dydt[:-len(cst_sp)] * 1e6
