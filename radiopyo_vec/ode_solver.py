#!/usr/bin/env python3
"""
The ode_solver submodule contains all functions used to resolve ODE system.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    List,
    Optional,
    Protocol,
    Tuple,
)

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo_vec.physics.beam import ParticleBeam

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "jacobian",
    "derive",
    "ScipyODEResult",
    "ODEResult",
]

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
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
    """ Computes the jacobian dy/ds_i where s_i are species[i]
    Returns:
        NDArray[np.float64]: jacobian (square) matrix
                 |   ds_0        ds_1    ....      ds_j
        ---------|---------------------------------------
            ds_0 | ds_0/ds_0   ds_0/ds_1 ....    ds_0/ds_j
            ds_1 | ds_1/ds_0   ds_1/ds_1 ....    ds_1/ds_j
            .... |    ....        ....   ....      ....
            ds_i | ds_i/ds_0   ds_i/ds_1 ....    ds_i/ds_j
    """
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
    """ Compute the time derivative of species concentration
    See the RPM_extended_doc.pdf file for a more detailed explanation.
    """
    # Force negative cc values to 0 + convert to [mol/l]
    y = y.clip(min=0.0,) * 1e-6
    y = np.append(y, cst_sp)

    # Compute the reaction rate of all k_reactions
    Rki = np.nan_to_num(np.einsum("j, jik -> ik ", y, Ki_mat[2]), nan=1)

    # Compute dydt, first only with k reactions term:
    dydt = np.matmul(Ki_mat[1], (np.prod(Rki, axis=1) * Ki_mat[0]).T)

    # Add the Time dependent Radiolytic term
    dydt += Gi_mat * beam.at(t).dose_rate()

    # Add Enzymatic term
    if Ez_mat is not None:
        # En_subs, En_enzy, En_kmki => Enzi[0], Enzi[1], Enzi[2]
        # Extract Enzyme cc of each enzymatic reactions.
        E = np.matmul(y, Ez_mat[1])
        # Extract substrate cc of each enzymatic reactions.
        S = np.matmul(y, Ez_mat[0])
        # Compute the reaction rate of all enzymatic reactions.
        MM = Ez_mat[2][1, :] * E * S / (Ez_mat[2][0, :] + S)
        # Add this term to the total
        dydt += np.matmul(Ez_mat[3], MM.T)

    # Last step ==> Handling acid/base equilibrium
    # => Force relationship between the total dydt's for A/B Partners
    # ABy, ABi, ABo => ab_mat
    dydt = AB_mat[0]*dydt + np.matmul(AB_mat[2], np.matmul(dydt, AB_mat[1]).T)
    # select only the DynSpecies and convert to [µmol/l]
    return dydt[:-len(cst_sp)] * 1e6

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class ScipyODEResult(Protocol):
    """
    Protocol mapping the output of scipy.integrate.solve_ivp
    """
    t: NDArray[np.float64]
    y: NDArray[np.float64]
    sol: Optional[Any]
    t_events: Optional[List[NDArray[np.float64]]]
    y_events: Optional[List[NDArray[np.float64]]]
    nfev: int
    njev: int
    nlu: int
    status: int
    message: str
    success: bool


class ODEResult(object):
    """
    Collection storing the results from scipy.integrate.solve_ivp in order to easily
    append together the results from several integration steps.
    """
    # Properties of the bunch object sent by scipy.integrate.solve_ivp
    t: NDArray[np.float64]
    y: NDArray[np.float64]
    sol: List[Any]
    t_events: List[Optional[List[NDArray[np.float64]]]]
    y_events: List[Optional[List[NDArray[np.float64]]]]
    nfev: List[int]
    njev: List[int]
    nlu: List[int]
    status: List[int]
    message: List[str]
    success: List[bool]

    # Add new useful properties
    n_points: List[int]
    exec_time: List[Optional[float]]

    def __init__(self, sol: Optional[ScipyODEResult] = None) -> None:
        self._default_init()
        if sol is not None:
            self.append_sol(sol)

    def __len__(self):
        """Number of simulation run"""
        return len(self.n_points)

    @property
    def shape(self):
        return self.y.shape

    def _default_init(self):
        self.t = np.array([], dtype=np.float64)
        self.y = np.array([], dtype=np.float64)
        self.sol = []
        self.t_events, self.y_events = [], []
        self.nfev, self.njev, self.nlu = [], [], []
        self.status, self.message, self.success = [], [], []
        self.n_points, self.exec_time = [], []

    def append_sol(self,
                   sol: ScipyODEResult,
                   exec_time: Optional[float] = None,
                   ) -> None:
        self.t = np.append(self.t, sol.t)
        if len(self.y) == 0:
            self.y = sol.y
        else:
            self.y = np.hstack([self.y, sol.y])
        self.sol.append(sol.sol)
        self.t_events.append(sol.t_events)
        self.y_events.append(sol.y_events)
        self.nfev.append(sol.nfev)
        self.njev.append(sol.njev)
        self.nlu.append(sol.nlu)
        self.status.append(sol.status)
        self.message.append(sol.message)
        self.success.append(sol.success)
        self.n_points.append(len(sol.t))
        self.exec_time.append(exec_time)
