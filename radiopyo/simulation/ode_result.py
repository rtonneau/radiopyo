#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 26/07/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import copy
import datetime
import typing as tp
from itertools import chain

import numpy as np
import pandas as pd
from scipy.integrate import simpson  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.simulation.ode_functions_solver import ScipyODEResult, jacobian

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

    from radiopyo.physics.beam import ParticleBeam
    from radiopyo.simulation.environment_build import SimEnv
    from radiopyo.utils.sim_types import SimMatrices

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = []

# -------------------------------------- Logging ------------------------------------- #

# ------------------------------------------------------------------------------------ #
#                                DECORATOR DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                 FUNCTION DEFINITIONS                                 #
# ------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class ODEResult(object):
    """
    Collection storing the results from scipy.integrate.solve_ivp in order to easily
    append together the results from several integration steps.
    """
    # Keep track of the original sim env
    env: SimEnv
    beam: ParticleBeam
    mat: SimMatrices
    LET: tp.Optional[float]
    date: datetime.datetime

    # Properties of the bunch object sent by scipy.integrate.solve_ivp
    t: NDArray[np.float64]
    y: NDArray[np.float64]
    sol: tp.List[tp.Any]
    t_events: tp.List[tp.Optional[tp.List[NDArray[np.float64]]]]
    y_events: tp.List[tp.Optional[tp.List[NDArray[np.float64]]]]
    nfev: tp.List[int]
    njev: tp.List[int]
    nlu: tp.List[int]
    status: tp.List[int]
    message: tp.List[str]
    success: tp.List[bool]

    # Add new useful properties
    n_points: tp.List[int]
    exec_time: tp.List[tp.Optional[float]]

    def __init__(self,
                 env: SimEnv,
                 beam: ParticleBeam,
                 mat: SimMatrices,
                 let: tp.Optional[float],
                 sol: tp.Optional[ScipyODEResult] = None,
                 ) -> None:
        self.date = datetime.datetime.now()
        self.env = copy.deepcopy(env)
        self.beam = copy.deepcopy(beam)
        self.mat = copy.deepcopy(mat)
        self.LET = let
        self._default_init()
        if sol is not None:
            self.append_sol(sol)

    def __len__(self) -> int:
        """Number of simulation run"""
        return len(self.n_points)

    def __getitem__(self, key: str) -> pd.Series:
        """ """
        if not self.env.species.is_dyn_species(key):
            raise KeyError(f"'{key}' is not a valid DynSpecies label")
        idx = self.env.species.dyn_sp[key].index
        return pd.Series(data=self.y[idx, :], index=self.t)

    def __add__(self, other: tp.Any) -> ODEResult:
        """ """
        # Sanity checks
        if not isinstance(other, ODEResult):
            raise TypeError(f"Cannot add type: {type(other)}, to ODEResult")
        if other.t[0] < self.t[-1]:
            raise TypeError("Trying to add ODEResult with starting "
                            "time inferior to previous results")
        if other.y.shape[0] != self.y.shape[0]:
            raise TypeError("Trying to add ODEResults with different number of tracked"
                            "species")

        # Combine all together
        self.t = np.append(self.t, other.t)
        self.y = np.hstack([self.y, other.y])
        self.sol = self.sol + other.sol
        self.t_events += other.t_events
        self.y_events += other.y_events
        self.nfev += other.nfev
        self.njev += other.njev
        self.nlu += other.nlu
        self.status += other.status
        self.message += other.message
        self.success += other.success
        self.n_points += other.n_points
        self.exec_time += other.exec_time
        return self

    @property
    def time(self) -> NDArray[np.float64]:
        return self.t

    @property
    def shape(self) -> tp.Tuple:
        return self.y.shape

    @property
    def final_cc(self) -> NDArray[np.float64]:
        """
        """
        return copy.copy(self.y[:, -1])

    @property
    def initial_cc(self) -> NDArray[np.float64]:
        """
        """
        return self.env.initial_values(all=True)

    def _default_init(self) -> None:
        self.t = np.array([], dtype=np.float64)
        self.y = np.array([], dtype=np.float64)
        self.sol = []
        self.t_events, self.y_events = [], []
        self.nfev, self.njev, self.nlu = [], [], []
        self.status, self.message, self.success = [], [], []
        self.n_points, self.exec_time = [], []

    def append_sol(self,
                   sol: ScipyODEResult,
                   exec_time: tp.Optional[float] = None,
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

    # --------------------------------- Utils methods -------------------------------- #
    def to_pandas(self, clip: bool = True) -> pd.DataFrame:
        """
        """
        labels = [(elt.label, elt.index)
                  for elt in self.env.species.iter_dyn_species()]
        labels = sorted(labels, key=lambda x: x[1])
        return pd.DataFrame(index=self.t,
                            columns=[elt[0] for elt in labels],
                            dtype=np.float64,
                            data=self.y.T,
                            ).clip(0.0 if clip else None)

    def reconstruct_jacobian(self) -> pd.DataFrame:
        if self.beam is None:
            raise TypeError("Beam is None, unable to compute jacobian")

        sp_label = [elt.label for elt in self.env.iter_DynSpecies()]
        dsp_label = ["".join(["d", elt]) for elt in sp_label]

        out = pd.DataFrame(index=self.t,
                           columns=pd.MultiIndex.from_product([sp_label, dsp_label],
                                                              names=["sp", "dsp"]),
                           dtype=np.float64,)

        for (idx, t_val) in enumerate(self.t):
            jac = jacobian(t_val, self.y[:, idx],
                           self.beam,
                           **self.mat,
                           )
            out.iloc[idx, :] = jac.flatten()
        return out

    def reconstruct_reactions(self) -> pd.DataFrame:
        """ 
        Rem: to integrate df
            >>> df.apply(lambda x: simpson(y=x.values, x=x.index))
        """
        # Compute column names and Pre compute slices to ease access in the dataframe
        # columns
        cols_kr = [kr.as_label() for kr in self.env.reactions.k_reaction]
        idx_kr = slice(0, len(cols_kr))

        cols_gi = [gr.as_label() for gr in self.env.reactions.radiolytic]
        idx_gi = slice(idx_kr.stop, idx_kr.stop+len(cols_gi))

        if self.mat["Ez_mat"] is not None:
            cols_ez = [er.as_label() for er in self.env.reactions.enzymatic]
            idx_ez = slice(idx_gi.stop, idx_gi.stop+len(cols_ez))
        else:
            cols_ez = []

        out = pd.DataFrame(columns=list(chain.from_iterable([cols_kr,
                                                             cols_gi,
                                                             cols_ez])),
                           index=self.t,
                           dtype=np.float64,
                           )

        # Prefill values for radiolytic reactions
        for label, reaction in zip(cols_gi, self.env.reactions.radiolytic):
            out.loc[:, label] = reaction.kr()

        for ((idx, t_val), y) in zip(enumerate(self.t), self.y.T):
            y = y.clip(min=0.0,) * 1e-6
            y = np.append(y, self.mat["cst_sp"])

            # K Reactions
            Rki = np.nan_to_num(
                np.einsum("j, jik -> ik ", y, self.mat["Ki_mat"][2]), nan=1)
            Rki = np.prod(Rki, axis=1) * self.mat["Ki_mat"][0]
            out.iloc[idx, idx_kr] = Rki

            # Radiolytic Reactions
            dr = 0.0
            if self.beam is not None:
                dr = self.beam.at(t_val).dose_rate()
            out.iloc[idx, idx_gi] *= dr

            # Enzimatic reactions
            if self.mat["Ez_mat"] is not None:
                E = np.matmul(y, self.mat["Ez_mat"][1])
                S = np.matmul(y, self.mat["Ez_mat"][0])
                MM = self.mat["Ez_mat"][2][1, :] * E * \
                    S / (self.mat["Ez_mat"][2][0, :] + S)
                out.iloc[idx, idx_ez] = MM
        return out

    def species_balance_evolution(self) -> pd.DataFrame:
        rr = self.reconstruct_reactions()
        out = pd.DataFrame(index=self.t,
                           columns=[elt.label
                                    for elt in self.env.iter_DynSpecies()],
                           dtype=np.float64,
                           )
        out.iloc[:, :] = 0

        for species in self.env.iter_DynSpecies():
            for reaction in self.env.reactions_involving_product(species):
                stoi = reaction.stoi_products[species.label]
                out.loc[:, species.label] += rr[reaction.as_label()] * stoi
            for reaction in self.env.reactions_involving_reactant(species):
                stoi = reaction.stoi_reactants[species.label]
                out.loc[:, species.label] -= rr[reaction.as_label()] * stoi
        return out

    def integrate_species(self,
                          start: float = 0,
                          stop: float = float("inf"),
                          ) -> pd.Series:
        """ 
        """
        labels = [(elt.label, elt.index)
                  for elt in self.env.species.iter_dyn_species()]
        out = pd.Series(index=[elt[0] for elt in sorted(labels, key=lambda x: x[1])],
                        dtype=np.float64)
        mask = np.where((self.t >= start) & (self.t < stop))
        for idx in range(len(self.y)):
            out.iloc[idx] = simpson(
                y=self.y[idx, mask].clip(min=0), x=self.t[mask])
        return out
