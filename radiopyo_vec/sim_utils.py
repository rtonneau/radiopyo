# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from itertools import chain
from typing import (
    TYPE_CHECKING,
    Any,
)

import numpy as np
import pandas as pd

# --------------------------------- Internal Imports --------------------------------- #
from radiopyo_vec.ode_solver import ODEResult, jacobian

if TYPE_CHECKING:
    from radiopyo_vec.sim_hand import UnitCell

__all__ = [
    "SimUtils",
]


class SimUtils(object):
    """Class gathering useful operations on sim results.
    This one should contains only static methods"""

    @staticmethod
    def reconstruct_reactions(sim: UnitCell, sol: Any) -> pd.DataFrame:
        # Compute column names and Pre compute slices to ease access in the dataframe
        # columns
        cols_kr = [kr.as_label() for kr in sim.env.reactions.k_reaction]
        idx_kr = slice(0, len(cols_kr))

        cols_gi = [gr.as_label() for gr in sim.env.reactions.radiolytic]
        idx_gi = slice(idx_kr.stop, idx_kr.stop+len(cols_gi))

        if sim.mat["Ez_mat"] is not None:
            cols_ez = [er.as_label() for er in sim.env.reactions.enzymatic]
            idx_ez = slice(idx_gi.stop, idx_gi.stop+len(cols_ez))
        else:
            cols_ez = []

        out = pd.DataFrame(columns=list(chain.from_iterable([cols_kr,
                                                             cols_gi,
                                                             cols_ez])),
                           index=sol.t,
                           dtype=np.float64,
                           )

        # Prefill values for radiolytic reactions
        for label, reaction in zip(cols_gi, sim.env.reactions.radiolytic):
            out.loc[:, label] = reaction.kr()

        for ((idx, t_val), y) in zip(enumerate(sol.t), sol.y.T):
            y = y.clip(min=0.0,) * 1e-6
            y = np.append(y, sim.mat["cst_sp"])

            # K Reactions
            Rki = np.nan_to_num(
                np.einsum("j, jik -> ik ", y, sim.mat["Ki_mat"][2]), nan=1)
            Rki = np.prod(Rki, axis=1) * sim.mat["Ki_mat"][0]
            out.iloc[idx, idx_kr] = Rki

            # Radiolytic Reactions
            dr = 0.0
            if sim.beam is not None:
                dr = sim.beam.at(t_val).dose_rate()
            out.iloc[idx, idx_gi] *= dr

            # Enzimatic reactions
            if sim.mat["Ez_mat"] is not None:
                E = np.matmul(y, sim.mat["Ez_mat"][1])
                S = np.matmul(y, sim.mat["Ez_mat"][0])
                MM = sim.mat["Ez_mat"][2][1, :] * E * \
                    S / (sim.mat["Ez_mat"][2][0, :] + S)
                out.iloc[idx, idx_ez] = MM
        return out

    @staticmethod
    def species_balance_evolution(sim: UnitCell, sol: Any) -> pd.DataFrame:
        rr = SimUtils.reconstruct_reactions(sim, sol)
        out = pd.DataFrame(index=sol.t,
                           columns=[elt.label
                                    for elt in sim.env.iter_DynSpecies()],
                           dtype=np.float64,
                           )
        out.iloc[:, :] = 0

        for species in sim.env.iter_DynSpecies():
            for reaction in sim.env.reactions_involving_product(species):
                stoi = reaction.stoi_products[species.label]
                out.loc[:, species.label] += rr[reaction.as_label()] * stoi
            for reaction in sim.env.reactions_involving_reactant(species):
                stoi = reaction.stoi_reactants[species.label]
                out.loc[:, species.label] -= rr[reaction.as_label()] * stoi
        return out

    @staticmethod
    def reconstruct_jacobian(sim: UnitCell, sol: Any) -> pd.DataFrame:
        if sim.beam is None:
            raise TypeError("Beam is None, unable to compute jacobian")
        sp_label = [elt.label for elt in sim.env.iter_DynSpecies()]
        dsp_label = ["".join(["d", elt]) for elt in sp_label]

        out = pd.DataFrame(index=sol.t,
                           columns=pd.MultiIndex.from_product([sp_label, dsp_label],
                                                              names=["sp", "dsp"]),
                           dtype=np.float64,)

        for (idx, t_val) in enumerate(sol.t):
            jac = jacobian(t_val, sol.y[:, idx],
                           sim.beam,
                           **sim.mat,
                           )
            out.iloc[idx, :] = jac.flatten()
        return out

    @staticmethod
    def to_pandas(sim: UnitCell, sol: ODEResult) -> pd.DataFrame:
        """
        """
        labels = [elt.label for elt in sim.env.species.iter_dyn_species()]
        return pd.DataFrame(index=sol.t,
                            columns=labels,
                            dtype=np.float64,
                            data=sol.y.T)
