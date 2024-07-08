#!/usr/bin/env python3
"""

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 26/07/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import copy
import datetime
import json
import logging
import typing as tp
from dataclasses import dataclass, field
from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import simpson  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.utils.units import Temperature
from radiopyo.physics.beam import BeamCollection
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.ode_functions_solver import ScipyODEResult, jacobian

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:

    from numpy.typing import NDArray

    from radiopyo.physics.beam import ParticleBeam
    from radiopyo.simulation.utils import H5Manager
    from radiopyo.utils.sim_types import SimMatrices

    T = tp.TypeVar("T", bound="ODEResult")

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__: tp.List = []

# -------------------------------------- Logging ------------------------------------- #
logger = logging.getLogger("radiopyo")

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
    T: Temperature

    # Properties of the bunch object sent by scipy.integrate.solve_ivp
    t: NDArray[np.float64]
    y: NDArray[np.float64]
    dydt: NDArray[np.float64]
    sol: tp.List[tp.Any]
    t_events: tp.List[tp.List[NDArray[np.float64]]]
    y_events: tp.List[tp.List[NDArray[np.float64]]]
    nfev: tp.List[int]
    njev: tp.List[int]
    nlu: tp.List[int]
    status: tp.List[int]
    message: tp.List[str]
    success: tp.List[bool]

    # Add new useful properties
    n_points: tp.List[int]
    exec_time: tp.List[tp.Optional[float]]

    def __init__(
        self,
        env: SimEnv,
        beam: ParticleBeam,
        mat: SimMatrices,
        let: tp.Optional[float],
        sol: tp.Optional[ScipyODEResult] = None,
        temperature: tp.Optional[Temperature] = None,
    ) -> None:
        self.date = datetime.datetime.now()
        self.env = copy.deepcopy(env)
        self.beam = copy.deepcopy(beam)
        self.mat = copy.deepcopy(mat)
        self.LET = let
        self.T = temperature if temperature is not None else Temperature()
        self._default_init()
        if sol is not None:
            self.append_sol(sol)

    def __len__(self) -> int:
        """Number of simulation run"""
        return len(self.n_points)

    def __getitem__(self, key: str) -> pd.Series:
        """ """
        if key not in self.env.species:
            raise KeyError(f"'{key}' is not a valid species for this env.")

        # Is it a constant species?
        if not self.env.species.is_dyn_species(key):
            return pd.Series(
                data=self.env.species[key].initial_cc() * 1e6, index=self.t, name=key
            )

        sp, idx = self.env.species.get_sim_species(key)
        return pd.Series(data=self.y[idx, :], index=self.t, name=sp.label)

    def __add__(self, other: tp.Any) -> ODEResult:
        """ """
        # Sanity checks
        if not isinstance(other, ODEResult):
            raise TypeError(f"Cannot add type: {type(other)}, to ODEResult")
        if other.t[0] < self.t[-1]:
            raise TypeError(
                "Trying to add ODEResult with starting "
                "time inferior to previous results"
            )
        if other.y.shape[0] != self.y.shape[0]:
            raise TypeError(
                "Trying to add ODEResults with different number of tracked species"
            )

        # Combine all together
        self.t = np.append(self.t, other.t)
        self.y = np.hstack([self.y, other.y])
        if other.sol is not None:
            self.sol = self.sol + other.sol
        if other.t_events is not None:
            self.t_events += other.t_events
        if other.y_events is not None:
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
        """ """
        return copy.copy(self.y[:, -1])

    @property
    def initial_cc(self) -> NDArray[np.float64]:
        """ """
        return self.env.initial_values(all=True)

    def to_hdf(
        self, file: str | Path, name: str, group: tp.Optional[str] = None
    ) -> H5Manager:
        """ """
        from .utils import H5Manager

        h5file = H5Manager.from_path(file)
        with h5file:
            grp = h5file.create_group(group)
            data = self.to_pandas().reset_index()
            ds = h5file.insert_dataset(grp, name, data.values)
            labels = [str(_) for _ in data.columns]
            labels[0] = "time"
            ds.attrs["columns"] = labels
            ds.attrs["sol"] = self.sol
            ds.attrs["t_events"] = self.t_events
            ds.attrs["y_events"] = self.y_events
            ds.attrs["nfev"] = self.nfev
            ds.attrs["njev"] = self.njev
            ds.attrs["nlu"] = self.nlu
            ds.attrs["status"] = self.status
            ds.attrs["message"] = self.message
            ds.attrs["success"] = self.success
            ds.attrs["n_points"] = self.n_points
            ds.attrs["exec_time"] = self.exec_time
        return h5file

    def _default_init(self) -> None:
        self.t = np.array([], dtype=np.float64)
        self.y = np.array([], dtype=np.float64)
        self.dydt = np.array([], dtype=np.float64)
        self.sol = []
        self.t_events, self.y_events = [], []
        self.nfev, self.njev, self.nlu = [], [], []
        self.status, self.message, self.success = [], [], []
        self.n_points, self.exec_time = [], []

    def append_sol(
        self,
        sol: ScipyODEResult,
        exec_time: tp.Optional[float] = None,
    ) -> None:
        self.t = np.append(self.t, sol.t)
        if len(self.y) == 0:
            self.y = sol.y
        else:
            self.y = np.hstack([self.y, sol.y])
        if sol.sol is not None:
            self.sol.append(sol.sol)
        if sol.t_events is not None:
            self.t_events.append(sol.t_events)
        if sol.y_events is not None:
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
        Provides simulation results in a nice pandas DataFrame.
        Units are µM.
        """
        labels = [elt.label for elt in self.env.species.sorted_iter_dyn_species()]
        labels_cst = [
            elt.label for elt in self.env.species.sorted_iter_cst_species()]

        df_sp = pd.DataFrame(
            index=self.t,
            columns=labels,
            dtype=np.float64,
            data=self.y.T,
        ).clip(0.0 if clip else None)

        df_cst = pd.DataFrame(
            index=self.t,
            columns=labels_cst,
            dtype=np.float64,
        )
        for sp_cst in self.env.species.iter_cst_species():
            df_cst.loc[:, sp_cst.label] = (
                sp_cst.initial_cc() * 1e6
            )  # TODO include temperature!

        return pd.concat([df_sp, df_cst], axis=1)

    def reconstruct_jacobian(self) -> pd.DataFrame:
        raise NotImplementedError
        if self.beam is None:
            raise TypeError("Beam is None, unable to compute jacobian")

        sp_label = [elt.label for elt in self.env.species.iter_dyn_species()]
        dsp_label = ["".join(["d", elt]) for elt in sp_label]

        out = pd.DataFrame(
            index=self.t,
            columns=pd.MultiIndex.from_product(
                [sp_label, dsp_label], names=["sp", "dsp"]
            ),
            dtype=np.float64,
        )

        for idx, t_val in enumerate(self.t):
            jac = jacobian(
                t_val,
                self.y[:, idx],
                self.beam,
                self.env.species["H_plus"].index,
                **self.mat,
            )
            out.iloc[idx, :] = jac.flatten()
        return out

    def reconstruct_reactions(self, recompute: bool = False) -> pd.DataFrame:
        """
        returns:
            dataframe [mol/l/s]
        Rem: to integrate df
            >>> df.apply(lambda x: simpson(y=x.values, x=x.index))
        """
        # Compute column names and Pre compute slices to ease access in the dataframe
        # columns
        cols_kr = [kr.as_label()
                   for kr in self.env.reactions.iter_all_kreactions()]
        idx_kr = slice(0, len(cols_kr))

        cols_gi = [gr.as_label() for gr in self.env.reactions.radiolytic]
        idx_gi = slice(idx_kr.stop, idx_kr.stop + len(cols_gi))

        if self.mat["Ez_mat"] is not None:
            cols_ez = [er.as_label() for er in self.env.reactions.enzymatic]
            idx_ez = slice(idx_gi.stop, idx_gi.stop + len(cols_ez))
        else:
            cols_ez = []

        if self.dydt.shape[-1] == self.y.shape[-1] > 0 and not recompute:
            return pd.DataFrame(
                columns=list(chain.from_iterable([cols_kr, cols_gi, cols_ez])),
                index=self.t,
                dtype=np.float64,
                data=self.dydt.T,
            )

        out = pd.DataFrame(
            columns=list(chain.from_iterable([cols_kr, cols_gi, cols_ez])),
            index=self.t,
            dtype=np.float64,
        )

        # Prefill values for radiolytic reactions
        for label, reaction in zip(cols_gi, self.env.reactions.radiolytic):
            out.loc[:, label] = reaction.kr()

        for (idx, t_val), y in zip(enumerate(self.t), self.y.T):
            y = y * 1.0e-6  # µmol/L => mol/l

            y = np.append(y, self.mat["cst_sp"])
            # K Reactions
            Rki = np.nan_to_num(
                np.einsum("j, jik -> ik ", y, self.mat["Ki_mat"][2]), nan=1
            )
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
                MM = (
                    self.mat["Ez_mat"][2][1, :]
                    * E
                    * S
                    / (self.mat["Ez_mat"][2][0, :] + S)
                )
                out.iloc[idx, idx_ez] = MM
        self.dydt = out.values.T
        return out

    def species_balance_evolution(self) -> pd.DataFrame:
        """ """
        rr = self.reconstruct_reactions()  # if stored then no "heavy calculation"
        out = pd.DataFrame(
            index=self.t,
            columns=[elt.label for elt in self.env.species.iter_dyn_species()],
            dtype=np.float64,
        )
        out.iloc[:, :] = 0

        for species in self.env.species.iter_dyn_species():
            for reaction in self.env.reactions_involving_product(species):
                stoi = reaction.stoi_products[species.label]
                out.loc[:, species.label] += rr[reaction.as_label()] * stoi
            for reaction in self.env.reactions_involving_reactant(species):
                stoi = reaction.stoi_reactants[species.label]
                out.loc[:, species.label] -= rr[reaction.as_label()] * stoi
        return out

    def integrate_species(
        self, start: float = 0, stop: float = float("inf"), compute_AB: bool = True
    ) -> pd.Series:
        """ """
        df = self.to_pandas(compute_AB)
        df = df[(df.index >= start) & (df.index < stop)]
        out = {}

        for col, data in df.items():
            out[col] = simpson(y=data.values, x=data.index)
        return pd.Series(out, dtype=np.float64)

    def to_json_config(self, file_path: Path | str, overwrite: bool = False) -> None:
        """ """
        file_path = Path(file_path)
        if not overwrite and file_path.exists():
            raise FileExistsError(file_path)
        out = self.env.convert_to_config_dict()
        out["parameters"]["timestamp"] = self.date.timestamp()
        out["parameters"]["temperature"] = self.T.kelvin
        out["parameters"]["species_order"] = {
            "dynamical": [
                sp.label for sp in self.env.species.sorted_iter_dyn_species()
            ],
            "constant": [sp.label for sp in self.env.species.sorted_iter_cst_species()],
        }
        kwargs = {}
        if self.LET is not None:
            kwargs["LET"] = self.LET
        out["beam"] = self.beam.serialize(**kwargs)
        with file_path.open("w") as fp:
            json.dump(out, fp, indent=4)

    @classmethod
    def from_json_config(cls: tp.Type[tp.Self], file_path: str | Path) -> tp.Self:
        """
            def __init__(
        self,
        env: SimEnv,
        beam: ParticleBeam,
        mat: SimMatrices,
        let: tp.Optional[float],
        sol: tp.Optional[ScipyODEResult] = None,
        """
        _prev_log_level = logger.level
        logger.setLevel(logging.ERROR)

        file_path = Path(file_path)
        with file_path.open("r") as fp:
            config = json.load(fp)
        env = SimEnv.from_dict(config)
        T = (
            Temperature(config["parameters"]["temperature"], unit="kelvin")
            if "temperature" in config["parameters"]
            else Temperature()
        )
        beam_type, beam_config = next(iter(config["beam"].items()))
        label, beam, LET = BeamCollection.build_beam_from_dict(
            beam_type, beam_config)
        mat = env.build_all_matrices()
        out = cls(env, beam, mat, LET, temperature=T)
        out.date = datetime.datetime.fromtimestamp(
            config["parameters"]["timestamp"])

        logger.setLevel(_prev_log_level)
        return out


@dataclass
class DummyScipyOdeResult:
    """ """

    t: NDArray[np.float64]
    y: NDArray[np.float64]
    sol: tp.Optional[tp.Any] = field(repr=False, default_factory=list)
    t_events: tp.Optional[tp.List[NDArray[np.float64]]] = field(
        repr=False,
        default=None,
    )
    y_events: tp.Optional[tp.List[NDArray[np.float64]]] = field(
        repr=False,
        default=None,
    )
    nfev: int = field(repr=False, default=0)
    njev: int = field(repr=False, default=0)
    nlu: int = field(repr=False, default=0)
    status: int = field(repr=False, default=0)
    message: str = field(repr=False, default="")
    success: bool = field(repr=False, default=True)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, order: tp.List[str]) -> tp.Self:
        df = df.reindex(columns=order)
        return cls(
            t=df.index.values,
            y=df.values.T,
        )
