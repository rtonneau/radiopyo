#!/usr/bin/env python3
"""
The sim_hand submodule contains class(es) to easily handle simulation. This one should
be the entry point for lambda user.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import collections
import copy
import logging
import time
from contextlib import suppress
from itertools import chain, repeat
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Iterable,
    List,
    NamedTuple,
    Optional,
    Tuple,
    Type,
    TypedDict,
)

import numpy as np
import pandas as pd
from scipy.integrate import simpson, solve_ivp  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
from . import exceptions as rp_excep
from .chemistry.species import CstSpecies
from .input.g_values import GValuesBoscolo
from .ode_solver import ScipyODEResult, derive, jacobian
from .physics.beam import NoBeam, beam_from_dict
from .sim_env import SimEnv
from .utils.ron_parser import RonFileParser
from .utils.toml_parser import TOMLFileParser

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from .physics.beam import ParticleBeam
    from .sim_types import ConfigDict, FileParser

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "UnitCell",
    "SimChunk",
    "SimChunksIter",
]

# ------------------------------------------------------------------------------------ #
logger = logging.getLogger("radiopyo")

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class SimMatrices(TypedDict):
    """Convenient TypedDict use to leverage type checking and store matrices used in ODE
    solver."""
    AB_mat: Tuple[NDArray[np.int8], NDArray[np.int8], NDArray[np.float64]]
    Ki_mat: Tuple[NDArray[np.float64], NDArray[np.int8],
                  NDArray[np.float32], NDArray[np.float32]]
    Gi_mat: NDArray[np.float64]
    Ez_mat: Optional[Tuple[NDArray[np.float32], NDArray[np.float32],
                           NDArray[np.float64], NDArray[np.int8],]]
    cst_sp: NDArray[np.float64]
    O2_mat: Optional[Tuple[NDArray[np.float64], NDArray[np.float64]]]


class UnitCell(object):
    """ 
    Class UnitCell represents the tiniest simulation volume. All concentrations are
    assumed to be homogenous in it. It also defines method to solve ODE's in order to
    compute the evolution of concentration over time.

    Attributes:
        env: SimEnv
            SimEnv build from the input file
        mat: SimMatrices
            Collections of matrices used to solve the ODE's
        beam: ParticleBeam
            Radiation source to use for the simulation
        last_exec_time: float
            Time duration of last simulation step
        g_val: GValuesBoscolo
            GValues from the literature to possibly update the LET
        FIRST_STEP: float
            Initial time step of all simulation (default=1e-12)

    """
    FIRST_STEP: float = 1e-12

    env: SimEnv
    mat: SimMatrices
    beam: ParticleBeam
    last_exec_time: float
    g_val: GValuesBoscolo
    LET: Optional[float]

    def __init__(self,
                 env: SimEnv,
                 beam: Optional[ParticleBeam] = None,
                 ) -> None:
        self.env = env
        self.mat = SimMatrices(
            AB_mat=self.env.make_acidBase_matrix(),
            Ki_mat=self.env.make_kreaction_matrix(),
            Gi_mat=self.env.make_radiolytic_vector(),
            Ez_mat=self.env.make_enzymatic_matrix(),
            cst_sp=self.env.make_cst_species_vector(),
            O2_mat=self.env.make_O2_intake_matrix(),
        )

        self.beam = beam or NoBeam()
        self.last_exec_time = -1
        self.g_val = GValuesBoscolo()
        self.LET = None

    @classmethod
    def from_config_dict(cls, config: ConfigDict) -> UnitCell:
        """Create a Simulation Environment from a the generic ConfigDict:

        Args:
            config: ConfigDict
                Dictionary (see radiopyo.sim_types.py) containing simulation parameters.

        Returns:
            UnitCell 
        """
        beam: Optional[ParticleBeam] = None

        if len(config["beam"]) > 0:
            logger.info("Beam configuration detected in config file")
            beam, LET = beam_from_dict(config["beam"])
            logger.info(f"Beam {beam} is setup")
        else:
            logger.info("No beam configuration detected")

        out = UnitCell(env=SimEnv.from_dict(config),
                       beam=beam)

        if LET is not None:
            out.change_LET(LET)
        return out

    @classmethod
    def from_parser(cls, file: str | Path, Parser: Type[FileParser]) -> UnitCell:
        """Create a Simulation Environment from a user defined FileParser.

        Args:
            file: str | pathlib.Path
                Path to the configuration file
            parser: Type[FileParser]
                USer defined parser to use to parse the provided file

        Returns:
            UnitCell
        """
        logger.info(f"Loading config from file: {file}")
        config = Parser(file).parse()
        logger.info("File is loaded")
        return UnitCell.from_config_dict(config)

    @classmethod
    def from_toml(cls, file: str | Path) -> UnitCell:
        """Create a Simulation Environment from a TOML file:

        Args:
            file: str | pathlib.Path
                Path to the TOML file

        Returns:
            UnitCell
        """
        logger.info(f"Loading config from TOML file: {file}")
        config = TOMLFileParser(file).parse()
        logger.info("File is loaded")
        return UnitCell.from_config_dict(config)

    @ classmethod
    def from_ron(cls, file: str | Path) -> UnitCell:
        """Create a Simulation Environment from a RON file:

        Args:
            file: str | pathlib.Path
                Path to the RON file

        Returns:
            UnitCell
        """
        logger.info(f"Loading config from RON file: {file}")
        config = RonFileParser(file).parse()
        logger.info("File is loaded")
        return UnitCell.from_config_dict(config)

    def set_beam(self, beam: ParticleBeam) -> None:
        """Define Radiation Source structure/parameters

        Args:
            beam: ParticleBeam

        Returns:
            None
        """
        self.beam = beam
        logger.info(f"Setting new beam definition: {beam}, done!")

    def change_LET(self, let: float) -> None:
        """Update the radiolytic yields from the simulation env according to literature
        values.

        Args:
            let: float
                LET in [kev/µm] and in range [0.14, 150]

        Returns:
            None
        """
        logger.info(
            f"Adapting (from Boscolo) Radiolytic yields to {let} keV/µm")
        for reaction in self.env.reactions.radiolytic:
            sp = reaction.species()
            if sp.label not in self.g_val or isinstance(sp, CstSpecies):
                continue
            reaction.update_from_ge(self.g_val.values_at(sp.label, let))

        # Recompute the Gi_mat accordingly
        self.mat["Gi_mat"] = self.env.make_radiolytic_vector()
        self.LET = let

        logger.debug("New G values [radical / 100eV / incident particle] are:")
        for reaction in self.env.reactions.radiolytic:
            logger.debug(f"-> {reaction.label()}: {reaction.ge():.3f}")

    def _raw_run(self,
                 t_span: Tuple[float, float],
                 use_jac: bool,
                 method: str,
                 atol: Optional[float] = None,
                 rtol: Optional[float] = None,
                 max_step: Optional[float] = None,
                 y0: Optional[NDArray[np.float64]] = None,
                 t_eval: Optional[NDArray[np.float64]] = None,
                 ) -> ScipyODEResult:

        if self.beam is NoBeam:
            rp_excep.NoBeamDefinedError("Simulation aborted")
        first_step = self.FIRST_STEP

        # Building kwargs for solve_ivp:
        kwargs: Dict[str, Any] = {}
        kwargs["first_step"] = self.FIRST_STEP
        kwargs["jac"] = jacobian if use_jac else None
        if rtol is not None:
            kwargs["rtol"] = rtol
        if atol is not None:
            kwargs["atol"] = atol

        kwargs["y0"] = y0
        if y0 is None:
            kwargs["y0"] = self.env.initial_values()

        if max_step is not None:
            kwargs["max_step"] = max_step
            # Check if first step is lower than the defined max_step
            if first_step >= max_step:
                kwargs["first_step"] = max_step / 10

        t0 = time.time()
        sol: ScipyODEResult = solve_ivp(derive,
                                        t_span,
                                        args=(self.beam,
                                              self.mat["cst_sp"],
                                              self.mat["AB_mat"],
                                              self.mat["Ki_mat"],
                                              self.mat["Gi_mat"],
                                              self.mat["Ez_mat"],
                                              self.mat["O2_mat"],
                                              ),
                                        t_eval=t_eval,
                                        method=method,
                                        **kwargs
                                        )
        t1 = time.time()
        self.last_exec_time = t1 - t0
        return sol

    def run(self,
            t_span: Tuple[float, float],
            t_eval: Optional[NDArray[np.float64]] = None,
            use_jac: bool = True,
            method: str = "LSODA",
            atol: Optional[float] = None,
            rtol: Optional[float] = None,
            max_step: Optional[float] = None,
            y0: Optional[NDArray[np.float64]] = None,
            ) -> ODEResult:
        """ Solve ODE over t_span
        Args:
            t_span: Tuple[float, float]
                t_span[0] -> start time, in seconds
                t_span[1] -> End time, in seconds
            t_eval: [NDArray[np.float64], optional
                Force function evaluation on t_eval (default: None)
            use_jac: bool
                Use the Jacobian ?
            method: str
                Which ODE solver to use? Define in scipy.integrate.solve_ivp
            atol: float, optional
                Absolute Tolerance of the ODE Solver (default: 1e-8)
            max_step: float, optional
                Larger authorized time step for the ODE solver (default: inf)
            y0: NDArray[np.float64], optional
                Initial concentration of DynSpecies (default: None). If None, calls 
                EnvSim.initial_values()

        Returns:
            ScipyODEResult: integration result provided by scipy.integrate.solve_ivp
        """
        logger.info(
            f"Starting simulation with t_span: {t_span} "
            f"(jacobian use: {use_jac})")
        sol = self._raw_run(t_span,
                            use_jac,
                            method,
                            atol,
                            rtol,
                            max_step,
                            y0,
                            t_eval,
                            )
        logger.info("Simulation done!")
        return ODEResult(self.env, self.beam, self.LET, sol)

    def run_chunked(self,
                    t_chunks: NDArray[np.float64],
                    t_eval: Optional[NDArray[np.float64]] = None,
                    use_jac: bool = True,
                    method: str = "LSODA",
                    atol: Optional[float] = None,
                    rtol: Optional[float] = None,
                    max_step: float | Iterable[float] = np.inf,
                    ) -> ODEResult:
        """TODO: Clean and remove"""
        sol = ODEResult(self.env, self.beam, self.LET)

        # Check up max_step
        if isinstance(max_step, float):
            max_step = repeat(max_step, len(t_chunks)-1)
        elif t_chunks.shape[0] != len(list(max_step)):
            raise ValueError(f"t_chunks (len={len(t_chunks)}) and max_step "
                             f"(len={len(list(max_step))}) must have the same "
                             f"length if max_step is a list or and array")

        # Checkup t_eval
        iter_t_eval: Iterable | List
        if t_eval is None:
            iter_t_eval = repeat(None, len(t_chunks)-1)
        else:
            t_eval = np.array(t_eval)
            iter_t_eval = []
            for idx, stop in enumerate(t_chunks[1:], 1):
                to_push = []
                with suppress(IndexError):
                    while t_eval[0] < stop:
                        to_push.append(t_eval[0])
                        # remove first element
                        t_eval = np.delete(t_eval, (0))

                # if nothing to add in this interval, at least at startup point
                if len(to_push) == 0:
                    to_push.append(t_chunks[idx-1])
                    to_push.append(t_chunks[idx-1]+(stop-t_chunks[idx-1])/2)
                iter_t_eval.append(np.array(to_push))

        # Loop through
        for idx, (stop, step_size, _t_eval) in enumerate(zip(t_chunks[1:],
                                                             max_step,
                                                             iter_t_eval,),
                                                         start=0):
            res = self._raw_run(y0=sol.y[:, -1] if len(sol) > 0 else None,
                                t_span=(t_chunks[idx], stop),
                                use_jac=use_jac,
                                method=method,
                                atol=atol,
                                rtol=rtol,
                                max_step=step_size,
                                t_eval=_t_eval,
                                )
            sol.append_sol(res, exec_time=self.last_exec_time)
        return sol

    def prepare_chunked_run(self,
                            t_span: Tuple[float, float],
                            use_jac: bool = True,
                            method: str = "LSODA",
                            atol: Optional[float] = None,
                            rtol: Optional[float] = None,
                            max_step_size_on: Optional[float] = None,
                            max_step_size_off: Optional[float] = None,
                            max_step_size_final: Optional[float] = None,
                            y0: Optional[NDArray[np.float64]] = None,
                            ) -> SimChunksIter:
        """ Method used to prepare a multi step ODE integration
        Args:
            t_span: Tuple[float, float]
                t_span[0] -> start time, in seconds
                t_span[1] -> End time, in seconds
            use_jac: bool, optional
                Whether to use the jacobian along with the ODE solver (default:False)
            method: str, optional
                Which ODE solver to use? Define in scipy.integrate.solve_ivp
            atol: float, optional
                Absolute Tolerance of the ODE Solver (default: 1e-10)
            max_step_size_on: float, optional
                Larger authorized time step for the ODE solver while beam is ON
                (default: inf)
            max_step_size_off: float, optional
                Larger authorized time step for the ODE solver while beam is OFF
                (default: inf)
            max_step_size_final: float, optional
                Larger authorized time step for the ODE solver for the very last
                integration step (can be the larger if long sim time and few pulses).
                If None, 'max_step_size_off|on' is used instead (default: None).
            y0: NDArray[np.float64], optional
                Initial concentration of DynSpecies (default: None). If None, calls 
                EnvSim.initial_values()

        Returns:
            SimChunksIter: Iterable object performing the ODE integration
        """
        if self.beam is NoBeam:
            logger.error("No beam is defined for the simulation")
            raise rp_excep.NoBeamDefinedError
        chunks = self.beam.make_time_intervals(start=t_span[0], stop=t_span[1])
        chunks[np.where(chunks[:, 1] > 0)[0], 1] = max_step_size_on
        chunks[np.where(chunks[:, 1] == 0)[0], 1] = max_step_size_off
        if max_step_size_final is not None:
            chunks[-2, 1] = max_step_size_final

        chunk_list: List[SimChunk] = []
        for idx in range(len(chunks[0:-1, :])):
            chunk = SimChunk(start=chunks[idx, 0],
                             stop=chunks[idx+1, 0],
                             max_step_size=chunks[idx, 1])
            chunk_list.append(chunk)
        logger.info(f"Simulation divided in {len(chunk_list)} chunks")
        return SimChunksIter(sim=self,
                             chunks=chunk_list,
                             atol=atol,
                             rtol=rtol,
                             use_jac=use_jac,
                             method=method,
                             y0=y0,)


class SimChunk(NamedTuple):
    """Convenient (immutable) container defining a chunk/ step of ODE integration.

    Attributes:
        start: float
            Start time, in seconds
        stop: float
            Stop time, in seconds
        max_step_size: float
            Max step size allowed over this integration step
    """
    start: float
    stop: float
    max_step_size: float


class SimChunksIter():
    """ Iterable class performing a multi-step integration.

    Attributes:
        sim: UniCell
            Simulation volume in which simulation takes place
        chucks: List[SimChunk]
            List of simulation steps/chunks
        sol: ODEResult
            Aggregate of ODE results from scipy.integrate.solve_ivp
        current: int
            Keeps tracks of the current integration step.
        atol: float
            Absolute tolerance used for scipy.integrate.solve_ivp
        use_jac: bool
            Whether jacobian is used with scipy.integrate.solve_ivp
        method: str
            Defines the ODE method used with scipy.integrate.solve_ivp
        y0: NDArray[np.float64], optional
            Initial concentration of DynSpecies (default: None). If None, calls 
            EnvSim.initial_values()
    Usage:
        >>> 
    """
    sim: UnitCell
    chunks: List[SimChunk]
    sol: ODEResult
    current: int
    atol: Optional[float]
    rtol: Optional[float]
    use_jac: bool
    method: str
    y0: Optional[NDArray[np.float64]]

    def __init__(self, sim: UnitCell,
                 chunks: List[SimChunk],
                 atol: Optional[float] = None,
                 rtol: Optional[float] = None,
                 use_jac: bool = True,
                 method: str = "LSODA",
                 y0: Optional[NDArray[np.float64]] = None,
                 ):
        self.sim = sim
        self.chunks = chunks
        self.sol = ODEResult(self.sim.env, self.sim.beam, self.sim.LET)
        self.atol = atol
        self.rtol = rtol
        self.use_jac = use_jac
        self.method = method
        self.current = -1
        self.y0 = y0

    def __iter__(self) -> SimChunksIter:
        return self

    def __next__(self) -> ODEResult:
        self.current += 1
        try:
            chunk = self.chunks[self.current]
        except IndexError:
            raise StopIteration from None
        res = self.sim._raw_run(
            t_span=(chunk.start, chunk.stop),
            use_jac=self.use_jac,
            method=self.method,
            atol=self.atol,
            rtol=self.rtol,
            max_step=chunk.max_step_size,
            y0=self.sol.y[:, -1] if self.current > 0 else self.y0,
        )
        self.sol.append_sol(res, exec_time=self.sim.last_exec_time)
        return self.sol

    def __len__(self) -> int:
        return len(self.chunks)

    def run(self) -> ODEResult:
        """Run all simulation chunks
        Tries to import tqdm to give a nice progressbar to the user. But it can safely 
        fail.
        get_ipython().__class__.__name__ --> detect if python was started from a 
        Jupyter Notebook (=="ZMQInteractiveShell").
        If ran from a regular python, the call to "get_ipython()" produces a NameError.
        """
        logger.info("Running all simulation chunks")
        try:
            from tqdm import tqdm  # type: ignore[import]
            with suppress(NameError):
                gip = get_ipython()  # type: ignore[name-defined] # noqa: F821
                if (gip.__class__.__name__ == "ZMQInteractiveShell"):
                    from tqdm.notebook import tqdm  # type: ignore[import] # noqa: F811
            # Perform Simulation with TQDM.
            sim = tqdm(self)
        except ImportError:
            sim = self

        # Efficient way to consume the iterator while returning the last element.
        res = collections.deque(sim, maxlen=1).pop()
        logger.info("Simulation done!")
        return res

    def subdivide_last_chunk(self,
                             n: int,
                             max_step_size: Optional[float] = None,
                             scale: str = "log",
                             ) -> SimChunksIter:
        last = self.chunks.pop()
        if max_step_size is None or max_step_size < last.max_step_size:
            max_step_size = last.max_step_size

        if scale == "linear":
            chunks = np.linspace(last.start, last.stop, num=n+1)
            steps = np.linspace(last.max_step_size, max_step_size, num=n)
        elif scale == "log":
            chunks = np.logspace(np.log10(last.start),
                                 np.log10(last.stop),
                                 num=n+1)
            steps = np.logspace(np.log10(last.max_step_size),
                                np.log10(max_step_size),
                                num=n)

        for idx in range(n):
            self.chunks.append(SimChunk(chunks[idx],
                                        chunks[idx+1],
                                        steps[idx],
                                        ))
        return self


class ODEResult(object):
    """
    Collection storing the results from scipy.integrate.solve_ivp in order to easily
    append together the results from several integration steps.
    """
    # Keep track of the original sim env
    env: SimEnv
    beam: ParticleBeam
    LET: Optional[float]

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

    def __init__(self,
                 env: SimEnv,
                 beam: ParticleBeam,
                 let: Optional[float],
                 sol: Optional[ScipyODEResult] = None,
                 ) -> None:
        self.env = copy.deepcopy(env)
        self.beam = copy.deepcopy(beam)
        self.LET = let
        self._default_init()
        if sol is not None:
            self.append_sol(sol)

    def __len__(self) -> int:
        """Number of simulation run"""
        return len(self.n_points)

    def __add__(self, other: Any) -> ODEResult:
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

    @classmethod
    def from_cell(cls,
                  cell: UnitCell,
                  sol: Optional[ScipyODEResult] = None,
                  ) -> ODEResult:
        return cls(cell.env, cell.beam, cell.LET, sol)

    @property
    def time(self) -> NDArray[np.float64]:
        return self.t

    @property
    def shape(self) -> Tuple:
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

        # (Re)Build used matrices
        mat = SimMatrices(
            AB_mat=self.env.make_acidBase_matrix(),
            Ki_mat=self.env.make_kreaction_matrix(),
            Gi_mat=self.env.make_radiolytic_vector(),
            Ez_mat=self.env.make_enzymatic_matrix(),
            cst_sp=self.env.make_cst_species_vector(),
            O2_mat=self.env.make_O2_intake_matrix(),
        )
        sp_label = [elt.label for elt in self.env.iter_DynSpecies()]
        dsp_label = ["".join(["d", elt]) for elt in sp_label]

        out = pd.DataFrame(index=self.t,
                           columns=pd.MultiIndex.from_product([sp_label, dsp_label],
                                                              names=["sp", "dsp"]),
                           dtype=np.float64,)

        for (idx, t_val) in enumerate(self.t):
            jac = jacobian(t_val, self.y[:, idx],
                           self.beam,
                           **mat,
                           )
            out.iloc[idx, :] = jac.flatten()
        return out

    def reconstruct_reactions(self) -> pd.DataFrame:
        """ 
        Rem: to integrate df
            >>> df.apply(lambda x: simpson(y=x.values, x=x.index))
        """
        # (Re)Build used matrices
        mat = SimMatrices(
            AB_mat=self.env.make_acidBase_matrix(),
            Ki_mat=self.env.make_kreaction_matrix(),
            Gi_mat=self.env.make_radiolytic_vector(),
            Ez_mat=self.env.make_enzymatic_matrix(),
            cst_sp=self.env.make_cst_species_vector(),
            O2_mat=self.env.make_O2_intake_matrix(),
        )

        # Compute column names and Pre compute slices to ease access in the dataframe
        # columns
        cols_kr = [kr.as_label() for kr in self.env.reactions.k_reaction]
        idx_kr = slice(0, len(cols_kr))

        cols_gi = [gr.as_label() for gr in self.env.reactions.radiolytic]
        idx_gi = slice(idx_kr.stop, idx_kr.stop+len(cols_gi))

        if mat["Ez_mat"] is not None:
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
            y = np.append(y, mat["cst_sp"])

            # K Reactions
            Rki = np.nan_to_num(
                np.einsum("j, jik -> ik ", y, mat["Ki_mat"][2]), nan=1)
            Rki = np.prod(Rki, axis=1) * mat["Ki_mat"][0]
            out.iloc[idx, idx_kr] = Rki

            # Radiolytic Reactions
            dr = 0.0
            if self.beam is not None:
                dr = self.beam.at(t_val).dose_rate()
            out.iloc[idx, idx_gi] *= dr

            # Enzimatic reactions
            if mat["Ez_mat"] is not None:
                E = np.matmul(y, mat["Ez_mat"][1])
                S = np.matmul(y, mat["Ez_mat"][0])
                MM = mat["Ez_mat"][2][1, :] * E * \
                    S / (mat["Ez_mat"][2][0, :] + S)
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
