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
import logging
import time
import typing as tp
from contextlib import suppress

import numpy as np
from scipy.integrate import solve_ivp  # type: ignore[import]

# --------------------------------- Internal imports --------------------------------- #
from radiopyo.chemistry.species import CstSpecies
from radiopyo.parser.file_parser import FileParser, resolve_includes
from radiopyo.parser.ron import RonFileParser
from radiopyo.parser.toml import TOMLFileParser
from radiopyo.physics.beam import BeamCollection, ParticleBeam
from radiopyo.simulation import exceptions as rp_excep
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.ode_functions_solver import ScipyODEResult, derive, jacobian
from radiopyo.simulation.ode_result import ODEResult
from radiopyo.utils.g_values import GValuesBoscolo
from radiopyo.utils.sim_types import SimMatrices

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from radiopyo.utils.sim_types import ConfigDict, KwargsRun

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = [
    "UnitCell",
]

# -------------------------------------- Logging ------------------------------------- #
logger = logging.getLogger("radiopyo")

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


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
    beam_list: BeamCollection
    last_exec_time: float
    g_val: GValuesBoscolo

    def __init__(self,
                 env: SimEnv,
                 beam: tp.Optional[ParticleBeam | BeamCollection] = None,
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

        self.g_val = GValuesBoscolo()  # To do before set_beam

        if isinstance(beam, BeamCollection):
            self.beam_list = beam
            self.set_beam(self.beam_list.current_name)
        else:
            self.beam_list = BeamCollection()
            if isinstance(beam, ParticleBeam):
                self.beam_list.add_beam(beam, beam.label)

        self.last_exec_time = -1

    @property
    def beam(self) -> ParticleBeam:
        return self.beam_list.current

    @classmethod
    def from_config_dict(cls, config: ConfigDict) -> UnitCell:
        """Create a Simulation Environment from a the generic ConfigDict:

        Args:
            config: ConfigDict
                Dictionary (see radiopyo.sim_types.py) containing simulation parameters.

        Returns:
            UnitCell 
        """
        beam_list: tp.Optional[BeamCollection] = None

        if len(config["beam"]) > 0:
            logger.info("Beam configuration detected in config file")
            beam_list = BeamCollection.from_config_dict(config["beam"])
        else:
            logger.info("No beam configuration detected")

        return UnitCell(env=SimEnv.from_dict(config),
                        beam=beam_list)

    @classmethod
    def from_parser(cls, file: str | Path, Parser: tp.Type[FileParser]) -> UnitCell:
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
        config = resolve_includes(Parser(file).parse())
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
        return cls.from_parser(file, TOMLFileParser)

    @ classmethod
    def from_ron(cls, file: str | Path) -> UnitCell:
        """Create a Simulation Environment from a RON file:

        Args:
            file: str | pathlib.Path
                Path to the RON file

        Returns:
            UnitCell
        """
        return cls.from_parser(file, RonFileParser)

    def use_acid_base(self, enable: bool = True) -> None:
        """ """
        if enable:
            self.mat["AB_mat"] = self.env.make_acidBase_matrix()
        else:
            ABi = np.zeros([len(self.env.species), len(self.env.reactions.acid_base)],
                           dtype=np.int8)
            ABy = np.ones(len(self.env.species), dtype=np.int8)
            self.mat["AB_mat"] = (ABy, ABi, self.mat["AB_mat"][2])

    def set_beam(self,
                 beam: ParticleBeam | str,
                 LET: tp.Optional[float] = None,
                 ) -> None:
        """Define Radiation Source structure/parameters

        Args:
            beam: ParticleBeam

        Returns:
            None
        """
        if LET is None:
            LET = self.beam_list.current_LET

        if isinstance(beam, str):
            self.beam_list.use_beam(beam)
        else:
            self.beam_list.add_beam(beam, beam.label, let=LET, use=True)
        logger.info(f"Setting new beam definition ({self.beam_list._current}): "
                    f"{self.beam_list.current}, done!")
        with suppress(ValueError):
            self.change_LET(LET)

    def change_LET(self, let: float) -> None:
        """Update the radiolytic yields from the simulation env according to literature
        values.

        Args:
            let: float
                LET in [kev/µm] and in range [0.14, 150]

        Returns:
            None
        """
        self.beam_list.update_let(let)
        logger.info(
            f"Adapting (from Boscolo) Radiolytic yields to {let} keV/µm")
        for reaction in self.env.reactions.radiolytic:
            sp = reaction.species()
            if sp.label not in self.g_val or isinstance(sp, CstSpecies):
                continue
            reaction.update_from_ge(self.g_val.values_at(sp.label, let))

        # Recompute the Gi_mat accordingly
        self.mat["Gi_mat"] = self.env.make_radiolytic_vector()

        logger.debug("New G values [radical / 100eV / incident particle] are:")
        for reaction in self.env.reactions.radiolytic:
            logger.debug(f"-> {reaction.label()}: {reaction.ge():.3f}")

    def _raw_run(self,
                 t_span: tp.Tuple[float, float],
                 use_jac: bool,
                 method: str,
                 atol: tp.Optional[float] = None,
                 rtol: tp.Optional[float] = None,
                 max_step: tp.Optional[float] = None,
                 y0: tp.Optional[NDArray[np.float64]] = None,
                 t_eval: tp.Optional[NDArray[np.float64]] = None,
                 ) -> ScipyODEResult:

        if len(self.beam_list) <= 0:
            rp_excep.NoBeamDefinedError("Simulation aborted")
        first_step = self.FIRST_STEP

        # Building kwargs for solve_ivp:
        kwargs: tp.Dict[str, tp.Any] = {}
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
                                        args=(self.beam_list.current,
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

    def run_till_max_dose(
            self,
            time_start: float = 1e-9,
            **kwargs: tp.Unpack[KwargsRun],  # type: ignore [misc]
    ) -> ODEResult:
        """ """
        if self.beam_list.current.time_at_max_dose >= np.inf:
            raise ValueError("Unable to use 'run_till_max_dose', time at max dose "
                             f"is inf (-> {self.beam_list.current.time_at_max_dose})")
        return self.run(t_span=(time_start, self.beam_list.current.time_at_max_dose),
                        **kwargs)

    def run(self,
            t_span: tp.Tuple[float, float],
            t_eval: tp.Optional[NDArray[np.float64]] = None,
            use_jac: bool = False,
            method: str = "LSODA",
            atol: tp.Optional[float] = None,
            rtol: tp.Optional[float] = None,
            max_step: tp.Optional[float] = None,
            y0: tp.Optional[NDArray[np.float64]] = None,
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
        logger.info(f"Simulation done! ({self.last_exec_time:.2g} s)")
        return ODEResult(self.env,
                         self.beam_list.current,
                         self.mat,
                         self.beam_list.current_LET,
                         sol)

    def prepare_chunked_run(self,
                            t_span: tp.Tuple[float, float],
                            use_jac: bool = False,
                            method: str = "LSODA",
                            atol: tp.Optional[float] = None,
                            rtol: tp.Optional[float] = None,
                            max_step_size_on: tp.Optional[float] = None,
                            max_step_size_off: tp.Optional[float] = None,
                            max_step_size_final: tp.Optional[float] = None,
                            y0: tp.Optional[NDArray[np.float64]] = None,
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
        if len(self.beam_list) <= 0:
            logger.error("No beam is defined for the simulation")
            raise rp_excep.NoBeamDefinedError
        chunks = self.beam_list.current.make_time_intervals(start=t_span[0],
                                                            stop=t_span[1])
        chunks[np.where(chunks[:, 1] > 0)[0], 1] = max_step_size_on
        chunks[np.where(chunks[:, 1] == 0)[0], 1] = max_step_size_off
        if max_step_size_final is not None:
            chunks[-2, 1] = max_step_size_final

        chunk_list: tp.List[SimChunk] = []
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


# ------------------------------ Chunked Simulation Run ------------------------------ #


class SimChunk(tp.NamedTuple):
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
    chunks: tp.List[SimChunk]
    sol: ODEResult
    current: int
    atol: tp.Optional[float]
    rtol: tp.Optional[float]
    use_jac: bool
    method: str
    y0: tp.Optional[NDArray[np.float64]]

    def __init__(self, sim: UnitCell,
                 chunks: tp.List[SimChunk],
                 atol: tp.Optional[float] = None,
                 rtol: tp.Optional[float] = None,
                 use_jac: bool = False,
                 method: str = "LSODA",
                 y0: tp.Optional[NDArray[np.float64]] = None,
                 ):
        self.sim = sim
        self.chunks = chunks
        self.sol = ODEResult(self.sim.env,
                             self.sim.beam_list.current,
                             self.sim.mat,
                             self.sim.beam_list.current_LET)
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
                             max_step_size: tp.Optional[float] = None,
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
