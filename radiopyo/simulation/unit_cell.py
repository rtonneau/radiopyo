#!/usr/bin/env python3

"""
The sim_hand submodule contains class(es) to easily handle simulation.
This one should be the entry point for lambda user.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import collections
import logging
import time
import typing as tp
from contextlib import suppress

import numpy as np
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

from radiopyo.chemistry.utils.units import Temperature

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.parser.file_parser import FileParser, resolve_includes
from radiopyo.parser.ron import RonFileParser
from radiopyo.parser.toml import TOMLFileParser
from radiopyo.physics.beam import BeamCollection, ParticleBeam
from radiopyo.simulation import exceptions as rp_excep
from radiopyo.simulation.environment_build import SimEnv
from radiopyo.simulation.ode_functions_solver import (
    ScipyODEResult,
    derive,
    jacobian,
)
from radiopyo.simulation.ode_result import ODEResult
from radiopyo.utils.g_values import GValuesBoscolo

# ------------------------- Type checking ONLY imports------------------------ #
if tp.TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

    from radiopyo.utils.sim_types import ConfigDict, KwargsRun, SimMatrices

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "UnitCell",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class UnitCell(object):
    """
    Class UnitCell represents the tiniest simulation volume. All
    concentrations are assumed to be homogenous in it. It also defines
    method to solve ODE's in order to compute the evolution of
    concentration over time.

    Attributes:
        T: Temperature
            Sets the values of chemical reaction constants (default=25).
            Unit = °C.
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
        T_ZERO: float
            Initial time after when switching to BEAM ON phase
            (default=1e-7)

    """

    FIRST_STEP: float = 1e-15
    T_ZERO: float = 1e-7

    _T: Temperature
    env: SimEnv
    mat: SimMatrices
    beam_list: BeamCollection
    last_exec_time: float
    g_val: GValuesBoscolo

    def __init__(
        self,
        env: SimEnv,
        beam: tp.Optional[ParticleBeam | BeamCollection] = None,
        init: bool = True,
    ) -> None:
        self.env = env
        self._T = Temperature()
        self.g_val = GValuesBoscolo()  # To do before set_beam

        if init:
            self._update_matrices()

        if isinstance(beam, BeamCollection):
            self.beam_list = beam
            self.set_beam(self.beam_list.current_name)
        else:
            self.beam_list = BeamCollection()
            if isinstance(beam, ParticleBeam):
                self.beam_list.add_beam(beam, beam.label)

        self.last_exec_time = -1

    @property
    def T(self) -> float:
        return self._T.celsius

    @property
    def beam(self) -> ParticleBeam:
        return self.beam_list.current

    @property
    def LET(self) -> float:
        if self.beam_list.current_LET is None:
            raise ValueError("No LET defined... o_O ?")
        return self.beam_list.current_LET

    @classmethod
    def from_config_dict(
        cls,
        config: ConfigDict,
        init: bool = True,
    ) -> tp.Self:
        """
        Create a Simulation Environment from a the generic ConfigDict:

        Args:
            config: ConfigDict
                Dictionary (see radiopyo.sim_types.py) containing
                simulation parameters.

        Returns:
            UnitCell
        """
        beam_list: tp.Optional[BeamCollection] = None

        if len(config["beam"]) > 0:
            logger.info("Beam configuration detected in config file")
            beam_list = BeamCollection.from_config_dict(config["beam"])
        else:
            logger.info("No beam configuration detected")

        return cls(env=SimEnv.from_dict(config), beam=beam_list, init=init)

    @classmethod
    def from_parser(
        cls,
        file: str | Path,
        Parser: tp.Type[FileParser],
        init: bool = True,
    ) -> tp.Self:
        """
        Create a Simulation Environment from a user defined FileParser.

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
        return cls.from_config_dict(config, init)

    @classmethod
    def from_toml(cls, file: str | Path, init: bool = True) -> tp.Self:
        """Create a Simulation Environment from a TOML file:

        Args:
            file: str | pathlib.Path
                Path to the TOML file

        Returns:
            UnitCell
        """
        return cls.from_parser(file, TOMLFileParser, init)

    @classmethod
    def from_ron(cls, file: str | Path) -> tp.Self:
        """Create a Simulation Environment from a RON file:

        Args:
            file: str | pathlib.Path
                Path to the RON file

        Returns:
            UnitCell
        """
        return cls.from_parser(file, RonFileParser)

    def _update_matrices(self) -> None:
        """ """
        self.mat = self.env.build_all_matrices(self._T)

    def set_temperature(self, T: float, unit: str = "celsius") -> tp.Self:
        """ """
        self._T = Temperature(T, unit)
        self._update_matrices()
        return self

    def set_beam(
        self,
        beam: ParticleBeam | str,
        LET: tp.Optional[float] = None,
    ) -> tp.Self:
        """Define Radiation Source structure/parameters

        Args:
            beam: ParticleBeam

        Returns:
            None
        """
        LET = LET if LET is not None else self.beam_list.current_LET
        if isinstance(beam, str):
            self.beam_list.use_beam(beam)
            logger.info(
                f"Setting beam definition ({self.beam_list._current}): "
                f"{self.beam_list.current}, done!"
            )
        else:
            self.beam_list.add_beam(beam, beam.label, let=LET, use=True)
            logger.info(
                f"Setting new beam definition ({self.beam_list._current}): "
                f"{self.beam_list.current}, done!"
            )
        self.change_LET(LET)  # Update GValues
        return self

    def change_LET(self, let: float) -> None:
        """
        Update the radiolytic yields from the simulation env according
        to literature values.

        Args:
            let: float
                LET in [kev/µm] and in range [0.14, 150]

        Returns:
            None
        """
        self.beam_list.update_let(let)
        logger.info(
            f"Adapting (from Boscolo) Radiolytic yields to {let} keV/µm for {self.beam.dose_rate_peak:.2g} Gy/s"
        )
        for reaction in self.env.reactions.radiolytic:
            sp = reaction.species()
            if sp.label not in self.g_val or (not self.env.species.is_dyn_species(sp)):
                continue
            reaction.update_from_ge(
                self.g_val.values_at(sp.label, let, self.beam.dose_rate_peak)
            )

        # Recompute the Gi_mat accordingly
        # self.mat["Gi_mat"] = self.env.make_radiolytic_vector()
        self._update_matrices()

        logger.debug(
            f"New G values [radical / 100eV / incident particle] at {self.beam.dose_rate_peak:.2g} Gy/s are:"
        )
        for reaction in self.env.reactions.radiolytic:
            logger.debug(f"-> {reaction.label()}: {reaction.ge():.3f}")

    def _raw_run(
        self,
        t_span: tp.Tuple[float, float],
        use_jac: bool,
        method: str,
        atol: tp.Optional[float] = None,
        rtol: tp.Optional[float] = None,
        max_step: tp.Optional[float] = None,
        y0: tp.Optional[NDArray[np.float64]] = None,
        t_eval: tp.Optional[NDArray[np.float64]] = None,
        use_beam: bool = True,
        beam_ON: bool = True,
        force_positive: bool = False,
    ) -> ScipyODEResult:
        if len(self.beam_list) <= 0:
            rp_excep.NoBeamDefinedError("Simulation aborted")
        first_step = self.FIRST_STEP

        if t_span[0] == t_span[1]:
            t_span = (t_span[0] / 1e2, t_span[1])

        # Building kwargs for solve_ivp:
        kwargs: tp.Dict[str, tp.Any] = {}
        kwargs["first_step"] = self.FIRST_STEP
        kwargs["jac"] = jacobian if use_jac else None
        kwargs["rtol"] = rtol if rtol is not None else 1e-10
        kwargs["atol"] = atol if atol is not None else 1e-12

        kwargs["y0"] = y0
        if y0 is None:
            kwargs["y0"] = self.env.initial_values()

        if max_step is not None:
            kwargs["max_step"] = max_step
            # Check if first step is lower than the defined max_step
            if first_step >= max_step:
                kwargs["first_step"] = max_step / 10

        if kwargs["first_step"] > t_span[1] - t_span[0]:
            kwargs["first_step"] = (t_span[1] - t_span[0]) / 1e3

        beam = self.beam_list.current if use_beam else None

        # Check if beam_on, then skip the next 1e-7 s of sim and compute
        # evolution of cc's according to Gvalues
        if beam_ON:  # assume beam is switched ON at t_span[0]
            t_span = (t_span[0] + self.T_ZERO, t_span[1])
            if t_span[0] > t_span[1]:
                raise ValueError(f"Too Short beam pulse {t_span}")
            for r in self.env.reactions.radiolytic:
                # If the species must be forced constant, skip it.
                if self.env.species.is_cst_species(r.species()):
                    continue

                # Keep initial value and increment
                kwargs["y0"][r.species().index] += (
                    r.compute_initial_cc_before(
                        t0=self.T_ZERO, dose_rate=self.beam.at(t_span[0]).dose_rate()
                    )
                    * 1e6
                )  # Units of initial cc are µM

        t0 = time.time()
        sol: ScipyODEResult = solve_ivp(
            derive,
            t_span,
            args=(
                beam,
                self.mat["cst_sp"],
                self.mat["Ki_mat"],
                self.mat["Gi_mat"],
                self.mat["Ez_mat"],
                force_positive,
            ),
            t_eval=t_eval,
            method=method,
            **kwargs,
        )
        t1 = time.time()
        self.last_exec_time = t1 - t0
        return sol

    def run_till_max_dose(
        self,
        time_start: float = 1e-9,
        **kwargs: tp.Unpack[KwargsRun],
    ) -> ODEResult:
        """ """
        if self.beam_list.current.time_at_max_dose >= np.inf:
            raise ValueError(
                "Unable to use 'run_till_max_dose', time at max dose "
                f"is inf (-> {self.beam_list.current.time_at_max_dose})"
            )
        kwargs["t_span"] = (time_start, self.beam_list.current.time_at_max_dose)
        return self.run(**kwargs)

    def run(
        self,
        t_span: tp.Tuple[float, float],
        t_eval: tp.Optional[NDArray[np.float64]] = None,
        use_jac: bool = False,
        method: str = "LSODA",
        atol: tp.Optional[float] = None,
        rtol: tp.Optional[float] = None,
        max_step: tp.Optional[float] = None,
        y0: tp.Optional[NDArray[np.float64]] = None,
        use_beam: bool = True,
        force_positive: bool = False,
    ) -> ODEResult:
        """Solve ODE over t_span
        Args:
            t_span: Tuple[float, float]
                t_span[0] -> start time, in seconds
                t_span[1] -> End time, in seconds
            t_eval: [NDArray[np.float64], optional
                Force function evaluation on t_eval (default: None)
            use_jac: bool
                Use the Jacobian ?
            method: str
                Which ODE solver to use? Define in
                scipy.integrate.solve_ivp
            atol: float, optional
                Absolute Tolerance of the ODE Solver (default: 1e-8)
            max_step: float, optional
                Larger authorized time step for the ODE solver
                (default: inf)
            y0: NDArray[np.float64], optional
                Initial concentration of DynSpecies (default: None).
                If None, calls EnvSim.initial_values()
            use_beam: bool, optional
                whether to use the beam or not (default: True).

        Returns:
            ScipyODEResult: integration result provided by
            scipy.integrate.solve_ivp
        """
        logger.info(
            f"Starting simulation with t_span: {t_span} (jacobian use: {use_jac})"
        )
        sol = self._raw_run(
            t_span,
            use_jac,
            method,
            atol,
            rtol,
            max_step,
            y0,
            t_eval,
            use_beam,
            force_positive,
        )
        logger.info(f"Simulation done! ({self.last_exec_time:.2g} s)")
        return ODEResult(
            self.env,
            self.beam_list.current,
            self.mat,
            self.beam_list.current_LET,
            self.T_ZERO,
            sol,
            force_positive=force_positive,
        )

    def prepare_chunked_run(
        self,
        t_span: tp.Tuple[float, float],
        use_jac: bool = False,
        method: str = "LSODA",
        atol: tp.Optional[float] = None,
        rtol: tp.Optional[float] = None,
        max_step_size_on: tp.Optional[float] = None,
        max_step_size_off: tp.Optional[float] = None,
        max_step_size_final: tp.Optional[float] = None,
        y0: tp.Optional[NDArray[np.float64]] = None,
        use_beam: bool = True,
        force_positive: bool = False,
    ) -> SimChunksIter:
        """Method used to prepare a multi step ODE integration
        Args:
            t_span: Tuple[float, float]
                t_span[0] -> start time, in seconds
                t_span[1] -> End time, in seconds
            use_jac: bool, optional
                Whether to use the jacobian along with the ODE solver
                (default:False)
            method: str, optional
                Which ODE solver to use? Define in
                scipy.integrate.solve_ivp
            atol: float, optional
                Absolute Tolerance of the ODE Solver (default: 1e-10)
            max_step_size_on: float, optional
                Larger authorized time step for the ODE solver while
                beam is ON (default: inf)
            max_step_size_off: float, optional
                Larger authorized time step for the ODE solver while
                beam is OFF (default: inf)
            max_step_size_final: float, optional
                Larger authorized time step for the ODE solver for the
                very last integration step (can be the larger if long
                sim time and few pulses). If None,
                'max_step_size_off|on' is used instead (default: None).
            y0: NDArray[np.float64], optional
                Initial concentration of DynSpecies (default: None).
                If None, calls EnvSim.initial_values()
            use_beam: bool, optional
                whether to use the beam or not (default: True).
            force_positive: bool, optional
                Whether to force the concentration to be positive
                (default: False).

        Returns:
            SimChunksIter: Iterable object performing the ODE
            integration.
        """
        if len(self.beam_list) <= 0:
            logger.error("No beam is defined for the simulation")
            raise rp_excep.NoBeamDefinedError
        chunks = self.beam_list.current.make_time_intervals(
            start=t_span[0], stop=t_span[1]
        )
        # chunks[np.where(chunks[:, 1] > 0)[0], 1] = max_step_size_on
        # chunks[np.where(chunks[:, 1] == 0)[0], 1] = max_step_size_off
        # if max_step_size_final is not None:
        #    chunks[-2, 1] = max_step_size_final

        chunk_list: tp.List[SimChunk] = []
        for idx in range(len(chunks[0:-2, :])):
            chunk = SimChunk(
                start=chunks[idx, 0],
                stop=chunks[idx + 1, 0],
                max_step_size=(
                    max_step_size_on if chunks[idx, 1] > 0 else max_step_size_off
                ),
                beam_on=chunks[idx, 1] > 0,
            )
            chunk_list.append(chunk)

        # Make final chunk
        max_step = max_step_size_final
        if max_step is None:
            max_step = max_step_size_on if chunks[-2, 1] > 0 else max_step_size_off
        chunk_list.append(
            SimChunk(
                start=chunks[-2, 0],
                stop=chunks[-1, 0],
                max_step_size=max_step,
                beam_on=chunks[-2, 1] > 0,
            )
        )

        logger.info(f"Beam span ==> [{t_span[0]} - {t_span[1]}]")
        logger.info(f"From beam chunks ==> {len(chunks)}")
        logger.info(f"Simulation divided in {len(chunk_list)} chunks")
        return SimChunksIter(
            sim=self,
            chunks=chunk_list,
            atol=atol,
            rtol=rtol,
            use_jac=use_jac,
            method=method,
            y0=y0,
            use_beam=use_beam,
            force_positive=force_positive,
        )

    def prepare_single_pulse(
        self,
        t_span: tp.Tuple[float, float],
        use_jac: bool = False,
        method: str = "LSODA",
        atol: tp.Optional[float] = None,
        rtol: tp.Optional[float] = None,
        max_step_size_on: tp.Optional[float] = None,
        max_step_size_trans: tp.Optional[float] = None,
        max_step_size_final: tp.Optional[float] = None,
        trans_time: tp.Optional[float] = None,
        y0: tp.Optional[NDArray[np.float64]] = None,
        use_beam: bool = True,
        force_positive: bool = False,
    ) -> SimChunksIter:
        """Method used to prepare a multi step ODE integration
        Args:
            t_span: Tuple[float, float]
                t_span[0] -> start time, in seconds
                t_span[1] -> End time, in seconds
            use_jac: bool, optional
                Whether to use the jacobian along with the ODE solver
                (default:False)
            method: str, optional
                Which ODE solver to use? Define in
                scipy.integrate.solve_ivp
            atol: float, optional
                Absolute Tolerance of the ODE Solver (default: 1e-10)
            max_step_size_on: float, optional
                Larger authorized time step for the ODE solver while
                beam is ON (default: inf)
            max_step_size_off: float, optional
                Larger authorized time step for the ODE solver while
                beam is OFF (default: inf)
            max_step_size_final: float, optional
                Larger authorized time step for the ODE solver for the
                very last integration step (can be the larger if long
                sim time and few pulses). If None,
                'max_step_size_off|on' is used instead (default: None).
            y0: NDArray[np.float64], optional
                Initial concentration of DynSpecies (default: None).
                If None, calls EnvSim.initial_values()
            use_beam: bool, optional
                whether to use the beam or not (default: True).
            force_positive: bool, optional
                Whether to force the concentration to be positive
                (default: False).

        Returns:
            SimChunksIter: Iterable object performing the ODE
            integration.
        """
        if len(self.beam_list) <= 0:
            logger.error("No beam is defined for the simulation")
            raise rp_excep.NoBeamDefinedError

        chunk_list: tp.List[SimChunk] = []

        trans_time = trans_time if trans_time is not None else 100e-6

        # Beam ON
        chunk_list.append(
            SimChunk(
                start=t_span[0],
                stop=self.beam_list.current.time_at_max_dose,
                max_step_size=max_step_size_on,
                beam_on=True,
            )
        )

        # Beam ON -> Beam OFF transition
        chunk_list.append(
            SimChunk(
                start=self.beam_list.current.time_at_max_dose,
                stop=self.beam_list.current.time_at_max_dose + trans_time,
                max_step_size=max_step_size_trans,
                beam_on=False,
            )
        )

        # Make final chunk
        chunk_list.append(
            SimChunk(
                start=self.beam_list.current.time_at_max_dose + trans_time,
                stop=t_span[1],
                max_step_size=max_step_size_final,
                beam_on=False,
            )
        )

        logger.info(f"Beam span ==> [{t_span[0]} - {t_span[1]}]")
        logger.info(f"Simulation divided in {len(chunk_list)} chunks")
        return SimChunksIter(
            sim=self,
            chunks=chunk_list,
            atol=atol,
            rtol=rtol,
            use_jac=use_jac,
            method=method,
            y0=y0,
            use_beam=use_beam,
            force_positive=force_positive,
        )


# -------------------------- Chunked Simulation Run ------------------------- #


class SimChunk(tp.NamedTuple):
    """
    Convenient (immutable) container defining a chunk/ step of
    ODE integration.

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
    max_step_size: float | None
    beam_on: bool


class SimChunksIter:
    """Iterable class performing a multi-step integration.

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
            Initial concentration of DynSpecies (default: None).
            If None, calls EnvSim.initial_values()
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
    use_beam: bool
    force_positive: bool

    def __init__(
        self,
        sim: UnitCell,
        chunks: tp.List[SimChunk],
        atol: tp.Optional[float] = None,
        rtol: tp.Optional[float] = None,
        use_jac: bool = False,
        method: str = "LSODA",
        y0: tp.Optional[NDArray[np.float64]] = None,
        use_beam: bool = True,
        force_positive: bool = False,
    ):
        self.sim = sim
        self.chunks = chunks
        self.sol = ODEResult(
            self.sim.env,
            self.sim.beam_list.current,
            self.sim.mat,
            self.sim.beam_list.current_LET,
            self.sim.T_ZERO,
            force_positive=force_positive,
        )
        self.atol = atol
        self.rtol = rtol
        self.use_jac = use_jac
        self.method = method
        self.current = -1
        self.y0 = y0
        self.use_beam = use_beam
        self.force_positive = force_positive

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
            use_beam=self.use_beam,
            beam_ON=chunk.beam_on,
            force_positive=self.force_positive,
        )
        self.sol.append_sol(res, exec_time=self.sim.last_exec_time)
        return self.sol

    def __len__(self) -> int:
        return len(self.chunks)

    def run(self) -> ODEResult:
        """Run all simulation chunks
        Tries to import tqdm to give a nice progressbar to the user.
        But it can safely fail.
        get_ipython().__class__.__name__ --> detect if python was
        started from a Jupyter Notebook (=="ZMQInteractiveShell").
        If ran from a regular python, the call to "get_ipython()"
        produces a NameError.
        """
        logger.info("Running all simulation chunks")
        try:
            from tqdm import tqdm  # type: ignore[import-untyped]

            with suppress(NameError):
                gip = get_ipython()  # type: ignore[name-defined] # noqa: F821
                if gip.__class__.__name__ == "ZMQInteractiveShell":
                    from tqdm.notebook import (  # type: ignore[import-untyped]
                        tqdm,
                    )
            # Perform Simulation with TQDM.
            sim = tqdm(self)
        except ImportError:
            sim = self

        # Efficient way to consume the iterator while returning the
        # last element.
        res = collections.deque(sim, maxlen=1).pop()
        logger.info("Simulation done!")
        return res

    # def subdivide_last_chunk(
    #     self,
    #     n: int,
    #     max_step_size: tp.Optional[float] = None,
    #     scale: str = "log",
    # ) -> SimChunksIter:
    #     last = self.chunks.pop()
    #     if max_step_size is None or max_step_size < last.max_step_size:
    #         max_step_size = last.max_step_size

    #     if scale == "linear":
    #         chunks = np.linspace(last.start, last.stop, num=n + 1)
    #         steps = np.linspace(last.max_step_size, max_step_size, num=n)
    #     elif scale == "log":
    #         chunks = np.logspace(np.log10(last.start), np.log10(last.stop), num=n + 1)
    #         steps = np.logspace(
    #             np.log10(last.max_step_size), np.log10(max_step_size), num=n
    #         )

    #     for idx in range(n):
    #         self.chunks.append(
    #             SimChunk(
    #                 chunks[idx],
    #                 chunks[idx + 1],
    #                 steps[idx],
    #             )
    #         )
    #     return self
