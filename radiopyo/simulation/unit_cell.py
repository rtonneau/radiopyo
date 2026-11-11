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
import dataclasses
import logging
import time
import typing as tp
from contextlib import suppress

import numpy as np
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.chemistry.utils.units import Temperature
from radiopyo.parser.file_parser import FileParser, resolve_includes
from radiopyo.parser.toml import TOMLFileParser
from radiopyo.physics.beam import ConstantDoseRate, ParticleBeam
from radiopyo.physics.beam_collection import BeamCollection
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

    from radiopyo.physics.beam import DoseRateSource
    from radiopyo.utils.sim_types import ConfigDict, KwargsRun, SimMatrices

# ------------------------- Expose importable stuffs ------------------------ #
__all__ = [
    "RunOptions",
    "UnitCell",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                            FUNCTION DEFINITIONS                             #
# --------------------------------------------------------------------------- #

# Tolerances passed to solve_ivp when the user gives none.
DEFAULT_ATOL: float = 1e-12
DEFAULT_RTOL: float = 1e-10
# solve_ivp methods that accept a Jacobian (the others would warn)
JAC_METHODS: tp.Final = ("LSODA", "BDF", "Radau")
# Duration [s] of the ON -> OFF transition chunk of prepare_single_pulse.
DEFAULT_TRANS_TIME: float = 100e-6


def widen_empty_span(t_span: tp.Tuple[float, float]) -> tp.Tuple[float, float]:
    """
    solve_ivp needs a non empty integration interval. An empty span
    [t, t] is therefore widened to [t/100, t]: the solver starts two
    decades earlier and still ends at the requested time t. Any other
    span is returned unchanged.

    The trick only works for t > 0: [0, 0] stays empty.
    """
    if t_span[0] == t_span[1]:
        return (t_span[0] / 1e2, t_span[1])
    return t_span


def chunk_dose_rate(beam: ParticleBeam, start: float, stop: float) -> float:
    """
    Dose rate [Gy/s] of the beam inside the chunk [start, stop], read at
    its middle. The chunks are cut at the beam transitions, so the dose
    rate is constant over a chunk. An unbounded chunk is read at its start.
    """
    time_in_chunk = 0.5 * (start + stop) if np.isfinite(stop) else start
    return float(beam.at(time_in_chunk).dose_rate())


def _make_chunks(
    boundaries: tp.Sequence[float] | NDArray[np.float64],
    beam_on: tp.Sequence[bool] | NDArray[np.bool_],
    max_steps: tp.Sequence[tp.Optional[float]],
    beam: tp.Optional[ParticleBeam] = None,
) -> tp.List[SimChunk]:
    """
    Build consecutive chunks: chunk i spans [boundaries[i],
    boundaries[i+1]], with beam state beam_on[i] and max step size
    max_steps[i] (len(boundaries) == len(beam_on) + 1 == len(max_steps) + 1).
    With a beam, every chunk holds the (constant) dose rate of the beam
    inside it (see chunk_dose_rate).
    """
    return [
        SimChunk(
            start=start,
            stop=stop,
            max_step_size=max_step,
            beam_on=on,
            dose_rate=(
                chunk_dose_rate(beam, start, stop) if beam is not None else None
            ),
        )
        for start, stop, on, max_step in zip(
            boundaries[:-1], boundaries[1:], beam_on, max_steps, strict=True
        )
    ]


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


@dataclasses.dataclass(frozen=True, kw_only=True)
class RunOptions:
    """
    Solver options shared by every ODE integration of a UnitCell (plain
    run and every chunk of a chunked run). Built internally from the
    keyword arguments of UnitCell.run / prepare_chunked_run /
    prepare_single_pulse: users never need to create one.

    Attributes:
        atol: float, optional
            Absolute tolerance of solve_ivp (None -> DEFAULT_ATOL = 1e-12)
        rtol: float, optional
            Relative tolerance of solve_ivp (None -> DEFAULT_RTOL = 1e-10)
        method: str
            solve_ivp integration method (default: "LSODA")
        max_step: float, optional
            Largest step allowed to solve_ivp (None -> inf). In a
            chunked run, it is replaced by the max step of each chunk.
        first_step: float, optional
            Initial step of solve_ivp (None -> UnitCell.FIRST_STEP,
            read on the cell at run time). It is reduced to max_step / 10
            when it is not lower than max_step, then to span / 1000 when
            it is larger than the integration span.
        use_beam: bool
            Whether the beam dose rate enters the ODE (default: True)
        use_jac: bool
            Pass the analytical jacobian() to the solver; used only by
            LSODA, BDF and Radau, ignored by other methods (default: True)
        t_eval: NDArray[np.float64], optional
            Times at which solve_ivp stores the solution (None -> solver
            steps)
    """

    atol: tp.Optional[float] = None
    rtol: tp.Optional[float] = None
    method: str = "LSODA"
    max_step: tp.Optional[float] = None
    first_step: tp.Optional[float] = None
    use_beam: bool = True
    use_jac: bool = True
    t_eval: tp.Optional[NDArray[np.float64]] = None


class UnitCell(object):
    """
    Class UnitCell represents the tiniest simulation volume. All
    concentrations are assumed to be homogenous in it. It also defines
    method to solve ODE's in order to compute the evolution of
    concentration over time.

    Attributes:
        T: float
            Simulation temperature in °C (default=25, or the
            "temperature" parameter of the config file). Read-only view
            of env.parameters["temperature"], the single source of truth;
            set it with set_temperature (or the T setter). It sets the
            rate constants, the pKa and the concentrations of the
            constant species (e.g. H2O).
        env: SimEnv
            SimEnv build from the input file
        mat: SimMatrices
            Collections of matrices used to solve the ODE's
        beam: ParticleBeam
            Radiation source to use for the simulation
        last_exec_time: float
            Time duration of last simulation step
        last_burst_dose_rate: float, optional
            Dose rate [Gy/s] used for the initial radiolytic burst of the
            last integration step (None if no burst was applied)
        g_val: GValuesBoscolo
            GValues from the literature to possibly update the LET
        FIRST_STEP: float
            Initial time step of all simulation (default=1e-15 s)
        T_ZERO: float
            Duration skipped at the start of a BEAM ON phase; the
            radiolytic production over it is added analytically to the
            initial concentrations (default=1e-7 s)

    """

    FIRST_STEP: float = 1e-15
    T_ZERO: float = 1e-7

    env: SimEnv
    mat: SimMatrices
    beam_list: BeamCollection
    last_exec_time: float
    last_burst_dose_rate: tp.Optional[float]
    g_val: GValuesBoscolo

    def __init__(
        self,
        env: SimEnv,
        beam: tp.Optional[ParticleBeam | BeamCollection] = None,
        init: bool = True,
    ) -> None:
        self.env = env
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
        self.last_burst_dose_rate = None

    @property
    def T(self) -> float:
        """Simulation temperature in °C (see set_temperature)."""
        return self.env.temperature.celsius

    @T.setter
    def T(self, value: float) -> None:
        self.set_temperature(value)

    @property
    def _T(self) -> Temperature:
        return self.env.temperature

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

    def update_env(self) -> None:
        """ """
        self._update_matrices()

    def _update_matrices(self) -> None:
        """ """
        self.mat = self.env.build_all_matrices(self.env.temperature)

    def set_temperature(self, T: float, unit: str = "celsius") -> tp.Self:
        """
        Set the simulation temperature (stored in
        env.parameters["temperature"], in °C) and rebuild the matrices.

        Units: T is given in `unit` ("celsius" default, or "kelvin").
        Internally the k polynomials (Arrhenius/equilibria) use 1/K,
        while the pKa and the initial concentrations (polynomials, e.g.
        H2O) use °C: the conversion is done for you.

        Args:
            T: float
                Temperature
            unit: str
                "celsius" (default) or "kelvin"

        Returns:
            Self
        """
        self.env.set_temperature(T, unit)
        self._update_matrices()
        return self

    def set_beam(
        self,
        beam: ParticleBeam | str,
        LET: tp.Optional[float] = None,
    ) -> tp.Self:
        """Define Radiation Source structure/parameters, then update the
        G values for the LET (see change_LET).

        Args:
            beam: ParticleBeam | str
                New beam (added to the beam list and selected) or label
                of a beam already in the beam list (selected).
            LET: float, optional
                LET [keV/µm] of the beam. If None, the LET of the beam
                selected *before* this call is silently reused (current
                behaviour, kept on purpose).

        Returns:
            Self
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
        if self.env.parameters["variable_GValues"]:
            _dr = self.beam.dose_rate_peak
            logger.info(
                f"Adapting (from Boscolo) Radiolytic yields to {let} keV/µm for {_dr:.2g} Gy/s"
            )
        else:
            _dr = None
            logger.info(
                f"Adapting (from Boscolo) Radiolytic yields to {let} keV/µm (No dose rate dep.)"
            )

        for reaction in self.env.reactions.radiolytic:
            sp = reaction.species()
            if sp.label not in self.g_val or (not self.env.species.is_dyn_species(sp)):
                continue
            reaction.update_from_ge(self.g_val.values_at(sp.label, let, _dr))

        # Recompute the radiolytic yields (mat.radiolytic) accordingly
        self._update_matrices()

        logger.debug("New G values [radical / 100eV / incident particle] are:")
        for reaction in self.env.reactions.radiolytic:
            logger.debug(f"-> {reaction.label()}: {reaction.ge():.3f}")

    # ----------------------------- ODE integration ------------------------- #
    def _raw_run(
        self,
        t_span: tp.Tuple[float, float],
        options: RunOptions,
        y0: tp.Optional[NDArray[np.float64]] = None,
        beam_on: bool = True,
        dose_rate: tp.Optional[float] = None,
    ) -> ScipyODEResult:
        """
        Integrate the ODE system over t_span with solve_ivp.

        Steps: check the options, widen an empty span (widen_empty_span),
        build the solve_ivp keyword arguments (_solver_kwargs), apply the
        initial radiolytic burst if the beam is switched ON at t_span[0]
        (_apply_initial_burst), then solve and time (_timed_solve).

        Args:
            t_span: Tuple[float, float]
                Start and stop times [s]
            options: RunOptions
                Solver options
            y0: NDArray[np.float64], optional
                Initial concentrations of the dynamic species [µmol/L].
                If None, SimEnv.initial_values() is used. Never modified:
                the integration starts from a copy.
            beam_on: bool
                Whether the beam is switched ON at t_span[0] (default:
                True). The initial burst is applied only if the beam is
                also used (options.use_beam), see _apply_initial_burst.
            dose_rate: float, optional
                Dose rate [Gy/s] held constant over t_span (chunk of a
                chunked run, and of its initial burst). None -> the
                solver reads the beam at every time (plain run).

        Returns:
            ScipyODEResult: raw solve_ivp result

        Raises:
            NoBeamDefinedError: no beam is defined
            SolverFailedError: solve_ivp did not reach t_span[1]
        """
        self._require_beam()

        t_span = widen_empty_span(t_span)
        # First step is clamped against the span *before* the burst shift.
        solver_kwargs = self._solver_kwargs(t_span, options)

        # Always work on a copy: the burst must not alter the caller's y0
        # (nor the last stored point of a previous chunk).
        y_start = (
            np.array(y0, dtype=np.float64, copy=True)
            if y0 is not None
            else self.env.initial_values()
        )
        beam: tp.Optional[DoseRateSource] = None
        if options.use_beam:
            beam = (
                self.beam_list.current
                if dose_rate is None
                else ConstantDoseRate(dose_rate)
            )

        self.last_burst_dose_rate = None
        if beam_on and options.use_beam:  # beam switched ON at t_span[0]
            t_span = self._apply_initial_burst(t_span, y_start, dose_rate)
            options = self._drop_skipped_t_eval(t_span, options)

        sol = self._timed_solve(t_span, y_start, beam, options, solver_kwargs)
        if not sol.success:
            t_last = float(sol.t[-1]) if len(sol.t) > 0 else t_span[0]
            raise rp_excep.SolverFailedError(
                f"solve_ivp ({options.method}) failed on t_span={t_span}: "
                f"last reached t={t_last:.6e} s; scipy message: {sol.message}"
            )
        return sol

    def _solver_kwargs(
        self,
        t_span: tp.Tuple[float, float],
        options: RunOptions,
    ) -> tp.Dict[str, tp.Any]:
        """
        Step tolerances of solve_ivp: rtol, atol (defaults DEFAULT_RTOL,
        DEFAULT_ATOL), first_step (see RunOptions.first_step) and
        max_step (omitted when None, i.e. solve_ivp's inf).
        """
        first_step = (
            options.first_step if options.first_step is not None else self.FIRST_STEP
        )
        kwargs: tp.Dict[str, tp.Any] = {}
        kwargs["first_step"] = first_step
        kwargs["rtol"] = options.rtol if options.rtol is not None else DEFAULT_RTOL
        kwargs["atol"] = options.atol if options.atol is not None else DEFAULT_ATOL

        if options.max_step is not None:
            kwargs["max_step"] = options.max_step
            # The first step must be lower than the max step
            if first_step >= options.max_step:
                kwargs["first_step"] = options.max_step / 10

        # The first step must fit in the integration span
        if kwargs["first_step"] > t_span[1] - t_span[0]:
            kwargs["first_step"] = (t_span[1] - t_span[0]) / 1e3
        return kwargs

    def _apply_initial_burst(
        self,
        t_span: tp.Tuple[float, float],
        y0: NDArray[np.float64],
        dose_rate: tp.Optional[float] = None,
    ) -> tp.Tuple[float, float]:
        """
        Beam switched ON at t_span[0]: skip the first T_ZERO seconds and
        add *in place* to y0 [µmol/L] the radiolytic production of the
        dynamic species over that interval (computed analytically from
        the G values, at `dose_rate` or, when None, at the dose rate of
        the beam at t_span[0] + T_ZERO; stored in last_burst_dose_rate).
        Constant species are left untouched.

        Returns:
            Tuple[float, float]: (t_span[0] + T_ZERO, t_span[1])

        Raises:
            ValueError: the span is shorter than T_ZERO
        """
        t_span = (t_span[0] + self.T_ZERO, t_span[1])
        if t_span[0] >= t_span[1]:
            raise ValueError(f"Too Short beam pulse {t_span}")
        if dose_rate is None:
            dose_rate = self.beam.at(t_span[0]).dose_rate()
        self.last_burst_dose_rate = dose_rate
        for r in self.env.reactions.radiolytic:
            # If the species must be forced constant, skip it.
            if self.env.species.is_cst_species(r.species()):
                continue

            # Keep initial value and increment
            y0[r.species().index] += (
                r.compute_initial_cc_before(t0=self.T_ZERO, dose_rate=dose_rate)
                * 1e6
            )  # Units of initial cc are µM
        return t_span

    @staticmethod
    def _drop_skipped_t_eval(
        t_span: tp.Tuple[float, float],
        options: RunOptions,
    ) -> RunOptions:
        """
        The initial burst moves the start of the integration to
        t_span[0] + T_ZERO: the t_eval times before it (inside the
        skipped interval) cannot be evaluated by solve_ivp and are
        dropped.
        """
        if options.t_eval is None:
            return options
        t_eval = np.asarray(options.t_eval)
        kept = t_eval[t_eval >= t_span[0]]
        if len(kept) == len(t_eval):
            return options
        logger.warning(
            f"{len(t_eval) - len(kept)} t_eval time(s) before "
            f"{t_span[0]} s (end of the skipped initial burst) dropped"
        )
        return dataclasses.replace(options, t_eval=kept)

    def _timed_solve(
        self,
        t_span: tp.Tuple[float, float],
        y0: NDArray[np.float64],
        beam: tp.Optional[DoseRateSource],
        options: RunOptions,
        solver_kwargs: tp.Dict[str, tp.Any],
    ) -> ScipyODEResult:
        """
        Call solve_ivp on derive() and store its wall time in last_exec_time.

        With options.use_jac, jacobian() (same arguments as derive) is
        passed as jac, only for the methods that use one (LSODA, BDF,
        Radau); other methods ignore use_jac without any warning.
        """
        if options.use_jac and options.method in JAC_METHODS:
            solver_kwargs = {**solver_kwargs, "jac": jacobian}
        t0 = time.time()
        sol: ScipyODEResult = solve_ivp(
            derive,
            t_span,
            y0,
            args=(beam, self.mat),
            t_eval=options.t_eval,
            method=options.method,
            **solver_kwargs,
        )
        t1 = time.time()
        self.last_exec_time = t1 - t0
        return sol

    # ------------------------------ Public runs ---------------------------- #
    def run_till_max_dose(
        self,
        time_start: float = 1e-9,
        **kwargs: tp.Unpack[KwargsRun],
    ) -> ODEResult:
        """
        Run from time_start up to the time at which the current beam
        reaches its max dose.

        Args:
            time_start: float
                Start time [s] (default: 1e-9)
            **kwargs: KwargsRun
                Other keyword arguments of run(). A 't_span' given here
                is overwritten.

        Returns:
            ODEResult

        Raises:
            ValueError: the beam has no max dose (time at max dose = inf)
        """
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
        use_jac: bool = True,
        method: str = "LSODA",
        atol: tp.Optional[float] = None,
        rtol: tp.Optional[float] = None,
        max_step: tp.Optional[float] = None,
        y0: tp.Optional[NDArray[np.float64]] = None,
        use_beam: bool = True,
    ) -> ODEResult:
        """Solve ODE over t_span

        Args:
            t_span: Tuple[float, float]
                t_span[0] -> start time, in seconds
                t_span[1] -> End time, in seconds
            t_eval: [NDArray[np.float64], optional
                Force function evaluation on t_eval (default: None)
            use_jac: bool
                Use the analytical Jacobian? Only for LSODA, BDF and
                Radau, ignored by other methods (default: True).
            method: str
                Which ODE solver to use? Define in
                scipy.integrate.solve_ivp (default: "LSODA")
            atol: float, optional
                Absolute Tolerance of the ODE Solver (default: 1e-12)
            rtol: float, optional
                Relative Tolerance of the ODE Solver (default: 1e-10)
            max_step: float, optional
                Larger authorized time step for the ODE solver
                (default: inf)
            y0: NDArray[np.float64], optional
                Initial concentration of DynSpecies [µmol/L] (default:
                None). If None, calls EnvSim.initial_values().
            use_beam: bool, optional
                whether to use the beam or not (default: True).

        The beam is assumed switched ON at t_span[0]: unless
        use_beam=False, the first T_ZERO seconds are skipped and the
        radiolytic production over them is added to y0 (a copy: the
        given y0 is never modified).

        Returns:
            ODEResult: integration result (solve_ivp solution wrapped
            with the environment, beam, matrices, LET and T_ZERO)

        Raises:
            NoBeamDefinedError: no beam is defined
            ValueError: beam ON at t_span[0] and t_span not longer than
                T_ZERO (the initial burst would cover the whole span)
        """
        logger.info(
            f"Starting simulation with t_span: {t_span} (jacobian use: {use_jac})"
        )
        options = RunOptions(
            atol=atol,
            rtol=rtol,
            method=method,
            max_step=max_step,
            use_beam=use_beam,
            use_jac=use_jac,
            t_eval=t_eval,
        )
        sol = self._raw_run(t_span, options, y0=y0, beam_on=True)
        logger.info(f"Simulation done! ({self.last_exec_time:.2g} s)")
        res = ODEResult(
            self.env,
            self.beam_list.current,
            self.mat,
            self.beam_list.current_LET,
            self.T_ZERO,
        )
        res.append_sol(
            sol,
            exec_time=self.last_exec_time,
            bursts=self._last_bursts(),
        )
        return res

    def prepare_chunked_run(
        self,
        t_span: tp.Tuple[float, float],
        use_jac: bool = True,
        method: str = "LSODA",
        atol: tp.Optional[float] = None,
        rtol: tp.Optional[float] = None,
        max_step_size_on: tp.Optional[float] = None,
        max_step_size_off: tp.Optional[float] = None,
        max_step_size_final: tp.Optional[float] = None,
        y0: tp.Optional[NDArray[np.float64]] = None,
        use_beam: bool = True,
    ) -> SimChunksIter:
        """Method used to prepare a multi step ODE integration, one chunk
        per beam ON / beam OFF interval of the current beam (from
        beam.make_time_intervals). The initial radiolytic burst is applied
        at the start of every beam ON chunk.

        Args:
            t_span: Tuple[float, float]
                t_span[0] -> start time, in seconds
                t_span[1] -> End time, in seconds
            use_jac: bool, optional
                Whether to use the analytical jacobian along with the ODE
                solver (default: True). Only for LSODA, BDF and Radau,
                ignored by other methods.
            method: str, optional
                Which ODE solver to use? Define in
                scipy.integrate.solve_ivp (default: "LSODA")
            atol: float, optional
                Absolute Tolerance of the ODE Solver (default: 1e-12)
            rtol: float, optional
                Relative Tolerance of the ODE Solver (default: 1e-10)
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
                Initial concentration of DynSpecies [µmol/L] (default:
                None). If None, calls EnvSim.initial_values()
            use_beam: bool, optional
                whether to use the beam or not (default: True).

        Returns:
            SimChunksIter: Iterable object performing the ODE
            integration.

        Raises:
            NoBeamDefinedError: no beam is defined
        """
        self._require_beam()
        intervals = self.beam_list.current.make_time_intervals(
            start=t_span[0], stop=t_span[1]
        )
        # intervals: one row per boundary, [time, 1 if beam ON else 0];
        # the last row only holds the stop time.
        beam_on = intervals[:-1, 1] > 0
        max_steps = [
            max_step_size_on if on else max_step_size_off for on in beam_on[:-1]
        ]
        # Final chunk: max_step_size_final, if any, overrides on/off.
        max_step = max_step_size_final
        if max_step is None:
            max_step = max_step_size_on if intervals[-2, 1] > 0 else max_step_size_off
        max_steps.append(max_step)

        chunk_list = _make_chunks(
            intervals[:, 0], beam_on, max_steps, self.beam_list.current
        )
        return self._chunks_iter(
            t_span,
            chunk_list,
            RunOptions(
                atol=atol,
                rtol=rtol,
                method=method,
                use_beam=use_beam,
                use_jac=use_jac,
            ),
            y0,
            n_beam_intervals=len(intervals),
        )

    def prepare_single_pulse(
        self,
        t_span: tp.Tuple[float, float],
        use_jac: bool = True,
        method: str = "LSODA",
        atol: tp.Optional[float] = None,
        rtol: tp.Optional[float] = None,
        max_step_size_on: tp.Optional[float] = None,
        max_step_size_trans: tp.Optional[float] = None,
        max_step_size_final: tp.Optional[float] = None,
        trans_time: tp.Optional[float] = None,
        y0: tp.Optional[NDArray[np.float64]] = None,
        use_beam: bool = True,
    ) -> SimChunksIter:
        """Method used to prepare a multi step ODE integration of a
        single pulse (tm = time at max dose of the current beam):
            1. the beam intervals over [t_span[0], tm] (beam.
               make_time_intervals, as prepare_chunked_run: a leading
               beam OFF chunk if the beam starts after t_span[0], one ON
               chunk per pulse, OFF chunks in between; the initial burst
               is applied at the start of every ON chunk);
            2. ON -> OFF transition over [tm, tm + trans_time];
            3. beam OFF over [tm + trans_time, t_span[1]].
        Every boundary is clipped to t_span[1] (empty chunks dropped).
        With a single pulse, 1. is one ON chunk.

        Args:
            t_span: Tuple[float, float]
                t_span[0] -> start time, in seconds
                t_span[1] -> End time, in seconds
            use_jac: bool, optional
                Whether to use the analytical jacobian along with the ODE
                solver (default: True). Only for LSODA, BDF and Radau,
                ignored by other methods.
            method: str, optional
                Which ODE solver to use? Define in
                scipy.integrate.solve_ivp (default: "LSODA")
            atol: float, optional
                Absolute Tolerance of the ODE Solver (default: 1e-12)
            rtol: float, optional
                Relative Tolerance of the ODE Solver (default: 1e-10)
            max_step_size_on: float, optional
                Larger authorized time step for the ODE solver while
                beam is ON before tm (default: inf)
            max_step_size_trans: float, optional
                Larger authorized time step for the ODE solver during
                the ON -> OFF transition (default: inf)
            max_step_size_final: float, optional
                Larger authorized time step for the ODE solver for the
                last (beam OFF) step, and for the beam OFF chunks before
                tm (default: inf, no fallback to the other max step
                sizes).
            trans_time: float, optional
                Duration of the ON -> OFF transition, in seconds
                (default: 100e-6)
            y0: NDArray[np.float64], optional
                Initial concentration of DynSpecies [µmol/L] (default:
                None). If None, calls EnvSim.initial_values()
            use_beam: bool, optional
                whether to use the beam or not (default: True).

        Returns:
            SimChunksIter: Iterable object performing the ODE
            integration.

        Raises:
            NoBeamDefinedError: no beam is defined
        """
        self._require_beam()
        trans_time = trans_time if trans_time is not None else DEFAULT_TRANS_TIME
        t_max_dose = self.beam_list.current.time_at_max_dose
        t_beam_off = t_max_dose + trans_time

        t_start, t_stop = t_span
        beam = self.beam_list.current

        # Before tm: the beam intervals up to tm (or t_span[1]), exactly as
        # prepare_chunked_run (leading OFF interval of a beam starting after
        # t_span[0], several pulses, gaps between them).
        boundaries = [t_start]
        beam_on: tp.List[bool] = []
        max_steps: tp.List[tp.Optional[float]] = []
        t_pre = min(t_max_dose, t_stop)
        if t_pre > t_start:
            intervals = beam.make_time_intervals(start=t_start, stop=t_pre)
            boundaries = intervals[:, 0].tolist()
            beam_on = (intervals[:-1, 1] > 0).tolist()
            max_steps = [
                max_step_size_on if on else max_step_size_final for on in beam_on
            ]
        # Then the ON -> OFF transition and the final OFF chunk, clipped to
        # t_span[1] (a chunk beyond it, or empty, is dropped).
        for stop, max_step in (
            (t_beam_off, max_step_size_trans),
            (t_stop, max_step_size_final),
        ):
            stop = min(stop, t_stop)
            if stop > boundaries[-1]:
                boundaries.append(stop)
                beam_on.append(False)
                max_steps.append(max_step)

        chunk_list = _make_chunks(
            boundaries=boundaries,
            beam_on=beam_on,
            max_steps=max_steps,
            beam=beam,
        )
        return self._chunks_iter(
            t_span,
            chunk_list,
            RunOptions(
                atol=atol,
                rtol=rtol,
                method=method,
                use_beam=use_beam,
                use_jac=use_jac,
            ),
            y0,
        )

    def _last_bursts(self) -> tp.List[float]:
        """Dose rates of the bursts applied by the last integration step."""
        rate = self.last_burst_dose_rate
        return [] if rate is None else [rate]

    def _require_beam(self) -> None:
        """Raise NoBeamDefinedError (and log it) when no beam is defined."""
        if len(self.beam_list) <= 0:
            logger.error("No beam is defined for the simulation")
            raise rp_excep.NoBeamDefinedError

    def _chunks_iter(
        self,
        t_span: tp.Tuple[float, float],
        chunk_list: tp.List[SimChunk],
        options: RunOptions,
        y0: tp.Optional[NDArray[np.float64]],
        n_beam_intervals: tp.Optional[int] = None,
    ) -> SimChunksIter:
        """Log the chunk layout and wrap the chunks in a SimChunksIter."""
        logger.info(f"Beam span ==> [{t_span[0]} - {t_span[1]}]")
        if n_beam_intervals is not None:
            logger.info(f"From beam chunks ==> {n_beam_intervals}")
        logger.info(f"Simulation divided in {len(chunk_list)} chunks")
        return SimChunksIter(sim=self, chunks=chunk_list, options=options, y0=y0)


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
        max_step_size: float, optional
            Max step size allowed over this integration step (None ->
            inf)
        beam_on: bool
            Whether the beam is switched ON at start (the initial
            radiolytic burst is then applied)
        dose_rate: float, optional
            Dose rate [Gy/s] of the beam inside the chunk, held constant
            by the solver (None -> the beam is read at every time)
    """

    start: float
    stop: float
    max_step_size: float | None
    beam_on: bool
    dose_rate: float | None = None


class SimChunksIter:
    """Iterable class performing a multi-step integration.

    Attributes:
        sim: UniCell
            Simulation volume in which simulation takes place
        chunks: List[SimChunk]
            List of simulation steps/chunks
        sol: ODEResult
            Aggregate of ODE results from scipy.integrate.solve_ivp
        current: int
            Keeps tracks of the current integration step.
        options: RunOptions
            Solver options shared by every chunk (the max step of each
            chunk replaces options.max_step).
        y0: NDArray[np.float64], optional
            Initial concentration of DynSpecies [µmol/L] (default:
            None). If None, calls EnvSim.initial_values(). Later chunks
            start from the last point of the previous one.
    Usage:
        >>> res = cell.prepare_chunked_run(t_span=(0, 1)).run()
    """

    sim: UnitCell
    chunks: tp.List[SimChunk]
    sol: ODEResult
    current: int
    options: RunOptions
    y0: tp.Optional[NDArray[np.float64]]

    def __init__(
        self,
        sim: UnitCell,
        chunks: tp.List[SimChunk],
        options: tp.Optional[RunOptions] = None,
        y0: tp.Optional[NDArray[np.float64]] = None,
    ):
        self.sim = sim
        self.chunks = chunks
        self.sol = ODEResult(
            self.sim.env,
            self.sim.beam_list.current,
            self.sim.mat,
            self.sim.beam_list.current_LET,
            self.sim.T_ZERO,
        )
        self.options = options if options is not None else RunOptions()
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
            (chunk.start, chunk.stop),
            dataclasses.replace(self.options, max_step=chunk.max_step_size),
            # _raw_run copies y0: the stored last point stays untouched
            y0=self.sol.y[:, -1] if self.current > 0 else self.y0,
            beam_on=chunk.beam_on,
            dose_rate=chunk.dose_rate,
        )
        self.sol.append_sol(
            res,
            exec_time=self.sim.last_exec_time,
            bursts=self.sim._last_bursts(),
            dose_rate=chunk.dose_rate if self.options.use_beam else None,
        )
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
