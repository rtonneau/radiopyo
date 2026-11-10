#!/usr/bin/env python3
"""
The beam submodule contains the radiation source (particle beam)
definitions. The timers live in `timer.py` and the beam container in
`beam_collection.py`; both are re-exported here for convenience.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging
import typing as tp
from dataclasses import dataclass, field
from math import floor, isclose, isinf

import numpy as np

# ----------------------------- Internal imports ---------------------------- #
from .exceptions import ParameterError
from .timer import (
    COUNT_RTOL,
    BaseTimer,
    MacroTimer,
    Timer,
    TimerState,
    split_periods,
)

# ------------------------- Type checking ONLY imports----------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

# -------------------------- Expose importable stuffs ----------------------- #
__all__ = [
    "BEAM_TYPES",
    "BaseTimer",
    "BeamCollection",
    "ConstantBeam",
    "MacroPulsedBeam",
    "MacroTimer",
    "ParticleBeam",
    "PulsedBeam",
    "SinglePulseBeam",
    "TimeMessage",
    "Timer",
    "TimerState",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class TimeMessage(tp.NamedTuple):
    """
    NamedTuple representing a dose rate at a specific time.

    Attributes:
        time (float): Time in seconds.
        current_dose_rate (float): Dose rate in Gy/s at the given time.

    This named tuple is used to encapsulate the dose rate and the time
    at which it was measured. It is typically returned by methods of
    radiation sources classes to provide dose rate information at a
    specific time.
    """

    time: float
    current_dose_rate: float

    def dose_rate(self) -> float:
        """
        Convenient method to shorten (+ encapsulate) the call to
        self.current_dose_rate.

        Returns:
            float: dose rate in Gy/s
        """
        return self.current_dose_rate


@dataclass
class ParticleBeam:
    """
    Generic definition of a Radiation Source via an instance of Timer
    and dose rate/dose max pieces of information. ParticleBeam should
    not be used as is but via class heritage instead.

    Attributes:
        dose_rate: float
            Average dose rate in Gy/s
        timer: Timer
            Instance of Timer NamedTuple
        max_dose: float
            integrated dose from which dose rate will drop to 0.0 Gy/s
            regardless of the Timer state.
        time_at_max_dose: float
            Time at which the maximal dose will be reached.
        dose_rate_peak: float
            Dose rate during on_time of the Timer.
    """

    dose_rate_peak: float = field(init=True)
    timer: BaseTimer = field(init=True)
    max_dose: float = field(default=float("inf"))
    label: str = field(default="default")
    time_at_max_dose: float = field(init=False)
    dose_rate: float = field(init=False)
    GROWTH_TIME: tp.ClassVar[float] = 1e-6

    def __post_init__(self) -> None:
        """
        Method called after ParticleBeam__init__() to
        finalize initialization
        """
        if self.dose_rate_peak <= 0:
            raise ValueError(
                f"Null or Negative inst. dose rate for ConstantBeam {self}"
            )
        self.time_at_max_dose = self._compute_time_at_max_dose()
        self.dose_rate = self._compute_average_dose_rate()

    def set_name(self, name: str) -> None:
        """Gives the beam a name"""
        self.label = name

    def beam_type(self) -> str:
        """
        Returns beam type as a string.
        Returns:
            str
                Beam type
        """
        raise NotImplementedError

    def update_beam_param(self, **kwargs: float) -> None:
        """
        Modify some beam parameter and recomputes what needs to be
        changed. This method should be overridden for new ParticleBeam
        inherited class. This one takes care of dose_rate, max_dose and
        label(name).

        Args:
            kwargs: float
                key words argument containing float value
                (except for 'label').
        Raises:
            AttributeError
                if kwargs contains unknown key values.
        """
        for key, value in kwargs.items():
            if key not in ["dose_rate_peak", "max_dose", "label"]:
                raise AttributeError(f"ParticleBeam has not attribute: {key}")
            setattr(self, key, value)
        self.__post_init__()

    def _compute_time_at_max_dose(self, max_dose: tp.Optional[float] = None) -> float:
        """Helper method to compute the time at which max dose rate will
        be reached with the peak dose rate and Timer defined.

        Returns:
            float: time in seconds.
        """
        max_dose = self.max_dose if max_dose is None else max_dose

        # if max_dose == inf, max_dose is never reached
        if max_dose >= float("inf") or self.dose_rate_peak <= 0:
            return float("inf")
        if max_dose == 0:
            return 0.0
        # if self.timer.period is inf -> no repetition: max_dose
        # delivered in ONE pulse. Continuous beam and mono_pulse beam
        # will end up here.In this case, (average) dose rate and peak
        # dose rate are the same.
        if self.timer.period >= float("inf"):
            return self.timer.start_time + (max_dose / self.dose_rate_peak)

        # From here, radiation source is pulsed with period < inf
        # How many pulses to reach max_dose?
        #   ==> (max dose) / (dose per pulse)

        # How many full pulses needed?
        n = max_dose // self.dose_per_pulse()

        # Fraction of the last pulse needed to reached max_dose
        last_frac = round((max_dose % self.dose_per_pulse()) / self.dose_per_pulse(), 6)

        # if n == 0, max dose reached in less than a single pulse:
        if n <= 0.0:
            return self.timer.start_time + (last_frac * self.timer.on_time)

        # From now on n >= 1.0 i.e. max_dose delivered in more than one
        #  pulse.If max_dose delivered with entire periods then last one
        # only last 'on_time'
        if last_frac <= 0.0:
            return (
                self.timer.start_time + (n - 1) * self.timer.period + self.timer.on_time
            )

        # From now on n >= 1.0 AND 0.0 < n_remainder < 1.0
        return (
            self.timer.start_time
            + n * self.timer.period
            + last_frac * self.timer.on_time
        )

    def _compute_average_dose_rate(self) -> float:
        """
        Average dose rate: delivered dose over the time elapsed between
        the beam start and `time_at_max_dose`. Without a max dose (never
        reached) it is the peak dose rate times the timer duty cycle.
        """
        if isinf(self.timer.period):
            return self.dose_rate_peak
        duration = self.time_at_max_dose - self.timer.start_time
        if isinf(duration) or duration <= 0:
            return self.dose_rate_peak * self.timer.duty_cycle()
        return self.dose_rate_peak * (
            self.tot_pulses() * self.timer.on_time / duration
        )

    def dose_per_pulse(self) -> float:
        """ """
        return self.dose_rate_peak * self.timer.on_time

    def average_dose_rate(self) -> float:
        """
        Method to make it clear we are dealing with average dose rate.

        Returns:
            float: average dose rate in Gy/s
        """
        return self.dose_rate

    def at(self, time: float) -> TimeMessage:
        """
        Computes instant dose rate at a given time.
        Args:
            time: float
                Time, in seconds, at which one want to know the dose
                rate in Gy/s

        Returns:
            TimeMessage
        """
        dr = 0.0
        if self.timer.start_time <= time < self.time_at_max_dose:
            match self.timer.state_at(time):
                case TimerState.isOff:
                    dr = 0.0
                case TimerState.isOn:
                    dr = self.dose_rate_peak
        return TimeMessage(
            time=time,
            current_dose_rate=dr,
        )

    def tot_pulses(self) -> float:
        """
        Computes the total number of pulses delivered between the beam
        start and `time_at_max_dose` (the last one may be partial).

        Returns:
            float
                Number of pulses (inf if the max dose is never reached)
        """
        if isinf(self.timer.period):
            return 1.0
        duration = self.time_at_max_dose - self.timer.start_time
        if isinf(duration):
            return float("inf")
        if duration <= 0:
            return 0.0
        n, remainder = split_periods(duration, self.timer.period)
        return n + min(remainder / self.timer.on_time, 1.0)

    def _first_on_timer(self) -> BaseTimer:
        """Timer whose start time is the earliest possible beam ON."""
        return self.timer

    def _transitions(self, start: float) -> tp.Iterator[tp.Tuple[float, bool]]:
        """
        Beam state transitions from `start`, as (time, is_on) pairs. The
        first one is the state at `start`; the following ones alternate
        ON/OFF with increasing time (times accumulate period by period).
        A timer that never switches OFF yields (inf, False) last.
        """
        timer = self.timer
        if isinf(timer.period):  # continuous / single pulse
            if start < timer.start_time:
                yield start, False
                yield timer.start_time, True
            else:
                yield start, True
            yield float("inf"), False
            return
        if start < timer.start_time:
            yield start, False
            t_on = timer.start_time
        elif timer.state_at(start) is TimerState.isOn:
            yield start, True
            yield start + timer.time_till_next_OFF(start), False
            t_on = start + timer.time_till_next_period(start)
        else:
            yield start, False
            t_on = start + timer.time_till_next_period(start)
        while True:
            yield t_on, True
            yield t_on + timer.on_time, False
            t_on = t_on + timer.period

    def make_matrix(self, stop: float, start: float = 0.0) -> NDArray[np.float64]:
        """
        Method computing the radiation source time structure in the
        range [start, stop] i.e. computes 2 columns matrix whose
        columns are:
            0) Time in seconds
            1) Related instant. dose rate
        For example:
            [   #Col0 #Col1
                [0.0, 0.5],
                [0.5, 0.0],
                [1.5, 0.5],
                [2.0, 0.0],
                [3.0, 0.5],
                ...
            ]
        The way to read it is the following:
            - From 0.0 s to 0.5 s: dose rate is 0.5 Gy/s
            - From 0.5 s to 1.5 s: dose rate is 0.0 Gy/s
            - From 1.5 s tp 2.0 s: dose rate is 0.5 Gy/s
            - ...
        -> Representation of a Pulsed beam with:
            period: 1.5s
            on_time: 0.5 s
            peek_dose_rate: 0.5 Gy/s (average dose rate: 0.166 Gy/s)

        This matrix definition of the beam structure is useful to
        predefine integration steps width (i.e. used by the
        "make_time_intervals" method). The matrix is built by walking
        the timer transitions up to min(stop, time_at_max_dose); the
        last row is always [stop, 0.0].

        Args:
            stop: float
                End time of the beam structure, in seconds
            start: float, optional
                Start time of the beam structure, in seconds
                (default: 0.0)

        Returns:
            NDArray[np.float64]: 2 cols matrix
        """
        # Always OFF cases:
        if start > self.time_at_max_dose or stop < self._first_on_timer().start_time:
            return np.array([[start, 0.0], [stop, 0.0]], dtype=np.float64)

        end = min(stop, self.time_at_max_dose)
        walk = self._transitions(start)
        start_time, start_on = next(walk)
        rows = [[start_time, self.dose_rate_peak if start_on else 0.0]]
        for time, is_on in walk:
            if time < end:
                rows.append([time, self.dose_rate_peak if is_on else 0.0])
            elif not is_on:
                # The beam is ON at `end` (max dose reached or stop): cut.
                rows.append([end, 0.0])
                break
            else:
                break
        if rows[-1][0] < stop:
            rows.append([stop, 0.0])
        else:
            rows[-1][1] = 0.0
        return np.array(rows, dtype=np.float64)

    def make_time_intervals(
        self,
        stop: float,
        start: float = 0.0,
    ) -> NDArray[np.float64]:
        """
        Create a 2D matrix similar to the one produced by the
        make_matrix(stop, start) method. The major differences are:
            - takes into account max_dose in order to trim the matrix
            if needed.
            - Symmetrically expand the 'on_time' so that it is possible
            to define a refined ODE integration time step for beam
            OFF-ON & ON-OFF transitions.
            - The first column is still time but the second is:
                * 0 for beam OFF
                * 1 for beam ON

        Returns:
            NDArray[np.float64]: 2 cols matrix
        """
        # Retrieve Time Structure Matrix
        out = self.make_matrix(stop, start)

        # Shift every OFF row that is neither the first nor the last one
        # (pulse ends) by GROWTH_TIME, unless the next row is closer:
        # the beam ON time is expanded to refine the integration step
        # at OFF transitions.
        beam_on = out[:, 1] > 0.0
        for idx in range(1, out.shape[0] - 1):
            if not beam_on[idx] and (out[idx + 1, 0] - out[idx, 0]) > self.GROWTH_TIME:
                out[idx, 0] += self.GROWTH_TIME

        out[:, 1] = beam_on
        return out

    def as_dict(self) -> tp.Dict:
        """ """
        raise NotImplementedError

    @classmethod
    def from_dict(cls, param: tp.Dict[str, float]) -> ParticleBeam:
        """Build a beam out of a dictionary (see subclasses)."""
        raise NotImplementedError

    def serialize(self, **kwargs: float) -> dict:
        """ """
        out = self.as_dict()
        out.update(kwargs)
        return {self.beam_type(): out}


class ConstantBeam(ParticleBeam):
    """Class representing a Constant Radiation Source
    Attributes:
        dose_rate: float
            Average dose rate in Gy/s
        max_dose: float
            Max dose to deliver in Gy
    """

    def __init__(
        self,
        dose_rate: float,
        max_dose: tp.Optional[float] = None,
        start_time: tp.Optional[float] = None,
    ) -> None:
        super().__init__(
            dose_rate,
            Timer.new_constant(start_time),
            max_dose=max_dose if max_dose is not None else float("inf"),
        )

    def __str__(self) -> str:
        return f"ConstantBeam(dose_rate={self.dose_rate:.2g}, max_dose={self.max_dose})"

    def update_beam_param(self, **kwargs: float) -> None:
        """
        Add start_time as possible keyword.
        """
        if "start_time" in kwargs:
            self.timer = Timer.new_constant(kwargs.pop("start_time"))
        super().update_beam_param(**kwargs)

    def beam_type(self) -> str:
        return "Constant"

    def as_dict(self) -> tp.Dict[str, float]:
        """ """
        return {
            "dose_rate": self.dose_rate,
            "max_dose": self.max_dose,
            "start_time": self.timer.start_time,
        }

    @classmethod
    def from_dict(cls, param: tp.Dict[str, float]) -> ConstantBeam:
        """Build a ConstantBeam out of a dictionary"""

        if "on_time" in param and "max_dose" in param:
            dose_rate = param["max_dose"] / param["on_time"]
        else:
            dose_rate = 0.0
            for key in ["dose_rate_peak", "dose_rate_avg", "dose_rate"]:
                if key in param:
                    dose_rate = param[key]
                    break

        return ConstantBeam(
            dose_rate=dose_rate,
            max_dose=param.get("max_dose"),
            start_time=param.get("start_time"),
        )


class SinglePulseBeam(ConstantBeam):
    """Class representing a Single Pulse Radiation Source
    Attributes:
        dose_rate: float
            Average dose rate in Gy/s
        max_dose: float
            Max dose to deliver in Gy
        start_time: float, optional
            Time to start beam ON (default: 0.0)

    TODO: remove this class as it just mimic ConstantBeam?
    /!\\ useful constructor from'on_time' though
    """

    def __init__(
        self,
        dose_rate: float,
        max_dose: float,
        start_time: tp.Optional[float] = None,
    ) -> None:
        super().__init__(
            dose_rate,
            max_dose=max_dose,
            start_time=start_time,
        )

    def __str__(self) -> str:
        return (
            f"SinglePulseBeam(dose_rate={self.dose_rate:.2g}, "
            f"max_dose={self.max_dose}, "
            f"start_time={self.timer.start_time})"
        )

    def beam_type(self) -> str:
        return "singlepulse"

    @classmethod
    def from_on_time(
        cls,
        dose_rate: float,
        on_time: float,
        start_time: tp.Optional[float] = None,
    ) -> SinglePulseBeam:
        return cls(
            dose_rate,
            max_dose=on_time * dose_rate,
            start_time=start_time,
        )

    @classmethod
    def from_dict(cls, param: tp.Dict[str, float]) -> SinglePulseBeam:
        """
        Create a SinglePulseBeam from a Dictionary.

        Args:
            param: Dict[str, float]
                Dictionary containing beam parameters

        Returns:
            SinglePulseBeam

        Raises:
            radiopyo.physics.exceptions.ParameterError
        """
        if "dose_rate" not in param:
            raise ParameterError(
                "SinglePulseBeam must be defined with a 'dose_rate' parameter"
            )

        if "max_dose" in param:
            return SinglePulseBeam(
                dose_rate=param["dose_rate"],
                max_dose=param["max_dose"],
                start_time=param.get("start_time"),
            )
        if "on_time" in param:
            return SinglePulseBeam.from_on_time(
                dose_rate=param["dose_rate"],
                on_time=param["on_time"],
                start_time=param.get("start_time"),
            )

        raise ParameterError(
            "SinglePulseBeam must be defined with either one of "
            "'on_time' or 'max_dose' parameter"
        )


class PulsedBeam(ParticleBeam):
    """Class representing a Pulsed Radiation Source. Helper functions
    are also available to easily create specific Pulsed Radiation
    Sources:
        - from_dose_rate_peak
        - from_dose_per_pulse

    Attributes:
        dose_rate: float
            Average dose rate in Gy/s
        max_dose: float
            Max dose to deliver in Gy
        period: float
            Repetition period of the source, in seconds
        on_time: float
            Time, within 1 period, during which irradiation is ON
        start_time: float, optional
            Time to start beam ON (default: 0.0)
    """

    def __init__(
        self,
        dose_rate_peak: float,
        period: float,
        on_time: float,
        max_dose: tp.Optional[float] = None,
        start_time: tp.Optional[float] = None,
    ):
        super().__init__(
            dose_rate_peak,
            Timer.new_pulsed(period, on_time, start_time),
            max_dose if max_dose is not None else float("inf"),
        )

    def _update_timers(self, kwargs: tp.Dict[str, float]) -> None:
        """
        Rebuild the timer(s) from the timing keywords found in `kwargs`
        (they are removed from it): start_time, on_time, period or
        frequency.
        """
        if (
            "start_time" in kwargs
            or "on_time" in kwargs
            or "period" in kwargs
            or "frequency" in kwargs
        ):
            start_time = kwargs.pop("start_time", self.timer.start_time)
            on_time = kwargs.pop("on_time", self.timer.on_time)
            if (freq := kwargs.pop("frequency", None)) is not None:
                period = 1.0 / freq
            else:
                period = kwargs.pop("period", self.timer.period)
            self.timer = Timer.new_pulsed(period, on_time, start_time)

    # @override waiting for this baby to show up (Python3.12 fall 2023!)
    # => https://peps.python.org/pep-0698/
    def update_beam_param(self, **kwargs: float) -> None:
        self._update_timers(kwargs)

        if "dose_per_pulse" in kwargs or "n_pulse" in kwargs:
            dpp = kwargs.pop("dose_per_pulse", self.dose_per_pulse())
            n_pulse = kwargs.pop("n_pulse", self.tot_pulses())
            kwargs["max_dose"] = n_pulse * dpp
            kwargs["dose_rate_peak"] = dpp / self.timer.on_time

        super().update_beam_param(**kwargs)

    def beam_type(self) -> str:
        return "pulsed"

    def as_dict(self) -> tp.Dict[str, float]:
        """ """
        return {
            "dose_rate_peak": self.dose_rate_peak,
            "period": self.timer.period,
            "on_time": self.timer.on_time,
            "max_dose": self.max_dose,
            "start_time": self.timer.start_time,
        }

    @classmethod
    def from_dose_rate_peak(
        cls,
        peak_dr: float,
        period: float,
        on_time: float,
        n_pulse: tp.Optional[float] = None,
        max_dose: tp.Optional[float] = None,
        start_time: tp.Optional[float] = None,
    ) -> PulsedBeam:
        """
        Class method building a PulsedBeam based on peak dose rate
          => dose rate during ON time.
        Args:
            peak_dr: float
                peak dose rate, in Gy/s
            period: float
                Pulse period, in seconds
            on_time: float
                Pulse ON time, in seconds
            n_pulse: float
                Total number of pulse
            max_dose: float, optional
                Max dose to deliver, in Gy (default is None -> inf)

        Returns:
            PulsedBeam
        """
        n_pulse = n_pulse if n_pulse is not None else float("inf")
        max_dose_n_pulse = on_time * n_pulse * peak_dr
        max_dose = max_dose if max_dose is not None else float("inf")

        if (
            not isinf(max_dose)
            and not isinf(max_dose_n_pulse)
            and not isclose(max_dose_n_pulse, max_dose, rel_tol=1e-9)
        ):
            raise ValueError(
                "Inconsistent max dose provided by both max dose and number of pulses"
                f" (max dose: {max_dose} || DPP: {max_dose_n_pulse})"
            )
        return cls(
            dose_rate_peak=peak_dr,
            period=period,
            on_time=on_time,
            max_dose=min(max_dose, max_dose_n_pulse),
            start_time=start_time,
        )

    @classmethod
    def from_dose_per_pulse(
        cls,
        pulse_dose: float,
        period: float,
        on_time: float,
        n_pulse: tp.Optional[float] = None,
        max_dose: tp.Optional[float] = None,
        start_time: tp.Optional[float] = None,
    ) -> PulsedBeam:
        """Class method building a PulsedBeam based on a dose per pulse
        Args:
            pulse_dose: float
                dose delivered during ONE single pulse
            period: float
                Pulse period, in seconds
            on_time: float
                Pulse ON time, in seconds
            n_pulse: float
                Total number of pulse
            max_dose: float, optional
                Max dose to deliver, in Gy (default is None -> inf)

        Returns:
            PulsedBeam
        """
        return PulsedBeam.from_dose_rate_peak(
            peak_dr=pulse_dose / on_time,
            period=period,
            on_time=on_time,
            n_pulse=n_pulse,
            max_dose=max_dose,
            start_time=start_time,
        )

    @classmethod
    def from_dict(cls, param: tp.Dict[str, float]) -> PulsedBeam:
        """
        Create a PulsedBeam from a Dictionary. The caller's dictionary
        is left untouched.

        Raises:
            radiopyo.physics.exceptions.ParameterError
        """
        param = dict(param)
        if "period" not in param and "frequency" not in param:
            raise ParameterError(
                "PulsedBeam must be defined with a 'period' or 'frequency' parameter"
            )

        if "frequency" in param:
            param["period"] = 1 / param["frequency"]

        if "pulse_dose" in param and "dose_rate_peak" in param:
            param["on_time"] = param["pulse_dose"] / param["dose_rate_peak"]
        elif "pulse_dose" in param and "on_time" in param:
            param["dose_rate_peak"] = param["pulse_dose"] / param["on_time"]
        elif "dose_rate_peak" in param and "on_time" in param:
            param["pulse_dose"] = param["dose_rate_peak"] * param["on_time"]
        elif "dose_rate" in param and "on_time" in param:
            timer = Timer.new_pulsed(
                period=param["period"],
                on_time=param["on_time"],
            )
            param["dose_rate_peak"] = param["dose_rate"] / timer.duty_cycle()
        else:
            raise ParameterError(
                "PulsedBeam must be defined with at least 2 parameters out of "
                "'dose_rate_peak', 'on_time' or 'pulse_dose' "
                f"Beam def was: {param.get('label', 'default')}"
            )

        if "dose_rate_peak" in param:
            return PulsedBeam.from_dose_rate_peak(
                param["dose_rate_peak"],
                period=param["period"],
                on_time=param["on_time"],
                n_pulse=param.get("n_pulse"),
                max_dose=param.get("max_dose"),
                start_time=param.get("start_time"),
            )
        if "pulse_dose" in param:
            return PulsedBeam.from_dose_per_pulse(
                param["pulse_dose"],
                period=param["period"],
                on_time=param["on_time"],
                n_pulse=param.get("n_pulse"),
                max_dose=param.get("max_dose"),
                start_time=param.get("start_time"),
            )
        raise NotImplementedError(f"Unknown Pulsed beam definition: {param}")


class MacroPulsedBeam(PulsedBeam):
    """
    Class representing a Macro/Micro pulsed Radiation Source: trains
    (macro pulses) of micro pulses. `timer` is the micro timer (its start
    time, `micro_start_time`, is relative to the start of a train) and
    `macro_timer` the train timer (`start_time`, `period`, `on_time`);
    the micro timer is also the `subtimer` of `macro_timer`. Both are
    always kept in sync.
    """

    macro_timer: MacroTimer

    def __init__(
        self,
        dose_rate_peak: float,
        period: float,
        on_time: float,
        micro_period: float,
        micro_on_time: float,
        max_dose: tp.Optional[float] = None,
        micro_start_time: tp.Optional[float] = None,
        start_time: tp.Optional[float] = None,
    ):
        self.macro_timer = MacroTimer.new_macropulse(
            period, on_time, micro_period, micro_on_time, start_time, micro_start_time
        )

        super().__init__(
            dose_rate_peak,
            micro_period,
            micro_on_time,
            max_dose if max_dose is not None else float("inf"),
            micro_start_time,
        )
        # One micro timer object, shared with the macro timer.
        self.timer = self.macro_timer.subtimer

    def beam_type(self) -> str:
        return "macro_pulsed"

    def _update_timers(self, kwargs: tp.Dict[str, float]) -> None:
        """
        Rebuild both timers from the timing keywords found in `kwargs`
        (they are removed from it). Same names as `as_dict`: period (or
        frequency), on_time and start_time are the macro ones;
        micro_period (or micro_frequency), micro_on_time and
        micro_start_time the micro ones.
        """
        names = (
            "start_time",
            "on_time",
            "period",
            "frequency",
            "micro_start_time",
            "micro_on_time",
            "micro_period",
            "micro_frequency",
        )
        if not any(name in kwargs for name in names):
            return
        macro, micro = self.macro_timer, self.timer
        if (freq := kwargs.pop("frequency", None)) is not None:
            period = 1.0 / freq
        else:
            period = kwargs.pop("period", macro.period)
        if (micro_freq := kwargs.pop("micro_frequency", None)) is not None:
            micro_period = 1.0 / micro_freq
        else:
            micro_period = kwargs.pop("micro_period", micro.period)
        self.macro_timer = MacroTimer.new_macropulse(
            period,
            kwargs.pop("on_time", macro.on_time),
            micro_period,
            kwargs.pop("micro_on_time", micro.on_time),
            kwargs.pop("start_time", macro.start_time),
            kwargs.pop("micro_start_time", micro.start_time),
        )
        self.timer = self.macro_timer.subtimer

    def tot_trains(self) -> float:
        """ """
        return self.max_dose / self.dose_per_train()

    def tot_pulses(self) -> float:
        """Computes the total number of pulses.
        Returns:
            float
                Number of pulses (.2 digits precision)
        """
        return round(self.max_dose / self.dose_per_pulse(), 2)

    def dose_per_train(self) -> float:
        return self.dose_per_pulse() * self.macro_timer.micro_pulses_per_train()

    def _compute_average_dose_rate(self) -> float:
        """
        Average dose rate: delivered dose over the time elapsed between
        the start of the first train and `time_at_max_dose`. Without a
        max dose (never reached) it is the peak dose rate times the
        macro timer duty cycle.
        """
        duration = self.time_at_max_dose - self.macro_timer.start_time
        if isinf(duration) or duration <= 0:
            return self.dose_rate_peak * self.macro_timer.duty_cycle()
        return self.dose_rate_peak * (
            self.tot_pulses() * self.timer.on_time / duration
        )

    def _compute_time_at_max_dose(self, max_dose: tp.Optional[float] = None) -> float:
        """
        Helper method to compute the time at which max dose rate will be
        reached with the average dose rate and Timer defined.

        Returns:
            float: time in seconds.
        """
        max_dose = self.max_dose if max_dose is None else max_dose

        if max_dose == 0:
            return 0

        # if max_dose == inf, max_dose is never reached
        if max_dose >= float("inf") or self.dose_rate_peak <= 0:
            return float("inf")

        # How many full pulse trains needed? (tolerance: 3 trains of
        # 1/3 Gy are exactly 1 Gy only up to rounding)
        dose_per_train = self.dose_per_train()
        n = floor(max_dose / dose_per_train + COUNT_RTOL)
        remaining = max_dose - n * dose_per_train
        if remaining <= COUNT_RTOL * dose_per_train:
            # Whole trains: the max dose is reached at the end of the
            # last micro pulse of the last train.
            return (
                self.macro_timer.start_time
                + (n - 1) * self.macro_timer.period
                + super()._compute_time_at_max_dose(max_dose=dose_per_train)
            )

        return (
            self.macro_timer.start_time
            + n * self.macro_timer.period
            + super()._compute_time_at_max_dose(max_dose=remaining)
        )

    def at(self, time: float) -> TimeMessage:
        """
        Computes instant dose rate at a given time.
        Args:
            time: float
                Time, in seconds, at which one want to know the
                dose rate in Gy/s

        Returns:
            TimeMessage
        """
        dr = 0.0
        if time < self.time_at_max_dose:
            match self.macro_timer.state_at(time):
                case TimerState.isOff:
                    dr = 0.0
                case TimerState.isOn:
                    dr = self.dose_rate_peak
        return TimeMessage(
            time=time,
            current_dose_rate=dr,
        )

    def as_dict(self) -> tp.Dict[str, float]:
        """ """
        return {
            "dose_rate_peak": self.dose_rate_peak,
            "period": self.macro_timer.period,
            "on_time": self.macro_timer.on_time,
            "micro_period": self.timer.period,
            "micro_on_time": self.timer.on_time,
            "max_dose": self.max_dose,
            "micro_start_time": self.timer.start_time,
            "start_time": self.macro_timer.start_time,
        }

    def _first_on_timer(self) -> BaseTimer:
        return self.macro_timer

    def _transitions(self, start: float) -> tp.Iterator[tp.Tuple[float, bool]]:
        """
        Beam state transitions from `start` (see `ParticleBeam`): trains
        start every `macro_timer.period` from `macro_timer.start_time`;
        micro pulses start every micro period, `micro_start_time` after
        the train start, and are cut by the end of the train.
        """
        macro = self.macro_timer
        micro = macro.subtimer
        is_on = macro.state_at(start) is TimerState.isOn
        yield start, is_on
        train_idx = max(0, floor((start - macro.start_time) / macro.period))
        tolerance = COUNT_RTOL * micro.on_time
        while True:
            train_start = macro.start_time + train_idx * macro.period
            train_end = train_start + macro.on_time
            pulse_idx = 0
            while True:
                t_on = train_start + micro.start_time + pulse_idx * micro.period
                if t_on >= train_end - tolerance:
                    break
                t_off = min(t_on + micro.on_time, train_end)
                for time, state in ((t_on, True), (t_off, False)):
                    if time >= start and state != is_on:
                        is_on = state
                        yield time, state
                pulse_idx += 1
            train_idx += 1

    @classmethod
    def from_dict(cls, param: tp.Dict[str, float]) -> MacroPulsedBeam:
        """
        Create a MacroPulsedBeam from a Dictionary. The caller's
        dictionary is left untouched.

        Raises:
            radiopyo.physics.exceptions.ParameterError
        """
        param = dict(param)
        # Check data for macro timer
        if "period" not in param and "frequency" not in param:
            raise ParameterError(
                "MacroPulsedBeam must be defined with a 'period' or 'frequency' parameter"
            )

        if "frequency" in param:
            param["period"] = 1 / param["frequency"]

        if "on_time" not in param:
            raise ParameterError(
                "MacroPulsedBeam must be defined with an 'on_time' parameter"
            )
        ##
        # Check data for micro timer
        if "micro_period" not in param and "micro_frequency" not in param:
            raise ParameterError(
                "MacroPulsedBeam must be defined with an 'micro_period' or 'micro_frequency'"
                " parameter"
            )
        if "micro_frequency" in param:
            param["micro_period"] = 1 / param["micro_frequency"]

        if "micro_on_time" not in param:
            raise ParameterError(
                "MacroPulsedBeam must be defined with an 'micro_on_time' parameter"
            )
        ##

        if "dose_rate_peak" not in param:
            raise ParameterError(
                "MacroPulsedBeam must be defined with a 'dose_rate_peak' parameter"
            )

        return MacroPulsedBeam(
            param["dose_rate_peak"],
            period=param["period"],
            on_time=param["on_time"],
            micro_period=param["micro_period"],
            micro_on_time=param["micro_on_time"],
            max_dose=param.get("max_dose"),
            start_time=param.get("start_time"),
            micro_start_time=param.get("micro_start_time"),
        )


# --------------------------------------------------------------------------- #
#                                BEAM REGISTRY                                #
# --------------------------------------------------------------------------- #
# Beam classes by (lower case) beam type name, as returned by
# `ParticleBeam.beam_type().lower()` and used as section name in beam
# configurations.
BEAM_TYPES: tp.Dict[str, tp.Type[ParticleBeam]] = {
    "constant": ConstantBeam,
    "pulsed": PulsedBeam,
    "singlepulse": SinglePulseBeam,
    "macro_pulsed": MacroPulsedBeam,
}

# Re-exported for backward compatibility. Imported last as beam_collection
# needs BEAM_TYPES (lazily) from this module.
from .beam_collection import BeamCollection
