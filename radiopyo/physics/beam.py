#!/usr/bin/env python3
"""
The beam submodule contains all class definitions related to radiation source
structures.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

import logging
import typing as tp
from contextlib import suppress
from dataclasses import dataclass, field
from enum import Enum
from math import isinf, isnan

import numpy as np

# --------------------------------- Internal imports --------------------------------- #
from .exceptions import ParameterError

# ----------------------------- Type checking ONLY imports---------------------------- #
if tp.TYPE_CHECKING:
    from numpy.typing import NDArray

# -------------------------------------- Logging ------------------------------------- #
logger = logging.getLogger("radiopyo")

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = ["TimeMessage",
           "TimerState",
           "Timer",
           "ParticleBeam",
           "ConstantBeam",
           "PulsedBeam",
           "BeamCollection"
           ]

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class TimerState(Enum):
    """ Convenient Enum use to define radiation source states, ON or OFF."""
    isOff = 0
    isOn = 1


class TimeMessage(tp.NamedTuple):
    """ NamedTuple containing the dose rate for a given time.
    Attributes:
        time: float
            Time expressed in second
        current_dose_rate: float
            Dose rate delivered at the given time
    """
    time: float
    current_dose_rate: float

    def dose_rate(self) -> float:
        """Convenient method to shorten (+ encapsulate) the call to
        self.current_dose_rate. 

        Returns:
            float: dose rate in Gy/s
        """
        return self.current_dose_rate


class Timer(tp.NamedTuple):
    """NamedTuple use as a binary timer for radiation source. Currently only two states
    are possible for a radiation source, either ON or OFF i.e. 0.0% or 100%
    This class holds the time structure of the radiation source, defined as pulses.
    Pulses are repeated with a given 'period' (seconds) and the timer keeps its ON state
    during 'on_time' seconds.
    For example:
        - period = 2 seconds
        - on_time = 0.5 second
    --> each 2 seconds, timer ON during 0.5 second (i.e. duty cycle of 25%)

    Timer should not be constructed from scratch but from helper functions (classmethod)
    => Either one of:
        - Timer.new_constant()
        - Timer.new_pulsed(period, on_time)

    Rem:
        NamedTuple is a convenient structure as it is immutable and lightweight.

    Usage:
        >>> timer = Timer.new_pulsed(period=2, on_time=0.5)
        >>> timer.state_at(0.25)      # -> TimerState.isOn
        >>> timer.state_at(0.50)      # -> TimerState.isOff
        >>> timer.state_at(1.50)      # -> TimerState.isOff
        >>> timer.state_at(2.00)      # -> TimerState.isOn
        >>> print(timer.duty_cycle()) # -> 0.25

        >>> timer = Timer.new_constant()
        >>> timer.state_at(0.25)      # -> TimerState.isOn
        >>> timer.state_at(0.50)      # -> TimerState.isOn
        >>> timer.state_at(1.50)      # -> TimerState.isOn
        >>> timer.state_at(2.00)      # -> TimerState.isOn
        >>> print(timer.duty_cycle()) # -> 1.00
    """
    period: float
    on_time: float
    start_time: float

    @staticmethod
    def _default_start_time() -> float:
        return 0.0

    @classmethod
    def new_constant(cls, start_time: tp.Optional[float] = None) -> Timer:
        """ Class method to generate an always ON timer"""
        stime = start_time if start_time is not None else cls._default_start_time()
        return cls(float("inf"), float("inf"), stime)

    @classmethod
    def new_pulsed(cls,
                   period: float,
                   on_time: float,
                   start_time: tp.Optional[float] = None) -> Timer:
        """ Class method to generate a pulsed timer
        Args:
            period: float
                Period (repetition) of the timer in seconds
            on_time: float
                Time during which timer is ON in seconds. on_time should be <= period.
        """
        stime = start_time if start_time is not None else cls._default_start_time()
        return cls(period, on_time, stime)

    def is_continuous(self) -> bool:
        """ Is timer continuous? 
        Returns:
            bool
        """
        return (self.period == float("inf") and self.on_time == float("inf"))\
            or (self.on_time >= self.period)

    def is_mono_pulse(self) -> bool:
        """ Is the timer defined as a single pulse? i.e. 
            - period = inf
            - 0 < on_time < period 
        Returns:
            bool
        """
        return self.period == float("inf") and 0 < self.on_time < float("inf")

    def duty_cycle(self) -> float:
        """Computes the duty cycle of the timer
        Returns:
            float
        """
        out = self.on_time / self.period
        # In case of "inf"/"inf" => nan
        if isnan(out) or isinf(out):
            return 1.0
        return out

    def state_at(self, time: float) -> TimerState:
        """
        Compute the state of the timer at a given time.
        Args:
            time: float
                Time in seconds
        Returns:
            TimerState: Enum being either isOn or isOff
        """
        if ((time >= self.start_time) and
                ((time-self.start_time) % self.period) <= self.on_time):
            return TimerState.isOn
        return TimerState.isOff

    def period_index(self, time: float) -> int:
        """
        Get the 'position' of the period containing 'time'.
        For example, for 3s period & start_time=0, t=8s then period_index = 2
        """
        time = time - self.start_time
        return int(time // self.period)

    def time_till_next_ON(self, time: float) -> float:
        """
        Computes time remaining before next pulse.
        """
        if self.state_at(time) is TimerState.isOn:
            return 0.0
        return self.time_till_next_period(time)

    def time_till_next_OFF(self, time: float) -> float:
        """
        Computes time remaining before beam is off. Returns 0 if beam already off.
        """
        if self.state_at(time) is TimerState.isOff:
            return 0.0
        time = time - self.start_time - self.period * self.period_index(time)
        return self.on_time - time

    def time_till_next_period(self, time: float) -> float:
        """
        Computes how much time is remaining before next period starts.
        """
        time = time - self.start_time - self.period * self.period_index(time)
        return self.period - time

    def copy(self, start_time: tp.Optional[float] = None) -> Timer:
        """Make a simple copy of the Timer"""
        stime = start_time if start_time is not None else Timer._default_start_time()
        return Timer(self.period,
                     self.on_time,
                     stime)

    def off_time(self) -> float:
        """Computes the OFF time"""
        return self.period - self.on_time


@dataclass
class ParticleBeam():
    """
    Generic definition of a Radiation Source via an instance of Timer and dose rate/
    dose max pieces of information. ParticleBeam should not be used as is but via class
    heritage instead.

    Attributes:
        dose_rate: float
            Average dose rate in Gy/s
        timer: Timer
            Instance of Timer NamedTuple
        max_dose: float
            integrated dose from which dose rate will drop to 0.0 Gy/s regardless of the
            Timer state.
        time_at_max_dose: float
            Time at which the maximal dose will be reached.
        dose_per_pulse: float
            Dose delivered during the on_time of the Timer.
        peak_dose_rate: float
            Dose rate during on_time of the Timer.
    """
    dose_rate: float = field(init=True)
    timer: Timer = field(init=True)
    max_dose: float = field(default=float("inf"))
    label: str = field(default="default")
    time_at_max_dose: float = field(init=False)
    dose_per_pulse: float = field(init=False)
    peak_dose_rate: float = field(init=False)
    GROWTH_FACTOR: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """ Method called after ParticleBeam__init__() to finalize initialization"""
        self.GROWTH_FACTOR = 0.1
        self.peak_dose_rate = self._compute_peak_dose_rate()
        self.dose_per_pulse = (self.peak_dose_rate * self.timer.on_time)
        self.time_at_max_dose = self.timer.start_time + self._compute_time_at_max_dose()

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
        Modify some beam parameter and recomputes what needs to be changed.
        This method should be overridden for new ParticleBeam inherited class.
        This one takes care of dose_rate, max_dose and label(name).

        Args:
            kwargs: float
                key words argument containing float value (except for 'label').
        Raises:
            AttributeError
                if kwargs contains unknown key values.
        """
        for key, value in kwargs.items():
            if key not in ["dose_rate", "max_dose", "label"]:
                raise AttributeError(f"ParticleBeam has not attribute: {key}")
            setattr(self, key, value)
        self.__post_init__()

    def _compute_time_at_max_dose(self) -> float:
        """Helper method to compute the time at which max dose rate will be reached with
        the average dose rate and Timer defined.

        Returns:
            float: time in seconds.
        """
        # if max_dose == inf, max_dose is never reached
        if self.max_dose >= float("inf") or self.peak_dose_rate <= 0:
            return float("inf")

        # if self.timer.period is inf -> no repetition: max_dose delivered in ONE pulse
        # Continuous beam and mon_pulse beam will end up here.
        if self.timer.period >= float("inf"):
            return self.max_dose / self.peak_dose_rate

        # From here, radiation source is pulsed with period < inf
        # How many pulses to reach max_dose?
        #   ==> (max dose) / (dose per pulse)

        # How many full pulses needed?
        n = self.max_dose // self.dose_per_pulse

        # Fraction of the last pulse needed to reached max_dose
        last_frac = round((self.max_dose % self.dose_per_pulse) /
                          self.dose_per_pulse, 6)

        # if n == 0, max dose reached in less than a single pulse:
        if n <= 0.0:
            return last_frac * self.timer.on_time

        # From now on n >= 1.0 i.e. max_dose delivered in more than one pulse.
        # if max_dose delivered with entire periods then last one only last 'on_time'
        if last_frac <= 0.0:
            return ((n-1)*self.timer.period + self.timer.on_time)

        # From now on n >= 1.0 AND 0.0 < n_remainder < 1.0
        return (n * self.timer.period + last_frac * self.timer.on_time)

    def _compute_peak_dose_rate(self) -> float:
        """ Compute the dose rate during ON_time

        Returns:
            float: peek dose rate in Gy/s
        """
        return (self.dose_rate / self.timer.duty_cycle())

    def average_dose_rate(self) -> float:
        """ Method to make it clear we are dealing with average dose rate.
        Returns:
            float: average dose rate in Gy/s
        """
        return self.dose_rate

    def at(self, time: float) -> TimeMessage:
        """ 
        Computes instant dose rate at a given time.
        Args:
            time: float
                Time, in seconds, at which one want to know the dose rate in Gy/s

        Returns:
            TimeMessage
        """
        dr = 0.0
        if self.timer.start_time <= time < self.time_at_max_dose:
            match self.timer.state_at(time):
                case TimerState.isOff:
                    dr = 0.0
                case TimerState.isOn:
                    dr = self.peak_dose_rate
        return TimeMessage(
            time=time,
            current_dose_rate=dr,
        )

    def make_beam_func(self) -> tp.Callable[[float], float]:
        """
        USELESS FOR NOW
        => TODO ? or to clean?
        """
        def beam_func(time: float) -> float:
            if (time >= self.timer.start_time and
                    ((time-self.timer.start_time) % self.timer.period)
                    <= self.timer.on_time):
                return self.peak_dose_rate
            return 0.0
        return beam_func

    def make_matrix(self, stop: float, start: float = 0.0) -> NDArray[np.float64]:
        """ Method computing the radiation source time structure in the range [start,
        stop] i.e. computes 2 columns matrix whose columns are:
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

        This matrix definition of the beam structure is useful to predefine integration
        steps width (i.e. used by the "make_time_intervals" method).
        No default implementation provided.

        Args:
            stop: float
                End time of the beam structure, in seconds
            start: float, optional
                Start time of the beam structure, in seconds (default: 0.0)

        Returns:
            NDArray[np.float64]: 2 cols matrix
        """
        raise NotImplementedError

    def make_time_intervals(self,
                            stop: float,
                            start: float = 0.0,
                            ) -> NDArray[np.float64]:
        """
        Create a 2D matrix similar to the one produced by the make_matrix(stop, start)
        method. The major differences are:
            - takes into account max_dose in order to trim the matrix if needed.
            - Symmetrically expand the 'on_time' so that it is possible to define a 
              refined ODE integration time step for beam OFF-ON & ON-OFF transitions.
            - The first column is still time but the second is:
                * 0 for beam OFF
                * 1 for beam ON

        Returns:
            NDArray[np.float64]: 2 cols matrix
        """
        # Retrieve Time Structure Matrix
        out = self.make_matrix(stop, start)

        # If 3 lines ==> line 0 is for ON, line 1 is for OFF, line 2 is for END TIME
        # ==> Only need to grow time interval def by line 1 - line 0
        if out.shape[0] == 3:
            out[0, 1] = 1.0
            out[1, 0] = out[1, 0] + (out[1, 0]-out[0, 0])*self.GROWTH_FACTOR
            if out[1, 0] >= out[2, 0]:
                out = out[0:2, :]
            return out

        # if 2 lines ==> line 0 is for ON, line 1 is for OFF AND End time.
        # ==> only ON, no need to expand anything.
        if out.shape[0] == 2:
            return out

        # From now on, at least 2 ON times.
        # Define the ideal growth factor of on_time (default is 10% -> 5% each sides)
        delta = self.timer.on_time*self.GROWTH_FACTOR

        # Check whether it is no bigger than on period.
        if delta >= (self.timer.period - self.timer.on_time):
            # If so, refine the computed delta
            # TO DO? OK? Or return a 2 lines matrix?
            delta = (self.timer.period - self.timer.on_time)*self.GROWTH_FACTOR

        # 1) Shift the pulse start to "pulse start"  - (delta/2)
        # -------------------------------------------------------
        # select indices of the lines to modify (where dr > 0.0)
        # indexes must be a tuple of len 2 as matrix of 2 cols.
        # Want to modify time -> col 0
        indexes = np.where(out[:, 1] > 0.0)[0], [0]

        # IF first item is index 0, remove it, it is the start point.
        if indexes[0][0] == 0:
            indexes = (indexes[0][1:], indexes[1])
        # If last item is corresponding to end time, do not modify it.
        if indexes[0][-1] == out.shape[0]-1:
            indexes = (indexes[0][0:-1], indexes[1])
        # Shift
        out[indexes] = out[indexes] - (delta/2)

        # 2) Shift the pulse end to "pulse end"  + (delta/2)
        # -------------------------------------------------------
        indexes = np.where(out[:, 1] == 0.0)[0], [0]
        if indexes[0][0] == 0:
            indexes = (indexes[0][1:], indexes[1])
        if indexes[0][-1] == out.shape[0]-1:
            indexes = (indexes[0][0:-1], indexes[1])
        out[indexes] = out[indexes] + (delta/2)

        # Set all Dose rate values (-> col 1) > 0 to 1
        out[:, 1] /= self.peak_dose_rate
        return out

    def as_dict(self) -> tp.Dict:
        """ 
        """
        raise NotImplementedError


class ConstantBeam(ParticleBeam):
    """ Class representing a Constant Radiation Source
    Attributes:
        dose_rate: float
            Average dose rate in Gy/s
        max_dose: float
            Max dose to deliver in Gy
    """

    def __init__(self,
                 dose_rate: float,
                 max_dose: tp.Optional[float] = None,
                 start_time: tp.Optional[float] = None,
                 ) -> None:
        super().__init__(dose_rate,
                         Timer.new_constant(start_time),
                         max_dose=max_dose if max_dose is not None else float(
                             "inf"),
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
        return "ConstantBeam"

    def make_matrix(self, stop: float, start: float = 0.0) -> NDArray[np.float64]:
        out = np.empty([0, 2])

        # Always OFF cases:
        if stop < self.timer.start_time or start > self.time_at_max_dose:
            return np.array([[start, 0.0], [stop, 0.0]])

        # Always ON cases:
        if start >= self.timer.start_time and stop < self.time_at_max_dose:
            dr = self.peak_dose_rate
            return np.array([[start, dr], [stop, dr]])

        # 2 possibilities for start:
        if start >= self.timer.start_time:
            out = np.vstack([out, [start, self.peak_dose_rate]])
        else:
            out = np.vstack([out,
                             [start, 0.0],
                             [self.timer.start_time, self.peak_dose_rate],
                             ])

        # 2 possibilities for stop:
        if stop < self.time_at_max_dose:
            return np.vstack([out, [stop, self.peak_dose_rate]])

        return np.vstack([out,
                          [self.time_at_max_dose, 0.0],
                          [stop, 0.0],
                          ])

    def as_dict(self) -> tp.Dict[str, float]:
        """ """
        return {"dose_rate": self.dose_rate,
                "max_dose": self.max_dose,
                "start_time": self.timer.start_time,
                }

    @classmethod
    def from_dict(cls, param: tp.Dict[str, float]) -> ConstantBeam:
        """Build a ConstantBeam out of a dictionary"""
        dose_rate = 0.0
        for key in ["dose_rate_peak", "dose_rate_avg", "dose_rate"]:
            with suppress(KeyError):
                dose_rate = param[key]
                break

        return ConstantBeam(dose_rate=dose_rate,
                            max_dose=param.get("max_dose", None),
                            start_time=param.get("start_time", None))


class SinglePulseBeam(ConstantBeam):
    """Class representing a Single Pulse Radiation Source
    Attributes:
        dose_rate: float
            Average dose rate in Gy/s
        max_dose: float
            Max dose to deliver in Gy
        start_time: float, optional
            Time to start beam ON (default: 0.0)

    TODO: remove this class as it just mimic ConstantBeam? /!\ useful constructor from 
          'on_time' though
    """

    def __init__(self,
                 dose_rate: float,
                 max_dose: float,
                 start_time: tp.Optional[float] = None,
                 ) -> None:
        super().__init__(dose_rate,
                         max_dose=max_dose,
                         start_time=start_time,
                         )

    def __str__(self) -> str:
        return (f"SinglePulseBeam(dose_rate={self.dose_rate:.2g}, "
                f"max_dose={self.max_dose}, "
                f"start_time={self.timer.start_time})")

    def beam_type(self) -> str:
        return "SinglePulseBeam"

    @classmethod
    def from_on_time(cls,
                     dose_rate: float,
                     on_time: float,
                     start_time: tp.Optional[float] = None,
                     ) -> SinglePulseBeam:
        return cls(dose_rate, max_dose=on_time*dose_rate, start_time=start_time)

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
            raise ParameterError("SinglePulseBeam must be defined with a 'dose_rate' "
                                 "parameter")

        if "max_dose" in param:
            return SinglePulseBeam(dose_rate=param["dose_rate"],
                                   max_dose=param["max_dose"],
                                   start_time=param.get("start_time", None))
        if "on_time" in param:
            return SinglePulseBeam.from_on_time(dose_rate=param["dose_rate"],
                                                on_time=param["on_time"],
                                                start_time=param.get(
                                                    "start_time", None)
                                                )

        raise ParameterError("SinglePulseBeam must be defined with either one of "
                             "'on_time' or 'max_dose' parameter")


class PulsedBeam(ParticleBeam):
    """ Class representing a Pulsed Radiation Source. Helper functions are also 
    available to easily create specific Pulsed Radiation Sources:
        - from_peak_dose_rate
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

    def __init__(self,
                 dose_rate: float,
                 period: float,
                 on_time: float,
                 max_dose: tp.Optional[float] = None,
                 start_time: tp.Optional[float] = None,
                 ):
        super().__init__(dose_rate,
                         Timer.new_pulsed(period, on_time, start_time),
                         max_dose if max_dose is not None else float("inf"),
                         )

    # @override waiting for this baby to show up (Python3.12 fall 2023!)
    # => https://peps.python.org/pep-0698/
    def update_beam_param(self, **kwargs: float) -> None:
        if "peak_dose_rate" in kwargs:
            self = PulsedBeam.from_peak_dose_rate(
                peak_dr=kwargs["peak_dose_rate"],
                period=kwargs.pop("period", self.timer.period),
                on_time=kwargs.pop("on_time", self.timer.on_time),
                n_pulse=kwargs.pop("n_pulse", self.tot_pulse()),
                max_dose=kwargs.pop("max_dose", self.max_dose),
                start_time=kwargs.pop("start_time", self.timer.start_time),
            )
            return

        if "dose_per_pulse" in kwargs:
            self = PulsedBeam.from_dose_per_pulse(
                pulse_dose=kwargs["dose_per_pulse"],
                period=kwargs.pop("period", self.timer.period),
                on_time=kwargs.pop("on_time", self.timer.on_time),
                n_pulse=kwargs.pop("n_pulse", self.tot_pulse()),
                max_dose=kwargs.pop("max_dose", self.max_dose),
                start_time=kwargs.pop("start_time", self.timer.start_time),
            )
            return

        if "start_time" in kwargs or "on_time" in kwargs or "period" in kwargs:
            start_time = kwargs.pop("start_time", self.timer.start_time)
            on_time = kwargs.pop("on_time", self.timer.on_time)
            period = kwargs.pop("period", self.timer.period)
            self.timer = Timer.new_pulsed(period, on_time, start_time)

        super().update_beam_param(**kwargs)

    def beam_type(self) -> str:
        return "PulsedBeam"

    def tot_pulse(self) -> float:
        """Computes the total number of pulses.
        Returns:
            float
                Number of pulses (.2 digits precision)
        """
        fp = self.time_at_max_dose // self.timer.period
        out = fp + ((self.time_at_max_dose-fp *
                    self.timer.period) / self.timer.on_time)
        return round(out, 2)

    def make_matrix(self,
                    stop: float,
                    start: float = 0.0,
                    ) -> NDArray[np.float64]:
        if start > self.time_at_max_dose or stop < self.timer.start_time:
            return np.array([[start, 0.0], [stop, 0.0]], dtype=np.float64)

        out = np.empty([0, 2])
        if start < self.timer.start_time:
            out = np.vstack([out,
                             [start, 0.0,],
                             [self.timer.start_time, self.peak_dose_rate],
                             ])
        else:
            if self.timer.state_at(start) is TimerState.isOn:
                out = np.vstack([out,
                                [start, self.peak_dose_rate],
                                [start +
                                    self.timer.time_till_next_OFF(start), 0.0],
                                [start+self.timer.time_till_next_period(start),
                                 self.peak_dose_rate],
                                 ])
            else:
                out = np.vstack([out,
                                [start, 0.0],
                                [start+self.timer.time_till_next_period(start),
                                 self.peak_dose_rate],
                                 ])

        end = min(stop, self.time_at_max_dose)
        while out[-1, 0] < end:
            out = np.vstack([out,
                             [out[-1, 0]+self.timer.on_time, 0.0],
                             [out[-1, 0]+self.timer.period, self.peak_dose_rate],
                             ])

        out = out[:-1, :]
        out[-1, 0] = end

        if end < stop:
            out = np.vstack([out, [stop, 0.0]])
        return out

    def as_dict(self) -> tp.Dict[str, float]:
        """ """
        return {"dose_rate": self.dose_rate,
                "period": self.timer.period,
                "on_time": self.timer.on_time,
                "max_dose": self.max_dose,
                "start_time": self.timer.start_time,
                }

    @ classmethod
    def from_peak_dose_rate(cls,
                            peak_dr: float,
                            period: float,
                            on_time: float,
                            n_pulse: tp.Optional[float] = None,
                            max_dose: tp.Optional[float] = None,
                            start_time: tp.Optional[float] = None,
                            ) -> PulsedBeam:
        """ Class method building a PulsedBeam based on peak dose rate => dose rate
        during ON time.
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
        timer = Timer.new_pulsed(period, on_time,)

        n_pulse = n_pulse if n_pulse is not None else float("inf")
        max_dose_n_pulse = on_time*n_pulse*peak_dr
        max_dose = max_dose if max_dose is not None else float("inf")
        return cls(
            dose_rate=timer.duty_cycle() * peak_dr,
            period=period,
            on_time=on_time,
            max_dose=min(max_dose, max_dose_n_pulse),
            start_time=start_time,
        )

    @ classmethod
    def from_dose_per_pulse(cls,
                            pulse_dose: float,
                            period: float,
                            on_time: float,
                            n_pulse: tp.Optional[float] = None,
                            max_dose: tp.Optional[float] = None,
                            start_time: tp.Optional[float] = None,
                            ) -> PulsedBeam:
        """ Class method building a PulsedBeam based on a dose per pulse
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
        return PulsedBeam.from_peak_dose_rate(
            peak_dr=pulse_dose / on_time,
            period=period,
            on_time=on_time,
            n_pulse=n_pulse,
            max_dose=max_dose,
            start_time=start_time,
        )

    @classmethod
    def from_dict(cls, param: tp.Dict[str, float]) -> PulsedBeam:
        """ """
        if "period" not in param and "frequency" not in param:
            raise ParameterError(
                "PulsedBeam must be defined with a 'period' or 'frequency' parameter")

        if "frequency" in param:
            param["period"] = 1 / param["frequency"]

        if "on_time" not in param:
            if "dose_rate_peak" in param and "dose_rate" in param:
                param["on_time"] = param["dose_rate"] * param["period"]
                param["on_time"] /= param["dose_rate_peak"]
            else:
                raise ParameterError("PulsedBeam must be defined with an "
                                     "'on_time' parameter")

        if "dose_rate_peak" in param:
            return PulsedBeam.from_peak_dose_rate(
                param["dose_rate_peak"],
                period=param["period"],
                on_time=param["on_time"],
                n_pulse=param.get("n_pulse", None),
                max_dose=param.get("max_dose", None),
                start_time=param.get("start_time", None),
            )
        if "pulse_dose" in param:
            return PulsedBeam.from_dose_per_pulse(
                param["pulse_dose"],
                period=param["period"],
                on_time=param["on_time"],
                n_pulse=param.get("n_pulse", None),
                max_dose=param.get("max_dose", None),
                start_time=param.get("start_time", None),
            )

        if "dose_rate_avg" in param or "dose_rate" in param:
            return PulsedBeam(
                dose_rate=param.get(
                    "dose_rate_avg", param.get("dose_rate", 0.0)),
                period=param["period"],
                on_time=param["on_time"],
                max_dose=param.get("max_dose", None),
                start_time=param.get("start_time", None),
            )
        raise NotImplementedError(f"Unknown Pulsed beam definition: {param}")


class BeamCollection(object):
    """
    Container class for efficient beam storage so that multiple beams can be define
    (even though only one at a time can be use). 
    """
    # Beam definitions are stored in this dict
    _beam: tp.Dict[str, ParticleBeam]
    # Keyword to access the currently used beam from _beam
    _current: str
    # LET to use with the corresponding beams
    _let: tp.Dict[str, float | None]

    def __init__(self) -> None:
        self._beam = {}
        self._let = {}
        self._current = ""

    def __contains__(self, key: str) -> bool:
        return key in self._beam

    def __iter__(self) -> tp.Iterator[tp.Tuple[str, ParticleBeam, float | None]]:
        for name, beam in self._beam.items():
            yield name, beam, self._let[name]

    def __len__(self) -> int:
        return len(self._beam)

    def __getitem__(self, key: str) -> tp.Tuple[ParticleBeam, float | None]:
        return self._beam[key], self._let[key]

    def __str__(self) -> str:
        out = ["BeamCollection:\n"]
        for name, beam, let in self:
            out.append(f"Beam {name} [{let}keV/µm] -> {beam}\n")
        return "".join(out)

    @property
    def current(self) -> ParticleBeam:
        return self._beam[self._current]

    @property
    def current_name(self) -> str:
        return self._current

    @property
    def current_LET(self) -> float:
        if self._let[self._current] is not None:
            return self._let[self._current]  # type: ignore [return-value]
        raise ValueError(f"No LET defined for beam: {self._current}")

    def remove_beam(self, key: str) -> None:
        """Remove a beam entry.
        Args:
            key: str
                name of the beam to delete
        Raises:
            KeyError
                No beam with name: 'key'
        """
        if key not in self:
            raise KeyError(f"No beam definition named: {key}")
        del self._let[key]
        del self._beam[key]

    def update_let(self, let: float, key: tp.Optional[str] = None) -> None:
        """Modifies the LET of a beam
        Args:
            let: float
                LET value
            key: str, optional (default=None)
                name of the beam to update. Default is the current one.
        Raises:
            KeyError
                Beam with name 'key' not found
        """
        if key is None:
            key = self._current
        elif key not in self:
            raise KeyError(f"No beam definition named: {key}")
        self._let[key] = let

    def add_beam(self,
                 beam: ParticleBeam,
                 key: str = "default",
                 let: tp.Optional[float] = None,
                 use: bool = False) -> None:
        """
        Pushes a new beam into the container.
        Args:
            beam: ParticleBeam
                beam to push in
            key: str, optional (default="default")
                Name of the beam
            let: float, optional (default=None)
                LET of the beam 
            use: bool, optional (default=False)
                Defines the newly pushed beam as the current one
        """
        if key is None or len(key) == 0:
            key = "default"
        # First remove any previous definition with the same 'key'
        if key in self:
            self.remove_beam(key)
        self._beam[key] = beam
        self._let[key] = let
        if use or len(self._current) == 0:
            self._current = key

    def use_beam(self, key: str) -> tp.Tuple[ParticleBeam, float]:
        """
        Defines which beam to use
        Args:
            key: str
                Beam name
        Raises:
            KeyError
                Beam with name 'key' not found
        Returns:
            ParticleBeam
            float
                The beam LET
        """
        if key not in self._beam:
            raise KeyError(f"No beam definition named: {key}")
        self._current = key
        return self.current, self.current_LET

    def add_beam_from_dict(self,
                           beam_type: str,
                           config: tp.Dict[str, float],
                           LET: tp.Optional[float] = None,
                           verbose: bool = True
                           ) -> None:
        """
        Add new beam defined in a dictionary.
        Args:
            beam_type: str
                type of the beam to add
            config: dict
                Dictionary containing values defining the beam
            LET: float, optional
                LET of the beam
            verbose: bool, optional (default=True)
                Enable or disable logging

        Raises:
            NotImplementedError
                beam_type not recognized
        """
        # if LET not provided, check if present in config.
        if LET is None and "LET" in config:
            LET = config["LET"]

        if LET is not None and "LET" not in config:
            config["LET"] = LET  # So that LET is displayed in the bottom log

        label: str
        if "label" not in config:
            msg = ("BEAM:: Assigning 'default' name to beam definition "
                   f"{beam_type} -> {config}")
            if verbose:
                logger.info(msg)
            label = "default"
        else:
            label = str(config["label"])
            msg = (f"BEAM:: Found beam '{label}' as {beam_type} -> {config}")
            if verbose:
                logger.info(msg)

        beam: ParticleBeam
        if beam_type.lower() == "constant":
            beam = ConstantBeam.from_dict(config)
        elif beam_type.lower() == "pulsed":
            beam = PulsedBeam.from_dict(config)
        elif beam_type.lower() == "singlepulse":
            beam = SinglePulseBeam.from_dict(config)
        else:
            raise NotImplementedError(f"Unknown beam type: {beam_type}")

        if verbose and label in self:
            msg = (f"BEAM:: Beam definition: {label} already defined. "
                   "Erasing previous definition")
            logger.info(msg)

        self.add_beam(beam=beam, key=label, let=LET)

    @classmethod
    def from_config_dict(cls,
                         param: tp.Dict[str, tp.List[tp.Dict]],
                         verbose: bool = True,
                         ) -> BeamCollection:
        """
        Create a BeamCollection out of a Dict containing multiple beam definitions.
        param dictionary keywords should be valid beam type names. Values can be either
        be a single Dict-like beam or a list of it.
        Args:
            param: Dict
                Dict structure whose keywords are valid beam type names.
            verbose: bool, optional (default=True)
        """
        beams = cls()

        for beam_type, beam_config in param.items():
            # If only one beam, not in a list, put it in a list!
            if isinstance(beam_config, dict):
                beam_config = [beam_config, ]
            for config in beam_config:
                LET: tp.Optional[float] = None
                with suppress(KeyError):
                    LET = config.pop("LET")

                beams.add_beam_from_dict(beam_type, config, LET, verbose)
        return beams

    def as_dict(self) -> tp.Dict[str, tp.List[tp.Dict]]:
        """
        Extract all the beam in a dictionary whose structure comply with the method
        'from_config_dict'
        """
        out: tp.Dict[str, tp.List[tp.Dict]] = {}
        for name, beam, LET in self:
            x = beam.as_dict()
            if name != "default":
                x["label"] = name
            if LET is not None:
                x["LET"] = self._let[name]
            beam_type = beam.beam_type().lower().split("beam")[0]
            if beam_type not in out:
                out[beam_type] = [x, ]
            else:
                out[beam_type].append(x)
        return out

    def merge(self, other: BeamCollection, verbose: bool = True) -> BeamCollection:
        """
        Merge two BeamCollection. Do not modify the current definition but return a new
        one.

        Args:
            other: BeamCollection
                collection to merge with
            verbose: bool, optional (default=True)
                Flag for logging

        Returns:
            BeamCollection
        """
        s: tp.Dict[str, tp.Tuple[tp.Dict, str]] = {}
        for name, _beam, let in self:
            d = _beam.as_dict()
            if let is not None:
                d["LET"] = let
            d["label"] = name
            # Bug Fix: important to use plit('beam') and not strip('beam')
            s[name] = (d, _beam.beam_type().lower().split("beam")[0])

        for beam_type, beams in other.as_dict().items():
            for beam in beams:
                name = beam.get("label", "default")
                if name not in s:
                    s[name] = (beam, beam_type)
                    continue
                if s[name][1] == beam_type:
                    s[name][0] |= beam
                    continue
                s[name] = (beam, beam_type)

        out = BeamCollection()
        for _, (beam, beam_type) in s.items():
            out.add_beam_from_dict(beam_type,
                                   beam,
                                   verbose=verbose
                                   )

        return out


def beam_from_dict(param: tp.Dict[str, tp.Dict[str, float]],
                   ) -> tp.Tuple[ParticleBeam, tp.Optional[float]]:
    """ 
    Construct beam from a dictionary.

    Parameters:
        param: Dict[str, Dict[str, float]]
            Parameters defining the beam to construct. The first key define the beam
            type.

    Returns:
        ParticleBeam:
            Beam class.
        float:
            LET.

    """
    if len(param) > 1:
        raise ParameterError(f"Found more than one beam definition ({param})")
    beam_type, config = param.popitem()

    LET: tp.Optional[float] = None
    with suppress(KeyError):
        LET = config.pop("LET")

    if beam_type.lower() == "constant":
        return ConstantBeam.from_dict(config), LET

    if beam_type.lower() == "pulsed":
        return PulsedBeam.from_dict(config), LET

    if beam_type.lower() == "singlepulse":
        return SinglePulseBeam.from_dict(config), LET

    raise NotImplementedError(f"Unknown beam type: {beam_type}")
