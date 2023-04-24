#!/usr/bin/env python3
"""
The beam submodule contains all class definitions related to radiation source
structures.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023 
"""
# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from math import isinf, isnan
from typing import (
    TYPE_CHECKING,
    Callable,
    NamedTuple,
)

import numpy as np

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #
if TYPE_CHECKING:
    from numpy.typing import NDArray

# ----------------------------- Expose importable stuffs ----------------------------- #
__all__ = ["TimeMessage",
           "TimerState",
           "Timer",
           "ParticleBeam",
           "ConstantBeam",
           "PulsedBeam",
           ]

# ------------------------------------------------------------------------------------ #
#                                  CLASS DEFINITIONS                                   #
# ------------------------------------------------------------------------------------ #


class TimerState(Enum):
    """ Convenient Enum use to define radiation source ON or OFF state."""
    isOff = 0
    isOn = 1


class TimeMessage(NamedTuple):
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


class Timer(NamedTuple):
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

    @classmethod
    def new_constant(cls) -> Timer:
        """ Class method to generate an always ON timer"""
        return cls(float("inf"), float("inf"))

    @classmethod
    def new_pulsed(cls, period: float, on_time: float) -> Timer:
        """ Class method to generate a pulsed timer
        Args:
            period: float
                Period (repetition) of the timer in seconds
            on_time: float
                Time during which timer is ON in seconds. on_time should be <= period.
        """
        return cls(period, on_time)

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
        if (time % self.period) <= self.on_time:
            return TimerState.isOn
        return TimerState.isOff


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
    dose_rate: float
    timer: Timer
    max_dose: float = float("inf")
    time_at_max_dose: float = field(init=False)
    dose_per_pulse: float = field(init=False)
    peak_dose_rate: float = field(init=False)
    GROWTH_FACTOR: float = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """ Method called after ParticleBeam__init__() to finalize initialization"""
        self.GROWTH_FACTOR = 0.1
        self.peak_dose_rate = self._compute_peak_dose_rate()
        self.dose_per_pulse = (self.peak_dose_rate * self.timer.on_time)
        self.time_at_max_dose = self._compute_time_at_max_dose()

    def _compute_time_at_max_dose(self) -> float:
        """Helper method to compute the time at which max dose rate will be reached with
        the average dose rate and Timer defined.

        Returns:
            float: time in seconds.
        """
        # if max_dose == inf, max_dose is never reached
        if self.max_dose >= float("inf"):
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
        last_frac = (self.max_dose % self.dose_per_pulse) / \
            self.dose_per_pulse

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
        if time < self.time_at_max_dose:
            match self.timer.state_at(time):
                case TimerState.isOff:
                    dr = 0.0
                case TimerState.isOn:
                    dr = self.peak_dose_rate
        return TimeMessage(
            time=time,
            current_dose_rate=dr,
        )

    def make_beam_func(self) -> Callable[[float], float]:
        """
        USELESS FOR NOW
        => TODO ? or to clean?
        """
        def beam_func(time: float) -> float:
            if (time % self.timer.period) <= self.timer.on_time:
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
            - increases the the on_time so that it is possible to define a refined ODE
              integration step for beam OFF-ON & ON-OFF transitions.
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
        # ==> only ON, no need to grow anything.
        if out.shape[0] == 2:
            return out

        # From now on, at least 2 ON times.
        # Define the ideal grow factor of on_time (default is 10% -> 5% each sides)
        delta = self.timer.on_time*self.GROWTH_FACTOR

        # Check whether it is no bigger than on period.
        if delta >= (self.timer.period - self.timer.on_time):
            # If so, refine the computed delta
            # TO DO? OK? Or return a 2 lines matrix?
            delta = (self.timer.period - self.timer.on_time)*self.GROWTH_FACTOR

        # 1) Shift the pulse start to "pulse start"  - (delta/2)
        # -------------------------------------------------------
        # select indices of the lines to modify (where dr > 0.0)
        # indexes mist be a tuple of len 2 as matrix of 2 cols.
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
                 max_dose: float = float("inf"),
                 ):
        super().__init__(dose_rate,
                         Timer.new_constant(),
                         max_dose,
                         )

    def make_matrix(self, stop: float, start: float = 0.0) -> NDArray[np.float64]:
        # If range after max dose was delivered.
        if start >= self.time_at_max_dose:
            return np.array([[start, 0.0],
                             [stop, 0.0]])
        # If always ON within [start, stop]
        if stop < self.time_at_max_dose:
            return np.array([[start, self.peak_dose_rate],
                             [stop, self.peak_dose_rate]])
        # If ON, then OFF after max dose is delivered
        return np.array([[start, self.peak_dose_rate],
                         [self.time_at_max_dose, 0.0],
                         [stop, 0.0]])


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
    """

    def __init__(self,
                 dose_rate: float,
                 period: float,
                 on_time: float,
                 max_dose: float = float("inf"),
                 ):
        super().__init__(dose_rate,
                         Timer.new_pulsed(period, on_time),
                         max_dose,
                         )

    def make_matrix(self, stop: float, start: float = 0.0) -> NDArray[np.float64]:
        # If sim starts after beam off:
        if start >= self.time_at_max_dose:
            return np.array([[start, 0.0], [stop, 0.0]], dtype=np.float64)

        # If sim ends in less than one on_time:
        if stop < self.timer.on_time:
            return np.array([[start, self.peak_dose_rate],
                             [stop, 0.0]],
                            dtype=np.float64)

        # If total dose is delivered in less than one pulse:
        if self.time_at_max_dose < self.timer.on_time:
            return np.array([[start, self.peak_dose_rate],
                             [self.time_at_max_dose, 0.0],
                             [stop, 0.0]],
                            dtype=np.float64)

        # Else ...
        start_period: int = int(start//self.timer.period)
        stop_period: int
        if stop > self.time_at_max_dose:
            stop_period = int(self.time_at_max_dose//self.timer.period)+1
        else:
            stop_period = int(stop//self.timer.period)+1
        times = np.empty([stop_period-start_period, 2],
                         dtype=np.float64)
        times[:, 0] = np.arange(start_period, stop_period) * self.timer.period
        times[:, 1] = times[:, 0] + self.timer.on_time

        dr = np.zeros_like(times)
        indexes = (np.where(times[:, 0] < self.time_at_max_dose)[0], [0])
        dr[indexes] = self.peak_dose_rate

        out = np.empty([dr.shape[0]*2, 2])
        out[:, 0] = times.flatten()
        out[:, 1] = dr.flatten()

        # Trim based on start/stop times:
        if out[1, 0] > start:
            out[0, 0] = start
        if out[-1, 0] < stop:
            out = np.vstack([out, [stop, 0.0]])

        return out

    @ classmethod
    def from_peak_dose_rate(cls,
                            peak_dr: float,
                            period: float,
                            on_time: float,
                            max_dose: float = float("inf")
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
            max_dose: float, optional
                Max dose to deliver, in Gy (default is inf)

        Returns:
            PulsedBeam
        """
        timer = Timer.new_pulsed(period, on_time)
        return cls(dose_rate=timer.duty_cycle() * peak_dr,
                   period=period,
                   on_time=on_time,
                   max_dose=max_dose,
                   )

    @ classmethod
    def from_dose_per_pulse(cls,
                            pulse_dose: float,
                            period: float,
                            on_time: float,
                            max_dose: float = float("inf")
                            ) -> PulsedBeam:
        """ Class method building a PulsedBeam based on a dose per pulse
        Args:
            pulse_dose: float
                dose delivered during ONE single pulse
            period: float
                Pulse period, in seconds
            on_time: float
                Pulse ON time, in seconds
            max_dose: float, optional
                Max dose to deliver, in Gy (default is inf)

        Returns:
            PulsedBeam
        """
        return PulsedBeam.from_peak_dose_rate(peak_dr=pulse_dose / on_time,
                                              period=period,
                                              on_time=on_time,
                                              max_dose=max_dose,)


class MonoPulsedBeam(PulsedBeam):
    """ 
    TODO
    """
