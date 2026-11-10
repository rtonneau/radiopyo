#!/usr/bin/env python3
"""
The timer submodule contains the time structure (ON/OFF state machines)
of the radiation sources.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 16/03/2023
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import logging
import math
import typing as tp
from abc import ABC, abstractmethod
from enum import Enum

# ----------------------------- Internal imports ---------------------------- #
# ------------------------- Type checking ONLY imports----------------------- #
# -------------------------- Expose importable stuffs ----------------------- #
__all__ = [
    "BaseTimer",
    "MacroTimer",
    "Timer",
    "TimerState",
]

# ---------------------------------- Logging -------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
# Relative tolerance used when counting how many periods fit in a duration
# (0.001 / 0.0001 is not exactly 10 in floating point).
COUNT_RTOL = 1e-9


def split_periods(duration: float, period: float) -> tp.Tuple[int, float]:
    """
    Split a duration in a whole number of periods and a remainder, with
    a relative tolerance on the period.

    Args:
        duration: float
            Duration to split, in seconds.
        period: float
            Period, in seconds.

    Returns:
        Tuple[int, float]: (number of full periods, remainder >= 0).
        A remainder within COUNT_RTOL * period of 0 is returned as 0.
    """
    n = math.floor(duration / period + COUNT_RTOL)
    remainder = duration - n * period
    if remainder <= COUNT_RTOL * period:
        remainder = 0.0
    return n, remainder


# --------------------------------------------------------------------------- #
#                              CLASS DEFINITIONS                              #
# --------------------------------------------------------------------------- #


class TimerState(Enum):
    """
    Convenient Enum use to define radiation source states, ON or OFF.
    """

    isOff = 0
    isOn = 1


class BaseTimer(ABC):
    """
    Abstract base class defining a radiation source timer.

    This class is designed to be inherited by concrete radiation source
    timer classes. It implements the basic logic to manage the state of
    a radiation source and provides a method to get the dose rate at a
    specific time.

    Attributes:
        period (float): The time (in seconds) between two consecutive
            ON periods.
        on_time (float): The duration (in seconds) of the ON period.
        start_time (float): The time (in seconds) at which the radiation
            source starts.

    """

    period: float
    on_time: float
    start_time: float

    def __init__(
        self, period: float, on_time: float, start_time: tp.Optional[float] = None
    ) -> None:
        self.period = period
        self.on_time = on_time
        self.start_time = (
            start_time if start_time is not None else self._default_start_time()
        )

    def __str__(self) -> str:
        return f"BaseTimer({self.period}, {self.on_time}, {self.start_time})"

    @staticmethod
    def _default_start_time() -> float:
        """
        Returns the default start time for a radiation source timer.

        This static method returns the default start time for a
        radiation source timer. The default start time is 0.0 seconds.

        Returns:
            float: The default start time in seconds.
        """
        return 0.0

    @abstractmethod
    def state_at(self, time: float) -> TimerState:
        """
        Computes the state of the timer at a given time.

        The ON window of a period is closed on both sides: a time
        exactly at the end of a pulse (``on_time`` after the pulse
        start) is still ON.

        Parameters:
            time (float): The time in seconds.

        Returns:
            TimerState: An enumeration representing the state of the
            timer. Possible values are:
                - isOn: The timer is on.
                - isOff: The timer is off.
        """

    @abstractmethod
    def time_till_next_ON(self, time: float) -> float:
        """
        Computes the time remaining before the next pulse of the timer.

        This method calculates and returns the time remaining before the
        next pulse of the timer. It takes a single parameter, `time`,
        which represents the current time in seconds. The method returns
        a float representing the time remaining in seconds.

        Returns:
            float: The time remaining before the next pulse in seconds.
        """

    @abstractmethod
    def time_till_next_OFF(self, time: float) -> float:
        """
        Computes the time remaining before the next pulse of the timer.

        This method calculates and returns the time remaining before the
        next pulse of the timer. If the timer is currently on, it returns
        the time remaining in the current pulse. If the timer is off, it
        returns 0.0. The method takes a single parameter, `time`, which
        represents the current time in seconds. The method returns a float
        representing the time remaining in seconds.

        Returns:
            float: Remaining time pulse in seconds. If the timer is off,
            it returns 0.0.
        """

    @abstractmethod
    def copy(self, start_time: tp.Optional[float] = None) -> tp.Self:
        """
        Create a copy of the timer.

        This method creates a copy of the timer and returns it. The copy
        is a new instance of the same class as the original timer. The
        start time of the copy can be optionally specified. If not
        provided, the start time of the copy will be the same as the
        original timer.

        Parameters:
            start_time (float, optional): The start time of the copy
                in seconds. If not provided, the start time of the copy
                will be the same as the original timer.

        Returns:
            Self: A new instance of the same class as the original timer.
        """

    def period_index(self, time: float) -> int:
        """
        Get the index of the period that contains a given time.

        Parameters:
            time (float): The time for which to find the period index.

        Returns:
            int: The index of the period that contains the given time.
                The period index is calculated by subtracting the start
                time from the given time and then dividing the result by
                the period duration. The integer division is used to
                determine the index.

        Example:
            For a period duration of 3 seconds and a start time of 0
            seconds, the period index of a time of 8 seconds would be 2.
        """
        # Calculate the time relative to the start time.
        time_relative = time - self.start_time

        # Calculate the period index by dividing the time relative to the
        # start time by the period duration and taking the integer part.
        return int(time_relative // self.period)

    def _phase(self, time: float) -> float:
        """
        Time elapsed since the beginning of the period containing `time`.

        No check is made on `time` being later than the timer start.
        """
        return time - self.start_time - self.period * self.period_index(time)

    def time_till_next_period(self, time: float) -> float:
        """
        Computes how much time is remaining before next period starts.

        Parameters:
            time (float): The current time in seconds.

        Returns:
            float: The remaining time before the next period starts in
                seconds.
        """
        return self.period - self._phase(time)

    def off_time(self) -> float:
        """
        Compute the length of the OFF time.

        This method calculates and returns the length of the OFF time
        of the timer. The OFF time is the duration of each period
        between the ON times of the timer. The method takes no
        parameters and returns a float representing the length of the
        OFF time in seconds.

        Returns:
            float: The length of the OFF time in seconds.
        """
        return self.period - self.on_time

    def _relative_time_pulse(self, time: float) -> float:
        """
        Calculate the relative time within a single pulse.

        This method calculates the relative time within a single pulse
        of the timer. The relative time is the time within the ON time
        of the timer. The method takes a time parameter and returns a
        float representing the relative time in seconds.

        Parameters:
            time (float): The time for which to calculate the relative time.

        Returns:
            float: The relative time within a single pulse in seconds.

        Raises:
            ValueError: If the time is earlier than the start time of
                the timer.

        Example:
            For a timer with a start time of 0 seconds, a period of 5
            seconds, and an ON time of 2 seconds, the relative time at
            time 4 seconds would be 2 seconds. The relative time at time
            6 seconds would be 1 second.
        """
        if time < self.start_time:
            raise ValueError("Time is earlier than timer.start_time")
        return self._phase(time)

    def duty_cycle(self) -> float:
        """
        Compute the duty cycle of the timer.

        The duty cycle is the ratio of the ON time to the period.
        If the timer is continuous, the duty cycle is 1.0.

        Returns:
            float: The duty cycle of the timer.

        Example:
            For a timer with a period of 5 seconds and an ON time of 2
            seconds, the duty cycle would be 0.4 (2/5).

            For a continuous timer, the duty cycle would be 1.0.
        """
        # If the timer is continuous, return 1.0.
        if self.is_continuous():
            return 1.0

        # Calculate the duty cycle.
        return self.on_time / self.period

    def is_continuous(self) -> bool:
        """Is timer continuous?
        Returns:
            bool
        """
        return (self.period == float("inf") and self.on_time == float("inf")) or (
            self.on_time >= self.period
        )

    def is_mono_pulse(self) -> bool:
        """Is the timer defined as a single pulse? i.e.
            - period = inf
            - 0 < on_time < period
        Returns:
            bool
        """
        return self.period == float("inf") and 0 < self.on_time < float("inf")


class Timer(BaseTimer):
    """Class used as a binary timer for radiation source. Currently only
    two states are possible for a radiation source, either ON or OFF
    i.e. 0.0% or 100%.
    his class holds the time structure of the radiation source,
    defined as pulses. Pulses are repeated with a given 'period'
    (seconds) and the timer keeps its ON state during 'on_time' seconds.
    For example:
        - period = 2 seconds
        - on_time = 0.5 second
    --> each 2 seconds, timer ON during 0.5 second
    (i.e. duty cycle of 25%)

    Timer should not be constructed from scratch but from helper
    functions (classmethod)
    => Either one of:
        - Timer.new_constant()
        - Timer.new_pulsed(period, on_time)

    Usage:
        >>> timer = Timer.new_pulsed(period=2, on_time=0.5)
        >>> timer.state_at(0.25)      # -> TimerState.isOn
        >>> timer.state_at(0.50)      # -> TimerState.isOn (pulse end is ON)
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

    @classmethod
    def new_constant(cls, start_time: tp.Optional[float] = None) -> tp.Self:
        """Class method to generate an always ON timer"""
        stime = start_time if start_time is not None else cls._default_start_time()
        return cls(float("inf"), float("inf"), stime)

    @classmethod
    def new_pulsed(
        cls, period: float, on_time: float, start_time: tp.Optional[float] = None
    ) -> tp.Self:
        """
        Class method to generate a pulsed timer

        Args:
            period: float
                Period (repetition) of the timer in seconds
            on_time: float
                Time during which timer is ON in seconds.
                on_time should be <= period.
        """
        stime = start_time if start_time is not None else cls._default_start_time()
        return cls(period, on_time, stime)

    def state_at(self, time: float) -> TimerState:
        """
        Compute the state of the timer at a given time.

        A time exactly at the end of a pulse (``start_time + k * period
        + on_time``) is still ON: the ON window is closed on both sides.

        Args:
            time: float
                Time in seconds
        Returns:
            TimerState: Enum being either isOn or isOff
        """
        if (time >= self.start_time) and (
            (time - self.start_time) % self.period
        ) <= self.on_time:
            return TimerState.isOn
        return TimerState.isOff

    def time_till_next_ON(self, time: float) -> float:
        """
        Computes time remaining before next pulse.
        """
        if self.state_at(time) is TimerState.isOn:
            return 0.0
        return self.time_till_next_period(time)

    def time_till_next_OFF(self, time: float) -> float:
        """
        Computes time remaining before beam is off.
        Returns 0 if beam already off.
        """
        if self.state_at(time) is TimerState.isOff:
            return 0.0
        return self.on_time - self._phase(time)

    def copy(self, start_time: tp.Optional[float] = None) -> Timer:
        """Make a simple copy of the Timer (same start time unless given)"""
        stime = start_time if start_time is not None else self.start_time
        return Timer(self.period, self.on_time, stime)


class MacroTimer(BaseTimer):
    """
    Macro timer enabling Micro/Macro Pulse.
    """

    subtimer: BaseTimer

    def __init__(
        self,
        period: float,
        on_time: float,
        subtimer: BaseTimer,
        start_time: tp.Optional[float] = None,
    ) -> None:
        self.subtimer = subtimer
        super().__init__(period, on_time, start_time)

    @classmethod
    def new_macropulse(
        cls,
        macro_period: float,
        macro_on_time: float,
        micro_period: float,
        micro_on_time: float,
        macro_start_time: tp.Optional[float] = None,
        micro_start_time: tp.Optional[float] = None,
    ) -> tp.Self:
        """Class method to generate a macro pulsed timer
        Args:
            macro_period: float
                Period (repetition) of the timer in seconds
            period: float
                Period (repetition) of the timer in seconds
            on_time: float
                Time during which timer is ON in seconds.
                on_time should be <= period.

        Raises:
            ValueError: if macro_period <= micro_period, or if
                micro_start_time >= macro_on_time (no micro pulse would
                fit in a train).
        """
        if macro_period <= micro_period:
            raise ValueError(
                "Macro Period is <= than micro period. "
                "Do you really need a macro timer?"
            )
        if micro_start_time is not None and micro_start_time >= macro_on_time:
            raise ValueError(
                f"micro_start_time ({micro_start_time}) must be < the train "
                f"on_time ({macro_on_time}): no micro pulse would fit in a train."
            )
        stime = (
            macro_start_time
            if macro_start_time is not None
            else cls._default_start_time()
        )
        return cls(
            macro_period,
            macro_on_time,
            Timer.new_pulsed(micro_period, micro_on_time, micro_start_time),
            stime,
        )

    def state_at(self, time: float) -> TimerState:
        """
        Compute the state of the timer at a given time.

        A time exactly at the end of a micro pulse (or of a train) is
        still ON: the ON windows are closed on both sides.

        Args:
            time: float
                Time in seconds
        Returns:
            TimerState: Enum being either isOn or isOff
        """
        if time < self.start_time:
            return TimerState.isOff
        time = self._relative_time_pulse(time)
        if time <= self.on_time:
            return self.subtimer.state_at(time)
        return TimerState.isOff

    def micro_pulses_per_train(self) -> float:
        """
        Number of micro pulses delivered during one train. Micro pulses
        start at `subtimer.start_time` after the train start; a pulse cut
        by the end of the train counts for the fraction of its on_time
        that is delivered (train 0.205 s, micro period 0.1 s, micro
        on_time 0.01 s gives 2.5).
        """
        available = self.on_time - self.subtimer.start_time
        if available <= 0:
            return 0.0
        n, remainder = split_periods(available, self.subtimer.period)
        out = float(n)
        if remainder > 0:
            out += min(remainder / self.subtimer.on_time, 1.0)
        return out

    def micro_period_index(self, time: float) -> int:
        """ """
        time = self._relative_time_pulse(time)
        if time > self.on_time:
            return -1
        return self.subtimer.period_index(time)

    def time_till_end_micro_pulse(self, time: float) -> float:
        """ """
        time = self._relative_time_pulse(time)
        if time < self.on_time:
            return self.on_time - time
        return 0.0

    def time_till_next_train(self, time: float) -> float:
        """
        Computes time remaining before next train starts.
        """
        if time < self.start_time:
            return self.start_time - time
        return self.period - self._relative_time_pulse(time)

    def time_till_next_ON(self, time: float) -> float:
        """
        Computes time remaining before next pulse (0.0 if the beam is
        ON at `time`). The next micro pulse is searched in the current
        train, then at the start of the next one.
        """
        if time < self.start_time:
            return self.start_time - time + self.subtimer.start_time
        if self.state_at(time) is TimerState.isOn:
            return 0.0
        rel = self._relative_time_pulse(time)
        if rel <= self.on_time:
            wait = self.subtimer.time_till_next_ON(rel)
            if rel + wait < self.on_time:
                return wait
        return self.period - rel + self.subtimer.start_time

    def time_till_next_OFF(self, time: float) -> float:
        """
        Computes time remaining before beam is off (the end of a train
        cuts the micro pulse). Returns 0 if beam already off.
        """
        if self.state_at(time) is TimerState.isOff:
            return 0.0
        rel = self._relative_time_pulse(time)
        return min(self.subtimer.time_till_next_OFF(rel), self.on_time - rel)

    def duty_cycle(self) -> float:
        """Computes the duty cycle of the MacroTimer.
        Rem:
            To avoid confusion for cases where the ON time is not a
            multiple of subperiods, the duty cycle is computed by
            calculating the total time ON during one macro pulse.
        Returns:
            float
        """
        return self.micro_pulses_per_train() * self.subtimer.on_time / self.period

    def copy(self, start_time: tp.Optional[float] = None) -> MacroTimer:
        """Make a simple copy of the Timer (same start time unless given)"""
        stime = start_time if start_time is not None else self.start_time
        return MacroTimer(self.period, self.on_time, self.subtimer.copy(), stime)
