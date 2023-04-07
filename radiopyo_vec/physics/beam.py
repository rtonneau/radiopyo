# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from math import isinf, isnan
from typing import (
    TYPE_CHECKING,
    Callable,
    List,
    NamedTuple,
    Tuple,
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
    isOff = 0
    isOn = 1


class TimeMessage(NamedTuple):
    time: float
    current_dose_rate: float

    def dose_rate(self) -> float:
        return self.current_dose_rate


class Timer(NamedTuple):
    period: float
    on_time: float

    @classmethod
    def new_constant(cls) -> Timer:
        return cls(float("inf"), float("inf"))

    @classmethod
    def new_pulsed(cls, period: float, on_time: float) -> Timer:
        return cls(period, on_time)

    def is_continuous(self) -> bool:
        return self.period == float("inf") and self.on_time == float("inf")

    def is_mono_pulse(self) -> bool:
        return self.period == float("inf") and 0 < self.on_time < float("inf")

    def duty_cycle(self) -> float:
        out = self.on_time / self.period
        # In case of "inf"/"inf" => nan
        if isnan(out) or isinf(out):
            return 1.0
        return out

    def state_at(self, time: float) -> TimerState:
        if (time % self.period) <= self.on_time:
            return TimerState.isOn
        return TimerState.isOff


@dataclass
class ParticleBeam():
    dose_rate: float
    timer: Timer
    max_dose: float = float("inf")
    time_at_max_dose: float = field(init=False)
    dose_per_pulse: float = field(init=False)
    peak_dose_rate: float = field(init=False)

    def __post_init__(self) -> None:
        self.peak_dose_rate = self._compute_peak_dose_rate()
        self.dose_per_pulse = (self.peak_dose_rate * self.timer.on_time)
        self.time_at_max_dose = self._compute_time_at_max_dose()

    def _compute_time_at_max_dose(self) -> float:
        if self.max_dose >= float("inf"):
            return float("inf")
        if self.timer.period >= float("inf"):
            return self.max_dose / self.peak_dose_rate
        # How many pulses to reach max_dose?
        #   ==> (max dose) / (dose per pulse)
        n = self.max_dose / self.dose_per_pulse
        return int(n) * self.timer.period + (n-int(n)) * self.timer.on_time

    def _compute_peak_dose_rate(self) -> float:
        return (self.dose_rate / self.timer.duty_cycle())

    def average_dose_rate(self) -> float:
        return self.dose_rate

    def at(self, time: float) -> TimeMessage:
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
        def beam_func(time: float) -> float:
            if (time % self.timer.period) <= self.timer.on_time:
                return self.peak_dose_rate
            return 0.0
        return beam_func

    def make_matrix(self, start: float, stop: float) -> NDArray[np.float64]:
        raise NotImplementedError

    def make_time_intervals(self,
                            start: float,
                            stop: float,
                            ) -> NDArray[np.float64]:
        out = self.make_matrix(start, stop)
        if out.shape[0] <= 2:
            return out[:, 0]
        return out


class ConstantBeam(ParticleBeam):
    def __init__(self,
                 dose_rate: float,
                 max_dose: float = float("inf"),
                 ):
        super().__init__(dose_rate,
                         Timer.new_constant(),
                         max_dose,
                         )

    def make_matrix(self, start: float, stop: float) -> NDArray[np.float64]:
        if start >= self.time_at_max_dose:
            return np.array([[start, 0.0],
                             [stop, 0.0]])
        if stop < self.time_at_max_dose:
            return np.array([[start, self.peak_dose_rate],
                             [stop, self.peak_dose_rate]])
        return np.array([[start, self.peak_dose_rate],
                         [self.timer.on_time, 0.0],
                         [stop, 0.0]])


class PulsedBeam(ParticleBeam):
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

    def make_matrix(self, start: float, stop: float) -> NDArray[np.float64]:
        times = np.empty([int((stop-start)/self.timer.period)+2, 2],
                         dtype=np.float64)
        times[:, 1] = np.arange(times.shape[0]-1, -1, -1) * self.timer.period
        times[:, 0] = times[:, 1] + self.timer.on_time

        dr = np.zeros_like(times)
        indexes = (np.where(times[:, 1] <= self.time_at_max_dose)[0], [1])
        dr[indexes] = self.peak_dose_rate

        out = np.empty([dr.shape[0]*2, 2])
        out[:, 0] = times.flatten()
        out[:, 1] = dr.flatten()

        while out[-2, 0] < start:
            out = np.delete(out, (-1), axis=0)
        while out[1, 1] <= 0.0:
            out = np.delete(out, (1), axis=0)
        return out

    def make_time_intervals(self,
                            start: float,
                            stop: float,
                            ) -> NDArray[np.float64]:
        out = self.make_matrix(start, stop)
        if out.shape[0] <= 2:
            return out[:, 0]
        delta = self.timer.on_time*0.05
        if delta >= (self.timer.period - self.timer.on_time):
            delta = (self.timer.period - self.timer.on_time)*0.05
        indexes = np.where(out[:, 1] > 0.0)[0], [0]
        if indexes[0][0] == 0:
            indexes = (indexes[0][1:], indexes[1])
        if indexes[0][-1] == out.shape[0]-1:
            indexes = (indexes[0][0:-1], indexes[1])
        out[indexes] = out[indexes] - delta

        indexes = np.where(out[:, 1] == 0.0)[0], [0]
        if indexes[0][0] == 0:
            indexes = (indexes[0][1:], indexes[1])
        if indexes[0][-1] == out.shape[0]-1:
            indexes = (indexes[0][0:-1], indexes[1])
        out[indexes] = out[indexes] + delta
        out[:, 1] /= self.peak_dose_rate
        return out

    @ classmethod
    def from_peak_dose_rate(cls,
                            peak_dr: float,
                            period: float,
                            on_time: float,
                            max_dose: float = float("inf")
                            ) -> PulsedBeam:
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
        return PulsedBeam.from_peak_dose_rate(peak_dr=pulse_dose / on_time,
                                              period=period,
                                              on_time=on_time,
                                              max_dose=max_dose,)


class MonoPulsedBeam(PulsedBeam):
    pass
