# --------------------------------- External imports --------------------------------- #
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from math import isinf, isnan
from typing import (
    NamedTuple,
)

# --------------------------------- Internal imports --------------------------------- #

# ----------------------------- Type checking ONLY imports---------------------------- #

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


class ConstantBeam(ParticleBeam):
    def __init__(self,
                 dose_rate: float,
                 max_dose: float = float("inf"),
                 ):
        super().__init__(dose_rate,
                         Timer.new_constant(),
                         max_dose,
                         )


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

    @classmethod
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

    @classmethod
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
