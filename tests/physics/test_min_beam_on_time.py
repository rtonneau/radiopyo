#!/usr/bin/env python3
"""
Tests of the shortest configured beam ON time (MIN_BEAM_ON_TIME = 1e-6 s):
a beam whose configured ON time is shorter raises ShortBeamOnError when it
is built (or its timing is updated); a pulse cut short by max_dose is not
a configured ON time and is accepted.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics.beam import (
    BEAM_TYPES,
    MIN_BEAM_ON_TIME,
    ConstantBeam,
    MacroPulsedBeam,
    ParticleBeam,
    PulsedBeam,
    SinglePulseBeam,
)
from radiopyo.physics.exceptions import BeamError, ShortBeamOnError
from radiopyo.simulation.unit_cell import UnitCell

SHORT = 5e-7


def test_threshold_value():
    assert MIN_BEAM_ON_TIME == 1e-6
    assert UnitCell.MIN_BEAM_ON_TIME == MIN_BEAM_ON_TIME
    assert UnitCell.T_ZERO < UnitCell.MIN_BEAM_ON_TIME


def test_exception_is_a_beam_error_and_a_value_error():
    assert issubclass(ShortBeamOnError, BeamError)
    assert issubclass(ShortBeamOnError, ValueError)


# ------------------------------- Too short --------------------------------- #
SHORT_BEAMS: tp.Dict[str, tp.Callable[[], ParticleBeam]] = {
    "pulsed": lambda: PulsedBeam(dose_rate_peak=1e6, period=1e-3, on_time=SHORT),
    "pulsed_dose_per_pulse": lambda: PulsedBeam.from_dose_per_pulse(
        pulse_dose=1.0, period=0.1, on_time=1e-7, n_pulse=3.0
    ),
    "pulsed_dose_rate_peak": lambda: PulsedBeam.from_dose_rate_peak(
        1e6, 1e-3, SHORT, max_dose=1.0
    ),
    "pulsed_from_dict": lambda: PulsedBeam.from_dict(
        {"pulse_dose": 1.0, "dose_rate_peak": 1e7, "period": 0.1, "max_dose": 2.0}
    ),
    "macro_micro": lambda: MacroPulsedBeam(
        dose_rate_peak=1e6,
        period=1e-3,
        on_time=3e-5,
        micro_period=1e-5,
        micro_on_time=SHORT,
    ),
    "macro_train": lambda: MacroPulsedBeam(
        dose_rate_peak=1e6,
        period=1e-3,
        on_time=SHORT,
        micro_period=1e-5,
        micro_on_time=SHORT / 2,
    ),
    "constant_max_dose": lambda: ConstantBeam(dose_rate=1e7, max_dose=1.0),
    "single_pulse": lambda: SinglePulseBeam(dose_rate=1e7, max_dose=1.0),
    "single_pulse_on_time": lambda: SinglePulseBeam.from_on_time(
        dose_rate=1e6, on_time=SHORT
    ),
    "single_pulse_start_time": lambda: SinglePulseBeam(
        dose_rate=1e7, max_dose=1.0, start_time=1e-3
    ),
}


@pytest.mark.parametrize("name", list(SHORT_BEAMS))
def test_short_configured_on_time_raises(name: str):
    with pytest.raises(ShortBeamOnError, match="MIN_BEAM_ON_TIME"):
        SHORT_BEAMS[name]()


def test_macro_train_shorter_than_threshold_raises_on_the_train():
    with pytest.raises(ShortBeamOnError, match=r"on_time=5e-07"):
        MacroPulsedBeam(
            dose_rate_peak=1e6,
            period=1e-3,
            on_time=SHORT,
            micro_period=2e-6,
            micro_on_time=1e-6,
            micro_start_time=0.0,
        )


def test_beam_types_registry_from_dict_raises():
    with pytest.raises(ShortBeamOnError):
        BEAM_TYPES["pulsed"].from_dict(
            {"dose_rate_peak": 1e6, "period": 1e-3, "on_time": SHORT}
        )


def test_value_error_still_catches_it():
    with pytest.raises(ValueError):
        PulsedBeam(dose_rate_peak=1e6, period=1e-3, on_time=SHORT)


@pytest.mark.parametrize(
    "update",
    [{"on_time": SHORT}, {"on_time": 1e-7, "period": 1e-4}],
    ids=["on_time", "on_time_period"],
)
def test_update_beam_param_to_a_short_on_time_raises(update: tp.Dict[str, float]):
    beam = PulsedBeam(dose_rate_peak=1e6, period=1e-3, on_time=1e-5)
    with pytest.raises(ShortBeamOnError):
        beam.update_beam_param(**update)


def test_update_macro_micro_on_time_raises():
    beam = MacroPulsedBeam(
        dose_rate_peak=1e6, period=1e-3, on_time=3e-5, micro_period=1e-5, micro_on_time=2e-6
    )
    with pytest.raises(ShortBeamOnError):
        beam.update_beam_param(micro_on_time=SHORT)


def test_update_single_pulse_max_dose_raises():
    beam = SinglePulseBeam(dose_rate=1e6, max_dose=2.0)
    with pytest.raises(ShortBeamOnError):
        beam.update_beam_param(max_dose=0.1)


# ------------------------------- Accepted ---------------------------------- #
ACCEPTED_BEAMS: tp.Dict[str, tp.Callable[[], ParticleBeam]] = {
    "pulsed_at_threshold": lambda: PulsedBeam(
        dose_rate_peak=1e6, period=1e-3, on_time=MIN_BEAM_ON_TIME
    ),
    # 1 Gy at 1e6 Gy/s: on_time computed as 1 / 1e6
    "pulsed_computed_at_threshold": lambda: PulsedBeam.from_dict(
        {"pulse_dose": 1.0, "dose_rate_peak": 1e6, "period": 0.1, "max_dose": 2.0}
    ),
    "macro_at_threshold": lambda: MacroPulsedBeam(
        dose_rate_peak=1e6,
        period=1e-3,
        on_time=3e-5,
        micro_period=1e-5,
        micro_on_time=MIN_BEAM_ON_TIME,
    ),
    "single_pulse_at_threshold": lambda: SinglePulseBeam(dose_rate=1e6, max_dose=1.0),
    "single_pulse_on_time_at_threshold": lambda: SinglePulseBeam.from_on_time(
        dose_rate=3e6, on_time=MIN_BEAM_ON_TIME
    ),
    "constant_without_max_dose": lambda: ConstantBeam(dose_rate=1e9),
    # the last pulse is cut by max_dose to 0.1 us: not a configured ON time
    "pulsed_partial_last_pulse": lambda: PulsedBeam.from_dose_rate_peak(
        1e6, 1e-3, 1e-5, max_dose=10.1
    ),
    # less than one pulse: the single (partial) pulse is cut by max_dose
    "pulsed_partial_first_pulse": lambda: PulsedBeam.from_dose_rate_peak(
        1e6, 1e-3, 1e-5, max_dose=0.1
    ),
}


@pytest.mark.parametrize("name", list(ACCEPTED_BEAMS))
def test_configured_on_time_at_or_above_threshold_is_accepted(name: str):
    beam = ACCEPTED_BEAMS[name]()
    assert beam.time_at_max_dose > 0
