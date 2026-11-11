#!/usr/bin/env python3
"""
Tests of ParticleBeam.make_time_intervals for every beam type: the rows
are the exact beam transitions (ADR 0007), an ON interval ends exactly at
the end of the pulse (or at the time of max dose) and the beam dose rate
is constant inside every interval.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics.beam import (
    ConstantBeam,
    ConstantDoseRate,
    MacroPulsedBeam,
    ParticleBeam,
    PulsedBeam,
    SinglePulseBeam,
)


def _macro() -> MacroPulsedBeam:
    return MacroPulsedBeam(
        dose_rate_peak=1e6,
        period=0.1,
        on_time=1e-3,
        micro_period=1e-4,
        micro_on_time=2e-6,
        start_time=0.01,
    )


# (id, beam, (start, stop), width of a pulse (None: single ON interval))
CASES: tp.List[tp.Tuple[str, ParticleBeam, tp.Tuple[float, float], tp.Optional[float]]] = [
    (
        "constant",
        ConstantBeam(dose_rate=1e3, max_dose=5.0, start_time=0.2),
        (0.0, 1.0),
        None,
    ),
    (
        "constant_stop_inside",
        ConstantBeam(dose_rate=1e3, max_dose=5.0, start_time=0.2),
        (0.201, 0.204),
        None,
    ),
    (
        "single_pulse",
        SinglePulseBeam(dose_rate=1e6, max_dose=2.0, start_time=1e-3),
        (0.0, 0.1),
        None,
    ),
    (
        "pulsed",
        PulsedBeam(1e6, 1e-2, 2e-6, max_dose=7.0),
        (0.0, 0.05),
        2e-6,
    ),
    (
        "pulsed_start_time",
        PulsedBeam(1e6, 1e-2, 2e-6, max_dose=7.0, start_time=1e-3),
        (0.0, 0.05),
        2e-6,
    ),
    ("macro_pulsed", _macro(), (0.0, 0.35), 2e-6),
]
IDS = [case[0] for case in CASES]


def _widths(intervals: np.ndarray) -> np.ndarray:
    return np.diff(intervals[:, 0])


@pytest.mark.parametrize(("_name", "beam", "span", "_width"), CASES, ids=IDS)
def test_intervals_cover_the_span(_name, beam, span, _width):
    start, stop = span
    intervals = beam.make_time_intervals(stop=stop, start=start)
    assert intervals[0, 0] == start
    assert intervals[-1, 0] == stop
    assert intervals[-1, 1] == 0.0
    assert np.all(_widths(intervals) > 0)
    assert set(intervals[:, 1].tolist()) <= {0.0, 1.0}


@pytest.mark.parametrize(("_name", "beam", "span", "width"), CASES, ids=IDS)
def test_on_intervals_end_exactly_at_the_pulse_end(_name, beam, span, width):
    start, stop = span
    intervals = beam.make_time_intervals(stop=stop, start=start)
    on = np.flatnonzero(intervals[:-1, 1] > 0)
    assert on.size > 0
    for i in on:
        t_on, t_off = intervals[i, 0], intervals[i + 1, 0]
        if width is None:  # one ON interval, ends at the max dose
            expected = min(beam.time_at_max_dose, stop)
        else:
            expected = min(t_on + width, beam.time_at_max_dose, stop)
        assert t_off == expected
        # and the next interval is OFF (never a refined tail after the pulse)
        assert intervals[i + 1, 1] == 0.0


@pytest.mark.parametrize(("_name", "beam", "span", "_width"), CASES, ids=IDS)
def test_dose_rate_is_constant_in_every_interval(_name, beam, span, _width):
    start, stop = span
    intervals = beam.make_time_intervals(stop=stop, start=start)
    for (t0, state), (t1, _) in zip(intervals[:-1], intervals[1:], strict=True):
        expected = beam.dose_rate_peak if state > 0 else 0.0
        samples = t0 + (t1 - t0) * np.array([1e-6, 0.25, 0.5, 0.75, 1 - 1e-6])
        for time in samples:
            assert beam.at(time).dose_rate() == expected, (t0, t1, time)


@pytest.mark.parametrize(("_name", "beam", "span", "_width"), CASES, ids=IDS)
def test_intervals_integrate_to_the_delivered_dose(_name, beam, span, _width):
    start, stop = span
    intervals = beam.make_time_intervals(stop=stop, start=start)
    dose = np.sum(intervals[:-1, 1] * _widths(intervals)) * beam.dose_rate_peak
    reached = min(stop, beam.time_at_max_dose) >= beam.time_at_max_dose
    if reached and start <= beam.timer.start_time:
        assert dose == pytest.approx(beam.max_dose, rel=1e-9)
    else:
        assert dose <= beam.max_dose * (1 + 1e-9)


def test_pulsed_intervals_are_exact_transitions():
    beam = PulsedBeam(1e6, 1e-2, 2e-6, max_dose=7.0)
    intervals = beam.make_time_intervals(stop=0.05)
    # 7 Gy = 3 pulses of 2 Gy + half of a 4th one, cut at the max dose
    expected_time = [0.0, 2e-6, 1e-2, 1e-2 + 2e-6, 2e-2, 2e-2 + 2e-6, 3e-2]
    expected_time += [beam.time_at_max_dose, 0.05]
    assert beam.time_at_max_dose == pytest.approx(3e-2 + 1e-6, rel=1e-12)
    np.testing.assert_allclose(intervals[:, 0], expected_time, rtol=1e-12)
    np.testing.assert_array_equal(intervals[:, 1], [1, 0, 1, 0, 1, 0, 1, 0, 0])


def test_growth_time_is_gone():
    assert not hasattr(ParticleBeam, "GROWTH_TIME")
    assert not hasattr(PulsedBeam, "GROWTH_TIME")


def test_constant_dose_rate_has_the_beam_interface():
    source = ConstantDoseRate(12.5)
    for time in (0.0, 1e-9, 3.0):
        message = source.at(time)
        assert message.time == time
        assert message.dose_rate() == 12.5
