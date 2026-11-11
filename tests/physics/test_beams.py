#!/usr/bin/env python3
"""
Tests of the beam module split (timer / beam / beam_collection): public
imports, beam type registry, constructor argument order, dictionary
builders (no mutation of caller data) and crash fixes.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import math

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics import beam as beam_module
from radiopyo.physics import beam_collection as beam_collection_module
from radiopyo.physics import timer as timer_module
from radiopyo.physics.beam import (
    BEAM_TYPES,
    BeamCollection,
    ConstantBeam,
    MacroPulsedBeam,
    PulsedBeam,
    SinglePulseBeam,
)
from radiopyo.physics.exceptions import BeamError, ParameterError
from radiopyo.physics.timer import Timer, TimerState


# --------------------------------------------------------------------------- #
#                                Module layout                                #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", ["BaseTimer", "Timer", "MacroTimer", "TimerState"])
def test_timers_reexported_from_beam(name):
    assert getattr(beam_module, name) is getattr(timer_module, name)


def test_beam_collection_reexported_from_beam():
    assert beam_module.BeamCollection is beam_collection_module.BeamCollection


def test_package_level_beams():
    import radiopyo

    assert radiopyo.PulsedBeam is PulsedBeam
    assert radiopyo.ConstantBeam is ConstantBeam


def test_parameter_error_is_beam_error():
    assert issubclass(ParameterError, BeamError)


# --------------------------------------------------------------------------- #
#                                  Registry                                   #
# --------------------------------------------------------------------------- #
def test_registry_contents():
    assert BEAM_TYPES == {
        "constant": ConstantBeam,
        "pulsed": PulsedBeam,
        "singlepulse": SinglePulseBeam,
        "macro_pulsed": MacroPulsedBeam,
    }


def test_registry_keys_match_beam_type():
    samples = [
        ConstantBeam(1.0),
        SinglePulseBeam(1.0, 1.0),
        PulsedBeam(1.0, 2.0, 1.0),
        MacroPulsedBeam(1.0, 10.0, 3.0, 1.0, 0.25),
    ]
    for beam in samples:
        assert BEAM_TYPES[beam.beam_type().lower()] is type(beam)


def test_unknown_beam_type():
    with pytest.raises(NotImplementedError):
        BeamCollection().add_beam_from_dict("nope", {"dose_rate": 1.0})


def test_beam_type_is_case_insensitive():
    beams = BeamCollection()
    beams.add_beam_from_dict("Constant", {"dose_rate": 1.0, "label": "c"})
    assert isinstance(beams["c"][0], ConstantBeam)


# --------------------------------------------------------------------------- #
#                         Constructor argument order                          #
# --------------------------------------------------------------------------- #
def test_from_dose_rate_peak_positional_order():
    beam = PulsedBeam.from_dose_rate_peak(1000, 1e-3, 5e-6, 10)
    # 4th positional is n_pulse: max dose = on_time * n_pulse * peak
    assert beam.max_dose == pytest.approx(5e-6 * 10 * 1000)


def test_from_dose_per_pulse_positional_order():
    beam = PulsedBeam.from_dose_per_pulse(0.005, 1e-3, 5e-6, 10)
    assert beam.max_dose == pytest.approx(0.05)


def test_pulsed_constructor_positional_order():
    beam = PulsedBeam(1000, 1e-3, 5e-6, 0.5, 2.0)
    assert (beam.dose_rate_peak, beam.timer.period, beam.timer.on_time) == (
        1000,
        1e-3,
        5e-6,
    )
    assert (beam.max_dose, beam.timer.start_time) == (0.5, 2.0)


def test_macro_constructor_positional_order():
    beam = MacroPulsedBeam(1000, 10.0, 3.0, 1.0, 0.25, 5.0, 0.5, 2.0)
    assert (beam.macro_timer.period, beam.macro_timer.on_time) == (10.0, 3.0)
    assert (beam.timer.period, beam.timer.on_time) == (1.0, 0.25)
    assert beam.max_dose == 5.0
    assert beam.timer.start_time == 0.5  # micro_start_time
    assert beam.macro_timer.start_time == 2.0  # start_time


# --------------------------------------------------------------------------- #
#                         from_dict never mutates input                       #
# --------------------------------------------------------------------------- #
PULSED_CONFIGS = [
    {"frequency": 1000.0, "dose_rate_peak": 1000.0, "on_time": 5e-6},
    {"period": 1e-3, "pulse_dose": 5e-3, "on_time": 5e-6},
    {"period": 1e-3, "pulse_dose": 5e-3, "dose_rate_peak": 1000.0},
    {"period": 1e-3, "dose_rate": 5.0, "on_time": 5e-6, "max_dose": 1.0},
]


@pytest.mark.parametrize("config", PULSED_CONFIGS)
def test_pulsed_from_dict_does_not_mutate(config):
    before = copy.deepcopy(config)
    PulsedBeam.from_dict(config)
    assert config == before


def test_macro_from_dict_does_not_mutate():
    config = {
        "frequency": 0.1,
        "on_time": 3.0,
        "micro_frequency": 1.0,
        "micro_on_time": 0.25,
        "dose_rate_peak": 100.0,
    }
    before = copy.deepcopy(config)
    MacroPulsedBeam.from_dict(config)
    assert config == before


@pytest.mark.parametrize(
    "beam_type, config",
    [
        ("constant", {"dose_rate": 1.0, "label": "a"}),
        ("singlepulse", {"dose_rate": 1.0, "max_dose": 2.0, "label": "a"}),
        (
            "pulsed",
            {"frequency": 1000.0, "dose_rate_peak": 1000.0, "on_time": 5e-6},
        ),
        (
            "macro_pulsed",
            {
                "frequency": 0.1,
                "on_time": 3.0,
                "micro_frequency": 1.0,
                "micro_on_time": 0.25,
                "dose_rate_peak": 100.0,
            },
        ),
    ],
)
def test_collection_builders_do_not_mutate(beam_type, config):
    before = copy.deepcopy(config)
    BeamCollection.build_beam_from_dict(beam_type, config, LET=2.0, verbose=False)
    assert config == before

    beams = BeamCollection()
    beams.add_beam_from_dict(beam_type, config, LET=2.0, verbose=False)
    assert config == before


def test_from_config_dict_does_not_mutate_and_keeps_let():
    config = {
        "constant": {"dose_rate": 1.0, "label": "c", "LET": 3.0},
        "pulsed": [
            {
                "period": 1e-3,
                "dose_rate_peak": 1000.0,
                "on_time": 5e-6,
                "label": "p",
                "LET": 5.0,
            }
        ],
    }
    before = copy.deepcopy(config)
    beams = BeamCollection.from_config_dict(config, verbose=False)
    assert config == before
    assert beams["c"][1] == 3.0
    assert beams["p"][1] == 5.0


def test_as_dict_round_trip():
    beams = BeamCollection()
    beams.add_beam(ConstantBeam(2.0), key="c", let=1.0)
    beams.add_beam(SinglePulseBeam(2.0, 1.0), key="s")
    beams.add_beam(PulsedBeam(1000.0, 1e-3, 5e-6, 0.5), key="p", let=4.0)
    beams.add_beam(MacroPulsedBeam(100.0, 10.0, 3.0, 1.0, 0.25, 5.0), key="m")
    out = beams.as_dict()
    assert set(out) == {"constant", "singlepulse", "pulsed", "macro_pulsed"}

    rebuilt = BeamCollection.from_config_dict(out, verbose=False)
    assert rebuilt.as_dict() == out


# --------------------------------------------------------------------------- #
#                                 Crash fixes                                 #
# --------------------------------------------------------------------------- #
def test_update_beam_param_period():
    beam = PulsedBeam(1000.0, 1e-3, 5e-6)
    beam.update_beam_param(period=2e-3)
    assert beam.timer.period == 2e-3
    assert beam.timer.on_time == 5e-6


def test_update_beam_param_frequency():
    beam = PulsedBeam(1000.0, 1e-3, 5e-6)
    beam.update_beam_param(frequency=500.0)
    assert beam.timer.period == pytest.approx(2e-3)


@pytest.mark.parametrize(
    "beam_class, config",
    [
        (PulsedBeam, {"period": 1e-3}),  # not enough parameters, no label
        (PulsedBeam, {"on_time": 5e-6}),  # no period
        (MacroPulsedBeam, {"period": 10.0}),
        (SinglePulseBeam, {"max_dose": 1.0}),
    ],
)
def test_unlabelled_bad_config_raises_parameter_error(beam_class, config):
    with pytest.raises(ParameterError):
        beam_class.from_dict(config)


def test_unlabelled_bad_config_through_collection():
    with pytest.raises(ParameterError):
        BeamCollection().add_beam_from_dict("pulsed", {"period": 1e-3}, verbose=False)


# --------------------------------------------------------------------------- #
#                              Removed duplicates                             #
# --------------------------------------------------------------------------- #
def test_timer_pulse_end_is_on():
    timer = Timer.new_pulsed(period=2, on_time=0.5)
    assert timer.state_at(0.5) is TimerState.isOn
    assert "ON" in (Timer.state_at.__doc__ or "")


@pytest.mark.parametrize(
    "owner, name",
    [
        ("PulsedBeam", "delivery"),
        ("PulsedBeam", "_compute_dose_rate_peak"),
        ("MacroTimer", "within_train"),
        ("MacroTimer", "micro_pulses_left_in_train"),
        ("MacroTimer", "micro_pulse_index"),
        ("MacroTimer", "time_left_in_train"),
    ],
)
def test_unused_methods_removed(owner, name):
    cls = getattr(beam_module, owner)
    assert not hasattr(cls, name)


def test_duplicated_timer_methods_defined_once():
    assert "time_till_next_period" not in vars(Timer)
    assert "period_index" not in vars(timer_module.MacroTimer)


# --------------------------------------------------------------------------- #
#                         Generic make_matrix (all beams)                      #
# --------------------------------------------------------------------------- #
def test_make_matrix_is_defined_once():
    # One generic implementation walking the timer transitions
    for cls in (ConstantBeam, SinglePulseBeam, PulsedBeam, MacroPulsedBeam):
        assert "make_matrix" not in vars(cls)
        assert "make_time_intervals" not in vars(cls)


def test_constant_beam_on_at_stop_closes_with_zero():
    beam = ConstantBeam(dose_rate=1e3, max_dose=5.0, start_time=0.2)
    matrix = beam.make_matrix(stop=0.204, start=0.201)
    np.testing.assert_array_equal(matrix, [[0.201, 1e3], [0.204, 0.0]])
    intervals = beam.make_time_intervals(stop=0.204, start=0.201)
    np.testing.assert_array_equal(intervals, [[0.201, 1.0], [0.204, 0.0]])


def test_single_pulse_matrix():
    beam = SinglePulseBeam(dose_rate=1e6, max_dose=2.0, start_time=1e-3)
    matrix = beam.make_matrix(stop=0.1)
    expected = [[0.0, 0.0], [1e-3, 1e6], [beam.time_at_max_dose, 0.0], [0.1, 0.0]]
    np.testing.assert_allclose(matrix, expected, rtol=1e-12)


def test_make_time_intervals_state_column_is_zero_or_one():
    # Start before the beam: the first row is OFF (it used to be flagged ON)
    beam = ConstantBeam(dose_rate=1e3, max_dose=5.0, start_time=0.2)
    intervals = beam.make_time_intervals(stop=0.204, start=0.0)
    np.testing.assert_array_equal(intervals[:, 1], [0.0, 1.0, 0.0])
    # The ON to OFF transition of a pulse is not expanded (ADR 0007)
    beam = PulsedBeam(1e6, 1e-2, 2e-6, max_dose=7.0)
    intervals = beam.make_time_intervals(stop=0.05)
    assert set(intervals[:, 1].tolist()) == {0.0, 1.0}
    assert intervals[1, 0] == 2e-6


# --------------------------------------------------------------------------- #
#                               Pulsed / averages                             #
# --------------------------------------------------------------------------- #
def test_pulsed_average_dose_rate_accounts_for_start_time():
    beam = PulsedBeam(1e6, 1e-2, 2e-6, max_dose=7.0, start_time=1e-3)
    assert beam.tot_pulses() == pytest.approx(3.5)
    # 7 Gy delivered between start_time and time_at_max_dose
    assert beam.dose_rate == pytest.approx(7.0 / (beam.time_at_max_dose - 1e-3))
    same_without_delay = PulsedBeam(1e6, 1e-2, 2e-6, max_dose=7.0)
    assert same_without_delay.tot_pulses() == pytest.approx(3.5)
    assert same_without_delay.dose_rate == pytest.approx(beam.dose_rate)


@pytest.mark.parametrize("start_time", [None, 0.5])
def test_pulsed_without_max_dose_has_no_nan(start_time):
    beam = PulsedBeam(1e3, 1e-2, 1e-3, start_time=start_time)
    assert beam.time_at_max_dose == math.inf
    assert beam.tot_pulses() == math.inf
    assert beam.dose_rate == pytest.approx(100.0)  # peak * duty cycle


def test_from_dose_rate_peak_float_tolerance():
    # 0.1 * 3 * 1.0 is 0.30000000000000004, not 0.3
    beam = PulsedBeam.from_dose_rate_peak(
        peak_dr=1.0, period=1.0, on_time=0.1, n_pulse=3.0, max_dose=0.3
    )
    assert beam.max_dose == pytest.approx(0.3)
    with pytest.raises(ValueError):
        PulsedBeam.from_dose_rate_peak(
            peak_dr=1.0, period=1.0, on_time=0.1, n_pulse=3.0, max_dose=0.4
        )


# --------------------------------------------------------------------------- #
#                                MacroPulsedBeam                              #
# --------------------------------------------------------------------------- #
def _macro(**kwargs):
    param = {
        "dose_rate_peak": 1e6,
        "period": 0.1,
        "on_time": 1e-3,
        "micro_period": 1e-4,
        "micro_on_time": 2e-6,
        "start_time": 0.01,
    }
    param.update(kwargs)
    return MacroPulsedBeam(**param)


def test_macro_without_max_dose_terminates_and_respects_stop():
    beam = _macro()
    assert beam.time_at_max_dose == math.inf
    matrix = beam.make_matrix(stop=0.35)
    assert matrix[-1].tolist() == [0.35, 0.0]
    assert matrix[:, 0].max() == 0.35
    assert np.all(np.diff(matrix[:, 0]) > 0)
    # 4 trains (0.01, 0.11, 0.21, 0.31) of 10 micro pulses
    assert np.count_nonzero(matrix[:, 1]) == 40
    intervals = beam.make_time_intervals(stop=0.35)
    assert intervals[-1].tolist() == [0.35, 0.0]


def test_macro_without_max_dose_no_nan():
    beam = _macro()
    assert beam.tot_pulses() == math.inf
    # peak * (10 micro pulses of 2 us every 0.1 s)
    assert beam.dose_rate == pytest.approx(1e6 * 10 * 2e-6 / 0.1)


def test_macro_stop_inside_the_beam():
    beam = _macro(max_dose=50.0)
    matrix = beam.make_matrix(stop=0.05)
    assert matrix[0].tolist() == [0.0, 0.0]
    assert matrix[-1].tolist() == [0.05, 0.0]
    assert matrix.shape[0] == 22  # start + 10 pulses (ON, OFF) + stop


def test_macro_no_duplicate_row_at_train_start():
    beam = _macro(max_dose=50.0)
    for start, stop in [(0.0, 1.0), (0.01, 0.5), (0.05, 1.0)]:
        matrix = beam.make_matrix(stop=stop, start=start)
        assert np.all(np.diff(matrix[:, 0]) > 0), (start, stop)
    matrix = beam.make_matrix(stop=1.0)
    # 2.5 trains of 10 micro pulses: 25 pulses, start row and stop row
    assert matrix.shape[0] == 2 * 25 + 2
    assert matrix[1].tolist() == [0.01, 1e6]


def test_macro_matrix_matches_instant_dose_rate():
    beam = _macro(max_dose=50.0, micro_start_time=2e-5)
    matrix = beam.make_matrix(stop=0.3)
    for (t0, dr), (t1, _) in zip(matrix[:-1], matrix[1:]):
        middle = 0.5 * (t0 + t1)
        assert beam.at(middle).dose_rate() == dr, middle


def test_macro_time_at_max_dose():
    # Dose per train: 10 pulses of 2 Gy. 2.5 trains: 5 pulses of the last
    beam = _macro(max_dose=50.0)
    assert beam.time_at_max_dose == pytest.approx(0.01 + 2 * 0.1 + 4e-4 + 2e-6)
    # whole trains: reached at the end of the last micro pulse
    beam = _macro(max_dose=40.0)
    assert beam.time_at_max_dose == pytest.approx(0.01 + 0.1 + 9e-4 + 2e-6)
    assert beam.at(0.1105).dose_rate() == 1e6
    assert beam.at(0.1110).dose_rate() == 0.0


def test_macro_time_at_max_dose_float_tolerance():
    # 3 * dose_per_train is not exactly representable the same way twice
    peak = 1e4 / 3
    train_dose = _macro(dose_rate_peak=peak).dose_per_train()
    beam = _macro(dose_rate_peak=peak, max_dose=3 * train_dose)
    assert beam.time_at_max_dose == pytest.approx(0.01 + 2 * 0.1 + 9e-4 + 2e-6)


def test_macro_average_dose_rate_accounts_for_start_time():
    beam = _macro(max_dose=50.0)
    assert beam.dose_rate == pytest.approx(50.0 / (beam.time_at_max_dose - 0.01))
    assert _macro(max_dose=50.0, start_time=0.0).dose_rate == pytest.approx(
        beam.dose_rate
    )


def test_macro_micro_start_time_reaches_the_macro_timer():
    beam = _macro(micro_start_time=2e-5, max_dose=50.0)
    assert beam.macro_timer.subtimer.start_time == 2e-5
    assert beam.macro_timer.subtimer is beam.timer
    # first micro pulse of a train: 2e-5 after the train start
    assert beam.at(0.01001).dose_rate() == 0.0
    assert beam.at(0.010021).dose_rate() == 1e6
    assert beam.make_matrix(stop=1.0)[1].tolist() == [0.01 + 2e-5, 1e6]


def test_macro_update_beam_param_keeps_timers_in_sync():
    beam = _macro(max_dose=50.0)
    beam.update_beam_param(micro_period=2e-4, micro_on_time=4e-6, start_time=0.02)
    assert beam.timer is beam.macro_timer.subtimer
    assert (beam.timer.period, beam.timer.on_time) == (2e-4, 4e-6)
    assert beam.macro_timer.start_time == 0.02
    assert beam.macro_timer.period == 0.1
    assert beam.as_dict()["micro_period"] == 2e-4
    assert beam.at(0.02).dose_rate() == 1e6

    beam.update_beam_param(frequency=5.0, micro_start_time=1e-5, max_dose=40.0)
    assert beam.macro_timer.period == pytest.approx(0.2)
    assert beam.timer is beam.macro_timer.subtimer
    assert beam.timer.start_time == 1e-5
    # 5 micro pulses of 4 Gy per train: 2 trains
    assert beam.time_at_max_dose == pytest.approx(0.02 + 0.2 + 8e-4 + 4e-6 + 1e-5)

    beam.update_beam_param(dose_rate_peak=2e6)
    assert beam.dose_rate_peak == 2e6
    assert beam.at(0.02 + 1e-5 + 1e-6).dose_rate() == 2e6


# --------------------------------------------------------------------------- #
#                               BeamCollection.merge                          #
# --------------------------------------------------------------------------- #
def test_merge_overrides_same_name_same_type():
    first = BeamCollection()
    first.add_beam(PulsedBeam(1e3, 1e-2, 1e-3, max_dose=3.0), "a", let=1.0)
    first.add_beam(ConstantBeam(5.0), "c")
    second = BeamCollection()
    second.add_beam(PulsedBeam(1e3, 1e-2, 1e-3, max_dose=5.0), "a")
    second.add_beam(PulsedBeam(2e3, 1e-2, 1e-3, max_dose=7.0), "b")
    second.add_beam(SinglePulseBeam(10.0, 1.0), "c")

    merged = first.merge(second, verbose=False)

    assert sorted(name for name, _, _ in merged) == ["a", "b", "c"]
    beam_a, let_a = merged["a"]
    assert beam_a.max_dose == 5.0
    assert let_a == 1.0  # kept: the second collection has no LET for "a"
    assert merged["b"][0].dose_rate_peak == 2e3
    assert isinstance(merged["c"][0], SinglePulseBeam)  # other type: replaced
    # the operands are left untouched
    assert first["a"][0].max_dose == 3.0
    assert isinstance(first["c"][0], ConstantBeam)
    assert len(second) == 3


def test_merge_macro_beam():
    first = BeamCollection()
    first.add_beam(_macro(max_dose=50.0), "m")
    second = BeamCollection()
    second.add_beam(_macro(max_dose=40.0), "m")
    merged = first.merge(second, verbose=False)
    assert merged["m"][0].max_dose == 40.0


def test_macro_micro_start_beyond_train_on_time_raises() -> None:
    """No micro pulse fits in a train: reject it instead of looping forever."""
    with pytest.raises(ValueError, match="micro_start_time"):
        MacroPulsedBeam(10.0, 1.0, 0.205, 0.1, 0.01, micro_start_time=0.3)
    beam = MacroPulsedBeam(10.0, 1.0, 0.205, 0.1, 0.01)
    with pytest.raises(ValueError, match="micro_start_time"):
        beam.update_beam_param(micro_start_time=0.205)
