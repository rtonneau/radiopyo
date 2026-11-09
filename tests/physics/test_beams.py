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
