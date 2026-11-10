#!/usr/bin/env python3
"""
Tests pinning the current behaviour of the beam timers
(radiopyo.physics.beam.Timer and MacroTimer).

Expected values were computed by hand from the timer definitions
(period / on_time / start_time). Behaviours that later tickets are
expected to change are flagged with a comment naming the ticket and are
NOT to be read as "correct".

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import math

import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics.beam import MacroTimer, Timer, TimerState

ON, OFF = TimerState.isOn, TimerState.isOff


# --------------------------------------------------------------------------- #
#                                    Timer                                    #
# --------------------------------------------------------------------------- #
class TestTimer:
    @pytest.fixture
    def timer(self):
        # ON during [0, 0.5] every 2 s
        return Timer.new_pulsed(period=2, on_time=0.5)

    @pytest.mark.parametrize(
        "time, state",
        [
            (-0.1, OFF),  # before start
            (0.0, ON),
            (0.25, ON),
            (0.5, ON),  # boundary is inclusive (on_time included)
            (0.75, OFF),
            (1.5, OFF),
            (2.0, ON),  # 2nd period
            (2.25, ON),
            (2.75, OFF),
            (6.25, ON),  # 4th period
        ],
    )
    def test_state_at(self, timer, time, state):
        assert timer.state_at(time) is state

    @pytest.mark.parametrize(
        "time, till_on, till_off",
        [
            (0.25, 0.0, 0.25),
            (1.0, 1.0, 0.0),
            (0.75, 1.25, 0.0),
            (2.25, 0.0, 0.25),
            (3.5, 0.5, 0.0),
        ],
    )
    def test_time_till_next(self, timer, time, till_on, till_off):
        assert timer.time_till_next_ON(time) == pytest.approx(till_on)
        assert timer.time_till_next_OFF(time) == pytest.approx(till_off)

    def test_period_helpers(self, timer):
        assert timer.period_index(0.0) == 0
        assert timer.period_index(1.99) == 0
        assert timer.period_index(2.0) == 1
        assert timer.period_index(8.5) == 4
        assert timer.time_till_next_period(0.5) == pytest.approx(1.5)
        assert timer.time_till_next_period(2.5) == pytest.approx(1.5)
        assert timer.off_time() == pytest.approx(1.5)

    def test_duty_cycle(self, timer):
        assert timer.duty_cycle() == pytest.approx(0.25)
        assert not timer.is_continuous()
        assert not timer.is_mono_pulse()

    def test_start_time(self):
        timer = Timer.new_pulsed(period=2, on_time=0.5, start_time=1.0)
        assert timer.state_at(0.5) is OFF  # before start
        assert timer.state_at(1.0) is ON
        assert timer.state_at(1.25) is ON
        assert timer.state_at(1.75) is OFF
        assert timer.state_at(3.25) is ON
        assert timer.time_till_next_ON(0.5) == pytest.approx(0.5)
        assert timer.time_till_next_ON(1.75) == pytest.approx(1.25)
        assert timer.time_till_next_OFF(3.25) == pytest.approx(0.25)

    def test_constant(self):
        timer = Timer.new_constant()
        assert timer.is_continuous()
        assert timer.duty_cycle() == 1.0
        for time in (0.0, 0.25, 1.5, 1e6):
            assert timer.state_at(time) is ON
            assert timer.time_till_next_ON(time) == 0.0
            # Quirk (no ticket assigned): period_index is 0 and
            # `inf * 0` makes time_till_next_OFF NaN rather than inf.
            assert math.isnan(timer.time_till_next_OFF(time))

    def test_mono_pulse(self):
        timer = Timer.new_pulsed(period=float("inf"), on_time=0.5)
        assert timer.is_mono_pulse()
        assert not timer.is_continuous()
        assert timer.state_at(0.25) is ON
        assert timer.state_at(1.0) is OFF
        assert timer.state_at(1e6) is OFF
        # Same `inf * 0 = nan` quirk as for the constant timer
        assert math.isnan(timer.time_till_next_OFF(0.25))

    def test_copy_keeps_start_time(self):
        timer = Timer.new_pulsed(period=2, on_time=0.5, start_time=1.0)
        copy = timer.copy()
        assert copy is not timer
        assert (copy.period, copy.on_time) == (2, 0.5)
        assert copy.start_time == 1.0
        assert timer.copy(start_time=3.0).start_time == 3.0


# --------------------------------------------------------------------------- #
#                                  MacroTimer                                 #
# --------------------------------------------------------------------------- #
class TestMacroTimer:
    @pytest.fixture
    def macro(self):
        # Trains of 3 s every 10 s, each made of 0.25 s pulses every 1 s
        return MacroTimer.new_macropulse(
            macro_period=10,
            macro_on_time=3,
            micro_period=1,
            micro_on_time=0.25,
        )

    def test_invalid_periods(self):
        with pytest.raises(ValueError):
            MacroTimer.new_macropulse(1, 0.5, 2, 0.25)

    @pytest.mark.parametrize(
        "time, state",
        [
            (-1.0, OFF),
            (0.0, ON),
            (0.1, ON),
            (0.5, OFF),
            (1.1, ON),
            (2.1, ON),
            (2.5, OFF),
            (3.0, ON),  # boundary of the train, inclusive
            (3.5, OFF),  # train is over
            (9.9, OFF),
            (10.1, ON),  # second train
            (10.5, OFF),
            (11.2, ON),
            (13.5, OFF),
            (20.05, ON),  # third train
        ],
    )
    def test_state_at(self, macro, time, state):
        assert macro.state_at(time) is state

    def test_period_indexes(self, macro):
        assert macro.period_index(0.0) == 0
        assert macro.period_index(9.99) == 0
        assert macro.period_index(10.0) == 1
        assert macro.period_index(25.0) == 2
        assert macro.micro_period_index(0.5) == 0
        assert macro.micro_period_index(2.5) == 2
        assert macro.micro_period_index(12.5) == 2
        assert macro.micro_period_index(5.0) == -1  # out of the train

    @pytest.mark.parametrize(
        "time, till_on, till_off",
        [
            (0.1, 0.0, 0.15),
            (0.5, 0.5, 0.0),
            (1.9, 0.1, 0.0),
            (11.5, 0.5, 0.0),
            (12.1, 0.0, 0.15),
            (5.0, 5.0, 0.0),  # between trains: next train in 5 s
            (15.0, 5.0, 0.0),
        ],
    )
    def test_time_till_next(self, macro, time, till_on, till_off):
        assert macro.time_till_next_ON(time) == pytest.approx(till_on)
        assert macro.time_till_next_OFF(time) == pytest.approx(till_off)

    def test_time_till_next_train(self, macro):
        assert macro.time_till_next_train(4.0) == pytest.approx(6.0)
        assert macro.time_till_next_train(14.0) == pytest.approx(6.0)

    def test_start_time(self):
        macro = MacroTimer.new_macropulse(
            macro_period=10,
            macro_on_time=3,
            micro_period=1,
            micro_on_time=0.25,
            macro_start_time=5.0,
        )
        assert macro.state_at(4.0) is OFF  # before start
        assert macro.state_at(5.1) is ON
        assert macro.state_at(5.5) is OFF
        assert macro.state_at(8.5) is OFF  # train is over
        assert macro.state_at(15.1) is ON  # second train
        assert macro.period_index(14.9) == 0
        assert macro.period_index(15.0) == 1

    def test_micro_pulses_per_train_exact_multiple(self, macro):
        assert macro.micro_pulses_per_train() == 3.0

    def test_micro_pulses_per_train_remainder(self):
        # Train of 3.5 s: 3 full pulses + remainder 0.5 s > on_time -> 4
        macro = MacroTimer.new_macropulse(10, 3.5, 1, 0.25)
        assert macro.micro_pulses_per_train() == 4.0

    def test_micro_pulses_per_train_partial_remainder(self):
        # 3 + 0.1 / 0.25 = 3.4 pulses
        macro = MacroTimer.new_macropulse(10, 3.1, 1, 0.25)
        assert macro.micro_pulses_per_train() == pytest.approx(3.4)

    def test_micro_pulses_per_train_half_pulse(self):
        # train 0.205 s, micro period 0.1 s, micro on_time 0.01 s
        macro = MacroTimer.new_macropulse(1.0, 0.205, 0.1, 0.01)
        assert macro.micro_pulses_per_train() == pytest.approx(2.5)
        assert macro.duty_cycle() == pytest.approx(0.025)

    def test_micro_pulses_per_train_float_noise(self):
        # 0.001 // 0.0001 is 9.0 in floating point
        macro = MacroTimer.new_macropulse(0.1, 1e-3, 1e-4, 2e-6)
        assert macro.micro_pulses_per_train() == 10.0

    def test_micro_pulses_per_train_micro_start(self):
        # micro pulses at 0.5, 1.5, 2.5 (train of 3 s)
        macro = MacroTimer.new_macropulse(10, 3, 1, 0.25, micro_start_time=0.5)
        assert macro.micro_pulses_per_train() == 3.0
        # ... and only 0.1 / 0.25 of the pulse at 2.5 in a train of 2.6 s
        macro = MacroTimer.new_macropulse(10, 2.6, 1, 0.25, micro_start_time=0.5)
        assert macro.micro_pulses_per_train() == pytest.approx(2.4)

    def test_duty_cycle(self, macro):
        # 3 pulses of 0.25 s over a 10 s period
        assert macro.duty_cycle() == pytest.approx(0.075)

    def test_duty_cycle_partial_remainder(self):
        # (3 * 0.25 + 0.1) / 10
        macro = MacroTimer.new_macropulse(10, 3.1, 1, 0.25)
        assert macro.duty_cycle() == pytest.approx(0.085)

    def test_time_till_next_ON_with_start_time(self):
        # trains at 5, 15, ...: [5, 8], micro pulses at 5, 6, 7 (0.25 s)
        macro = MacroTimer.new_macropulse(
            macro_period=10, macro_on_time=3, micro_period=1, micro_on_time=0.25,
            macro_start_time=5.0,
        )
        assert macro.time_till_next_ON(2.0) == pytest.approx(3.0)  # before start
        assert macro.time_till_next_ON(5.1) == 0.0  # ON
        assert macro.time_till_next_ON(5.5) == pytest.approx(0.5)  # next micro
        assert macro.time_till_next_ON(7.5) == pytest.approx(7.5)  # next train
        assert macro.time_till_next_ON(12.0) == pytest.approx(3.0)
        assert macro.time_till_next_train(12.0) == pytest.approx(3.0)
        assert macro.time_till_next_train(2.0) == pytest.approx(3.0)

    def test_time_till_next_ON_micro_start(self):
        macro = MacroTimer.new_macropulse(10, 3, 1, 0.25, micro_start_time=0.5)
        assert macro.time_till_next_ON(0.0) == pytest.approx(0.5)
        assert macro.time_till_next_ON(9.0) == pytest.approx(1.5)

    def test_time_till_next_OFF_cut_by_train_end(self):
        # micro pulse [2.9, 3.4] cut by the end of the train at 3.0
        macro = MacroTimer.new_macropulse(10, 3, 1, 0.5, micro_start_time=0.9)
        assert macro.time_till_next_OFF(2.95) == pytest.approx(0.05)

    def test_copy_keeps_start_time(self):
        macro = MacroTimer.new_macropulse(
            10, 3, 1, 0.25, macro_start_time=5.0, micro_start_time=0.0
        )
        copy = macro.copy()
        assert copy is not macro
        assert copy.subtimer is not macro.subtimer
        assert (copy.period, copy.on_time) == (10, 3)
        assert copy.start_time == 5.0
        assert macro.copy(start_time=7.0).start_time == 7.0
