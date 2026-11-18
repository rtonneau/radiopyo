#!/usr/bin/env python3
"""
Tests of the global min beam ON time setting (rp.set_min_beam_on_time /
rp.min_beam_on_time) and of its read-only view UnitCell.MIN_BEAM_ON_TIME.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import logging
import math
import typing as tp

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.physics.beam import MIN_BEAM_ON_TIME, T_ZERO, PulsedBeam
from radiopyo.physics.exceptions import ShortBeamOnError
from radiopyo.simulation.unit_cell import UnitCell, is_short_on_span


@pytest.fixture(autouse=True)
def _restore_default() -> tp.Iterator[None]:
    yield
    rp.set_min_beam_on_time(None)


@pytest.fixture(scope="module")
def base_cell() -> UnitCell:
    return UnitCell.from_toml(rp.PURE_WATER_CONFIG_FILE)


def test_default_and_exports() -> None:
    assert rp.min_beam_on_time() == MIN_BEAM_ON_TIME == 0.4e-6
    assert UnitCell.MIN_BEAM_ON_TIME == MIN_BEAM_ON_TIME
    assert UnitCell.T_ZERO == T_ZERO


def test_set_and_reset() -> None:
    rp.set_min_beam_on_time(2e-7)
    assert rp.min_beam_on_time() == 2e-7
    assert UnitCell.MIN_BEAM_ON_TIME == 2e-7
    rp.set_min_beam_on_time(None)
    assert rp.min_beam_on_time() == MIN_BEAM_ON_TIME


@pytest.mark.parametrize("bad", [0.0, -1e-6, math.inf, math.nan])
def test_invalid_values_raise(bad: float) -> None:
    with pytest.raises(ValueError):
        rp.set_min_beam_on_time(bad)
    assert rp.min_beam_on_time() == MIN_BEAM_ON_TIME


def test_below_t_zero_warns(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        rp.set_min_beam_on_time(5e-8)
    assert rp.min_beam_on_time() == 5e-8
    assert any("below T_ZERO" in m for m in caplog.messages)


def test_view_is_read_only(base_cell: UnitCell) -> None:
    with pytest.raises(AttributeError):
        base_cell.MIN_BEAM_ON_TIME = 1e-6  # type: ignore[misc]
    rp.set_min_beam_on_time(1e-6)
    assert base_cell.MIN_BEAM_ON_TIME == 1e-6


def test_beam_check_follows_setting() -> None:
    rp.set_min_beam_on_time(None)
    with pytest.raises(ShortBeamOnError):
        PulsedBeam(dose_rate_peak=1e6, period=1e-3, on_time=2e-7)
    rp.set_min_beam_on_time(1e-7)
    PulsedBeam(dose_rate_peak=1e6, period=1e-3, on_time=2e-7)
    rp.set_min_beam_on_time(1e-5)
    with pytest.raises(ShortBeamOnError, match="set_min_beam_on_time"):
        PulsedBeam(dose_rate_peak=1e6, period=1e-3, on_time=2e-6)


def test_is_short_on_span_default_reads_setting() -> None:
    assert is_short_on_span((0.0, 3e-7))
    rp.set_min_beam_on_time(2e-7)
    assert not is_short_on_span((0.0, 3e-7))


def test_lower_threshold_integrates_short_pulse(base_cell: UnitCell) -> None:
    rp.set_min_beam_on_time(2e-7)
    cell = copy.deepcopy(base_cell)
    cell.set_beam(
        PulsedBeam.from_dose_rate_peak(1e6, 1e-3, 1e-5, max_dose=10.3), LET=10.0
    )
    res = cell.prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    assert res.short_on_spans == []  # 0.3 us last pulse now integrated


def test_spans_up_to_t_zero_stay_analytic(base_cell: UnitCell) -> None:
    rp.set_min_beam_on_time(1e-8)
    cell = copy.deepcopy(base_cell)
    cell.set_beam(
        PulsedBeam.from_dose_rate_peak(1e6, 1e-3, 1e-5, max_dose=10.1), LET=10.0
    )
    res = cell.prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    assert len(res.short_on_spans) == 1  # 0.1 us = T_ZERO: still analytic
