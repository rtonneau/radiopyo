#!/usr/bin/env python3
"""
Tests of the opt-in point thinning of the chunked runs
(points_per_chunk_on / points_per_chunk_off): N log-spaced points stored
per chunk, the solver steps (hence the chained state) unchanged, the
default run untouched, and one warning when integrals are read from a
thinned result.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import logging

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.physics.beam import PulsedBeam, SinglePulseBeam
from radiopyo.simulation.unit_cell import (
    UnitCell,
    check_points_per_chunk,
    thinned_times,
)

T_SPAN = (0.0, 3e-4)


# ------------------------------- Fixtures ---------------------------------- #
@pytest.fixture(scope="module")
def _template() -> UnitCell:
    cell = UnitCell.from_toml(rp.PURE_WATER_CONFIG_FILE)
    cell.set_beam(
        PulsedBeam.from_dose_per_pulse(
            pulse_dose=1.0, period=1e-4, on_time=2e-6, n_pulse=3.0
        ),
        LET=10.0,
    )
    return cell


@pytest.fixture
def cell(_template: UnitCell) -> UnitCell:
    return copy.deepcopy(_template)


@pytest.fixture(scope="module")
def full_run(_template: UnitCell) -> tuple:
    cell = copy.deepcopy(_template)
    it = cell.prepare_chunked_run(t_span=T_SPAN)
    return it, it.run()


# --------------------------------- Tests ----------------------------------- #
def test_thinned_times_log_spaced() -> None:
    t = thinned_times(1.0, 2.0, 10)
    assert t.size == 10
    assert t[0] == 1.0 and t[-1] == 2.0
    off = np.diff(t[1:])
    assert np.all(off > 0)
    np.testing.assert_allclose(off[1:] / off[:-1], off[1] / off[0], rtol=1e-6)
    assert t[1] - 1.0 == pytest.approx(1e-6, rel=1e-6)
    np.testing.assert_array_equal(thinned_times(1.0, 1.0, 5), [1.0])
    assert thinned_times(0.0, 1.0, 2).tolist() == [0.0, 1.0]


@pytest.mark.parametrize("bad", [0, 1, -3, 2.5, "4", True, False])
def test_check_points_per_chunk_rejects(bad: object) -> None:
    with pytest.raises(ValueError):
        check_points_per_chunk(bad, "x")  # type: ignore[arg-type]


def test_check_points_per_chunk_accepts() -> None:
    check_points_per_chunk(None, "x")
    check_points_per_chunk(2, "x")
    check_points_per_chunk(np.int64(7), "x")


def test_prepare_rejects_bad_value(cell: UnitCell) -> None:
    with pytest.raises(ValueError, match="points_per_chunk_on"):
        cell.prepare_chunked_run(t_span=T_SPAN, points_per_chunk_on=1)
    with pytest.raises(ValueError, match="points_per_chunk_off"):
        cell.prepare_single_pulse(t_span=(0.0, 1e-3), points_per_chunk_off=0)


def test_exact_points_per_chunk(cell: UnitCell) -> None:
    it = cell.prepare_chunked_run(
        t_span=T_SPAN, points_per_chunk_on=7, points_per_chunk_off=12
    )
    res = it.run()
    assert len(res.n_points) == len(it.chunks)
    for n, chunk in zip(res.n_points, it.chunks):
        assert n == (7 if chunk.beam_on else 12)
        assert chunk.n_points == n
    assert res.thinned
    assert res.t[-1] == pytest.approx(T_SPAN[1])


def test_only_off_thinned(cell: UnitCell, full_run: tuple) -> None:
    it_full, res_full = full_run
    res = cell.prepare_chunked_run(t_span=T_SPAN, points_per_chunk_off=5).run()
    for n, n_full, chunk in zip(res.n_points, res_full.n_points, it_full.chunks):
        if chunk.beam_on:
            assert n == n_full
        else:
            assert n == 5


def test_default_unchanged(cell: UnitCell, full_run: tuple) -> None:
    _, res_full = full_run
    res = cell.prepare_chunked_run(t_span=T_SPAN).run()
    assert not res.thinned
    np.testing.assert_array_equal(res.t, res_full.t)
    np.testing.assert_array_equal(res.y, res_full.y)


def test_chunk_ends_match_full_run(cell: UnitCell, full_run: tuple) -> None:
    _, res_full = full_run
    res = cell.prepare_chunked_run(
        t_span=T_SPAN, points_per_chunk_on=4, points_per_chunk_off=4
    ).run()
    ends_full = np.cumsum(res_full.n_points) - 1
    ends = np.cumsum(res.n_points) - 1
    np.testing.assert_allclose(res.t[ends], res_full.t[ends_full], rtol=1e-12)
    # The state at the chunk ends (hence the chained y0) is the same up to
    # the dense-output interpolation of the last step.
    scale = np.max(np.abs(res_full.y), axis=1, keepdims=True) + 1e-30
    np.testing.assert_allclose(
        res.y[:, ends] / scale, res_full.y[:, ends_full] / scale, atol=1e-6
    )


def test_single_pulse_thinning(cell: UnitCell) -> None:
    cell.set_beam(
        SinglePulseBeam.from_on_time(dose_rate=5e5, on_time=2e-6, start_time=0.0),
        LET=10.0,
    )
    it = cell.prepare_single_pulse(
        t_span=(0.0, 1e-3), points_per_chunk_on=6, points_per_chunk_off=9
    )
    res = it.run()
    for n, chunk in zip(res.n_points, it.chunks):
        assert n == (6 if chunk.beam_on else 9)
    assert res.thinned


def test_integrals_warn_once(cell: UnitCell, caplog: pytest.LogCaptureFixture) -> None:
    res = cell.prepare_chunked_run(t_span=T_SPAN, points_per_chunk_off=5).run()
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        res.integrate_species()
        res.integrate_species(stop=1e-4)
    msgs = [r for r in caplog.records if "thinn" in r.getMessage().lower()]
    assert len(msgs) == 1


def test_no_warning_on_full_run(
    full_run: tuple, caplog: pytest.LogCaptureFixture
) -> None:
    _, res_full = full_run
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        res_full.integrate_species()
    assert not [r for r in caplog.records if "thinn" in r.getMessage().lower()]


def test_short_on_span_with_thinning(cell: UnitCell) -> None:
    # A max dose cutting the last pulse leaves a short ON span: it is not
    # integrated (analytic yield, no t_eval), the other chunks are thinned.
    peak, period, on_time = 1e6, 1e-3, 1e-5
    cell.set_beam(
        PulsedBeam.from_dose_rate_peak(peak, period, on_time, max_dose=10.1),
        LET=10.0,
    )
    it = cell.prepare_chunked_run(
        t_span=(0.0, 3e-3), points_per_chunk_on=5, points_per_chunk_off=5
    )
    res = it.run()
    assert len(res.short_on_spans) == 1
    assert res.thinned
    assert res.t[-1] == pytest.approx(3e-3)
    for n, chunk in zip(res.n_points, it.chunks):
        if not chunk.beam_on:
            assert n == 5
