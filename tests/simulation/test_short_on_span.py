#!/usr/bin/env python3
"""
Tests of the short ON spans: a beam ON span shorter than MIN_BEAM_ON_TIME
at run time (last pulse cut by max_dose, pulse clipped by t_span) is not
integrated; its radiolytic yield is added analytically, a warning is
logged and the span is recorded in ODEResult.short_on_spans.

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
from radiopyo.physics.beam import MIN_BEAM_ON_TIME, PulsedBeam
from radiopyo.simulation.ode_result import ShortOnSpan
from radiopyo.simulation.unit_cell import UnitCell, is_short_on_span

PEAK = 1e6
PERIOD = 1e-3
ON_TIME = 1e-5


@pytest.fixture(scope="module")
def base_cell() -> UnitCell:
    return UnitCell.from_toml(rp.PURE_WATER_CONFIG_FILE)


def _cell(base: UnitCell, max_dose: float) -> UnitCell:
    cell = copy.deepcopy(base)
    beam = PulsedBeam.from_dose_rate_peak(PEAK, PERIOD, ON_TIME, max_dose=max_dose)
    cell.set_beam(beam, LET=10.0)
    return cell


# ------------------------------ is_short_on_span ---------------------------- #
def test_is_short_on_span_threshold():
    assert is_short_on_span((0.0, 3e-7))
    assert not is_short_on_span((0.0, MIN_BEAM_ON_TIME))
    # rounding of a chunk boundary computed at a large start time
    start = 1234.5
    assert not is_short_on_span((start, start + MIN_BEAM_ON_TIME))
    assert is_short_on_span((start, start + 0.5 * MIN_BEAM_ON_TIME))


# ------------------------------ chunked runs ------------------------------- #
@pytest.mark.parametrize("max_dose", [10.1, 0.1], ids=["last_pulse", "first_pulse"])
def test_partial_pulse_is_a_short_on_span(base_cell, caplog, max_dose):
    cell = _cell(base_cell, max_dose)
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        res = cell.prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    assert len(res.short_on_spans) == 1
    span = res.short_on_spans[0]
    assert isinstance(span, ShortOnSpan)
    assert span.stop - span.start == pytest.approx(1e-7, rel=1e-6)
    assert span.dose_rate == pytest.approx(PEAK)
    assert any("shorter than MIN_BEAM_ON_TIME" in m for m in caplog.messages)
    assert res.t[-1] == pytest.approx(3e-3)


@pytest.mark.parametrize("max_dose", [10.1, 10.3, 0.1])
def test_production_integral_counts_the_short_span_once(base_cell, max_dose):
    cell = _cell(base_cell, max_dose)
    res = cell.prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    prod = res.integrate_production_of("e_aq")
    (r,) = cell.env.radiolytic_reaction_involving("e_aq")
    expected = r.compute_initial_cc_before(max_dose / PEAK, PEAK) * 1e6
    assert prod[r.as_label()] == pytest.approx(expected, rel=1e-9)


def test_full_pulses_have_no_short_on_span(base_cell):
    res = _cell(base_cell, 20.0).prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    assert res.short_on_spans == []


def test_reconstruct_and_pandas_work_with_a_short_span(base_cell):
    res = _cell(base_cell, 10.1).prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    df = res.reconstruct_reactions()
    assert len(df) == len(res.t)
    assert np.all(np.isfinite(df.to_numpy(dtype=float)))
    assert len(res.to_pandas()) == len(res.t)


def test_span_clipped_by_t_span_is_a_short_on_span(base_cell):
    cell = _cell(base_cell, 50.0)
    res = cell.prepare_chunked_run(t_span=(9.7e-6, 2e-3)).run()
    assert len(res.short_on_spans) == 1
    assert res.short_on_spans[0].start == pytest.approx(9.7e-6)


def test_plain_run_short_span(base_cell, caplog):
    cell = _cell(base_cell, 50.0)
    with caplog.at_level(logging.WARNING, logger="radiopyo"):
        res = cell.run(t_span=(0.0, 3e-7))
    assert len(res.short_on_spans) == 1
    assert res.t[-1] == pytest.approx(3e-7)
    assert not [dr for dr in (res.burst_dose_rates or []) if dr is not None]
    (r,) = cell.env.radiolytic_reaction_involving("e_aq")
    idx = r.species().index
    gained = res.y[idx, -1] - cell.env.initial_values()[idx]
    expected = r.compute_initial_cc_before(3e-7, PEAK) * 1e6
    assert gained == pytest.approx(expected, rel=1e-9)
    assert any("shorter than MIN_BEAM_ON_TIME" in m for m in caplog.messages)


def test_single_pulse_clipped_by_t_span(base_cell):
    cell = _cell(base_cell, 50.0)
    res = cell.prepare_single_pulse(t_span=(0.0, PERIOD + 3e-7)).run()
    assert len(res.short_on_spans) == 1
    assert res.short_on_spans[0].start == pytest.approx(PERIOD)


def test_add_concatenates_short_on_spans(base_cell):
    cell = _cell(base_cell, 50.0)
    first = cell.run(t_span=(0.0, 3e-7))
    second = cell.run(t_span=(3e-7, 4e-7), y0=first.y[:, -1])
    total = first + second
    assert total.short_on_spans == [*first.short_on_spans, *second.short_on_spans]
    assert len(total.short_on_spans) == 2


def test_micro_pulse_cut_by_its_train_end_is_a_short_on_span(base_cell):
    from radiopyo.physics.beam import MacroPulsedBeam

    cell = copy.deepcopy(base_cell)
    beam = MacroPulsedBeam(
        dose_rate_peak=PEAK,
        period=1e-3,
        on_time=2.03e-5,
        micro_period=1e-5,
        micro_on_time=5e-6,
    )
    cell.set_beam(beam, LET=10.0)
    res = cell.prepare_chunked_run(t_span=(0.0, 5e-4)).run()
    assert len(res.short_on_spans) == 1
    span = res.short_on_spans[0]
    assert span.start == pytest.approx(2e-5)
    assert span.stop == pytest.approx(2.03e-5)


def test_result_loaded_from_hdf5_has_no_short_on_span(base_cell, tmp_path):
    import dataclasses as dc

    from radiopyo.simulation.persistence import H5Manager, SimIdentifier

    @dc.dataclass(kw_only=True)
    class _ID(SimIdentifier):
        case: str = "short"

    res = _cell(base_cell, 10.1).prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    assert res.short_on_spans
    hdf = H5Manager.from_path(tmp_path / "short.h5")
    hdf.push_ODEResult(_ID(), res)
    loaded = hdf.read_ODEResult(_ID())
    assert loaded.short_on_spans == []
    np.testing.assert_allclose(loaded.y[:, -1], res.y[:, -1], rtol=1e-12)
