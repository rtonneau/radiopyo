#!/usr/bin/env python3
"""
Rust backend (radiopyo.rust, ADR 0010), chunked and single-pulse modes:
SimChunksIter.run(backend="rust") gives the same results as the Python
backend (chunks, thinning, short ON spans, reconstructed rates). Skipped
without radiopyo-rs.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("radiopyo_rs")

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics.beam import MacroPulsedBeam, PulsedBeam, SinglePulseBeam
from radiopyo.simulation.exceptions import RustBackendError
from radiopyo.simulation.unit_cell import SimChunksIter

from .conftest import template

RTOL = 1e-6
ATOL = 1e-10
LET = 0.201


def _beam(name: str):
    if name == "pulsed":
        return PulsedBeam.from_dose_rate_peak(1e6, 1e-3, 1e-5, max_dose=50.0)
    if name == "macro":
        return MacroPulsedBeam(
            dose_rate_peak=1e6,
            period=1e-3,
            on_time=2e-5,
            micro_period=1e-5,
            micro_on_time=5e-6,
        )
    return SinglePulseBeam(dose_rate=1e6, max_dose=10.0)


def _pair(config, beam, prepare, **kw):
    cell = template(config)
    cell.set_beam(_beam(beam), LET=LET)
    py = getattr(cell, prepare)(**kw).run()
    rs = getattr(cell, prepare)(**kw).run(backend="rust")
    return py, rs


def _same(py, rs) -> None:
    """Same final state and time span (the two LSODA codes choose their
    own steps, so the stored grids differ unless thinned)."""
    assert rs.t[0] == pytest.approx(py.t[0])
    assert rs.t[-1] == pytest.approx(py.t[-1])
    np.testing.assert_allclose(rs.final_cc, py.final_cc, rtol=RTOL, atol=ATOL)


@pytest.mark.parametrize("config", ["pure_water", "labarbe"])
@pytest.mark.parametrize("beam", ["pulsed", "macro"])
def test_chunked_matches_python(config, beam):
    py, rs = _pair(config, beam, "prepare_chunked_run", t_span=(0.0, 3e-3))
    _same(py, rs)

    np.testing.assert_allclose(rs.step_dose_rates, py.step_dose_rates)


@pytest.mark.parametrize("config", ["pure_water", "labarbe"])
def test_single_pulse_matches_python(config):
    py, rs = _pair(config, "single", "prepare_single_pulse", t_span=(0.0, 1.0))
    _same(py, rs)


def test_thinning_matches_python():
    py, rs = _pair(
        "pure_water",
        "pulsed",
        "prepare_chunked_run",
        t_span=(0.0, 3e-3),
        points_per_chunk_on=5,
        points_per_chunk_off=7,
    )
    assert rs.thinned and py.thinned
    np.testing.assert_allclose(rs.t, py.t, rtol=1e-12)
    np.testing.assert_allclose(rs.y, py.y, rtol=RTOL, atol=ATOL)


def test_short_on_span_from_max_dose_cut():
    cell = template("pure_water")
    cell.set_beam(
        PulsedBeam.from_dose_rate_peak(1e6, 1e-3, 1e-5, max_dose=10.1), LET=LET
    )
    py = cell.prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    rs = cell.prepare_chunked_run(t_span=(0.0, 3e-3)).run(backend="rust")
    assert len(rs.short_on_spans) == len(py.short_on_spans) == 1
    a, b = rs.short_on_spans[0], py.short_on_spans[0]
    assert a.start == pytest.approx(b.start)
    assert a.stop == pytest.approx(b.stop)
    assert a.dose_rate == pytest.approx(b.dose_rate)
    _same(py, rs)


def test_reconstruct_reactions_matches_python():
    py, rs = _pair(
        "labarbe",
        "pulsed",
        "prepare_chunked_run",
        t_span=(0.0, 3e-3),
        points_per_chunk_on=50,
        points_per_chunk_off=50,
    )
    dpy = py.reconstruct_reactions()
    drs = rs.reconstruct_reactions()
    assert list(drs.columns) == list(dpy.columns)
    np.testing.assert_allclose(
        drs.to_numpy(dtype=float), dpy.to_numpy(dtype=float), rtol=RTOL, atol=1e-8
    )


def test_hand_built_iterator_refuses_rust():
    cell = template("pure_water")
    cell.set_beam(_beam("pulsed"), LET=LET)
    it = cell.prepare_chunked_run(t_span=(0.0, 3e-3))
    with pytest.raises(RustBackendError, match="prepare_chunked_run"):
        SimChunksIter(cell, it.chunks).run(backend="rust")
    next(it)
    with pytest.raises(RustBackendError, match="fresh"):
        it.run(backend="rust")
