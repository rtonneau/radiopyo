#!/usr/bin/env python3
"""
Rust backend batch API (radiopyo.rust.run_batch, ADR 0010): parallel runs
equal the single-cell Rust runs, input order is kept, errors are raised
after all cells finish or returned in their slots. Skipped without
radiopyo-rs.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import os
import time

import numpy as np
import pytest

pytest.importorskip("radiopyo_rs")

# ----------------------------- Internal imports ---------------------------- #
from radiopyo import rust
from radiopyo.physics.beam import ConstantBeam, PulsedBeam
from radiopyo.simulation.exceptions import NoBeamDefinedError, ParityError
from radiopyo.simulation.ode_result import ODEResult
from radiopyo.simulation.unit_cell import UnitCell

from .conftest import template

T_SPAN = (0.0, 3e-3)


def _cells(n: int = 4, config: str = "pure_water") -> list[UnitCell]:
    cells = []
    for i in range(n):
        cell = template(config)
        cell.set_beam(
            PulsedBeam.from_dose_rate_peak(1e6 * (i + 1), 1e-3, 1e-5, max_dose=50.0),
            LET=0.201,
        )
        cells.append(cell)
    return cells


def _same(a: ODEResult, b: ODEResult) -> None:
    np.testing.assert_array_equal(a.t, b.t)
    np.testing.assert_array_equal(a.y, b.y)


@pytest.mark.parametrize("mode", ["chunked", "single-pulse"])
def test_batch_equals_single_runs(mode):
    cells = _cells()
    res = rust.run_batch(cells, T_SPAN, mode=mode, n_threads=2)
    prep = "prepare_chunked_run" if mode == "chunked" else "prepare_single_pulse"
    for cell, r in zip(cells, res, strict=True):
        assert isinstance(r, ODEResult)
        _same(r, getattr(cell, prep)(t_span=T_SPAN).run(backend="rust"))


def test_plain_run_mode_and_options():
    cells = [template("labarbe") for _ in range(2)]
    for cell in cells:
        cell.set_beam(ConstantBeam(dose_rate=0.1, max_dose=1e-3), LET=0.201)
    t_eval = np.linspace(0.0, 0.1, 11)
    res = rust.run_batch(cells, (0.0, 0.1), t_eval=t_eval, rtol=1e-6)
    single = cells[0].run(t_span=(0.0, 0.1), t_eval=t_eval, rtol=1e-6, backend="rust")
    _same(res[0], single)
    np.testing.assert_allclose(res[1].t, single.t)


def test_order_and_mixed_configs():
    cells = _cells(2) + _cells(2, config="labarbe")
    res = rust.run_batch(cells, T_SPAN, mode="chunked")
    sizes = [r.y.shape[0] for r in res]
    assert sizes[0] == sizes[1] != sizes[2] == sizes[3]
    # higher peak dose rate -> more e_aq at the end, slot by slot
    idx = cells[0].env.species["e_aq"].index
    assert res[0].y[idx, -1] != res[1].y[idx, -1]


def _with_failures() -> list[UnitCell]:
    cells = _cells(4)
    cells[1] = template("pure_water")  # no beam
    # k edited after the build: radkin would solve another system
    reaction = next(iter(cells[2].env.reactions.k_reaction))
    reaction.set_k_value(reaction.k_value(cells[2]._T) * 2.0)
    cells[2].update_env()
    return cells


def test_failures_raise_after_all():
    cells = _with_failures()
    with pytest.raises(NoBeamDefinedError) as info:
        rust.run_batch(cells, T_SPAN, mode="chunked")
    assert any("cell 1" in n for n in info.value.__notes__)


def test_failures_returned_in_their_slots():
    cells = _with_failures()
    res = rust.run_batch(cells, T_SPAN, mode="chunked", return_exceptions=True)
    assert isinstance(res[0], ODEResult)
    assert isinstance(res[1], NoBeamDefinedError)
    assert isinstance(res[2], ParityError)
    assert isinstance(res[3], ODEResult)


def test_bad_inputs():
    with pytest.raises(ValueError, match="mode"):
        rust.run_batch(_cells(1), T_SPAN, mode="other")
    with pytest.raises(ValueError, match="unknown"):
        rust.run_batch(_cells(1), T_SPAN, bogus=1)


@pytest.mark.skipif((os.cpu_count() or 1) < 2, reason="needs 2 cores")
def test_batch_is_faster_than_serial():
    cells = [template("labarbe") for _ in range(4)]
    for cell in cells:
        cell.set_beam(
            PulsedBeam.from_dose_rate_peak(1e6, 1e-3, 1e-5, max_dose=50.0), LET=0.201
        )
    span = (0.0, 2e-2)
    rust.run_batch(cells[:1], span, mode="chunked")  # warm up
    t0 = time.perf_counter()
    for cell in cells:
        cell.prepare_chunked_run(t_span=span).run(backend="rust")
    serial = time.perf_counter() - t0
    t0 = time.perf_counter()
    rust.run_batch(cells, span, mode="chunked")
    batch = time.perf_counter() - t0
    assert batch < serial
