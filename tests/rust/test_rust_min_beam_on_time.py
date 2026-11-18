#!/usr/bin/env python3
"""
Rust backend (radiopyo.rust, ADR 0010) with a non-default min beam ON time
(rp.set_min_beam_on_time): radkin receives the value, the parity guard
compares it, and the results match the Python backend. Skipped without
radiopyo-rs.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import numpy as np
import pytest

pytest.importorskip("radiopyo_rs")

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo import rust
from radiopyo.physics.beam import PulsedBeam
from radiopyo.physics.exceptions import ShortBeamOnError
from radiopyo.simulation.exceptions import ParityError

from .conftest import template

LET = 0.201


@pytest.fixture(autouse=True)
def _restore_default() -> tp.Iterator[None]:
    yield
    rp.set_min_beam_on_time(None)


def _cell(on_time: float, max_dose: float):
    cell = template("pure_water")
    cell.set_beam(
        PulsedBeam.from_dose_rate_peak(1e6, 1e-3, on_time, max_dose=max_dose),
        LET=LET,
    )
    return cell


# 10.3 Gy: 0.3 us last pulse; 10.1 Gy: 0.1 us (= T_ZERO) last pulse
@pytest.mark.parametrize(
    ("value", "max_dose", "n_short"),
    [(None, 10.3, 1), (2e-7, 10.3, 0), (5e-8, 10.1, 1)],
    ids=["default", "0.2us", "below_T_ZERO"],
)
def test_last_pulse_matches_python(value, max_dose, n_short):
    rp.set_min_beam_on_time(value)
    cell = _cell(1e-5, max_dose)
    py = cell.prepare_chunked_run(t_span=(0.0, 3e-3)).run()
    rs = cell.prepare_chunked_run(t_span=(0.0, 3e-3)).run(backend="rust")
    assert len(py.short_on_spans) == len(rs.short_on_spans) == n_short
    np.testing.assert_allclose(rs.final_cc, py.final_cc, rtol=1e-6, atol=1e-10)


def test_short_configured_pulse_runs_with_lower_threshold():
    rp.set_min_beam_on_time(2e-7)
    cell = _cell(3e-7, 5.0)
    py = cell.prepare_chunked_run(t_span=(0.0, 3e-2)).run()
    rs = cell.prepare_chunked_run(t_span=(0.0, 3e-2)).run(backend="rust")
    assert rs.short_on_spans == py.short_on_spans == []
    np.testing.assert_allclose(rs.final_cc, py.final_cc, rtol=1e-6, atol=1e-10)


def test_parity_guard_compares_the_value():
    cell = _cell(1e-5, 50.0)
    system = dict(rust.build_cell(cell).system())
    assert system["min_beam_on_time"] == rp.min_beam_on_time()
    rust.check_parity(cell, system)
    system["min_beam_on_time"] = 1e-6
    with pytest.raises(ParityError, match="min beam ON time"):
        rust.check_parity(cell, system)


def test_radkin_rejects_with_the_higher_value():
    cell = _cell(1e-5, 50.0)
    rp.set_min_beam_on_time(2e-5)  # set after the Python beam was built
    with pytest.raises(ShortBeamOnError):
        cell.run(t_span=(0.0, 1e-3), backend="rust")
