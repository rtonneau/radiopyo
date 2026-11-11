#!/usr/bin/env python3
"""
Regression test of the LSODA stall of the Fenton reproduction (ADR 0007).

A pulsed beam of 0.1 Gy pulses at 1e5 Gy/s on Labarbe + Fe2+ (2.2e-8 M) +
RH (1 M) used to hold LSODA at steps of ~1e-8 s for minutes (about 80 s for
the three pulses below, 1 s now). The full grid is in
``examples/bench_fenton_stall.py``.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 11/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import time
import typing as tp
from pathlib import Path

import pytest

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics.beam import PulsedBeam
from radiopyo.simulation.unit_cell import UnitCell

# --------------------------------------------------------------------------- #
#                                  Constants                                  #
# --------------------------------------------------------------------------- #
CONFIG = """\
[includes]
toml = ["radiopyo/config_Labarbe.toml"]

[concentrations.fixed]
Fe_plus_plus = 2.2e-8
RH = 1

[concentrations.initial]
O2 = 50e-6
"""
PULSE_DOSE = 0.1  # Gy
DOSE_RATE = 1e5  # Gy/s
PERIOD = 0.1  # s
N_PULSE = 3
T_END = N_PULSE * PERIOD + 0.3
WALL_LIMIT = 30.0  # s (about 1 s now; the stall took more than 80 s)

# --------------------------------------------------------------------------- #
#                               FUNCTION / TESTS                              #
# --------------------------------------------------------------------------- #


def _run(config: Path, **solver_options: tp.Any) -> tp.Tuple[bool, float]:
    """Chunked run of the stalling point: (success, final ROOH [mol/L])."""
    cell = UnitCell.from_toml(str(config))
    on_time = PULSE_DOSE / DOSE_RATE
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=PULSE_DOSE, period=PERIOD, on_time=on_time, n_pulse=N_PULSE
    )
    cell.set_beam(beam, LET=0.201)
    res = cell.prepare_chunked_run(
        t_span=(0.0, T_END),
        max_step_size_on=on_time / 10,
        max_step_size_off=PERIOD / 10,
        **solver_options,
    ).run()
    return bool(all(res.success)), float(res.to_pandas()["ROOH"].iloc[-1])


def test_fenton_pulses_do_not_stall(tmp_path: Path) -> None:
    config = tmp_path / "fenton.toml"
    config.write_text(CONFIG)

    t0 = time.perf_counter()
    ok, rooh = _run(config)
    wall = time.perf_counter() - t0

    assert ok
    assert wall < WALL_LIMIT, f"run took {wall:.1f} s (stall?)"
    assert rooh > 0.0
    # Same answer with a tighter absolute tolerance
    ok_ref, rooh_ref = _run(config, atol=1e-10)
    assert ok_ref
    assert rooh == pytest.approx(rooh_ref, rel=1e-6)
