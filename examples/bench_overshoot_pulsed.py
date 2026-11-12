#!/usr/bin/env python3
"""
Benchmark of the overshoot relaxation path (ADR 0007) in a pulsed chunked run.

Scenario (scout report 2026-11-11): bundled config_Labarbe, PulsedBeam of peak dose rate 1e4 Gy/s,
period 10 ms, ON time 10 us, 10 Gy in total (100 pulses, 200 chunks), run
with ``prepare_chunked_run`` and the default solver options.

Prints the wall time, the number of chunks and, from a cProfile run, the
cumulative time spent in the overshoot path (``_overshoot``) and in the
derivative helpers, with their share of the profiled time.

Usage:
    python examples/bench_overshoot_pulsed.py [--repeat N]

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 12/11/2026
"""
# ruff: noqa: T201
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import argparse
import cProfile
import logging
import pstats
import sys
import time
import typing as tp

# ----------------------------- Internal imports ---------------------------- #
from radiopyo import LABARBE_CONFIG_FILE
from radiopyo.physics.beam import PulsedBeam
from radiopyo.simulation.unit_cell import UnitCell

# --------------------------------------------------------------------------- #
#                                 Expose stuff                                #
# --------------------------------------------------------------------------- #
__all__ = ["main", "run_once"]

# --------------------------------------------------------------------------- #
#                                   Logging                                   #
# --------------------------------------------------------------------------- #
logging.getLogger("radiopyo").setLevel(logging.WARNING)

# --------------------------------------------------------------------------- #
#                                  Constants                                  #
# --------------------------------------------------------------------------- #
PEAK_DOSE_RATE = 1e4  # Gy/s
PERIOD = 10e-3  # s
ON_TIME = 10e-6  # s
TOTAL_DOSE = 10.0  # Gy
PROFILED = (
    "_overshoot",
    "_k_relaxation_terms",
    "_ez_relaxation_terms",
    "depletion_switch",
    "depletion_switch_derivatives",
    "mass_action_derivatives",
    "jacobian",
    "derive",
)

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


def run_once(config: str) -> tp.Tuple[float, int]:
    """One chunked run. Returns (wall time of the run [s], number of chunks)."""
    cell = UnitCell.from_toml(config)
    beam = PulsedBeam.from_dose_rate_peak(
        peak_dr=PEAK_DOSE_RATE, period=PERIOD, on_time=ON_TIME, max_dose=TOTAL_DOSE
    )
    cell.set_beam(beam)
    n_pulse = TOTAL_DOSE / (PEAK_DOSE_RATE * ON_TIME)
    chunks = cell.prepare_chunked_run(t_span=(0.0, n_pulse * PERIOD))
    t0 = time.perf_counter()
    res = chunks.run()
    wall = time.perf_counter() - t0
    if not all(res.success):
        raise RuntimeError("a chunk failed")
    return wall, len(res.success)


def main(argv: tp.Optional[tp.Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--repeat", type=int, default=3, help="timed runs (best kept)")
    args = parser.parse_args(argv)
    config = str(LABARBE_CONFIG_FILE)

    walls = []
    n_chunks = 0
    for _ in range(args.repeat):
        wall, n_chunks = run_once(config)
        walls.append(wall)
    print(f"config: {config}")
    print(f"chunks: {n_chunks}")
    print(f"wall (best of {args.repeat}): {min(walls):.2f} s  all: "
          + " ".join(f"{w:.2f}" for w in walls))

    profiler = cProfile.Profile()
    profiler.enable()
    run_once(config)
    profiler.disable()
    stats = pstats.Stats(profiler)
    total = stats.total_tt  # type: ignore[attr-defined]
    print(f"profiled total: {total:.2f} s")
    for (_file, _line, name), (_cc, ncalls, _tt, cumtime, _callers) in sorted(
        stats.stats.items(),  # type: ignore[attr-defined]
        key=lambda item: -item[1][3],
    ):
        if name in PROFILED and "ode_functions_solver" in _file:
            print(f"  {name:32s} {ncalls:8d} calls {cumtime:7.2f} s "
                  f"({100 * cumtime / total:4.1f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
