#!/usr/bin/env python3
"""
Benchmark of the cost of a chunked pulsed run versus the number of pulses.

Scenario: bundled config_Labarbe, PulsedBeam of peak dose rate 1e4 Gy/s,
period 10 ms, ON time 10 us (0.1 Gy per pulse), run with
``prepare_chunked_run`` and the default solver options; one ON and one OFF
chunk per pulse. Optional point thinning (``--points-on`` /
``--points-off``) stores N log-spaced points per chunk.

For every pulse count it prints the number of chunks, the stored points,
the wall time, the mean cost per chunk, and the mean cost per chunk over
the first and the last 10% of the chunks: with a linear storage, the last
equals the first (no growth with the stored history).

Usage:
    python examples/bench_chunked_scaling.py [--pulses 100 400 1600]
        [--points-on N] [--points-off N]

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 12/11/2026
"""
# ruff: noqa: T201
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import argparse
import logging
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
DEFAULT_PULSES = (100, 400, 1600)

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


def run_once(
    n_pulse: int,
    points_on: tp.Optional[int] = None,
    points_off: tp.Optional[int] = None,
) -> tp.Dict[str, float]:
    """One chunked run of n_pulse pulses. Returns its figures."""
    cell = UnitCell.from_toml(str(LABARBE_CONFIG_FILE))
    beam = PulsedBeam.from_dose_rate_peak(
        peak_dr=PEAK_DOSE_RATE,
        period=PERIOD,
        on_time=ON_TIME,
        max_dose=n_pulse * PEAK_DOSE_RATE * ON_TIME,
    )
    cell.set_beam(beam)
    chunks = cell.prepare_chunked_run(
        t_span=(0.0, n_pulse * PERIOD),
        points_per_chunk_on=points_on,
        points_per_chunk_off=points_off,
    )
    durations: tp.List[float] = []
    t0 = time.perf_counter()
    last = t0
    res = None
    for res in chunks:
        now = time.perf_counter()
        durations.append(now - last)
        last = now
    wall = last - t0
    if res is None or not all(res.success):
        raise RuntimeError("a chunk failed")
    tenth = max(1, len(durations) // 10)
    return {
        "chunks": len(durations),
        "points": res.t.size,
        "wall": wall,
        "per_chunk": wall / len(durations),
        "first": sum(durations[:tenth]) / tenth,
        "last": sum(durations[-tenth:]) / tenth,
    }


def main(argv: tp.Optional[tp.Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--pulses", type=int, nargs="+", default=list(DEFAULT_PULSES),
        help="pulse counts to run (default: 100 400 1600)",
    )
    parser.add_argument("--points-on", type=int, default=None,
                        help="stored points per ON chunk (default: every step)")
    parser.add_argument("--points-off", type=int, default=None,
                        help="stored points per OFF chunk (default: every step)")
    args = parser.parse_args(argv)

    print(f"config: {LABARBE_CONFIG_FILE}")
    print(f"thinning: on={args.points_on} off={args.points_off}")
    print(f"{'pulses':>7} {'chunks':>7} {'points':>9} {'wall [s]':>9} "
          f"{'ms/chunk':>9} {'first 10%':>10} {'last 10%':>9}")
    for n_pulse in args.pulses:
        r = run_once(n_pulse, args.points_on, args.points_off)
        print(f"{n_pulse:7d} {r['chunks']:7.0f} {r['points']:9.0f} "
              f"{r['wall']:9.1f} {1e3 * r['per_chunk']:9.1f} "
              f"{1e3 * r['first']:10.1f} {1e3 * r['last']:9.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
