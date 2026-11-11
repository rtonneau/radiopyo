#!/usr/bin/env python3
"""
Benchmark of the Fenton reproduction grid (LSODA stalls, ADR 0007).

Chemistry: config_Labarbe + Fe2+ fixed (2.2e-8 M), O2 50 umol/L, RH = 1 M
(reactions_biology); 10 Gy delivered as pulses of dose D every 0.1 s, each
at the peak dose rate DR (LET 0.201 keV/um), then 2 s of relaxation. Grid:
D in {0.1, 1, 2} Gy x DR in {1e3, 1e5, 1e7} Gy/s (points whose pulse is shorter than
UnitCell.MIN_BEAM_ON_TIME = 1e-6 s, i.e. every point at 1e7 Gy/s, are skipped), plus Fe2+ = 8.9e-7 M at
D = 1 Gy and DR in {1e3, 1e5}. Before the fix some of these points stalled
for tens of minutes.

Every point must succeed in less than TIME_LIMIT seconds with the default
solver options. Where a BDF run with atol 1e-10 (and an LSODA run with atol
1e-10) succeeds in time, the final ROOH must agree within 1e-6 (relative).

Usage:
    python examples/bench_fenton_stall.py [--no-reference]

Exit status 0 when every point is ok, 1 otherwise.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 11/11/2026
"""
# ruff: noqa: T201
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import argparse
import dataclasses as dc
import logging
import sys
import tempfile
import time
import typing as tp
from pathlib import Path

# ----------------------------- Internal imports ---------------------------- #
from radiopyo.physics.beam import PulsedBeam
from radiopyo.simulation.unit_cell import UnitCell

# --------------------------------------------------------------------------- #
#                                 Expose stuff                                #
# --------------------------------------------------------------------------- #
__all__ = ["BenchPoint", "main", "run_point"]

# --------------------------------------------------------------------------- #
#                                   Logging                                   #
# --------------------------------------------------------------------------- #
logging.getLogger("radiopyo").setLevel(logging.WARNING)

# --------------------------------------------------------------------------- #
#                                  Constants                                  #
# --------------------------------------------------------------------------- #
TIME_LIMIT = 60.0  # s, per grid point
ROOH_RTOL = 1e-6
LET = 0.201  # keV/um
TOTAL_DOSE = 10.0  # Gy
PERIOD = 0.1  # s
RELAX_TIME = 2.0  # s
DOSES = (0.1, 1.0, 2.0)
DOSE_RATES = (1e3, 1e5, 1e7)

CONFIG_TEMPLATE = """\
[includes]
toml = ["radiopyo/config_Labarbe.toml"]

[concentrations.fixed]
Fe_plus_plus = {fe}
RH = 1

[concentrations.initial]
O2 = 50e-6
"""

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


@dc.dataclass(frozen=True)
class BenchPoint:
    """One point of the grid."""

    pulse_dose: float  # Gy
    dose_rate: float  # Gy/s
    fe: float = 2.2e-8  # mol/L

    @property
    def on_time(self) -> float:
        return self.pulse_dose / self.dose_rate

    def __str__(self) -> str:
        return f"D={self.pulse_dose:g} Gy DR={self.dose_rate:.0e} Gy/s Fe={self.fe:g}"


def grid() -> tp.List[BenchPoint]:
    """The reproduction grid (pulses shorter than MIN_BEAM_ON_TIME are left out)."""
    points = [
        BenchPoint(d, dr)
        for d in DOSES
        for dr in DOSE_RATES
        if d / dr >= UnitCell.MIN_BEAM_ON_TIME
    ]
    points += [BenchPoint(1.0, dr, fe=8.9e-7) for dr in (1e3, 1e5)]
    return points


def run_point(
    point: BenchPoint, config_dir: Path, **solver_options: tp.Any
) -> tp.Tuple[bool, float, float]:
    """Run one point. Returns (success, wall time [s], final ROOH [mol/L])."""
    cfg = config_dir / f"fe_{point.fe:g}.toml"
    if not cfg.exists():
        cfg.write_text(CONFIG_TEMPLATE.format(fe=point.fe))
    cell = UnitCell.from_toml(str(cfg))
    n_pulse = TOTAL_DOSE / point.pulse_dose
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=point.pulse_dose,
        period=PERIOD,
        on_time=point.on_time,
        n_pulse=n_pulse,
    )
    cell.set_beam(beam, LET=LET)
    t0 = time.perf_counter()
    try:
        res = cell.prepare_chunked_run(
            t_span=(0.0, n_pulse * PERIOD + RELAX_TIME),
            max_step_size_on=point.on_time / 10,
            max_step_size_off=PERIOD / 10,
            **solver_options,
        ).run()
    except Exception as err:  # noqa: BLE001  SolverFailedError & co: report, do not abort
        print(f"    failed: {type(err).__name__}: {str(err)[:100]}")
        return False, time.perf_counter() - t0, float("nan")
    wall = time.perf_counter() - t0
    return bool(all(res.success)), wall, float(res.to_pandas()["ROOH"].iloc[-1])


def main(argv: tp.Optional[tp.Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--no-reference",
        action="store_true",
        help="skip the BDF / LSODA atol=1e-10 reference runs",
    )
    args = parser.parse_args(argv)

    n_bad = 0
    with tempfile.TemporaryDirectory() as tmp:
        config_dir = Path(tmp)
        for point in grid():
            ok, wall, rooh = run_point(point, config_dir)
            good = ok and wall < TIME_LIMIT
            msg = f"{point}: ok={ok} wall={wall:.1f}s ROOH={rooh:.6e}"
            if good and not args.no_reference:
                for name, opts in (
                    ("BDF", {"method": "BDF", "atol": 1e-10}),
                    ("LSODA", {"atol": 1e-10}),
                ):
                    r_ok, r_wall, r_rooh = run_point(point, config_dir, **opts)
                    if r_ok and r_wall < TIME_LIMIT:
                        rel = abs(rooh - r_rooh) / max(abs(r_rooh), 1e-300)
                        close = rel <= ROOH_RTOL or rooh == r_rooh
                        msg += f" | {name} ref rel={rel:.1e} ({r_wall:.1f}s)"
                        good = good and close
                    else:
                        msg += f" | {name} ref n/a"
            if not good:
                n_bad += 1
                msg += "  <-- FAIL"
            print(msg, flush=True)
    print(f"{n_bad} failing point(s)")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
