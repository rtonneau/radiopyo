#!/usr/bin/env python3
"""
Runs and analysis reproducing the main findings of Tonneau et al. 2025,
"Decoding the chemistry of the FLASH effect: a physicochemical model of
dose-rate, pH and oxygen-dependent H2O2 production", Phys. Med. Biol. 70
235021 (doi:10.1088/1361-6560/ae1fcc).

Shared by the notebook ``tonneau_2025_pmb.ipynb`` (figures) and by the
golden test ``tests/golden/test_paper_tonneau_2025.py`` (the golden data
``tests/golden/data/paper_tonneau_2025.npz`` are written by
``tests/golden/record_golden.py --only paper``), so that both run the
very same code.

Scenario (paper tables 3 and 5): bundled pure-water config, LET 1.17
keV/µm, dose-rate dependent G-values (``variable_GValues``, the GF(DR) fit
shipped in radiopyo, see ADR 0009), one single pulse per run delivered from
t = 0, run over (0, 500 s) with ``prepare_single_pulse``; maximum step of
pulse/500 with the beam ON and 0.1 s for the final beam OFF chunk; initial
[O2] = pO2 [%] / 100 x 1300 µM (Henry); pH 7; 25 degC.

Units: concentrations in µM, dose rates in Gy/s, doses in Gy, pO2 in %,
net yields in #/100 eV, times in s.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 12/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import collections
import contextlib
import copy
import dataclasses
import functools
import logging
import typing as tp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

# ----------------------------- Internal imports ---------------------------- #
from radiopyo import PURE_WATER_CONFIG_FILE
from radiopyo.chemistry.radiolytic import kr_to_ge
from radiopyo.physics.beam import SinglePulseBeam
from radiopyo.simulation.unit_cell import UnitCell

# ----------------------------- Type checking ONLY -------------------------- #
if tp.TYPE_CHECKING:
    from radiopyo.simulation.ode_result import ODEResult

# --------------------------------------------------------------------------- #
#                                 Expose stuff                                #
# --------------------------------------------------------------------------- #
__all__ = [
    "BAR_PO2",
    "DOSES",
    "GRID_KEYS",
    "MAP_DOSE_RATES",
    "MAP_PO2",
    "PAPER_DOSE_RATES",
    "TIME_GRID",
    "TS_PO2",
    "PointResult",
    "delta_o2",
    "index_of",
    "load_grid",
    "net_yield",
    "net_yield_map",
    "o2_from_po2",
    "run_grid",
    "run_point",
    "save_grid",
]

# --------------------------------------------------------------------------- #
#                                   Logging                                   #
# --------------------------------------------------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                                  Constants                                  #
# --------------------------------------------------------------------------- #
LET = 1.17  # keV/µm (paper table 3)
T_END = 500.0  # s
HENRY = 1.3e3  # µM/atm, [O2] = pO2 [atm] x HENRY
ON_STEPS = 500  # minimum number of solver steps over the pulse
MAX_STEP_OFF = 0.1  # s, final beam OFF chunk

# Paper table 5 (dose rates of the time-series and bar figures)
PAPER_DOSE_RATES: tp.Tuple[float, ...] = (
    0.1, 1.0, 10.0, 50.0, 1e2, 1e3, 1e4, 1e5, 1e6, 5e6, 1e7,
)
# Map grid: 0.1 to 1e7 Gy/s at two per decade, plus the table 5 values
MAP_DOSE_RATES: tp.Tuple[float, ...] = tuple(
    sorted({float(10.0 ** (k / 2)) for k in range(-2, 15)} | set(PAPER_DOSE_RATES))
)
# Initial pO2 [%] of the bar figures (paper figs 6 and 9)
BAR_PO2: tp.Tuple[float, ...] = (0.0, 0.1, 0.5, 2.0, 4.0, 6.0, 12.0, 16.0, 21.0)
# Initial pO2 [%] of the time-series figures: hypoxia, physioxia
TS_PO2: tp.Tuple[float, ...] = (1.0, 5.0)
MAP_PO2: tp.Tuple[float, ...] = tuple(
    sorted(
        set(BAR_PO2)
        | set(TS_PO2)
        | {0.25, 0.75, 1.5, 3.0, 5.0, 7.5, 10.0, 14.0, 18.0}
    )
)
# Doses of the net-yield fit; the last one is the paper's 30 Gy delivery.
# At 1e7 Gy/s, 10 Gy is a 1 µs pulse (MIN_BEAM_ON_TIME).
DOSES: tp.Tuple[float, ...] = (10.0, 20.0, 30.0)
# Time grid on which O2, H2O2 and pH are stored (first stored point: 1e-7 s)
TIME_GRID: NDArray[np.float64] = np.logspace(-7.0, np.log10(T_END), 300)

# Keys of the npz file written by save_grid
GRID_KEYS = (
    "dose_rates", "po2", "doses", "o2_initial",
    "h2o2_final", "o2_final", "ph_min",
    "o2_production", "o2_consumption", "h2o2_production", "h2o2_consumption",
    "time", "ts_dose_rates", "ts_po2", "ts_o2", "ts_h2o2", "ts_ph",
)

# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #


def o2_from_po2(po2: float) -> float:
    """Initial [O2] in µM for an initial pO2 in %."""
    return po2 / 100.0 * HENRY


@contextlib.contextmanager
def _quiet_radiopyo() -> tp.Iterator[None]:
    """Silence the INFO logs of radiopyo for one run, then restore."""
    level = logger.level
    logger.setLevel(logging.WARNING)
    try:
        yield
    finally:
        logger.setLevel(level)


@functools.lru_cache(maxsize=1)
def _template_cell() -> UnitCell:
    with _quiet_radiopyo():
        cell = UnitCell.from_toml(PURE_WATER_CONFIG_FILE)
    cell.env.update_flag_GValues(True)
    return cell


def _solve(dose_rate: float, po2: float, dose: float) -> ODEResult:
    """One paper run: single pulse of `dose` at `dose_rate` from t = 0,
    initial pO2 `po2` [%], integrated over (0, T_END)."""
    cell = copy.deepcopy(_template_cell())
    cell.env.species["O2"].set_initial_cc(o2_from_po2(po2) * 1e-6)  # mol/L
    with _quiet_radiopyo():
        cell.set_beam(SinglePulseBeam(dose_rate=dose_rate, max_dose=dose), LET=LET)
        chunks = cell.prepare_single_pulse(
            t_span=(0.0, T_END),
            max_step_size_on=dose / dose_rate / ON_STEPS,
            max_step_size_final=MAX_STEP_OFF,
        )
        # Consume the chunks without the tqdm bar of SimChunksIter.run
        return collections.deque(chunks, maxlen=1).pop()


def _on_grid(t: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.interp(TIME_GRID, t, y)


def run_point(
    dose_rate: float,
    po2: float,
    dose: float,
    integrals: bool = False,
) -> PointResult:
    """Run one grid point and analyse it.

    Args:
        dose_rate: float
            Dose rate of the single pulse [Gy/s].
        po2: float
            Initial pO2 [%].
        dose: float
            Delivered dose [Gy].
        integrals: bool, optional
            Also compute the gross production and consumption of O2 and
            H2O2 over the run (default: False; they rebuild every reaction
            rate at every stored point).

    Returns:
        PointResult
    """
    res = _solve(dose_rate, po2, dose)
    df = res.to_pandas()
    t = df.index.to_numpy(dtype=np.float64)
    o2 = df["O2"].to_numpy(dtype=np.float64)
    h2o2 = df["H2O2"].to_numpy(dtype=np.float64)
    ph = -np.log10(df["H3O_plus"].to_numpy(dtype=np.float64) * 1e-6)
    out = PointResult(
        dose_rate=dose_rate,
        po2=po2,
        dose=dose,
        o2_initial=float(o2[0]),
        o2_final=float(o2[-1]),
        h2o2_final=float(h2o2[-1]),
        ph_min=float(ph.min()),
        ts_o2=_on_grid(t, o2),
        ts_h2o2=_on_grid(t, h2o2),
        ts_ph=_on_grid(t, ph),
    )
    if integrals:
        with _quiet_radiopyo():
            out.o2_production = float(res.integrate_production_of("O2").sum())
            out.o2_consumption = -float(res.integrate_consumption_of("O2").sum())
            out.h2o2_production = float(res.integrate_production_of("H2O2").sum())
            out.h2o2_consumption = -float(res.integrate_consumption_of("H2O2").sum())
    return out


def net_yield(
    doses: tp.Sequence[float] | NDArray[np.float64],
    concentrations: tp.Sequence[float] | NDArray[np.float64],
) -> float:
    """Net yield [#/100 eV]: slope of a linear fit (with intercept) of the
    final concentrations [µM] versus the delivered doses [Gy]."""
    slope = np.polyfit(np.asarray(doses, float), np.asarray(concentrations, float), 1)[0]
    return float(kr_to_ge(slope * 1e-6))  # µM/Gy -> mol/L/Gy -> #/100 eV


def index_of(values: tp.Sequence[float] | NDArray[np.float64], value: float) -> int:
    """Index of `value` in `values` (relative match, raises if absent)."""
    arr = np.asarray(values, dtype=np.float64)
    hits = np.flatnonzero(np.isclose(arr, value, rtol=1e-9, atol=1e-12))
    if hits.size != 1:
        raise KeyError(f"{value} not found in {arr}")
    return int(hits[0])


def _run_job(job: tp.Tuple[float, float, float, bool]) -> PointResult:
    return run_point(*job)


def run_grid(
    dose_rates: tp.Sequence[float] = MAP_DOSE_RATES,
    po2: tp.Sequence[float] = MAP_PO2,
    doses: tp.Sequence[float] = DOSES,
    ts_dose_rates: tp.Sequence[float] = PAPER_DOSE_RATES,
    ts_po2: tp.Sequence[float] = TS_PO2,
    workers: int = 1,
    progress: tp.Optional[tp.Callable[[int, int], None]] = None,
) -> tp.Dict[str, NDArray[np.float64]]:
    """Run every (dose rate, pO2, dose) of the grid.

    The gross production / consumption integrals are computed at the last
    (largest) dose only. Time series are kept for the (ts_dose_rates x
    ts_po2) points of the grid, at the last dose.

    Args:
        dose_rates, po2, doses: Sequence[float]
            Grid axes [Gy/s], [%], [Gy].
        ts_dose_rates, ts_po2: Sequence[float]
            Points whose time series are kept (those in the grid only).
        workers: int, optional
            Number of worker processes (default: 1, serial). With more
            than one, this module must be importable by name (on sys.path).
        progress: Callable[[int, int], None], optional
            Called with (done, total) after every run.

    Returns:
        Dict[str, NDArray]: arrays keyed by GRID_KEYS. Axes are
        (dose rate, pO2, dose) for h2o2_final, o2_final and ph_min,
        (dose rate, pO2) for the integrals (at the last dose) and
        (ts dose rate, ts pO2, time) for the time series.
    """
    n_dr, n_po2, n_d = len(dose_rates), len(po2), len(doses)
    last = n_d - 1
    jobs = [
        (float(dr), float(p), float(d), k == last)
        for dr in dose_rates
        for p in po2
        for k, d in enumerate(doses)
    ]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = _collect(pool.map(_run_job, jobs, chunksize=4), len(jobs), progress)
    else:
        results = _collect(map(_run_job, jobs), len(jobs), progress)

    shape = (n_dr, n_po2, n_d)
    out: tp.Dict[str, NDArray[np.float64]] = {
        "dose_rates": np.asarray(dose_rates, dtype=np.float64),
        "po2": np.asarray(po2, dtype=np.float64),
        "doses": np.asarray(doses, dtype=np.float64),
        "o2_initial": np.asarray([o2_from_po2(p) for p in po2], dtype=np.float64),
        "h2o2_final": np.empty(shape),
        "o2_final": np.empty(shape),
        "ph_min": np.empty(shape),
        "time": TIME_GRID.copy(),
    }
    for key in ("o2_production", "o2_consumption", "h2o2_production", "h2o2_consumption"):
        out[key] = np.empty((n_dr, n_po2))

    ts_dr = [dr for dr in ts_dose_rates if np.isclose(dose_rates, dr).any()]
    ts_p = [p for p in ts_po2 if np.isclose(po2, p).any()]
    out["ts_dose_rates"] = np.asarray(ts_dr, dtype=np.float64)
    out["ts_po2"] = np.asarray(ts_p, dtype=np.float64)
    for key in ("ts_o2", "ts_h2o2", "ts_ph"):
        out[key] = np.empty((len(ts_dr), len(ts_p), TIME_GRID.size))

    it_results = iter(results)
    for i in range(n_dr):
        for j in range(n_po2):
            for k in range(n_d):
                r = next(it_results)
                out["h2o2_final"][i, j, k] = r.h2o2_final
                out["o2_final"][i, j, k] = r.o2_final
                out["ph_min"][i, j, k] = r.ph_min
                if k != last:
                    continue
                out["o2_production"][i, j] = r.o2_production
                out["o2_consumption"][i, j] = r.o2_consumption
                out["h2o2_production"][i, j] = r.h2o2_production
                out["h2o2_consumption"][i, j] = r.h2o2_consumption
                if r.dose_rate in ts_dr and r.po2 in ts_p:
                    a, b = ts_dr.index(r.dose_rate), ts_p.index(r.po2)
                    out["ts_o2"][a, b] = r.ts_o2
                    out["ts_h2o2"][a, b] = r.ts_h2o2
                    out["ts_ph"][a, b] = r.ts_ph
    return out


def _collect(
    results: tp.Iterable[PointResult],
    total: int,
    progress: tp.Optional[tp.Callable[[int, int], None]],
) -> tp.List[PointResult]:
    out = []
    for n, r in enumerate(results, start=1):
        out.append(r)
        if progress is not None:
            progress(n, total)
    return out


def delta_o2(grid: tp.Mapping[str, NDArray[np.float64]]) -> NDArray[np.float64]:
    """O2 depletion [µM] at T_END after the last dose, (dose rate, pO2)."""
    return grid["o2_initial"][None, :] - grid["o2_final"][:, :, -1]


def net_yield_map(grid: tp.Mapping[str, NDArray[np.float64]]) -> NDArray[np.float64]:
    """Net H2O2 yield [#/100 eV] over the doses of the grid, (dose rate, pO2)."""
    doses = grid["doses"]
    h2o2 = grid["h2o2_final"]
    return np.asarray(
        [
            [net_yield(doses, h2o2[i, j]) for j in range(h2o2.shape[1])]
            for i in range(h2o2.shape[0])
        ]
    )


def save_grid(grid: tp.Mapping[str, NDArray[np.float64]], path: str | Path) -> None:
    """Write the grid arrays (GRID_KEYS) to a compressed npz file."""
    missing = set(GRID_KEYS) - set(grid)
    if missing:
        raise KeyError(f"missing grid arrays: {sorted(missing)}")
    arrays: tp.Dict[str, tp.Any] = {k: np.asarray(grid[k]) for k in GRID_KEYS}
    np.savez_compressed(Path(path), **arrays)


def load_grid(path: str | Path) -> tp.Dict[str, NDArray[np.float64]]:
    """Read a grid written by save_grid."""
    with np.load(Path(path)) as data:
        return {k: data[k] for k in GRID_KEYS}


# --------------------------------------------------------------------------- #
#                               CLASS DEFINITIONS                             #
# --------------------------------------------------------------------------- #
@dataclasses.dataclass
class PointResult:
    """Analysis of one run (concentrations in µM).

    The gross production and consumption are positive amounts [µM]
    integrated over (0, T_END), radiolytic yield included in the
    production; NaN when not computed (run_point(integrals=False)).
    Time series are sampled on TIME_GRID.
    """

    dose_rate: float
    po2: float
    dose: float
    o2_initial: float
    o2_final: float
    h2o2_final: float
    ph_min: float
    ts_o2: NDArray[np.float64]
    ts_h2o2: NDArray[np.float64]
    ts_ph: NDArray[np.float64]
    o2_production: float = float("nan")
    o2_consumption: float = float("nan")
    h2o2_production: float = float("nan")
    h2o2_consumption: float = float("nan")

    @property
    def delta_o2(self) -> float:
        """O2 depletion [µM] over the run."""
        return self.o2_initial - self.o2_final
