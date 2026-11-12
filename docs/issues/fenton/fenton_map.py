#!/usr/bin/env python3
"""
Map of the ROOH yield versus dose per pulse and average dose rate
(Fenton runaway case, Configuration_Bordet_Emballement.toml).

Author:
    radiopyo contributors
"""

# --------------------------------------------------------------------------- #
#                               External imports                              #
# --------------------------------------------------------------------------- #
from __future__ import annotations

import hashlib
import json
import logging
import sys
import time
import typing as tp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

# radiopyo is not pip-installed everywhere: fall back on the repo checkout
_REPO = Path(__file__).resolve().parents[3]
if (_REPO / "radiopyo" / "__init__.py").exists() and str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

# --------------------------------------------------------------------------- #
#                               Internal imports                              #
# --------------------------------------------------------------------------- #
import radiopyo as rp
from radiopyo.parser.file_parser import resolve_includes
from radiopyo.parser.toml import TOMLFileParser
from radiopyo.parser.utils import combine_config
from radiopyo.physics.beam import PulsedBeam

# --------------------------------------------------------------------------- #
#                         Expose importable stuffs                            #
# --------------------------------------------------------------------------- #
__all__ = ["compute_map", "compute_point", "load_config", "CONFIG_FILE", "CACHE_DIR"]

# --------------------------------------------------------------------------- #
#                                   Logging                                   #
# --------------------------------------------------------------------------- #
logger = logging.getLogger("radiopyo")

# --------------------------------------------------------------------------- #
#                                  Constants                                  #
# --------------------------------------------------------------------------- #
HERE = Path(__file__).resolve().parent
CONFIG_FILE = HERE / "Configuration_Bordet_Emballement.toml"
CACHE_DIR = HERE / ".map_cache"

ON_TIME = 1.2e-6  # s
LET = 0.2116  # keV/um
MAX_DOSE = 10.0  # Gy
DPP_RANGE = (0.03, 15.0)  # Gy per pulse
DR_RANGE = (0.2, 1000.0)  # average dose rate, Gy/s
MAX_STEP_ON = 1e-4
MAX_STEP_OFF = 0.5

# One UnitCell per worker process and per Fe2+ value (building is slow)
_CELLS: tp.Dict[tp.Optional[float], tp.Any] = {}

# --------------------------------------------------------------------------- #
#                              FUNCTION DEFINITIONS                           #
# --------------------------------------------------------------------------- #


def load_config(fe: tp.Optional[float] = None) -> tp.Any:
    """ConfigDict of the Fenton case: the water radiolysis, THIS folder's
    reactions_biology.toml and the beams of Configuration_Bordet_Emballement.toml
    (both files unchanged; the Configuration file includes the packaged
    biology, which is not the runaway case). The folder's biology file also
    sets H2O = 55, which the water radiolysis already defines: it is dropped
    in memory. Fe_plus_plus is overridden if `fe` is not None (mol/L)."""
    own = TOMLFileParser(CONFIG_FILE).parse()
    water = resolve_includes({**own, "includes": {"toml": ["radiopyo/advanced_water_radiolysis.toml"]}})
    bio = TOMLFileParser(HERE / "reactions_biology.toml").parse()
    bio["concentrations"]["fixed"].pop("H2O", None)
    own["includes"] = {}
    config = combine_config(combine_config(water, bio), own)
    if fe is not None:
        config["concentrations"]["fixed"]["Fe_plus_plus"] = fe
    return config


def _get_cell(fe: tp.Optional[float]) -> tp.Any:
    """UnitCell of the Fenton case (see load_config)."""
    if fe not in _CELLS:
        _CELLS[fe] = rp.UnitCell.from_config_dict(load_config(fe))
    return _CELLS[fe]


def compute_point(
    dose_per_pulse: float,
    avg_dose_rate: float,
    period: float = 1200.0,
    fe: tp.Optional[float] = None,
    on_time: float = ON_TIME,
    max_dose: float = MAX_DOSE,
    let: float = LET,
) -> tp.Tuple[float, float]:
    """Run one point of the map.

    Args:
        dose_per_pulse: Gy
        avg_dose_rate: Gy/s (frequency = avg_dose_rate / dose_per_pulse)
        period: simulated time, s (default 1200)
        fe: Fe_plus_plus fixed concentration in mol/L (None: config value)

    Returns:
        (ROOH at max dose, ROOH at end time) in umol/L; (nan, nan) if the
        point is infeasible (pulse period shorter than the pulse, max dose
        not reached within `period`) or if the solver failed.
    """
    nan = float("nan")
    frequency = avg_dose_rate / dose_per_pulse
    pulse_period = 1.0 / frequency
    if pulse_period <= on_time:
        return nan, nan
    t0 = time.perf_counter()
    try:
        beam = PulsedBeam.from_dose_per_pulse(
            dose_per_pulse, pulse_period, on_time, max_dose=max_dose
        )
        if beam.time_at_max_dose >= period:
            return nan, nan
        uc = _get_cell(fe)
        uc.set_beam(beam, LET=let)
        res = uc.prepare_chunked_run(
            [1e-9, period],
            max_step_size_on=MAX_STEP_ON,
            max_step_size_off=MAX_STEP_OFF,
        ).run()
        df = res.to_pandas()
        rooh = df["ROOH"].to_numpy()
        t = df.index.to_numpy()
        # last sample at (or just after) the end of the last pulse
        i_max = max(int(np.searchsorted(t, beam.time_at_max_dose * (1 + 1e-9), side="right")) - 1, 0)
        out = float(rooh[i_max]), float(rooh[-1])
    except Exception:  # noqa: BLE001 - a failed point must not kill the map
        logger.exception(
            f"Point failed (dpp={dose_per_pulse:.3g} Gy, dr={avg_dose_rate:.3g} Gy/s)"
        )
        return nan, nan
    logger.info(
        f"point dpp={dose_per_pulse:.3g} Gy dr={avg_dose_rate:.3g} Gy/s "
        f"done in {time.perf_counter() - t0:.1f} s"
    )
    return out


def _point_task(args: tp.Tuple[float, float, float, tp.Optional[float]]) -> tp.Tuple[float, float]:
    return compute_point(*args)


def _cache_path(cache: tp.Optional[str | Path], params: tp.Dict[str, tp.Any]) -> tp.Optional[Path]:
    if cache is None:
        return None
    key = hashlib.sha1(json.dumps(params, sort_keys=True).encode()).hexdigest()[:12]
    return Path(cache) / f"fenton_map_{key}.npz"


def compute_map(
    n: int = 12,
    dpp_range: tp.Tuple[float, float] = DPP_RANGE,
    dr_range: tp.Tuple[float, float] = DR_RANGE,
    period: float = 1200.0,
    fe: tp.Optional[float] = None,
    cache: tp.Optional[str | Path] = CACHE_DIR,
    max_workers: tp.Optional[int] = None,
) -> tp.Dict[str, np.ndarray]:
    """Compute ROOH on a log-spaced n x n grid.

    Returns:
        dict with `dose_per_pulse` (n,), `avg_dose_rate` (n,),
        `rooh_max_dose` and `rooh_end` of shape (n_dr, n_dpp) in umol/L
        (row = dose rate, column = dose per pulse), NaN where infeasible
        or failed.
    Cached in an .npz keyed by the grid parameters (cache=None disables).
    """
    dpp = np.geomspace(*dpp_range, n)
    dr = np.geomspace(*dr_range, n)
    params = {
        "n": n, "dpp": list(dpp_range), "dr": list(dr_range), "period": period,
        "fe": fe, "on_time": ON_TIME, "let": LET, "max_dose": MAX_DOSE,
        "step_on": MAX_STEP_ON, "step_off": MAX_STEP_OFF,
    }  # fmt: skip
    path = _cache_path(cache, params)
    if path is not None and path.exists():
        logger.info(f"Loading cached map {path}")
        with np.load(path) as f:
            return {k: f[k] for k in f.files}

    tasks = [(float(p), float(r), period, fe) for r in dr for p in dpp]
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        results = list(pool.map(_point_task, tasks))
    logger.info(f"Map of {len(tasks)} points done in {time.perf_counter() - t0:.0f} s")

    arr = np.array(results, dtype=np.float64).reshape(n, n, 2)
    out = {
        "dose_per_pulse": dpp,
        "avg_dose_rate": dr,
        "rooh_max_dose": arr[:, :, 0],
        "rooh_end": arr[:, :, 1],
    }
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, **out)
        logger.info(f"Saved map to {path}")
    return out
