#!/usr/bin/env python3
"""
Golden test of the paper reproduction (Tonneau et al. 2025, PMB 70
235021): ``examples/tonneau_2025_pmb/paper_runs.py`` must still produce the
grid recorded in ``tests/golden/data/paper_tonneau_2025.npz``.

The full grid (19 dose rates x 19 pO2 x 3 doses) takes minutes, so only a
few points are re-run here: two net yields (three doses each), the gross
production / consumption integrals of four bar points, and two time
series. The npz is written by ``tests/golden/record_golden.py --only
paper`` from the same module; re-record it only for an approved
result-changing fix.

The values are what radiopyo produces with the GF(DR) fit it ships, not
the numbers printed in the paper (ADR 0009); the notebook
``examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb`` compares both.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 12/11/2026
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import importlib.util
import sys
import types
import typing as tp
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------- #
#                                  CONSTANTS                                  #
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parents[2]
PAPER_RUNS_FILE = ROOT / "examples" / "tonneau_2025_pmb" / "paper_runs.py"
PAPER_FILE = Path(__file__).parent / "data" / "paper_tonneau_2025.npz"
# Same grid with the GF(DR) fit printed in the paper (record_golden.py
# --only paper-gf), for examples/tonneau_2025_pmb/gf_fit_comparison.ipynb
PAPER_GF_FILE = Path(__file__).parent / "data" / "paper_tonneau_2025_paper_gf.npz"
# Same grid with the inputs of the original runs quoted in the paper's text
# (paper_runs.SCENARIOS["original"]; record_golden.py --only paper-original)
PAPER_ORIGINAL_FILE = Path(__file__).parent / "data" / "paper_tonneau_2025_original.npz"

# Same code, same inputs: only the summation order of a refactor may move
# the solver path, so 1e-6 is far below any physical change.
RTOL = 1e-6
ATOL = 1e-9  # µM, for the near-zero integrals (H2O2 consumption at UHDR)
# Time series: np.interp of the solver points onto TIME_GRID. The LSODA step
# sequence moves with the floating-point rounding of the machine (BLAS
# kernels differ between CI runners), and a linear interpolation between
# other solver points shifts the sampled values by up to ~1e-5 relative
# (seen on GitHub runners: 7e-6 on H2O2 at 0.1 Gy/s), while the final
# values still agree to 1e-6.
TS_RTOL = 1e-4

# (dose rate [Gy/s], pO2 [%]) re-run at every dose: net yield, integrals and
# time series at the last dose (both are time-series points of the grid)
FULL_POINTS = ((0.1, 1.0), (1e7, 5.0))
# (dose rate, pO2) re-run at the last dose only: bar integrals
BAR_POINTS = ((1e3, 21.0), (1e7, 0.0))


# --------------------------------------------------------------------------- #
#                             FUNCTION DEFINITIONS                            #
# --------------------------------------------------------------------------- #
def load_paper_runs() -> types.ModuleType:
    """Load paper_runs.py by path, registered as ``paper_runs`` and with its
    directory on sys.path (so that run_grid(workers > 1) can pickle it)."""
    name = "paper_runs"
    if name in sys.modules:
        return sys.modules[name]
    folder = str(PAPER_RUNS_FILE.parent)
    if folder not in sys.path:
        sys.path.insert(0, folder)
    spec = importlib.util.spec_from_file_location(name, PAPER_RUNS_FILE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # required by dataclasses
    spec.loader.exec_module(module)
    return module


# ------------------------------- Fixtures ---------------------------------- #
@pytest.fixture(scope="module")
def pr() -> types.ModuleType:
    return load_paper_runs()


@pytest.fixture(scope="module")
def golden() -> tp.Dict[str, np.ndarray]:
    if not PAPER_FILE.exists():
        pytest.fail(f"{PAPER_FILE} is missing: run record_golden.py --only paper")
    return load_paper_runs().load_grid(PAPER_FILE)


@pytest.fixture(scope="module")
def full_runs(pr: types.ModuleType) -> tp.Dict[tp.Tuple[float, float], tp.List[tp.Any]]:
    last = len(pr.DOSES) - 1
    return {
        (dr, po2): [
            pr.run_point(dr, po2, dose, integrals=(k == last))
            for k, dose in enumerate(pr.DOSES)
        ]
        for dr, po2 in FULL_POINTS
    }


# --------------------------------- Tests ----------------------------------- #
def test_grid_layout(pr: types.ModuleType, golden: tp.Dict[str, np.ndarray]) -> None:
    """The npz matches the grid of paper_runs (axes, shapes, time grid)."""
    assert set(golden) == set(pr.GRID_KEYS)
    np.testing.assert_array_equal(golden["dose_rates"], pr.MAP_DOSE_RATES)
    np.testing.assert_array_equal(golden["po2"], pr.MAP_PO2)
    np.testing.assert_array_equal(golden["doses"], pr.DOSES)
    np.testing.assert_array_equal(golden["ts_dose_rates"], pr.PAPER_DOSE_RATES)
    np.testing.assert_array_equal(golden["ts_po2"], pr.TS_PO2)
    np.testing.assert_allclose(golden["time"], pr.TIME_GRID, rtol=1e-15)
    shape = (len(pr.MAP_DOSE_RATES), len(pr.MAP_PO2), len(pr.DOSES))
    for key in ("h2o2_final", "o2_final", "ph_min"):
        assert golden[key].shape == shape, key
        assert np.all(np.isfinite(golden[key])), key
    for key in ("o2_production", "o2_consumption", "h2o2_production", "h2o2_consumption"):
        assert golden[key].shape == shape[:2], key
        assert np.all(golden[key] >= 0.0), key
    ts_shape = (len(pr.PAPER_DOSE_RATES), len(pr.TS_PO2), pr.TIME_GRID.size)
    for key in ("ts_o2", "ts_h2o2", "ts_ph"):
        assert golden[key].shape == ts_shape, key


@pytest.mark.parametrize("point", FULL_POINTS, ids=lambda p: f"{p[0]:g}Gy/s-{p[1]:g}%")
def test_net_yield(
    pr: types.ModuleType,
    golden: tp.Dict[str, np.ndarray],
    full_runs: tp.Dict[tp.Tuple[float, float], tp.List[tp.Any]],
    point: tp.Tuple[float, float],
) -> None:
    i = pr.index_of(golden["dose_rates"], point[0])
    j = pr.index_of(golden["po2"], point[1])
    runs = full_runs[point]
    h2o2 = [r.h2o2_final for r in runs]
    np.testing.assert_allclose(h2o2, golden["h2o2_final"][i, j], rtol=RTOL)
    np.testing.assert_allclose(
        [r.o2_final for r in runs], golden["o2_final"][i, j], rtol=RTOL
    )
    np.testing.assert_allclose(
        [r.ph_min for r in runs], golden["ph_min"][i, j], rtol=RTOL
    )
    expected = pr.net_yield(golden["doses"], golden["h2o2_final"][i, j])
    assert pr.net_yield(pr.DOSES, h2o2) == pytest.approx(expected, rel=RTOL)
    assert pr.net_yield_map(golden)[i, j] == pytest.approx(expected, rel=1e-12)


INTEGRALS = ("o2_production", "o2_consumption", "h2o2_production", "h2o2_consumption")


@pytest.mark.parametrize(
    "point", FULL_POINTS + BAR_POINTS, ids=lambda p: f"{p[0]:g}Gy/s-{p[1]:g}%"
)
def test_bar_integrals(
    pr: types.ModuleType,
    golden: tp.Dict[str, np.ndarray],
    full_runs: tp.Dict[tp.Tuple[float, float], tp.List[tp.Any]],
    point: tp.Tuple[float, float],
) -> None:
    i = pr.index_of(golden["dose_rates"], point[0])
    j = pr.index_of(golden["po2"], point[1])
    if point in full_runs:
        r = full_runs[point][-1]
    else:
        r = pr.run_point(point[0], point[1], pr.DOSES[-1], integrals=True)
        assert r.h2o2_final == pytest.approx(golden["h2o2_final"][i, j, -1], rel=RTOL)
    got = [getattr(r, key) for key in INTEGRALS]
    np.testing.assert_allclose(
        got, [golden[key][i, j] for key in INTEGRALS], rtol=RTOL, atol=ATOL
    )


@pytest.mark.parametrize("point", FULL_POINTS, ids=lambda p: f"{p[0]:g}Gy/s-{p[1]:g}%")
def test_time_series(
    pr: types.ModuleType,
    golden: tp.Dict[str, np.ndarray],
    full_runs: tp.Dict[tp.Tuple[float, float], tp.List[tp.Any]],
    point: tp.Tuple[float, float],
) -> None:
    a = pr.index_of(golden["ts_dose_rates"], point[0])
    b = pr.index_of(golden["ts_po2"], point[1])
    r = full_runs[point][-1]
    for key in ("ts_o2", "ts_h2o2", "ts_ph"):
        np.testing.assert_allclose(
            getattr(r, key), golden[key][a, b], rtol=TS_RTOL, atol=ATOL, err_msg=key
        )


# ------------------------- GF fit of the paper ----------------------------- #
# (dose rate, pO2) re-run at the last dose: the two fits differ most at UHDR
GF_POINTS = ((0.1, 1.0), (1e7, 5.0))


@pytest.fixture(scope="module")
def golden_paper_gf(pr: types.ModuleType) -> tp.Dict[str, np.ndarray]:
    if not PAPER_GF_FILE.exists():
        pytest.fail(f"{PAPER_GF_FILE} is missing: run record_golden.py --only paper-gf")
    return pr.load_grid(PAPER_GF_FILE)


def test_paper_gf_layout(
    golden: tp.Dict[str, np.ndarray], golden_paper_gf: tp.Dict[str, np.ndarray]
) -> None:
    """Same axes as the grid with the shipped fit."""
    for key in ("dose_rates", "po2", "doses", "ts_dose_rates", "ts_po2", "time"):
        np.testing.assert_array_equal(golden_paper_gf[key], golden[key], err_msg=key)
    for key in golden:
        assert golden_paper_gf[key].shape == golden[key].shape, key


@pytest.mark.parametrize("point", GF_POINTS, ids=lambda p: f"{p[0]:g}Gy/s-{p[1]:g}%")
def test_paper_gf_point(
    pr: types.ModuleType,
    golden_paper_gf: tp.Dict[str, np.ndarray],
    point: tp.Tuple[float, float],
) -> None:
    i = pr.index_of(golden_paper_gf["dose_rates"], point[0])
    j = pr.index_of(golden_paper_gf["po2"], point[1])
    r = pr.run_point(point[0], point[1], pr.DOSES[-1], integrals=True, gf_fit="paper")
    got = [r.h2o2_final, r.o2_final, r.ph_min] + [getattr(r, k) for k in INTEGRALS]
    expected = [golden_paper_gf[k][i, j, -1] for k in ("h2o2_final", "o2_final", "ph_min")]
    expected += [golden_paper_gf[k][i, j] for k in INTEGRALS]
    np.testing.assert_allclose(got, expected, rtol=RTOL, atol=ATOL)


def test_gf_fit_is_restored(pr: types.ModuleType) -> None:
    """A run with the paper fit leaves the shipped GF(DR) in place."""
    from radiopyo.utils.g_values import GValuesBoscolo

    before = GValuesBoscolo._scale_factor(1e7)
    pr.run_point(1e7, 0.0, pr.DOSES[-1], gf_fit="paper")
    assert GValuesBoscolo._scale_factor(1e7) == before
    assert pr.gf_factor(1e7, "paper") == pytest.approx(0.982 * 1e7**-0.09)

# ------------------------ Inputs of the original runs ---------------------- #
# Final H2O2 [µM] at 30 Gy of the original runs whose values are quoted in the
# paper's text (author's sim_GVal_Dr_pH3.h5, 2025-02-28, written by the code
# of that time with CONFIG_PAPER_01.toml, LET 0.56 keV/µm, shipped GF fit),
# rounded to 4 digits. The current code reproduces them to about 3e-4.
ORIGINAL_RUNS_H2O2 = {(1e2, 1.0): 2.185, (1e7, 5.0): 1.103}


@pytest.fixture(scope="module")
def golden_original(pr: types.ModuleType) -> tp.Dict[str, np.ndarray]:
    if not PAPER_ORIGINAL_FILE.exists():
        pytest.fail(
            f"{PAPER_ORIGINAL_FILE} is missing: run record_golden.py --only paper-original"
        )
    return pr.load_grid(PAPER_ORIGINAL_FILE)


def test_original_layout(
    golden: tp.Dict[str, np.ndarray], golden_original: tp.Dict[str, np.ndarray]
) -> None:
    """Same axes as the grid of the paper methods."""
    for key in ("dose_rates", "po2", "doses", "ts_dose_rates", "ts_po2", "time"):
        np.testing.assert_array_equal(golden_original[key], golden[key], err_msg=key)
    for key in golden:
        assert golden_original[key].shape == golden[key].shape, key


@pytest.mark.parametrize(
    "point", tuple(ORIGINAL_RUNS_H2O2), ids=lambda p: f"{p[0]:g}Gy/s-{p[1]:g}%"
)
def test_original_point(
    pr: types.ModuleType,
    golden_original: tp.Dict[str, np.ndarray],
    point: tp.Tuple[float, float],
) -> None:
    """The original scenario matches its golden grid and the original runs."""
    i = pr.index_of(golden_original["dose_rates"], point[0])
    j = pr.index_of(golden_original["po2"], point[1])
    r = pr.run_point(point[0], point[1], pr.DOSES[-1], integrals=True, scenario="original")
    got = [r.h2o2_final, r.o2_final, r.ph_min] + [getattr(r, k) for k in INTEGRALS]
    expected = [golden_original[k][i, j, -1] for k in ("h2o2_final", "o2_final", "ph_min")]
    expected += [golden_original[k][i, j] for k in INTEGRALS]
    np.testing.assert_allclose(got, expected, rtol=RTOL, atol=ATOL)
    assert r.h2o2_final == pytest.approx(ORIGINAL_RUNS_H2O2[point], rel=1e-3)


def test_unknown_scenario(pr: types.ModuleType) -> None:
    with pytest.raises(KeyError, match="unknown scenario"):
        pr.run_grid(dose_rates=(1.0,), po2=(1.0,), doses=(1.0,), scenario="nope")
