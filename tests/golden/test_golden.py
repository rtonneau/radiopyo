#!/usr/bin/env python3
"""
Golden-output (characterization) tests.

The reference data in ``tests/golden/data`` were produced by
``tests/golden/record_golden.py`` running the package code of commit
2261b92, *before* any refactoring. These tests pin the observable
behaviour of the package (ODE right-hand side, full solver runs, per
reaction rates, beam time structure, HDF5 persistence) so that later
behaviour-preserving refactors can be checked against it.

Design rules:

* Only stable public entry points are used (``UnitCell.from_toml``,
  ``set_beam``, ``run``, ``prepare_chunked_run(...).run()``,
  ``ODEResult.to_pandas``, ``ODEResult.reconstruct_reactions``,
  ``beam.make_matrix``, ``H5Manager`` / ``SimCollection``). The only
  exception is ``derive()``, which is isolated in ``_derive_adapter``.
* Species indices depend on ``set`` iteration order (and therefore on
  ``PYTHONHASHSEED``). Every golden is keyed by species label, and
  reactions by ``str(reaction)`` (``as_label()`` of a k-reaction is
  ``KR_<index>``, i.e. index based, so it is NOT used for k-reactions).
* ``from __future__ import annotations`` is deliberately NOT used in this
  module: ``SimIdentifier`` subclasses need real (non string) type
  annotations because ``SimIdentifier.to_string`` calls ``field.type(...)``.

This module also hosts the ``compute_*`` functions that the recorder
re-uses, so that the recorded data and the tested data are produced by
exactly the same code path.

Author:
    Romain Tonneau (romain.tonneau@unamur.be) - 09/11/2026
"""
# ----------------------------- External imports ---------------------------- #
import dataclasses as dc
import json
import logging
import shutil
import typing as tp
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.physics.beam import (
    ConstantBeam,
    MacroPulsedBeam,
    ParticleBeam,
    PulsedBeam,
    SinglePulseBeam,
)
from radiopyo.simulation.ode_functions_solver import derive
from radiopyo.simulation.ode_result import ODEResult
from radiopyo.simulation.unit_cell import UnitCell
from radiopyo.simulation.utils import H5Manager, SimCollection, SimIdentifier

# ---------------------------------- Logging -------------------------------- #
logging.getLogger("radiopyo").setLevel(logging.ERROR)

# --------------------------------------------------------------------------- #
#                                  CONSTANTS                                  #
# --------------------------------------------------------------------------- #
DATA_DIR = Path(__file__).parent / "data"
DERIVE_FILE = DATA_DIR / "derive.json"
RUNS_FILE = DATA_DIR / "runs.npz"
BEAMS_FILE = DATA_DIR / "beams.npz"
HDF_FILE = DATA_DIR / "golden_results.h5"
HDF_EXPECTED_FILE = DATA_DIR / "hdf_expected.npz"

# ------------------------------- Tolerances -------------------------------- #
# derive(), beam matrices: pure floating point algebra on identical inputs.
# The only differences a behaviour-preserving refactor can introduce are
# rounding differences from a different summation order (a few ulps,
# ~1e-15 relative), so 1e-12 leaves ~3 orders of magnitude of headroom
# while still catching any real change (a wrong coefficient is >> 1e-6).
RTOL_ALGEBRA = 1e-12
# Entries of dydt that are the result of strong cancellation can have a
# large *relative* rounding error while being negligible with respect to
# the rest of the vector. They are therefore also allowed an absolute error
# of RTOL_ALGEBRA * max|golden vector|.
ATOL_ALGEBRA_FRAC_OF_MAX = 1e-12

# Full solver runs (LSODA, solver rtol=1e-10 and atol=1e-12 uM by default).
# The adaptive solver amplifies last-bit rounding differences (different
# PYTHONHASHSEED = different species order = different summation order)
# into different step sequences: even the *number* of solver points of an
# unchanged code varies by several percent with the hash seed. What stays
# comparable are values taken at fixed times where the solver output is
# its own high order dense output (t_eval), or at chunk boundaries.
# Measured on the unmodified package over many hash seeds:
#   * relative differences are <~1e-9 for every quantity far above the
#     solver absolute tolerance;
#   * radicals living at 1e-8..1e-13 uM, i.e. at the level of the solver
#     absolute tolerance (1e-12 uM), differ by up to a few 1e-12 uM
#     absolute (up to ~1e-3 relative).
# Since ADR 0003 runs use the analytic Jacobian by default (use_jac=True)
# while the data were recorded with finite differences; the benchmark
# (final states, analytic vs FD, LSODA and BDF, bundled configs) showed
# differences <= 2e-9 uM absolute and <= 8e-9 relative, so the tolerances
# below were left unchanged.
# RTOL_RUN is thus ~3 orders of magnitude above the relative spread...
RTOL_RUN = 1e-6
# ... and the absolute floor for species (uM) is ~100x the solver atol.
ATOL_SPECIES_RUN = 1e-10
# Reaction rates built from those noisy low-level radicals (e.g. the
# HO3 -> OH_r + O2 rate, first order in a ~1e-13 uM species) get an absolute
# floor expressed as a fraction of the largest |rate| of ALL reactions over
# the compared points (measured spread over 11 hash seeds: <= 3.5e-10, so
# ~300x headroom).
ATOL_RATES_RUN_FRAC_OF_MAX = 1e-7
# Number of solver points of the plain run. Large variations mean that the
# solver is driven differently (measured spread with hash seed: ~+-4%).
N_POINTS_RELATIVE_SPREAD = 0.15

# Data read back from HDF5 are float64 stored losslessly: must be (almost)
# bit-identical to the recorded expectation.
RTOL_STORED = 1e-12

# ------------------------------- Scenarios --------------------------------- #
CONFIGS: tp.Dict[str, Path] = {
    "labarbe": rp.LABARBE_CONFIG_FILE,
    "pure_water": rp.PURE_WATER_CONFIG_FILE,
}
LET_DERIVE = 0.201  # keV/µm
T_BEAM_ON = 5e-7  # s, inside the 1e-6 s pulse (1e7 Gy/s, 10 Gy)
T_BEAM_OFF = 1e-3  # s, after max_dose has been delivered
RUN_T_SPAN = (1e-7, 10.0)
N_CURVE = 100  # t_eval points of the dense-output comparison of the Labarbe run

# Pulsed beam used by the chunked run (3 pulses of 2 Gy, 2 µs, 100 Hz)
CHUNK_T_SPAN = (0.0, 0.05)
CHUNK_LET = 2.5  # keV/µm

HDF_KEYS = ["golden"]
N_HDF_POINTS = 40


def _beam_cases() -> tp.Dict[str, tp.Tuple[ParticleBeam, tp.List[tp.Tuple[float, float]]]]:
    """Beams (finite max_dose) and the (start, stop) windows to evaluate."""
    return {
        "constant": (
            ConstantBeam(dose_rate=1e3, max_dose=5.0, start_time=0.2),
            [
                (0.0, 0.1),
                (0.0, 1.0),
                (0.201, 0.204),
                (0.201, 1.0),
                (0.3, 1.0),
                (0.0, 0.204),
                (0.25, 0.3),
            ],
        ),
        "singlepulse": (
            SinglePulseBeam(dose_rate=1e6, max_dose=2.0, start_time=1e-3),
            [
                (0.0, 0.1),
                (0.0, 1.5e-3),
                (1.0005e-3, 0.1),
                (0.5, 1.0),
                (0.0, 1.0005e-3),
            ],
        ),
        "pulsed": (
            PulsedBeam(
                dose_rate_peak=1e6,
                period=1e-2,
                on_time=2e-6,
                max_dose=7.0,
                start_time=1e-3,
            ),
            [
                (0.0, 0.1),
                (0.0, 0.02),
                (1e-3, 0.05),
                (0.0105, 0.1),
                (0.011 + 1e-6, 0.1),
                (0.5, 1.0),
            ],
        ),
        "pulsed_whole_pulses": (
            PulsedBeam(
                dose_rate_peak=1e6,
                period=1e-2,
                on_time=2e-6,
                max_dose=6.0,
                start_time=1e-3,
            ),
            [(0.0, 0.1), (0.0, 0.02), (1e-3, 0.05)],
        ),
        "macropulsed": (
            MacroPulsedBeam(
                dose_rate_peak=1e6,
                period=0.1,
                on_time=1e-3,
                micro_period=1e-4,
                micro_on_time=2e-6,
                max_dose=50.0,
                micro_start_time=0.0,
                start_time=0.01,
            ),
            [(0.0, 1.0), (0.01, 0.5), (0.05, 1.0), (0.0, 0.05)],
        ),
    }


# --------------------------------------------------------------------------- #
#                                  ADAPTERS                                   #
# --------------------------------------------------------------------------- #
# Everything that touches a part of the package that is expected to change
# signature/layout during the refactor lives in this section. When a ticket
# changes these internals, ONLY these functions need updating.


def _derive_adapter(
    cell: UnitCell,
    t: float,
    y: np.ndarray,
    use_beam: bool,
    force_positive: bool,
) -> np.ndarray:
    """
    Single place calling ``derive()``. Today its signature is
    ``derive(t, y, beam, mat, force_positive)`` with ``mat`` the
    ``SimMatrices`` of ``cell.mat``. Later tickets restructure these
    matrices: update only this function.

    Args:
        y: state vector in dynamical-species index order, [µmol/L]

    Returns:
        dydt in dynamical-species index order, [µmol/L/s]
    """
    return derive(
        t,
        y,
        cell.beam if use_beam else None,
        cell.mat,
        force_positive,
    )


def _dyn_species_labels(cell: UnitCell) -> tp.List[str]:
    """Labels of the dynamical species, ordered by their state index."""
    sp = sorted(cell.env.species.iter_dyn_species(), key=lambda s: s.index)
    assert [s.index for s in sp] == list(range(len(sp)))
    return [s.label for s in sp]


def _reaction_keys(res: ODEResult) -> tp.List[str]:
    """
    Index independent keys (``str(reaction)``) of the columns of
    ``res.reconstruct_reactions()``: k-reactions (acid/base, equilibrium
    and plain k-reactions), then radiolytic, then enzymatic.
    """
    reactions = res.env.reactions
    keys = [str(r) for r in reactions.iter_all_kreactions()]
    keys += [str(r) for r in reactions.radiolytic]
    keys += [str(r) for r in reactions.enzymatic]
    assert len(set(keys)) == len(keys), "str(reaction) is not unique"
    return keys


def _g_values(cell: UnitCell) -> tp.Dict[str, float]:
    """Radiolytic yields currently in the env, keyed by 'G_<species>'."""
    return {r.as_label(): float(r.ge()) for r in cell.env.reactions.radiolytic}


# --------------------------------------------------------------------------- #
#                         COMPUTE FUNCTIONS (shared)                          #
# --------------------------------------------------------------------------- #


def build_cell(config: str) -> UnitCell:
    """UnitCell of a bundled config, with a beam and LET defined."""
    cell = UnitCell.from_toml(CONFIGS[config])
    if len(cell.beam_list) == 0:  # pure water has no beam in its file
        cell.set_beam(ConstantBeam(dose_rate=1e7, max_dose=10.0), LET=LET_DERIVE)
    return cell


def _synthetic_value(label: str) -> float:
    """Deterministic (hash seed independent) concentration in [1e-4, 10] µM."""
    frac = (zlib.crc32(label.encode()) % 10007) / 10007.0
    return float(10 ** (-4.0 + 5.0 * frac))


def _derive_states(cell: UnitCell) -> tp.Dict[str, tp.Dict[str, float]]:
    """Fixed states, as {state name: {species label: µmol/L}}."""
    labels = _dyn_species_labels(cell)
    initial = dict(zip(labels, cell.env.initial_values()))
    synthetic = {lb: _synthetic_value(lb) + initial[lb] for lb in labels}
    # Exactly-zero radicals exercise the "negative product" guard of the
    # k-reactions and the force_positive branch.
    zeros = {
        lb: (0.0 if lb in ("e_aq", "H_r", "OH_r") else val)
        for lb, val in synthetic.items()
    }
    return {"initial": initial, "synthetic": synthetic, "zeros": zeros}


# (name, time, use_beam, force_positive)
DERIVE_CASES: tp.List[tp.Tuple[str, float, bool, bool]] = [
    ("beam_on", T_BEAM_ON, True, False),
    ("beam_off", T_BEAM_OFF, True, False),
    ("no_beam", T_BEAM_ON, False, False),
    ("beam_on_force_positive", T_BEAM_ON, True, True),
]


def compute_derive() -> tp.Dict[str, tp.Any]:
    """dydt (keyed by species label) at several fixed states, both configs."""
    out: tp.Dict[str, tp.Any] = {}
    for config in CONFIGS:
        cell = build_cell(config)
        # The scenario must not be vacuous: beam ON at T_BEAM_ON, OFF later
        assert cell.beam.at(T_BEAM_ON).dose_rate() > 0.0
        assert cell.beam.at(T_BEAM_OFF).dose_rate() == 0.0
        labels = _dyn_species_labels(cell)
        entry: tp.Dict[str, tp.Any] = {"g_values": _g_values(cell), "dydt": {}}
        for state_name, state in _derive_states(cell).items():
            y = np.array([state[lb] for lb in labels], dtype=np.float64)
            for case, t, use_beam, force_positive in DERIVE_CASES:
                dydt = _derive_adapter(cell, t, y.copy(), use_beam, force_positive)
                entry["dydt"][f"{state_name}__{case}"] = dict(
                    zip(labels, (float(v) for v in dydt))
                )
        out[config] = entry
    return out


def _result_to_arrays(
    prefix: str,
    res: ODEResult,
    indices: tp.Optional[tp.Sequence[int]] = None,
    with_rates: bool = True,
) -> tp.Dict[str, np.ndarray]:
    """
    Flatten (some points of) a run into label keyed arrays:
        <prefix>__time, <prefix>__n_points,
        <prefix>__species_labels, <prefix>__species,
        <prefix>__rate_labels, <prefix>__rate
    Species are in µM, rates in the units of ``reconstruct_reactions()``.

    Args:
        indices: positions of the solver points to keep (default: all)
    """
    species = res.to_pandas()
    species = species[sorted(species.columns)]
    rows = slice(None) if indices is None else list(indices)

    arrays: tp.Dict[str, np.ndarray] = {}
    arrays[f"{prefix}__time"] = species.index.to_numpy()[rows]
    arrays[f"{prefix}__n_points"] = np.array(res.n_points)
    arrays[f"{prefix}__species_labels"] = np.array(list(species.columns))
    arrays[f"{prefix}__species"] = species.to_numpy()[rows].T

    if with_rates:
        rates = res.reconstruct_reactions()
        rates.columns = pd.Index(_reaction_keys(res))
        rates = rates[sorted(rates.columns)]
        arrays[f"{prefix}__rate_labels"] = np.array(list(rates.columns))
        arrays[f"{prefix}__rate"] = rates.to_numpy()[rows].T
    return arrays


def compute_labarbe_run() -> tp.Dict[str, np.ndarray]:
    """
    ``run(t_span=[1e-7, 10])`` on LABARBE_CONFIG_FILE (default beam).

    Two runs are stored:
    * 'labarbe_run': the plain run, only first and last solver points are
      kept (solver steps in between are not reproducible across hash
      seeds), plus its number of points;
    * 'labarbe_curve': same run with ``t_eval`` (N_CURVE log spaced times)
      giving a reproducible description of the whole trajectory.
    Note: ``run()`` currently skips the initial radiolytic burst (it passes
    ``force_positive`` in the ``beam_ON`` slot of ``_raw_run``); this is the
    behaviour being recorded.
    """
    cell = build_cell("labarbe")
    res = cell.run(t_span=RUN_T_SPAN)
    out = _result_to_arrays("labarbe_run", res, indices=[0, -1])

    grid = np.geomspace(RUN_T_SPAN[0], RUN_T_SPAN[1], N_CURVE)
    cell = build_cell("labarbe")
    res = cell.run(t_span=RUN_T_SPAN, t_eval=grid)
    out.update(_result_to_arrays("labarbe_curve", res))
    return out


def compute_chunked_run() -> tp.Dict[str, np.ndarray]:
    """
    Chunked pulsed run on PURE_WATER_CONFIG_FILE, LET set. Only the first
    and last solver point of every chunk are kept: they sit at the exact
    chunk boundary times, whatever the solver steps in between.
    """
    cell = UnitCell.from_toml(CONFIGS["pure_water"])
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=2.0, period=1e-2, on_time=2e-6, n_pulse=3.0
    )
    cell.set_beam(beam, LET=CHUNK_LET)
    res = cell.prepare_chunked_run(
        t_span=CHUNK_T_SPAN,
        max_step_size_on=1e-7,
        max_step_size_off=1e-3,
    ).run()
    ends = np.cumsum(res.n_points) - 1
    starts = np.concatenate([[0], ends[:-1] + 1])
    edges = sorted(set(starts.tolist()) | set(ends.tolist()))
    arrays = _result_to_arrays("chunked_run", res, indices=edges)

    gv = _g_values(cell)
    arrays["chunked_run__g_labels"] = np.array(sorted(gv))
    arrays["chunked_run__g_values"] = np.array([gv[k] for k in sorted(gv)])
    return arrays


def compute_beams() -> tp.Dict[str, np.ndarray]:
    """make_matrix / make_time_intervals of every beam type and window."""
    out: tp.Dict[str, np.ndarray] = {}
    for name, (beam, windows) in _beam_cases().items():
        out[f"{name}__scalars"] = np.array(
            [beam.time_at_max_dose, beam.dose_rate, beam.dose_rate_peak]
        )
        for idx, (start, stop) in enumerate(windows):
            out[f"{name}__{idx}__matrix"] = beam.make_matrix(stop=stop, start=start)
            out[f"{name}__{idx}__intervals"] = beam.make_time_intervals(
                stop=stop, start=start
            )
    return out


# ------------------------------ HDF5 scenario ------------------------------ #
@dc.dataclass(kw_only=True)
class GoldenSimID(SimIdentifier):
    """Identifier of the stored results (real annotations on purpose)."""

    config: str = "labarbe"
    LET: float = LET_DERIVE
    t_end: float = 10.0
    with_rates: bool = False


def _hdf_runs() -> tp.Dict[str, tp.Tuple[GoldenSimID, UnitCell, float]]:
    """{expectation prefix: (identifier, cell, t_end)}"""
    out = {}
    for config, t_end, with_rates in (
        ("labarbe", 10.0, True),
        ("pure_water", 1e-3, False),
    ):
        cell = build_cell(config)
        ident = GoldenSimID(
            config=config, LET=cell.LET, t_end=t_end, with_rates=with_rates
        )
        ident.set_keys(list(HDF_KEYS))
        out[f"hdf_{config}"] = (ident, cell, t_end)
    return out


def _read_collection(path: Path) -> tp.Dict[str, ODEResult]:
    """Read every stored result through SimCollection, by identifier string."""
    collection = SimCollection(path, GoldenSimID(), keys=list(HDF_KEYS))
    return {sim_id.to_string(): res for _label, sim_id, res in collection}


def _rebuilt_rates_arrays(prefix: str, res: ODEResult) -> tp.Dict[str, np.ndarray]:
    """Rates recomputed (recompute=True) on a result read back from HDF5."""
    rates = res.reconstruct_reactions(recompute=True)
    rates.columns = pd.Index(_reaction_keys(res))
    rates = rates[sorted(rates.columns)]
    return {
        f"{prefix}__rebuilt_rate_labels": np.array(list(rates.columns)),
        f"{prefix}__rebuilt_rate": rates.to_numpy().T,
    }


def compute_hdf(path: Path) -> tp.Dict[str, np.ndarray]:
    """
    Write the HDF5 file at ``path`` (and its sibling json store folder)
    through SimCollection/H5Manager, read it back, and return the
    expectation arrays:
    * '<prefix>__species'/'__rate': the original in-memory results;
    * '<prefix>__rebuilt_rate': rates recomputed on the result read back.
      The committed json store was written by the old code, which doubled
      the k_value of the enzymatic (catalase) reaction (1.324e8 -> 2.648e8):
      reading it back still gives these doubled rates (pinned on purpose).
      Since MichaelisMenten.to_dict was fixed, a freshly written store
      round-trips the k_value and gives rates equal to the original ones.
    """
    path = Path(path)
    collection = SimCollection(path, GoldenSimID(), keys=list(HDF_KEYS), new=True)
    expected: tp.Dict[str, np.ndarray] = {}
    for prefix, (ident, cell, t_end) in _hdf_runs().items():
        t_eval = np.geomspace(RUN_T_SPAN[0], t_end, N_HDF_POINTS)
        res = cell.run(t_span=(RUN_T_SPAN[0], t_end), t_eval=t_eval)
        collection.hdf.push_ODEResult(
            ident, res, keys=list(HDF_KEYS), production_rate=ident.with_rates
        )
        expected.update(_result_to_arrays(prefix, res, with_rates=ident.with_rates))
        expected[f"{prefix}__identifier"] = np.array(ident.to_string())

    stored = _read_collection(path)
    for prefix, (ident, _cell, _t_end) in _hdf_runs().items():
        if ident.with_rates:
            expected.update(_rebuilt_rates_arrays(prefix, stored[ident.to_string()]))
    return expected


# --------------------------------------------------------------------------- #
#                               COMPARE HELPERS                               #
# --------------------------------------------------------------------------- #


def _assert_rows_close(
    actual: tp.Dict[str, np.ndarray],
    golden: tp.Dict[str, np.ndarray],
    rtol: float,
    atol: float = 0.0,
    atol_frac_of_max: float = 0.0,
    what: str = "",
) -> None:
    """
    Compare label keyed rows with
    |actual - golden| <= atol + atol_frac_of_max * max|all golden rows| + rtol * |golden|
    and report every offending label.
    """
    assert set(actual) == set(golden), (
        f"{what}: label mismatch; only in actual: "
        f"{sorted(set(actual) - set(golden))}; only in golden: "
        f"{sorted(set(golden) - set(actual))}"
    )
    peak = max(
        (float(np.max(np.abs(v))) for v in golden.values() if v.size), default=0.0
    )
    failures = []
    for label, ref in golden.items():
        val = actual[label]
        if val.shape != ref.shape:
            failures.append(f"{label}: shape {val.shape} != {ref.shape}")
            continue
        floor = atol + atol_frac_of_max * peak
        bad = ~np.isclose(val, ref, rtol=rtol, atol=floor)
        if bad.any():
            worst = int(np.argmax(np.abs(val - ref) - floor - rtol * np.abs(ref)))
            failures.append(
                f"{label}: {int(bad.sum())} point(s) differ, worst at "
                f"#{worst}: actual={val[worst]!r} golden={ref[worst]!r}"
            )
    assert not failures, f"{what}: {len(failures)} label(s) differ:\n" + "\n".join(
        failures[:15]
    )


def _rows(
    arrays: tp.Dict[str, np.ndarray], prefix: str, kind: str
) -> tp.Dict[str, np.ndarray]:
    """{label: row} of '<prefix>__<kind>' using '<prefix>__<kind>_labels'."""
    labels = [str(lb) for lb in arrays[f"{prefix}__{kind}_labels"]]
    values = arrays[f"{prefix}__{kind}"]
    assert len(labels) == values.shape[0]
    return dict(zip(labels, values))


def _load_npz(path: Path) -> tp.Dict[str, np.ndarray]:
    if not path.exists():
        pytest.fail(f"Missing golden file {path}: run tests/golden/record_golden.py")
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def _assert_run_matches(
    actual: tp.Dict[str, np.ndarray],
    golden: tp.Dict[str, np.ndarray],
    prefix: str,
) -> None:
    """Compare the species and the reaction rates of a (sub-sampled) run."""
    np.testing.assert_allclose(
        actual[f"{prefix}__time"], golden[f"{prefix}__time"], rtol=1e-12
    )
    _assert_rows_close(
        _rows(actual, prefix, "species"),
        _rows(golden, prefix, "species"),
        rtol=RTOL_RUN,
        atol=ATOL_SPECIES_RUN,
        what=f"{prefix} species",
    )
    _assert_rows_close(
        _rows(actual, prefix, "rate"),
        _rows(golden, prefix, "rate"),
        rtol=RTOL_RUN,
        atol_frac_of_max=ATOL_RATES_RUN_FRAC_OF_MAX,
        what=f"{prefix} rates",
    )


def _assert_n_points(
    actual: tp.Dict[str, np.ndarray], golden: tp.Dict[str, np.ndarray], prefix: str
) -> None:
    n_new = int(np.sum(actual[f"{prefix}__n_points"]))
    n_ref = int(np.sum(golden[f"{prefix}__n_points"]))
    assert abs(n_new - n_ref) <= N_POINTS_RELATIVE_SPREAD * n_ref, (
        f"{prefix}: number of solver points changed a lot ({n_ref} -> {n_new})"
    )


def _assert_vector_close(
    actual: tp.Dict[str, float], golden: tp.Dict[str, float], what: str
) -> None:
    """dydt vectors: the absolute floor scales with the whole vector."""
    assert set(actual) == set(golden), f"{what}: species label mismatch"
    labels = sorted(golden)
    ref = np.array([golden[lb] for lb in labels])
    val = np.array([actual[lb] for lb in labels])
    atol = ATOL_ALGEBRA_FRAC_OF_MAX * float(np.max(np.abs(ref)))
    bad = ~np.isclose(val, ref, rtol=RTOL_ALGEBRA, atol=atol)
    assert not bad.any(), f"{what}: " + "; ".join(
        f"{lb}: actual={a!r} golden={g!r}"
        for lb, a, g in zip(np.array(labels)[bad], val[bad], ref[bad])
    )


# --------------------------------------------------------------------------- #
#                                    TESTS                                    #
# --------------------------------------------------------------------------- #


def test_derive_matches_golden() -> None:
    if not DERIVE_FILE.exists():
        pytest.fail(f"Missing golden file {DERIVE_FILE}: run record_golden.py")
    golden = json.loads(DERIVE_FILE.read_text())
    actual = compute_derive()
    assert set(actual) == set(golden)
    for config in golden:
        assert actual[config]["g_values"] == pytest.approx(
            golden[config]["g_values"], rel=RTOL_ALGEBRA
        ), f"{config}: radiolytic yields changed"
        assert set(actual[config]["dydt"]) == set(golden[config]["dydt"])
        for case, ref in golden[config]["dydt"].items():
            _assert_vector_close(
                actual[config]["dydt"][case], ref, f"derive {config}/{case}"
            )


def test_labarbe_run_matches_golden() -> None:
    golden = _load_npz(RUNS_FILE)
    actual = compute_labarbe_run()
    _assert_n_points(actual, golden, "labarbe_run")
    _assert_run_matches(actual, golden, "labarbe_run")  # first and last point
    _assert_run_matches(actual, golden, "labarbe_curve")  # whole trajectory


def test_chunked_pulsed_run_matches_golden() -> None:
    golden = _load_npz(RUNS_FILE)
    actual = compute_chunked_run()
    assert len(actual["chunked_run__n_points"]) == len(golden["chunked_run__n_points"])
    _assert_n_points(actual, golden, "chunked_run")
    _assert_run_matches(actual, golden, "chunked_run")
    # Radiolytic yields at the requested LET (Boscolo fit), label keyed
    ref_g = dict(zip(golden["chunked_run__g_labels"], golden["chunked_run__g_values"]))
    new_g = dict(zip(actual["chunked_run__g_labels"], actual["chunked_run__g_values"]))
    assert set(ref_g) == set(new_g)
    for key, val in ref_g.items():
        assert new_g[key] == pytest.approx(val, rel=RTOL_ALGEBRA), key


def test_beam_matrices_match_golden() -> None:
    golden = _load_npz(BEAMS_FILE)
    actual = compute_beams()
    assert set(actual) == set(golden)
    for key, ref in golden.items():
        np.testing.assert_allclose(
            actual[key], ref, rtol=RTOL_ALGEBRA, atol=0.0, err_msg=key
        )


def test_hdf5_written_by_old_code_is_readable(tmp_path: Path) -> None:
    """
    The committed .h5 (and its json store folder) were written by the
    pre-refactor code. Read them with the current code and compare to the
    recorded expectation. A copy is used because H5Manager opens the file
    in 'r+' mode and creates folders.
    """
    for path in (HDF_FILE, HDF_EXPECTED_FILE):
        if not path.exists():
            pytest.fail(f"Missing golden file {path}: run record_golden.py")
    expected = _load_npz(HDF_EXPECTED_FILE)

    work = tmp_path / "data"
    shutil.copytree(DATA_DIR, work)
    stored = _read_collection(work / HDF_FILE.name)

    runs = _hdf_runs()
    assert sorted(stored) == sorted(
        str(expected[f"{prefix}__identifier"]) for prefix in runs
    )

    for prefix, (ident, _cell, _t_end) in runs.items():
        res = stored[str(expected[f"{prefix}__identifier"])]

        # Stored time series (lossless float64 round trip)
        species = res.to_pandas()
        species = species[sorted(species.columns)]
        np.testing.assert_allclose(
            species.index.to_numpy(), expected[f"{prefix}__time"], rtol=RTOL_STORED
        )
        _assert_rows_close(
            dict(zip(species.columns, species.to_numpy().T)),
            _rows(expected, prefix, "species"),
            rtol=RTOL_STORED,
            what=f"{prefix} stored species",
        )
        if not ident.with_rates:
            continue

        # Stored production rates (pushed with production_rate=True)
        rates = res.reconstruct_reactions()
        rates.columns = pd.Index(_reaction_keys(res))
        rates = rates[sorted(rates.columns)]
        _assert_rows_close(
            dict(zip(rates.columns, rates.to_numpy().T)),
            _rows(expected, prefix, "rate"),
            rtol=RTOL_STORED,
            what=f"{prefix} stored rates",
        )

        # Environment rebuilt from the json store: rates recomputed on it
        rebuilt = _rebuilt_rates_arrays(prefix, res)
        _assert_rows_close(
            _rows(rebuilt, prefix, "rebuilt_rate"),
            _rows(expected, prefix, "rebuilt_rate"),
            rtol=RTOL_ALGEBRA,
            atol_frac_of_max=ATOL_ALGEBRA_FRAC_OF_MAX,
            what=f"{prefix} rebuilt rates",
        )
