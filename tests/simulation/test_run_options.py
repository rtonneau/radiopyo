#!/usr/bin/env python3
"""
Tests of the UnitCell run flow: RunOptions, the named steps of _raw_run
(span widening, solver kwargs, initial burst), the chunk building of
prepare_chunked_run / prepare_single_pulse, SimChunksIter, and the
public keyword arguments and defaults of the run methods.

run() applies the initial radiolytic burst on beam ON (unless the beam
is not used), forwards force_positive to the solver, never modifies the
y0 it is given and raises NoBeamDefinedError without beam.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import dataclasses
import inspect
import typing as tp

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
import radiopyo.simulation.unit_cell as unit_cell_module
from radiopyo.physics.beam import PulsedBeam, SinglePulseBeam
from radiopyo.simulation.exceptions import NoBeamDefinedError
from radiopyo.simulation.ode_result import ODEResult
from radiopyo.simulation.unit_cell import (
    DEFAULT_ATOL,
    DEFAULT_RTOL,
    DEFAULT_TRANS_TIME,
    RunOptions,
    SimChunk,
    SimChunksIter,
    UnitCell,
    widen_empty_span,
)


# ------------------------------- Fixtures ---------------------------------- #
# Parsing a config takes ~1.5 s: parse once per module, hand out deep copies
# (some tests change the beam, hence the G values of the env).
@pytest.fixture(scope="module")
def _labarbe_template() -> UnitCell:
    return UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)


@pytest.fixture(scope="module")
def _pulsed_template() -> UnitCell:
    cell = UnitCell.from_toml(rp.PURE_WATER_CONFIG_FILE)
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=2.0, period=1e-2, on_time=2e-6, n_pulse=3.0
    )
    cell.set_beam(beam, LET=10.0)
    return cell


@pytest.fixture
def labarbe_cell(_labarbe_template: UnitCell) -> UnitCell:
    """Labarbe config: constant beam of 1e7 Gy/s up to 10 Gy (ON [0, 1e-6])."""
    return copy.deepcopy(_labarbe_template)


@pytest.fixture
def pulsed_cell(_pulsed_template: UnitCell) -> UnitCell:
    """Pure water with 3 pulses of 2 Gy (on_time 2e-6 s, period 1e-2 s)."""
    return copy.deepcopy(_pulsed_template)


@pytest.fixture
def no_beam_cell(labarbe_cell: UnitCell) -> UnitCell:
    return UnitCell(labarbe_cell.env)


class _Captured(Exception):
    """Raised by the _raw_run spy to stop run() once its arguments are known."""


def _spy_raw_run(
    monkeypatch: pytest.MonkeyPatch, cell: UnitCell
) -> tp.List[tp.Dict[str, tp.Any]]:
    """Replace cell._raw_run by a spy recording its arguments and stopping."""
    calls: tp.List[tp.Dict[str, tp.Any]] = []

    def spy(
        t_span: tp.Tuple[float, float],
        options: RunOptions,
        y0: tp.Any = None,
        beam_on: bool = True,
    ) -> tp.Any:
        calls.append(
            {"t_span": t_span, "options": options, "y0": y0, "beam_on": beam_on}
        )
        raise _Captured

    monkeypatch.setattr(cell, "_raw_run", spy)
    return calls


# --------------------------------------------------------------------------- #
#                                  RunOptions                                 #
# --------------------------------------------------------------------------- #
def test_run_options_defaults() -> None:
    opts = RunOptions()
    assert opts.atol is None
    assert opts.rtol is None
    assert opts.method == "LSODA"
    assert opts.max_step is None
    assert opts.first_step is None
    assert opts.use_beam is True
    assert opts.force_positive is False
    assert opts.use_jac is False
    assert opts.t_eval is None


def test_run_options_frozen_and_keyword_only() -> None:
    opts = RunOptions(atol=1e-9)
    with pytest.raises(dataclasses.FrozenInstanceError):
        opts.atol = 1.0  # type: ignore[misc]
    with pytest.raises(TypeError):
        RunOptions(1e-9)  # type: ignore[misc]
    assert dataclasses.replace(opts, max_step=1e-3).max_step == 1e-3
    assert opts.max_step is None


def test_default_constants() -> None:
    assert DEFAULT_ATOL == 1e-12
    assert DEFAULT_RTOL == 1e-10
    assert DEFAULT_TRANS_TIME == 100e-6
    assert UnitCell.FIRST_STEP == 1e-15
    assert UnitCell.T_ZERO == 1e-7


# --------------------------------------------------------------------------- #
#                               _raw_run steps                                #
# --------------------------------------------------------------------------- #
def test_widen_empty_span() -> None:
    assert widen_empty_span((1e-6, 1e-6)) == (1e-8, 1e-6)
    assert widen_empty_span((1e-7, 1.0)) == (1e-7, 1.0)
    assert widen_empty_span((0.0, 0.0)) == (0.0, 0.0)


def test_solver_kwargs_defaults(labarbe_cell: UnitCell) -> None:
    kwargs = labarbe_cell._solver_kwargs((1e-7, 1.0), RunOptions())
    assert kwargs == {
        "first_step": UnitCell.FIRST_STEP,
        "rtol": DEFAULT_RTOL,
        "atol": DEFAULT_ATOL,
    }


def test_solver_kwargs_user_values(labarbe_cell: UnitCell) -> None:
    opts = RunOptions(atol=1e-9, rtol=1e-7, max_step=1e-3, first_step=1e-9)
    kwargs = labarbe_cell._solver_kwargs((1e-7, 1.0), opts)
    assert kwargs == {"first_step": 1e-9, "rtol": 1e-7, "atol": 1e-9, "max_step": 1e-3}


def test_solver_kwargs_first_step_clamped_by_max_step(
    labarbe_cell: UnitCell,
) -> None:
    labarbe_cell.FIRST_STEP = 1e-4  # instance override is honoured
    kwargs = labarbe_cell._solver_kwargs((1e-7, 1.0), RunOptions(max_step=1e-5))
    assert kwargs["first_step"] == 1e-5 / 10
    # first_step == max_step is clamped too
    kwargs = labarbe_cell._solver_kwargs((1e-7, 1.0), RunOptions(max_step=1e-4))
    assert kwargs["first_step"] == 1e-4 / 10


def test_solver_kwargs_first_step_clamped_by_span(labarbe_cell: UnitCell) -> None:
    span = (1e-7, 1e-3)
    kwargs = labarbe_cell._solver_kwargs(span, RunOptions(first_step=1.0))
    assert kwargs["first_step"] == (span[1] - span[0]) / 1e3


def test_apply_initial_burst(labarbe_cell: UnitCell) -> None:
    cell = labarbe_cell
    y0 = cell.env.initial_values()
    before = y0.copy()
    span = cell._apply_initial_burst((0.0, 1.0), y0)
    assert span == (0.0 + cell.T_ZERO, 1.0)

    expected = before.copy()
    for r in cell.env.reactions.radiolytic:
        if cell.env.species.is_cst_species(r.species()):
            continue
        expected[r.species().index] += (
            r.compute_initial_cc_before(
                t0=cell.T_ZERO, dose_rate=cell.beam.at(span[0]).dose_rate()
            )
            * 1e6
        )
    np.testing.assert_array_equal(y0, expected)  # modified in place
    assert np.any(y0 != before)


def test_apply_initial_burst_too_short(labarbe_cell: UnitCell) -> None:
    with pytest.raises(ValueError, match="Too Short"):
        labarbe_cell._apply_initial_burst(
            (1e-7, 1.5e-7), labarbe_cell.env.initial_values()
        )


def test_raw_run_uses_options(
    labarbe_cell: UnitCell, monkeypatch: pytest.MonkeyPatch
) -> None:
    """_raw_run hands every option to its named steps and to solve_ivp."""
    seen: tp.Dict[str, tp.Any] = {}

    def fake_solve(
        t_span: tp.Tuple[float, float],
        y0: np.ndarray,
        beam: tp.Any,
        options: RunOptions,
        solver_kwargs: tp.Dict[str, tp.Any],
    ) -> str:
        seen.update(
            t_span=t_span, y0=y0, beam=beam, options=options, kwargs=solver_kwargs
        )
        return "sol"

    monkeypatch.setattr(labarbe_cell, "_timed_solve", fake_solve)
    opts = RunOptions(max_step=1e-3)
    y0 = labarbe_cell.env.initial_values()
    y0_ref = y0.copy()
    assert labarbe_cell._raw_run((1e-6, 1e-6), opts, y0=y0, beam_on=True) == "sol"
    # empty span widened, then shifted by T_ZERO (burst)
    assert seen["t_span"] == (1e-8 + labarbe_cell.T_ZERO, 1e-6)
    assert seen["y0"] is not y0  # a copy carries the burst
    np.testing.assert_array_equal(y0, y0_ref)
    assert np.any(seen["y0"] != y0_ref)
    assert seen["beam"] is labarbe_cell.beam
    assert seen["options"] is opts
    assert seen["kwargs"]["max_step"] == 1e-3

    # beam not used: no burst, no beam
    labarbe_cell._raw_run((1e-7, 1.0), RunOptions(use_beam=False), y0=y0, beam_on=True)
    assert seen["t_span"] == (1e-7, 1.0)
    assert seen["beam"] is None
    np.testing.assert_array_equal(seen["y0"], y0_ref)
    assert labarbe_cell.last_burst_dose_rate is None

    labarbe_cell._raw_run((1e-7, 1.0), RunOptions(), beam_on=False)
    assert seen["t_span"] == (1e-7, 1.0)
    assert seen["beam"] is labarbe_cell.beam
    np.testing.assert_array_equal(seen["y0"], labarbe_cell.env.initial_values())


# --------------------------------------------------------------------------- #
#                                  UnitCell.run                               #
# --------------------------------------------------------------------------- #
def test_run_forwards_options(
    labarbe_cell: UnitCell, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = _spy_raw_run(monkeypatch, labarbe_cell)
    t_eval = np.array([1e-6, 1e-3])
    y0 = labarbe_cell.env.initial_values()
    with pytest.raises(_Captured):
        labarbe_cell.run(
            t_span=(1e-7, 1.0),
            t_eval=t_eval,
            method="BDF",
            atol=1e-9,
            rtol=1e-7,
            max_step=1e-3,
            y0=y0,
            use_beam=False,
        )
    (call,) = calls
    assert call["t_span"] == (1e-7, 1.0)
    assert call["y0"] is y0
    opts = call["options"]
    assert (opts.atol, opts.rtol, opts.method, opts.max_step) == (
        1e-9,
        1e-7,
        "BDF",
        1e-3,
    )
    assert opts.use_beam is False
    assert opts.use_jac is False
    assert opts.first_step is None
    assert opts.t_eval is t_eval


@pytest.mark.parametrize("force_positive", [False, True])
def test_run_forwards_force_positive_to_options(
    labarbe_cell: UnitCell, monkeypatch: pytest.MonkeyPatch, force_positive: bool
) -> None:
    """force_positive lands in RunOptions, whatever the burst does."""
    calls = _spy_raw_run(monkeypatch, labarbe_cell)
    with pytest.raises(_Captured):
        labarbe_cell.run(t_span=(1e-7, 1.0), force_positive=force_positive)
    (call,) = calls
    assert call["options"].force_positive is force_positive
    assert call["beam_on"] is True  # independent of force_positive


@pytest.mark.parametrize("force_positive", [False, True])
def test_force_positive_reaches_derive(
    labarbe_cell: UnitCell, monkeypatch: pytest.MonkeyPatch, force_positive: bool
) -> None:
    """The args given by solve_ivp to derive() carry force_positive."""
    seen: tp.Dict[str, tp.Any] = {}

    def fake_solve_ivp(
        fun: tp.Any, t_span: tp.Any, y0: tp.Any, **kwargs: tp.Any
    ) -> str:
        seen["args"] = kwargs["args"]
        return "sol"

    monkeypatch.setattr(unit_cell_module, "solve_ivp", fake_solve_ivp)
    labarbe_cell._raw_run(
        (1e-7, 1.0), RunOptions(force_positive=force_positive), beam_on=True
    )
    beam, mat, forced = seen["args"]
    assert beam is labarbe_cell.beam
    assert mat is labarbe_cell.mat
    assert forced is force_positive


def test_run_applies_burst_by_default(labarbe_cell: UnitCell) -> None:
    """End to end: the first stored point is the burst shifted state."""
    cell = labarbe_cell
    res = cell.run(t_span=(1e-7, 1e-6))
    y_burst = cell.env.initial_values()
    cell._apply_initial_burst((1e-7, 1e-6), y_burst)
    np.testing.assert_array_equal(res.y[:, 0], y_burst)
    assert res.t[0] == 1e-7 + cell.T_ZERO
    assert res.burst_dose_rates == [cell.beam.dose_rate_peak]


def test_run_without_beam_use_has_no_burst(labarbe_cell: UnitCell) -> None:
    cell = labarbe_cell
    res = cell.run(t_span=(1e-7, 1e-6), use_beam=False)
    np.testing.assert_array_equal(res.y[:, 0], cell.env.initial_values())
    assert res.t[0] == 1e-7
    assert res.burst_dose_rates == []


def test_run_never_mutates_y0(labarbe_cell: UnitCell) -> None:
    y0 = labarbe_cell.env.initial_values()
    y0_ref = y0.copy()
    res = labarbe_cell.run(t_span=(1e-7, 1e-6), y0=y0)
    np.testing.assert_array_equal(y0, y0_ref)
    assert np.any(res.y[:, 0] != y0_ref)  # the burst was applied to a copy


def test_run_drops_t_eval_inside_skipped_interval(labarbe_cell: UnitCell) -> None:
    t_eval = np.array([1e-7, 1.5e-7, 2e-7, 5e-7])
    res = labarbe_cell.run(t_span=(1e-7, 1e-6), t_eval=t_eval)
    np.testing.assert_array_equal(res.t, t_eval[2:])
    # beam not used: nothing is skipped
    res = labarbe_cell.run(t_span=(1e-7, 1e-6), t_eval=t_eval, use_beam=False)
    np.testing.assert_array_equal(res.t, t_eval)


def test_run_without_beam_raises(no_beam_cell: UnitCell) -> None:
    with pytest.raises(NoBeamDefinedError):
        no_beam_cell.run(t_span=(1e-7, 1.0))
    with pytest.raises(NoBeamDefinedError):
        no_beam_cell.run(t_span=(1e-7, 1.0), use_beam=False)
    with pytest.raises(NoBeamDefinedError):
        no_beam_cell._raw_run((1e-7, 1.0), RunOptions())


def test_integrate_production_follows_applied_burst(labarbe_cell: UnitCell) -> None:
    """The burst term of the production integrals exists only if the burst
    was applied to the initial concentrations."""
    cell = labarbe_cell
    with_burst = cell.run(t_span=(1e-7, 1e-6))
    without = cell.run(t_span=(1e-7, 1e-6), use_beam=False)
    (reaction,) = cell.env.radiolytic_reaction_involving("e_aq")
    label = reaction.as_label()
    burst = (
        reaction.compute_initial_cc_before(cell.T_ZERO, cell.beam.dose_rate_peak) * 1e6
    )

    def simulated(res: ODEResult) -> float:
        rate = res.reconstruct_reactions()[label]
        return float(np.trapezoid(rate.to_numpy() * 1e6, x=rate.index))

    # Burst applied: integral of the rate over the simulated part + burst
    prod = with_burst.integrate_production_of("e_aq")[label]
    assert prod == pytest.approx(simulated(with_burst) + burst)
    # Not applied (use_beam=False): the rate integral only
    prod_without = without.integrate_production_of("e_aq")[label]
    assert prod_without == pytest.approx(simulated(without))

    # Unknown burst information (result read from a file): one burst
    # per beam pulse at the peak dose rate.
    with_burst.burst_dose_rates = None
    legacy = with_burst.integrate_production_of("e_aq")[label]
    assert legacy == pytest.approx(
        simulated(with_burst) + burst * cell.beam.tot_pulses()
    )


def test_run_till_max_dose_forwards_kwargs(
    labarbe_cell: UnitCell, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = _spy_raw_run(monkeypatch, labarbe_cell)
    with pytest.raises(_Captured):
        labarbe_cell.run_till_max_dose(
            time_start=1e-8,
            t_span=(5.0, 6.0),  # overwritten
            max_step=1e-9,
            force_positive=True,
        )
    (call,) = calls
    assert call["t_span"] == (1e-8, labarbe_cell.beam.time_at_max_dose)
    assert call["options"].max_step == 1e-9
    assert call["options"].force_positive is True
    assert call["beam_on"] is True


# --------------------------------------------------------------------------- #
#                              Public signatures                              #
# --------------------------------------------------------------------------- #
_EMPTY = inspect.Parameter.empty
_COMMON_CHUNKED = [
    ("use_jac", False),
    ("method", "LSODA"),
    ("atol", None),
    ("rtol", None),
    ("max_step_size_on", None),
]


@pytest.mark.parametrize(
    "method, expected",
    [
        (
            "run",
            [
                ("t_span", _EMPTY),
                ("t_eval", None),
                ("use_jac", False),
                ("method", "LSODA"),
                ("atol", None),
                ("rtol", None),
                ("max_step", None),
                ("y0", None),
                ("use_beam", True),
                ("force_positive", False),
            ],
        ),
        (
            "prepare_chunked_run",
            [
                ("t_span", _EMPTY),
                *_COMMON_CHUNKED,
                ("max_step_size_off", None),
                ("max_step_size_final", None),
                ("y0", None),
                ("use_beam", True),
                ("force_positive", False),
            ],
        ),
        (
            "prepare_single_pulse",
            [
                ("t_span", _EMPTY),
                *_COMMON_CHUNKED,
                ("max_step_size_trans", None),
                ("max_step_size_final", None),
                ("trans_time", None),
                ("y0", None),
                ("use_beam", True),
                ("force_positive", False),
            ],
        ),
        ("run_till_max_dose", [("time_start", 1e-9), ("kwargs", _EMPTY)]),
    ],
)
def test_public_signatures_unchanged(
    method: str, expected: tp.List[tp.Tuple[str, tp.Any]]
) -> None:
    params = list(inspect.signature(getattr(UnitCell, method)).parameters.values())
    assert params[0].name == "self"
    assert [(p.name, p.default) for p in params[1:]] == expected
    for p in params[1:]:
        assert p.kind in (
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.VAR_KEYWORD,
        )


# --------------------------------------------------------------------------- #
#                                Chunked runs                                 #
# --------------------------------------------------------------------------- #
def _chunk_table(chunks: tp.List[SimChunk]) -> tp.List[tp.Tuple[tp.Any, ...]]:
    return [(c.start, c.stop, c.max_step_size, bool(c.beam_on)) for c in chunks]


def test_prepare_chunked_run_chunks(pulsed_cell: UnitCell) -> None:
    t_span = (0.0, 5e-2)
    intervals = pulsed_cell.beam.make_time_intervals(start=0.0, stop=5e-2)
    it = pulsed_cell.prepare_chunked_run(
        t_span=t_span, max_step_size_on=1e-7, max_step_size_off=1e-3
    )
    on, off = 1e-7, 1e-3
    expected = [
        (
            intervals[i, 0],
            intervals[i + 1, 0],
            on if intervals[i, 1] > 0 else off,
            bool(intervals[i, 1] > 0),
        )
        for i in range(len(intervals) - 1)
    ]
    assert _chunk_table(it.chunks) == expected
    assert len(it) == len(intervals) - 1
    # 3 pulses: ON, OFF, ON, OFF, ON, OFF
    assert [c.beam_on for c in it.chunks] == [True, False] * 3


def test_prepare_chunked_run_final_max_step(pulsed_cell: UnitCell) -> None:
    it = pulsed_cell.prepare_chunked_run(
        t_span=(0.0, 5e-2),
        max_step_size_on=1e-7,
        max_step_size_off=1e-3,
        max_step_size_final=1e-2,
    )
    assert [c.max_step_size for c in it.chunks] == [1e-7, 1e-3] * 2 + [1e-7, 1e-2]


def test_prepare_chunked_run_options(pulsed_cell: UnitCell) -> None:
    y0 = pulsed_cell.env.initial_values()
    it = pulsed_cell.prepare_chunked_run(
        t_span=(0.0, 5e-2),
        method="BDF",
        atol=1e-9,
        rtol=1e-7,
        y0=y0,
        use_beam=False,
        force_positive=True,
    )
    assert it.options == RunOptions(
        atol=1e-9, rtol=1e-7, method="BDF", use_beam=False, force_positive=True
    )
    assert it.y0 is y0


def test_prepare_single_pulse_chunks(labarbe_cell: UnitCell) -> None:
    t_max = labarbe_cell.beam.time_at_max_dose
    it = labarbe_cell.prepare_single_pulse(
        t_span=(0.0, 1.0), max_step_size_on=1e-8, max_step_size_trans=1e-6
    )
    assert _chunk_table(it.chunks) == [
        (0.0, t_max, 1e-8, True),
        (t_max, t_max + DEFAULT_TRANS_TIME, 1e-6, False),
        (t_max + DEFAULT_TRANS_TIME, 1.0, None, False),  # no fallback
    ]
    assert it.options == RunOptions()


def test_prepare_single_pulse_trans_time_and_options(labarbe_cell: UnitCell) -> None:
    cell = labarbe_cell
    cell.set_beam(SinglePulseBeam(dose_rate=1e6, max_dose=10.0), LET=0.2)
    t_max = cell.beam.time_at_max_dose
    it = cell.prepare_single_pulse(
        t_span=(0.0, 1.0),
        max_step_size_final=1e-3,
        trans_time=1e-5,
        atol=1e-9,
        use_jac=True,
        force_positive=True,
    )
    assert _chunk_table(it.chunks) == [
        (0.0, t_max, None, True),
        (t_max, t_max + 1e-5, None, False),
        (t_max + 1e-5, 1.0, 1e-3, False),
    ]
    assert it.options == RunOptions(atol=1e-9, use_jac=True, force_positive=True)


@pytest.mark.parametrize("method", ["prepare_chunked_run", "prepare_single_pulse"])
def test_prepare_without_beam_raises(no_beam_cell: UnitCell, method: str) -> None:
    with pytest.raises(NoBeamDefinedError):
        getattr(no_beam_cell, method)(t_span=(0.0, 1.0))


def test_sim_chunks_iter_calls(
    labarbe_cell: UnitCell, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Each chunk: its span, its max step and beam state, the shared
    options, and the last point of the previous chunk as y0."""
    cell = labarbe_cell
    calls: tp.List[tp.Dict[str, tp.Any]] = []
    real_raw_run = cell._raw_run

    def spy(
        t_span: tp.Tuple[float, float],
        options: RunOptions,
        y0: tp.Any = None,
        beam_on: bool = True,
    ) -> tp.Any:
        calls.append(
            {
                "t_span": t_span,
                "options": options,
                "y0": None if y0 is None else y0.copy(),
                "beam_on": beam_on,
            }
        )
        return real_raw_run(t_span, options, y0=y0, beam_on=beam_on)

    monkeypatch.setattr(cell, "_raw_run", spy)
    chunks = [
        SimChunk(start=0.0, stop=1e-6, max_step_size=1e-8, beam_on=True),
        SimChunk(start=1e-6, stop=1e-4, max_step_size=None, beam_on=False),
    ]
    opts = RunOptions(rtol=1e-8, max_step=5.0)
    it = SimChunksIter(sim=cell, chunks=chunks, options=opts)

    first = next(it)
    n_first = first.y.shape[1]
    last_of_first = first.y[:, -1].copy()
    res = next(it)
    with pytest.raises(StopIteration):
        next(it)

    assert [c["t_span"] for c in calls] == [(0.0, 1e-6), (1e-6, 1e-4)]
    assert [c["beam_on"] for c in calls] == [True, False]
    assert [c["options"] for c in calls] == [
        dataclasses.replace(opts, max_step=1e-8),
        dataclasses.replace(opts, max_step=None),
    ]
    assert calls[0]["y0"] is None
    np.testing.assert_array_equal(calls[1]["y0"], last_of_first)
    assert res.n_points == [n_first, res.y.shape[1] - n_first]
    assert it.options is opts  # never replaced


def test_sim_chunks_iter_default_options(labarbe_cell: UnitCell) -> None:
    it = SimChunksIter(sim=labarbe_cell, chunks=[])
    assert it.options == RunOptions()
    assert it.y0 is None
    assert len(it) == 0


def test_sim_chunks_iter_keeps_previous_last_point(pulsed_cell: UnitCell) -> None:
    """The burst of an ON chunk starts from a copy of the previous last
    point: the stored value is not overwritten. Chunks are ON, OFF, ON."""
    cell = pulsed_cell
    it = cell.prepare_chunked_run(
        t_span=(0.0, 1e-2 + 2e-6),
        max_step_size_on=1e-7,
        max_step_size_off=1e-3,
    )
    assert [c.beam_on for c in it.chunks] == [True, False, True]
    next(it)
    res = next(it)
    n_stored = res.y.shape[1]
    last_off = res.y[:, -1].copy()
    res = next(it)
    np.testing.assert_array_equal(res.y[:, n_stored - 1], last_off)
    # ... and the burst really was added to the next chunk's first point
    assert cell.last_burst_dose_rate == cell.beam.dose_rate_peak
    assert np.any(res.y[:, n_stored] != last_off)
    assert res.burst_dose_rates == [cell.beam.dose_rate_peak] * 2


def test_sim_chunks_iter_never_mutates_y0(pulsed_cell: UnitCell) -> None:
    y0 = pulsed_cell.env.initial_values()
    y0_ref = y0.copy()
    pulsed_cell.prepare_chunked_run(t_span=(0.0, 3e-6), y0=y0).run()
    np.testing.assert_array_equal(y0, y0_ref)


def test_chunked_run_without_beam_use_has_no_burst(pulsed_cell: UnitCell) -> None:
    res = pulsed_cell.prepare_chunked_run(t_span=(0.0, 3e-6), use_beam=False).run()
    assert res.burst_dose_rates == []
    assert res.t[0] == 0.0
    np.testing.assert_array_equal(res.y[:, 0], pulsed_cell.env.initial_values())
