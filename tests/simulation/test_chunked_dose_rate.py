#!/usr/bin/env python3
"""
Tests of the dose rate seen by the solver in a chunked run (ADR 0007): the
chunks are cut at the beam transitions, every chunk holds the (constant)
dose rate of the beam inside it, run() keeps reading the beam, and the
ODEResult reconstructs the rates with the dose rate that was integrated.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy
import typing as tp

import numpy as np
import pytest
from scipy.integrate import trapezoid

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
import radiopyo.simulation.unit_cell as unit_cell_module
from radiopyo.physics.beam import (
    ConstantBeam,
    ConstantDoseRate,
    MacroPulsedBeam,
    ParticleBeam,
    PulsedBeam,
    SinglePulseBeam,
)
from radiopyo.simulation.unit_cell import UnitCell


# ------------------------------- Fixtures ---------------------------------- #
@pytest.fixture(scope="module")
def _template() -> UnitCell:
    return UnitCell.from_toml(rp.PURE_WATER_CONFIG_FILE)


@pytest.fixture
def cell(_template: UnitCell) -> UnitCell:
    return copy.deepcopy(_template)


def _beams() -> tp.Dict[str, tp.Tuple[ParticleBeam, tp.Tuple[float, float]]]:
    """Every beam type with the time span to chunk."""
    return {
        "constant_max_dose": (
            ConstantBeam(dose_rate=1e6, max_dose=2.0, start_time=1e-6),
            (0.0, 1e-4),
        ),
        "single_pulse": (
            SinglePulseBeam(dose_rate=1e6, max_dose=2.0, start_time=1e-6),
            (0.0, 1e-4),
        ),
        "pulsed": (
            PulsedBeam.from_dose_per_pulse(
                pulse_dose=2.0, period=1e-5, on_time=2e-6, n_pulse=3.0
            ),
            (0.0, 5e-5),
        ),
        "macro_pulsed": (
            MacroPulsedBeam(
                dose_rate_peak=1e6,
                period=1e-3,
                on_time=3e-5,
                micro_period=1e-5,
                micro_on_time=2e-6,
                start_time=1e-5,
            ),
            (0.0, 1.5e-3),
        ),
    }


BEAM_IDS = list(_beams())


def _cell_with(cell: UnitCell, beam: ParticleBeam) -> UnitCell:
    cell.set_beam(beam, LET=10.0)
    return cell


# --------------------------------------------------------------------------- #
#                         Chunks follow the beam transitions                  #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", BEAM_IDS)
def test_chunk_boundaries_and_dose_rate_follow_the_beam(cell: UnitCell, name: str):
    beam, span = _beams()[name]
    cell = _cell_with(cell, beam)
    intervals = beam.make_time_intervals(start=span[0], stop=span[1])
    chunks = cell.prepare_chunked_run(t_span=span).chunks

    # boundaries are the beam transitions (make_matrix rows, ON/OFF alike)
    matrix = beam.make_matrix(start=span[0], stop=span[1])
    np.testing.assert_array_equal(intervals[:, 0], matrix[:, 0])
    assert [c.start for c in chunks] == intervals[:-1, 0].tolist()
    assert [c.stop for c in chunks] == intervals[1:, 0].tolist()
    assert [c.beam_on for c in chunks] == (intervals[:-1, 1] > 0).tolist()

    # the dose rate of a chunk is the beam dose rate, at any time inside it
    for chunk in chunks:
        assert chunk.dose_rate is not None
        for frac in (1e-6, 0.3, 0.5, 0.9, 1 - 1e-6):
            time = chunk.start + frac * (chunk.stop - chunk.start)
            assert beam.at(time).dose_rate() == chunk.dose_rate
        assert chunk.dose_rate == (beam.dose_rate_peak if chunk.beam_on else 0.0)


def test_single_pulse_chunks_dose_rate(cell: UnitCell):
    beam = SinglePulseBeam(dose_rate=1e6, max_dose=2.0)
    cell = _cell_with(cell, beam)
    chunks = cell.prepare_single_pulse(t_span=(0.0, 1e-3)).chunks
    t_max = beam.time_at_max_dose
    assert [(c.start, c.stop) for c in chunks] == [
        (0.0, t_max),
        (t_max, t_max + unit_cell_module.DEFAULT_TRANS_TIME),
        (t_max + unit_cell_module.DEFAULT_TRANS_TIME, 1e-3),
    ]
    assert [c.dose_rate for c in chunks] == [beam.dose_rate_peak, 0.0, 0.0]


# --------------------------------------------------------------------------- #
#                         What the solver is given                           #
# --------------------------------------------------------------------------- #
def _record_derive(monkeypatch: pytest.MonkeyPatch) -> tp.List[tp.Tuple[float, tp.Any]]:
    """Spy on the right hand side: records (t, beam) of every evaluation."""
    calls: tp.List[tp.Tuple[float, tp.Any]] = []
    real = unit_cell_module.derive

    def spy(t: float, y: np.ndarray, beam: tp.Any, mat: tp.Any) -> np.ndarray:
        calls.append((t, beam))
        return real(t, y, beam, mat)

    monkeypatch.setattr(unit_cell_module, "derive", spy)
    return calls


def test_chunked_solver_sees_a_constant_dose_rate_per_chunk(
    cell: UnitCell, monkeypatch: pytest.MonkeyPatch
):
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=2.0, period=1e-5, on_time=2e-6, n_pulse=2.0
    )
    cell = _cell_with(cell, beam)
    calls = _record_derive(monkeypatch)
    it = cell.prepare_chunked_run(t_span=(0.0, 2e-5), max_step_size_on=1e-7)

    for chunk in it.chunks:
        calls.clear()
        next(it)
        assert calls
        for time, source in calls:
            assert isinstance(source, ConstantDoseRate)
            assert source.at(time).dose_rate() == chunk.dose_rate
            # the evaluation times stay inside the chunk (burst shifts the start)
            assert chunk.start <= time <= chunk.stop


def test_plain_run_still_reads_the_beam(
    cell: UnitCell, monkeypatch: pytest.MonkeyPatch
):
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=2.0, period=1e-5, on_time=2e-6, n_pulse=2.0
    )
    cell = _cell_with(cell, beam)
    calls = _record_derive(monkeypatch)
    res = cell.run(t_span=(0.0, 3e-6))
    assert calls
    assert all(source is cell.beam for _, source in calls)
    assert res.step_dose_rates == [None]


def test_chunked_run_without_beam_use_passes_no_beam(
    cell: UnitCell, monkeypatch: pytest.MonkeyPatch
):
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=2.0, period=1e-5, on_time=2e-6, n_pulse=2.0
    )
    cell = _cell_with(cell, beam)
    calls = _record_derive(monkeypatch)
    res = cell.prepare_chunked_run(t_span=(0.0, 1e-5), use_beam=False).run()
    assert calls
    assert all(source is None for _, source in calls)
    assert set(res.step_dose_rates) == {None}


def test_initial_burst_uses_the_chunk_dose_rate(cell: UnitCell):
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=2.0, period=1e-5, on_time=2e-6, n_pulse=2.0
    )
    cell = _cell_with(cell, beam)
    res = cell.prepare_chunked_run(t_span=(0.0, 2e-5)).run()
    # chunks: ON, OFF, ON, OFF; one burst per ON chunk, at the chunk dose rate
    assert res.burst_dose_rates == [beam.dose_rate_peak] * 2
    assert res.step_dose_rates == [beam.dose_rate_peak, 0.0] * 2


# --------------------------------------------------------------------------- #
#                  The ODEResult reconstructs what was integrated            #
# --------------------------------------------------------------------------- #
def test_reconstructed_rates_use_the_integrated_dose_rate(cell: UnitCell):
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=2.0, period=1e-5, on_time=2e-6, n_pulse=2.0
    )
    cell = _cell_with(cell, beam)
    res = cell.prepare_chunked_run(t_span=(0.0, 2e-5), max_step_size_on=1e-7).run()
    rates = res.reconstruct_reactions(recompute=True)
    label = next(
        r.as_label() for r in res.env.reactions.radiolytic if r.species().label == "e_aq"
    )
    rad = rates[label]
    rate_on = rad.max()
    assert rate_on > 0

    ends = np.cumsum(res.n_points) - 1
    starts = np.concatenate([[0], ends[:-1] + 1])
    for i, (first, last) in enumerate(zip(starts, ends, strict=True)):
        on = i % 2 == 0  # chunks: ON, OFF, ON, OFF
        expected = rate_on if on else 0.0
        # every point of the chunk, boundaries included (the last point of an
        # ON chunk lies on the pulse end where the beam is already OFF)
        np.testing.assert_allclose(rad.to_numpy()[first : last + 1], expected, rtol=1e-12)

    # integrated radiolytic production of an ON chunk = rate * its duration
    # (the first T_ZERO of the chunk is skipped, added analytically as the
    # burst); a chunk cut short of its last point would lose half a step
    values, times = rad.to_numpy(), rad.index.to_numpy()
    for i in (0, 2):
        chunk = slice(starts[i], ends[i] + 1)
        integral = trapezoid(values[chunk], x=times[chunk])
        assert integral == pytest.approx(rate_on * (2e-6 - cell.T_ZERO), rel=1e-9)


def test_step_dose_rates_follow_concatenation(cell: UnitCell):
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=2.0, period=1e-5, on_time=2e-6, n_pulse=2.0
    )
    cell = _cell_with(cell, beam)
    first = cell.prepare_chunked_run(t_span=(0.0, 1e-5)).run()
    second = cell.run(t_span=(1e-5, 1.1e-5))
    total = first + second
    assert total.step_dose_rates == [*first.step_dose_rates, None]
    assert len(total.step_dose_rates) == len(total.n_points)


@pytest.mark.parametrize(
    "beam",
    [
        SinglePulseBeam(dose_rate=1e6, max_dose=2.0, start_time=1e-3),
        ConstantBeam(dose_rate=1e6, max_dose=2.0, start_time=1e-3),
    ],
    ids=["single_pulse", "constant"],
)
def test_single_pulse_chunks_with_beam_starting_after_t_span(
    cell: UnitCell, beam: ParticleBeam
):
    cell = _cell_with(cell, beam)
    t_span = (0.0, 1e-2)
    chunks = cell.prepare_single_pulse(t_span=t_span).chunks
    t_start, t_max = beam.timer.start_time, beam.time_at_max_dose
    t_off = t_max + unit_cell_module.DEFAULT_TRANS_TIME
    assert [(c.start, c.stop) for c in chunks] == [
        (0.0, t_start),
        (t_start, t_max),
        (t_max, t_off),
        (t_off, t_span[1]),
    ]
    assert [c.beam_on for c in chunks] == [False, True, False, False]
    assert [c.dose_rate for c in chunks] == [0.0, beam.dose_rate_peak, 0.0, 0.0]
    # every chunk sits in one beam interval: the dose rate is constant in it
    for chunk in chunks:
        for frac in (1e-6, 0.5, 1 - 1e-6):
            time = chunk.start + frac * (chunk.stop - chunk.start)
            assert beam.at(time).dose_rate() == chunk.dose_rate

    # the pulse delivers its dose: the ON chunk is integrated with the beam
    # dose rate and a burst at that rate
    res = cell.prepare_single_pulse(t_span=t_span).run()
    assert res.burst_dose_rates == [beam.dose_rate_peak]
    assert res.step_dose_rates == [0.0, beam.dose_rate_peak, 0.0, 0.0]


# --------------------------------------------------------------------------- #
#             prepare_single_pulse: pre-tm chunks are the beam intervals       #
# --------------------------------------------------------------------------- #
def _multi_pulse_beams() -> tp.Dict[str, ParticleBeam]:
    return {
        "pulsed": PulsedBeam.from_dose_per_pulse(
            pulse_dose=2.0, period=1e-5, on_time=2e-6, n_pulse=3.0
        ),
        "pulsed_start_time": PulsedBeam.from_dose_per_pulse(
            pulse_dose=2.0,
            period=1e-5,
            on_time=2e-6,
            n_pulse=3.0,
            start_time=1e-4,
        ),
        "macro_start_time": MacroPulsedBeam(
            dose_rate_peak=1e6,
            period=1e-3,
            on_time=3e-5,
            micro_period=1e-5,
            micro_on_time=2e-6,
            start_time=1e-4,
            max_dose=10.0,
        ),
    }


@pytest.mark.parametrize("name", list(_multi_pulse_beams()))
def test_single_pulse_chunks_follow_a_multi_pulse_beam(cell: UnitCell, name: str):
    beam = _multi_pulse_beams()[name]
    cell = _cell_with(cell, beam)
    t_span = (0.0, 3e-3)
    chunks = cell.prepare_single_pulse(t_span=t_span).chunks
    t_max = beam.time_at_max_dose
    assert t_max < t_span[1]

    # contiguous, inside t_span, and the pre-tm chunks are the beam intervals
    assert chunks[0].start == t_span[0]
    assert chunks[-1].stop == t_span[1]
    for left, right in zip(chunks[:-1], chunks[1:], strict=True):
        assert left.stop == right.start
    assert all(c.stop > c.start for c in chunks)
    intervals = beam.make_time_intervals(start=t_span[0], stop=t_max)
    n_pre = len(intervals) - 1
    assert [c.start for c in chunks[:n_pre]] == intervals[:-1, 0].tolist()
    assert [c.stop for c in chunks[:n_pre]] == intervals[1:, 0].tolist()
    assert [c.beam_on for c in chunks[:n_pre]] == (intervals[:-1, 1] > 0).tolist()
    assert sum(c.beam_on for c in chunks) > 1 or name == "macro_start_time"

    # the dose rate of a chunk is the beam one, at any time inside it
    for chunk in chunks:
        assert chunk.dose_rate is not None
        for frac in (1e-6, 0.3, 0.5, 0.9, 1 - 1e-6):
            time = chunk.start + frac * (chunk.stop - chunk.start)
            assert beam.at(time).dose_rate() == chunk.dose_rate, (chunk, time)


@pytest.mark.parametrize("name", list(_multi_pulse_beams()))
def test_single_pulse_delivers_the_dose_of_the_chunked_run(
    _template: UnitCell, name: str
):
    beam = _multi_pulse_beams()[name]
    t_span = (0.0, 3e-3)
    finals = []
    for method in ("prepare_chunked_run", "prepare_single_pulse"):
        cell = _cell_with(copy.deepcopy(_template), copy.deepcopy(beam))
        finals.append(getattr(cell, method)(t_span=t_span).run().y[:, -1])
    np.testing.assert_allclose(finals[1], finals[0], rtol=1e-5, atol=1e-9)
    # something was delivered (the species are not left at their initial value)
    initial = _cell_with(copy.deepcopy(_template), beam).env.initial_values()
    assert np.any(np.abs(finals[1] - initial) > 1e-3)


def test_single_pulse_chunks_are_clipped_to_t_span(cell: UnitCell):
    beam = SinglePulseBeam(dose_rate=1e3, max_dose=1.0)  # tm = 1e-3 s
    cell = _cell_with(cell, beam)
    trans = unit_cell_module.DEFAULT_TRANS_TIME

    # t_span ends inside the pulse: one (clipped) ON chunk
    chunks = cell.prepare_single_pulse(t_span=(0.0, 5e-4)).chunks
    assert [(c.start, c.stop, c.beam_on) for c in chunks] == [(0.0, 5e-4, True)]
    res = cell.prepare_single_pulse(t_span=(0.0, 5e-4)).run()
    assert res.t[-1] == 5e-4
    assert np.all(np.diff(res.t) >= 0)

    # t_span ends inside the transition
    chunks = cell.prepare_single_pulse(t_span=(0.0, 1e-3 + trans / 2)).chunks
    assert [(c.start, c.stop) for c in chunks] == [
        (0.0, 1e-3),
        (1e-3, 1e-3 + trans / 2),
    ]

    # t_span ends exactly at tm: no empty transition chunk
    chunks = cell.prepare_single_pulse(t_span=(0.0, 1e-3)).chunks
    assert [(c.start, c.stop) for c in chunks] == [(0.0, 1e-3)]


def test_single_pulse_off_chunks_before_tm_use_the_transition_max_step(
    cell: UnitCell,
):
    beam = _multi_pulse_beams()["pulsed_start_time"]
    cell = _cell_with(cell, beam)
    on, trans, final = 1e-8, 2e-7, 3e-6
    chunks = cell.prepare_single_pulse(
        t_span=(0.0, 3e-3),
        max_step_size_on=on,
        max_step_size_trans=trans,
        max_step_size_final=final,
    ).chunks
    t_max = beam.time_at_max_dose
    pre = [c for c in chunks if c.stop <= t_max]
    # leading OFF chunk, ON pulses and the OFF gaps between them
    assert not pre[0].beam_on and pre[0].start == 0.0
    assert any(not c.beam_on for c in pre[1:])
    for c in pre:
        assert c.max_step_size == (on if c.beam_on else trans), c
    # ON -> OFF transition, then the final OFF chunk
    assert [c.max_step_size for c in chunks[len(pre) :]] == [trans, final]
