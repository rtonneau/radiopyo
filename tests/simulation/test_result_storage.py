#!/usr/bin/env python3
"""
Tests of the per-step storage of ODEResult: appending a step never copies
the earlier ones (the cost of a chunked run stays linear in the number of
chunks), t / y join the steps lazily, and the results are unchanged.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""
# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import copy

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.physics.beam import PulsedBeam
from radiopyo.simulation.ode_result import DummyScipyOdeResult, ODEResult
from radiopyo.simulation.unit_cell import UnitCell


# ------------------------------- Fixtures ---------------------------------- #
@pytest.fixture(scope="module")
def _template() -> UnitCell:
    cell = UnitCell.from_toml(rp.PURE_WATER_CONFIG_FILE)
    cell.set_beam(
        PulsedBeam.from_dose_per_pulse(
            pulse_dose=1.0, period=1e-4, on_time=2e-6, n_pulse=3.0
        ),
        LET=10.0,
    )
    return cell


@pytest.fixture
def cell(_template: UnitCell) -> UnitCell:
    return copy.deepcopy(_template)


def _empty(cell: UnitCell) -> ODEResult:
    return ODEResult(
        cell.env, cell.beam_list.current, cell.mat, None, cell.T_ZERO
    )


def _step(t0: float, n: int, n_sp: int) -> DummyScipyOdeResult:
    t = np.linspace(t0, t0 + 1.0, n)
    y = np.arange(n_sp * n, dtype=np.float64).reshape(n_sp, n) + t0
    return DummyScipyOdeResult(t=t, y=y)


# --------------------------------- Tests ----------------------------------- #
def test_append_does_not_copy_earlier_steps(cell: UnitCell) -> None:
    res = _empty(cell)
    n_sp = len(cell.env.initial_values())
    first = _step(0.0, 5, n_sp)
    res.append_sol(first)
    for i in range(1, 50):
        res.append_sol(_step(float(i), 5, n_sp))
        # The first step is still the very array that was appended.
        assert res._y_blocks[0] is first.y
    assert len(res._y_blocks) == 50
    assert res.y.shape == (n_sp, 250)
    assert res.t.shape == (250,)


def test_read_joins_and_append_invalidates(cell: UnitCell) -> None:
    res = _empty(cell)
    n_sp = len(cell.env.initial_values())
    steps = [_step(float(i), 4, n_sp) for i in range(3)]
    for s in steps[:2]:
        res.append_sol(s)
    y = res.y
    assert res.y is y  # cached: second read is free
    np.testing.assert_array_equal(y, np.hstack([s.y for s in steps[:2]]))
    res.append_sol(steps[2])
    np.testing.assert_array_equal(res.y, np.hstack([s.y for s in steps]))
    np.testing.assert_array_equal(res.t, np.concatenate([s.t for s in steps]))


def test_in_place_edit_persists_across_append(cell: UnitCell) -> None:
    res = _empty(cell)
    n_sp = len(cell.env.initial_values())
    res.append_sol(_step(0.0, 3, n_sp))
    res.append_sol(_step(1.0, 3, n_sp))
    res.y[0, 0] = -42.0
    res.append_sol(_step(2.0, 3, n_sp))
    assert res.y[0, 0] == -42.0


def test_empty_state_and_assignment(cell: UnitCell) -> None:
    res = _empty(cell)
    assert len(res.y) == 0 and res.t.size == 0
    n_sp = len(cell.env.initial_values())
    s = _step(0.0, 6, n_sp)
    res.t = s.t
    res.y = s.y
    np.testing.assert_array_equal(res.t, s.t)
    np.testing.assert_array_equal(res.y, s.y)
    res.y = np.array([], dtype=np.float64)
    assert len(res.y) == 0


def test_final_cc_without_join(cell: UnitCell) -> None:
    res = _empty(cell)
    n_sp = len(cell.env.initial_values())
    res.append_sol(_step(0.0, 3, n_sp))
    last = _step(1.0, 3, n_sp)
    res.append_sol(last)
    res.append_sol(DummyScipyOdeResult(t=np.array([]), y=np.empty((n_sp, 0))))
    out = res.final_cc
    np.testing.assert_array_equal(out, last.y[:, -1])
    assert len(res._y_blocks) == 3  # not joined
    out[0] = 1e9  # a copy
    assert res.y[0, 5] != 1e9


def test_add_leaves_operands_unchanged(cell: UnitCell) -> None:
    n_sp = len(cell.env.initial_values())
    a, b = _empty(cell), _empty(cell)
    a.append_sol(_step(0.0, 3, n_sp))
    b.append_sol(_step(2.0, 4, n_sp))
    ya, yb = a.y.copy(), b.y.copy()
    c = a + b
    assert c.y.shape == (n_sp, 7)
    c.append_sol(_step(10.0, 2, n_sp))
    np.testing.assert_array_equal(a.y, ya)
    np.testing.assert_array_equal(b.y, yb)
    assert len(a.n_points) == 1


def test_shallow_copy_owns_its_blocks(cell: UnitCell) -> None:
    n_sp = len(cell.env.initial_values())
    a = _empty(cell)
    a.append_sol(_step(0.0, 3, n_sp))
    b = copy.copy(a)
    b.append_sol(_step(1.0, 3, n_sp))
    assert a.y.shape == (n_sp, 3)
    assert b.y.shape == (n_sp, 6)


def test_chunked_run_matches_joined_history(cell: UnitCell) -> None:
    it = cell.prepare_chunked_run(t_span=(0.0, 3e-4))
    res = it.run()
    assert len(res._y_blocks) == len(it.chunks)
    # Every chunk starts from the last point of the previous one (up to
    # the initial burst of the ON chunks, which only adds).
    first = 0
    joined = res.y
    for n_prev, chunk in zip(res.n_points[:-1], it.chunks[1:]):
        first += n_prev
        if not chunk.beam_on:
            np.testing.assert_array_equal(joined[:, first], joined[:, first - 1])
    np.testing.assert_array_equal(res.final_cc, joined[:, -1])
