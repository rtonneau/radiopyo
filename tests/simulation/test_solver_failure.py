#!/usr/bin/env python3
"""
Tests of ``SolverFailedError``: a solve_ivp result with success=False is
turned into an exception by UnitCell runs; successful runs are unchanged.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import typing as tp

import pytest
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.exceptions import RadioPyoError
from radiopyo.physics.beam import ConstantBeam
from radiopyo.simulation import unit_cell as uc_mod
from radiopyo.simulation.exceptions import SimulationError, SolverFailedError
from radiopyo.simulation.unit_cell import UnitCell

# ------------------------------- Fixtures ---------------------------------- #
T_SPAN = (0.0, 1e-5)
FAIL_MESSAGE = "Required step size is less than spacing between numbers."


@pytest.fixture(scope="module")
def cell() -> UnitCell:
    out = UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    if len(out.beam_list) == 0:
        out.set_beam(ConstantBeam(dose_rate=1e7, max_dose=10.0), LET=0.201)
    return out


@pytest.fixture
def failing_solver(monkeypatch: pytest.MonkeyPatch) -> None:
    """Run the real solver, then report a failure as scipy does."""

    def fake(*args: tp.Any, **kwargs: tp.Any) -> tp.Any:
        sol = solve_ivp(*args, **kwargs)
        sol.success = False
        sol.status = -1
        sol.message = FAIL_MESSAGE
        return sol

    monkeypatch.setattr(uc_mod, "solve_ivp", fake)


# --------------------------------- Tests ----------------------------------- #
def test_hierarchy() -> None:
    assert issubclass(SolverFailedError, SimulationError)
    assert issubclass(SolverFailedError, RadioPyoError)


@pytest.mark.usefixtures("failing_solver")
def test_run_raises(cell: UnitCell) -> None:
    with pytest.raises(SolverFailedError) as info:
        cell.run(t_span=T_SPAN)
    msg = str(info.value)
    assert FAIL_MESSAGE in msg
    assert "t_span=" in msg
    assert "last reached t=" in msg


@pytest.mark.usefixtures("failing_solver")
def test_chunked_run_raises(cell: UnitCell) -> None:
    with pytest.raises(SolverFailedError):
        cell.prepare_chunked_run(t_span=T_SPAN).run()


def test_successful_run_does_not_raise(cell: UnitCell) -> None:
    res = cell.run(t_span=T_SPAN)
    assert res.y.shape[1] > 0
