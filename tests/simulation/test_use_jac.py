#!/usr/bin/env python3
"""
Tests of ``use_jac``: the analytic Jacobian is passed to solve_ivp for
LSODA, BDF and Radau only, and the results agree with those obtained
without it.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import warnings

import numpy as np
import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.physics.beam import ConstantBeam
from radiopyo.simulation.unit_cell import RunOptions, UnitCell

# ------------------------------- Fixtures ---------------------------------- #
T_SPAN = (0.0, 1e-5)
RTOL = 1e-4


@pytest.fixture(scope="module")
def cell() -> UnitCell:
    out = UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    if len(out.beam_list) == 0:
        out.set_beam(ConstantBeam(dose_rate=1e7, max_dose=10.0), LET=0.201)
    return out


def assert_close(a: np.ndarray, b: np.ndarray) -> None:
    scale = max(float(np.max(np.abs(a))), 1e-30)
    np.testing.assert_allclose(a, b, rtol=RTOL, atol=1e-6 * scale)


# --------------------------------- Tests ----------------------------------- #
def test_use_jac_default_is_false() -> None:
    assert RunOptions().use_jac is False


@pytest.mark.parametrize("method", ["LSODA", "BDF", "Radau"])
def test_use_jac_agrees_with_no_jac(cell: UnitCell, method: str) -> None:
    ref = cell.run(t_span=T_SPAN, method=method, use_jac=False)
    res = cell.run(t_span=T_SPAN, method=method, use_jac=True)
    assert_close(res.y[:, -1], ref.y[:, -1])


def test_use_jac_chunked_agrees_with_no_jac(cell: UnitCell) -> None:
    ref = cell.prepare_chunked_run(t_span=T_SPAN, use_jac=False).run()
    res = cell.prepare_chunked_run(t_span=T_SPAN, use_jac=True).run()
    assert_close(res.y[:, -1], ref.y[:, -1])


def test_use_jac_rk45_no_warning(cell: UnitCell) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = cell.run(t_span=T_SPAN, method="RK45", use_jac=True)
    assert all(res.success)
