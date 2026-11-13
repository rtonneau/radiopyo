#!/usr/bin/env python3
"""
Rust backend (backend="rust", radiopyo.rust, ADR 0010), plain run mode:
same results as the Python backend, parity guard, error mapping, and
the ODEResult features on a Rust result. Skipped without radiopyo-rs.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import dataclasses as dc
import typing as tp

import numpy as np
import pytest

pytest.importorskip("radiopyo_rs")

# ----------------------------- Internal imports ---------------------------- #
from radiopyo import rust
from radiopyo.physics.beam import (
    ConstantBeam,
    MacroPulsedBeam,
    PulsedBeam,
    SinglePulseBeam,
)
from radiopyo.simulation.exceptions import (
    NoBeamDefinedError,
    ParityError,
    RustBackendError,
)
from radiopyo.simulation.persistence import H5Manager, SimIdentifier
from radiopyo.simulation.unit_cell import UnitCell

if tp.TYPE_CHECKING:
    from radiopyo.physics.beam import ParticleBeam
    from radiopyo.simulation.ode_result import ODEResult

# Python vs radkin tolerance on concentrations [µmol/L].
RTOL = 1e-6
ATOL = 1e-10
LET = 0.201  # keV/µm


_Case = tp.Tuple["ParticleBeam", tp.Tuple[float, float], tp.Optional[float]]


def _beams() -> tp.Dict[str, _Case]:
    """name -> (beam, t_span, max_step). A plain run reads the beam at
    every evaluation: max_step below the ON time keeps LSODA from
    stepping over a pulse (otherwise both backends miss pulses, each in
    its own way)."""
    return {
        "constant": (ConstantBeam(dose_rate=0.1, max_dose=1e-3), (0.0, 0.1), None),
        "single_pulse": (
            SinglePulseBeam(dose_rate=1e6, max_dose=10.0),
            (0.0, 1e-3),
            None,
        ),
        "pulsed": (
            PulsedBeam(dose_rate_peak=1e6, period=1e-3, on_time=2e-6, max_dose=6.0),
            (0.0, 3e-3),
            1e-6,
        ),
        "macro_pulsed": (
            MacroPulsedBeam(
                dose_rate_peak=1e6,
                period=1e-3,
                on_time=2e-5,
                micro_period=1e-5,
                micro_on_time=5e-6,
                max_dose=40.0,
            ),
            (0.0, 2e-3),
            2.5e-6,
        ),
        "from_dose_rate_peak": (
            PulsedBeam.from_dose_rate_peak(1e6, 1e-3, 2e-6, 3, start_time=0.0),
            (0.0, 3e-3),
            1e-6,
        ),
    }


def _assert_same(py: ODEResult, rs: ODEResult) -> None:
    np.testing.assert_allclose(rs.final_cc, py.final_cc, rtol=RTOL, atol=ATOL)


def _both(cell: UnitCell, **kwargs: tp.Any) -> tp.Tuple[ODEResult, ODEResult]:
    return cell.run(**kwargs), cell.run(backend="rust", **kwargs)


@pytest.mark.parametrize("config", ["pure_water", "labarbe"])
@pytest.mark.parametrize("beam_name", list(_beams()))
def test_run_matches_python(config: str, beam_name: str, request) -> None:  # type: ignore[no-untyped-def]
    cell = request.getfixturevalue(config)
    beam, t_span, max_step = _beams()[beam_name]
    cell.set_beam(beam, LET=LET)
    py, rs = _both(cell, t_span=t_span, max_step=max_step)
    _assert_same(py, rs)
    assert rs.burst_dose_rates == pytest.approx(py.burst_dose_rates)


def test_paper_setup_matches_python(pure_water: UnitCell) -> None:
    """GF(DR) on and O2 edited on the env after the build (paper notebook)."""
    pure_water.env.update_flag_GValues(True)
    pure_water.env.species["O2"].set_initial_cc(50e-6)
    pure_water.set_beam(SinglePulseBeam(dose_rate=1e6, max_dose=10.0), LET=1.17)
    _assert_same(*_both(pure_water, t_span=(0.0, 1e-3)))


def test_t_eval(labarbe: UnitCell) -> None:
    t_eval = np.geomspace(1e-6, 1e-2, 30)
    py, rs = _both(labarbe, t_span=(0.0, 1e-2), t_eval=t_eval)
    np.testing.assert_allclose(rs.t, py.t, rtol=1e-15)
    np.testing.assert_allclose(rs.y, py.y, rtol=RTOL, atol=ATOL)


def test_temperature_and_let_changes(labarbe: UnitCell) -> None:
    labarbe.set_temperature(37.0)
    labarbe.change_LET(2.5)
    _assert_same(*_both(labarbe, t_span=(0.0, 1e-3)))


def test_y0_and_no_beam_use(labarbe: UnitCell) -> None:
    y0 = labarbe.env.initial_values() * 1.5
    y0_copy = y0.copy()
    _assert_same(*_both(labarbe, t_span=(0.0, 1e-3), y0=y0))
    np.testing.assert_array_equal(y0, y0_copy)
    _assert_same(*_both(labarbe, t_span=(0.0, 1e-3), use_beam=False))


def test_run_till_max_dose(labarbe: UnitCell) -> None:
    py = labarbe.run_till_max_dose()
    rs = labarbe.run_till_max_dose(backend="rust")
    _assert_same(py, rs)
    assert rs.short_on_spans == py.short_on_spans


def test_parity_guard_raises_after_editing_a_k(labarbe: UnitCell) -> None:
    reaction = next(iter(labarbe.env.reactions.k_reaction))
    reaction.set_k_value(reaction.k_value(labarbe._T) * 2.0)
    labarbe.update_env()
    with pytest.raises(ParityError, match="k values"):
        labarbe.run((0.0, 1e-3), backend="rust")


def test_cell_without_config_raises(labarbe: UnitCell) -> None:
    cell = UnitCell(env=labarbe.env, beam=labarbe.beam)
    with pytest.raises(RustBackendError, match="config"):
        cell.run((0.0, 1e-3), backend="rust")


def test_other_method_raises(labarbe: UnitCell) -> None:
    with pytest.raises(RustBackendError, match="LSODA"):
        labarbe.run((0.0, 1e-3), method="BDF", backend="rust")


def test_no_beam_raises(labarbe: UnitCell) -> None:
    cell = UnitCell(env=labarbe.env)
    cell.config = labarbe.config
    with pytest.raises(NoBeamDefinedError):
        cell.run((0.0, 1e-3), backend="rust")


def test_result_features(labarbe: UnitCell, tmp_path) -> None:  # type: ignore[no-untyped-def]
    py, rs = _both(labarbe, t_span=(0.0, 1e-3))
    df = rs.to_pandas()
    assert list(df.columns) == list(py.to_pandas().columns)
    rates = rs.reconstruct_reactions()
    assert rates is not None
    assert rs.nfev and rs.success == [True]

    @dc.dataclass(kw_only=True)
    class _ID(SimIdentifier):
        case: str = "rust"

    hdf = H5Manager.from_path(tmp_path / "rust.h5")
    hdf.push_ODEResult(_ID(), rs)
    loaded = hdf.read_ODEResult(_ID())
    np.testing.assert_allclose(loaded.y[:, -1], rs.y[:, -1], rtol=1e-12)


def test_rust_module_helpers() -> None:
    assert rust.is_available()
