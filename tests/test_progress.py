#!/usr/bin/env python3
"""
Tests of the optional tqdm progress bar (radiopyo.utils.progress) and of
the ``progress=`` keyword of SimChunksIter.run.

A fake ``tqdm`` module stands in for the real one, so the tests neither
depend on tqdm being installed nor print a bar; ``sys.modules["tqdm"] =
None`` makes the import fail as when tqdm is not installed.

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import sys
import types
import typing as tp

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo.physics.beam import PulsedBeam
from radiopyo.simulation.unit_cell import UnitCell
from radiopyo.utils import progress


# ------------------------------- Fixtures ---------------------------------- #
@pytest.fixture(autouse=True)
def _restore_global_setting() -> tp.Iterator[None]:
    saved = progress.progress_bar_enabled()
    yield
    progress.set_progress_bar(saved)


@pytest.fixture
def fake_tqdm(monkeypatch: pytest.MonkeyPatch) -> list[int | None]:
    """Install a fake tqdm module; returns the list of `total` it got."""
    calls: list[int | None] = []

    def tqdm(iterable: tp.Iterable[tp.Any], total: int | None = None) -> tp.Any:
        calls.append(total)
        return iterable

    module = types.ModuleType("tqdm")
    module.tqdm = tqdm  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "tqdm", module)
    return calls


@pytest.fixture
def no_tqdm(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "tqdm", None)  # import raises ImportError


@pytest.fixture(scope="module")
def _pulsed_template() -> UnitCell:
    cell = UnitCell.from_toml(rp.PURE_WATER_CONFIG_FILE)
    beam = PulsedBeam.from_dose_per_pulse(
        pulse_dose=1.0, period=1e-3, on_time=1e-6, n_pulse=2.0
    )
    cell.set_beam(beam, LET=10.0)
    return cell


# --------------------------------- Tests ----------------------------------- #
def test_exported_from_radiopyo() -> None:
    assert rp.set_progress_bar is progress.set_progress_bar
    assert rp.progress_bar_enabled is progress.progress_bar_enabled
    assert progress.progress_bar_enabled() is True  # on by default


def test_default_uses_tqdm(fake_tqdm: list[int | None]) -> None:
    data = [1, 2, 3]
    assert list(progress.with_progress_bar(data, total=3)) == data
    assert fake_tqdm == [3]


def test_global_switch_off(fake_tqdm: list[int | None]) -> None:
    rp.set_progress_bar(False)
    assert rp.progress_bar_enabled() is False
    data = [1, 2]
    assert progress.with_progress_bar(data) is data
    assert fake_tqdm == []


@pytest.mark.parametrize("global_on", [True, False])
def test_keyword_overrides_global(
    fake_tqdm: list[int | None], global_on: bool
) -> None:
    rp.set_progress_bar(global_on)
    data = [1]
    progress.with_progress_bar(data, progress=True)
    assert len(fake_tqdm) == 1
    assert progress.with_progress_bar(data, progress=False) is data
    assert len(fake_tqdm) == 1


def test_off_does_not_import_tqdm(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(*args: tp.Any, **kwargs: tp.Any) -> tp.Any:
        raise AssertionError("tqdm imported while the bar is off")

    monkeypatch.setattr("builtins.__import__", fail)
    data = [1]
    assert progress.with_progress_bar(data, progress=False) is data


@pytest.mark.parametrize("flag", [None, True, False])
def test_without_tqdm(no_tqdm: None, flag: bool | None) -> None:
    data = [1, 2]
    assert progress.with_progress_bar(data, progress=flag) is data


def test_jupyter_uses_tqdm_notebook(monkeypatch: pytest.MonkeyPatch) -> None:
    used: list[str] = []
    std = types.ModuleType("tqdm")
    std.tqdm = lambda it, total=None: used.append("std") or it  # type: ignore[attr-defined]
    nb = types.ModuleType("tqdm.notebook")
    nb.tqdm = lambda it, total=None: used.append("notebook") or it  # type: ignore[attr-defined]
    std.notebook = nb  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "tqdm", std)
    monkeypatch.setitem(sys.modules, "tqdm.notebook", nb)
    monkeypatch.setattr(progress, "_in_jupyter", lambda: True)
    progress.with_progress_bar([1])
    assert used == ["notebook"]


def test_in_jupyter_false_outside_ipython() -> None:
    assert progress._in_jupyter() is False


@pytest.mark.parametrize(
    ("global_on", "flag", "expected"),
    [
        (True, None, True),
        (False, None, False),
        (False, True, True),
        (True, False, False),
    ],
)
def test_sim_chunks_iter_run(
    _pulsed_template: UnitCell,
    fake_tqdm: list[int | None],
    global_on: bool,
    flag: bool | None,
    expected: bool,
) -> None:
    rp.set_progress_bar(global_on)
    chunks = _pulsed_template.prepare_chunked_run(t_span=(0.0, 2e-3))
    res = chunks.run(progress=flag)
    assert res.t[-1] == pytest.approx(2e-3)
    assert fake_tqdm == ([len(chunks)] if expected else [])


def test_sim_chunks_iter_run_without_tqdm(
    _pulsed_template: UnitCell, no_tqdm: None
) -> None:
    res = _pulsed_template.prepare_chunked_run(t_span=(0.0, 2e-3)).run()
    assert res.t[-1] == pytest.approx(2e-3)
