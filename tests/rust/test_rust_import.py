#!/usr/bin/env python3
"""
The Rust backend without the radiopyo-rs extension: ImportError with the
install hint (runs whether or not the extension is installed).

Author:
    Romain Tonneau (romain.tonneau@unamur.be)
"""

# ----------------------------- External imports ---------------------------- #
from __future__ import annotations

import sys

import pytest

# ----------------------------- Internal imports ---------------------------- #
import radiopyo as rp
from radiopyo import rust


def test_missing_extension_raises_import_error_with_hint(monkeypatch) -> None:  # type: ignore[no-untyped-def]
    monkeypatch.setitem(sys.modules, "radiopyo_rs", None)
    monkeypatch.setitem(sys.modules, "radiopyo_rs._core", None)
    assert not rust.is_available()
    cell = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    with pytest.raises(ImportError, match=r"pip install radiopyo\[rust\]") as info:
        cell.run((0.0, 1e-3), backend="rust")
    assert "pip install ./rust" in str(info.value)


def test_unknown_backend_raises() -> None:
    cell = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    with pytest.raises(ValueError, match="backend"):
        cell.run((0.0, 1e-3), backend="fortran")


def test_config_is_stored_and_deepcopied() -> None:
    import copy

    cell = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    assert cell.config is not None and "reactions" in cell.config
    other = copy.deepcopy(cell)
    assert other.config == cell.config and other.config is not cell.config
    assert rp.UnitCell(env=cell.env).config is None
