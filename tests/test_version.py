"""The version is exposed and logged before a config file is loaded."""
from __future__ import annotations

import logging
import tomllib
from pathlib import Path

import pytest

import radiopyo as rp


def test_version_matches_pyproject() -> None:
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    version = tomllib.loads(pyproject.read_text(encoding="utf-8"))["project"]["version"]
    assert rp.__version__ == version


def test_version_logged_before_config(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO, logger="radiopyo"):
        rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    messages = [r.getMessage() for r in caplog.records]
    i_version = messages.index(f"radiopyo version {rp.__version__}")
    i_loading = next(i for i, m in enumerate(messages) if m.startswith("Loading config"))
    assert i_version == i_loading - 1
