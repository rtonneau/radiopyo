"""The version is exposed and logged before a config file is loaded."""
from __future__ import annotations

import logging
from importlib import metadata

import pytest

import radiopyo as rp


def test_version_matches_metadata() -> None:
    assert rp.__version__ == metadata.version("radiopyo")


def test_version_logged_before_config(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.INFO, logger="radiopyo"):
        rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    messages = [r.getMessage() for r in caplog.records]
    i_version = messages.index(f"radiopyo version {rp.__version__}")
    i_loading = next(i for i, m in enumerate(messages) if m.startswith("Loading config"))
    assert i_version == i_loading - 1
